import os
import glob
from typing import List, Optional
from Bio import SeqIO
import os, glob
from .utils import run_cmd, timing
from .fasta_utils import FastaFile
from .ms_postprocess import (
    filter_by_fdr,
    blastp_remove_canonical,
    predict_binding_affinity_mixmhcpred,
)


class DBSearchPipeline:
    """
    Pipeline steps:
      1) Build a search FASTA for each engine:
         - Human proteome (PE=0) + TE subdatabase (PE=1) + contaminants (PE=0) [+ decoy]
         - Prefer a TE subdatabase selected by 'de novo tag → TE reference (blastp-short) → extract hits'
           If there are no hits, fall back to the de novo merged FASTA itself as the TE candidate set.
      2) Run Comet / MSFragger / MS-GF+ (if enabled)
      3) (Optional) iProphet integration, FDR filtering, homology removal, and HLA binding prediction
    """

    def __init__(self, cfg):
        self.cfg = cfg
        self.sample_path = cfg.sample_path
        self.sample_name = cfg.raw["sample"]["name"]

        self.db_root = os.path.join(self.sample_path, "DB_search_iProphet")
        os.makedirs(self.db_root, exist_ok=True)

    # -------------------------- helpers --------------------------

    def _denovo_merged_fasta(self) -> str:
        """Get the de novo merged FASTA (written by stage_denovo)"""
        path = self.cfg.raw["database_build"].get("denovo_te_tagged_fasta", "AUTO")
        if path == "AUTO":
            return os.path.join(
                self.sample_path, "Denovo", "Denovo_TE_SoftMerged_InstaNovo.fasta"
            )
        return path

    @staticmethod
    def _write_concat(fastas: List[str], out_path: str) -> None:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as oh:
            for fa in fastas:
                with open(fa, "r") as ih:
                    oh.write(ih.read())

    def _normalize_fasta(self, in_fa: str, pe_tag: str, extra_tag: Optional[str] = None) -> str:
        """Single-line entries + unify PE= tag + optional source tag; return a temporary normalized path"""
        tmp_dir = os.path.join(self.db_root, "tmp_norm")
        os.makedirs(tmp_dir, exist_ok=True)
        base = os.path.basename(in_fa)
        out_path = os.path.join(tmp_dir, f"PE{pe_tag}__{base}")
        ff = FastaFile(in_fa, uniform_pe_tag=pe_tag, extra_header_tag=extra_tag)
        ff.save(out_path)
        return out_path

    def _make_decoy_inline(self, target_fa: str) -> None:
        """Generate reverse decoys from target_fa and append them to the same file"""
        decoy_path = target_fa + ".decoy.tmp.fasta"
        FastaFile.generate_decoy(target_fa, decoy_path)
        with open(target_fa, "a") as oh, open(decoy_path, "r") as ih:
            oh.write(ih.read())
        os.remove(decoy_path)

    # ---------------------- TE selection (blast) ----------------------

    def _build_te_candidates_via_blast(self, raw_tag_fasta: str) -> Optional[str]:
        """
        Use de novo tags (short peptides) to run blastp-short against the TE reference DB and extract hit TE entries to build a sample-level TE subdatabase.
        Return the subdatabase path; return None if there are no hits.
        """
        db_cfg = self.cfg.raw.get("database_build", {})
        te_db = db_cfg.get("te_db_fasta")
        sel = db_cfg.get("te_selection", {})
        if not te_db or sel.get("method", "direct").lower() != "blast":
            return None

        blast_cfg = self.cfg.raw.get("blast", {})
        blastp_bin = blast_cfg.get("blastp_bin")
        if not blastp_bin:
            raise RuntimeError("Configuration is missing blastp_bin (under .blast.blastp_bin).")

        chunks = int(sel.get("chunks", 15))
        threads = int(sel.get("threads_per_chunk", 4))

        work_dir = os.path.join(self.sample_path, "Denovo", "TE_blast_work")
        os.makedirs(work_dir, exist_ok=True)

        # 1) [[EN REQUIRED]] FASTA [[EN REQUIRED]] blastp-short
        q = FastaFile(raw_tag_fasta)
        q.split_for_parallel(work_dir, num_files=chunks)

        filtered_txt = FastaFile.blastp_short_parallel(
            fasta_chunks_dir=work_dir,
            blastp_bin=blastp_bin,
            blast_db=te_db,
            threads_per_chunk=threads,
            poll_interval=5,
        )

        # 2) [[EN REQUIRED]] I/L [[EN REQUIRED]] + [[EN REQUIRED]] + gapopen=0 [[EN REQUIRED]]
        df = FastaFile.load_and_filter_blastp_results(filtered_txt)
        if df is None or len(df) == 0:
            return None

        # [[EN REQUIRED]] I/L [[EN REQUIRED]]
        if "IL_identity" in df.columns:
            df = df[df["IL_identity"] == True]
        if len(df) == 0:
            return None

        te_ids = sorted(set(df["saccver"].tolist()))
        if not te_ids:
            return None

        # 3) [[EN REQUIRED]] ID [[EN REQUIRED]] TE [[EN REQUIRED]]
        out_te_fa = os.path.join(self.sample_path, "Denovo", "TE_candidates_from_DB.fasta")

        # [[EN REQUIRED]] FASTA [[EN REQUIRED]] BLAST [[EN REQUIRED]]
        def _is_blast_db(base: str) -> bool:
            # , .00.pin/.01.pin... [[EN REQUIRED]] .pal/.pin
            if os.path.exists(base + ".pal") or os.path.exists(base + ".pin"):
                return True
            if glob.glob(base + ".[0-9][0-9].pin"):
                return True
            return False

        out_te_fa = os.path.join(self.sample_path, "Denovo", "TE_candidates_from_DB.fasta")
        if _is_blast_db(te_db):
            # [[EN REQUIRED]] blastdbcmd [[EN REQUIRED]] BLAST [[EN REQUIRED]] ID [[EN REQUIRED]] FASTA
            ids_txt = os.path.join(work_dir, "te_ids.txt")
            with open(ids_txt, "w") as fh:
                fh.write("\n".join(te_ids) + "\n")

            blastdbcmd_bin = blast_cfg.get("blastdbcmd_bin", "blastdbcmd")
            cmd = (
                f'"{blastdbcmd_bin}" -db "{te_db}" -dbtype prot '
                f'-entry_batch "{ids_txt}" -out "{out_te_fa}" -outfmt %f'
            )
            run_cmd(cmd, env=os.environ)
        else:
            # [[EN REQUIRED]] FASTA [[EN REQUIRED]]
            with open(out_te_fa, "w") as oh:
                for rec in SeqIO.parse(te_db, "fasta"):
                    if rec.id in te_ids:
                        oh.write(f">{rec.id}\n{str(rec.seq)}\n")

        print(f"[TE] Extracted {len(te_ids)} TE entries → {out_te_fa}")
        return out_te_fa

    # ---------------------- Stage: build DB ----------------------

    @timing
    def build_search_database(self) -> None:
        """
        Build each engine's DB FASTA:
          human(PE=0) + TE(PE=1) + contaminants(PE=0) [+ decoy]
        Prefer the TE subdatabase selected by 'tag→TE DB (blast)'; if no hits, fall back to the de novo tag FASTA.
        """
        db_cfg = self.cfg.raw["database_build"]

        human_fa = db_cfg["human_fasta"]
        cont_fa = db_cfg["contaminants_fasta"]
        denovo_tag_fa = self._denovo_merged_fasta()

        # 1) [[EN REQUIRED]] TE [[EN REQUIRED]]
        te_subdb = self._build_te_candidates_via_blast(denovo_tag_fa)

        # 2) ( + [[EN REQUIRED]] PE + )
        human_norm = self._normalize_fasta(human_fa, pe_tag="0", extra_tag="SRC=HUMAN")
        cont_norm  = self._normalize_fasta(cont_fa,  pe_tag="0", extra_tag="SRC=CONTAM")
        te_source  = te_subdb if te_subdb else denovo_tag_fa
        te_norm    = self._normalize_fasta(te_source, pe_tag="1", extra_tag="SRC=TE")

        # 3) [[EN REQUIRED]] DB [[EN REQUIRED]] decoy
        add_decoy = bool(db_cfg.get("add_decoy", True))
        engines = self.cfg.raw.get("search_engines", {})

        for eng in ("comet", "msfragger", "msgfplus"):
            eng_cfg = engines.get(eng, {})
            if not eng_cfg or not eng_cfg.get("enable", False):
                continue

            eng_dir = os.path.join(self.db_root, eng.upper())
            os.makedirs(eng_dir, exist_ok=True)
            db_fa = os.path.join(eng_dir, f"{eng.upper()}_DBsearch.fasta")

            # [[EN REQUIRED]]
            self._write_concat([human_norm, te_norm, cont_norm], db_fa)

            # [[EN REQUIRED]] decoy
            if add_decoy:
                self._make_decoy_inline(db_fa)

            print(f"[DB] {eng} database ready:{db_fa}")

    # ---------------------- Stage: run engines ----------------------

    def _mzml_files(self) -> List[str]:
        return sorted(glob.glob(os.path.join(self.sample_path, "mzML", "*.mzML")))

    @timing
    def run_comet(self) -> None:
        if not self.cfg.raw["search_engines"]["comet"]["enable"]:
            return
        eng_dir = os.path.join(self.db_root, "COMET")
        db_fa   = os.path.join(eng_dir, "COMET_DBsearch.fasta")
        params  = self.cfg.raw["search_engines"]["comet"]["params_file"]
        comet   = self.cfg.raw["search_engines"]["comet"]["binary"]

        for mzml in self._mzml_files():
            cmd = f'{comet} -P{params} -D{db_fa} "{mzml}"'
            run_cmd(cmd, env=os.environ)

    @timing
    def run_msfragger(self) -> None:
        if not self.cfg.raw["search_engines"]["msfragger"]["enable"]:
            return
        eng_dir = os.path.join(self.db_root, "MSFRAGGER")
        db_fa   = os.path.join(eng_dir, "MSFRAGGER_DBsearch.fasta")
        java    = self.cfg.raw["search_engines"]["msfragger"]["java_bin"]
        jar     = self.cfg.raw["search_engines"]["msfragger"]["jar"]
        basepar = self.cfg.raw["search_engines"]["msfragger"]["base_params"]

        for mzml in self._mzml_files():
            cmd = f'{java} -Xmx16G -jar "{jar}" "{basepar}" "{mzml}" --database "{db_fa}"'
            run_cmd(cmd, env=os.environ)

    @timing
    def run_msgfplus(self) -> None:
        if not self.cfg.raw["search_engines"]["msgfplus"]["enable"]:
            return
        eng_dir = os.path.join(self.db_root, "MSGFPLUS")
        db_fa   = os.path.join(eng_dir, "MSGFPLUS_DBsearch.fasta")
        java    = self.cfg.raw["search_engines"]["msgfplus"]["java_bin"]
        jar     = self.cfg.raw["search_engines"]["msgfplus"]["jar"]
        params  = self.cfg.raw["search_engines"]["msgfplus"]["params_file"]

        for mzml in self._mzml_files():
            out_mzid = os.path.splitext(os.path.basename(mzml))[0] + ".mzid"
            cmd = (
                f'{java} -Xmx8G -jar "{jar}" '
                f'-d "{db_fa}" -s "{mzml}" -o "{os.path.join(eng_dir, out_mzid)}" '
                f'-conf "{params}"'
            )
            run_cmd(cmd, env=os.environ)

    # ---------------------- Stage: (optional) postprocess ----------------------

    @timing
    def run_peptideprophet_and_iprophet(self, fdr_threshold: float = 0.03) -> None:
        """
            Placeholder implementation: integrate parsing and combination according to your local TPP/Philosopher installation.
            Below is a typical order of subsequent calls:
              1) filter_by_fdr(...)
              2) blastp_remove_canonical(...)
              3) predict_binding_affinity_mixmhcpred(...)
        """
        print(
            "[run_peptideprophet_and_iprophet] [[EN REQUIRED]] TPP/Philosopher,.\n"
            "[[EN REQUIRED]] iProphet , ms_postprocess.filter_by_fdr() [[EN REQUIRED]] FDR ,\n"
            "[[EN REQUIRED]] blastp_remove_canonical() , predict_binding_affinity_mixmhcpred() [[EN REQUIRED]] HLA ."
        )
