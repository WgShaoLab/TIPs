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
      1) 构建每个引擎的搜库 FASTA：
         - 人类蛋白组(PE=0) + TE 子库(PE=1) + 污染库(PE=0) [+ decoy]
         - TE 子库优先来自 “de novo tag → TE 大库(blastp-short) → 抽取命中条目”
           若没有命中，则退回使用 de novo 合并 FASTA 本体作为 TE 候选。
      2) 运行 Comet / MSFragger / MS-GF+（如启用）
      3) （可选）iProphet 整合、FDR 过滤、去同源、HLA 结合预测
    """

    def __init__(self, cfg):
        self.cfg = cfg
        self.sample_path = cfg.sample_path
        self.sample_name = cfg.raw["sample"]["name"]

        self.db_root = os.path.join(self.sample_path, "DB_search_iProphet")
        os.makedirs(self.db_root, exist_ok=True)

    # -------------------------- helpers --------------------------

    def _denovo_merged_fasta(self) -> str:
        """获取 de novo 合并 FASTA（由 stage_denovo 写入）"""
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
        """单行化 + 统一 PE=标签 + 可选追加来源标签，返回临时标准化路径"""
        tmp_dir = os.path.join(self.db_root, "tmp_norm")
        os.makedirs(tmp_dir, exist_ok=True)
        base = os.path.basename(in_fa)
        out_path = os.path.join(tmp_dir, f"PE{pe_tag}__{base}")
        ff = FastaFile(in_fa, uniform_pe_tag=pe_tag, extra_header_tag=extra_tag)
        ff.save(out_path)
        return out_path

    def _make_decoy_inline(self, target_fa: str) -> None:
        """基于 target_fa 生成反向 decoy 并追加到同一文件尾部"""
        decoy_path = target_fa + ".decoy.tmp.fasta"
        FastaFile.generate_decoy(target_fa, decoy_path)
        with open(target_fa, "a") as oh, open(decoy_path, "r") as ih:
            oh.write(ih.read())
        os.remove(decoy_path)

    # ---------------------- TE selection (blast) ----------------------

    def _build_te_candidates_via_blast(self, raw_tag_fasta: str) -> Optional[str]:
        """
        用 de novo tag（短肽）在 TE 大库上跑 blastp-short，抽取命中的 TE 条目，生成“样本级 TE 子库”。
        返回子库路径；如果没有命中则返回 None。
        """
        db_cfg = self.cfg.raw.get("database_build", {})
        te_db = db_cfg.get("te_db_fasta")
        sel = db_cfg.get("te_selection", {})
        if not te_db or sel.get("method", "direct").lower() != "blast":
            return None

        blast_cfg = self.cfg.raw.get("blast", {})
        blastp_bin = blast_cfg.get("blastp_bin")
        if not blastp_bin:
            raise RuntimeError("配置缺少 blastp_bin（在 .blast.blastp_bin 下）")

        chunks = int(sel.get("chunks", 15))
        threads = int(sel.get("threads_per_chunk", 4))

        work_dir = os.path.join(self.sample_path, "Denovo", "TE_blast_work")
        os.makedirs(work_dir, exist_ok=True)

        # 1) 拆分短肽 FASTA 并行跑 blastp-short
        q = FastaFile(raw_tag_fasta)
        q.split_for_parallel(work_dir, num_files=chunks)

        filtered_txt = FastaFile.blastp_short_parallel(
            fasta_chunks_dir=work_dir,
            blastp_bin=blastp_bin,
            blast_db=te_db,
            threads_per_chunk=threads,
            poll_interval=5,
        )

        # 2) 载入并按 I/L 等价 + 全长 + gapopen=0 过滤
        df = FastaFile.load_and_filter_blastp_results(filtered_txt)
        if df is None or len(df) == 0:
            return None

        # 仅保留 I/L 等价一致
        if "IL_identity" in df.columns:
            df = df[df["IL_identity"] == True]
        if len(df) == 0:
            return None

        te_ids = sorted(set(df["saccver"].tolist()))
        if not te_ids:
            return None

        # 3) 按 ID 从 TE 大库中抽取条目
        out_te_fa = os.path.join(self.sample_path, "Denovo", "TE_candidates_from_DB.fasta")

        # 根据路径判断是 FASTA 还是 BLAST 数据库
        def _is_blast_db(base: str) -> bool:
            # 兼容分卷索引，如 .00.pin/.01.pin... 或存在 .pal/.pin
            if os.path.exists(base + ".pal") or os.path.exists(base + ".pin"):
                return True
            if glob.glob(base + ".[0-9][0-9].pin"):
                return True
            return False

        out_te_fa = os.path.join(self.sample_path, "Denovo", "TE_candidates_from_DB.fasta")
        if _is_blast_db(te_db):
            # 用 blastdbcmd 从 BLAST 库按 ID 批量导出为 FASTA
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
            # 仍然支持直接给 FASTA 的场景
            with open(out_te_fa, "w") as oh:
                for rec in SeqIO.parse(te_db, "fasta"):
                    if rec.id in te_ids:
                        oh.write(f">{rec.id}\n{str(rec.seq)}\n")

        print(f"[TE] 抽取到 {len(te_ids)} 条 TE 条目 → {out_te_fa}")
        return out_te_fa

    # ---------------------- Stage: build DB ----------------------

    @timing
    def build_search_database(self) -> None:
        """
        构建每个引擎的 DB FASTA：
          human(PE=0) + TE(PE=1) + contaminants(PE=0) [+ decoy]
        TE 优先采用 “tag→TE库(blast)” 的子库；无命中则退回 de novo tag FASTA。
        """
        db_cfg = self.cfg.raw["database_build"]

        human_fa = db_cfg["human_fasta"]
        cont_fa = db_cfg["contaminants_fasta"]
        denovo_tag_fa = self._denovo_merged_fasta()

        # 1) 先尝试从 TE 大库中抽取样本级子库
        te_subdb = self._build_te_candidates_via_blast(denovo_tag_fa)

        # 2) 归一化（单行 + 统一 PE + 追加来源标记）
        human_norm = self._normalize_fasta(human_fa, pe_tag="0", extra_tag="SRC=HUMAN")
        cont_norm  = self._normalize_fasta(cont_fa,  pe_tag="0", extra_tag="SRC=CONTAM")
        te_source  = te_subdb if te_subdb else denovo_tag_fa
        te_norm    = self._normalize_fasta(te_source, pe_tag="1", extra_tag="SRC=TE")

        # 3) 合成每个引擎的 DB 并加 decoy
        add_decoy = bool(db_cfg.get("add_decoy", True))
        engines = self.cfg.raw.get("search_engines", {})

        for eng in ("comet", "msfragger", "msgfplus"):
            eng_cfg = engines.get(eng, {})
            if not eng_cfg or not eng_cfg.get("enable", False):
                continue

            eng_dir = os.path.join(self.db_root, eng.upper())
            os.makedirs(eng_dir, exist_ok=True)
            db_fa = os.path.join(eng_dir, f"{eng.upper()}_DBsearch.fasta")

            # 组合目标库
            self._write_concat([human_norm, te_norm, cont_norm], db_fa)

            # 追加 decoy
            if add_decoy:
                self._make_decoy_inline(db_fa)

            print(f"[DB] {eng} 数据库已就绪：{db_fa}")

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
            占位实现：你可以根据本地 TPP/Philosopher 安装把解析与整合接上。
            这里保留典型的后续调用顺序示例：
              1) filter_by_fdr(...)
              2) blastp_remove_canonical(...)
              3) predict_binding_affinity_mixmhcpred(...)
        """
        print(
            "[run_peptideprophet_and_iprophet] 该步骤依赖 TPP/Philosopher，请按你的安装在此调用。\n"
            "完成 iProphet 后，可用 ms_postprocess.filter_by_fdr() 进一步按 FDR 分组过滤，\n"
            "再用 blastp_remove_canonical() 去同源，最后用 predict_binding_affinity_mixmhcpred() 做 HLA 结合预测。"
        )
