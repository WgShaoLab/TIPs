import os
import re
import time
import pandas as pd
from typing import List
from Bio import SeqIO
from .utils import run_cmd


class FastaFile:
    """
    Utility class to:
      - normalize FASTA into single-line sequences
      - enforce PE= tags
      - append custom tags in headers
      - split into sub-FASTAs for parallel blastp
      - generate decoys
    """

    def __init__(
        self,
        fasta_path: str,
        uniform_pe_tag: str | None = None,
        extra_header_tag: str | None = None,
    ):
        self.fasta_path = fasta_path
        with open(self.fasta_path, "r") as handle:
            self.fasta_str = handle.read()

        # Count records
        self.record_count = len(self.fasta_str.split(">")[1:])

        # If multiline FASTA (seq wrapped across many lines),
        # convert to single-line-per-entry
        if len(self.fasta_str.split("\n")) > self.record_count * 2 + 4:
            self._to_single_line()

        # Normalize PE tag if requested
        if uniform_pe_tag:
            self._uniform_pe_tag(uniform_pe_tag)

        # Append extra header tag if requested
        if extra_header_tag:
            self._add_header_tag(extra_header_tag)

        # Collect PE levels present
        pe_list = re.findall(r"PE=(\d)", self.fasta_str)
        self.pe_levels = list(sorted(set(pe_list)))

    def _to_single_line(self):
        """
        Force each FASTA entry to be of the form:
            >header
            SEQUENCE(no newlines)
        """
        out_lines = []
        for record in self.fasta_str.split(">")[1:]:
            parts = record.strip().split("\n")
            header = parts[0]
            seq = "".join(parts[1:])
            out_lines.append(f">{header}\n{seq}\n")
        self.fasta_str = "".join(out_lines)
        print("[FastaFile] Converted multi-line FASTA to single-line entries")

    def _uniform_pe_tag(self, tag: str):
        """
        Ensure every header has a 'PE={tag}' annotation.
        If PE exists, replace it; else append it.
        """
        new_records = []
        for record in self.fasta_str.split(">")[1:]:
            header, seq = record.strip().split("\n", 1)
            if "PE=" in header:
                # Replace: PE=digit (surrounded by spaces)
                header = re.sub(r"PE=(\d)", f"PE={tag}", header)
            else:
                header = header + f" PE={tag}"
            # keep sequence as a single block
            seq = seq.replace("\n", "")
            new_records.append(f">{header}\n{seq}\n")
        self.fasta_str = "".join(new_records)
        print("[FastaFile] Unified PE tag")

    def _add_header_tag(self, tag: str):
        """
        Append arbitrary tag text to each FASTA header.
        e.g. to label TE-derived source, contaminants, etc.
        """
        new_records = []
        for record in self.fasta_str.split(">")[1:]:
            header, seq = record.strip().split("\n", 1)
            header = header + f" {tag}"
            seq = seq.replace("\n", "")
            new_records.append(f">{header}\n{seq}\n")
        self.fasta_str = "".join(new_records)
        print(f"[FastaFile] Added header tag '{tag}'")

    def save(self, out_fasta: str):
        """
        Write current FASTA content to file.
        """
        with open(out_fasta, "w") as handle:
            handle.write(self.fasta_str)

    def split_for_parallel(self, out_dir: str, num_files: int = 15) -> List[str]:
        """
        Split current FASTA into N chunks for parallel blastp on short peptides.
        Returns the list of created FASTA paths.
        """
        os.makedirs(out_dir, exist_ok=True)
        records = self.fasta_str.split(">")[1:]
        per_file = len(records) // num_files
        remainder = len(records) % num_files

        paths = []
        start = 0
        for i in range(num_files):
            end = start + per_file + (1 if i < remainder else 0)
            sub_records = records[start:end]
            start = end

            sub_path = os.path.join(out_dir, f"sub_fasta_{i+1}.fasta")
            with open(sub_path, "w") as sub_handle:
                sub_handle.write("".join([">" + r for r in sub_records]))
            paths.append(sub_path)

        print(f"[FastaFile] Split into {len(paths)} FASTA chunks under {out_dir}")
        return paths

    @staticmethod
    def generate_decoy(input_fasta: str, out_fasta: str):
        """
        Reverse each sequence to produce a decoy entry (Target-Decoy strategy).
        Decoys are prefixed with 'rev_'.
        """
        with open(out_fasta, "w") as out_handle:
            for rec in SeqIO.parse(input_fasta, "fasta"):
                decoy_seq = str(rec.seq)[::-1]
                out_handle.write(f">rev_{rec.id}\n{decoy_seq}\n")

    @staticmethod
    def blastp_short_parallel(
        fasta_chunks_dir: str,
        blastp_bin: str,
        blast_db: str,
        threads_per_chunk: int = 4,
        poll_interval: int = 5,
    ) -> str:
        """
        Run blastp-short on each split FASTA file in parallel (background &),
        then merge and filter.
        Returns path to merged+filtered TSV.
        NOTE: This function assumes GNU tools (nohup, awk, etc.).
        """

        # launch blastp for each sub_fasta
        for fname in os.listdir(fasta_chunks_dir):
            if not fname.startswith("sub_fasta_") or not fname.endswith(".fasta"):
                continue
            sub_q = os.path.join(fasta_chunks_dir, fname)
            sub_out = sub_q.replace(".fasta", ".txt")

            cmd = (
                f"nohup {blastp_bin} -task blastp-short "
                f"-query {sub_q} "
                f"-db {blast_db} "
                f"-out {sub_out} "
                f"-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend "
                f"sstart send evalue bitscore qseq sseq' "
                f"-evalue 20000 -num_threads {threads_per_chunk} &"
            )
            run_cmd(cmd)

        # wait for all blastp to stop
        while True:
            try:
                out = int(
                    subprocess.check_output(
                        "ps aux | grep -v grep | grep -c 'blastp'",
                        shell=True,
                    ).strip()
                )
            except Exception:
                out = 0
            if out <= 0:
                break
            time.sleep(poll_interval)

        print("[FastaFile] blastp-short done for all chunks")

        merged_txt = os.path.join(fasta_chunks_dir, "merged_blastp.txt")
        filtered_txt = os.path.join(fasta_chunks_dir, "merged_blastp_filter.txt")

        # merge
        run_cmd(f"cat {fasta_chunks_dir}/*.txt > {merged_txt}")

        # filter: qstart==1, gapopen==0, pident>=75
        # (mirrors your awk logic)
        run_cmd(
            "awk '$7 == 1 && $6 == 0 && $3 >= 75' "
            f"{merged_txt} > {filtered_txt}"
        )

        return filtered_txt

    @staticmethod
    def load_and_filter_blastp_results(filtered_txt: str) -> pd.DataFrame:
        """
        Convert merged_blastp_filter.txt into a dataframe,
        compute I/L-equivalent identity, and retain full-length perfect matches.
        """
        cols = [
            "qaccver", "saccver", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "qseq", "sseq"
        ]
        df = pd.read_csv(filtered_txt, sep="\t", names=cols)
        df["q_length"] = df["qaccver"].str.len()
        df["qseq_I2L"] = df["qseq"].str.replace("I", "L")
        df["sseq_I2L"] = df["sseq"].str.replace("I", "L")

        df = df[
            (df["length"] == df["q_length"]) &
            (df["qseq"].str.len() == df["sseq"].str.len()) &
            (df["gapopen"] == 0)
        ]

        df["IL_identity"] = (df["qseq_I2L"] == df["sseq_I2L"])
        return df
