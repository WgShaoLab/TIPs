import os
import re
import glob
import pandas as pd
import numpy as np
from io import StringIO
from typing import List
from .utils import run_cmd, timing
from .sequence_filters import extract_confident_subsequences


class DeNovoEngine:
    """
    Base class for de novo sequencing engines.
    Each engine:
      - run()   -> generate raw prediction files in its own work dir
      - parse() -> load & normalize results into a common dataframe:
                   columns: ['Sequence', 'ScoreList'] where ScoreList is List[float]
    """

    def __init__(self, sample_path: str, engine_name: str, gpu_id: str):
        self.sample_path = sample_path
        self.engine_name = engine_name
        self.gpu_id = gpu_id
        self.work_dir = os.path.join(sample_path, "Denovo", engine_name)
        os.makedirs(self.work_dir, exist_ok=True)
        self.mgf_dir = os.path.join(sample_path, "mgf")
        self.mgf_list = sorted(glob.glob(os.path.join(self.mgf_dir, "*.mgf")))

    def run(self):
        raise NotImplementedError

    def parse(self) -> pd.DataFrame:
        raise NotImplementedError

    def filter_subseq(self, min_score: float, min_len: int) -> pd.DataFrame:
        """
        Apply extract_confident_subsequences row-wise.
        Return df with 'Filtered_Subsequences' list for each spectrum.
        """
        df = self.parse()
        df["Filtered_Subsequences"] = df.apply(
            lambda row: extract_confident_subsequences(
                row["Sequence"], row["ScoreList"],
                min_score=min_score,
                min_length=min_len
            ),
            axis=1
        )
        df = df[df["Filtered_Subsequences"].str.len() > 0]
        return df


class CasanovoEngine(DeNovoEngine):
    """
    Wrapper for Casanovo predictions.
    Requires:
      - casanovo binary
      - model checkpoint
      - config yaml
    """

    def __init__(self, sample_path, gpu_id, binary, model_ckpt, config_yaml):
        super().__init__(sample_path, "Casanovo", gpu_id)
        self.binary = binary
        self.model_ckpt = model_ckpt
        self.config_yaml = config_yaml

    @timing
    def run(self):
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = self.gpu_id
        for mgf_file in self.mgf_list:
            base = os.path.basename(mgf_file).split(".")[0]
            out_txt = os.path.join(self.work_dir, f"{base}_result.txt")
            cmd = (
                f"{self.binary} sequence {mgf_file} "
                f"-m {self.model_ckpt} "
                f"-o {out_txt} "
                f"--config {self.config_yaml}"
            )
            run_cmd(cmd, env=env)

    def parse(self) -> pd.DataFrame:
        # Casanovo writes mztab-like TSVs; filter out "MTD" header lines.
        # We'll concatenate all such outputs.
        all_rows = []
        for path in glob.glob(os.path.join(self.work_dir, "*.mztab*")):
            with open(path, "r") as fh:
                lines = [ln for ln in fh if not ln.startswith("MTD")]
            df_single = pd.read_csv(StringIO("".join(lines)), sep="\t")
            all_rows.append(df_single)
        if not all_rows:
            return pd.DataFrame(columns=["Sequence", "ScoreList"])

        df = pd.concat(all_rows, ignore_index=True)
        df = df.rename(columns={
            "sequence": "Sequence",
            "opt_ms_run[1]_aa_scores": "Score"
        })
        df["Sequence"] = df["Sequence"].apply(lambda s: re.sub("[^A-Z]", "", s))
        df["ScoreList"] = df["Score"].apply(
            lambda x: list(map(float, str(x).split(",")))
        )
        return df[["Sequence", "ScoreList"]]


class PepNetEngine(DeNovoEngine):
    """
    Wrapper for PepNet.
    Requires:
      - Python interpreter for PepNet env
      - pepnet denovo script
      - model.h5
    """

    def __init__(self, sample_path, gpu_id, python_bin, script, model_h5):
        super().__init__(sample_path, "PepNet", gpu_id)
        self.python_bin = python_bin
        self.script = script
        self.model_h5 = model_h5

    @timing
    def run(self):
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = self.gpu_id  # adjust GPU use
        for mgf_file in self.mgf_list:
            base = os.path.basename(mgf_file).split(".")[0]
            out_file = os.path.join(self.work_dir, f"{base}.result")
            cmd = (
                f"{self.python_bin} {self.script} "
                f"--input {mgf_file} "
                f"--model {self.model_h5} "
                f"--output {out_file}"
            )
            run_cmd(cmd, env=env)

    def parse(self) -> pd.DataFrame:
        dfs = []
        for res_file in glob.glob(os.path.join(self.work_dir, "*.result")):
            df_single = pd.read_csv(res_file, sep="\t")
            dfs.append(df_single)
        if not dfs:
            return pd.DataFrame(columns=["Sequence", "ScoreList"])

        df = pd.concat(dfs, ignore_index=True)
        df = df.rename(columns={
            "DENOVO": "Sequence",
            "Positional Score": "Score"
        })
        # Clean "Score": "[0.1,0.2,...]" → list[float]
        df["Score"] = df["Score"].apply(
            lambda x: re.sub(r"[\[\]\s]", "", str(x))
        )
        df["ScoreList"] = df["Score"].apply(
            lambda s: [float(v) for v in s.split(",") if v != ""]
        )

        df["Sequence"] = df["Sequence"].apply(lambda s: re.sub("[^A-Z]", "", s))
        df = df.dropna()
        return df[["Sequence", "ScoreList"]]


class InstaNovoEngine(DeNovoEngine):
    """
    Wrapper for InstaNovo.
    Requires:
      - python env with instanovo installed
      - model checkpoint
    """

    def __init__(self, sample_path, gpu_id, python_bin, module_call, model_ckpt, extra_env=None):
        super().__init__(sample_path, "InstaNovo", gpu_id)
        self.python_bin = python_bin
        self.module_call = module_call
        self.model_ckpt = model_ckpt
        self.extra_env = extra_env or {}

    @timing
    def run(self):
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = self.gpu_id
        env.update(self.extra_env)

        # InstaNovo can process *.mgf/*.mzML glob in one go
        mgf_glob = os.path.join(self.mgf_dir, "*.mgf")
        out_file = os.path.join(self.work_dir, "InstaNovo.result")

        cmd = (
            f'{self.python_bin} -m {self.module_call} '
            f'data_path={mgf_glob} '
            f'model_path={self.model_ckpt} '
            f'output_path={out_file} '
            f'denovo=True'
        )
        run_cmd(cmd, env=env)

    def parse(self) -> pd.DataFrame:
        csv_path = os.path.join(self.work_dir, "InstaNovo.result")
        if not os.path.exists(csv_path):
            return pd.DataFrame(columns=["Sequence", "ScoreList"])

        df = pd.read_csv(csv_path)
        df = df.rename(columns={
            "preds": "Sequence",
            "token_log_probs": "RawScores"
        })
        # "RawScores" like "[logp1,logp2,...]". Convert each logp to prob via exp().
        def to_probs(x):
            x = x.strip()[1:-1]  # remove outer brackets
            return [float(np.exp(float(v))) for v in x.split(",")]

        df["ScoreList"] = df["RawScores"].apply(to_probs)
        df["Sequence"] = df["Sequence"].apply(lambda s: re.sub("[^A-Z]", "", s))
        return df[["Sequence", "ScoreList"]]


@timing
def write_denovo_merged_fasta(
    sample_path: str,
    engine_filtered_tables: List[pd.DataFrame],
    out_fasta: str,
):
    """
    Combine all confident subsequences from all engines,
    deduplicate, and write a FASTA where header==sequence.
    This FASTA will later be PE=1-tagged and merged into DB.
    """
    peptides = set()
    for df in engine_filtered_tables:
        if df is None or len(df) == 0:
            continue
        for _, row in df.iterrows():
            for pep in row["Filtered_Subsequences"]:
                peptides.add(pep)

    os.makedirs(os.path.dirname(out_fasta), exist_ok=True)
    with open(out_fasta, "w") as fh:
        for pep in sorted(peptides):
            fh.write(f">{pep}\n{pep}\n")

    print(f"[write_denovo_merged_fasta] wrote {out_fasta} with {len(peptides)} peptides")

