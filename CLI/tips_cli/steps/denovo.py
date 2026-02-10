from __future__ import annotations
from importlib import resources

import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

def default_denovo_script_path() -> Path:
    # This returns a real filesystem path in editable install and normal install.
    return Path(resources.files("tips_cli.tools").joinpath("run_denovo_tool_batch_final.py"))

@dataclass
class DenovoOutputs:
    sample_path: Path
    mgf_dir: Path
    out_dir: Path          # typically: sample_path/denovo
    logs_dir: Path         # typically: sample_path/logs
    commands_sh: Path      # logs/denovo_<tool>_commands.sh
    log_file: Path         # logs/denovo_<tool>.log


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _list_mgfs(mgf_dir: Path) -> List[Path]:
    if not mgf_dir.exists():
        return []
    return sorted([p for p in mgf_dir.iterdir() if p.is_file() and p.suffix.lower() == ".mgf"])


def _instanovo_expected_files(out_dir: Path, sample: Optional[str], mgfs: List[Path]) -> List[Path]:
    inst_dir = out_dir / "instanovo"
    if sample is None:
        return [inst_dir / f"{m.stem}.result.csv" for m in mgfs]
    return [inst_dir / f"{sample}__{m.stem}.result.csv" for m in mgfs]


def _nonempty(p: Path) -> bool:
    return p.exists() and p.is_file() and p.stat().st_size > 0


def _expected_files_for_tool(out_dir: Path, tool: str, sample: Optional[str], mgfs: List[Path]) -> List[Path]:
    """
    v0.2 strategy:
      - instanovo: coarse check per-mgf expected outputs (csv)
      - casanovo: coarse check per-mgf expected outputs (mztab)
      - pepnet: coarse check (result filename may vary; require at least one non-empty output)
    """
    tool_dir = out_dir / tool

    if tool == "instanovo":
        return _instanovo_expected_files(out_dir, sample, mgfs)

    if tool == "casanovo":
        # Casanovo outputs: <sample>__<mgf_basename>_result.mztab
        if sample is None:
            return [tool_dir / f"{m.stem}_result.mztab" for m in mgfs]
        return [tool_dir / f"{sample}__{m.stem}_result.mztab" for m in mgfs]

    if tool == "pepnet":
        # v0.2 coarse: any non-empty result file under pepnet dir
        if not tool_dir.exists():
            return []
        return sorted([p for p in tool_dir.rglob("*") if p.is_file()])

    if tool == "all":
        # Combine expected outputs from all supported tools.
        expected: List[Path] = []
        expected += _expected_files_for_tool(out_dir, "instanovo", sample, mgfs)
        expected += _expected_files_for_tool(out_dir, "casanovo", sample, mgfs)
        expected += _expected_files_for_tool(out_dir, "pepnet", sample, mgfs)
        return expected

    return []


def denovo_step_batch_script(
    outputs: DenovoOutputs,
    script_path: Path,
    tool: str,
    sample: Optional[str],
    gpu: str = "auto",
    param_file: Optional[Path] = None,
    dry_run: bool = False,
    resume: bool = False,
) -> None:
    """
    v0.2: Call the external batch runner script directly.

    The script is responsible for docker execution.
    CLI responsibilities:
      - validate mgf inputs
      - decide outdir/log paths
      - provide reproducible command script
      - implement resume based on expected final outputs
    """
    script_path = script_path.resolve()
    if not script_path.is_file():
        raise FileNotFoundError(f"denovo batch script not found: {script_path}")

    _ensure_dir(outputs.out_dir)
    _ensure_dir(outputs.logs_dir)
    _ensure_dir(outputs.commands_sh.parent)

    mgfs = _list_mgfs(outputs.mgf_dir)
    if not mgfs:
        raise FileNotFoundError(f"No .mgf files found under: {outputs.mgf_dir}")

    allowed = {"instanovo", "casanovo", "pepnet", "all"}
    if tool not in allowed:
        raise ValueError(f"Unsupported tool: {tool}. Choose from: {sorted(allowed)}")

    if resume:
        expected = _expected_files_for_tool(outputs.out_dir, tool, sample, mgfs)

        if tool == "instanovo":
            # Strict per-mgf check
            if expected and all(_nonempty(p) for p in expected):
                return
        else:
            # Coarse check for casanovo/pepnet: if there is any non-empty result, assume done
            if expected and any(_nonempty(p) for p in expected):
                return

    # Use glob pattern for the batch script (it supports glob and avoids repeated model loading)
    mgf_glob = str((outputs.mgf_dir / "*.mgf").resolve())

    cmd: List[str] = [
        "python",
        str(script_path),
        "--tool",
        tool,
        "--input",
        mgf_glob,
        "--outdir",
        str(outputs.out_dir.resolve()),
        "--gpu",
        gpu,
    ]

    if sample is not None and sample.strip():
        cmd += ["--sample", sample]

    if param_file is not None:
        cmd += ["--param-file", str(param_file.resolve())]

    # Write reproducible commands.sh
    with outputs.commands_sh.open("w", encoding="utf-8") as f:
        f.write("#!/usr/bin/env bash\nset -euo pipefail\n\n")
        f.write(" ".join(shlex.quote(x) for x in cmd) + "\n")

    if dry_run:
        return

    # Stream stdout/stderr to log file
    with outputs.log_file.open("w", encoding="utf-8") as logf:
        p = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"denovo batch script failed (rc={p.returncode}). See log: {outputs.log_file}")

    # Post-check expected outputs for InstaNovo
    expected = _expected_files_for_tool(outputs.out_dir, tool, sample, mgfs)

    if tool == "instanovo":
        missing = [p for p in expected if not _nonempty(p)]
        if missing:
            raise RuntimeError("InstaNovo outputs missing/empty:\n" + "\n".join(str(x) for x in missing))
    else:
        # Coarse post-check: require at least one result for casanovo/pepnet
        if not expected or not any(_nonempty(p) for p in expected):
            raise RuntimeError(f"{tool} outputs missing/empty under: {outputs.out_dir / tool}")

