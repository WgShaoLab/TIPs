from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from tips_cli.runner import DockerRunner


@dataclass
class PrepareOutputs:
    sample_path: Path
    mzml_dir: Path
    mgf_dir: Path
    raw_dir: Path
    log_dir: Path
    commands_sh: Path
    log_file: Path


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _find_raw_files(d: Path) -> List[Path]:
    if not d.exists():
        return []
    return sorted([p for p in d.iterdir() if p.is_file() and p.name.lower().endswith(".raw")])


def _stage_raw_files(sample_path: Path, raw_dir: Path, stage_mode: str) -> List[Path]:
    """
    Stage RAW files into sample_path/raw to match TIPs folder convention.

    stage_mode:
      - move: move *.raw from sample_path/ into sample_path/raw/
      - copy: copy *.raw from sample_path/ into sample_path/raw/
      - none: do not move/copy; if raw_dir doesn't exist, use sample_path as raw source
    """
    _ensure_dir(sample_path)
    raw_in_rawdir = _find_raw_files(raw_dir)
    if raw_in_rawdir:
        return raw_in_rawdir

    raw_in_root = _find_raw_files(sample_path)
    if not raw_in_root:
        raise FileNotFoundError(f"No .raw found in {sample_path} or {raw_dir}")

    if stage_mode == "none":
        return raw_in_root

    _ensure_dir(raw_dir)
    staged: List[Path] = []
    for rf in raw_in_root:
        dst = raw_dir / rf.name
        if stage_mode == "move":
            shutil.move(str(rf), str(dst))
        elif stage_mode == "copy":
            shutil.copy2(str(rf), str(dst))
        else:
            raise ValueError("stage must be one of: move, copy, none")
        staged.append(dst)

    return staged


def _expected_out(out_dir: Path, raw_file: Path, ext: str) -> Path:
    return out_dir / (raw_file.stem + ext)


def _write_cmd(commands_sh: Path, cmd_line: str) -> None:
    with commands_sh.open("a", encoding="utf-8") as f:
        f.write(cmd_line + "\n")


def _convert_one(
    runner: DockerRunner,
    raw_file: Path,
    out_dir: Path,
    fmt_flag: int,
    expected: Path,
    dry_run: bool,
    resume: bool,
    commands_sh: Path,
    log_file: Path,
) -> Tuple[Path, int]:
    """
    Convert a single RAW to mzML/mgf.
    """
    if resume and expected.exists():
        return expected, 0

    inner_cmd = [
        "mono",
        "/opt/thermorawfileparser/ThermoRawFileParser.exe",
        "-i",
        str(raw_file),
        "-o",
        str(out_dir),
        f"--format={fmt_flag}",
    ]

    full_cmd = runner.build_run_cmd(inner_cmd)
    _write_cmd(commands_sh, runner.cmd_to_shell(full_cmd))

    rc = runner.run(inner_cmd, dry_run=dry_run, log_file=log_file)
    return expected, rc


def prepare_step(
    outputs: PrepareOutputs,
    image: str,
    do_mzml: bool,
    do_mgf: bool,
    stage_mode: str,
    jobs: int,
    dry_run: bool,
    resume: bool,
) -> None:
    """
    v0.1 prepare:
      - ensure TIPs folder layout: raw/ mzML/ mgf/ logs/
      - stage RAW files into raw/ (default move)
      - run ThermoRawFileParser in docker (tips-core) to generate mzML/mgf
      - write a reproducible commands.sh and a log file
    """
    sample_path = outputs.sample_path.resolve()

    _ensure_dir(outputs.log_dir)
    _ensure_dir(outputs.mzml_dir)
    _ensure_dir(outputs.mgf_dir)

    # Always start a fresh commands.sh for easier debugging
    outputs.commands_sh.write_text("#!/usr/bin/env bash\nset -euo pipefail\n\n", encoding="utf-8")

    raw_files = _stage_raw_files(sample_path, outputs.raw_dir, stage_mode)

    runner = DockerRunner(
        image=image,
        mount_root=sample_path,
        user_map=True,
        workdir_in_container="/tmp",
    )

    tasks = []
    with ThreadPoolExecutor(max_workers=max(1, int(jobs))) as ex:
        for rf in raw_files:
            if do_mzml:
                expected = _expected_out(outputs.mzml_dir, rf, ".mzML")
                tasks.append(
                    ex.submit(
                        _convert_one,
                        runner,
                        rf,
                        outputs.mzml_dir,
                        1,
                        expected,
                        dry_run,
                        resume,
                        outputs.commands_sh,
                        outputs.log_file,
                    )
                )
            if do_mgf:
                expected = _expected_out(outputs.mgf_dir, rf, ".mgf")
                tasks.append(
                    ex.submit(
                        _convert_one,
                        runner,
                        rf,
                        outputs.mgf_dir,
                        0,
                        expected,
                        dry_run,
                        resume,
                        outputs.commands_sh,
                        outputs.log_file,
                    )
                )

        for fut in as_completed(tasks):
            expected, rc = fut.result()
            if rc != 0:
                raise RuntimeError(f"RAW conversion failed: expected={expected} rc={rc}")
