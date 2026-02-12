#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import glob
import os
import re
import shlex
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import List, Optional, Tuple


# ----------------------------
# Defaults (images and runtime knobs)
# ----------------------------
DEFAULT_IMG_CASANOVO = os.environ.get("IMG_CASANOVO", "tips-denovo-casanovo:0.2")
DEFAULT_IMG_PEPNET = os.environ.get("IMG_PEPNET", "tips-denovo-pepnet:0.2")
DEFAULT_IMG_INSTANOVO = os.environ.get("IMG_INSTANOVO", "tips-denovo-instanovo:0.2")

# GPU_FLAG is used only when we decide to expose GPU to container.
# Keep it as "--gpus all" for auto GPU selection to work correctly.
GPU_FLAG = os.environ.get("GPU_FLAG", "--gpus all")

# Per-tool GPU selection (optional): "auto" | "-1" | "0" | "1" ...
# If set, these override --gpu for the corresponding tool.
ENV_CASA_GPU = os.environ.get("CASA_GPU", "").strip()
ENV_PEPNET_GPU = os.environ.get("PEPNET_GPU", "").strip()
ENV_INSTA_GPU = os.environ.get("INSTA_GPU", "").strip()

# Optional persistent CUDA cache on host (helpful for TF/PTX JIT)
CUDA_CACHE_HOST = os.environ.get("CUDA_CACHE_HOST", "").strip()
CUDA_CACHE_CONT = "/tmp/nv_cache"

# Only support RTX 40 series and below => compute capability <= 8.9 (SM89).
# Use MAX_CC env to override, e.g. MAX_CC=86 for Ampere-only.
MAX_COMPUTE_CAP_INT = int(os.environ.get("MAX_CC", "89"))  # "8.9" -> 89

DOCKER_COMMON_OPTS = [
    "--rm",
    "-e",
    "HOME=/tmp",
    "-e",
    "TF_FORCE_GPU_ALLOW_GROWTH=true",
    "-e",
    "PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:512",
    "--shm-size=2g",
]


# ----------------------------
# Helper utilities
# ----------------------------

def eprint(*args) -> None:
    print(*args, file=sys.stderr)


def run_cmd(cmd: List[str], check: bool = True, capture: bool = False) -> Tuple[int, str]:
    """Run a command. Return (return_code, combined_stdout_stderr_if_capture_else_empty)."""
    try:
        if capture:
            p = subprocess.run(
                cmd,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            out = p.stdout or ""
        else:
            p = subprocess.run(cmd, check=False)
            out = ""
        if check and p.returncode != 0:
            raise subprocess.CalledProcessError(p.returncode, cmd, output=out)
        return p.returncode, out
    except FileNotFoundError:
        return 127, ""


def docker_available() -> bool:
    rc, _ = run_cmd(["docker", "version"], check=False, capture=True)
    return rc == 0


def abs_path(p: str) -> str:
    return str(Path(p).expanduser().resolve())


def expand_mgf_inputs(items: List[str]) -> List[Path]:
    """Expand user inputs into a sorted list of existing .mgf files.

    Supported inputs:
      - A file path (must exist)
      - A directory path (all *.mgf files in that directory)
      - A glob pattern (e.g., /data/mgf/*.mgf)

    Notes:
      - This function does not recurse into subdirectories.
      - Returned paths are absolute and de-duplicated.
    """
    mgfs: List[Path] = []
    for item in items:
        has_glob = any(ch in item for ch in ["*", "?", "["])
        p = Path(item).expanduser()

        if has_glob:
            candidates = sorted(glob.glob(str(p)))
        elif p.is_dir():
            candidates = sorted(glob.glob(str(p / "*.mgf")))
        else:
            candidates = [str(p)]

        for c in candidates:
            cp = Path(c).expanduser().resolve()
            if cp.is_file() and cp.suffix.lower() == ".mgf":
                mgfs.append(cp)

    # De-duplicate while preserving order
    seen = set()
    out: List[Path] = []
    for p in mgfs:
        if p not in seen:
            seen.add(p)
            out.append(p)

    return sorted(out)


def ensure_dirs(out_abs: str) -> None:
    Path(out_abs, "casanovo").mkdir(parents=True, exist_ok=True)
    Path(out_abs, "pepnet").mkdir(parents=True, exist_ok=True)
    Path(out_abs, "instanovo").mkdir(parents=True, exist_ok=True)
    Path(out_abs, "logs").mkdir(parents=True, exist_ok=True)
    Path(out_abs, ".tips_tmp").mkdir(parents=True, exist_ok=True)


def list_gpu_indices() -> List[str]:
    """Return a list of GPU indices from nvidia-smi. If nvidia-smi is unavailable, return []."""
    rc, out = run_cmd(
        ["bash", "-lc", "command -v nvidia-smi >/dev/null 2>&1 && echo OK || echo NO"],
        check=False,
        capture=True,
    )
    if "OK" not in out:
        return []

    rc, out = run_cmd(
        ["nvidia-smi", "--query-gpu=index", "--format=csv,noheader"],
        check=False,
        capture=True,
    )
    if rc == 0 and out.strip():
        return [line.strip() for line in out.strip().splitlines() if line.strip()]

    # Fallback: parse "nvidia-smi -L"
    rc, out = run_cmd(["nvidia-smi", "-L"], check=False, capture=True)
    idxs: List[str] = []
    if rc == 0 and out.strip():
        for line in out.splitlines():
            m = re.match(r"GPU\s+(\d+):", line.strip())
            if m:
                idxs.append(m.group(1))
    return idxs


def docker_gpu_runtime_usable(img_for_check: str) -> bool:
    """Check whether docker GPU runtime is usable at all."""
    if not GPU_FLAG.strip():
        return False

    script = 'python -c "import torch; print(\\"1\\" if torch.cuda.is_available() else \\"0\\")"'
    cmd = ["docker", "run", "--rm"] + shlex.split(GPU_FLAG) + [img_for_check, "/bin/bash", "-lc", script]
    rc, out = run_cmd(cmd, check=False, capture=True)
    return (rc == 0) and ("1" in out)


def container_gpu_ok_for_index(img_for_check: str, idx: str) -> bool:
    """Verify GPU usability INSIDE container and enforce compute capability <= MAX_COMPUTE_CAP_INT."""
    script = textwrap.dedent(
        f"""\
        python - <<'PY'
        import torch

        ok_alloc = False
        cap_int = 999

        if torch.cuda.is_available():
            try:
                cap_tuple = torch.cuda.get_device_capability(0)
                cap_int = cap_tuple[0] * 10 + cap_tuple[1]
                x = torch.ones(1, device="cuda")
                torch.cuda.synchronize()
                ok_alloc = True
            except Exception:
                ok_alloc = False

        print(f"CAP={{cap_int}}")
        if ok_alloc and cap_int <= {MAX_COMPUTE_CAP_INT}:
            print("GPU_OK")
        else:
            print("GPU_BAD")
        PY
        """
    ).strip()

    cmd = (
        ["docker", "run", "--rm"]
        + shlex.split(GPU_FLAG)
        + ["-e", f"CUDA_VISIBLE_DEVICES={idx}", img_for_check, "/bin/bash", "-lc", script]
    )
    rc, out = run_cmd(cmd, check=False, capture=True)
    return (rc == 0) and ("GPU_OK" in out)


def pick_compatible_gpu_index(img_for_check: str) -> Optional[str]:
    """Auto-pick the first GPU that passes container capability check."""
    idxs = list_gpu_indices()
    if not idxs:
        return None

    for idx in idxs:
        if container_gpu_ok_for_index(img_for_check, idx):
            return idx
    return None


def build_mounts(mgf_abs: str, out_abs: str, param_file: Optional[str]) -> List[str]:
    """Build docker -v mounts."""
    in_dir = str(Path(mgf_abs).parent)
    mounts = [
        "-v",
        f"{in_dir}:/in:ro",
        "-v",
        f"{out_abs}:/out:rw",
    ]

    if CUDA_CACHE_HOST:
        Path(CUDA_CACHE_HOST).mkdir(parents=True, exist_ok=True)
        mounts += ["-v", f"{CUDA_CACHE_HOST}:{CUDA_CACHE_CONT}:rw"]

    if param_file:
        p = Path(param_file)
        mounts += ["-v", f"{str(p.parent.resolve())}:/user_params:ro"]

    return mounts


def tool_gpu_choice(tool: str, cli_gpu: str) -> str:
    """Determine GPU selection string for a tool."""
    if tool == "casanovo" and ENV_CASA_GPU:
        return ENV_CASA_GPU
    if tool == "pepnet" and ENV_PEPNET_GPU:
        return ENV_PEPNET_GPU
    if tool == "instanovo" and ENV_INSTA_GPU:
        return ENV_INSTA_GPU
    return cli_gpu


def resolve_gpu_index_for_tool(
    gpu_runtime_ok: bool,
    img_for_check: str,
    selection: str,
) -> Optional[str]:
    """Resolve GPU index based on selection."""
    s = selection.strip().lower()
    if not gpu_runtime_ok:
        return None
    if s == "-1":
        return None
    if s == "" or s == "auto":
        return pick_compatible_gpu_index(img_for_check)
    if not s.isdigit():
        return None
    idx = s
    if container_gpu_ok_for_index(img_for_check, idx):
        return idx
    return None


def docker_run(
    img: str,
    script: str,
    mgf_abs: str,
    out_abs: str,
    sample: str,
    gpu_index: Optional[str],
    param_file: Optional[str],
) -> None:
    """Run docker container with proper mounts and optional GPU exposure."""
    envs = ["-e", f"SAMPLE={sample}"]

    if CUDA_CACHE_HOST:
        envs += [
            "-e",
            f"CUDA_CACHE_PATH={CUDA_CACHE_CONT}",
            "-e",
            "CUDA_CACHE_MAXSIZE=2147483648",
        ]

    if param_file:
        envs += ["-e", f"USER_PARAM_NAME={Path(param_file).name}"]

    gpu_opts: List[str] = []
    if gpu_index is not None:
        gpu_opts += shlex.split(GPU_FLAG)
        envs += ["-e", f"CUDA_VISIBLE_DEVICES={gpu_index}"]
        print(f"[RUN] image={img}  GPU=ON  CUDA_VISIBLE_DEVICES={gpu_index}")
    else:
        envs += ["-e", "CUDA_VISIBLE_DEVICES="]
        print(f"[RUN] image={img}  GPU=OFF -> CPU")

    user_flag = ["--user", f"{os.getuid()}:{os.getgid()}"]
    mounts = build_mounts(mgf_abs, out_abs, param_file)

    cmd = (
        ["docker", "run"]
        + DOCKER_COMMON_OPTS
        + gpu_opts
        + user_flag
        + envs
        + ["-w", "/out"]
        + mounts
        + [img, "/bin/bash", "-lc", script]
    )

    run_cmd(cmd, check=True, capture=False)


# ----------------------------
# Container scripts (bash)
# ----------------------------

SCRIPT_CASANOVO = r"""
set -euo pipefail
: "${SAMPLE:?missing SAMPLE}"

LIST="/out/.tips_tmp/mgf_list.txt"
[[ -f "$LIST" ]] || { echo "[Casanovo][ERROR] mgf list not found: $LIST" >&2; exit 2; }

mkdir -p /out/casanovo /out/logs /tmp/casa_in
rm -f /tmp/casa_in/*.mgf 2>/dev/null || true
while IFS= read -r f; do
  [[ -n "$f" ]] || continue
  ln -sf "/in/$f" "/tmp/casa_in/$f"
done < "$LIST"

shopt -s nullglob
IN_FILES=(/tmp/casa_in/*.mgf)
[[ ${#IN_FILES[@]} -gt 0 ]] || { echo "[Casanovo][ERROR] no .mgf files found under /tmp/casa_in" >&2; exit 2; }

OUT_MERGED="/out/.tips_tmp/${SAMPLE}__casanovo_merged_result.mztab"
LOG_MERGED="/out/.tips_tmp/${SAMPLE}.casanovo.batch.log"

MODEL_CAND=(/opt/models/casanovo/*.{ckpt,pt,pth} /opt/models/casanovo/**/*.{ckpt,pt,pth})
CFG_CAND=(/opt/models/casanovo/*.{yaml,yml} /opt/models/casanovo/**/*.{yaml,yml})
MODEL="${MODEL_CAND[0]:-}"
CFG="${CFG_CAND[0]:-}"

if [[ -n "${USER_PARAM_NAME:-}" ]] && [[ -f "/user_params/${USER_PARAM_NAME}" ]]; then
  CFG="/user_params/${USER_PARAM_NAME}"
fi

echo "[Casanovo] inputs=${#IN_FILES[@]} files"
echo "[Casanovo] model=$MODEL"
echo "[Casanovo] config=${CFG:-<none>}"
[[ -n "$MODEL" ]] || { echo "[Casanovo][ERROR] cannot find model under /opt/models/casanovo" >&2; exit 2; }

USE_GPU=0
set +e
python - <<'PY' >/tmp/casa_gpu_check.txt 2>&1
import torch
ok = False
if torch.cuda.is_available():
    try:
        x = torch.ones(1, device="cuda")
        torch.cuda.synchronize()
        ok = True
    except Exception:
        ok = False
print("GPU_OK" if ok else "GPU_BAD")
PY
rc=$?
set -e
if [[ $rc -eq 0 ]] && grep -q "GPU_OK" /tmp/casa_gpu_check.txt; then
  USE_GPU=1
else
  USE_GPU=0
fi

if [[ "$USE_GPU" -eq 1 ]]; then
  echo "[Casanovo] GPU usable: yes"
else
  echo "[Casanovo] GPU usable: no -> run on CPU"
  export CUDA_VISIBLE_DEVICES=""

  if [[ -n "${CFG:-}" ]] && [[ -f "$CFG" ]]; then
    export CFG_SRC="$CFG"
    set +e
    python - <<'PY' >/tmp/casa_cfg_patch.txt 2>&1
import os, sys
from pathlib import Path

src = os.environ.get("CFG_SRC","")
dst = "/tmp/casanovo.cpu.yaml"

try:
    import yaml
    cfg = yaml.safe_load(open(src))
    def patch(obj):
        if isinstance(obj, dict):
            if "accelerator" in obj: obj["accelerator"] = "cpu"
            if "devices" in obj: obj["devices"] = 1
            if "gpus" in obj: obj["gpus"] = 0
            if "precision" in obj and str(obj["precision"]).startswith("16"):
                obj["precision"] = 32
            for v in obj.values():
                patch(v)
        elif isinstance(obj, list):
            for v in obj:
                patch(v)
    patch(cfg)
    yaml.safe_dump(cfg, open(dst, "w"), sort_keys=False)
    print("WROTE", dst)
    sys.exit(0)
except Exception:
    try:
        txt = Path(src).read_text()
        import re
        if re.search(r'^\s*accelerator\s*:', txt, flags=re.M):
            txt = re.sub(r'^\s*accelerator\s*:\s*.*$', 'accelerator: "cpu"', txt, flags=re.M)
        else:
            txt = 'accelerator: "cpu"\n' + txt

        if re.search(r'^\s*devices\s*:', txt, flags=re.M):
            txt = re.sub(r'^\s*devices\s*:\s*.*$', 'devices: 1', txt, flags=re.M)
        else:
            txt = 'devices: 1\n' + txt

        Path(dst).write_text(txt)
        print("WROTE", dst)
        sys.exit(0)
    except Exception as e2:
        print("PATCH_FAILED", e2)
        sys.exit(2)
PY
    prc=$?
    set -e
    if [[ $prc -eq 0 ]] && grep -q "WROTE" /tmp/casa_cfg_patch.txt; then
      CFG="/tmp/casanovo.cpu.yaml"
      echo "[Casanovo] patched cpu config: $CFG"
    else
      echo "[Casanovo][WARN] could not patch config to CPU. Will try with original config."
    fi
  fi
fi

if [[ -n "${CFG:-}" ]]; then
  casanovo sequence /tmp/casa_in/*.mgf -m "$MODEL" -o "$OUT_MERGED" --config "$CFG" 2>&1 | tee "$LOG_MERGED"
else
  casanovo sequence /tmp/casa_in/*.mgf -m "$MODEL" -o "$OUT_MERGED" 2>&1 | tee "$LOG_MERGED"
fi

[[ -s "$OUT_MERGED" ]] || { echo "[Casanovo][ERROR] merged output not created or empty: $OUT_MERGED" >&2; exit 3; }

export OUT_MERGED="$OUT_MERGED"
export LOG_MERGED="$LOG_MERGED"
python - <<'PY'
import os, re
from pathlib import Path

merged = Path(os.environ.get("OUT_MERGED", "/out/.tips_tmp/merged.mztab"))
sample = os.environ.get("SAMPLE", "sample")
out_dir = Path("/out/casanovo")
log_dir = Path("/out/logs")
log_merged = Path(os.environ.get("LOG_MERGED", "/out/.tips_tmp/merged.log"))

text = merged.read_text().splitlines(True)

spectra_ref_col = None
for line in text:
    if line.startswith("PSH\t"):
        cols = line.rstrip("\n").split("\t")[1:]
        if "spectra_ref" in cols:
            spectra_ref_col = cols.index("spectra_ref")
        break

msrun_loc_re = re.compile(r"^MTD\tms_run\[(\d+)\]-(\w+)\t(.+?)\s*$")
spectra_ref_re = re.compile(r"ms_run\[(\d+)\]:")

run_to_base = {}
for line in text:
    m = msrun_loc_re.match(line)
    if not m:
        continue
    idx = int(m.group(1))
    key = m.group(2)
    val = m.group(3)
    if key != "location":
        continue
    val2 = val[7:] if val.startswith("file://") else val
    val2 = val2.replace("\\", "/")
    fname = val2.rsplit("/", 1)[-1]
    run_to_base[idx] = Path(fname).stem

if not run_to_base:
    list_path = Path("/out/.tips_tmp/mgf_list.txt")
    if list_path.exists():
        first = [x.strip() for x in list_path.read_text().splitlines() if x.strip()][0]
        run_to_base = {1: Path(first).stem}
    else:
        run_to_base = {1: "input"}

def rewrite_ms_run(s: str, src_idx: int) -> str:
    return re.sub(rf"ms_run\[{src_idx}\]", "ms_run[1]", s)

for run_idx in sorted(run_to_base):
    out_lines = []
    for line in text:
        if line.startswith("COM") or line.strip() == "":
            out_lines.append(line)
            continue

        if line.startswith("MTD\tms_run["):
            m = msrun_loc_re.match(line)
            if m:
                idx = int(m.group(1))
                if idx != run_idx:
                    continue
                out_lines.append(rewrite_ms_run(line, run_idx))
                continue
            out_lines.append(line)
            continue

        if line.startswith("PSM\t"):
            if spectra_ref_col is not None:
                fields = line.rstrip("\n").split("\t")
                if 1 + spectra_ref_col >= len(fields):
                    continue
                sref = fields[1 + spectra_ref_col]
                m2 = spectra_ref_re.search(sref)
                if not m2:
                    continue
                idx = int(m2.group(1))
                if idx != run_idx:
                    continue
                fields[1 + spectra_ref_col] = rewrite_ms_run(sref, run_idx)
                out_lines.append("\t".join(fields) + "\n")
            else:
                m2 = spectra_ref_re.search(line)
                if not m2:
                    continue
                idx = int(m2.group(1))
                if idx != run_idx:
                    continue
                out_lines.append(rewrite_ms_run(line, run_idx))
            continue

        out_lines.append(line)

    base = run_to_base[run_idx]
    out_path = out_dir / f"{sample}__{base}_result.mztab"
    out_path.write_text("".join(out_lines))

    if log_merged.exists():
        (log_dir / f"{sample}__{base}.casanovo.log").write_text(log_merged.read_text())

for run_idx, base in run_to_base.items():
    out_path = out_dir / f"{sample}__{base}_result.mztab"
    if not out_path.exists() or out_path.stat().st_size == 0:
        raise SystemExit(f"[Casanovo][ERROR] split output missing/empty: {out_path}")
print(f"[Casanovo][OK] wrote {len(run_to_base)} files into {out_dir}")
PY
""".strip()


SCRIPT_PEPNET = r"""
set -euo pipefail
: "${SAMPLE:?missing SAMPLE}"

LIST="/out/.tips_tmp/mgf_list.txt"
[[ -f "$LIST" ]] || { echo "[PepNet][ERROR] mgf list not found: $LIST" >&2; exit 2; }

mkdir -p /out/pepnet /out/logs /tmp/pep_in
rm -f /tmp/pep_in/*.mgf 2>/dev/null || true
while IFS= read -r f; do
  [[ -n "$f" ]] || continue
  ln -sf "/in/$f" "/tmp/pep_in/$f"
done < "$LIST"

shopt -s nullglob
IN_FILES=(/tmp/pep_in/*.mgf)
[[ ${#IN_FILES[@]} -gt 0 ]] || { echo "[PepNet][ERROR] no .mgf files found under /tmp/pep_in" >&2; exit 2; }

MODEL_CAND=(/opt/models/pepnet/*.{h5,keras} /opt/models/pepnet/**/*.{h5,keras})
MODEL="${MODEL_CAND[0]:-}"
[[ -n "$MODEL" ]] || { echo "[PepNet][ERROR] cannot find model under /opt/models/pepnet" >&2; exit 2; }

OUT_MERGED="/out/.tips_tmp/${SAMPLE}__pepnet_merged.result"
LOG_MERGED="/out/.tips_tmp/${SAMPLE}.pepnet.batch.log"
MERGED_MGF="/tmp/pepnet_merged.mgf"

echo "[PepNet] inputs=${#IN_FILES[@]} files"
echo "[PepNet] model=$MODEL"

USE_GPU=0
set +e
python - <<'PY' >/tmp/pep_gpu_check.txt 2>&1
import tensorflow as tf
ok = False
gpus = tf.config.list_physical_devices("GPU")
if gpus:
    try:
        with tf.device("/GPU:0"):
            x = tf.constant([1.0, 2.0, 3.0])
            y = x * 2.0
        _ = y.numpy()
        ok = True
    except Exception:
        ok = False
print("GPU_OK" if ok else "GPU_BAD")
print("tf:", tf.__version__)
print("gpus:", gpus)
PY
rc=$?
set -e
if [[ $rc -eq 0 ]] && grep -q "GPU_OK" /tmp/pep_gpu_check.txt; then
  USE_GPU=1
else
  USE_GPU=0
fi

cat /tmp/pep_gpu_check.txt || true
if [[ "$USE_GPU" -eq 1 ]]; then
  echo "[PepNet] GPU usable: yes"
else
  echo "[PepNet] GPU usable: no -> run on CPU"
  export CUDA_VISIBLE_DEVICES=""
fi

python - <<'PY'
from pathlib import Path
import re

in_dir = Path("/tmp/pep_in")
out_mgf = Path("/tmp/pepnet_merged.mgf")

begin_re = re.compile(r"^\s*BEGIN\s+IONS\s*$", re.IGNORECASE)
end_re = re.compile(r"^\s*END\s+IONS\s*$", re.IGNORECASE)
title_re = re.compile(r"^\s*TITLE\s*=\s*(.*)\s*$", re.IGNORECASE)

with out_mgf.open("w") as out:
    total = 0
    for f in sorted(in_dir.glob("*.mgf")):
        base = f.stem
        spec_i = 0
        in_block = False
        buf = []
        has_title = False

        with f.open("r", errors="ignore") as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                if begin_re.match(line):
                    in_block = True
                    buf = [line]
                    has_title = False
                    continue

                if not in_block:
                    continue

                m = title_re.match(line)
                if m:
                    has_title = True
                    orig_title = m.group(1).strip()
                    buf.append(f"TITLE={base}||{orig_title}")
                    continue

                if end_re.match(line):
                    spec_i += 1
                    if not has_title:
                        buf.insert(1, f"TITLE={base}||spec_{spec_i}")
                    buf.append(line)
                    out.write("\n".join(buf) + "\n")
                    total += 1
                    in_block = False
                    buf = []
                    continue

                buf.append(line)

print(f"WROTE {out_mgf} spectra={total}")
PY

python /opt/pepnet/denovo.py --input "$MERGED_MGF" --model "$MODEL" --output "$OUT_MERGED" 2>&1 | tee "$LOG_MERGED"
[[ -s "$OUT_MERGED" ]] || { echo "[PepNet][ERROR] merged output not created or empty: $OUT_MERGED" >&2; exit 3; }

python - <<'PY'
import os
from pathlib import Path

sample = os.environ.get("SAMPLE", "sample")
merged = Path("/out/.tips_tmp") / f"{sample}__pepnet_merged.result"
list_path = Path("/out/.tips_tmp/mgf_list.txt")
out_dir = Path("/out/pepnet")

bases = [Path(x.strip()).stem for x in list_path.read_text().splitlines() if x.strip()]

handles = {}
try:
    with merged.open("r", errors="ignore") as fin:
        header = fin.readline()
        if not header:
            raise SystemExit("[PepNet][ERROR] merged result is empty")

        for base in bases:
            p = out_dir / f"{sample}__{base}.result"
            h = p.open("w")
            h.write(header)
            handles[base] = h

        for line in fin:
            if not line.strip():
                continue
            first, rest = (line.split("\t", 1) + [""])[:2]
            if "||" not in first:
                continue
            base, orig = first.split("||", 1)
            h = handles.get(base)
            if h is None:
                continue
            if rest:
                h.write(orig + "\t" + rest)
            else:
                h.write(orig + "\n")
finally:
    for h in handles.values():
        try:
            h.close()
        except Exception:
            pass

for base in bases:
    p = out_dir / f"{sample}__{base}.result"
    if not p.exists() or p.stat().st_size == 0:
        raise SystemExit(f"[PepNet][ERROR] split output missing/empty: {p}")
print(f"[PepNet][OK] wrote {len(bases)} files into {out_dir}")
PY

while IFS= read -r f; do
  [[ -n "$f" ]] || continue
  base="$(basename "$f" .mgf)"
  cp -f "$LOG_MERGED" "/out/logs/${SAMPLE}__${base}.pepnet.log"
done < "$LIST"
""".strip()


SCRIPT_INSTANOVO = r"""
set -euo pipefail
: "${SAMPLE:?missing SAMPLE}"

LIST="/out/.tips_tmp/mgf_list.txt"
[[ -f "$LIST" ]] || { echo "[InstaNovo][ERROR] mgf list not found: $LIST" >&2; exit 2; }

mkdir -p /out/instanovo /out/logs /tmp/insta_in
rm -f /tmp/insta_in/*.mgf 2>/dev/null || true
while IFS= read -r f; do
  [[ -n "$f" ]] || continue
  ln -sf "/in/$f" "/tmp/insta_in/$f"
done < "$LIST"

shopt -s nullglob
IN_FILES=(/tmp/insta_in/*.mgf)
[[ ${#IN_FILES[@]} -gt 0 ]] || { echo "[InstaNovo][ERROR] no .mgf files found under /tmp/insta_in" >&2; exit 2; }

MODEL="/opt/models/instanovo/instanovo_extended.ckpt"
[[ -f "$MODEL" ]] || { echo "[InstaNovo][ERROR] required model not found: $MODEL" >&2; exit 2; }

OUT_MERGED="/out/.tips_tmp/${SAMPLE}__instanovo_merged.result.csv"
LOG_MERGED="/out/.tips_tmp/${SAMPLE}.instanovo.batch.log"

echo "[InstaNovo] inputs=${#IN_FILES[@]} files"
echo "[InstaNovo] model=$MODEL"

DEV="cpu"
FP16="false"

set +e
python - <<'PY' >/tmp/insta_gpu_check.txt 2>&1
import torch
ok = False
if torch.cuda.is_available():
    try:
        x = torch.ones(1, device="cuda")
        torch.cuda.synchronize()
        ok = True
    except Exception:
        ok = False
print("GPU_OK" if ok else "GPU_BAD")
print("torch:", torch.__version__)
print("cuda_available:", torch.cuda.is_available())
print("device_count:", torch.cuda.device_count())
PY
rc=$?
set -e

if [[ $rc -eq 0 ]] && grep -q "GPU_OK" /tmp/insta_gpu_check.txt; then
  DEV="cuda"
  FP16="true"
else
  DEV="cpu"
  FP16="false"
  export CUDA_VISIBLE_DEVICES=""
fi

cat /tmp/insta_gpu_check.txt || true
echo "[InstaNovo] device=$DEV  fp16=$FP16"

python -m instanovo.transformer.predict \
  data_path="/tmp/insta_in/*.mgf" \
  data_type=mgf \
  model_path="$MODEL" \
  output_path="$OUT_MERGED" \
  denovo=True \
  device="$DEV" \
  fp16="$FP16" \
  2>&1 | tee "$LOG_MERGED"

[[ -s "$OUT_MERGED" ]] || { echo "[InstaNovo][ERROR] merged output not created or empty: $OUT_MERGED" >&2; exit 3; }

python - <<'PY'
import csv, os
from pathlib import Path

sample = os.environ.get("SAMPLE", "sample")
merged = merged = Path("/out/.tips_tmp") / f"{sample}__instanovo_merged.result.csv"
list_path = Path("/out/.tips_tmp/mgf_list.txt")
out_dir = Path("/out/instanovo")

bases = [Path(x.strip()).stem for x in list_path.read_text().splitlines() if x.strip()]

with merged.open(newline="") as fin:
    reader = csv.DictReader(fin)
    header = reader.fieldnames
    if not header:
        raise SystemExit("[InstaNovo][ERROR] CSV header missing")

    handles = {}
    writers = {}
    try:
        for exp in bases:
            p = out_dir / f"{sample}__{exp}.result.csv"
            h = p.open("w", newline="")
            w = csv.DictWriter(h, fieldnames=header)
            w.writeheader()
            handles[exp] = h
            writers[exp] = w

        for row in reader:
            exp = row.get("experiment_name", "")
            w = writers.get(exp)
            if w is not None:
                w.writerow(row)
    finally:
        for h in handles.values():
            try:
                h.close()
            except Exception:
                pass

for exp in bases:
    p = out_dir / f"{sample}__{exp}.result.csv"
    if not p.exists() or p.stat().st_size == 0:
        raise SystemExit(f"[InstaNovo][ERROR] split output missing/empty: {p}")
print(f"[InstaNovo][OK] wrote {len(bases)} files into {out_dir}")
PY

while IFS= read -r f; do
  [[ -n "$f" ]] || continue
  base="$(basename "$f" .mgf)"
  cp -f "$LOG_MERGED" "/out/logs/${SAMPLE}__${base}.instanovo.log"
done < "$LIST"
""".strip()


# ----------------------------
# Main
# ----------------------------

def main() -> None:
    class RawTextDefaultsHelpFormatter(
        argparse.ArgumentDefaultsHelpFormatter,
        argparse.RawTextHelpFormatter,
    ):
        pass

    description = (
        "Run de novo docker tool (pepnet/casanovo/instanovo) with auto GPU selection.\n"
        "Batch mode: multiple MGF files in the same folder, single model load per tool.\n"
        "GPU auto-selection enforces compute capability <= 8.9 (RTX40 and below) by default.\n"
        "Use MAX_CC to override (e.g. MAX_CC=86).\n"
    )

    epilog = textwrap.dedent(
        """\
        GPU selection logic (priority high -> low):
          1) Per-tool ENV overrides: CASA_GPU / PEPNET_GPU / INSTA_GPU
          2) CLI flag: --gpu
          3) Default: auto

        GPU selection values:
          - auto : pick the first compatible GPU (cap checked inside container)
          - -1   : force CPU
          - 0/1..: specify GPU index (validated inside container)

        Environment variables:
          Images:
            IMG_CASANOVO / IMG_PEPNET / IMG_INSTANOVO
          GPU runtime:
            GPU_FLAG='--gpus all'
            MAX_CC=89
          Optional CUDA cache:
            CUDA_CACHE_HOST=/path/to/nv_cache

        Output layout:
          <outdir>/
            casanovo/   <sample>__<mgf_basename>_result.mztab
            pepnet/     <sample>__<mgf_basename>.result
            instanovo/  <sample>__<mgf_basename>.result.csv
            logs/       <sample>__<mgf_basename>.*.log
          <outdir>/.tips_tmp/   (temporary batch artifacts; safe to delete)

        Notes:
          - Casanovo supports --param-file (mounted and passed to casanovo as --config).
            If CPU fallback happens, the config may be patched to ensure CPU execution.
          - InstaNovo always uses: /opt/models/instanovo/instanovo_extended.ckpt (cannot be overridden here).
          - All input .mgf files must be in the same directory for a single run.

        Examples:
          ./run_denovo_tool.py --tool all --input /data/mgf/*.mgf --outdir out --sample S1
          ./run_denovo_tool.py --tool casanovo --input a.mgf b.mgf --outdir out --sample S1
          ./run_denovo_tool.py --tool all --input /data/mgf/*.mgf --outdir out --sample S1 --gpu -1
        """
    )

    parser = argparse.ArgumentParser(
        formatter_class=RawTextDefaultsHelpFormatter,
        description=description,
        epilog=epilog,
    )

    req = parser.add_argument_group("Required arguments")
    opt = parser.add_argument_group("Optional arguments")

    req.add_argument(
        "--tool",
        required=True,
        choices=["pepnet", "casanovo", "instanovo", "all"],
        help="Which tool to run. Use 'all' to run casanovo+pepnet+instanovo sequentially.",
    )
    req.add_argument(
        "--input",
        "-i",
        nargs="+",
        required=True,
        help=(
            "One or more .mgf files, directories, or glob patterns (same folder). "
            "Examples: --input a.mgf  |  --input /data/mgf/*.mgf"
        ),
    )
    req.add_argument(
        "--outdir",
        required=True,
        help="Output directory on host (will be created if missing).",
    )
    req.add_argument(
        "--sample",
        required=True,
        help="Sample ID prefix used in output filenames (e.g. S1).",
    )

    opt.add_argument(
        "--param-file",
        default="",
        help="Optional external config file for Casanovo (passed as casanovo --config).",
    )
    opt.add_argument(
        "--gpu",
        default="auto",
        help="GPU selection for tools without per-tool ENV override: auto | -1 | 0,1,2,...",
    )

    args = parser.parse_args()

    if not docker_available():
        eprint("[ERROR] docker not available. Please ensure docker is installed and running.")
        sys.exit(1)

    mgf_paths = expand_mgf_inputs(args.input)
    if not mgf_paths:
        eprint("[ERROR] no .mgf inputs found after expanding --input")
        sys.exit(1)

    mgf_dir = mgf_paths[0].parent
    if any(p.parent != mgf_dir for p in mgf_paths):
        eprint("[ERROR] all input .mgf files must be in the same directory for a single run.")
        eprint("        Please move/collect them into one folder, or run the script separately per folder.")
        sys.exit(1)

    mgf_abs = str(mgf_paths[0])
    out_abs = abs_path(args.outdir)
    param_file = abs_path(args.param_file) if args.param_file else None

    if param_file and (not Path(param_file).is_file()):
        eprint(f"[ERROR] param-file not found: {param_file}")
        sys.exit(1)

    ensure_dirs(out_abs)

    mgf_list_path = Path(out_abs, ".tips_tmp", "mgf_list.txt")
    mgf_list_path.write_text("\n".join([p.name for p in mgf_paths]) + "\n")

    print(f"[INFO] input dir : {mgf_dir}")
    print(f"[INFO] mgf count  : {len(mgf_paths)}")
    print(f"[INFO] first mgf : {mgf_abs}")
    print(f"[INFO] outdir    : {out_abs}")
    print(f"[INFO] sample    : {args.sample}")
    print("[INFO] images    :")
    print(f"  - {DEFAULT_IMG_CASANOVO}")
    print(f"  - {DEFAULT_IMG_PEPNET}")
    print(f"  - {DEFAULT_IMG_INSTANOVO}")
    print("")

    gpu_runtime_ok = docker_gpu_runtime_usable(DEFAULT_IMG_INSTANOVO)
    if gpu_runtime_ok:
        print("[INFO] GPU runtime: available ✅")
    else:
        print("[INFO] GPU runtime: NOT available -> will run CPU only ⚠️")
    print("")

    tools = [args.tool] if args.tool != "all" else ["casanovo", "pepnet", "instanovo"]

    for tool in tools:
        print("================================")
        print(f"[START] tool={tool}")
        print("================================")

        selection = tool_gpu_choice(tool, args.gpu)
        gpu_index = resolve_gpu_index_for_tool(
            gpu_runtime_ok=gpu_runtime_ok,
            img_for_check=DEFAULT_IMG_INSTANOVO,
            selection=selection,
        )

        if gpu_index is not None:
            print(f"[INFO] {tool} GPU selection: {selection} -> index={gpu_index} (cap<=8.9 verified in container)")
        else:
            print(f"[INFO] {tool} GPU selection: {selection} -> CPU")

        if tool == "casanovo":
            docker_run(
                img=DEFAULT_IMG_CASANOVO,
                script=SCRIPT_CASANOVO,
                mgf_abs=mgf_abs,
                out_abs=out_abs,
                sample=args.sample,
                gpu_index=gpu_index,
                param_file=param_file,
            )
        elif tool == "pepnet":
            docker_run(
                img=DEFAULT_IMG_PEPNET,
                script=SCRIPT_PEPNET,
                mgf_abs=mgf_abs,
                out_abs=out_abs,
                sample=args.sample,
                gpu_index=gpu_index,
                param_file=None,
            )
        elif tool == "instanovo":
            docker_run(
                img=DEFAULT_IMG_INSTANOVO,
                script=SCRIPT_INSTANOVO,
                mgf_abs=mgf_abs,
                out_abs=out_abs,
                sample=args.sample,
                gpu_index=gpu_index,
                param_file=None,
            )
        else:
            raise ValueError(f"Unknown tool: {tool}")

    print("\n==============================")
    print("[DONE] De novo docker run finished")
    for sub in ["casanovo", "pepnet", "instanovo", "logs"]:
        p = Path(out_abs, sub)
        if p.exists():
            print(f" - {p}")
    print("==============================\n")


if __name__ == "__main__":
    main()

