from __future__ import annotations

import os
import time
from pathlib import Path


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S", time.localtime())


def is_abs_path(p: Path) -> bool:
    return p.is_absolute()


def rel_to(child: Path, parent: Path) -> str:
    """
    Return a human-friendly relative path if possible, otherwise absolute.
    """
    try:
        return str(child.relative_to(parent))
    except Exception:
        return str(child)


def env_flag(name: str, default: bool = False) -> bool:
    v = os.environ.get(name, "")
    if v == "":
        return default
    return v.lower() in {"1", "true", "yes", "y", "on"}
