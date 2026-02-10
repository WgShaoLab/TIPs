from __future__ import annotations

import shlex
import subprocess
from pathlib import Path
from typing import List, Optional
from dataclasses import dataclass, field


@dataclass
class DockerRunner:
    """
    Minimal docker runner that enforces:
      - same-path bind mount (-v S:S)
      - user mapping (-u uid:gid) to avoid permission issues
      - optional working directory (-w)
    """
    image: str
    mount_root: Path
    user_map: bool = True
    workdir_in_container: Optional[str] = "/tmp"
    extra_mounts: List[Path] = field(default_factory=list)

    def build_run_cmd(self, inner_cmd: List[str]) -> List[str]:
        cmd: List[str] = ["docker", "run", "--rm"]

        if self.user_map:
            uid = str(subprocess.check_output(["id", "-u"]).decode().strip())
            gid = str(subprocess.check_output(["id", "-g"]).decode().strip())
            cmd += ["-u", f"{uid}:{gid}"]

        # Bind mount with identical host/container path to simplify path handling
        m = str(self.mount_root)
        cmd += ["-v", f"{m}:{m}"]
        # Extra bind mounts (same-path host:container)
        for p in self.extra_mounts:
            pp = str(p)
            cmd += ["-v", f"{pp}:{pp}"]

        if self.workdir_in_container:
            cmd += ["-w", self.workdir_in_container]

        cmd += [self.image]
        cmd += inner_cmd
        return cmd

    @staticmethod
    def cmd_to_shell(cmd: List[str]) -> str:
        """
        Convert command list to a shell-safe string.
        """
        return " ".join(shlex.quote(x) for x in cmd)

    def run(self, inner_cmd: List[str], dry_run: bool = False, log_file: Optional[Path] = None) -> int:
        full_cmd = self.build_run_cmd(inner_cmd)
        shell_line = self.cmd_to_shell(full_cmd)

        if log_file is not None:
            log_file.parent.mkdir(parents=True, exist_ok=True)
            with log_file.open("a", encoding="utf-8") as f:
                f.write(shell_line + "\n")

        if dry_run:
            return 0

        proc = subprocess.run(full_cmd)
        return int(proc.returncode)
