# utils.py
import os
import time
import subprocess
import functools
from typing import Callable

def run_cmd(cmd: str, env=None, shell: bool = True) -> None:
    """
    Run a shell command and raise if it fails.
    """
    if env is None:
        env = os.environ
    print(f"[CMD] {cmd}")
    result = subprocess.run(cmd, shell=shell, env=env)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with code {result.returncode}: {cmd}")

def wait_while_running(pattern: str, poll_interval: int = 5) -> None:
    """
    Block until no running processes match 'pattern' (used to wait for comet/blastp/etc).

    WARNING:
    - Uses 'ps aux | grep {pattern}' style inspection.
    - This is heuristic and assumes a Unix-like system.

    pattern: substring that identifies the external process (e.g. path to comet binary).
    """
    while True:
        try:
            out = subprocess.check_output(
                f"ps aux | grep -v grep | grep -c '{pattern}'",
                shell=True
            )
            n = int(out.strip())
        except Exception:
            # If ps or grep fails, assume process finished.
            n = 0

        if n <= 0:
            break
        time.sleep(poll_interval)


def timing(func: Callable) -> Callable:
    """
    Decorator to report execution time of a method, and also
    log the pipeline stage for easier debugging.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        self = args[0] if args else None
        class_name = self.__class__.__name__ if self else "<no_self>"
        soft_name = getattr(self, "Soft_name", None)
        if soft_name:
            print(f"[{class_name}] Running {func.__name__} for {soft_name}")

        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()

        print(f"[{class_name}] {func.__name__} finished in {end - start:.2f} sec")
        return result

    return wrapper
