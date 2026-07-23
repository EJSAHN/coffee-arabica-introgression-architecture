#!/usr/bin/env python3
"""Run a subprocess, mirror combined output to the console, and save a UTF-8 log."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) < 4 or sys.argv[1] != "--log":
        print("Usage: run_and_tee.py --log LOGFILE COMMAND [ARGS...]", file=sys.stderr)
        return 2

    log_path = Path(sys.argv[2]).resolve()
    command = sys.argv[3:]
    log_path.parent.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env.setdefault("PYTHONUTF8", "1")
    env.setdefault("PYTHONIOENCODING", "utf-8")

    print("RUN:", subprocess.list2cmdline(command), flush=True)
    print("LOG:", str(log_path), flush=True)

    with log_path.open("w", encoding="utf-8", newline="") as log:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=env,
            bufsize=0,
        )
        assert process.stdout is not None
        while True:
            chunk = process.stdout.readline()
            if not chunk:
                break
            text = chunk.decode("utf-8", errors="replace")
            sys.stdout.write(text)
            sys.stdout.flush()
            log.write(text)
            log.flush()
        return process.wait()


if __name__ == "__main__":
    raise SystemExit(main())
