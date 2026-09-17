"""Subprocess helpers.

Small on purpose. The one opinion encoded here is that a pipeline stage should
fail loudly with the tool's own stderr attached, because the failure modes in
this project — a CRAM whose reference is wrong, a remote read that times out, a
GATK jar that cannot find a `python` on PATH — all announce themselves clearly
in stderr and not at all in the exit code alone.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass


class ToolError(RuntimeError):
    """A tool exited non-zero, or a pipeline stage did."""


@dataclass
class Result:
    stdout: str
    stderr: str
    returncode: int


def require(*tools: str) -> None:
    """Fail early, with every missing tool named at once.

    One missing binary reported per run turns a cluster job into a guessing game;
    listing them all means one round trip.
    """
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise ToolError(
            "not on PATH: " + ", ".join(missing)
            + "\nLoad the environment first (see README: environment.yml, or the "
              "`prelude` setting for module-based sites)."
        )


def run(cmd: list[str], *, timeout: float | None = None,
        text_input: str | None = None, check: bool = True,
        env: dict | None = None, cwd: str | None = None) -> Result:
    """Run one command to completion, capturing both streams."""
    proc = subprocess.run(cmd, capture_output=True, text=True, input=text_input,
                          timeout=timeout, env=env, cwd=cwd)
    if check and proc.returncode != 0:
        raise ToolError(
            f"{cmd[0]} exited {proc.returncode}\n"
            f"  command: {' '.join(cmd)}\n"
            f"  stderr:  {proc.stderr.strip()[:2000]}"
        )
    return Result(proc.stdout, proc.stderr, proc.returncode)


def pipeline(stages: list[list[str]], *, env: dict | None = None,
             cwd: str | None = None, timeout: float | None = None,
             stdout_path: str | None = None) -> Result:
    """Run ``a | b | c``, failing if *any* stage fails.

    ``shell=True`` with ``set -o pipefail`` would be shorter, but it also makes
    every argument a quoting problem — and these arguments include URLs and
    region strings. Wiring the pipes explicitly keeps arguments as a list, and
    checking every stage's return code gives pipefail's semantics: a bowtie2
    killed by the scheduler must not be masked by a samtools that cheerfully
    writes a valid, truncated, empty BAM.

    Each stage's stderr goes to its own temporary file rather than a pipe. A
    pipe would deadlock the moment a chatty upstream stage — bowtie2 reports
    progress on stderr — filled its buffer while we were blocked reading the
    downstream one's stdout.
    """
    if not stages:
        raise ValueError("pipeline() needs at least one stage")

    with tempfile.TemporaryDirectory(prefix="lilrwgs_pipe_") as tmp:
        err_files = [open(os.path.join(tmp, f"stderr.{i}"), "w+") for i in range(len(stages))]
        sink = open(stdout_path, "wb") if stdout_path else None
        procs: list[subprocess.Popen] = []
        try:
            prev_stdout = None
            for i, cmd in enumerate(stages):
                last = i == len(stages) - 1
                proc = subprocess.Popen(
                    cmd,
                    stdin=prev_stdout,
                    stdout=(sink or subprocess.PIPE) if last else subprocess.PIPE,
                    stderr=err_files[i],
                    text=not last or sink is None,
                    env=env,
                    cwd=cwd,
                )
                # Hand the read end to the child and close ours, so an upstream
                # stage sees EOF when the downstream one exits rather than
                # blocking forever on a pipe nobody is draining.
                if prev_stdout is not None:
                    prev_stdout.close()
                prev_stdout = proc.stdout
                procs.append(proc)

            out = ""
            if sink is None:
                out, _ = procs[-1].communicate(timeout=timeout)
            else:
                procs[-1].wait(timeout=timeout)
            for proc in procs[:-1]:
                proc.wait()
        finally:
            if sink is not None:
                sink.close()

        errs = []
        for fh in err_files:
            fh.seek(0)
            errs.append(fh.read())
            fh.close()

        failed = [i for i, p in enumerate(procs) if p.returncode not in (0, None)]
        if failed:
            detail = "\n".join(
                f"  stage {i} ({' '.join(stages[i][:3])} ...) exited "
                f"{procs[i].returncode}: {errs[i].strip()[:800]}"
                for i in failed
            )
            raise ToolError("pipeline failed:\n" + detail)

        return Result(out, "\n".join(e for e in errs if e.strip()), 0)
