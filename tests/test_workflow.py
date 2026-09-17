"""Tests for the workflow text itself.

Not a Snakemake run — a handful of assertions about the shell the rules emit,
because that shell is the one part of the pipeline no other test executes and
the cheapest place to lose a whole cohort. Snakemake prefixes every shell block
with `set -euo pipefail`, so a habit that is harmless in an interactive shell is
fatal in a job, and the failure arrives as a one-line message inside a traceback
about the rule.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

SNAKEFILE = Path(__file__).resolve().parent.parent / "workflow" / "Snakefile"


class TestShellIsSafeUnderNounset:
    def snakefile(self) -> str:
        return SNAKEFILE.read_text()

    def test_pythonpath_tolerates_being_unset(self):
        """`export PYTHONPATH=...:$PYTHONPATH` under `set -u` is an error.

        PYTHONPATH is unset on a fresh compute node and set in the shell where
        anyone would test this by hand, so the bug hides exactly until it is
        expensive: all 101 jobs of a cohort exited in under a second with
        "PYTHONPATH: unbound variable".
        """
        text = self.snakefile()
        assert "$PYTHONPATH" in text, "the PYTHONPATH export has moved; update this test"
        assert ":+" in text or ":-" in text, (
            "PYTHONPATH is expanded without a default. Under `set -u` that is "
            "fatal whenever PYTHONPATH is unset, which is the normal state in a "
            "cluster job. Use ${PYTHONPATH:+:$PYTHONPATH}."
        )

    def test_the_expansion_actually_survives_nounset(self):
        """Assert the behaviour, not the spelling."""
        line = next(ln for ln in self.snakefile().splitlines()
                    if ln.strip().startswith("PYENV ="))
        namespace: dict = {"SRC": "/repo/src"}
        exec(line.strip(), namespace)          # noqa: S102 — our own file
        shell = namespace["PYENV"] + "echo ok"

        for env_desc, env in (("unset", {"PATH": "/usr/bin:/bin"}),
                              ("set", {"PATH": "/usr/bin:/bin",
                                       "PYTHONPATH": "/already/here"})):
            done = subprocess.run(["bash", "-c", "set -euo pipefail; " + shell],
                                  capture_output=True, text=True, env=env)
            assert done.returncode == 0, (
                f"PYENV fails under `set -u` with PYTHONPATH {env_desc}: "
                f"{done.stderr.strip()}"
            )
