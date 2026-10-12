"""Unit tests for dna_te_windows.py failure parsing and the TIR-Learner shim
(compat/edta_python/sitecustomize.py). Run: python3 -m unittest discover tests"""

import os
import subprocess
import sys
import unittest

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts")
sys.path.insert(0, SCRIPTS)

import dna_te_windows as dw  # noqa: E402

LOG = """Module 4, Step 8: Get FASTA sequences from CNN prediction
Module 4, Step 9: Check TIR and TSD
  Step 1/6: Checking TIR
Traceback (most recent call last):
  File "/usr/local/share/TIR-Learner3/TIR-Learner.py", line 119, in <module>
    main()
  File "/usr/local/share/TIR-Learner3/app/post_processing.py", line 22, in _combine_all
    df = df.sort_values(["seqid", "source"])
  File "/usr/local/lib/python3.12/site-packages/pandas/core/frame.py", line 7000, in sort_values
    raise KeyError(key)
KeyError: 'source'
Error: Could not load sequence. Empty file or bad format.
Error: TIR results not found!
""".splitlines(keepends=True)


class CrashSite(unittest.TestCase):
    def test_step_frame_error(self):
        self.assertEqual(dw.crash_site(LOG),
                         "Module 4, Step 9: Check TIR and TSD; post_processing.py:22; KeyError: 'source'")

    def test_block_stops_at_error_line(self):
        tb = dw.traceback_block(LOG)
        self.assertTrue(tb[0].startswith("Traceback"))
        self.assertEqual(tb[-1].strip(), "KeyError: 'source'")

    def test_no_traceback(self):
        self.assertEqual(dw.crash_site(["Error: TIR results not found!\n"]), "")
        self.assertEqual(dw.traceback_block([]), [])


SHIM_CHECK = """
import pandas as pd
df = pd.DataFrame({"id": pd.Series([], dtype=object), "TIR1": pd.Series([], dtype=object)})
df["check"] = df.apply(lambda x: "N" in x["TIR1"], axis=1)   # TIR-Learner's empty-frame pattern
print(len(df), list(df.columns))
full = pd.DataFrame({"a": [1, 2]})
print(type(full.apply(lambda r: r, axis=1)).__name__)       # non-empty frames unchanged
"""


@unittest.skipIf(subprocess.run([sys.executable, "-c", "import pandas"], capture_output=True).returncode,
                 "pandas not installed")
class Shim(unittest.TestCase):
    def test_empty_apply(self):
        env = dict(os.environ, PYTHONPATH=os.path.join(SCRIPTS, "compat", "edta_python"))
        out = subprocess.run([sys.executable, "-c", SHIM_CHECK], capture_output=True, text=True, env=env)
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(out.stdout.split("\n")[:2], ["0 ['id', 'TIR1', 'check']", "DataFrame"])


if __name__ == "__main__":
    unittest.main()
