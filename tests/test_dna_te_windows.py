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


TIRL_CHECK = r"""
import os, sys, tempfile
d = os.path.join(tempfile.mkdtemp(), "TIR-Learner3.0", "bin")
os.makedirs(d)
open(os.path.join(d, "post_processing.py"), "w").write(
    "import pandas as pd\n"
    "def combine_all(df_list):\n"
    "    df = pd.concat(df_list) if len(df_list) > 1 else df_list[0].copy()\n"
    "    return df.sort_values(['seqid', 'source'])\n")
open(os.path.join(d, "get_fasta_sequence.py"), "w").write(
    "import pandas as pd\n"
    "def get_fasta_pieces_SeqIO(genome_file, df_in, cpu_cores, flag_verbose):\n"
    "    return pd.concat([None])\n")
sys.path.insert(0, d)
import pandas as pd
import post_processing, get_fasta_sequence
empty_no_source = pd.DataFrame({"seqid": pd.Series([], dtype=object)})   # check_TIR_TSD's early return
print(len(post_processing.combine_all([empty_no_source])))
hit = pd.DataFrame({"seqid": ["c1"], "source": ["Module4"]})
print(len(post_processing.combine_all([empty_no_source, hit])))
out = get_fasta_sequence.get_fasta_pieces_SeqIO("g.fa", pd.DataFrame({"seqid": pd.Series([], dtype=object)}), 1, False)
print(list(out.columns), len(out))
"""


@unittest.skipIf(subprocess.run([sys.executable, "-c", "import pandas"], capture_output=True).returncode,
                 "pandas not installed")
class TirLearnerEmptySteps(unittest.TestCase):
    """TIR-Learner 3.0.x modules (loaded from a TIR-Learner path) get empty-step guards."""

    def test_combine_all_and_fasta_pieces(self):
        env = dict(os.environ, PYTHONPATH=os.path.join(SCRIPTS, "compat", "edta_python"))
        out = subprocess.run([sys.executable, "-c", TIRL_CHECK], capture_output=True, text=True, env=env)
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(out.stdout.split("\n")[:3], ["0", "1", "['seqid', 'seq'] 0"])


if __name__ == "__main__":
    unittest.main()
