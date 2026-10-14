"""Unit tests for workflow/scripts/family_tandem.py: long-period arrays and
the copies-by-mass rule, on synthetic .out files. Run:
python3 -m unittest discover tests"""

import csv
import os
import random
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(__file__)
SCRIPTS = os.path.join(HERE, "..", "workflow", "scripts")
sys.path.insert(0, SCRIPTS)

import family_tandem as ft  # noqa: E402


def out_line(contig, begin, end, family, cons_len, strand="+", cls="Unknown"):
    n = end - begin + 1
    if strand == "+":
        rep = f"1 {n} ({max(cons_len - n, 0)})"
    else:
        rep = f"({max(cons_len - n, 0)}) {n} 1"
    return f"1000 10.0 0.0 0.0 {contig} {begin} {end} (0) {strand} {family} {cls} {rep} 1\n"


def random_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def run_family_tandem(families, hits, extra=()):
    """families: {name: cons_len}; hits: .out lines. Returns {family: row}."""
    with tempfile.TemporaryDirectory() as d:
        lib = os.path.join(d, "lib.fa")
        with open(lib, "w") as fh:
            for i, (name, n) in enumerate(families.items()):
                fh.write(f">{name}#Unknown\n{random_seq(n, i)}\n")
        out = os.path.join(d, "x.fa.out")
        with open(out, "w") as fh:
            fh.write("   SW  perc perc perc  query\n\n")
            fh.writelines(hits)
        stats = os.path.join(d, "stats.tsv")
        with open(stats, "w") as fh:
            fh.write("non_n_bp\n1000000000\n")
        fam_tsv, cls_tsv = os.path.join(d, "f.tsv"), os.path.join(d, "c.tsv")
        subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "family_tandem.py"), "--out-file", out,
             "--library", lib, "--assembly-stats", stats, "--arm", "shared", "--sample", "s",
             "--out", fam_tsv, "--class-out", cls_tsv, *extra],
            check=True, capture_output=True)
        with open(fam_tsv) as fh:
            return {r["family"]: r for r in csv.DictReader(fh, delimiter="\t")}


class PeriodicRuns(unittest.TestCase):
    def test_regular_spacing_is_a_run(self):
        arrays = [("c", 1 + i * 20000, 4000 + i * 20000, 1, "+") for i in range(6)]
        runs = ft.periodic_runs(arrays, 50000, 0.1, 4)
        self.assertEqual(len(runs), 1)
        self.assertEqual(len(runs[0][1]), 6)
        self.assertEqual(runs[0][2], 20000)

    def test_irregular_spacing_is_not(self):
        starts = [1, 7000, 31000, 36000, 80000, 93000]
        arrays = [("c", s, s + 3999, 1, "+") for s in starts]
        self.assertEqual(ft.periodic_runs(arrays, 50000, 0.1, 4), [])

    def test_spacing_above_max_period_is_not(self):
        arrays = [("c", 1 + i * 80000, 4000 + i * 80000, 1, "+") for i in range(6)]
        self.assertEqual(ft.periodic_runs(arrays, 50000, 0.1, 4), [])

    def test_strands_kept_apart(self):
        arrays = [("c", 1 + i * 20000, 4000 + i * 20000, 1, "+-"[i % 2]) for i in range(6)]
        self.assertEqual(ft.periodic_runs(arrays, 50000, 0.1, 4), [])


class FamilyTandem(unittest.TestCase):
    def test_partial_unit_family_is_tandem(self):
        # 4 kb consensus, one hit per 20 kb unit, 100 units on 2 contigs:
        # no chain spans 3 consensus lengths, but the spacing is regular
        hits = [out_line(c, 1 + i * 20000, 4000 + i * 20000, "long", 4000)
                for c in ("c1", "c2") for i in range(50)]
        rows = run_family_tandem({"long": 4000}, hits)
        r = rows["long"]
        self.assertEqual(r["short_array_bp"], "0")
        self.assertEqual(r["tandem_family"], "True")
        self.assertEqual(r["long_period"], "20000")

    def test_off_switch(self):
        hits = [out_line(c, 1 + i * 20000, 4000 + i * 20000, "long", 4000)
                for c in ("c1", "c2") for i in range(50)]
        rows = run_family_tandem({"long": 4000}, hits, ["--max-period", "0"])
        self.assertEqual(rows["long"]["tandem_family"], "False")

    def test_dispersed_family_is_not(self):
        rng = random.Random(1)
        hits = []
        for c in range(100):
            pos = 1
            for _ in range(3):
                pos += rng.randint(5000, 200000)
                hits.append(out_line(f"c{c}", pos, pos + 3999, "te", 4000))
        rows = run_family_tandem({"te": 4000}, hits)
        self.assertEqual(rows["te"]["tandem_family"], "False")
        self.assertEqual(rows["te"]["periodic_bp"], "0")

    def test_long_unit_counts_copies_by_mass(self):
        # 4 kb unit, one contiguous 60-unit array: 240 kb = 60 "copies"
        # (< min_copies 100) but >= 100 x 2000 bp of tandem mass
        hits = [out_line("c1", 1 + i * 4000, 4000 + i * 4000, "unit", 4000) for i in range(60)]
        r = run_family_tandem({"unit": 4000}, hits)["unit"]
        self.assertEqual(r["pass_copies"], "False")
        self.assertEqual(r["tandem_family"], "True")
        self.assertEqual(r["satellite_like"], "False")

    def test_small_long_unit_array_still_fails(self):
        hits = [out_line("c1", 1 + i * 4000, 4000 + i * 4000, "unit", 4000) for i in range(20)]
        r = run_family_tandem({"unit": 4000}, hits)["unit"]
        self.assertEqual(r["tandem_family"], "False")


if __name__ == "__main__":
    unittest.main()
