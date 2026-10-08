"""Unit tests for workflow/scripts/numt_calls.py: origin-wrap folding and the
hit merge logic, on synthetic hits. Run: python3 -m unittest discover tests"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts"))

import numt_calls as nc  # noqa: E402

L = 16000


def hit(contig, start, end, strand, ms, me, wraps=False, pident=90.0, bitscore=None, mito="mt"):
    h = nc.Hit()
    h.mito, h.contig, h.start, h.end, h.strand = mito, contig, start, end, strand
    h.pident, h.aln_len, h.mito_len = pident, end - start, L
    h.mito_start, h.mito_end, h.wraps = ms, me, wraps
    h.bitscore = bitscore if bitscore is not None else float(end - start)
    return h


class TestFold(unittest.TestCase):
    def test_first_copy(self):
        self.assertEqual(nc.fold(101, 600, L), (101, 600, False))

    def test_second_copy(self):
        self.assertEqual(nc.fold(L + 101, L + 600, L), (101, 600, False))

    def test_across_origin(self):
        self.assertEqual(nc.fold(L - 99, L + 200, L), (L - 99, 200, True))

    def test_ends_exactly_at_origin(self):
        self.assertEqual(nc.fold(L - 99, L, L), (L - 99, L, False))

    def test_starts_exactly_at_origin(self):
        self.assertEqual(nc.fold(L + 1, L + 50, L), (1, 50, False))

    def test_longer_than_mitogenome(self):
        ms, me, wraps = nc.fold(1, L + 10, L)
        self.assertTrue(wraps)
        self.assertEqual((ms, me), (1, 10))


class TestMitoStep(unittest.TestCase):
    def test_plus_adjacent_gap_overlap(self):
        a = hit("c", 0, 500, "+", 1001, 1500)
        self.assertEqual(nc.mito_step(a, hit("c", 600, 900, "+", 1501, 1800)), 0)
        self.assertEqual(nc.mito_step(a, hit("c", 600, 900, "+", 1551, 1850)), 50)
        self.assertEqual(nc.mito_step(a, hit("c", 600, 900, "+", 1481, 1780)), -20)

    def test_minus_strand_runs_backwards(self):
        # genome order a, b; on '-' the mitogenome runs b.end -> a.start
        a = hit("c", 0, 500, "-", 1501, 2000)
        b = hit("c", 3500, 4000, "-", 1001, 1490)
        self.assertEqual(nc.mito_step(a, b), 10)

    def test_across_origin(self):
        a = hit("c", 0, 300, "+", L - 299, L)
        b = hit("c", 400, 700, "+", 21, 320)
        self.assertEqual(nc.mito_step(a, b), 20)
        a = hit("c", 0, 300, "-", 1, 300)
        b = hit("c", 400, 700, "-", L - 299, L - 10)
        self.assertEqual(nc.mito_step(a, b), 10)


class TestDuplicates(unittest.TestCase):
    def test_doubled_query_copies_collapse(self):
        full = hit("c", 1000, 1600, "+", L - 299, 300, wraps=True, bitscore=900)
        edge1 = hit("c", 1300, 1600, "+", 1, 300, bitscore=450)      # query copy 1, start
        edge2 = hit("c", 1000, 1300, "+", L - 299, L, bitscore=450)  # query copy 2, end
        same = hit("c", 1000, 1600, "+", L - 299, 300, wraps=True, bitscore=900)
        kept = nc.drop_duplicates([edge1, full, edge2, same], 0.5)
        self.assertEqual(len(kept), 1)
        self.assertTrue(kept[0].wraps)

    def test_neighbours_and_other_strand_survive(self):
        a = hit("c", 0, 500, "+", 1, 500)
        b = hit("c", 450, 1000, "+", 501, 1050)  # 50 bp overlap only
        c = hit("c", 0, 500, "-", 5001, 5500)
        self.assertEqual(len(nc.drop_duplicates([a, b, c], 0.5)), 3)


class TestMerge(unittest.TestCase):
    def chains(self, hits, merge_gap=500, mito_gap=100):
        return sorted((sorted((h.start for h in ch)) for ch in nc.merge(hits, merge_gap, mito_gap)))

    def test_split_by_insertion(self):
        a = hit("c", 0, 1000, "+", 2001, 3000)
        b = hit("c", 4000, 5000, "+", 3001, 4000)  # 3 kb TE between
        self.assertEqual(self.chains([a, b]), [[0], [4000]])
        self.assertEqual(self.chains([a, b], merge_gap=5000), [[0, 4000]])

    def test_not_colinear(self):
        a = hit("c", 0, 1000, "+", 2001, 3000)
        b = hit("c", 1100, 2000, "+", 8001, 9000)
        self.assertEqual(self.chains([a, b]), [[0], [1100]])

    def test_mito_gap_limit(self):
        a = hit("c", 0, 1000, "+", 2001, 3000)
        self.assertEqual(self.chains([a, hit("c", 1100, 2000, "+", 3101, 4000)]), [[0, 1100]])
        self.assertEqual(self.chains([a, hit("c", 1100, 2000, "+", 3102, 4000)]), [[0], [1100]])
        self.assertEqual(self.chains([a, hit("c", 1100, 2000, "+", 2950, 3800)]), [[0, 1100]])  # 51 bp overlap

    def test_strand_and_contig_separate(self):
        a = hit("c", 0, 1000, "+", 2001, 3000)
        self.assertEqual(len(self.chains([a, hit("c", 1100, 2000, "-", 3001, 4000)])), 2)
        self.assertEqual(len(self.chains([a, hit("d", 1100, 2000, "+", 3001, 4000)])), 2)

    def test_minus_strand(self):
        a = hit("c", 0, 1000, "-", 3001, 4000)
        b = hit("c", 1200, 2000, "-", 2001, 3000)
        self.assertEqual(self.chains([a, b]), [[0, 1200]])
        comp = nc.Compound(nc.merge([a, b], 500, 100)[0])
        self.assertEqual((comp.mito_start, comp.mito_end, comp.wraps), (2001, 4000, False))
        self.assertEqual(comp.max_gap, (200, 1000, 1200))

    def test_origin_join(self):
        a = hit("c", 0, 300, "+", L - 299, L)
        b = hit("c", 350, 700, "+", 1, 350)
        chains = nc.merge([a, b], 500, 100)
        self.assertEqual(len(chains), 1)
        comp = nc.Compound(chains[0])
        self.assertEqual((comp.mito_start, comp.mito_end, comp.wraps), (L - 299, 350, True))
        self.assertEqual(comp.aligned_bp, 650)

    def test_interleaved_unrelated_hit(self):
        # an unrelated hit between two colinear pieces doesn't break the chain
        a = hit("c", 0, 1000, "+", 2001, 3000)
        x = hit("c", 1050, 1150, "+", 9001, 9100)
        b = hit("c", 1200, 2000, "+", 3001, 3800)
        self.assertEqual(self.chains([a, x, b]), [[0, 1200], [1050]])

    def test_compound_identity_is_bp_weighted(self):
        a = hit("c", 0, 1000, "+", 1, 1000, pident=90.0)
        b = hit("c", 1000, 1100, "+", 1001, 1100, pident=80.0)
        comp = nc.Compound(nc.merge([a, b], 500, 100)[0])
        self.assertAlmostEqual(comp.identity, (90 * 1000 + 80 * 100) / 1100)
        self.assertEqual(comp.max_gap[0], 0)


class TestNeighbourPairs(unittest.TestCase):
    def test_any_gap_reported(self):
        a = hit("c", 0, 1000, "+", 2001, 3000)
        b = hit("c", 50000, 51000, "+", 3001, 4000)
        self.assertEqual([g for _a, _b, g in nc.neighbour_pairs([a, b], 100)], [49000])


class TestMitoContigs(unittest.TestCase):
    def test_threshold(self):
        hits = [hit("m", 0, 9000, "+", 1, 9000, pident=99.5), hit("m", 9000, 16000, "+", 9001, L, pident=99.0),
                hit("n", 0, 9000, "+", 1, 9000, pident=99.5)]
        mito = nc.mito_contigs(hits, {"m": 16500, "n": 2_000_000}, 0.8, 98)
        self.assertEqual(set(mito), {"m"})


class TestBlastNames(unittest.TestCase):
    def test_accession_tags_resolved(self):
        from fasta_utils import blast_subject_resolver
        r = blast_subject_resolver(["NC_002639.1", "ptg001005l", "scaf_1"])
        self.assertEqual(r("ref|NC_002639.1|"), "NC_002639.1")
        self.assertEqual(r("lcl|ptg001005l"), "ptg001005l")
        self.assertEqual(r("scaf_1"), "scaf_1")
        with self.assertRaises(ValueError):
            r("ref|NC_999999.1|")


if __name__ == "__main__":
    unittest.main()
