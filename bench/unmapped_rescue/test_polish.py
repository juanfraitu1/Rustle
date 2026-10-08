#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_polish.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polish as P  # noqa: E402


def sam(name, pos1, cigar, seq, flag=0):
    return "\t".join([name, str(flag), "cons", str(pos1), "60", cigar, "*", "0", "0", seq, "*"])


class Binomial(unittest.TestCase):
    def test_upper_tail_known_values(self):
        self.assertAlmostEqual(P.binom_tail(8, 10, 0.5), 0.0546875, places=9)
        self.assertAlmostEqual(P.binom_tail(0, 10, 0.3), 1.0, places=12)
        self.assertAlmostEqual(P.binom_tail(10, 10, 0.5), 1 / 1024, places=12)

    def test_degenerate_error_rates(self):
        self.assertEqual(P.binom_tail(1, 5, 0.0), 0.0)
        self.assertEqual(P.binom_tail(0, 5, 0.0), 1.0)
        self.assertEqual(P.binom_tail(5, 5, 1.0), 1.0)


class Pileup(unittest.TestCase):
    def test_counts_matches_mismatch_deletion_and_insertion(self):
        cons = "ACGTACGT"
        # read: ACG [T->C mismatch at 3] A [CG deleted at 5,6] ... then T; with an insertion "GG" between positions 0 and 1
        recs = [sam("r1", 1, "1=2I2=1X1=", "AGGCGCA")]   # A | GG inserted before consensus position 1 | CG | T->C | A
        cnt, ins, ngap = P.pileup(recs, cons)
        self.assertEqual(cnt[0], [1, 0, 0, 0, 0])      # A
        self.assertEqual(cnt[1], [0, 1, 0, 0, 0])      # C
        self.assertEqual(cnt[2], [0, 0, 1, 0, 0])      # G
        self.assertEqual(cnt[3], [0, 1, 0, 0, 0])      # T replaced by C
        self.assertEqual(cnt[4], [1, 0, 0, 0, 0])      # A
        self.assertEqual(ins[1], {"GG": 1})
        self.assertEqual(sum(cnt[5]), 0)               # not covered

    def test_deletion_counts_a_dash(self):
        cons = "ACGTA"
        cnt, ins, ngap = P.pileup([sam("r1", 1, "2=1D2=", "ACTA")], cons)
        self.assertEqual(cnt[2], [0, 0, 0, 0, 1])
        self.assertEqual(cnt[3][3], 1)

    def test_secondary_and_unmapped_records_are_skipped(self):
        cons = "ACGT"
        cnt, ins, ngap = P.pileup([sam("r1", 1, "4=", "ACGT", flag=256), sam("r2", 1, "4=", "ACGT", flag=4)], cons)
        self.assertTrue(all(sum(c) == 0 for c in cnt))


def stack(cons, reads):
    """a pileup of identical full-length reads built from 'sub' lists: reads = [(count, {pos: base})]"""
    recs = []
    for k, (n, subs) in enumerate(reads):
        s = "".join(subs.get(i, b) for i, b in enumerate(cons))
        cig = "".join(("1X" if i in subs and subs[i] != b else "1=") for i, b in enumerate(cons))
        recs += [sam(f"r{k}_{j}", 1, cig, s) for j in range(n)]
    return recs


class Corrections(unittest.TestCase):
    def test_a_consensus_error_held_by_all_reads_is_corrected(self):
        cons = "ACGTACGTACGTACGTACGT"
        recs = stack(cons, [(30, {5: "G"})])
        pile = P.pileup(recs, cons)
        applied, minority, e = P.corrections(cons, pile)
        self.assertEqual(applied, [("sub", 5, "G")])
        self.assertEqual(P.apply(cons, applied), "ACGTAGGTACGTACGTACGT")

    def test_a_significant_minority_variant_is_recorded_not_applied(self):
        cons = "ACGTACGTACGTACGTACGT"
        recs = stack(cons, [(24, {}), (6, {5: "G"})])
        applied, minority, e = P.corrections(cons, P.pileup(recs, cons))
        self.assertEqual(applied, [])
        self.assertEqual([m[:2] for m in minority], [("sub", 5)])

    def test_scattered_single_read_errors_change_nothing(self):
        cons = "ACGTACGTACGTACGTACGT"
        recs = stack(cons, [(27, {}), (1, {2: "T"}), (1, {9: "A"}), (1, {14: "T"})])
        applied, minority, e = P.corrections(cons, P.pileup(recs, cons))
        self.assertEqual((applied, minority), ([], []))

    def test_a_missing_base_is_inserted_and_an_extra_base_is_deleted(self):
        cons = "ACGTACGTACGT"
        # 26 of 30 reads carry an extra 'A' before position 4 (the consensus lacks it); the other 4 do not
        with_ins = [sam(f"a{j}", 1, "4=1I8=", "ACGT" + "A" + "ACGTACGT") for j in range(26)]
        without = [sam(f"b{j}", 1, "12=", cons) for j in range(4)]
        applied, minority, e = P.corrections(cons, P.pileup(with_ins + without, cons))
        self.assertEqual(applied, [("ins", 4, "A")])
        self.assertEqual(P.apply(cons, applied), "ACGTAACGTACGT")
        # 28 of 30 reads delete position 6
        dele = [sam(f"d{j}", 1, "6=1D5=", cons[:6] + cons[7:]) for j in range(28)]
        keep = [sam(f"k{j}", 1, "12=", cons) for j in range(2)]
        applied, minority, e = P.corrections(cons, P.pileup(dele + keep, cons))
        self.assertEqual(applied, [("del", 6, None)])
        self.assertEqual(P.apply(cons, applied), cons[:6] + cons[7:])

    def test_apply_several_edits_keeps_coordinates(self):
        cons = "ACGTACGT"
        self.assertEqual(P.apply(cons, [("sub", 0, "T"), ("ins", 3, "GG"), ("del", 6, None)]), "TCGGGTACT")


class Loop(unittest.TestCase):
    def test_polish_stops_when_a_round_changes_nothing_and_reports_rounds(self):
        truth = "ACGTACGTACGTACGTACGT"
        wrong = "ACGTAGGTACGTACGTACGT"
        reads = [truth] * 30

        def align(c):
            out = []
            for j, r in enumerate(reads):
                cig = "".join("1=" if a == b else "1X" for a, b in zip(c, r))
                out.append(sam(f"r{j}", 1, cig, r))
            return out

        cons, rounds, applied, minority = P.polish(wrong, align)
        self.assertEqual(cons, truth)
        self.assertEqual(rounds, 2)          # round 1 corrects, round 2 finds nothing
        self.assertEqual(len(applied), 1)

    def test_polish_caps_the_rounds(self):
        flip = {"A": "C", "C": "A"}
        state = {"n": 0}

        def align(c):
            state["n"] += 1
            r = flip[c[0]] + c[1:]       # every round the reads disagree with the consensus at its first base, whatever it is
            cig = "1X" + "".join("1=" for _ in c[1:])
            return [sam(f"r{j}", 1, cig, r) for j in range(30)]

        cons, rounds, applied, minority = P.polish("ACACACACAC", align, max_rounds=3)
        self.assertEqual((rounds, state["n"]), (3, 3))


class Cs(unittest.TestCase):
    def test_edit_columns_forward_strand(self):
        # 5 matches, a mismatch, 2 matches, an insertion of 'gg' in the query, 1 match, a deletion of 'tt' from the reference, 3 matches
        cols = P.cs_edit_columns(":5*ag:2+gg:1-tt:3", qstart=10, qend=10 + 5 + 1 + 2 + 2 + 1 + 3, qlen=100, strand="+")
        self.assertEqual(cols, {15, 18, 19, 21})

    def test_edit_columns_reverse_strand_are_reported_in_consensus_coordinates(self):
        cols = P.cs_edit_columns(":5*ag:2", qstart=10, qend=18, qlen=100, strand="-")
        # the alignment covers consensus positions 10..17; on the reverse strand it starts at the end: the mismatch is the 6th aligned base
        self.assertEqual(cols, {17 - 5})

    def test_splice_introns_do_not_advance_the_query(self):
        cols = P.cs_edit_columns(":3~gt100ag:2*ag", qstart=0, qend=6, qlen=50, strand="+")
        self.assertEqual(cols, {5})


class Allelic(unittest.TestCase):
    def test_share_of_edit_columns_that_are_minority_columns_within_one_base(self):
        self.assertAlmostEqual(P.allele_like_share({10, 20, 30, 40}, [("sub", 11, "A"), ("sub", 29, "C")]), 0.5)
        self.assertIsNone(P.allele_like_share(set(), [("sub", 11, "A")]))


class SplicedGaps(unittest.TestCase):
    def test_n_gaps_count_as_deletions_only_on_request(self):
        cons = "A" * 20
        line = sam("r1", 1, "5=10N5=", "AAAAAAAAAA")
        cnt, _ins, _ng = P.pileup([line], cons)
        self.assertTrue(all(c[4] == 0 for c in cnt))
        cnt, _ins, _ng = P.pileup([line], cons, n_as_del=True)
        self.assertEqual([c[4] for c in cnt[5:15]], [1] * 10)
        self.assertEqual(sum(cnt[3]), 1)


if __name__ == "__main__":
    unittest.main()
