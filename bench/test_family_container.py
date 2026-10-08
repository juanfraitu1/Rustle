#!/usr/bin/env python3
"""Unit tests for bench/family_container.py (docs/archive/2026-09/PREREG_fusion_container_sim_2026-09-28.md §5 step 1).

    python3 bench/test_family_container.py        (stdlib unittest; fixtures are written to a temporary directory)

The fixture (GFF 1-based closed; PAF offsets are 0-based within each record, record key START = offset 0):

  family MCL0: A = c1:1001-2000  gA  tA1 1001-1100,1301-1400,1901-2000 ; tA2 1051-1150,1901-1950
                 + folded record A2 = c1:1951-2100 (gA2, one exon 1951-2100)  => blocks 1001-1150, 1301-1400, 1901-2100
               B = c1:5001-6000  gB  tB1 5001-5150, 5301-5400, 5701-6000
  family MCL1: D = c2:1001-2000  gD  tD1 1001-1200, 1601-2000
               E = c2:5001-6000  gE  tE1 5001-5200, 5501-5600, 5901-6000
  unclustered: U = c3:1001-2000  gU  tU1 1001-1100, 1901-2000

  R1  A->B  +  q[0,150)    t[0,150)   50M2I48M2D50M     A.b1/B.b1 exon-exon, 148 aligned columns each side
  R2  A->B  +  q[399,900)  t[399,900) 501M              touches A.b2/B.b2 by exactly 1 base; ends 0 bases before A.b3
  R3  A->D  -  q[700,1000) t[100,400) 100M5D95M5I100M   reverse strand: only the first 100M is exon-exon, it joins
                                                        D[100,200) to A query [900,1000) (a forward walk would put
                                                        it on A's intron [700,800))
  R4  B->A  +  q[700,760)  t[450,510) 60M               B.b3 onto A's INTRON: not evidence
  R5  B->U  +  q[800,900)  t[0,100)   100M              U unclustered: skipped
  R6  A2->D +  q[100,150)  t[650,700) 50M               the folded record counts for A (A.b3 genome 2051-2100)
  R7  A->A2 +  q[950,1000) t[0,50)    50M               same locus: skipped
  R8  D->E  +  q[600,650)  t[0,50)    50M               D.b2 / E.b1 core
  R9  A->U  +  q[900,1000) t[900,1000) 100M             U unclustered: skipped
  R10 E->B  +  q[500,600)  t[0,100)   100M              E.b2 accessory related to MCL0 (B.b1 is core already)
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family_container as fc  # noqa: E402

GENES = {  # gene_id: (chrom, strand, {transcript: [(start1, end1), ...]})
    'gA': ('c1', '+', {'tA1': [(1001, 1100), (1301, 1400), (1901, 2000)], 'tA2': [(1051, 1150), (1901, 1950)]}),
    'gA2': ('c1', '+', {'tA2x': [(1951, 2100)]}),
    'gB': ('c1', '+', {'tB1': [(5001, 5150), (5301, 5400), (5701, 6000)]}),
    'gD': ('c2', '-', {'tD1': [(1001, 1200), (1601, 2000)]}),
    'gE': ('c2', '-', {'tE1': [(5001, 5200), (5501, 5600), (5901, 6000)]}),
    'gU': ('c3', '+', {'tU1': [(1001, 1100), (1901, 2000)]}),
}
KEY = {'gA': 'c1:1001-2000', 'gA2': 'c1:1951-2100', 'gB': 'c1:5001-6000', 'gD': 'c2:1001-2000',
       'gE': 'c2:5001-6000', 'gU': 'c3:1001-2000'}
CLUSTERS = [('MCL0', 'gA'), ('MCL0', 'gB'), ('MCL1', 'gD'), ('MCL1', 'gE')]
FOLDS = [('gA2', 'gA')]
PAF = [  # qgene, qs, qe, strand, tgene, ts, te, cigar
    ('gA', 0, 150, '+', 'gB', 0, 150, '50M2I48M2D50M'),
    ('gA', 399, 900, '+', 'gB', 399, 900, '501M'),
    ('gA', 700, 1000, '-', 'gD', 100, 400, '100M5D95M5I100M'),
    ('gB', 700, 760, '+', 'gA', 450, 510, '60M'),
    ('gB', 800, 900, '+', 'gU', 0, 100, '100M'),
    ('gA2', 100, 150, '+', 'gD', 650, 700, '50M'),
    ('gA', 950, 1000, '+', 'gA2', 0, 50, '50M'),
    ('gD', 600, 650, '+', 'gE', 0, 50, '50M'),
    ('gA', 900, 1000, '+', 'gU', 900, 1000, '100M'),
    ('gE', 500, 600, '+', 'gB', 0, 100, '100M'),
]


def key_len(g):
    k = fc.parse_key(KEY[g])
    return k[2] - k[1] + 1


def write_fixture(d, paf=PAF, clusters=CLUSTERS, folds=FOLDS, with_cg=True):
    gtf = os.path.join(d, 'x.gtf')
    with open(gtf, 'w') as fh:
        for g, (c, st, txs) in GENES.items():
            for t, exs in txs.items():
                fh.write(f'{c}\trustle\ttranscript\t{exs[0][0]}\t{exs[-1][1]}\t.\t{st}\t.\tgene_id "{g}"; '
                         f'transcript_id "{t}"; reads "5";\n')
                for i, (s, e) in enumerate(exs):
                    fh.write(f'{c}\trustle\texon\t{s}\t{e}\t.\t{st}\t.\tgene_id "{g}"; transcript_id "{t}"; '
                             f'exon_number "{i + 1}";\n')
    fam = os.path.join(d, 'x.fam')
    with open(fam + '.loci.gff3', 'w') as fh:
        fh.write('##gff-version 3\n')
        for g, (c, st, txs) in GENES.items():
            k = fc.parse_key(KEY[g])
            fh.write(f'{c}\t.\tgene\t{k[1]}\t{k[2]}\t.\t{st}\t.\tID=gene-{g};Name={g}\n')
            for s, e in list(txs.values())[0]:
                fh.write(f'{c}\t.\texon\t{s}\t{e}\t.\t{st}\t.\tParent=gene-{g};gene={g}\n')
    with open(fam + '.clusters.tsv', 'w') as fh:
        fh.write('cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n')
        for fid, g in clusters:
            k = fc.parse_key(KEY[g])
            n = sum(1 for f, _ in clusters if f == fid)
            fh.write(f'{fid}\t{n}\t1.0000\t1.0000\tNA\t{k[0]}\t{k[1]}\t{k[2]}\n')
    with open(fam + '.loci.tsv', 'w') as fh:
        fh.write('annotation\trepresentative\n')
        for a, r in folds:
            fh.write(f'{KEY[a]}\t{KEY[r]}\n')
    with open(fam + '.loci.paf', 'w') as fh:
        for qg, qs, qe, st, tg, ts, te, cg in PAF if paf is None else paf:
            tags = f'\tNM:i:0\ttp:A:P\tcg:Z:{cg}' if with_cg else '\tNM:i:0\ttp:A:P'
            fh.write(f'{KEY[qg]}\t{key_len(qg)}\t{qs}\t{qe}\t{st}\t{KEY[tg]}\t{key_len(tg)}\t{ts}\t{te}\t'
                     f'{qe - qs}\t{max(qe - qs, te - ts)}\t60{tags}\n')
    return gtf, fam


class Primitives(unittest.TestCase):
    def test_merge_blocks_overlap_merges_abutting_stays_separate(self):
        self.assertEqual(fc.merge_blocks([(10, 20), (0, 10), (15, 30), (40, 50), (45, 46)]),
                         [(0, 10), (10, 30), (40, 50)])

    def test_union_len(self):
        self.assertEqual(fc.union_len([(0, 10), (10, 20), (15, 25), (30, 31)]), 26)
        self.assertEqual(fc.union_len([]), 0)

    def test_aligned_runs_forward_with_eq_x_and_n(self):
        runs = list(fc.aligned_runs('10=2X3N4I5M1D6M', '+', 100, 127, 50, 77))
        self.assertEqual(runs, [(50, 100, 10), (60, 110, 2), (65, 116, 5), (71, 121, 6)])

    def test_aligned_runs_reverse(self):
        # '-': query consumed from qe downwards; run query interval [qe - consumed - n, qe - consumed)
        runs = list(fc.aligned_runs('100M5D95M5I100M', '-', 700, 1000, 100, 400))
        self.assertEqual(runs, [(100, 900, 100), (205, 805, 95), (300, 700, 100)])

    def test_cigar_length_mismatch_and_bad_ops_raise(self):
        with self.assertRaises(ValueError):
            list(fc.aligned_runs('10M', '+', 0, 11, 0, 10))
        with self.assertRaises(ValueError):
            list(fc.aligned_runs('5S10M', '+', 0, 10, 0, 10))
        with self.assertRaises(ValueError):
            list(fc.aligned_runs('10M', '.', 0, 10, 0, 10))

    def test_exon_columns_reverse_maps_each_column_exactly(self):
        # query record offset 0 = genome 0; target offset 0 = genome 1000. Query exon [3,5), target exon [1000,1002):
        # '-' 10M on q[0,10) / t[0,10): column i joins t 1000+i with q 9-i, so t [1000,1002) <-> q [8,10) -- not an
        # exon on the query -- and q [3,5) <-> t [1005,1007) -- not an exon on the target. No exon-exon column.
        q_blocks, t_blocks = ([3], [5]), ([1000], [1002])
        self.assertEqual(list(fc.exon_columns('10M', '-', 0, 10, 0, 10, 0, 1000, q_blocks, t_blocks)), [])
        # a target exon at [1005,1007) is exactly the mirror of the query exon: two columns
        t_blocks = ([1005], [1007])
        self.assertEqual(list(fc.exon_columns('10M', '-', 0, 10, 0, 10, 0, 1000, q_blocks, t_blocks)),
                         [(3, 5, 0, 1005, 1007, 0)])
        # and on '+' the same exons do not meet (q [3,5) <-> t [1003,1005))
        self.assertEqual(list(fc.exon_columns('10M', '+', 0, 10, 0, 10, 0, 1000, q_blocks, t_blocks)), [])


class Container(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR'))
        cls.gtf, cls.fam = write_fixture(cls.tmp.name)
        cls.out = os.path.join(cls.tmp.name, 'out')
        rc = fc.main(['--gtf', cls.gtf, '--fam', cls.fam, '--out', cls.out])
        assert rc == 0
        cls.rows = {}
        with open(cls.out + '.blocks.tsv') as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
            cls.header = hdr
            for line in fh:
                r = dict(zip(hdr, line.rstrip('\n').split('\t')))
                cls.rows[(r['locus'], int(r['block']))] = r
        with open(cls.out + '.relations.tsv') as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
            cls.rel = [dict(zip(hdr, line.rstrip('\n').split('\t'))) for line in fh]
        with open(cls.out + '.summary.tsv') as fh:
            fh.readline()
            cls.summary = {k: int(v) for k, v in (line.rstrip('\n').split('\t') for line in fh)}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def row(self, g, b):
        return self.rows[(KEY[g], b)]

    def test_header(self):
        self.assertEqual(self.header, fc.BLOCK_HEADER)

    def test_multi_transcript_and_folded_record_union(self):
        a = [self.row('gA', b) for b in (1, 2, 3)]
        self.assertEqual([(int(r['start']), int(r['end'])) for r in a], [(1001, 1150), (1301, 1400), (1901, 2100)])
        self.assertEqual(a[0]['gene_ids'], 'gA,gA2')
        self.assertEqual(a[0]['n_records'], '2')
        self.assertEqual(a[0]['n_blocks'], '3')
        self.assertNotIn((KEY['gA2'], 1), self.rows, 'a folded record is part of its locus, not a row of its own')

    def test_forward_cigar_with_insertion_and_deletion(self):
        a1, b1 = self.row('gA', 1), self.row('gB', 1)
        self.assertEqual((a1['class'], a1['core_bp'], a1['core_partners']), ('core', '148', f'{KEY["gB"]}=148'))
        self.assertEqual((b1['class'], b1['core_bp']), ('core', '148'))
        self.assertEqual(a1['rel_families'], '.')

    def test_one_base_touch_is_core_zero_is_accessory(self):
        self.assertEqual((self.row('gA', 2)['class'], self.row('gA', 2)['core_bp']), ('core', '1'))
        self.assertEqual((self.row('gB', 2)['class'], self.row('gB', 2)['core_bp']), ('core', '1'))
        a3 = self.row('gA', 3)  # R2 ends one base before it; its partner base there IS a B exon base
        self.assertEqual(a3['class'], 'accessory')
        self.assertEqual(a3['core_partners'], '.')

    def test_reverse_strand_projection_and_relation(self):
        a3 = self.row('gA', 3)
        self.assertEqual(a3['rel_families'], 'MCL1')
        self.assertEqual(a3['rel_bp'], '150')  # R3's 100 (genome 1901-2000) + R6's 50 via the folded record
        self.assertEqual(a3['rel_partners'], f'MCL1|{KEY["gD"]}=150')
        d1 = self.row('gD', 1)
        self.assertEqual((d1['class'], d1['rel_families'], d1['rel_bp']), ('accessory', 'MCL0', '100'))
        self.assertEqual(d1['rel_partners'], f'MCL0|{KEY["gA"]}=100')

    def test_partner_intron_is_not_evidence_and_unclustered_is_ignored(self):
        b3 = self.row('gB', 3)
        self.assertEqual((b3['class'], b3['rel_families'], b3['rel_bp']), ('accessory', '.', '0'))

    def test_core_blocks_are_not_tested_for_relations(self):
        d2, e1, b1 = self.row('gD', 2), self.row('gE', 1), self.row('gB', 1)
        self.assertEqual((d2['class'], d2['core_bp'], d2['rel_families']), ('core', '50', '.'))
        self.assertEqual((e1['class'], e1['core_bp']), ('core', '50'))
        self.assertEqual(b1['rel_families'], '.')  # R10 hits B.b1 from MCL1 but B.b1 is core

    def test_accessory_without_any_alignment(self):
        e3 = self.row('gE', 3)
        self.assertEqual((e3['class'], e3['rel_families'], e3['bp']), ('accessory', '.', '100'))
        e2 = self.row('gE', 2)
        self.assertEqual((e2['class'], e2['rel_families'], e2['rel_bp']), ('accessory', 'MCL0', '100'))

    def test_unclustered_locus_has_no_rows(self):
        self.assertFalse(any(k[0] == KEY['gU'] for k in self.rows))
        self.assertEqual(len(self.rows), 3 + 3 + 2 + 3)

    def test_family_relation_table(self):
        got = [(r['family_id'], r['related_family'], r['n_members'], r['n_blocks'], r['bp'], r['reciprocal'],
                r['members']) for r in self.rel]
        self.assertEqual(got, [('MCL0', 'MCL1', '1', '1', '150', 'yes', KEY['gA']),
                               ('MCL1', 'MCL0', '2', '2', '200', 'yes', f'{KEY["gD"]},{KEY["gE"]}')])

    def test_summary_counts(self):
        s = self.summary
        self.assertEqual((s['families'], s['loci'], s['folded_records'], s['blocks']), (2, 4, 1, 11))
        self.assertEqual((s['core_blocks'], s['accessory_blocks'], s['accessory_blocks_related']), (6, 5, 3))
        self.assertEqual((s['paf_records'], s['paf_unclustered'], s['paf_same_locus'], s['paf_no_exon_interval']),
                         (10, 2, 1, 1))
        self.assertEqual((s['paf_projected'], s['paf_exon_exon']), (6, 6))
        self.assertEqual((s['loci_all_core'], s['loci_with_accessory']), (0, 4))
        self.assertEqual(s['records_span_mismatch'], 0)


class Guards(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR'))

    def tearDown(self):
        self.tmp.cleanup()

    def test_record_without_cigar_raises(self):
        gtf, fam = write_fixture(self.tmp.name, with_cg=False)
        with self.assertRaises(ValueError):
            fc.run(gtf, fam + '.clusters.tsv', fam + '.loci.gff3', fam + '.loci.tsv', fam + '.loci.paf')

    def test_member_in_two_families_raises(self):
        gtf, fam = write_fixture(self.tmp.name, clusters=CLUSTERS + [('MCL1', 'gA')])
        with self.assertRaises(ValueError):
            fc.run(gtf, fam + '.clusters.tsv', fam + '.loci.gff3', fam + '.loci.tsv', fam + '.loci.paf')

    def test_without_the_fold_the_folded_record_is_unclustered(self):
        gtf, fam = write_fixture(self.tmp.name, folds=[])
        rows, rel, cnt = fc.run(gtf, fam + '.clusters.tsv', fam + '.loci.gff3', fam + '.loci.tsv',
                                fam + '.loci.paf')
        a = [r for r in rows if r['locus'] == KEY['gA']]
        self.assertEqual([(r['start'], r['end']) for r in a], [(1001, 1150), (1301, 1400), (1901, 2000)])
        self.assertEqual(a[2]['rel_bp'], 100)  # R6 (from the no-longer-folded record) no longer counts
        self.assertEqual(cnt['paf_same_locus'], 0)
        self.assertEqual(cnt['paf_unclustered'], 4)  # R5, R6, R7, R9


if __name__ == '__main__':
    unittest.main(verbosity=2)
