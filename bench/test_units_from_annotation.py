#!/usr/bin/env python3
"""Tests of units_from_annotation.py on hand-made GFF / GTF files (python3 bench/test_units_from_annotation.py). Stdlib only."""
import gzip
import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stdout, redirect_stderr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import units_from_annotation as U  # noqa: E402


def gff(genes):
    """genes: (name, type, strand, [exons 1-based closed], description) -> GFF text; each gene has one mRNA child."""
    out = ['##gff-version 3']
    for name, typ, strand, exons, desc in genes:
        lo, hi = min(a for a, _ in exons) if exons else 1, max(b for _, b in exons) if exons else 2
        d = f';description={desc}' if desc else ''
        out.append(f'c1\tRefSeq\t{typ}\t{lo}\t{hi}\t.\t{strand}\t.\tID=gene-{name};Name={name}{d}')
        if exons:
            out.append(f'c1\tRefSeq\tmRNA\t{lo}\t{hi}\t.\t{strand}\t.\tID=rna-{name};Parent=gene-{name};Name={name}')
            for k, (a, b) in enumerate(exons):
                out.append(f'c1\tRefSeq\texon\t{a}\t{b}\t.\t{strand}\t.\tID=exon-{name}-{k};Parent=rna-{name}')
    return '\n'.join(out) + '\n'


def gtf(txs):
    """txs: (tid, gene, strand, [exons]) -> GTF text."""
    out = []
    for tid, gene, strand, ex in txs:
        out.append(f'c1\trustle\ttranscript\t{ex[0][0]}\t{ex[-1][1]}\t.\t{strand}\t.\tgene_id "{gene}"; transcript_id "{tid}"; reads "5";')
        for k, (a, b) in enumerate(ex, 1):
            out.append(f'c1\trustle\texon\t{a}\t{b}\t.\t{strand}\t.\tgene_id "{gene}"; transcript_id "{tid}"; exon_number "{k}";')
    return '\n'.join(out) + '\n'


A = [(100, 200), (300, 400), (500, 600)]
B = [(2000, 2100), (2200, 2300)]
C = [(5000, 5100), (5200, 5300)]


def write(path, text):
    """a text file, gzipped when the name ends in `.gz`."""
    with (gzip.open(path, 'wt') if path.endswith('.gz') else open(path, 'w')) as fh:
        fh.write(text)


def run(genes, txs, **kw):
    with tempfile.TemporaryDirectory() as d:
        write(f'{d}/a.gff', gff(genes))
        write(f'{d}/t.gtf', gtf(txs))
        return U.run(f'{d}/a.gff', f'{d}/t.gtf', **kw)


class T(unittest.TestCase):
    def test_a_fusion_of_two_genes_is_cut_at_the_intron_between_them(self):
        rows, unc, n = run([('A', 'gene', '+', A, ''), ('B', 'gene', '+', B, '')],
                           [('F', 'g1', '+', A + B), ('A1', 'g1', '+', A), ('S', 'g1', '+', [(100, 200)])])
        self.assertEqual(rows, [('F', 'c1', '+', '601-1999', 'ann:A|B', 'g1')])
        self.assertEqual((unc, n['transcripts'], n['cuts'], n['spliced']), ([], 1, 1, 2), 'the single-gene and the unspliced transcripts are not listed')

    def test_the_widest_intron_between_the_blocks_ties_to_the_first(self):
        up = [(100, 200), (300, 400)]
        down = [(1500, 1600), (1700, 1800)]
        mid = [(500, 600), (800, 900), (1200, 1300)]      # exons over no gene: introns 401-499, 601-799, 901-1199, 1301-1499
        rows, _, _ = run([('U', 'gene', '+', up, ''), ('D', 'gene', '+', down, '')], [('F', 'g', '+', up + mid + down)])
        self.assertEqual(rows[0][3], '901-1199')
        # two introns of 199 bp between the blocks: the first, in genomic order
        rows, _, _ = run([('U', 'gene', '+', [(100, 200)], ''), ('D', 'gene', '+', [(700, 800)], '')],
                         [('F', 'g', '+', [(100, 200), (400, 500), (700, 800)])])
        self.assertEqual(rows[0][3], '201-399')

    def test_three_genes_two_cuts_in_genomic_order_on_either_strand(self):
        for strand in '+-':
            rows, _, _ = run([('A', 'gene', strand, A, ''), ('B', 'gene', strand, B, ''), ('C', 'pseudogene', strand, C, '')],
                             [('F', 'g', strand, A + B + C)])
            self.assertEqual(rows, [('F', 'c1', strand, '601-1999;2301-4999', 'ann:A|B|C', 'g')])

    def test_genes_of_the_other_strand_are_ignored(self):
        rows, _, n = run([('A', 'gene', '+', A, ''), ('B', 'gene', '-', B, '')], [('F', 'g', '+', A + B)])
        self.assertEqual((rows, n['ge2_genes']), ([], 0))

    def test_shared_exons_and_nested_genes_are_uncuttable_and_reported(self):
        # B2 is a second model over B's exons: one block, never cut
        rows, unc, n = run([('A', 'gene', '+', A, ''), ('B', 'gene', '+', B, ''), ('B2', 'pseudogene', '+', B, '')],
                           [('F', 'g', '+', A + B)])
        self.assertEqual(len(rows), 1, 'A | {B, B2}: still two blocks')
        self.assertEqual(rows[0][4], 'ann:A|B2|B', 'genes of one block sort by their id `name|start`, as the oracle did')
        # two genes over the SAME exons: no cut
        rows, unc, n = run([('A', 'gene', '+', A, ''), ('A2', 'pseudogene', '+', A, '')], [('F', 'g', '+', A)])
        self.assertEqual((rows, unc, n['uncuttable_shared_exon']), ([], [('F', 'c1', '+', 'A2|A', 'shared_exon')], 1))
        # nested: X overlaps exons 0 and 3, Y exons 1 and 2, no exon is in both
        t = [(100, 200), (300, 400), (500, 600), (700, 800)]
        rows, unc, n = run([('X', 'gene', '+', [t[0], t[3]], ''), ('Y', 'gene', '+', [t[1], t[2]], '')], [('F', 'g', '+', t)])
        self.assertEqual((rows, unc, n['uncuttable_nested_genes']), ([], [('F', 'c1', '+', 'X|Y', 'nested_genes')], 1))

    def test_readthrough_records_and_exonless_genes_are_not_genes(self):
        # RT is RefSeq's own fused model over A and B: excluded, so a transcript over A and B is still a fusion of A and B
        rt = A + B
        rows, _, _ = run([('A', 'gene', '+', A, ''), ('B', 'gene', '+', B, ''), ('A-B', 'gene', '+', rt, 'A-B readthrough%2C x')],
                         [('F', 'g', '+', A + B)])
        self.assertEqual((len(rows), rows[0][4]), (1, 'ann:A|B'))
        # a transcript over RT and A only: with RT excluded it touches one gene
        rows, _, n = run([('A', 'gene', '+', A, ''), ('A-B', 'gene', '+', rt, 'ReadThrough')], [('F', 'g', '+', A + [(700, 800)])])
        self.assertEqual((rows, n['ge2_genes']), ([], 0))
        # a gene record with no exon is dropped
        rows, _, n = run([('A', 'gene', '+', A, ''), ('E', 'gene', '+', [], '')], [('F', 'g', '+', A + B)])
        self.assertEqual((rows, n['ge2_genes']), ([], 0))

    def test_replacement_records_replace_a_same_named_record_or_are_added(self):
        with tempfile.TemporaryDirectory() as d:
            write(f'{d}/a.gff', gff([('A', 'gene', '+', A, ''), ('B', 'gene', '+', [(9000, 9100)], '')]))
            write(f'{d}/t.gtf', gtf([('F', 'g', '+', A + B)]))
            write(f'{d}/r.tsv', 'name\tchrom\tstrand\texons\nB\tc1\t+\t2000-2100,2200-2300\nC\tc1\t+\t5000-5100\n')
            self.assertEqual(U.run(f'{d}/a.gff', f'{d}/t.gtf')[0], [], 'B is elsewhere in the annotation')
            rows, _, _ = U.run(f'{d}/a.gff', f'{d}/t.gtf', replace=U.read_replacements(f'{d}/r.tsv'), label_prefix='H:')
            self.assertEqual(rows, [('F', 'c1', '+', '601-1999', 'H:A|B', 'g')])

    def test_cli_writes_the_list_and_the_report(self):
        with tempfile.TemporaryDirectory() as d:
            write(f'{d}/a.gff', gff([('A', 'gene', '+', A, ''), ('B', 'gene', '+', B, ''), ('A2', 'gene', '+', A, '')]))
            write(f'{d}/t.gtf', gtf([('F', 'g', '+', A + B), ('G', 'g', '+', A)]))
            out, err = io.StringIO(), io.StringIO()
            with redirect_stdout(out), redirect_stderr(err):
                U.main([f'{d}/a.gff', f'{d}/t.gtf', '--report', f'{d}/unc.tsv'])
            self.assertEqual(out.getvalue().splitlines()[0], 'tid\tchrom\tstrand\tjunctions\tlabel\tgene')
            self.assertEqual(out.getvalue().splitlines()[1], 'F\tc1\t+\t601-1999\tann:A2|A|B\tg')
            with open(f'{d}/unc.tsv') as fh:
                self.assertEqual(fh.read().splitlines()[1:], ['G\tc1\t+\tA2|A\tshared_exon'])
            self.assertIn('1 cut at 1 introns', err.getvalue())
            out = io.StringIO()
            with redirect_stdout(out), redirect_stderr(io.StringIO()):
                U.main([f'{d}/a.gff', f'{d}/t.gtf', '--contigs', 'c2'])
            self.assertEqual(out.getvalue().strip(), 'tid\tchrom\tstrand\tjunctions\tlabel\tgene', '--contigs filters both files')


    def test_gzipped_annotation_and_gtf_give_the_same_list(self):
        genes = [('A', 'gene', '+', A, ''), ('B', 'gene', '+', B, '')]
        txs = [('F', 'g1', '+', A + B), ('A1', 'g1', '+', A)]
        with tempfile.TemporaryDirectory() as d:
            for name in ('a.gff', 'a.gff.gz'):
                write(f'{d}/{name}', gff(genes))
            for name in ('t.gtf', 't.gtf.gz'):
                write(f'{d}/{name}', gtf(txs))
            plain = U.run(f'{d}/a.gff', f'{d}/t.gtf')
            self.assertEqual(plain[0], [('F', 'c1', '+', '601-1999', 'ann:A|B', 'g1')])
            for g, t in (('a.gff.gz', 't.gtf'), ('a.gff', 't.gtf.gz'), ('a.gff.gz', 't.gtf.gz')):
                self.assertEqual(U.run(f'{d}/{g}', f'{d}/{t}'), plain, (g, t))
            with open(f'{d}/a.gff.gz', 'rb') as fh:
                self.assertEqual(fh.read(2), b'\x1f\x8b', 'the fixture really is gzip')

    def test_an_empty_result_is_a_header_only_list_and_a_note_not_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            write(f'{d}/a.gff', gff([('A', 'gene', '+', A, ''), ('B', 'gene', '+', B, '')]))
            write(f'{d}/t.gtf', gtf([('A1', 'g1', '+', A), ('B1', 'g2', '+', B)]))     # no fusion in the sample
            out, err = io.StringIO(), io.StringIO()
            with redirect_stdout(out), redirect_stderr(err):
                U.main([f'{d}/a.gff', f'{d}/t.gtf'])                                    # exits normally: no SystemExit
            self.assertEqual(out.getvalue(), 'tid\tchrom\tstrand\tjunctions\tlabel\tgene\n')
            self.assertIn('0 cut at 0 introns', err.getvalue())
            self.assertIn('HEADER-ONLY', err.getvalue())
            self.assertIn('run without it', err.getvalue())
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                U.main([f'{d}/a.gff', f'{d}/t.gtf', '--out', f'{d}/list.tsv'])
            with open(f'{d}/list.tsv') as fh:
                self.assertEqual(fh.read(), 'tid\tchrom\tstrand\tjunctions\tlabel\tgene\n')

if __name__ == '__main__':
    unittest.main()
