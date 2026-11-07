"""famtruth.load_genes on any species: GTF, Ensembl GFF3 and RefSeq-style GFF3 annotations, own ids, taxon-list levels."""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F  # noqa: E402

GTF = '''c1\tx\texon\t101\t200\t.\t+\t.\tgene_id "ENSMUSG00000000001.4"; transcript_id "t1"; gene_name "Abc1"; gene_biotype "protein_coding";
c1\tx\texon\t151\t300\t.\t+\t.\tgene_id "ENSMUSG00000000001.4"; transcript_id "t2"; gene_name "Abc1"; gene_biotype "protein_coding";
c1\tx\texon\t1001\t1100\t.\t-\t.\tgene_id "ENSMUSG00000000002.1"; transcript_id "t3"; gene_type "lncRNA";
'''
ENS = '''##gff-version 3
c1\tx\tgene\t101\t300\t.\t+\t.\tID=gene:ENSMUSG00000000001;Name=Abc1;biotype=protein_coding
c1\tx\tmRNA\t101\t300\t.\t+\t.\tID=transcript:T1;Parent=gene:ENSMUSG00000000001
c1\tx\texon\t101\t300\t.\t+\t.\tParent=transcript:T1
c1\tx\tgene\t501\t600\t.\t+\t.\tID=gene:ENSMUSG00000000003;biotype=protein_coding
c1\tx\tmRNA\t501\t600\t.\t+\t.\tID=transcript:T3;Parent=gene:ENSMUSG00000000003
c1\tx\texon\t501\t600\t.\t+\t.\tParent=transcript:T3
'''
REFSEQ = '''##gff-version 3
c1\tx\tgene\t101\t300\t.\t+\t.\tID=gene-Abc1;gene=Abc1;gene_biotype=protein_coding
c1\tx\tmRNA\t101\t300\t.\t+\t.\tID=rna-XM_1;Parent=gene-Abc1
c1\tx\texon\t101\t200\t.\t+\t.\tParent=rna-XM_1
c1\tx\texon\t251\t300\t.\t+\t.\tParent=rna-XM_1
'''


def write(text, suffix):
    f = tempfile.NamedTemporaryFile('w', suffix=suffix, delete=False)
    f.write(text)
    f.close()
    return f.name


class LoadGenes(unittest.TestCase):
    def test_sniff_tells_the_three_formats_apart(self):
        self.assertEqual(F.sniff(write(GTF, '.gtf')), 'gtf')
        self.assertEqual(F.sniff(write(ENS, '.gff3')), 'ensembl')
        self.assertEqual(F.sniff(write(REFSEQ, '.gff3')), 'gff3')

    def test_gtf_exon_union_and_own_ids_drop_the_version(self):
        G = F.load_genes('mouse', write(GTF, '.gtf'), ids='own')
        d = G['ENSMUSG00000000001.4']
        self.assertEqual(d['exons'], [(100, 300)])
        self.assertEqual((d['name'], d['biotype']), ('Abc1', 'protein_coding'))
        self.assertEqual(d['ens'], {'ENSMUSG00000000001'})
        self.assertEqual(G['ENSMUSG00000000002.1']['biotype'], 'lncRNA')

    def test_gtf_without_biotype_is_coding_by_its_refseq_transcript(self):
        g = 'c1\tx\texon\t1\t50\t.\t+\t.\ttranscript_id "rna-XM_9.1"; gene_id "gene-Abc1"; gene_name "Abc1";\n'
        g += 'c1\tx\texon\t901\t950\t.\t+\t.\ttranscript_id "rna-NR_9.1"; gene_id "gene-Nc1"; gene_name "Nc1";\n'
        G = F.load_genes('mouse', write(g, '.gtf'), ids='own')
        self.assertEqual(G['gene-Abc1']['biotype'], 'protein_coding')
        self.assertEqual(G['gene-Nc1']['biotype'], 'other')

    def test_ensembl_names_from_the_name_attribute_without_xref(self):
        G = F.load_genes('mouse', write(ENS, '.gff3'), ids='own')
        self.assertEqual(G['gene:ENSMUSG00000000001']['name'], 'Abc1')
        self.assertEqual(G['gene:ENSMUSG00000000003']['name'], 'ENSMUSG00000000003')
        self.assertEqual(G['gene:ENSMUSG00000000003']['ens'], {'ENSMUSG00000000003'})
        self.assertNotIn('symbol', G['gene:ENSMUSG00000000001'])

    def test_refseq_gff3_own_ids_need_no_cat(self):
        saved = F.CAT_GFF, F.CAT_TSV
        F.CAT_GFF = F.CAT_TSV = None
        try:
            G = F.load_genes('mouse', write(REFSEQ, '.gff3'), ids='own')
        finally:
            F.CAT_GFF, F.CAT_TSV = saved
        self.assertEqual(G['gene-Abc1']['exons'], [(100, 200), (250, 300)])
        self.assertEqual(G['gene-Abc1']['ens'], {'Abc1'})

    def test_unnamed_is_loc_or_a_bare_ensembl_id_not_a_symbol(self):
        self.assertTrue(F.unnamed({'name': 'LOC101124002'}))
        self.assertTrue(F.unnamed({'name': 'ENSMUSG00000000003'}))
        self.assertFalse(F.unnamed({'name': 'ENSAP1'}))

    def test_levels_take_a_key_or_a_taxon_list(self):
        self.assertEqual(F.levels('recent')[0], 'Homo sapiens')
        self.assertEqual(F.levels('Mus musculus,Murinae'), ['Mus musculus', 'Murinae'])

    def test_missing_annotation_names_the_variable(self):
        with self.assertRaises(SystemExit) as e:
            F.load_genes('mouse')
        self.assertIn('--annotation', str(e.exception))


    def test_truth_pairs_drop_the_version_of_the_paralogue_table(self):
        G = {'a': dict(biotype='protein_coding', ens={'ENSG1'}), 'b': dict(biotype='protein_coding', ens={'ENSG2'})}
        with tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False) as fh:
            fh.write('Gene stable ID\tpar\ttype\tlca\nENSG1.3\tENSG2.7\tparalog\tHominidae\n')
        try:
            self.assertEqual(F.truth_pairs(G, fh.name, ['Hominidae']), ({('a', 'b')}, {'a'}))
        finally:
            os.unlink(fh.name)

if __name__ == '__main__':
    unittest.main()
