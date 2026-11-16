"""tools/gtf_names.py: names are transferred by intron chain, else by exon overlap on the strand; lines stay otherwise equal."""
import os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
GFF = """c1\tx\tgene\t100\t900\t.\t+\t.\tID=gene-A;Name=GENEA
c1\tx\tmRNA\t100\t900\t.\t+\t.\tID=rna-A1;Parent=gene-A;transcript_id=NM_1
c1\tx\texon\t100\t200\t.\t+\t.\tParent=rna-A1
c1\tx\texon\t300\t400\t.\t+\t.\tParent=rna-A1
c1\tx\texon\t800\t900\t.\t+\t.\tParent=rna-A1
c1\tx\tgene\t850\t1500\t.\t-\t.\tID=gene-B;Name=GENEB
c1\tx\tmRNA\t850\t1500\t.\t-\t.\tID=rna-B1;Parent=gene-B
c1\tx\texon\t850\t1500\t.\t-\t.\tParent=rna-B1
"""
GTF_ANNOT = 'c1\tx\texon\t100\t200\t.\t+\t.\tgene_id "gA"; transcript_id "tA"; gene_name "GENEA";\n' \
            'c1\tx\texon\t300\t400\t.\t+\t.\tgene_id "gA"; transcript_id "tA"; gene_name "GENEA";\n'

def tx(t, strand, *exons):
    a = f'gene_id "{t}"; transcript_id "{t}";'
    return "".join(f"c1\tr\t{k}\t{s}\t{e}\t.\t{strand}\t.\t{a}\n" for k, s, e in
                   [("transcript", exons[0][0], exons[-1][1])] + [("exon", s, e) for s, e in exons])

def run(annot, gtf, name="a.gff3"):
    with tempfile.TemporaryDirectory() as d:
        a, g, o = os.path.join(d, name), os.path.join(d, "in.gtf"), os.path.join(d, "out.gtf")
        open(a, "w").write(annot); open(g, "w").write(gtf)
        subprocess.run([sys.executable, "-I", os.path.join(HERE, "gtf_names.py"), a, g, o], check=True, capture_output=True)
        return open(o).read().splitlines()

def test_names():
    src = (tx("exact", "+", (50, 200), (300, 400), (800, 950)) + tx("part", "+", (100, 200), (300, 380))
           + tx("mono_minus", "-", (860, 1000)) + tx("unstranded", ".", (1200, 1300)) + tx("far", "+", (5000, 5100)))
    out = run(GFF, src)
    by = {l.split('transcript_id "')[1].split('"')[0]: l for l in out if "\ttranscript\t" in l}
    assert 'ref_gene_name "GENEA"; ref_gene_id "A"; ref_transcript_id "A1"; ref_match "exact"; ref_genes "GENEA,GENEB' not in by["exact"]
    assert 'ref_gene_name "GENEA"; ref_gene_id "A"; ref_transcript_id "A1"; ref_match "exact";' in by["exact"]  # other strand ignored
    assert 'ref_gene_name "GENEA"' in by["part"] and 'ref_match "overlap"' in by["part"] and "ref_transcript_id" not in by["part"]
    assert 'ref_gene_name "GENEB"' in by["mono_minus"] and 'ref_gene_name "GENEB"' in by["unstranded"]
    assert "ref_" not in by["far"]
    strip = [l.split(" ref_gene_name")[0] for l in out]
    assert strip == src.splitlines()  # nothing else changes
    assert sum("ref_gene_name" in l for l in out) == 4 + 3 + 2 + 2  # exon lines carry the name too

def test_gtf_annotation():
    out = run(GTF_ANNOT, tx("m", "+", (100, 200), (300, 400)), "a.gtf")
    assert 'ref_gene_name "GENEA"; ref_gene_id "gA"; ref_transcript_id "tA"; ref_match "exact";' in out[0]
