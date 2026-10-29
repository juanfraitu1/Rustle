#!/bin/bash
# PREREG U1-FRESH-HG002-INTRON scoring: intron filter -> MCL arms -> famtruth/pairjudge/protjudge/readouts. Run on a compute node or the frontend after job_hg002.sbatch.
set -euo pipefail
export PATH=$HOME/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin:$PATH
export DENOVO_CAT_GFF=$HOME/data/denovo_eval/chm13v2.0_gencode.genes.gff3 DENOVO_CAT_TSV=$HOME/work/jfim/rustle_data/annotation/chm13v2.0_CAT_Liftoff.genes.tsv
export DENOVO_PARALOGS=$HOME/data/denovo_eval/paralogs.tsv DENOVO_HGNC=$HOME/data/denovo_eval/hgnc_complete_set.txt
export DENOVO_HG002_GFF=$HOME/data/hg002/genes.gff3 DENOVO_HG002_XREF=$HOME/data/hg002/xref.tsv
RUN=$HOME/runs/hg002_kinnex; BIN=$HOME/rustle_target/release; B=$HOME/Rustle/bench/denovo_eval; FASTA=$HOME/data/hg002/hg002.pat.fa; T=${THREADS:-16}
cd $RUN
# 1. intron-homology evidence per U1 edge (species hg002 for the verdict column)
python3 $B/edge_dna.py hg.fam.loci.gff3 hg.fam.loci.gated.paf edgeh/hg/edges.u1.tsv --species hg002 > edgeh/edge_dna.tsv 2> edgeh/edge_dna.err
python3 - <<'PY'
import csv
rows = list(csv.DictReader(open('edgeh/edge_dna.tsv'), delimiter='\t'))
for name, T in (('intron', 1), ('intron500', 500)):
    keep = {(r['a'], r['b']) for r in rows if int(r['ii_cols']) >= T}
    n = 0
    with open(f'edgeh/edges.{name}.tsv', 'w') as w:
        for l in open('edgeh/hg/edges.u1.tsv'):
            a, b = l.split('\t')[:2]
            if (min(a, b), max(a, b)) in keep: w.write(l); n += 1
    print(name, 'edges kept', n, 'of', len(rows))
PY
# 2. MCL arms on the cached all-vs-all
for arm in u1:edgeh/hg/edges.u1.tsv intron:edgeh/edges.intron.tsv intron500:edgeh/edges.intron500.tsv; do
  n=${arm%%:*}; e=${arm#*:}
  RUSTLE_CACHE_DIR=$RUN/hg.cache $BIN/mcl_families --from-gtf hg.families.gtf --fasta $FASTA --threads $T --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --extra-edges $e --out hg$n.fam > hg$n.families.log 2>&1
  grep -h "extra-edges\|cluster(s)" hg$n.families.log | tail -2
done
# 3. truth
for p in hg hgu1 hgintron hgintron500; do
  python3 $B/famtruth.py score $p.fam.copies.tsv $p.famtruth.gtf.json --species hg002 --ggo-gff $DENOVO_HG002_GFF --xref $DENOVO_HG002_XREF --gtf hg.families.gtf --levels simiiformes --show all 2>&1 | tail -3
done
PAIRJUDGE_SPECIES=hg002 python3 $B/pairjudge.py hg.fam.copies.tsv hgu1.fam.copies.tsv hgintron.fam.copies.tsv hgintron500.fam.copies.tsv
# 4. protein mutual hit from HG002's own CDS
mkdir -p edgeh/prot
for n in u1 intron intron500; do python3 $B/protjudge.py pairs hg002 hg.fam.copies.tsv hg$n.fam.copies.tsv edgeh/prot/$n; done
cat edgeh/prot/*/genes.txt | sort -u > edgeh/prot/genes_all.txt
THREADS=$T python3 $B/protjudge.py prep hg002 $FASTA edgeh/prot/db edgeh/prot/genes_all.txt
python3 $B/protjudge.py rate edgeh/prot/db edgeh/prot/u1/base.pairs.tsv edgeh/prot/u1/added.pairs.tsv edgeh/prot/intron/added.pairs.tsv edgeh/prot/intron500/added.pairs.tsv
# 5. readouts
python3 $B/lift_misses.py hg edgeh/hg hgintron --species hg002 > edgeh/lift_misses.intron.tsv 2> edgeh/lift_misses.intron.err; tail -1 edgeh/lift_misses.intron.err
python3 $B/edge_verdicts.py hg.fam.loci.gff3 - edgeh/hg/edges.u1.tsv edgeh/edges.intron.tsv --species hg002 > edgeh/edge_verdicts.txt 2> edgeh/edge_verdicts.err || true
echo done
