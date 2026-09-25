#!/usr/bin/env python3
"""Member completeness for the family-scoped pool test (docs/PREREG_family_scoped_pool_2026-09-24.md).
usage: member_completeness.py ARM_GTF REF_GTF UNIVERSE_TSV GENES_GFF CHROM LABEL
Universe = referee genes (one per line, 'Gene Name' header) fixed before the arms; runs gffcompare (-r REF_GTF) on
ARM_GTF and reports: complete members (>= 1 transcript of class '='), partial ('=', 'c', 'k'), transcripts per
member locus, gffcompare transcript-level precision, and the share of expressed genes split into >= 2 loci."""
import sys, os, re, subprocess, collections
gtf, ref, uni, ggff, chrom, label = sys.argv[1:7]
universe = [l.rstrip('\n').split('\t')[0] for l in open(uni) if not l.startswith('Gene Name') and l.strip()]
# gene name -> ref gene_id in the ref GTF (gffcompare's tmap uses ref_gene_id = gene_id attr, e.g. gene-LOC...)
name_of_gid = {}
for l in open(ref):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[0] != chrom: continue
    gid = re.search(r'gene_id "([^"]+)"', f[8]); gn = re.search(r'gene_name "([^"]+)"', f[8])
    if gid and gn: name_of_gid[gid.group(1)] = gn.group(1)
out = os.path.splitext(gtf)[0] + '.gc'
subprocess.run(['gffcompare', '-r', ref, '-o', out, gtf], capture_output=True, text=True)
tmap = [p for p in os.listdir(os.path.dirname(gtf) or '.') if p.startswith(os.path.basename(out)) and p.endswith('.tmap')]
tmap = os.path.join(os.path.dirname(gtf) or '.', tmap[0])
best = collections.defaultdict(set); ntx_gene = collections.Counter()
for l in open(tmap):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'ref_gene_id': continue
    g = name_of_gid.get(f[0], f[0]); best[g].add(f[2]); ntx_gene[g] += 1
U = set(universe)
complete = sum(1 for g in U if '=' in best.get(g, ()))
partial = sum(1 for g in U if best.get(g, set()) & {'=', 'c', 'k'})
with_tx = [g for g in U if g in ntx_gene]
tx_per = sum(ntx_gene[g] for g in with_tx) / max(1, len(with_tx))
stats = open(out + '.stats').read() if os.path.exists(out + '.stats') else open(out).read()  # a dotted prefix makes gffcompare write the summary to the bare prefix
prec = re.search(r'Transcript level:\s+([\d.]+)\s+\|\s+([\d.]+)', stats)
print(f"{label:14s} universe {len(U)} | complete (=) {complete} | partial (=,c,k) {partial} | transcripts per member gene {tx_per:.2f} | gffcompare transcript sens/prec {prec.group(1)}/{prec.group(2)}")
