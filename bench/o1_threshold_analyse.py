#!/usr/bin/env python3
"""G3 analysis: anchor cohesion on gorilla + Soto scores per grid point.  usage: o1_threshold_analyse.py <g3_dir>"""
import sys, os, csv, glob, collections, subprocess
G = sys.argv[1]; M = '/mnt/linuxdisk/home/juanfraitu/mcl_ann'; S = '/mnt/linuxdisk/home/juanfraitu/soto_mcl'
R = '/mnt/c/Users/jfris/Desktop/Rustle'
def load(p): return [r for r in csv.DictReader(open(p), delimiter='\t')]
bp1 = load(f'{M}/rna_bp1_p9.clusters.tsv')
anchors = {'NPIP': 'MCL2', 'MCL1': 'MCL1', 'MCL3': 'MCL3', 'MCL4': 'MCL4', 'LCR16u': 'MCL7', 'L1blob': 'MCL0'}
A = {k: [r for r in bp1 if r['cluster_id'] == v] for k, v in anchors.items()}
def ov(a, b): return a['chrom'] == b['chrom'] and int(a['start']) < int(b['end']) and int(b['start']) < int(a['end'])
def modal(mem, rows):
    c = collections.Counter()
    for m in mem:
        hit = sorted({r['cluster_id'] for r in rows if ov(m, r)})
        c[hit[0] if hit else '-'] += 1
    top = c.most_common(1)[0]
    return top[0], top[1] / len(mem)
tags = sorted({os.path.basename(p)[4:-13] for p in glob.glob(f'{G}/ggo_*.clusters.tsv')})
print('tag\tNPIP\tMCL1\tMCL3\tMCL4\tNPIP||LCR16u\tL1blob\tSoto_det\tSoto_bandP\tSoto_Rboth\tSoto_fam')
for tag in tags:
    rows = load(f'{G}/ggo_{tag}.clusters.tsv')
    coh = {k: modal(A[k], rows) for k in anchors}
    sep = 'yes' if coh['NPIP'][0] != coh['LCR16u'][0] else 'NO'
    hs = f'{G}/hsa_{tag}'
    subprocess.run(['python3', f'{R}/bench/mcl_to_cat_copies.py', hs, f'{S}/allgenes.asm20.paf', f'{G}/arm_{tag}_u', f'{G}/arm_{tag}'], capture_output=True)
    subprocess.run(['python3', f'{R}/bench/mcl_edge_dump.py', hs, f'{S}/allgenes.asm20.paf', f'{G}/arm_{tag}'], capture_output=True)
    out = subprocess.run(['python3', f'{R}/bench/soto_adjudicate.py', f'{R}/bench/soto/80_fams.chr.bed', f'{G}/arm_{tag}'], capture_output=True, text=True).stdout
    sec = out.split('floor 50%')[-1]
    def grab(key):
        for l in sec.splitlines():
            if key in l: return l
        return ''
    det = grab('DETECTION').split('=')[1].split()[0] if grab('DETECTION') else 'NA'
    band = grab('[0.90,1.00)').split()[-1] if grab('[0.90,1.00)') else 'NA'
    rb = grab('both detected').split('=')[1].split()[0] if grab('both detected') else 'NA'
    fam = grab('FAMILY').split(':')[1].split()[0] if grab('FAMILY') else 'NA'
    print(f"{tag}\t{coh['NPIP'][1]:.2f}\t{coh['MCL1'][1]:.2f}\t{coh['MCL3'][1]:.2f}\t{coh['MCL4'][1]:.2f}\t{sep}\t{coh['L1blob'][1]:.2f}\t{det}\t{band}\t{rb}\t{fam}")
