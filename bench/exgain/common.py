"""Shared by the r1694 rebuild scripts (PREREG_o3_exon_gain_detectors_2026-11-13.md, Amendment 1)."""
import collections
C, T = 0.8, 0.90
def sexset(path):
    """SEX.tsv: genome (ref|mat|pat) <tab> contig. Replaces the workstation's chrX/chrY/chrM name test."""
    return {tuple(l.split()[:2]) for l in open(path) if l.strip() and not l.startswith('#')}
def load_hits(d, sex):
    """ex_sr_{ref,mat,pat}.paf -> exon -> [(target, start, end, strand)] at coverage >= C, identity >= T, sex contigs out."""
    H = {g: collections.defaultdict(list) for g in ('ref', 'mat', 'pat')}
    for g in H:
        for l in open(f'{d}/ex_sr_{g}.paf'):
            f = l.split('\t')
            if (int(f[3]) - int(f[2])) / int(f[1]) >= C and int(f[9]) / int(f[10]) >= T and (g, f[5]) not in sex:
                H[g][f[0]].append((f[5], int(f[7]), int(f[8]), f[4]))
    return H
def merged(xs):
    out = []
    for t, s, e, st in sorted(xs):
        if out and out[-1][0] == t and s < out[-1][2]: out[-1][2] = max(out[-1][2], e)
        else: out.append([t, s, e, st])
    return out
