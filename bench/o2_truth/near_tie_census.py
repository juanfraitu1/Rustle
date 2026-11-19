#!/usr/bin/env python3
"""PREREG r1829 (docs/PREREG_near_tie_census_2026-11-19.md), layer 1: near-tie tiers per sample from a BAM's as_table.

  near_tie_census.py SAMPLE AS_TABLE [--sd SEDEF.bed --exclude EXCL.bed --genes CAT.tsv] [--family-names NAMES.txt]
                     [--hash-pct 25 --hash-offset 0] > SAMPLE.census.tsv

as_table columns: name best_as second_as n_records primary_chrom primary_as best_ms second_ms primary_ms primary_mapq
primary_spliced rival_as rival_ms rival_spliced rival_chrom primary_pos rival_pos (rival = best-AS non-primary record).
Tiers (a rival at another place): exact = rival ms == primary ms and rival AS >= primary AS; ms10 = |primary ms - rival ms|
<= 10; good = rival AS >= 0.98 x best AS. Same place = rival on the primary's contig within 5,000 of its position.
Location (with --sd / --exclude): excluded / family / SD / SD one side / other, first that applies. `family` = the read's
primary position inside a truth-family gene span (--genes; famtruth; an approximation: span, not exons). --family-names
(HG002: neutral_truth's selection of residues 0-24) and --hash-* count that family subset apart (calibration K1)."""
import argparse, bisect, collections, hashlib, os, sys

def intervals(path, pairs=False):
    iv = collections.defaultdict(list)
    for l in open(path):
        if l[0] == '#' or not l.strip():
            continue
        f = l.split('\t')
        iv[f[0]].append((int(f[1]), int(f[2])))
        if pairs:                                           # SEDEF: both copies of a pair
            iv[f[3]].append((int(f[4]), int(f[5])))
    out = {}
    for c, v in iv.items():
        v.sort(); m = []
        for s, e in v:
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        out[c] = ([x[0] for x in m], [x[1] for x in m])
    return out

def inside(iv, c, p):
    if c not in iv:
        return False
    s, e = iv[c]
    i = bisect.bisect_right(s, p) - 1
    return i >= 0 and p < e[i]

def selected(name, pct, offset):
    return (int(hashlib.sha256(name.encode()).hexdigest()[:8], 16) - offset) % 100 < pct

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('sample'); ap.add_argument('table')
    ap.add_argument('--sd'); ap.add_argument('--exclude'); ap.add_argument('--genes'); ap.add_argument('--family-names')
    ap.add_argument('--hash-pct', type=int); ap.add_argument('--hash-offset', type=int, default=0)
    a = ap.parse_args()
    sd = intervals(a.sd, pairs=True) if a.sd else None
    ex = intervals(a.exclude) if a.exclude else None
    fam_iv = None
    if a.genes:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
        import famtruth as F
        G = F.load_genes('human')
        pairs, _ = F.truth_pairs(G, F.PARALOGS, set(F.LEVELS['simiiformes']))
        tmp = collections.defaultdict(list)
        for f in F.components(pairs):
            for g in f:
                v = G[g]
                if v['exons']:
                    tmp[v['chrom']].append((v['exons'][0][0], v['exons'][-1][1]))
        fam_iv = {}
        for c, v in tmp.items():
            v.sort(); m = []
            for s, e in v:
                if m and s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
                else: m.append([s, e])
            fam_iv[c] = ([x[0] for x in m], [x[1] for x in m])
    fam_names = set(l.split()[0].lstrip('>') for l in open(a.family_names) if l[0] == '>') if a.family_names else None
    C = collections.Counter()
    L = collections.defaultdict(collections.Counter)     # tier -> location counts
    for l in open(a.table):
        if l[0] == '#':
            continue
        f = l.rstrip('\n').split('\t')
        if f[4] == 'NA' or f[4] == '*':
            C['no_primary'] += 1; continue
        C['molecules'] += 1
        name = f[0]
        fam = None
        if fam_names is not None:
            fam = name in fam_names
        if a.hash_pct is not None and fam and selected(name, a.hash_pct, a.hash_offset):
            C['hash_family'] += 1
        if f[11] == 'NA':
            C['one_place'] += 1; continue
        pc, rc = f[4], f[14]
        pp, rp = int(f[15]), int(f[16])
        if rc == pc and abs(rp - pp) < 5000:
            C['same_place'] += 1; continue
        C['two_places'] += 1
        pms, rms, pas, ras, bas = int(f[8]), int(f[12]), int(f[5]), int(f[11]), int(f[1])
        tiers = []
        if rms == pms and ras >= pas: tiers.append('exact')
        if abs(pms - rms) <= 10: tiers.append('ms10')
        if ras >= 0.98 * bas: tiers.append('good')
        if 'good' in tiers and 'ms10' not in tiers: tiers.append('good_not_ms10')
        struct = f[10] != f[13]
        for t in tiers:
            C[t] += 1
            if struct: C[t + '_structure_differs'] += 1
            if a.hash_pct is not None and fam and selected(name, a.hash_pct, a.hash_offset):
                C['hash_family_' + t] += 1
            if sd is not None:
                fam_loc = fam_iv is not None and inside(fam_iv, pc, pp)
                if ex and (inside(ex, pc, pp) or inside(ex, rc, rp)): loc = 'excluded'
                elif fam_loc: loc = 'family'
                else:
                    n = inside(sd, pc, pp) + inside(sd, rc, rp)
                    loc = 'SD' if n == 2 else 'SD_one_side' if n == 1 else 'other'
                L[t][loc] += 1
                if inside(sd, pc, pp) and inside(sd, rc, rp): L[t]['both_in_SD_any_location'] += 1
            if t == 'ms10' and pc != rc: C['ms10_other_contig'] += 1
    w = sys.stdout
    w.write('sample\tkey\tcount\tshare_of_molecules\n')
    for k in sorted(C):
        w.write(f"{a.sample}\t{k}\t{C[k]}\t{C[k] / max(1, C['molecules']):.5f}\n")
    for t in sorted(L):
        for k in sorted(L[t]):
            w.write(f"{a.sample}\tloc:{t}:{k}\t{L[t][k]}\t{L[t][k] / max(1, C[t]):.5f}\n")

if __name__ == '__main__':
    main()
