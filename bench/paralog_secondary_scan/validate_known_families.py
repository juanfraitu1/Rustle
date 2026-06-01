#!/usr/bin/env python3
"""Validate the multi-copy-family definition against curated KNOWN families.

For each curated family (members assembled by GFF description keyword), evaluate
the M/H/X/I criteria against the read data and report whether the definition
holds:
  M (multiplicity)   = number of annotated members (>=2 required)
  H (homology/share) = reads that chain-match >=2 members (multi-map within family)
  X (expression)     = members with >=k reads that decisively own them: among the
                       family members a read chain-matches, best-member NM beats
                       2nd-best member by raw dNM>=T
  I (identifiability) = of reads chain-matching >=1 member, fraction that are
                       within-family TIES (no distinguishing position among members)
Verdict: 'family' if M>=2 and X>=2 (>=2 expressed members) and H>0; else note why.
Identifiability label: full / partial / none by tie fraction.

Usage: validate_known_families.py GENE_INTRONS.tsv GENE_META.tsv MM_PLACEMENTS.tsv
"""
import sys, re, json, bisect, collections, statistics as st

T = 2          # distinguishing-mismatch threshold
KEXPR = 3      # reads to call a member expressed

FAMILIES = {
    'DAZ (azoospermia, Y)':            'deleted in azoospermia',
    'NBPF (neuroblastoma bp)':         'neuroblastoma breakpoint',
    'TSPY (testis-specific Y)':        'testis-specific',
    'RBMY (RNA-binding Y)':            'RNA-binding motif protein, Y',
    'CDY (chromodomain Y)':            'chromodomain protein, Y',
    'GOLGA (golgin A)':                'golgin subfamily A',
    'MAGEA (melanoma antigen)':        'melanoma-associated antigen',
    'PRAME-like':                      'PRAME family',
    'tubulin beta':                    'tubulin beta',
    'beta-defensin':                   'beta-defensin',
    'histone H2/H3/H4':                'histone',
    'protocadherin':                   'protocadherin',
    'zinc finger (sample)':            'zinc finger protein',
    'amylase':                         'amylase',
}

def parse_cig(c): return re.findall(r'(\d+)([MIDNSH=X])', c)
def junctions(pos, cig):
    cur = pos; js = []
    for ln, op in parse_cig(cig):
        ln = int(ln)
        if op == 'N': js.append((cur - 1, cur + ln)); cur += ln
        elif op in 'MD=X': cur += ln
    return js
def alen(c): return sum(int(l) for l, o in parse_cig(c) if o in 'M=XD') or 1

def main():
    gi, meta, mm = sys.argv[1:4]
    genes = {}
    for line in open(gi):
        if line.startswith('gene_id'): continue
        gid, chrom, s, e, strand, ntx, chain = line.rstrip('\n').split('\t')
        ints = set()
        if chain:
            for pr in chain.split(','):
                d, a = pr.split(':'); ints.add((int(d), int(a)))
        genes[gid] = dict(chrom=chrom, start=int(s), end=int(e), strand=strand,
                          introns=ints, single_exon=(not chain))
    desc = {}
    for line in open(meta):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 4: desc[f[0]] = (f[2], f[3])    # name, description
    # placements
    reads = collections.defaultdict(list); idx = collections.defaultdict(list)
    for line in open(mm):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6 or f[5] == '' or f[4] == '*': continue
        qn, flag, chrom, pos, cig, nm = f
        p = int(pos); reads[qn].append((chrom, p, cig, int(nm))); idx[chrom].append((p, qn))
    for c in idx: idx[c].sort()
    def reads_in(chrom, a, b):
        arr = idx.get(chrom, []); lo = bisect.bisect_left(arr, (a - 10, '')); hi = bisect.bisect_right(arr, (b + 10, chr(0x10ffff)))
        return set(q for _, q in arr[lo:hi])

    print(f"{'family':28s} {'M':>3} {'mexon':>5} {'expr':>4} {'H_share':>7} {'tieFrac':>7}  verdict")
    print('-' * 92)
    for fam, kw in FAMILIES.items():
        kwl = kw.lower()
        members = [g for g, (nm, d) in desc.items() if kwl in (d or '').lower() and g in genes]
        if not members:
            print(f"{fam:28s} {'--':>3}  (no annotated members matched)"); continue
        M = len(members)
        mexon = sum(1 for g in members if not genes[g]['single_exon'])
        # gather reads chain-matching each member
        member_reads = collections.defaultdict(dict)   # gid -> {qn: (nm, alen)}
        for g in members:
            G = genes[g]
            if G['single_exon']:
                K = 0  # single-exon: count overlap only (no chain) — flag separately
            K = max(2, min(8, (len(G['introns']) + 2) // 3)) if G['introns'] else None
            for qn in reads_in(G['chrom'], G['start'], G['end']):
                best = None
                for (chrom, pos, cig, nm) in reads[qn]:
                    if chrom != G['chrom'] or not (G['start'] - 10 <= pos <= G['end']): continue
                    if K is not None and sum(1 for j in junctions(pos, cig) if j in G['introns']) < K: continue
                    if best is None or nm < best[0]: best = (nm, alen(cig))
                if best is not None:
                    member_reads[g][qn] = best
        # per read: assign to best member vs 2nd best member (within-family identifiability)
        read_members = collections.defaultdict(list)   # qn -> [(gid, nm, alen)]
        for g, rs in member_reads.items():
            for qn, (nm, al) in rs.items():
                read_members[qn].append((g, nm, al))
        assigned = collections.Counter(); ties = 0; total = 0; shared = 0
        for qn, lst in read_members.items():
            total += 1
            if len(lst) >= 2: shared += 1
            lst.sort(key=lambda x: x[1])      # by NM
            if len(lst) == 1:
                assigned[lst[0][0]] += 1       # only fits one member -> owns it
            else:
                dnm = lst[1][1] - lst[0][1]
                if dnm >= T: assigned[lst[0][0]] += 1
                else: ties += 1
        expr = sum(1 for g in members if assigned[g] >= KEXPR)
        tief = ties / total if total else 0.0
        Hsh = shared / total if total else 0.0
        connected = Hsh >= 0.10
        if M < 2: verdict = 'NOT family (M<2 annotated)'
        elif total < KEXPR: verdict = 'NOT expressed in data — correctly excluded'
        elif not connected: verdict = 'NOT one VG family: members do not share backbone (H~0)'
        elif expr >= 2: verdict = 'FAMILY ✓ (>=2 resolvable copies)'
        elif total >= 10 and tief >= 0.6: verdict = 'FAMILY, NON-IDENTIFIABLE (expressed; copies indistinguishable)'
        elif expr == 1: verdict = 'gene + unexpressed/non-identifiable paralog'
        else: verdict = 'weak / ambiguous'
        idlabel = 'full' if tief < 0.15 else ('partial' if tief < 0.6 else 'NONE')
        print(f"{fam:28s} {M:>3} {mexon:>5} {expr:>4} {Hsh:>7.2f} {tief:>7.2f}  {verdict}  [I:{idlabel}]")

if __name__ == '__main__':
    main()
