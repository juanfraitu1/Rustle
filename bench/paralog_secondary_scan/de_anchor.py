#!/usr/bin/env python3
"""de-based per-read copy anchoring — identical to nm_anchor.py but the
discriminant is the gap-compressed event count (de * aligned_len) instead of raw
NM. minimap2's `de` counts each indel as ONE event (not per-base) and is the
fairer 'how well does this read fit this copy' measure in indel-rich HiFi
alignments. Falls back to raw NM when `de` is absent.

Usage: de_anchor.py MM_PLACEMENTS_DE.tsv GENE_INTRONS.tsv ENRICHED.json OUT.json
  (MM_PLACEMENTS_DE.tsv columns: qname flag chrom pos cigar nm de AS ts)
Env: DE_T(2)  MIN_ANCHOR(3)  ANCHOR_K(0=adaptive)  EXTENT_FRAC(0.7)
"""
import sys, os, re, json, bisect, collections

DE_T        = float(os.environ.get('DE_T', '2'))
MIN_ANCHOR  = int(os.environ.get('MIN_ANCHOR', '3'))
K_FIXED     = int(os.environ.get('ANCHOR_K', '0'))
EXTENT_FRAC = float(os.environ.get('EXTENT_FRAC', '0.7'))

def parse_cig(cig): return re.findall(r'(\d+)([MIDNSH=X])', cig)
def junctions(pos, cig):
    cur = pos; js = []
    for ln, op in parse_cig(cig):
        ln = int(ln)
        if op == 'N': js.append((cur - 1, cur + ln)); cur += ln
        elif op in 'MD=X': cur += ln
    return js
def alen(cig): return sum(int(l) for l, o in parse_cig(cig) if o in 'M=XD') or 1

def main():
    placements_f, gi, enriched_f, out = sys.argv[1:5]
    genes = {}; by_chrom = collections.defaultdict(list)
    for line in open(gi):
        if line.startswith('gene_id'): continue
        gid, chrom, s, e, strand, ntx, chain = line.rstrip('\n').split('\t')
        ints = set()
        if chain:
            for pr in chain.split(','):
                d, a = pr.split(':'); ints.add((int(d), int(a)))
        genes[gid] = dict(chrom=chrom, start=int(s), end=int(e), strand=strand, introns=ints)
        by_chrom[chrom].append((int(s), int(e), gid))
    for c in by_chrom: by_chrom[c].sort()

    def gene_at(chrom, pos):
        arr = by_chrom.get(chrom)
        if not arr: return None
        i = bisect.bisect_right(arr, (pos, float('inf'), '')) - 1
        for j in range(i, max(-1, i - 6), -1):
            if j < 0: break
            ss, ee, gid = arr[j]
            if ss <= pos <= ee: return gid
        return None

    reads = collections.defaultdict(list); pos_index = collections.defaultdict(list)
    for line in open(placements_f):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6 or f[4] == '*' or f[5] == '': continue
        qn, flag, chrom, pos, cig, nm = f[:6]
        de = f[6] if len(f) > 6 else ''
        p = int(pos); al = alen(cig)
        metric = float(de) * al if de else float(nm)   # gap-compressed event count
        reads[qn].append((chrom, p, cig, metric))
        pos_index[chrom].append((p, qn))
    for c in pos_index: pos_index[c].sort()
    print(f"loaded placements for {len(reads)} reads", file=sys.stderr)

    def reads_in_locus(chrom, start, end):
        arr = pos_index.get(chrom)
        if not arr: return []
        lo = bisect.bisect_left(arr, (start - 10, '')); hi = bisect.bisect_right(arr, (end + 10, chr(0x10ffff)))
        return set(qn for _, qn in arr[lo:hi])

    hits = json.load(open(enriched_f))
    def kmin(g):
        n = len(g['introns'])
        return K_FIXED if K_FIXED > 0 else max(2, min(8, (n + 2) // 3))

    out_rows = []; tax = collections.Counter()
    for h in hits:
        gid = h['gene']; g = genes[gid]; K = kmin(g)
        ga = sa = amb = 0; sib_anchor = collections.Counter()
        for qn in reads_in_locus(g['chrom'], g['start'], g['end']):
            placs = reads[qn]
            g_m = None; g_al = None
            for (chrom, pos, cig, m) in placs:
                if chrom != g['chrom'] or not (g['start'] - 10 <= pos <= g['end']): continue
                if sum(1 for j in junctions(pos, cig) if j in g['introns']) >= K:
                    if g_m is None or m < g_m: g_m = m; g_al = alen(cig)
            if g_m is None: continue
            best_other = None; best_other_gene = None
            for (chrom, pos, cig, m) in placs:
                if chrom == g['chrom'] and g['start'] - 10 <= pos <= g['end']: continue
                if alen(cig) < EXTENT_FRAC * g_al: continue
                if best_other is None or m < best_other:
                    best_other = m; best_other_gene = gene_at(chrom, pos)
            if best_other is None: ga += 1; continue
            dM = best_other - g_m                  # gap-compressed event difference
            if dM >= DE_T: ga += 1
            elif dM <= -DE_T: sa += 1; sib_anchor[best_other_gene or 'intergenic'] += 1
            else: amb += 1
        total = ga + sa + amb
        if ga >= MIN_ANCHOR: t = 'expressed_real_copy'
        elif ga == 0 and sa > 0: t = 'pure_spillover'
        elif amb >= max(sa, ga) and amb >= 3: t = 'ambiguous_dominated'
        else: t = 'weak'
        tax[(h['klass'], t)] += 1
        dom = sib_anchor.most_common(1)[0] if sib_anchor else (None, 0)
        o = dict(h); o.update(de_anchor_G=ga, de_anchor_sibling=sa, de_anchor_ambiguous=amb,
                              de_anchor_total=total, de_taxonomy=t, de_dom_sibling=dom[0])
        out_rows.append(o)
    json.dump(out_rows, open(out, 'w'), indent=1)
    print("de taxonomy (klass, de_taxonomy):", file=sys.stderr)
    for k, v in sorted(tax.items()): print(f"   {k}: {v}", file=sys.stderr)

if __name__ == '__main__':
    main()
