#!/usr/bin/env python3
"""NM-based per-read copy anchoring (flag-independent).

For each candidate gene G and each multimapping read whose placement at G
matches G's annotated splice chain, compare the RAW edit distance (NM) at G vs
the read's best ALTERNATIVE placement of COMPARABLE alignment extent. The
discriminant is the number of distinguishing mismatches dNM = NM_other - NM_G
(each mismatch is a near-independent observation; the fingerprint-EM log-LR is
~ dNM * ln((1-eps)/eps), so even dNM>=2 is decisive). Per-base rate is the WRONG
discriminant: NM 14 vs 26 (12 distinguishing positions) is decisive, not a tie.

  G_anchored      : dNM >= T          (sequence says the read is G's)
  sibling_anchored: dNM <= -T         (spillover from a sibling copy)
  ambiguous       : |dNM| < T         (no distinguishing position -> the
                                       boundary-theorem non-identifiable case)

Extent guard: an alternative placement only competes if its aligned length is
>= EXTENT_FRAC of G's placement (a short partial hit elsewhere is not a real
competitor). A read with no comparable alternative is G_anchored (uniquely G's).

Sequence-grounded taxonomy per gene:
  expressed_real_copy : n_G_anchored >= MIN_ANCHOR        (recovery/quant win)
  pure_spillover      : n_G_anchored == 0                 (not expressed here;
                                                           assembling = fabrication)
  ambiguous_dominated : ambiguous is the plurality        (identifiability floor)
  weak                : otherwise

Usage: nm_anchor.py MM_PLACEMENTS.tsv GENE_INTRONS.tsv ENRICHED.json OUT.json
Env: DNM_T(2)  MIN_ANCHOR(3)  ANCHOR_K(0=adaptive)  EXTENT_FRAC(0.7)
"""
import sys, os, re, json, bisect, collections, statistics as st

DNM_T       = int(os.environ.get('DNM_T', '2'))
MIN_ANCHOR  = int(os.environ.get('MIN_ANCHOR', '3'))
K_FIXED     = int(os.environ.get('ANCHOR_K', '0'))
EXTENT_FRAC = float(os.environ.get('EXTENT_FRAC', '0.7'))

def parse_cig(cig):
    return re.findall(r'(\d+)([MIDNSH=X])', cig)

def junctions(pos, cig):
    cur = pos; js = []
    for ln, op in parse_cig(cig):
        ln = int(ln)
        if op == 'N':
            js.append((cur - 1, cur + ln)); cur += ln
        elif op in 'MD=X':
            cur += ln
    return js

def alen(cig):
    return sum(int(l) for l, o in parse_cig(cig) if o in 'M=XD') or 1

def main():
    placements_f, gi, enriched_f, out = sys.argv[1:5]
    genes = {}; by_chrom = collections.defaultdict(list)
    for line in open(gi):
        if line.startswith('gene_id'):
            continue
        gid, chrom, s, e, strand, ntx, chain = line.rstrip('\n').split('\t')
        ints = set()
        if chain:
            for pr in chain.split(','):
                d, a = pr.split(':'); ints.add((int(d), int(a)))
        genes[gid] = dict(chrom=chrom, start=int(s), end=int(e),
                          strand=strand, introns=ints)
        by_chrom[chrom].append((int(s), int(e), gid))
    for c in by_chrom:
        by_chrom[c].sort()

    def gene_at(chrom, pos):
        arr = by_chrom.get(chrom)
        if not arr:
            return None
        i = bisect.bisect_right(arr, (pos, float('inf'), '')) - 1
        for j in range(i, max(-1, i - 6), -1):
            if j < 0:
                break
            ss, ee, gid = arr[j]
            if ss <= pos <= ee:
                return gid
        return None

    # load placements per read + a locus index (chrom -> sorted [(pos,qname)])
    reads = collections.defaultdict(list)
    pos_index = collections.defaultdict(list)
    for line in open(placements_f):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6:
            continue
        qn, flag, chrom, pos, cig, nm = f
        if nm == '' or cig == '*':
            continue
        p = int(pos)
        reads[qn].append((chrom, p, cig, int(nm)))
        pos_index[chrom].append((p, qn))
    for c in pos_index:
        pos_index[c].sort()
    print(f"loaded placements for {len(reads)} reads", file=sys.stderr)

    def reads_in_locus(chrom, start, end):
        arr = pos_index.get(chrom)
        if not arr:
            return []
        lo = bisect.bisect_left(arr, (start - 10, ''))
        hi = bisect.bisect_right(arr, (end + 10, chr(0x10ffff)))
        return set(qn for _, qn in arr[lo:hi])

    hits = json.load(open(enriched_f))
    def kmin(g):
        n = len(g['introns'])
        return K_FIXED if K_FIXED > 0 else max(2, min(8, (n + 2) // 3))

    out_rows = []
    tax = collections.Counter()
    for h in hits:
        gid = h['gene']; g = genes[gid]; K = kmin(g)
        ga = sa = amb = 0
        sib_anchor = collections.Counter()
        margins = []
        for qn in reads_in_locus(g['chrom'], g['start'], g['end']):
            placs = reads[qn]
            # this read's chain-matching placement at G (lowest NM); keep raw NM + extent
            g_nm = None; g_alen = None
            for (chrom, pos, cig, nm) in placs:
                if chrom != g['chrom'] or not (g['start'] - 10 <= pos <= g['end']):
                    continue
                js = junctions(pos, cig)
                if sum(1 for j in js if j in g['introns']) >= K:
                    if g_nm is None or nm < g_nm:
                        g_nm = nm; g_alen = alen(cig)
            if g_nm is None:
                continue
            # best alternative placement of COMPARABLE extent (raw NM)
            best_other = None; best_other_gene = None
            for (chrom, pos, cig, nm) in placs:
                if chrom == g['chrom'] and g['start'] - 10 <= pos <= g['end']:
                    continue
                if alen(cig) < EXTENT_FRAC * g_alen:        # partial hit, not a competitor
                    continue
                if best_other is None or nm < best_other:
                    best_other = nm
                    best_other_gene = gene_at(chrom, pos)
            if best_other is None:
                ga += 1; continue                            # uniquely G's
            dnm = best_other - g_nm                           # distinguishing mismatches
            margins.append(dnm)
            if dnm >= DNM_T:
                ga += 1
            elif dnm <= -DNM_T:
                sa += 1
                sib_anchor[best_other_gene or 'intergenic'] += 1
            else:
                amb += 1
        total = ga + sa + amb
        if ga >= MIN_ANCHOR:
            t = 'expressed_real_copy'
        elif ga == 0 and sa > 0:
            t = 'pure_spillover'
        elif amb >= max(sa, ga) and amb >= 3:
            t = 'ambiguous_dominated'
        else:
            t = 'weak'
        tax[(h['klass'], t)] += 1
        dom_sib = sib_anchor.most_common(1)[0] if sib_anchor else (None, 0)
        o = dict(h)
        o.update(nm_anchor_G=ga, nm_anchor_sibling=sa, nm_anchor_ambiguous=amb,
                 nm_anchor_total=total, nm_taxonomy=t,
                 nm_dom_sibling=dom_sib[0], nm_dom_sibling_n=dom_sib[1],
                 nm_margin_median=round(st.median(margins), 4) if margins else None)
        out_rows.append(o)
    json.dump(out_rows, open(out, 'w'), indent=1)
    print("NM taxonomy (klass, nm_taxonomy):", file=sys.stderr)
    for k, v in sorted(tax.items()):
        print(f"   {k}: {v}", file=sys.stderr)

if __name__ == '__main__':
    main()
