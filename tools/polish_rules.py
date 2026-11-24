#!/usr/bin/env python3
"""The assembler's polishing rules, in Python: which transcripts are removed and by which rule.

    polish_rules.py UNPOLISHED.gtf OUT.gtf OUT.reasons.tsv [--mono-quantile 0.82] [--isoform-fraction 0.02] [--ism-ratio 0.3]
                    [--retained-ratio 10] [--no-mono-shadow] [--mode full|light]

A line-by-line port of `polish_drop_set_with` in src/copy_assign_bin.rs (the Rust assembler applies it per contig to the
transcripts it has just built, each with its read count in the attribute `reads`). It is here to be READ: every rule is
a short block below, in the order the assembler applies them, and OUT.reasons.tsv names the rule that removed each
transcript. Defaults are the driver's (tools/rustle_pipeline.sh). Parity with the Rust output is checked by
`polish_rules.py --parity UNPOLISHED.gtf RUST_POLISHED.gtf`.

The rules, in order (a transcript removed by one rule is not looked at by the next):
  fuzzy     (only with --fuzzy N > 0; off by default) two chains of the same length whose junctions all lie within N
            bases: the one with fewer reads goes.
  ism       sub-chain collapse: a spliced transcript whose intron chain is a contiguous piece of a longer chain of the
            same contig and strand is removed unless it has at least --ism-ratio x the longer one's reads.
  ism_mono  a single-exon transcript lying inside the span of a spliced transcript (the one with the longest chain that
            contains it) is removed unless it has at least --ism-ratio x that transcript's reads.
  floor     single-exon read floor: the --mono-quantile quantile of the read counts of the contig's surviving spliced
            transcripts; a single-exon transcript with fewer reads goes.
  shadow    a surviving single-exon transcript that overlaps an exon of any spliced transcript (either strand) or lies
            in the span of one on its own strand goes.
  frac      isoform fraction: a transcript with fewer reads than --isoform-fraction x the best transcript of its gene.
  ret       retained intron: a transcript one of whose exons contains an intron that other transcripts of its gene
            splice, when those have at least --retained-ratio x its reads.
Not in this file (they run before or after it inside the assembler and need the reads): junction snapping, the spliced
share test, the 5' end rescue, the internal-priming test of single-exon transcripts, junction replacement by support ratio.
Parity (2026-11-24, those steps switched off; workspace polish_parity_2026-11-24): two planted sets 267 = 267 and 301 = 301
transcripts kept; gorilla OR6737 every rule's count equal to the assembler's (ism 58,333, floor + shadow 28,772, frac 28,907,
ret 2,392), kept sets equal but for the 92 single-exon transcripts the internal-priming test removes."""
import sys, re, argparse, collections, bisect

# The assembler regroups genes AFTER polishing (bridge regrouping: gene ids get a suffix .rgN or .fusN). A file written with
# polishing off already carries those suffixes, so they are removed here to judge every transcript in the gene it had when
# the rules ran. Checked on gorilla OR6737: with this, every rule removes exactly the assembler's number of transcripts.
REGROUP = re.compile(r'\.(rg|fus)\d+$')

def attr(s, k):
    m = re.search(k + r' "([^"]*)"', s)
    return m.group(1) if m else None

def drop_set(lines, mode='full', mono_quantile=0.82, isoform_fraction=0.02, mono_shadow=True, ism_escape=False, ism_3p=False,
             ism_ratio=0.3, fraction_exempt=False, fuzzy=0, fuzzy_ism=False, fraction_min_reads=0, retained_ratio=10.0, shadow_ratio=0.0):
    """lines: the GTF lines of ONE contig. -> {transcript: rule that removed it}, the single-exon floor"""
    if mode == 'none': return {}, 0
    exons = collections.defaultdict(list); key = {}; reads = {}; gene = {}
    for line in lines:
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        tid = attr(f[8], 'transcript_id')
        if tid is None: continue
        r = attr(f[8], 'reads')
        if r is not None and r.isdigit(): reads[tid] = max(reads.get(tid, 0), int(r))
        g = attr(f[8], 'gene_id')
        if g is not None: gene.setdefault(tid, REGROUP.sub('', g))     # the gene as it was when the assembler polished (see REGROUP)
        if f[2] != 'exon': continue
        exons[tid].append((int(f[3]) - 1, int(f[4]))); key.setdefault(tid, (f[0], f[6]))
    chain = {}; span = {}
    for tid, ex in exons.items():
        ex.sort(); chain[tid] = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]; span[tid] = (ex[0][0], ex[-1][1])
    R = lambda t: reads.get(t, 0)
    def quantile(v):
        v = sorted(v)
        return 0 if not v or mono_quantile <= 0 else v[min(int(mono_quantile * len(v)), len(v) - 1)]
    support_level = quantile([R(t) for t, c in chain.items() if c])
    why = {}
    def mark(ts, rule):
        for t in ts: why.setdefault(t, rule)
    drop = set()

    # ---- fuzzy: near-identical chains (off unless --fuzzy)
    if fuzzy > 0:
        buckets = collections.defaultdict(list)
        for t, c in chain.items():
            if c and t in key: buckets[(key[t][0], key[t][1], len(c))].append(t)
        for bk in sorted(buckets):
            ts = sorted(buckets[bk], key=lambda t: (-R(t), t))
            for i in range(len(ts)):
                if ts[i] in drop: continue
                for j in range(i + 1, len(ts)):
                    if ts[j] not in drop and all(abs(a[0] - b[0]) <= fuzzy and abs(a[1] - b[1]) <= fuzzy for a, b in zip(chain[ts[i]], chain[ts[j]])): drop.add(ts[j])
    mark(drop, 'fuzzy')

    if mode == 'full':
        def supported(frag, cont):
            rf = R(frag)
            if ism_escape and support_level > 0 and rf >= support_level: return True
            rc = R(cont)
            return rc > 0 and rf >= ism_ratio * rc
        groups = collections.defaultdict(list)
        for tid, k in key.items(): groups[k].append(tid)
        for k, tids in groups.items():
            # ---- ism: a chain that is a contiguous piece of a longer one
            multi = sorted((t for t in tids if chain[t]), key=lambda t: (-len(chain[t]), t))
            exact_index = not (fuzzy_ism and fuzzy > 0)
            by_first = collections.defaultdict(list)
            if exact_index:
                for y in multi: by_first[chain[y][0]].append(y)
            tol = fuzzy if fuzzy_ism else 0
            for x in multi:
                if x in drop: continue
                cx = chain[x]
                cand = [y for j in cx for y in by_first.get(j, ())] if exact_index else multi
                for y in cand:
                    if y == x or y in drop: continue
                    cy = chain[y]
                    if len(cy) >= len(cx): continue
                    offs = ([0] if k[1] == '-' else [len(cx) - len(cy)]) if ism_3p else range(len(cx) - len(cy) + 1)
                    if any(all(abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol for a, b in zip(cx[o:o + len(cy)], cy)) for o in offs) and not supported(y, x): drop.add(y)
            # ---- ism_mono: a single-exon transcript inside a spliced one's span, judged against the longest such chain
            alive = [(r, span[m]) for r, m in enumerate(multi) if m not in drop]
            starts = sorted({s[0] for _, s in alive}); n = len(starts); fen = [None] * (n + 1)
            by_end = sorted(((s[1], r, s[0]) for r, s in alive), key=lambda z: (-z[0], z[1]))
            monos = sorted(((t, span[t]) for t in sorted(tids) if not chain[t] and t not in drop), key=lambda z: (-z[1][1], z[0]))
            p = 0
            for t, s in monos:
                while p < len(by_end) and by_end[p][0] >= s[1]:           # spliced transcripts that end at or after this one
                    i = bisect.bisect_left(starts, by_end[p][2]) + 1
                    while i <= n:
                        fen[i] = by_end[p][1] if fen[i] is None else min(fen[i], by_end[p][1]); i += i & -i
                    p += 1
                i = bisect.bisect_right(starts, s[0]); best = None      # ... and start at or before it: the longest chain among them
                while i > 0:
                    if fen[i] is not None: best = fen[i] if best is None else min(best, fen[i])
                    i -= i & -i
                if best is not None and not supported(t, multi[best]): drop.add(t)
    for t in drop: why.setdefault(t, 'ism' if chain.get(t) else 'ism_mono')

    # ---- floor: single-exon transcripts under the quantile of the surviving spliced ones
    floor = quantile([R(t) for t, c in chain.items() if c and t not in drop])
    if floor > 0:
        for t, c in chain.items():
            if not c and t not in drop and R(t) < floor: drop.add(t)
    mark(drop, 'floor')

    # ---- shadow: a single-exon transcript under an exon of any spliced transcript, or in a spliced span of its strand
    if mono_shadow:
        exons_any = collections.defaultdict(list); spans_same = collections.defaultdict(list)
        for t, c in chain.items():
            if not c or t in drop or t not in key: continue
            k = key[t]; r = R(t)
            spans_same[k].append([span[t][0], span[t][1], r]); exons_any[k[0]] += [[a, b, r] for a, b in exons[t]]
        def merge(v):
            out = []
            for a, b, r in sorted(v):
                if out and a <= out[-1][1]: out[-1][1] = max(out[-1][1], b); out[-1][2] = max(out[-1][2], r)
                else: out.append([a, b, r])
            return out
        exons_any = {k: merge(v) for k, v in exons_any.items()}; spans_same = {k: merge(v) for k, v in spans_same.items()}
        def hits(v, s, e):
            hi = bisect.bisect_right([x[0] for x in v], e); top = None; i = hi - 1
            while i >= 0 and v[i][1] >= s:
                top = v[i][2] if top is None else max(top, v[i][2]); i -= 1
            return top
        for t in sorted(t for t, c in chain.items() if not c and t not in drop):
            if t not in key: continue
            k = key[t]; s, e = span[t]
            a = hits(exons_any[k[0]], s, e) if k[0] in exons_any else None; b = hits(spans_same[k], s, e) if k in spans_same else None
            top = a if b is None else b if a is None else max(a, b)
            if top is not None and not (shadow_ratio > 0 and R(t) >= shadow_ratio * top): drop.add(t)
    mark(drop, 'shadow')

    # ---- frac: far fewer reads than the best transcript of the gene
    if isoform_fraction > 0:
        best = collections.Counter()
        for t, g in gene.items():
            if t not in drop: best[g] = max(best[g], R(t))
        for t in sorted(t for t in gene if t not in drop):
            b = best[gene[t]]; r = R(t)
            if b == 0 or (fraction_exempt and support_level > 0 and r >= support_level) or (fraction_min_reads > 0 and r >= fraction_min_reads): continue
            if r < b and r < isoform_fraction * b: drop.add(t)
    mark(drop, 'frac')

    # ---- ret: an exon that holds an intron the gene's other transcripts splice, with far more reads
    if retained_ratio > 0:
        support = collections.Counter(); gj = collections.defaultdict(set)
        for t, c in chain.items():
            if t in drop or t not in key: continue
            k = key[t]
            for j in c: support[(k[0], k[1], j)] += R(t)
            if t in gene: gj[(gene[t], k[1])].update(c)
        newly = set()
        for t in sorted(t for t in chain if t not in drop):
            if t not in gene or t not in key: continue
            k = key[t]; js = gj.get((gene[t], k[1]))
            if not js: continue
            own = set(chain[t]); rt = R(t)
            for j in sorted(js):
                if j in own: continue
                if any(a < j[0] and j[1] < b for a, b in exons[t]) and support[(k[0], k[1], j)] >= retained_ratio * rt: newly.add(t); break
        drop |= newly
    mark(drop, 'ret')
    return why, floor

def by_contig(path):
    out = collections.defaultdict(list); order = []
    for l in open(path):
        if l.startswith('#') or not l.strip(): continue
        c = l.split('\t', 1)[0]
        if c not in out: order.append(c)
        out[c].append(l)
    return out, order

def apply(path, **kw):
    """-> kept lines, {transcript: (contig, rule)}, {contig: floor}"""
    parts, order = by_contig(path); kept = []; why = {}; floors = {}
    for c in order:
        w, floors[c] = drop_set(parts[c], **kw)
        for t, r in w.items(): why[t] = (c, r)
        kept += [l for l in parts[c] if attr(l.split('\t')[8] if l.count('\t') >= 8 else '', 'transcript_id') not in w]
    return kept, why, floors

def chains_of(lines):
    ex = collections.defaultdict(list)
    for l in lines:
        f = l.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[attr(f[8], 'transcript_id')].append((f[0], f[6], int(f[3]), int(f[4])))
    return {t: tuple(sorted(e)) for t, e in ex.items()}

def selftest():
    def tr(t, g, ex, n, s='+'):
        o = f'c\tx\ttranscript\t{ex[0][0]}\t{ex[-1][1]}\t.\t{s}\t.\tgene_id "{g}"; transcript_id "{t}"; reads "{n}";\n'
        return [o] + [f'c\tx\texon\t{a}\t{b}\t.\t{s}\t.\tgene_id "{g}"; transcript_id "{t}";\n' for a, b in ex]
    L = (tr('full', 'G', [(1, 100), (201, 300), (401, 500)], 100) + tr('piece', 'G', [(201, 300), (401, 500)], 10) + tr('strong_piece', 'G', [(1, 100), (201, 310)], 40)
         + tr('minor', 'G', [(1, 100), (201, 300), (451, 500)], 1) + tr('retained', 'G', [(1, 300), (401, 450), (461, 500)], 5)
         + tr('mono_frag', 'M0', [(221, 280)], 20) + tr('mono_in', 'M1', [(211, 290)], 150) + tr('mono_weak', 'M2', [(2000, 2500)], 2) + tr('mono_ok', 'M3', [(5000, 5500)], 200)
         + tr('other', 'H', [(9001, 9100), (9201, 9300)], 30))
    why, floor = drop_set(L)
    assert why == {'piece': 'ism', 'mono_frag': 'ism_mono', 'mono_in': 'shadow', 'mono_weak': 'floor', 'minor': 'frac', 'retained': 'ret'}, why
    assert floor == 100 and 'strong_piece' not in why and 'mono_ok' not in why
    assert drop_set(L, mode='none') == ({}, 0)
    w2, _ = drop_set(L, ism_ratio=0.0, isoform_fraction=0, retained_ratio=0, mono_quantile=0, mono_shadow=False); assert w2 == {}, w2
    print('polish_rules selftest ok')

def main():
    a = sys.argv[1:]
    if a == ['--selftest']: return selftest()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gtf'); ap.add_argument('out'); ap.add_argument('reasons', nargs='?')
    ap.add_argument('--parity', action='store_true', help='OUT is the Rust-polished GTF of the same run: compare, write nothing')
    ap.add_argument('--mode', default='full'); ap.add_argument('--mono-quantile', type=float, default=0.82); ap.add_argument('--isoform-fraction', type=float, default=0.02)
    ap.add_argument('--ism-ratio', type=float, default=0.3); ap.add_argument('--retained-ratio', type=float, default=10.0); ap.add_argument('--no-mono-shadow', action='store_true')
    ap.add_argument('--fuzzy', type=int, default=0)
    o = ap.parse_args()
    kept, why, floors = apply(o.gtf, mode=o.mode, mono_quantile=o.mono_quantile, isoform_fraction=o.isoform_fraction, ism_ratio=o.ism_ratio,
                              retained_ratio=o.retained_ratio, mono_shadow=not o.no_mono_shadow, fuzzy=o.fuzzy)
    n = collections.Counter(r for _, r in why.values())
    if o.parity:
        mine = set(chains_of(kept).values()); rust = set(chains_of(open(o.out)).values())
        only = mine - rust; mono = sum(1 for c in only if len(c) == 1)
        print(f'polish_rules parity: python keeps {len(mine)} transcripts, rust {len(rust)}; shared {len(mine & rust)}; only python {len(only)} ({mono} of them single-exon: the assembler also removes internally primed single-exon transcripts, a rule that needs the reads and is not in this file); only rust {len(rust - mine)}; removed by rule ' + ', '.join(f'{k} {v}' for k, v in sorted(n.items())))
        return
    with open(o.out, 'w') as w: w.writelines(kept)
    if o.reasons:
        with open(o.reasons, 'w') as w:
            w.write('transcript\tcontig\trule\n'); w.writelines(f'{t}\t{c}\t{r}\n' for t, (c, r) in sorted(why.items(), key=lambda z: (z[1][0], z[0])))
    print(f'polish_rules: {len(why)} transcripts removed (' + ', '.join(f'{k} {v}' for k, v in sorted(n.items())) + f') -> {o.out}')

if __name__ == '__main__':
    main()
