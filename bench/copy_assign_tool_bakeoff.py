#!/usr/bin/env python3
"""O2 hard-locus tool bakeoff (docs/PREREG_tool_bakeoff_2026-09-08.md, PREREG hard_locus_bakeoff 5ca5c7e4), one script,
two modes.

  calls    Derive a per-molecule COPY CALL from any isoform tool's GTF, by one identical rule.
           usage: copy_assign_tool_bakeoff.py calls <tool.gtf> <bam> <copies.tsv> --label NAME [--own assignments.tsv] [--out PREFIX] [--restrict names.txt] [--fuzz N]
  compare  PREREG hard_locus_bakeoff (5ca5c7e4): compare per-tool derived calls (tool_bakeoff.py --out *.calls.tsv)
           usage: copy_assign_tool_bakeoff.py compare --assign ours.assignments.tsv [--gtf ours.gtf] [--bam B --min-mult N] ours=a.calls.tsv flair=b.calls.tsv ...

Derive a per-molecule COPY CALL from any isoform tool's GTF, by one identical rule.

Why this exists: scoring tools on "copy attribution accuracy" is vacuous — StringTie, flair and
isoseq collapse emit no copy attribute, so they score 0 by construction (docs/PREREG_tool_bakeoff_2026-09-08.md).
Every tool DOES emit transcripts with genomic coordinates, and a molecule's intron chain either matches a
transcript or does not. So the copy call is DERIVABLE for all of them, by the same code path.

  molecule --(exact intron chain)--> transcript(s) --(containment in a copy interval)--> copy

States, identical for every arm:
  derived_one    every matching transcript sits in ONE copy
  derived_multi  matching transcripts span >= 2 copies  (the tool conflated copies)
  derived_none   no transcript in the set carries this molecule's chain

Nothing here is tool-specific. A tool's own declared state (ours has assigned/undecidable/unadjudicated)
is reported separately via --own and never enters the derived columns.

usage: tool_bakeoff.py <tool.gtf> <bam> <copies.tsv> --label NAME [--own assignments.tsv] [--out PREFIX] [--restrict names.txt]
"""
from collections import Counter, defaultdict
import argparse, csv, re
import sys, re, csv, collections, subprocess


def calls(argv):

    ap = argv
    gtf, bam, copies_p = ap[0], ap[1], ap[2]
    def opt(k, d=None):
        return ap[ap.index(k) + 1] if k in ap else d
    label = opt('--label', 'tool')
    # Junction tolerance. isoseq collapse runs --max-fuzzy-junction 5 by default, so ITS transcript
    # junctions may sit up to 5 bp off the read's, while our GTF is an exact intron-chain collapse and
    # matches by construction. Scoring at fuzz 0 therefore biases `derived_none` in OUR favour. Applied
    # symmetrically to every arm; report both 0 and 5.
    FUZZ = int(opt('--fuzz', '0'))
    own_p = opt('--own')
    out_p = opt('--out')
    # PREREG hard_locus_bakeoff: score only the listed molecules (one read name per line)
    restrict = set(l.strip() for l in open(opt('--restrict'))) if opt('--restrict') else None

    # ---------------------------------------------------------------- copies
    cop = list(csv.DictReader(open(copies_p), delimiter='\t'))
    copies = [(r['chrom'], int(r['start']), int(r['end']), int(r['copy_idx'])) for r in cop]
    by_chrom = collections.defaultdict(list)
    for c, s, e, i in copies:
        by_chrom[c].append((s, e, i))
    for c in by_chrom:
        by_chrom[c].sort()

    def copy_of(chrom, s, e):
        """The copy a transcript belongs to: max reciprocal overlap, ties -> lowest copy_idx.
        Returns None when the transcript touches no copy at all."""
        best, best_ov = None, 0
        for cs, ce, ci in by_chrom.get(chrom, ()):
            ov = min(e, ce) - max(s, cs)
            if ov > best_ov:
                best, best_ov = ci, ov
        return best

    # ---------------------------------------------------------------- transcripts
    ex = collections.defaultdict(list)
    meta = {}
    for line in open(gtf):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        m = re.search(r'transcript_id[ =]"?([^";]*)"?', f[8])
        if not m:
            continue
        t = m.group(1)
        if f[2] == 'exon':
            ex[t].append((int(f[3]) - 1, int(f[4]), f[0]))

    # chain -> set of copies asserting it; also keep unspliced transcripts by span
    chain_to_copies = collections.defaultdict(set)
    unspliced = []          # (chrom, start, end, copy)
    n_tx = n_tx_in_copy = 0
    tx_multi = []
    copies_hit = set()
    fuzzy_chains = {}
    for t, v in ex.items():
        v.sort()
        chrom = v[0][2]
        s, e = v[0][0], v[-1][1]
        n_tx += 1
        ci = copy_of(chrom, s, e)
        if ci is None:
            continue
        n_tx_in_copy += 1
        # CONFLATION, measured on the transcript rather than on the chain. An intron chain carries
        # genomic coordinates, so two copies can NEVER assert the same chain and a chain-level
        # "spans >= 2 copies" test can never fire (found 2026-09-08, before the flair arm ran;
        # PREREG amendment 1). A transcript that OVERLAPS >= 2 copy intervals is the real thing.
        n_ov = sum(1 for cs, ce, _ in by_chrom.get(chrom, ()) if min(e, ce) - max(s, cs) > 0)
        if n_ov >= 2:
            tx_multi.append((t, chrom, s, e, n_ov))
        copies_hit.add(ci)
        chain = tuple((a[1], b[0]) for a, b in zip(v, v[1:]))
        if chain:
            # index every junction under a rounded key so a fuzzy lookup is O(1) per read
            key = (chrom,) + tuple((a // (2 * FUZZ + 1), b // (2 * FUZZ + 1)) for a, b in chain) if FUZZ else (chrom,) + chain
            chain_to_copies[key].add(ci)
            if FUZZ:
                fuzzy_chains.setdefault(key, []).append((chain, ci))
        else:
            unspliced.append((chrom, s, e, ci))

    # ---------------------------------------------------------------- molecules
    def introns(pos, cig):
        o, p = [], pos
        for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cig):
            n = int(n)
            if op in 'M=XD':
                p += n
            elif op == 'N':
                o.append((p, p + n))
                p += n
        return tuple(o)

    regions = []
    for c in by_chrom:
        cur = None
        for s, e, _ in by_chrom[c]:
            if cur and s <= cur[1]:
                cur[1] = max(cur[1], e)
            else:
                cur = [s, e]
                regions.append((c, cur))

    state = collections.Counter()
    mol_call = {}
    seen = set()
    for c, (lo, hi) in regions:
        # -F 2308: primary, mapped, non-supplementary -- the invariant for per-read statistics
        p = subprocess.run(['samtools', 'view', '-F', '2308', bam, f'{c}:{lo+1}-{hi}'],
                           capture_output=True, text=True)
        for ln in p.stdout.splitlines():
            f = ln.split('\t', 6)
            name, pos, cig = f[0], int(f[3]) - 1, f[5]
            if name in seen or (restrict is not None and name not in restrict):
                continue
            seen.add(name)
            ich = introns(pos, cig)
            if ich:
                if FUZZ:
                    cps = set()
                    # a read matches a transcript when every junction is within FUZZ bp. The bucket key
                    # can straddle a boundary, so probe the neighbouring bucket on each coordinate too.
                    base = 2 * FUZZ + 1
                    seenk = set()
                    for d1 in (0, -1, 1):
                        for d2 in (0, -1, 1):
                            k = (c,) + tuple(((a // base) + d1, (b // base) + d2) for a, b in ich)
                            if k in seenk:
                                continue
                            seenk.add(k)
                            for cand, ci2 in fuzzy_chains.get(k, ()):
                                if len(cand) == len(ich) and all(
                                        abs(x[0] - y[0]) <= FUZZ and abs(x[1] - y[1]) <= FUZZ
                                        for x, y in zip(cand, ich)):
                                    cps.add(ci2)
                else:
                    cps = chain_to_copies.get((c,) + ich, set())
            else:
                # empty-chain trap (register 757): an unspliced read matches every unspliced
                # transcript's empty chain, so it must be resolved by SPAN CONTAINMENT instead.
                end = pos
                for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cig):
                    if op in 'M=XDN':
                        end += int(n)
                cps = {ci for tc, ts, te, ci in unspliced if tc == c and ts <= pos and end <= te}
            if not cps:
                st = 'derived_none'
            elif len(cps) == 1:
                st = 'derived_one'
            else:
                st = 'derived_multi'
            state[st] += 1
            mol_call[name] = (st, sorted(cps))

    tot = sum(state.values())
    print(f'== {label}')
    print(f'   transcripts {n_tx}  in a copy {n_tx_in_copy}')
    print(f'   transcripts spanning >=2 copies {len(tx_multi)}  ({len(tx_multi)/n_tx_in_copy:.3f} of in-copy)'
          if n_tx_in_copy else '   transcripts spanning >=2 copies 0')
    print(f'   copies with >=1 transcript      {len(copies_hit)} / {len(copies)}')
    print(f'   molecules   {tot}')
    for k in ('derived_one', 'derived_multi', 'derived_none'):
        v = state[k]
        print(f'   {k:14s} {v:7d}  {v/tot:6.3f}' if tot else f'   {k:14s} {v:7d}')

    # ---------------------------------------------------------------- the tool's OWN state, if any
    if own_p:
        own = collections.Counter()
        own_call = {}
        with open(own_p) as fh:
            rd = csv.DictReader(fh, delimiter='\t')
            cols = rd.fieldnames or []
            sc = next((c for c in ('status', 'decision', 'call') if c in cols), None)
            nc = next((c for c in ('read', 'read_name', 'molecule', 'name') if c in cols), None)
            # catalog_copy_idx FIRST: `assigned_copy` is the SWEEP index and does not address
            # copies.tsv (register 756 -- the two were once printed side by side and disagreed).
            cc = next((c for c in ('catalog_copy_idx', 'copy_idx', 'copy') if c in cols), None)
            for r in rd:
                if not sc or not nc:
                    break
                own[r[sc]] += 1
                if cc:
                    own_call[r[nc]] = (r[sc], r[cc])
        if own:
            print(f'   -- own declared states ({sc}):')
            for k, v in own.most_common():
                print(f'      {k:20s} {v:7d}')
        if own_call:
            # Only rows the tool itself CALLED count. Including its abstentions would compare the
            # derived rule against a non-decision and inflate the denominator.
            DECIDED = {'assigned'}
            agree = dis = 0
            for n, (st, cps) in mol_call.items():
                o = own_call.get(n)
                if o and o[0] in DECIDED and st == 'derived_one' and o[1].isdigit():
                    if int(o[1]) == cps[0]:
                        agree += 1
                    else:
                        dis += 1
            if agree + dis:
                print(f'   -- derived vs own copy, where the TOOL decided: {agree}/{agree+dis} = {agree/(agree+dis):.3f}')

    if out_p:
        with open(out_p + '.calls.tsv', 'w') as fh:
            fh.write('molecule\tstate\tcopies\n')
            for n, (st, cps) in sorted(mol_call.items()):
                fh.write(f'{n}\t{st}\t{",".join(map(str, cps))}\n')
        print(f'   wrote {out_p}.calls.tsv')



def load_calls(p):
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["molecule"]] = (r["state"], [int(x) for x in r["copies"].split(",") if x])
    return out


def compare(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--assign", required=True)
    ap.add_argument("--gtf", help="our GTF, for copy_index of the transcript carrying each molecule (P6)")
    ap.add_argument("--bam", help="with --min-mult: primaries (-F 2308) give each molecule's intron chain")
    ap.add_argument("--min-mult", type=int, default=0, help="keep only molecules whose exact chain is carried by >= N molecules (support-policy control)")
    ap.add_argument("tools", nargs="+", help="label=calls.tsv; the first is ours")
    a = ap.parse_args(argv)
    tools = [(t.split("=")[0], load_calls(t.split("=")[1])) for t in a.tools]
    ours_l, ours = tools[0]
    assign = {r["read_name"]: r for r in csv.DictReader(open(a.assign), delimiter="\t")}
    allm = set(ours)
    for _, c in tools:
        allm |= set(c)
    if a.min_mult and a.bam:
        import re as _re, subprocess
        def introns(pos, cig):
            o, p = [], pos
            for n, op in _re.findall(r"(\d+)([MIDNSHP=X])", cig):
                n = int(n)
                if op in "M=XD":
                    p += n
                elif op == "N":
                    o.append((p, p + n)); p += n
            return tuple(o)
        chain = {}
        out = subprocess.run(["samtools", "view", "-F", "2308", a.bam], capture_output=True, text=True).stdout
        for ln in out.splitlines():
            f = ln.split("\t", 6)
            chain.setdefault(f[0], (f[2],) + introns(int(f[3]) - 1, f[5]))
        mult = Counter(chain.values())
        keep = {m for m in allm if mult.get(chain.get(m), 0) >= a.min_mult}
        print(f"--min-mult {a.min_mult}: keeping {len(keep)} of {len(allm)} molecules whose exact chain has >= {a.min_mult} molecules")
        allm = keep
    hard = {m for m in allm if m in assign}
    easy = allm - hard
    carried = lambda c, m: c.get(m, ("derived_none", []))[0] != "derived_none"
    strata = {
        "hard (all gate rows)": hard,
        "  contested": {m for m in hard if assign[m]["origin_rejected"] == "0" and int(assign[m]["n_candidates"]) >= 2},
        "    assigned": {m for m in hard if assign[m]["status"] == "assigned" and assign[m]["origin_rejected"] == "0" and int(assign[m]["n_candidates"]) >= 2},
        "    tied": {m for m in hard if assign[m]["status"] == "tied" and assign[m]["origin_rejected"] == "0" and int(assign[m]["n_candidates"]) >= 2},
        "    ambiguous": {m for m in hard if assign[m]["status"] == "ambiguous" and assign[m]["origin_rejected"] == "0" and int(assign[m]["n_candidates"]) >= 2},
        "  tie outside catalog": {m for m in hard if assign[m].get("tie_outside_catalog") == "1"},
        "  origin-rejected": {m for m in hard if assign[m]["origin_rejected"] == "1"},
        "easy (not admitted by the gate)": easy,
    }
    print(f"molecules: {len(allm)} total, hard {len(hard)}, easy {len(easy)}")
    print(f"\n{'stratum':34s} {'n':>6} " + " ".join(f"{l:>10}" for l, _ in tools) + "   (fraction carried = derived_one|multi)")
    rates = {}
    for name, S in strata.items():
        row = []
        for l, c in tools:
            k = sum(1 for m in S if carried(c, m))
            rates[(name, l)] = k / len(S) if S else float("nan")
            row.append(f"{k/len(S):10.3f}" if S else f"{'-':>10}")
        print(f"{name:34s} {len(S):>6} " + " ".join(row))
    # P2 discordance on hard and contested
    for name in ("hard (all gate rows)", "  contested"):
        S = strata[name]
        print(f"\n{name.strip()}: discordance vs ours")
        for l, c in tools[1:]:
            on = sum(1 for m in S if carried(ours, m) and not carried(c, m))
            xn = sum(1 for m in S if carried(c, m) and not carried(ours, m))
            both = sum(1 for m in S if carried(ours, m) and carried(c, m))
            nei = len(S) - on - xn - both
            print(f"   {l:10s} ours-not-{l}: {on:5d}   {l}-not-ours: {xn:5d}   ratio {on/max(1,xn):5.2f}   both {both}  neither {nei}")
    # P3 hard vs easy derived_none
    print("\nP3 derived_none hard / easy:")
    for l, c in tools:
        h = 1 - rates[("hard (all gate rows)", l)]; e = 1 - rates[("easy (not admitted by the gate)", l)]
        print(f"   {l:10s} hard {h:.3f}  easy {e:.3f}  ratio {h/e if e else float('nan'):.2f}")
    # P4 per copy: copies where ours carries >=1 hard molecule and X carries none (by the derived copy)
    print("\nP4 per-copy coverage on the hard set (copies where the tool carries >=1 hard molecule, by derived copy):")
    cov = {}
    for l, c in tools:
        cs = set()
        for m in hard:
            st, cps = c.get(m, ("derived_none", []))
            if st != "derived_none":
                cs.update(cps)
        cov[l] = cs
    for l, _ in tools[1:]:
        print(f"   {l:10s} copies {len(cov[l]):2d}; ours-only {sorted(cov[ours_l]-cov[l])}  {l}-only {sorted(cov[l]-cov[ours_l])}")
    print(f"   {ours_l:10s} copies {len(cov[ours_l])}: {sorted(cov[ours_l])}")
    # P6 copy attribution of the O2-assigned
    asg = strata["    assigned"]
    if a.gtf:
        tx_copy = {}
        for line in open(a.gtf):
            if "\ttranscript\t" not in line:
                continue
            t = re.search(r'transcript_id "([^"]+)"', line); ci = re.search(r'copy_index "([^"]+)"', line)
            if t and ci:
                tx_copy[t.group(1)] = ci.group(1)
    print(f"\nP6 the {len(asg)} O2-assigned molecules (report):")
    for l, c in tools:
        one = [m for m in asg if c.get(m, ("derived_none", []))[0] == "derived_one"]
        agree = sum(1 for m in one if str(c[m][1][0]) == assign[m]["catalog_copy_idx"])
        print(f"   {l:10s} carried {sum(1 for m in asg if carried(c, m)):3d}/{len(asg)}; derived_one {len(one):3d}, of which derived copy == O2 copy: {agree} ({100*agree/max(1,len(one)):.0f}%)")


if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in ('calls', 'compare'):
        sys.exit(__doc__)
    {'calls': calls, 'compare': compare}[sys.argv[1]](sys.argv[2:])
