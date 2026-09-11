#!/usr/bin/env python3
"""PREREG isoform_bakeoff (84e8ed50): lift-aware hard-locus bakeoff. A tool carries a molecule if a transcript's
chain equals the read's chain at its primary copy, or its LIFT to another copy (+-tol per boundary). The tool's
address set for that isoform = the copies where it emits the chain (+ our `copies` attribute).

  python3 bench/isoform_bakeoff.py --assign ours_final2.assignments.tsv --copies copies16.tsv --paf human_gspans.paf \
      --bam hsa16.bam ours=ours_copyset.gtf ours_raw=ours_final2.gtf flair=flair_family.gtf stringtie=stringtie_family.gtf isoseq=isoseq_family.gff
"""
import argparse, csv, os, re, subprocess, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from isoform_copy_lift import build_lifts, make_lift_pos, introns_of


def main():
    ap = argparse.ArgumentParser()
    for k in ("--assign", "--copies", "--paf", "--bam"):
        ap.add_argument(k, required=True)
    ap.add_argument("--tol", type=int, default=5); ap.add_argument("--min-mult", type=int, default=0)
    ap.add_argument("--summary", action="store_true",
                     help="B6 (docs/OPEN_ITEMS_2026-09-09.md): print the ONE composite headline per tool "
                          "instead of 'carried' alone (register row 803: lift-aware carried rewards a tool "
                          "that duplicates one isoform across several copies, e.g. isoseq's 45 phantom "
                          "transcripts / 200 arbitrary addresses). The composite is carried AND "
                          "address-correct (address subset of the molecule's own lift-reachable copy set) "
                          "AND reported beside the tool's own phantom rate (transcripts with zero backing "
                          "molecule) — never quote 'carried' by itself again.")
    ap.add_argument("tools", nargs="+")
    a = ap.parse_args()
    cop = {r["copy_idx"]: (r["chrom"], int(r["start"]), int(r["end"])) for r in csv.DictReader(open(a.copies), delimiter="\t")}
    span2idx = {f"{c}:{s+1}-{e}": i for i, (c, s, e) in cop.items()}
    lifts, ident = build_lifts(a.paf, span2idx); lift_pos = make_lift_pos(lifts, cop)
    def copy_of(chrom, s, e):
        best, bo = None, 0
        for i, (c, cs, ce) in cop.items():
            if c != chrom: continue
            ov = min(e, ce) - max(s, cs)
            if ov > bo: best, bo = i, ov
        return best
    # tools: (chrom, n_introns) -> list of (chain, address_set)
    tools = []
    for spec in a.tools:
        label, path = spec.split("=", 1)
        ex = defaultdict(list); attrs = {}
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or line.startswith("#"): continue
            t = re.search(r'transcript_id[ =]"?([^";]*)"?', f[8])
            if not t: continue
            t = t.group(1)
            if f[2] == "exon": ex[t].append((int(f[3]) - 1, int(f[4]), f[0]))
            elif f[2] == "transcript": attrs[t] = f[8]
        idx = defaultdict(list)
        for t, v in ex.items():
            v.sort(); chrom = v[0][2]; c = copy_of(chrom, v[0][0], v[-1][1])
            if c is None: continue
            chain = tuple((p[1], q[0]) for p, q in zip(v, v[1:]))
            if not chain: continue
            addr = {c}
            m = re.search(r' copies "([^"]*)"', attrs.get(t, ""))
            # A6 (`--name-outside-tie`) can print the outside partner as `outside:chrom:start-end` instead
            # of the bare `outside` this script already special-cases everywhere below (`- {"outside"}`,
            # `| {"outside"}`) and `sorted(..., key=int)` assumes every other token is a catalog index —
            # normalize any `outside*` token to the bare word here so both stay correct either way.
            if m: addr |= {("outside" if x.startswith("outside") else x) for x in m.group(1).split(",") if x}
            idx[(chrom, len(chain))].append((chain, addr))
        tools.append((label, idx))
    def match(idx, chrom, chain):
        out = set()
        for ch, addr in idx.get((chrom, len(chain)), ()):
            if all(abs(x0 - y0) <= a.tol and abs(x1 - y1) <= a.tol for (x0, x1), (y0, y1) in zip(chain, ch)):
                out |= addr
        return out
    # hard molecules with a primary in a copy span
    assign = {r["read_name"]: r for r in csv.DictReader(open(a.assign), delimiter="\t")}
    reads = {}; mult = Counter()
    out = subprocess.run(["samtools", "view", "-F", "2308", a.bam], capture_output=True, text=True).stdout
    for ln in out.splitlines():
        f = ln.split("\t", 6); name, chrom, pos, cig = f[0], f[2], int(f[3]) - 1, f[5]
        if name in reads: continue
        ich = introns_of(pos, cig)
        mult[(chrom,) + ich] += 1
        if name not in assign or not ich: continue
        end = pos + sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig) if op in "M=XDN")
        P = copy_of(chrom, pos, end)
        if P is None: continue
        reads[name] = (chrom, ich, P)
    if a.min_mult:
        reads = {n: v for n, v in reads.items() if mult[(v[0],) + v[1]] >= a.min_mult}
        print(f"--min-mult {a.min_mult}: {len(reads)} hard spliced molecules kept")
    # lifted chains per read
    def lifted_chains(chrom, ich, P):
        res = {P: ich}
        for c in cop:
            if c == P or (P, c) not in lifts: continue
            lc = []; ok = True
            for (x0, x1) in ich:
                l0, d0 = lift_pos(x0, P, c); l1, d1 = lift_pos(x1, P, c)
                if l0 is None or l1 is None or d0 > a.tol or d1 > a.tol: ok = False; break
                lc.append((min(l0, l1), max(l0, l1)))
            if ok: res[c] = tuple(sorted(lc))
        return res
    carried = {label: {} for label, _ in tools}
    touched = {label: set() for label, _ in tools}  # (chrom, chain) keys of idx actually matched by >=1 molecule
    lifted_of = {}  # name -> lc, kept for --summary's address-correctness check
    for name, (chrom, ich, P) in reads.items():
        lc = lifted_chains(chrom, ich, P)
        lifted_of[name] = lc
        for label, idx in tools:
            addr = set()
            for c, ch in lc.items():
                for tch, taddr in idx.get((chrom, len(ch)), ()):
                    if all(abs(x0 - y0) <= a.tol and abs(x1 - y1) <= a.tol for (x0, x1), (y0, y1) in zip(ch, tch)):
                        addr |= taddr
                        touched[label].add((chrom, tch))
            if addr: carried[label][name] = addr
    def is_con(r): return r["origin_rejected"] == "0" and int(r["n_candidates"]) >= 2
    strata = {"hard (spliced, primary in a copy)": set(reads),
              "  contested": {n for n in reads if is_con(assign[n])},
              "    assigned": {n for n in reads if is_con(assign[n]) and assign[n]["status"] == "assigned"},
              "    undecided (tied+ambiguous)": {n for n in reads if is_con(assign[n]) and assign[n]["status"] != "assigned"}}
    print(f"hard spliced molecules: {len(reads)}")
    print(f"{'stratum':36s} {'n':>5} " + " ".join(f"{l:>10}" for l, _ in tools) + "   (lift-aware carried)")
    for s, S in strata.items():
        print(f"{s:36s} {len(S):>5} " + " ".join(f"{sum(1 for n in S if n in carried[l])/max(1,len(S)):10.3f}" for l, _ in tools))
    asg = strata["    assigned"]; und = strata["    undecided (tied+ambiguous)"]
    print("\nP2 O2-assigned molecules carried: address set contains O2's copy")
    for l, _ in tools:
        c = [n for n in asg if n in carried[l]]; k = sum(1 for n in c if assign[n]["catalog_copy_idx"] in carried[l][n])
        one = sum(1 for n in c if len(carried[l][n]) == 1)
        print(f"   {l:10s} carried {len(c):3d}/{len(asg)}; contains O2 copy {k} ({100*k/max(1,len(c)):.0f}%); single-address {one}")
    print("\nP3 O2-undecided molecules carried: single address vs a set")
    for l, _ in tools:
        c = [n for n in und if n in carried[l]]; one = sum(1 for n in c if len(carried[l][n]) == 1)
        print(f"   {l:10s} carried {len(c):3d}/{len(und)}; ONE address {one} ({100*one/max(1,len(c)):.0f}%); set >= 2 {len(c)-one} ({100*(len(c)-one)/max(1,len(c)):.0f}%)")
    print("\nP4 per copy (addresses over the hard set):")
    cov = {l: set().union(*carried[l].values()) - {"outside"} if carried[l] else set() for l, _ in tools}
    ours = tools[0][0]
    for l, _ in tools[1:]:
        print(f"   {l:10s} copies {len(cov[l]):2d}; {ours}-only {sorted(cov[ours]-cov[l], key=int)}  {l}-only {sorted(cov[l]-cov[ours], key=int)}")
    print(f"   discordance on the hard set vs {ours}: " + ", ".join(f"{l}: ours-not-{l} {sum(1 for n in reads if n in carried[ours] and n not in carried[l])} / {l}-not-ours {sum(1 for n in reads if n in carried[l] and n not in carried[ours])}" for l, _ in tools[1:]))

    if a.summary:
        # B6: the ONE composite to quote — never "carried" alone (register row 803). A molecule counts only
        # if carried AND its claimed address never reaches beyond copies the molecule's OWN lifted chain
        # could legitimately sit at (an over-claimed address is exactly the isoseq-duplication failure mode
        # row 803 found: 45 phantom transcripts, 200 arbitrary addresses, that inflated "carried" alone).
        print(f"\nP7 (B6) THE composite headline — carried AND address-correct AND the tool's own phantom rate:")
        n_hard = len(reads)
        for label, idx in tools:
            composite = sum(
                1 for name in reads
                if name in carried[label] and carried[label][name] <= (set(lifted_of[name].keys()) | {"outside"})
            )
            total_tx = sum(len(v) for v in idx.values())
            n_touched = len(touched[label])
            phantom = total_tx - n_touched
            print(f"   {label:10s} composite {composite:4d}/{n_hard} = {composite/max(1,n_hard):.3f}   "
                  f"(carried {sum(1 for n in reads if n in carried[label])}/{n_hard} = "
                  f"{sum(1 for n in reads if n in carried[label])/max(1,n_hard):.3f}, "
                  f"phantom {phantom}/{total_tx} = {phantom/max(1,total_tx):.3f} of its own transcripts)")


if __name__ == "__main__":
    main()
