#!/usr/bin/env python3
"""PREREG isoform_pool (12576433): pool a chain's AS-tied, origin-pass reads and run the pairwise certificate
on the pooled trials (the binary's arithmetic: n = distinguishing columns, k = matches to the target,
p = P(Bin(n, e/3) >= k), LLR = lr * (matches A - matches B), bk = maximin, assigned iff every rival rejected
at alpha/(|C|-1) with LLR > 0; tied if some rival shares no distinguishing column).

  python3 bench/isoform_pool.py --dump ours_final2_dump.star_reads.tsv --assign ours_final2.assignments.tsv \
      --bam hsa16.bam --copies copies16.tsv [--gtf ours_final2.gtf --phantoms isoform_groups.tsv] [--out isoforms.tsv]
"""
import argparse, csv, math, re, subprocess
from collections import Counter, defaultdict

E = 0.003; EPS = E / 3.0; LR = math.log((1 - E) / EPS); ALPHA = 1e-3


def binom_tail(n, k):
    """P(X >= k), X ~ Bin(n, EPS), in log space."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    lp, l1p = math.log(EPS), math.log1p(-EPS)
    tot = 0.0
    for i in range(k, n + 1):
        tot += math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * lp + (n - i) * l1p)
    return min(1.0, tot)


def introns_of(pos, cig):
    o, p = [], pos
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
        n = int(n)
        if op in "M=XD":
            p += n
        elif op == "N":
            o.append((p, p + n)); p += n
    return tuple(o)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True); ap.add_argument("--assign", required=True)
    ap.add_argument("--bam", required=True); ap.add_argument("--copies", required=True)
    ap.add_argument("--gtf"); ap.add_argument("--phantoms"); ap.add_argument("--out")
    ap.add_argument("--min-reads", type=int, default=2)
    a = ap.parse_args()
    cop = {r["copy_idx"]: (r["chrom"], int(r["start"]), int(r["end"])) for r in csv.DictReader(open(a.copies), delimiter="\t")}
    def copy_at(chrom, s, e):
        return next((i for i, (c, cs, ce) in cop.items() if c == chrom and s < ce and e > cs), None)
    assign = {r["read_name"]: r for r in csv.DictReader(open(a.assign), delimiter="\t")}
    # per-read pairwise evidence from the dump: (n, kA, kB) for every candidate pair
    ev = {}
    for r in csv.DictReader(open(a.dump), delimiter="\t"):
        st = assign.get(r["read_name"])
        if st is None or st["origin_rejected"] != "0" or int(st["n_candidates"]) < 2:
            continue
        cands = r["candidates"].split(",")
        cols = [c.split(":") for c in r["columns"].split(",") if c] if r["columns"] else []
        pe = {}
        for i in range(len(cands)):
            for j in range(i + 1, len(cands)):
                n = ka = kb = 0
                for pos, o, al in cols:
                    x, y = al[i], al[j]
                    if x == "." or y == "." or x == y:
                        continue
                    n += 1
                    if o == x: ka += 1
                    elif o == y: kb += 1
                if n:
                    pe[(cands[i], cands[j])] = (n, ka, kb)
        ev[r["read_name"]] = (cands, pe, st["status"], st["catalog_copy_idx"])
    # chains from the BAM primaries; unique mappers per chain
    chain_of = {}; uniq_at = defaultdict(Counter)
    out = subprocess.run(["samtools", "view", "-F", "2308", a.bam], capture_output=True, text=True).stdout
    for ln in out.splitlines():
        f = ln.split("\t", 6)
        name, chrom, pos, cig = f[0], f[2], int(f[3]) - 1, f[5]
        if name in chain_of:
            continue
        key = (chrom,) + introns_of(pos, cig)
        chain_of[name] = key
        if name not in assign:
            end = pos + sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig) if op in "M=XDN")
            c = copy_at(chrom, pos, end)
            if c is not None:
                uniq_at[key][c] += 1
    groups = defaultdict(list)
    for name in ev:
        if name in chain_of and len(chain_of[name]) >= 2:   # chrom + >= 1 intron
            groups[chain_of[name]].append(name)
    abst = [n for n, (_, _, st, _) in ev.items() if st != "assigned"]
    in_multi = [n for n in abst if n in chain_of and len(groups.get(chain_of[n], [])) >= a.min_reads]
    print(f"AS-tied origin-pass contested reads with columns: {len(ev)}; abstaining: {len(abst)}")
    print(f"P1 abstaining reads in a chain group with >= {a.min_reads} contested reads: {len(in_multi)}/{len(abst)} = {100*len(in_multi)/max(1,len(abst)):.1f}%  (pass >=40)")
    # pooled verdict per group
    verdicts = {}; rows = []
    for key, names in groups.items():
        if len(names) < a.min_reads:
            continue
        C = sorted({c for n in names for c in ev[n][0]}, key=int)
        pool = defaultdict(lambda: [0, 0, 0])
        for n in names:
            for (A, B), (nn, ka, kb) in ev[n][1].items():
                p = pool[(A, B)]; p[0] += nn; p[1] += ka; p[2] += kb
        def stats(A, B):
            if (A, B) in pool:
                n, ka, kb = pool[(A, B)]
            elif (B, A) in pool:
                n, kb, ka = pool[(B, A)]
            else:
                return (0, 0, 0)
            return (n, ka, kb)
        def llr(A, B):
            n, ka, kb = stats(A, B); return LR * (ka - kb)
        def worst(A):
            return min((llr(A, B) for B in C if B != A), default=float("inf"))
        bk = max(C, key=lambda A: (worst(A), sum(llr(A, B) for B in C if B != A), -int(A)))
        thr = ALPHA / max(1, len(C) - 1)
        k0 = False; p_read = 0.0; margin = float("inf"); ndec = 10**9
        for B in C:
            if B == bk: continue
            n, ka, kb = stats(bk, B)
            if n == 0: k0 = True
            p = binom_tail(n, ka)
            p_read = max(p_read, p); margin = min(margin, llr(bk, B)); ndec = min(ndec, n)
        status = "tied" if k0 or len(C) < 2 else ("assigned" if p_read < thr and margin > 0 else "ambiguous")
        read_st = Counter(ev[n][2] for n in names)
        asg_copies = {ev[n][3] for n in names if ev[n][2] == "assigned"}
        verdicts[key] = (status, bk, margin, ndec, len(names), read_st, asg_copies, uniq_at.get(key, Counter()))
        rows.append((key[0], ";".join(f"{s}-{e}" for s, e in key[1:]), len(names), status, bk, f"{margin:.1f}", ndec, dict(read_st), ",".join(sorted(asg_copies)), ",".join(f"{c}:{k}" for c, k in uniq_at.get(key, Counter()).most_common())))
    st = Counter(v[0] for v in verdicts.values())
    print(f"P2 groups (>= {a.min_reads} contested reads): {len(verdicts)}; pooled verdicts: {dict(st)}; isoform-assigned {st['assigned']}/{len(verdicts)} = {100*st['assigned']/max(1,len(verdicts)):.1f}%  (pass 20..60)")
    # reads covered by an isoform assignment that were abstaining
    reads_res = sum(v[4] - v[5].get("assigned", 0) for v in verdicts.values() if v[0] == "assigned")
    print(f"   abstaining reads now covered by an isoform assignment: {reads_res}/{len(abst)} = {100*reads_res/max(1,len(abst)):.1f}%")
    # P3
    with_asg = [(k, v) for k, v in verdicts.items() if v[6]]
    contra = [(k, v) for k, v in with_asg if v[0] == "assigned" and v[1] not in v[6] or len(v[6]) > 1]
    print(f"P3 groups containing an assigned read: {len(with_asg)}; isoform copy != read copy (or reads disagree): {len(contra)}  (pass iff 0); isoform status among them: {Counter(v[0] for _, v in with_asg)}")
    # P5
    both = [(k, v) for k, v in verdicts.items() if v[0] == "assigned" and v[7]]
    agree = sum(1 for k, v in both if v[1] in v[7])
    print(f"P5 isoform-assigned groups with unique mappers: {len(both)}; certified copy has unique mappers too: {agree} ({100*agree/max(1,len(both)):.0f}%, report); certified copy is the unique-majority copy: {sum(1 for k, v in both if v[7].most_common(1)[0][0] == v[1])}")
    # P6
    if a.gtf and a.phantoms:
        # evidence-less transcripts = phantom copies + single-copy arbitrary; identify by chain from the GTF
        ex = defaultdict(list); tmeta = {}
        for line in open(a.gtf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            t = re.search(r'transcript_id "([^"]+)"', f[8])
            if not t: continue
            if f[2] == "exon": ex[t.group(1)].append((int(f[3]) - 1, int(f[4]), f[0]))
            elif f[2] == "transcript":
                ci = re.search(r'copy_index "([^"]+)"', f[8]); tmeta[t.group(1)] = ci.group(1) if ci else None
        tchain = {}
        for t, v in ex.items():
            v.sort(); tchain[t] = (v[0][2],) + tuple((p[1], q[0]) for p, q in zip(v, v[1:]))
        # evidence-less = chains whose reads (per §6hl) are abstaining only: recompute here = groups with no unique mappers and no assigned read
        evless = [k for k, v in verdicts.items() if not v[7] and not v[6]]
        cert = [k for k in evless if verdicts[k][0] == "assigned"]
        print(f"P6 pooled groups with NO unique mapper and NO assigned read (the evidence-less class): {len(evless)}; certified by pooling: {len(cert)}  (pass >=10); their copies: {Counter(verdicts[k][1] for k in cert).most_common()}")
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("chrom\tchain\tn_reads\tisoform_status\tisoform_copy\tmargin\tn_decisive_min\tread_statuses\tread_assigned_copies\tunique_mappers\n")
            for r in sorted(rows, key=lambda r: (-r[2])):
                fh.write("\t".join(map(str, r)) + "\n")
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
