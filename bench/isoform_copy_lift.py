#!/usr/bin/env python3
"""PREREG isoform_phantoms (62f34c64): lift each family transcript's intron chain through the copy-to-copy
alignments and find the SAME isoform emitted at several copies; classify each emitted copy by the evidence
its member reads carry (unique mapper / O2-assigned / abstaining only).

  python3 bench/isoform_copy_lift.py --gtf ours_final2.gtf --copies copies16.tsv --paf human_gspans.paf \
      --bam hsa16.bam --assign ours_final2.assignments.tsv --family MCL0 [--tol 5] [--out groups.tsv]
"""
import argparse, bisect, csv, re, subprocess
from collections import Counter, defaultdict


def introns_of(pos, cig):
    o, p = [], pos
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
        n = int(n)
        if op in "M=XD":
            p += n
        elif op == "N":
            o.append((p, p + n)); p += n
    return tuple(o)


class Lift:
    """Blocks of one PAF fragment: query position -> target position (both 0-based, span-relative)."""
    def __init__(self, qs, qe, ts, strand, cg):
        self.blocks = []  # (q_lo, q_hi, t_lo, sign)
        q, t = (qe if strand == "-" else qs), ts
        for n, op in re.findall(r"(\d+)([=XIDM])", cg):
            n = int(n)
            if op in "=XM":
                if strand == "-":
                    self.blocks.append((q - n, q, t, -1)); q -= n
                else:
                    self.blocks.append((q, q + n, t, +1)); q += n
                t += n
            elif op == "I":
                q = q - n if strand == "-" else q + n
            else:
                t += n
        self.blocks.sort()
        self.lo = [b[0] for b in self.blocks]

    def map(self, qpos):
        i = bisect.bisect_right(self.lo, qpos) - 1
        if i >= 0:
            q_lo, q_hi, t_lo, s = self.blocks[i]
            if qpos < q_hi:
                return (t_lo + (qpos - q_lo) if s > 0 else t_lo + (q_hi - 1 - qpos)), 0
        # not inside a block: nearest block edge
        best = None
        for j in (i, i + 1):
            if 0 <= j < len(self.blocks):
                q_lo, q_hi, t_lo, s = self.blocks[j]
                for qq in (q_lo, q_hi - 1):
                    d = abs(qq - qpos)
                    if best is None or d < best[1]:
                        tt = t_lo + (qq - q_lo) if s > 0 else t_lo + (q_hi - 1 - qq)
                        best = (tt + (qpos - qq) * (s), d)
        return (best[0], best[1]) if best else (None, None)


def build_lifts(paf_path, span2idx):
    """Copy-to-copy lifts in both directions (all fragments) and pairwise identity, from a minimap2 --eqx PAF
    over the genomic unit spans (query/target names = "chrom:start1-end")."""
    lifts = defaultdict(list); acc = defaultdict(lambda: [0, 0])
    for l in open(paf_path):
        f = l.rstrip().split("\t")
        if f[0] not in span2idx or f[5] not in span2idx or f[0] == f[5]:
            continue
        qi, ti = span2idx[f[0]], span2idx[f[5]]
        cg = next(x for x in f[12:] if x.startswith("cg:Z:"))[5:]
        qs, qe, ts = int(f[2]), int(f[3]), int(f[7])
        fwd = Lift(qs, qe, ts, f[4], cg)
        lifts[(qi, ti)].append(fwd)
        inv = Lift.__new__(Lift)
        inv.blocks = [(tl, tl + (qh - ql), (ql if s > 0 else ql), s) for (ql, qh, tl, s) in fwd.blocks]
        # for s < 0 the forward block maps q in [ql, qh) to t = tl + (qh-1-q); the inverse maps t in [tl, th) to
        # q = qh-1-(t-tl) = (ql) + (th-1-t) with th = tl+(qh-ql): store q_lo slot = ql so map() gives ql + (th-1-t)
        inv.blocks.sort(); inv.lo = [b[0] for b in inv.blocks]
        lifts[(ti, qi)].append(inv)
        k = tuple(sorted((qi, ti)))
        acc[k][0] += int(f[10])
        acc[k][1] += sum(int(n) for n, op in re.findall(r"(\d+)([=XID])", cg) if op == "X")
    ident = {k: 1 - x / al for k, (al, x) in acc.items()}
    return lifts, ident


def make_lift_pos(lifts, cop):
    def lift_pos(gpos, A, B):
        """genomic position in copy A -> genomic position in copy B (best fragment), with distance-to-block."""
        ca, sa, ea = cop[A]; cb, sb, eb = cop[B]
        rel = gpos - sa
        best = None
        for L in lifts.get((A, B), []):
            t, d = L.map(rel)
            if t is not None and (best is None or d < best[1]):
                best = (t, d)
        return (best[0] + sb, best[1]) if best else (None, None)
    return lift_pos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True); ap.add_argument("--copies", required=True); ap.add_argument("--paf", required=True)
    ap.add_argument("--bam", required=True); ap.add_argument("--assign", required=True); ap.add_argument("--family", default="MCL0")
    ap.add_argument("--tol", type=int, default=5); ap.add_argument("--out")
    a = ap.parse_args()
    cop = {}
    for r in csv.DictReader(open(a.copies), delimiter="\t"):
        cop[r["copy_idx"]] = (r["chrom"], int(r["start"]), int(r["end"]))
    span2idx = {f"{c}:{s+1}-{e}": i for i, (c, s, e) in cop.items()}
    lifts, ident = build_lifts(a.paf, span2idx)
    lift_pos = make_lift_pos(lifts, cop)

    # transcripts of the family
    tx = {}  # tid -> dict(copy, strand, chain, chrom)
    ex = defaultdict(list); meta = {}
    for line in open(a.gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        if a.family and f[2] == "transcript" and f'family_id "{a.family}"' not in f[8]:
            continue
        tm = re.search(r'transcript_id[ =]"?([^";]*)"?', f[8])
        if not tm:
            continue
        t = tm.group(1)
        if f[2] == "transcript":
            ci = re.search(r'copy_index "([^"]+)"', f[8]); st = re.search(r'copy_status "([^"]+)"', f[8])
            # B4 (docs/OPEN_ITEMS_2026-09-09.md, register/6hn P1): the binary's OWN evidence vocabulary is
            # `placed_by` (unique_mapper / assigned_read / aligner_primary), not `copy_status` (B3's old
            # vocabulary). A transcript lifted here from another copy's certified read (`placed_by
            # "assigned_read"`) is legitimate evidence this script's independent re-derivation didn't know
            # about, and called a phantom instead.
            pb = re.search(r'placed_by "([^"]+)"', f[8])
            meta[t] = (ci.group(1) if ci else None, f[6], f[0], st.group(1) if st else None, pb.group(1) if pb else None)
        elif f[2] == "exon":
            ex[t].append((int(f[3]) - 1, int(f[4])))
            meta.setdefault(t, (None, f[6], f[0], None, None))   # tools without transcript rows
    def copy_by_overlap(chrom, s, e):
        best, bo = None, 0
        for i, (c, cs, ce) in cop.items():
            if c != chrom: continue
            ov = min(e, ce) - max(s, cs)
            if ov > bo: best, bo = i, ov
        return best
    for t, (ci, strand, chrom, st, pb) in meta.items():
        if t not in ex:
            continue
        if ci is None:   # no copy_index attribute (other tools): max-overlap copy, None when outside every copy
            v = sorted(ex[t]); ci = copy_by_overlap(chrom, v[0][0], v[-1][1])
            if ci is None:
                continue
        v = sorted(ex[t]); chain = tuple((p[1], q[0]) for p, q in zip(v, v[1:]))
        tx[t] = dict(copy=ci, strand=strand, chrom=chrom, chain=chain, status=st, placed_by=pb)
    multi = {t: d for t, d in tx.items() if len(d["chain"]) >= 2}
    print(f"family {a.family}: {len(tx)} transcripts with a copy, {len(multi)} with >= 2 introns (the denominator)")

    # same-isoform pairs across copies by lifting
    by_n = defaultdict(list)
    for t, d in multi.items():
        by_n[(len(d["chain"]), d["strand"])].append(t)
    parent = {t: t for t in multi}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    exact_pairs = same_pairs = 0
    for key, ts in by_n.items():
        for i in range(len(ts)):
            for j in range(i + 1, len(ts)):
                t, u = ts[i], ts[j]
                A, B = multi[t]["copy"], multi[u]["copy"]
                if A == B or (A, B) not in lifts:
                    continue
                ok = True; exact = True
                for (x0, x1), (y0, y1) in zip(multi[t]["chain"], multi[u]["chain"]):
                    lx0, d0 = lift_pos(x0, A, B); lx1, d1 = lift_pos(x1, A, B)
                    if lx0 is None or lx1 is None or abs(lx0 - y0) > a.tol or abs(lx1 - y1) > a.tol:
                        ok = False; break
                    if lx0 != y0 or lx1 != y1:
                        exact = False
                if ok:
                    same_pairs += 1; exact_pairs += exact
                    parent[find(t)] = find(u)
    groups = defaultdict(list)
    for t in multi:
        groups[find(t)].append(t)
    print(f"same-isoform transcript pairs across copies: {same_pairs} (exact after lift: {exact_pairs}); isoform groups: {len(groups)}")

    # member reads and their evidence
    assign = {r["read_name"]: r for r in csv.DictReader(open(a.assign), delimiter="\t")}
    chain_to_tx = defaultdict(list)
    for t, d in multi.items():
        chain_to_tx[(d["chrom"],) + d["chain"]].append(t)
    def copy_at(chrom, s, e):
        return next((i for i, (c, cs, ce) in cop.items() if c == chrom and s < ce and e > cs), None)
    members = defaultdict(list)  # tid -> [(read, evidence_copy or None, kind)]
    out = subprocess.run(["samtools", "view", "-F", "2308", a.bam], capture_output=True, text=True).stdout
    seen = set()
    for ln in out.splitlines():
        f = ln.split("\t", 6)
        name, chrom, pos, cig = f[0], f[2], int(f[3]) - 1, f[5]
        if name in seen:
            continue
        seen.add(name)
        ich = introns_of(pos, cig)
        for t in chain_to_tx.get((chrom,) + ich, ()):
            end = pos + sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig) if op in "M=XDN")
            pc = copy_at(chrom, pos, end)
            r = assign.get(name)
            if r is None:
                members[t].append((name, pc, "unique"))
            elif r["status"] == "assigned" and r["origin_rejected"] == "0":
                members[t].append((name, r["catalog_copy_idx"], "assigned"))
            else:
                members[t].append((name, None, "abstain"))
    # classify emitted copies per group
    n_multi_copy = 0; groups_multi = []; phantom_groups = 0; shared_groups = 0
    phantom_tx = []; rows = []
    for g, ts in groups.items():
        copies = defaultdict(list)
        for t in ts:
            copies[multi[t]["copy"]].append(t)
        if len(copies) >= 2:
            n_multi_copy += len(ts)
            ev = {}
            for c, tl in copies.items():
                kinds = Counter(k for t in tl for (_, ec, k) in members[t] if (k != "assigned" or ec == c))
                # evidence at c: a unique mapper whose primary is at c, or a read assigned to c
                uniq_at = sum(1 for t in tl for (_, ec, k) in members[t] if k == "unique" and ec == c)
                asg_at = sum(1 for t in tl for (_, ec, k) in members[t] if k == "assigned" and ec == c)
                abst = sum(1 for t in tl for (_, ec, k) in members[t] if k == "abstain")
                # B4: a lifted placement (placed_by "assigned_read") IS evidence — the certified read lives
                # at another copy in this same evidenced group, and the binary deliberately placed its
                # isoform here too. Only "aligner_primary" (no certificate anywhere) is a true phantom.
                lifted = any(multi[t].get("placed_by") == "assigned_read" for t in tl)
                ev[c] = (uniq_at, asg_at, abst, [multi[t]["status"] for t in tl], lifted)
            with_ev = [c for c, (u, s, _, _, lifted) in ev.items() if u + s > 0 or lifted]
            no_ev = [c for c, (u, s, _, _, lifted) in ev.items() if u + s == 0 and not lifted]
            groups_multi.append((g, copies, ev, with_ev, no_ev))
            if no_ev: phantom_groups += 1
            if len(with_ev) >= 2: shared_groups += 1
            for c in no_ev:
                sib = max((ident.get(tuple(sorted((c, c2))), 0.0) for c2 in copies if c2 != c), default=0.0)
                phantom_tx.append((c, sib, ev[c][3]))
            rows.append((g, len(ts), ",".join(sorted(copies)), ",".join(sorted(with_ev)), ",".join(sorted(no_ev)),
                         ";".join(f"{c}:u{u}/a{s}/x{x}" for c, (u, s, x, _, _) in sorted(ev.items()))))
    print(f"\nP1 multi-intron isoforms emitted at >= 2 copies: {n_multi_copy}/{len(multi)} = {100*n_multi_copy/max(1,len(multi)):.1f}%  (pass 10..40)")
    print(f"P2 isoform groups at >= 2 copies: {len(groups_multi)}; with >= 1 copy lacking evidence (phantom): {phantom_groups} = {100*phantom_groups/max(1,len(groups_multi)):.0f}%  (pass >=50)")
    print(f"P3 genuinely shared (evidence at >= 2 copies): {shared_groups}  (pass >=5)")
    if phantom_tx:
        k = sum(1 for _, sib, _ in phantom_tx if sib >= 0.985)
        print(f"P4 phantom copies: {len(phantom_tx)}; with a sibling >= 0.985 in the group: {k} = {100*k/len(phantom_tx):.0f}%  (pass >=80); statuses of phantom transcripts: {Counter(s for _, _, sl in phantom_tx for s in sl)}")
        print(f"   phantom copies by copy: {Counter(c for c, _, _ in phantom_tx).most_common()}")
    asg_ph = sum(1 for _, _, sl in phantom_tx if "assigned" in sl)
    print(f"P5 `assigned` transcripts that are phantoms: {asg_ph}  (pass iff 0)")
    # extras
    single = [g for g, ts in groups.items() if len({multi[t]['copy'] for t in ts}) == 1]
    arb = 0
    for g in single:
        t = groups[g][0]
        kinds = Counter(k for (_, _, k) in members[t])
        if kinds and kinds.get("unique", 0) + kinds.get("assigned", 0) == 0:
            arb += 1
    print(f"extra: single-copy isoform groups {len(single)}, of which placed on abstaining reads only (arbitrary address): {arb}")
    print(f"extra: member-read kinds over multi-intron transcripts: {Counter(k for t in multi for (_, _, k) in members[t])}")
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("group\tn_tx\tcopies\tcopies_with_evidence\tcopies_without\tper_copy(u=unique,a=assigned,x=abstain)\n")
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
