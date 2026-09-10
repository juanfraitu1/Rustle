#!/usr/bin/env python3
"""PREREG gtf_copy_set (95409846): rewrite the family transcripts of our GTF so that every emitted address is
backed by evidence (a unique mapper or a read-level certified read), phantoms are dropped, certified reads
whose transcript sits at another copy get a LIFTED placement, and isoforms with no evidence are emitted ONCE
with a copy SET. Non-family transcripts and single-intron/unspliced family transcripts pass through.

  python3 bench/gtf_copy_set.py --gtf ours_final2.gtf --copies copies16.tsv --paf human_gspans.paf --bam hsa16.bam \
      --assign ours_final2.assignments.tsv --dump ours_final2_dump.star_reads.tsv [--pooled isoforms_pooled.tsv] \
      --family MCL0 --out ours_copyset.gtf
"""
import argparse, csv, os, re, subprocess, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from isoform_copy_lift import Lift, build_lifts, make_lift_pos, introns_of


def main():
    ap = argparse.ArgumentParser()
    for k in ("--gtf", "--copies", "--paf", "--bam", "--assign", "--dump", "--out"):
        ap.add_argument(k, required=True)
    ap.add_argument("--pooled"); ap.add_argument("--family", default="MCL0"); ap.add_argument("--tol", type=int, default=5)
    a = ap.parse_args()
    cop = {}; cstrand = {}
    for r in csv.DictReader(open(a.copies), delimiter="\t"):
        cop[r["copy_idx"]] = (r["chrom"], int(r["start"]), int(r["end"])); cstrand[r["copy_idx"]] = r["strand"]
    span2idx = {f"{c}:{s+1}-{e}": i for i, (c, s, e) in cop.items()}
    lifts, ident = build_lifts(a.paf, span2idx)
    lift_pos = make_lift_pos(lifts, cop)
    def copy_at(chrom, s, e):
        return next((i for i, (c, cs, ce) in cop.items() if c == chrom and s < ce and e > cs), None)

    # ---- GTF
    lines = open(a.gtf).read().split("\n")
    tx_lines = defaultdict(list)   # tid -> line indices
    meta = {}; ex = defaultdict(list)
    for i, line in enumerate(lines):
        f = line.split("\t")
        if len(f) < 9: continue
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        if not t: continue
        t = t.group(1); tx_lines[t].append(i)
        if f[2] == "transcript":
            fam = re.search(r'family_id "([^"]+)"', f[8]); ci = re.search(r'copy_index "([^"]+)"', f[8])
            meta[t] = dict(fam=fam.group(1) if fam else None, copy=ci.group(1) if ci else None, strand=f[6], chrom=f[0], attrs=f[8], gene=re.search(r'gene_id "([^"]+)"', f[8]).group(1), line=i)
        elif f[2] == "exon":
            ex[t].append((int(f[3]) - 1, int(f[4])))
    fam_tx = {t: m for t, m in meta.items() if m["fam"] == a.family and m["copy"] is not None and t in ex}
    for t in fam_tx:
        v = sorted(ex[t]); fam_tx[t]["exons"] = v; fam_tx[t]["chain"] = tuple((p[1], q[0]) for p, q in zip(v, v[1:]))
    multi = {t: m for t, m in fam_tx.items() if len(m["chain"]) >= 2}

    # ---- groups by lift (as in isoform_copy_lift.py)
    by_n = defaultdict(list)
    for t, m in multi.items(): by_n[(len(m["chain"]), m["strand"])].append(t)
    parent = {t: t for t in multi}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for key, ts in by_n.items():
        for i in range(len(ts)):
            for j in range(i + 1, len(ts)):
                t, u = ts[i], ts[j]; A, B = multi[t]["copy"], multi[u]["copy"]
                if A == B or (A, B) not in lifts: continue
                ok = True
                for (x0, x1), (y0, y1) in zip(multi[t]["chain"], multi[u]["chain"]):
                    lx0, d0 = lift_pos(x0, A, B); lx1, d1 = lift_pos(x1, A, B)
                    if lx0 is None or lx1 is None or abs(lx0 - y0) > a.tol or abs(lx1 - y1) > a.tol: ok = False; break
                if ok: parent[find(t)] = find(u)
    groups = defaultdict(list)
    for t in multi: groups[find(t)].append(t)

    # ---- reads: chains, evidence, candidates
    assign = {r["read_name"]: r for r in csv.DictReader(open(a.assign), delimiter="\t")}
    # the copy SET of an abstaining read = the catalog units at its AS-TIED placements (the gate's tie set), plus
    # "outside" when a tied placement overlaps no unit (tie_outside_catalog) — not the read-star candidate list
    # (every copy with a -p 0.3 hit), which is far wider than the tie
    import pysam
    tie_units = defaultdict(set); best_as = {}
    abst_names = {n for n, r in assign.items() if r["status"] != "assigned"}
    with pysam.AlignmentFile(a.bam) as fh:
        for al in fh:
            if al.query_name not in abst_names or al.is_unmapped or al.is_supplementary: continue
            AS = al.get_tag("AS"); u = copy_at(al.reference_name, al.reference_start, al.reference_end)
            b = best_as.get(al.query_name)
            if b is None or AS > b: best_as[al.query_name] = AS; tie_units[al.query_name] = {u if u is not None else "outside"}
            elif AS == b: tie_units[al.query_name].add(u if u is not None else "outside")
    cands = {n: sorted(str(x) for x in v) for n, v in tie_units.items()}
    for n in abst_names:
        if assign[n].get("tie_outside_catalog") == "1": cands.setdefault(n, []).append("outside")
    chain_to_tx = defaultdict(list)
    for t, m in multi.items(): chain_to_tx[(m["chrom"],) + m["chain"]].append(t)
    members = defaultdict(list)   # tid -> (name, kind, copy or candidate set)
    out = subprocess.run(["samtools", "view", "-F", "2308", a.bam], capture_output=True, text=True).stdout
    seen = set()
    for ln in out.splitlines():
        f = ln.split("\t", 6); name, chrom, pos, cig = f[0], f[2], int(f[3]) - 1, f[5]
        if name in seen: continue
        seen.add(name)
        ts = chain_to_tx.get((chrom,) + introns_of(pos, cig))
        if not ts: continue
        end = pos + sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig) if op in "M=XDN")
        pc = copy_at(chrom, pos, end); r = assign.get(name)
        if r is None:
            kind, val = ("unique", pc) if pc is not None else ("outside", None)
        elif r["status"] == "assigned" and r["origin_rejected"] == "0":
            kind, val = "assigned", r["catalog_copy_idx"]
        else:
            kind, val = "abstain", set(cands.get(name, []))
        for t in ts: members[t].append((name, kind, val))
    pooled = {}
    if a.pooled:
        for r in csv.DictReader(open(a.pooled), delimiter="\t"):
            key = (r["chrom"],) + tuple(tuple(int(x) for x in s.split("-")) for s in r["chain"].split(";") if s)
            pooled[key] = (r["isoform_status"], r["isoform_copy"], r["margin"])

    # ---- rewrite
    drop = set(); add = []; attr_extra = {}; stats = Counter()
    def fmt(c): return ",".join(f"{k}:{v}" for k, v in sorted(c.items(), key=lambda kv: int(kv[0])))
    for g, ts in groups.items():
        uniq = Counter(); asg = Counter(); und = set(); n_abst = 0
        for t in ts:
            for name, kind, val in members[t]:
                if kind == "unique": uniq[val] += 1
                elif kind == "assigned": asg[val] += 1
                elif kind == "abstain": und |= set(val); n_abst += 1
        E = set(uniq) | set(asg)
        key0 = (multi[ts[0]]["chrom"],) + multi[ts[0]]["chain"]
        pl = pooled.get(key0)
        base_attr = f' evidence_unique "{fmt(uniq)}"; evidence_assigned "{fmt(asg)}"; reads_undecided "{n_abst}";'
        if pl: base_attr += f' pooled_status "{pl[0]}"; pooled_copy "{pl[1]}"; pooled_margin "{pl[2]}";'
        have = {multi[t]["copy"]: t for t in ts}
        if E:
            for c, t in have.items():
                if c in E:
                    attr_extra[t] = base_attr + f' copies "{",".join(sorted(E, key=int))}"; copies_undecided "{",".join(sorted(und - E, key=lambda x: (x == "outside", int(x) if x != "outside" else 0)))}"; placed_by "{"unique_mapper" if uniq.get(c) else "assigned_read"}";'
                    stats["kept"] += 1
                else:
                    drop.add(t); stats["dropped_phantom"] += 1
            # certified copies with no transcript in the group: lift the best-supported transcript there
            rep = max(ts, key=lambda t: len(members[t]))
            for c in E - set(have):
                A = multi[rep]["copy"]
                if (A, c) not in lifts: stats["lift_failed_noalign"] += 1; continue
                lifted = []; ok = True
                for (s0, e0) in multi[rep]["exons"]:
                    ls, ds = lift_pos(s0, A, c); le, de = lift_pos(e0, A, c)
                    if ls is None or le is None or ds > a.tol or de > a.tol: ok = False; break
                    lifted.append((min(ls, le), max(ls, le)))
                lifted.sort()
                if not ok or any(p[1] > q[0] for p, q in zip(lifted, lifted[1:])): stats["lift_failed"] += 1; continue
                add.append((rep, c, lifted, base_attr + f' copies "{",".join(sorted(E, key=int))}"; copies_undecided "{",".join(sorted(und - E, key=lambda x: (x == "outside", int(x) if x != "outside" else 0)))}"; placed_by "assigned_read"; lifted_from "{rep}";'))
                stats["lifted_added"] += 1
        else:
            # no evidence anywhere: keep ONE transcript (the best-supported), drop the rest, copy SET
            rep = max(ts, key=lambda t: len(members[t]))
            for t in ts:
                if t != rep: drop.add(t); stats["dropped_duplicate_undecided"] += 1
            attr_extra[rep] = base_attr + f' copies "{",".join(sorted(und, key=lambda x: (x == "outside", int(x) if x != "outside" else 0)))}"; placed_by "aligner_primary"; copy_status_final "undecidable";'
            stats["undecided_single"] += 1
            if len(und) >= 2: stats["undecided_with_set>=2"] += 1
            if "outside" in und: stats["undecided_with_outside"] += 1
    # emit
    with open(a.out, "w") as fh:
        for t, idxs in tx_lines.items():
            if t in drop: continue
            for i in idxs:
                line = lines[i]
                if t in attr_extra and "\ttranscript\t" in line:
                    line = line.rstrip(";") + ";" + attr_extra[t]
                fh.write(line + "\n")
        for rep, c, lifted, attr in add:
            m = multi[rep]; chrom = cop[c][0]; strand = cstrand[c]; tid = f"{rep}_lift{c}"; gene = f"{m['gene']}_copy{c}"
            base = re.sub(r'copy_index "[^"]*"', f'copy_index "{c}"', m["attrs"]).replace(f'transcript_id "{rep}"', f'transcript_id "{tid}"').replace(f'gene_id "{m["gene"]}"', f'gene_id "{gene}"')
            fh.write("\t".join([chrom, "rustle", "transcript", str(lifted[0][0] + 1), str(lifted[-1][1]), ".", strand, ".", base.rstrip(";") + ";" + attr]) + "\n")
            order = lifted if strand == "+" else lifted[::-1]
            for k, (s0, e0) in enumerate(order, 1):
                fh.write("\t".join([chrom, "rustle", "exon", str(s0 + 1), str(e0), ".", strand, ".", f'gene_id "{gene}"; transcript_id "{tid}"; exon_number "{k}";']) + "\n")
        # pass-through lines without a transcript_id (headers)
        for i, line in enumerate(lines):
            if line and not re.search(r'transcript_id "', line): fh.write(line + "\n")
    n_fam = sum(1 for t in fam_tx if t not in drop) + len(add)
    print(f"groups {len(groups)}; {dict(stats)}; family transcripts {len(fam_tx)} -> {n_fam}; wrote {a.out}")


if __name__ == "__main__":
    main()
