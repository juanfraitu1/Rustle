#!/usr/bin/env python3
"""Resolve, EXACTLY, why each SD-predicted copy location has no catalog entry.

THE QUESTION. Segmental duplications predict 3,928 locations where the 627-family catalog has nothing.
Earlier work could not tell whether such a location NEVER BECAME A NODE or became a node that joined no
family, because `copies.tsv` lists family members only. Inspection of three well-supported examples
suggested the latter, but that is not a measurement.

THE DISCRIMINATOR. `RUSTLE_ER_EDGE_DUMP` emits one row PER REP with coordinates and a `degree` column
(`write_er_edge_dump`, denovo_pipeline.rs:5642-5662), from the SAME run and with no flag or definition
change. That splits the question three ways rather than two:

  (a)  absent from nodes.tsv        -> NEVER A NODE
  (b1) degree == 0                  -> a node with NO E_r edge at all
  (b2) degree > 0, not in copies.tsv -> a node that HAD an edge and was lost at the γ split
                                        or the >=2-distinct-loci gate

(b2) is a class no previous analysis could see, and it is the one that would indict the partition
rather than the relation.

⚠ REJECTED ALTERNATIVE, recorded so it is not retried: `--single-copy-baseline` returns before the
homology dispatch and builds families from the read-conflict graph (E_c), whose `colocated_families`
partitions by chromosome and drops spans > 5 Mb — so 1,681/2,019 = 83.3% of the catalog's own copies
would be emitted as "single-copy". Contaminated in both directions.
"""
import bisect, collections, csv, os, sys

R = "/mnt/linuxdisk/home/juanfraitu/o1_reps"
G = "/mnt/linuxdisk/home/juanfraitu/o1_gw"
SDBED = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/931c208e-8acb-4dd2-aacb-cf92d5ad051f/scratchpad/sd_absent_intervals.bed"


def find_nodes_tsv():
    import glob
    c = sorted(glob.glob(f"{R}/dump/*.nodes.tsv"), key=os.path.getsize, reverse=True)
    return c[0] if c else None


def main():
    nt = find_nodes_tsv()
    if not nt:
        sys.exit(f"no nodes.tsv under {R}/dump — has the run finished?")
    rows = list(csv.DictReader(open(nt), delimiter="\t"))
    print(f"rep table: {nt}\n  rows {len(rows)}   columns {list(rows[0].keys())}", flush=True)

    # ---- determinism gate: this run must reproduce the shipped catalog ----
    def sig(p):
        return sorted((r["chrom"], r["start"], r["end"])
                      for r in csv.DictReader(open(p), delimiter="\t"))
    a, b = f"{G}/ggo_gw.copies.tsv", f"{R}/ggo_reps.copies.tsv"
    if os.path.exists(b):
        same = sig(a) == sig(b)
        print(f"  DETERMINISM GATE — copies identical to the 2026-08-20 catalog: {same}"
              f"  ({len(sig(a))} vs {len(sig(b))} rows)")
        if not same:
            print("  ⚠ the run did NOT reproduce the shipped catalog; every number below describes "
                  "THIS run, not the shipped one.")
    print(flush=True)

    key = {k.lower(): k for k in rows[0]}
    C, S, E = key.get("chrom"), key.get("start"), key.get("end")
    DEG = key.get("degree")
    if not all((C, S, E, DEG)):
        sys.exit(f"missing expected columns; have {list(rows[0].keys())}")

    reps = collections.defaultdict(list)
    for r in rows:
        reps[r[C]].append((int(r[S]), int(r[E]), int(r[DEG])))
    for c in reps:
        reps[c].sort()
    starts = {c: [x[0] for x in v] for c, v in reps.items()}

    infam = set()
    for r in csv.DictReader(open(a), delimiter="\t"):
        infam.add((r["chrom"], int(r["start"]), int(r["end"])))

    def overlapping(c, s, e):
        v = reps.get(c)
        if not v:
            return []
        i = bisect.bisect_left(starts[c], e)
        return [v[j] for j in range(max(0, i - 800), i) if v[j][1] > s]

    cls = collections.Counter(); deg_when_node = []
    for line in open(SDBED):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        c, s, e = f[0], int(f[1]), int(f[2])
        ov = overlapping(c, s, e)
        if not ov:
            cls["(a) NEVER A NODE"] += 1
            continue
        best = max(ov, key=lambda x: x[2])
        deg_when_node.append(best[2])
        if best[2] == 0:
            cls["(b1) NODE, degree 0 — no E_r edge at all"] += 1
        elif (c, best[0], best[1]) in infam:
            cls["(c) node IS a catalog copy — SD interval merely abuts it"] += 1
        else:
            cls["(b2) NODE, degree>0 — lost at γ split / ≥2-distinct-loci gate"] += 1

    n = sum(cls.values())
    print("=== why does each SD-predicted location have no catalog entry? ===")
    for k, v in cls.most_common():
        print(f"  {k:<58} {v:>5}/{n}  = {v/n:.4f}")
    if deg_when_node:
        import statistics as st
        print(f"\n  where a node exists, median E_r degree: {st.median(deg_when_node):.0f}"
              f"   (degree 0: {sum(1 for d in deg_when_node if d==0)})")
    print(f"\n  rep set: {len(rows)} reps; catalog copies: {len(infam)}; "
          f"reps not in any family: {len(rows)-len(infam)}")


if __name__ == "__main__":
    main()
