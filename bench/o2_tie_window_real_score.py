# REPORT (LEDGER_O2): D vs T O2 on real reads, truth = labels.tsv (KB3781 haplotypes)
import csv, collections, json, pysam, sys
for p in ("p5", "p10"):
    W = f"/mnt/linuxdisk/tmp/rna_allele/{p}"
    lab = {r["read"]: r for r in csv.DictReader(open(f"{W}/labels.tsv"), delimiter="\t")}
    span = {r["tid"]: (r["chrom"], int(r["start"]), int(r["end"])) for r in csv.DictReader(open(f"{W}/c.copies.tsv"), delimiter="\t")}
    fam_copies = collections.defaultdict(list)
    for t, (c, s, e) in span.items(): fam_copies[t.split(":")[0]].append(t)
    rows = {a: list(csv.DictReader(open(f"{p}_{a}.o2.tsv"), delimiter="\t")) for a in ("D", "T")}
    need = {r["read"] for a in rows for r in rows[a]}
    prim = {}
    for rec in pysam.AlignmentFile(f"{W}/R0.bam").fetch(until_eof=True):
        if rec.query_name in need and not (rec.is_secondary or rec.is_supplementary or rec.is_unmapped):
            prim[rec.query_name] = (rec.reference_name, rec.reference_start, rec.reference_end)
    def prim_copy(n, f):
        if n not in prim: return None
        c, s, e = prim[n]
        hits = [t for t in fam_copies[f] if span[t][0] == c and s < span[t][2] and e > span[t][1]]
        return hits[0] if len(hits) == 1 else ("multi" if hits else "outside")
    print(f"===== {p}")
    asg_by = {}
    for a in ("D", "T"):
        c = collections.Counter(); asg = {}
        for r in rows[a]:
            n, f, k = r["read"], r["family"], r["assigned"]
            c["tied rows"] += 1
            L = lab.get(n); truth = L["copy"] if L and L["family"] == f else None
            c["labelled"] += truth is not None
            pc = prim_copy(n, f)
            if truth: c["primary right (labelled)"] += pc == truth
            if k in ("ABSTAIN", "HOPELESS"): c[k] += 1; continue
            asg[n] = k; c["assigned"] += 1
            if truth:
                c["assigned labelled"] += 1; c["right"] += k == truth; c["wrong"] += k != truth
                c["  fixed (primary wrong)"] += k == truth and pc != truth
                c["  broken (primary right)"] += k != truth and pc == truth
        asg_by[a] = asg
        print(f"  {a}: " + "; ".join(f"{k} {v}" for k, v in c.items()))
    sD, sT = {(r["family"], r["read"]) for r in rows["D"]}, {(r["family"], r["read"]) for r in rows["T"]}
    print(f"  tied rows only in D {len(sD - sT)}, only in T {len(sT - sD)}")
    diff = collections.Counter()
    for n in set(asg_by["D"]) | set(asg_by["T"]):
        if asg_by["D"].get(n) != asg_by["T"].get(n):
            L = lab.get(n)
            diff[(asg_by["D"].get(n, "-") == (L or {}).get("copy"), asg_by["T"].get(n, "-") == (L or {}).get("copy"), L is not None)] += 1
    print("  assignment changes (D right, T right, labelled):", dict(diff))
