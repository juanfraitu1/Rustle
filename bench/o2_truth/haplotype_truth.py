#!/usr/bin/env python3
"""Arm 4 of docs/archive/2026-10/merged/PREREG_o2_anchor_2026-10-13.md (Amendment 4): truth from the animal's own haplotype assemblies.

    haplotype_truth.py RUN_PREFIX CATALOG_PREFIX BAM WINDOWS.bed MAT.mmi PAT.mmi OUTDIR

A KB3781 read must match its true locus on mat or pat up to sequencing error. Reads (original orientation) and every
involved catalog copy's spliced sequence are aligned to both haplotypes (minimap2 -ax splice:hq -uf --eqx -N 10). A read's
truth = the catalog copies whose haplotype loci overlap one of its NEAR hits (NM <= best + 2), provided its best hit has
NM <= max(2, 0.003 x aligned bases). Scored per (read, family) row: O2's call is right iff its copy is in the truth set and
the truth set is one locus (its copies overlap in primary coordinates). Also scores the aligner's primary.
"""
import collections, csv, hashlib, json, os, re, subprocess, sys

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
SAMPLE = 3000


def parse_sam(lines):
    """{qname: [(chrom, start0, end, NM, aligned bases)]} over mapped records."""
    out = collections.defaultdict(list)
    for l in lines:
        if not l or l.startswith("@"):
            continue
        f = l.rstrip("\n").split("\t")
        if int(f[1]) & 4:
            continue
        ops = CIG.findall(f[5])
        end = int(f[3]) - 1 + sum(int(n) for n, op in ops if op in "MDN=X")
        aligned = sum(int(n) for n, op in ops if op in "M=X")
        nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), 10 ** 6)
        out[f[0]].append((f[2], int(f[3]) - 1, end, nm, aligned))
    return out


def near_hits(hits):
    """hits = [(hap, chrom, s, e, NM, aligned)] -> the hits within NM + 2 of the best, or None if the best is not a
    near-perfect match (NM <= max(2, 0.003 x aligned))."""
    if not hits:
        return None
    best = min(hits, key=lambda h: (h[4], -h[5]))
    if best[4] > max(2, 0.003 * best[5]):
        return None
    return [h for h in hits if h[4] <= best[4] + 2]


def copy_loci(hits):
    """A copy's loci: per haplotype, its hits within NM + 2 of its best there -> {(hap, chrom, s, e)}."""
    best = {}
    for h in hits:
        best[h[0]] = min(best.get(h[0], 10 ** 9), h[4])
    return {(h[0], h[1], h[2], h[3]) for h in hits if h[4] <= best[h[0]] + 2}


def truth_copies(near, loci):
    """near = a read's near hits [(hap, chrom, s, e, NM, aligned)]; loci = {copy: {(hap, chrom, s, e)}}."""
    out = set()
    for hap, ch, s, e, *_ in near:
        for k, ls in loci.items():
            if any(lh == hap and lc == ch and ls0 < e and s < le for lh, lc, ls0, le in ls):
                out.add(k)
    return out


def overlaps(a, b):
    return a[0] == b[0] and a[1] < b[2] and b[1] < a[2]


def one_locus(spans):
    return bool(spans) and all(overlaps(a, b) for i, a in enumerate(spans) for b in spans[i + 1:])


def tsv(path):
    csv.field_size_limit(sys.maxsize)
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main(prefix, catalog, bam, windows, mat, pat, out):
    os.makedirs(out, exist_ok=True)
    span = {(r["family_id"], r["copy_idx"]): (r["chrom"], int(r["start"]), int(r["end"]))
            for r in tsv(catalog + ".fam.copies.tsv")}
    rows = [r for r in tsv(prefix + ".assignments.tsv") if r["contested"] == "1"]
    assigned = {r["read_name"] for r in rows if r["status"] == "assigned"}
    names = sorted({r["read_name"] for r in rows}, key=lambda n: hashlib.md5(f"7:{n}".encode()).hexdigest())
    sample = set(names[:SAMPLE])
    want = assigned | sample
    fams = sorted({r["family_id"] for r in rows if r["read_name"] in want})
    with open(f"{out}/names.txt", "w") as fh:
        fh.write("\n".join(sorted(want)) + "\n")
    # primary records: sequence (original orientation) + primary span
    prim = {}
    res = subprocess.run(["samtools", "view", "-F", "2308", "-L", windows, "-N", f"{out}/names.txt", bam],
                         capture_output=True, text=True, check=True).stdout
    for l in res.splitlines():
        f = l.split("\t")
        if f[0] in prim or f[9] == "*":
            continue
        seq = f[9].translate(COMP)[::-1] if int(f[1]) & 16 else f[9]
        end = int(f[3]) - 1 + sum(int(n) for n, op in CIG.findall(f[5]) if op in "MDN=X")
        prim[f[0]] = (seq, (f[2], int(f[3]) - 1, end))
    # query = reads + the involved families' copies (spliced, `FAM|idx|...` headers)
    with open(f"{out}/query.fa", "w") as fh:
        fh.writelines(f">R|{q}\n{s}\n" for q, (s, _) in prim.items())
        on = False
        for l in open(catalog + ".fam.copies.fa"):
            if l.startswith(">"):
                h = l[1:].split("|")
                on = h[0] in set(fams)
                if on:
                    fh.write(f">C|{h[0]}|{h[1]}\n")
                    continue
            if on:
                fh.write(l)
    hits = collections.defaultdict(list)
    for hap, idx in (("mat", mat), ("pat", pat)):
        sam = subprocess.run(["minimap2", "-ax", "splice:hq", "-uf", "--eqx", "-N", "10", "-t", "4", idx,
                              f"{out}/query.fa"], capture_output=True, text=True, check=True).stdout
        for q, hs in parse_sam(sam.splitlines()).items():
            hits[q].extend((hap,) + h for h in hs)
    loci = {(k.split("|")[1], k.split("|")[2]): copy_loci(v) for k, v in hits.items() if k.startswith("C|")}
    c = collections.Counter(rows_considered=0)
    ex = {"wrong": [], "right": [], "outside": []}
    for r in rows:
        q, fam = r["read_name"], r["family_id"]
        if q not in want:
            continue
        c["rows_considered"] += 1
        if q not in prim:
            c["no_primary"] += 1
            continue
        near = near_hits(hits.get("R|" + q, []))
        if near is None:
            c["no_near_perfect_hit"] += 1
            continue
        T = truth_copies(near, loci)
        is_assigned = r["status"] == "assigned"
        tag = "assigned" if is_assigned else "abstained"
        if not T:
            c[f"truth_outside_catalog|{tag}"] += 1
            if is_assigned and len(ex["outside"]) < 20:
                ex["outside"].append({"read": q, "family": fam, "o2_copy": r["catalog_copy_idx"],
                                      "best": min(near, key=lambda h: h[4])})
            continue
        if not one_locus([span[k] for k in T]):
            c[f"truth_ambiguous|{tag}"] += 1
            continue
        in_fam = any(k[0] == fam for k in T)
        c[f"one_locus|{tag}"] += 1
        p = prim[q][1]
        al_right = any(overlaps(p, span[k]) for k in T)
        if is_assigned:
            ok = (fam, r["catalog_copy_idx"]) in T
            c["H1_assigned"] += 1
            c["H1_right" if ok else "H1_wrong"] += 1
            c["H2_aligner_right"] += al_right
            rec = {"read": q, "family": fam, "o2_copy": r["catalog_copy_idx"], "truth": sorted(T), "aligner_right": al_right,
                   "near": near[:4], "margin": r["margin"], "p": r["p_value"]}
            if not ok:
                ex["wrong"].append(rec)
            elif len(ex["right"]) < 20:
                ex["right"].append(rec)
        if q in sample:
            c["H3_one_locus"] += 1
            c["H3_truth_in_family"] += in_fam
            c["H3_aligner_right"] += al_right
            if is_assigned:
                c["H3_assigned_right" if (fam, r["catalog_copy_idx"]) in T else "H3_assigned_wrong"] += 1
    n1 = c["H1_assigned"]
    res = {"counts": dict(c),
           "H1_precision": {"assigned_one_locus": n1, "wrong": c["H1_wrong"],
                            "verdict": "UNDERPOWERED" if n1 < 30 else ("PASS" if c["H1_wrong"] / n1 <= 0.01 else "FAIL")},
           "H2_aligner_on_same": {"n": n1, "aligner_right": c["H2_aligner_right"]},
           "H3_sample": {"one_locus": c["H3_one_locus"], "o2_assigned_right": c["H3_assigned_right"],
                         "o2_assigned_wrong": c["H3_assigned_wrong"],
                         "o2_recall": round(c["H3_assigned_right"] / c["H3_one_locus"], 4) if c["H3_one_locus"] else None,
                         "aligner_right": c["H3_aligner_right"]},
           "examples": ex}
    json.dump(res, open(f"{out}/haplotype_truth.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "examples"}, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:8])
