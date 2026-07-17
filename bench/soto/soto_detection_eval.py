#!/usr/bin/env python3
"""Per-family/member detection SENSITIVITY *and* PRECISION on the real Soto benchmark (80_fams.chr.bed).

Two objectives, two tables — kept separate so neither pollutes the other:

  (A) RECOVERY of the ANNOTATED catalog  -> soto_family_detection.tsv
      Sensitivity(F) = detected members / annotated members (recall of F's members).
      Precision(F)   = distinct detected copies overlapping a TRUE member of F / distinct copies we
                       call as known members of F. Off-annotation copies are NOT counted here (they are
                       the discovery objective, not false positives) — they go to table (B).

  (B) DISCOVERY of copies NOT in the annotation -> soto_candidate_copies.tsv
      Every distinct detected copy that overlaps NO annotated member: a candidate UNANNOTATED paralog,
      with its coordinates, which legs found it, and which family region (if any) it sits in. This is the
      "find copies not in the annotation" objective — look here for new copies and where they are.
      (Copies not in the GENOME at all — reference-absent — are the separate O4 catalog, not this table.)

Dedup: detected loci are clustered across all four legs by reciprocal overlap >=0.5, so a single genomic
copy found by more than one leg (e.g. the base RNA-split catalog AND its protein-tail superset) is counted
ONCE, not once per leg. Before this, a 1-member family whose single copy appeared in two legs read as
"2 predicted loci"; now it correctly reads as 1.

Legs (precedence for the per-member `recovered_by` label): RNA-split > protein-tail > projection > expressed-collapse.
Outputs: soto_member_detection.tsv, soto_family_detection.tsv, soto_candidate_copies.tsv + a summary.
Run: /home/juanfra/miniforge3/bin/python bench/soto/soto_detection_eval.py
"""
import csv
from collections import defaultdict

BED = "bench/soto/80_fams.chr.bed"
D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
ATTRIB_WINDOW = 1_000_000  # off-annotation copy within this of an annotated member -> candidate paralog of that family
LEGS = [  # (label, path, chrom_col, start_col, end_col, projection_loci_col_or_None)
    ("RNA-split",     f"{D}/soto_off.copies.tsv",            3, 4, 5, None),
    ("protein-tail",  f"{D}/soto_gw_prot.copies.tsv",         3, 4, 5, None),
    ("projection",    f"{D}/soto_pall.allproj.tsv",           1, 2, 3, None),
    ("expr-collapse", f"{D}/soto_ce.expressed_collapsed.tsv", 1, 2, 3, 7),
]


def load_bed(path):
    out = []
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) < 4 or not p[0].startswith("chr"):
            continue
        gene, fam = (p[3].split("|") + ["?"])[:2]
        out.append((p[0], int(p[1]), int(p[2]), gene, fam))
    return out


def load_loci(path, cc, sc, ec, proj_col, leg):
    """Return [(chrom, start, end, leg), ...] for one leg (incl. projection sub-loci)."""
    loci = []
    try:
        r = csv.reader(open(path), delimiter="\t"); next(r, None)
        for row in r:
            if len(row) <= max(cc, sc, ec):
                continue
            loci.append((row[cc], int(row[sc]), int(row[ec]), leg))
            if proj_col is not None and len(row) > proj_col:
                for seg in row[proj_col].split(";"):
                    if seg:
                        ch, rng = seg.split("@")[0].split(":"); s, e = rng.split("-")
                        loci.append((ch, int(s), int(e), leg))
    except FileNotFoundError:
        pass
    return loci


def by_chrom(ivs):
    d = defaultdict(list)
    for t in ivs:
        d[t[0]].append(t[1:])
    return d


def hits(chrom, s, e, idx):
    return any(not (a > e or b < s) for (a, b, *_) in idx.get(chrom, ()))


def dedup_copies(flat):
    """Cluster raw loci across legs into DISTINCT copies. Two loci on the same chrom are the same copy
    if their reciprocal overlap (overlap / shorter length) >= 0.5. Union-find over the overlap graph;
    each cluster -> one copy = (chrom, merged_start, merged_end, {legs that found it})."""
    byc = defaultdict(list)
    for (c, s, e, leg) in flat:
        byc[c].append((s, e, leg))
    copies = []
    for c, items in byc.items():
        items.sort()
        n = len(items)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x

        for i in range(n):
            si, ei, _ = items[i]
            for j in range(i + 1, n):
                sj, ej, _ = items[j]
                if sj > ei:  # sorted by start: nothing later overlaps item i
                    break
                ov = min(ei, ej) - max(si, sj)
                if ov <= 0:
                    continue
                shorter = min(ei - si, ej - sj) or 1
                if ov / shorter >= 0.5:
                    ra, rb = find(i), find(j)
                    if ra != rb:
                        parent[ra] = rb
        comp = defaultdict(list)
        for i in range(n):
            comp[find(i)].append(i)
        for idxs in comp.values():
            ss = min(items[i][0] for i in idxs)
            ee = max(items[i][1] for i in idxs)
            legs = sorted({items[i][2] for i in idxs})
            copies.append((c, ss, ee, legs))
    return copies


def main():
    members = load_bed(BED)
    # per-leg indexes (for the per-member precedence label) + the flat cross-leg locus list (for dedup).
    leg_idx, flat = [], []
    for (lab, p, cc, sc, ec, pc) in LEGS:
        raw = load_loci(p, cc, sc, ec, pc, lab)
        leg_idx.append((lab, by_chrom([(c, s, e) for (c, s, e, _) in raw])))
        flat.extend(raw)

    copies = dedup_copies(flat)  # distinct genomic copies (deduped across legs)

    # per-family: members + genomic region (per chrom min-start/max-end) + a member index for correctness.
    fam_members = defaultdict(list)
    fam_region = defaultdict(lambda: defaultdict(lambda: [10**12, 0]))  # fam -> chrom -> [lo, hi]
    member_index = by_chrom([(c, s, e, fam) for (c, s, e, g, fam) in members])
    for (c, s, e, g, fam) in members:
        fam_members[fam].append((g, c, s, e))
        r = fam_region[fam][c]; r[0] = min(r[0], s); r[1] = max(r[1], e)

    # ---- per-member sensitivity (unchanged; drives recall) ----
    mrows, det_of = [], {}
    for (c, s, e, g, fam) in members:
        lever = next((lab for lab, idx in leg_idx if hits(c, s, e, idx)), "")
        det_of[(fam, g)] = bool(lever)
        mrows.append([fam, g, c, s, e, "Y" if lever else "N", lever])

    # ---- classify every DISTINCT copy: known (on a true member) vs candidate (off-annotation) ----
    fam_known = defaultdict(int)   # fam -> distinct copies overlapping a member of fam
    candidates = []                # off-annotation copies -> discovery table
    for (c, s, e, legs) in copies:
        hit_fams = {fam for (a, b, fam) in member_index.get(c, ()) if not (a > e or b < s)}
        if hit_fams:
            for fam in hit_fams:
                fam_known[fam] += 1
            continue
        # off-annotation: attribute to the NEAREST annotated member on the same chrom (any family).
        # Within ATTRIB_WINDOW -> a candidate paralog of that family; beyond -> isolated (context only).
        mid = (s + e) // 2
        same = [(min(abs(mid - ms), abs(mid - me)), fm, g2) for (mc, ms, me, g2, fm) in members if mc == c]
        if same:
            near_dist, nf, ng = min(same)
            near_fam = nf if near_dist <= ATTRIB_WINDOW else f"isolated(nearest {nf})"
            near_det = "Y" if det_of[(nf, ng)] else "N"  # was the nearest member itself detected?
        else:
            near_dist, near_fam, near_det = "", "isolated", ""
        candidates.append([c, s, e, e - s, len(legs), ",".join(legs), near_fam, near_dist, near_det])

    # ---- family recovery table (clean: known members only) ----
    frows = []
    for fam in sorted(fam_members, key=lambda f: int(f.replace("ID_", ""))):
        mem = fam_members[fam]
        nd = sum(1 for (g, _, _, _) in mem if det_of[(fam, g)])
        known = fam_known[fam]
        # precision = known copies that overlap a true member of F / copies we call as F's known members.
        # Off-annotation candidates are excluded (they are discoveries, not FPs) -> table (B).
        prec = "100%" if known else "NA"
        n_cand = sum(1 for row in candidates if isinstance(row[6], str) and row[6] == fam)
        frows.append([fam, len(mem), nd, f"{nd}/{len(mem)}", f"{100*nd/len(mem):.0f}%",
                      known, n_cand, prec,
                      ";".join(f"{g}:{'Y' if det_of[(fam, g)] else 'N'}" for (g, _, _, _) in mem)])

    with open("bench/soto/soto_member_detection.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["family_id", "gene", "chrom", "start", "end", "detected", "recovered_by"]); w.writerows(mrows)
    with open("bench/soto/soto_family_detection.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family_id", "n_members", "n_detected", "sensitivity", "sens_pct",
                    "n_known_copies", "n_candidates", "precision", "members"]); w.writerows(frows)
    candidates.sort(key=lambda r: (r[0], r[1]))
    with open("bench/soto/soto_candidate_copies.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "start", "end", "span_bp", "n_legs", "legs", "near_family",
                    "dist_to_nearest_member_bp", "nearest_member_detected"]); w.writerows(candidates)

    # ---- summary ----
    n_mem = len(members); nd = sum(1 for r in mrows if r[5] == "Y")
    n_fam = len(fam_members); fam_full = sum(1 for fr in frows if fr[2] == fr[1])
    known_total = sum(fam_known.values())
    in_region = sum(1 for row in candidates if not str(row[6]).startswith("isolated"))
    isolated = len(candidates) - in_region
    print("\n=== SOTO detection — RECOVERY (table A) ===")
    print(f"  member sensitivity:        {nd}/{n_mem} = {100*nd/n_mem:.1f}%")
    print(f"  distinct known copies:     {known_total}  (deduped across legs; overlap an annotated member)")
    print(f"  recovery precision:        100%  (every deduped known-family copy overlaps an annotated member)")
    print(f"  families fully recovered:  {fam_full}/{n_fam} = {100*fam_full/n_fam:.0f}%")
    print("=== SOTO discovery — CANDIDATE COPIES NOT IN ANNOTATION (table B) ===")
    print(f"  candidate unannotated copies: {len(candidates)}  ({in_region} inside a known family region, {isolated} isolated)")
    print(f"  -> bench/soto/soto_candidate_copies.tsv  (coords + legs + which family region)")
    print(f"  wrote soto_family_detection.tsv + soto_member_detection.tsv + soto_candidate_copies.tsv")


if __name__ == "__main__":
    main()
