#!/usr/bin/env python3
"""Held-out PSV concordance: does real RBMY read data carry a CONSISTENT per-copy
fingerprint, and for which copies? (Real-data validation of copy attribution.)

The problem: real RBMY reads have NO per-read source-copy label, and there is no
external per-copy annotation — so we cannot measure attribution accuracy directly.
The leakage-free substitute: a read that genuinely came from copy X carries copy X's
paralogous sequence variants (PSVs) CONSISTENTLY along its whole length. So:

  1. Align each read to EACH of the 6 RBMY1 reference copies separately (minimap2
     spliced), and in the READ's own coordinate frame record match/mismatch vs every
     copy at each read position (pysam aligned_pairs with_seq — no MSA, no genome
     coordinate joins).
  2. DISCRIMINATING read-positions = positions where the copies disagree (the read
     matches some copies, mismatches others) — i.e. PSV columns as seen by this read.
  3. best-copy(full) = the copy the read matches at the most discriminating positions;
     decisiveness gap = (#matches best) − (#matches 2nd-best).
  4. HELD-OUT test: split a read's discriminating positions into ODD and EVEN halves
     (disjoint, no leakage). best-copy(odd) vs best-copy(even). CONCORDANT iff equal.
     High concordance ⇒ the read carries a real, consistent copy signal. Chance-level
     ⇒ the "attribution" is noise / the copies are non-identifiable for that read —
     exactly where the EM SHOULD abstain.

Honest outputs: overall concordance vs the chance baseline; concordance + decisiveness
stratified by best-copy and by #discriminating-positions (the identifiability spectrum
on REAL data). A null result (concordance ≈ chance) is a valid, informative finding.
"""
import sys, os, subprocess, json, random
from collections import defaultdict
import pysam

COPIES_FA = "/tmp/rbmy_psv/copies.fa"
BAM = "/tmp/tspy.bam"
WORK = "/tmp/rbmy_psv"
# RBMY1 tandem array, ordered by genomic start. c5 = LOC129530242 (divergent, shortest).
COPY_ORDER = ["19602754-19616644", "19625715-19639621", "19648670-19662577",
              "19671638-19685525", "19694606-19708531", "19717578-19730926"]
COPY_LABEL = {c: f"c{i}" for i, c in enumerate(COPY_ORDER)}  # c0..c5 by position

def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)

def main():
    os.makedirs(WORK, exist_ok=True)
    # 1. distinct read sequences (primary, mapped, reasonable length) from the RBMY bam
    reads = {}
    with pysam.AlignmentFile(BAM, "rb") as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if r.query_sequence and len(r.query_sequence) >= 300 and r.query_name not in reads:
                reads[r.query_name] = r.query_sequence
    with open(f"{WORK}/reads.fa", "w") as fh:
        for n, s in reads.items():
            fh.write(f">{n}\n{s}\n")
    print(f"[rbmy] {len(reads)} distinct primary reads (>=300bp)")

    # split copies.fa into one file per copy, labelled c0..c5 by genomic order
    seqs, name, buf = {}, None, []
    for line in open(COPIES_FA):
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name = line[1:].strip().split(":")[1] if ":" in line else line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    # 2. align reads to EACH copy separately; record per-read-position match vs that copy
    #    match_vec[read][copy_label] = dict(read_pos -> 1 match / 0 mismatch)
    match_vec = defaultdict(dict)
    for span in COPY_ORDER:
        lab = COPY_LABEL[span]
        cfa = f"{WORK}/{lab}.fa"
        open(cfa, "w").write(f">{lab}\n{seqs[span]}\n")
        sh(f"minimap2 -ax splice:hq -uf --MD {cfa} {WORK}/reads.fa 2>/dev/null "
           f"| samtools sort -o {WORK}/{lab}.bam - 2>/dev/null && samtools index {WORK}/{lab}.bam")
        with pysam.AlignmentFile(f"{WORK}/{lab}.bam", "rb") as bam:
            for r in bam:
                if r.is_secondary or r.is_supplementary or r.is_unmapped:
                    continue
                d = {}
                for qpos, rpos, rbase in r.get_aligned_pairs(with_seq=True):
                    if qpos is None or rpos is None or rbase is None:
                        continue  # indel / softclip
                    qbase = r.query_sequence[qpos]
                    # pysam with_seq: ref base lowercase => mismatch, uppercase => match
                    d[qpos] = 1 if rbase.upper() == qbase.upper() and rbase.isupper() else 0
                if d:
                    # keep the read's BEST alignment to this copy (most matches) if multiple
                    prev = match_vec[r.query_name].get(lab)
                    if prev is None or sum(d.values()) > sum(prev.values()):
                        match_vec[r.query_name][lab] = d

    # 3-4. per read: discriminating positions, best copy, decisiveness gap, held-out concordance
    labels = [COPY_LABEL[c] for c in COPY_ORDER]
    def best_at(read_d, positions):
        """argmax copy by #matches over `positions`; returns (best_label, gap, scores)."""
        scores = {}
        for lab in labels:
            mv = read_d.get(lab)
            if mv is None:
                scores[lab] = None; continue
            scores[lab] = sum(mv.get(p, 0) for p in positions)
        present = {l: s for l, s in scores.items() if s is not None}
        if len(present) < 2:
            return None, 0, scores
        ordered = sorted(present.items(), key=lambda kv: -kv[1])
        gap = ordered[0][1] - ordered[1][1]
        return ordered[0][0], gap, scores

    # Held-out via AVERAGED RANDOM-HALF splits (not odd/even — a deterministic odd/even
    # split manufactures a parity correlation: it spuriously inflates one copy of a near-tie
    # to 1.0 and pushes its sister BELOW chance. Audited 2026-06-04). Per read: K random
    # disjoint halves, concordance = fraction where the two halves pick the same best copy.
    # Identifiability is gated on the NORMALIZED MARGIN (gap / n_disc) — net discriminating
    # matches per position — which cleanly separates a real signal (c0 ~0.44) from a coin-flip
    # near-tie (~0.03-0.04, ~1.5 net matches over ~48 positions).
    rng = random.Random(0)
    K = 50
    rows = []
    for rn, read_d in match_vec.items():
        covered = [l for l in labels if read_d.get(l)]
        if len(covered) < 2:
            continue
        common = set.intersection(*[set(read_d[l].keys()) for l in covered])
        disc = [p for p in common if len({read_d[l][p] for l in covered}) > 1]
        if len(disc) < 2:
            rows.append({"read": rn, "n_disc": len(disc), "best": None, "gap": 0,
                         "norm_margin": 0.0, "concordance": None, "n_copies": len(covered)})
            continue
        best, gap, _ = best_at(read_d, disc)
        nmargin = gap / len(disc)
        agree = 0
        for _ in range(K):
            shuf = disc[:]; rng.shuffle(shuf); h = len(shuf) // 2
            ba, _, _ = best_at(read_d, shuf[:h])
            bb, _, _ = best_at(read_d, shuf[h:])
            if ba is not None and ba == bb:
                agree += 1
        rows.append({"read": rn, "n_disc": len(disc), "best": best, "gap": gap,
                     "norm_margin": round(nmargin, 3), "concordance": agree / K,
                     "n_copies": len(covered)})

    scored = [r for r in rows if r["concordance"] is not None]
    n = len(scored)
    overall = sum(r["concordance"] for r in scored) / n if n else 0.0
    chance = sum(1.0 / r["n_copies"] for r in scored) / n if n else 0.0
    NM_GATE = 0.10  # normalized-margin identifiability threshold

    per_copy = defaultdict(lambda: {"n": 0, "conc": [], "nm": []})
    for r in scored:
        pc = per_copy[r["best"]]
        pc["n"] += 1; pc["conc"].append(r["concordance"]); pc["nm"].append(r["norm_margin"])

    out = {
        "n_reads_scored": n,
        "heldout_concordance_randomhalf": round(overall, 4),
        "chance_baseline": round(chance, 4),
        "norm_margin_gate": NM_GATE,
        "per_copy": {k: {"n": v["n"],
                          "concordance": round(sum(v["conc"]) / v["n"], 4) if v["n"] else None,
                          "mean_norm_margin": round(sum(v["nm"]) / v["n"], 3) if v["n"] else 0,
                          "identifiable": (sum(v["nm"]) / v["n"] > NM_GATE) if v["n"] else False}
                      for k, v in sorted(per_copy.items())},
        "note": "c0 = LOC129530243, c5 = LOC129530242 (divergent/starved). IDENTIFIABLE = mean "
                "normalized margin > gate: only a copy with a real distributed per-read signal "
                "clears it. Near-ties (~chance concordance, tiny margin) are NON-identifiable — "
                "the EM must abstain there. Random-half avoids the odd/even parity artifact.",
    }
    print(json.dumps(out, indent=2))
    json.dump(out, open(f"{WORK}/heldout_results.json", "w"), indent=2)

if __name__ == "__main__":
    main()
