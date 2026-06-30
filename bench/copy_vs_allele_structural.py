#!/usr/bin/env python
"""STRUCTURAL copy-vs-allele resolver on the current haploid T2T genome (the O4 close, WITHOUT DNA reads).

True read-depth parCN (Soto/QuicK-mer2) needs gorilla DNA WGS reads, which we do not have (only RNA IsoSeq).
But a COMPLETE, HAPLOID T2T assembly answers copy-vs-allele structurally:
  - a real COPY (paralog) is resolved as a SECOND distinct locus in the assembly;
  - an ALLELE (het) is the OTHER haplotype, ABSENT from the primary — so it has NO second locus.

So per candidate we combine two assembly-structure signals (the second is INDEPENDENT of the reads):
  (1) n_loci  — how many genome loci the candidate's CONSENSUS maps to (already computed, read-derived).
  (2) SEDEF   — does the candidate's locus sit in a segmental-duplication pair (final.bed = genome self-align,
                never sees the reads)? a distinct SEDEF partner = a second copy in the assembly.

Call:
  COPY            n_loci >= 2  OR  a distinct SEDEF segdup partner  → the assembly resolves >= 2 loci.
  ALLELE/HET      single locus, no SEDEF partner, divergence in the het window (2-20%) → one-copy gene, het variant.
  DIVERGENT-NOVEL single locus, no SEDEF, divergence > 20% → too divergent for a het AND no second locus
                  (a paralog the T2T truly lacks, or an artifact) — not callable without DNA / a 2nd haplotype.
  ~REF            divergence < 2% → essentially the reference allele.
A locus that is BOTH multi-copy AND het-range is flagged AMBIGUOUS (e.g. MHC: many paralogs AND high het).

Run: /home/juanfra/miniforge3/bin/python bench/copy_vs_allele_structural.py
"""
import json
import os
from bisect import bisect_right
from collections import defaultdict

R = "/home/juanfra/winloci_scratch/refabsent"
GW = f"{R}/gw_promoted"
SEDEF = next((p for p in ["/mnt/c/Users/jfris/Desktop/final.bed",
                          "/mnt/c/Users/jfris/Desktop/Rustle/final.bed"] if os.path.exists(p)), None)
HET_LO, HET_HI = 2.0, 20.0   # divergence window (%) where a variant could plausibly be a het allele


def load_sedef(path):
    """per-chrom sorted segdup intervals -> (start, end, partner_chrom, partner_start, partner_end). Each
    final.bed line is ONE segdup pair (region A cols 0-2, partner B cols 3-5); index BOTH sides."""
    by = defaultdict(list)
    with open(path) as fh:
        for ln in fh:
            f = ln.split("\t")
            if len(f) < 6:
                continue
            try:
                ca, sa, ea, cb, sb, eb = f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])
            except ValueError:
                continue  # header / malformed line
            by[ca].append((sa, ea, cb, sb, eb))
            by[cb].append((sb, eb, ca, sa, ea))
    for c in by:
        by[c].sort()
    return by


def sedef_partners(by, chrom, start, end):
    """distinct partner loci whose segdup region overlaps [start,end] on `chrom`, at a DIFFERENT position."""
    ivs = by.get(chrom, [])
    if not ivs:
        return []
    starts = [x[0] for x in ivs]
    i = bisect_right(starts, end)
    out = []
    for j in range(i - 1, -1, -1):
        s, e, pc, ps, pe = ivs[j]
        if e < start:
            if start - s > 5_000_000:   # intervals are start-sorted; stop after a long non-overlapping run
                break
            continue
        if s <= end and e >= start:     # overlap
            if not (pc == chrom and abs(ps - start) < 1000):   # partner is a genuinely different locus
                out.append((pc, ps, pe))
    # dedupe partners by (chrom, ~position)
    seen, uniq = set(), []
    for pc, ps, pe in out:
        key = (pc, ps // 10000)
        if key not in seen:
            seen.add(key); uniq.append((pc, ps, pe))
    return uniq


def classify(div, n_loci, n_sedef):
    multicopy = (n_loci is not None and n_loci >= 2) or n_sedef >= 1
    het_range = div is not None and HET_LO <= div <= HET_HI
    if div is not None and div < HET_LO:
        return "~REF", "divergence <2% — essentially the reference allele"
    if multicopy and het_range:
        return "AMBIGUOUS", f"multi-copy (n_loci={n_loci}, sedef_partners={n_sedef}) AND het-range div={div}% — could be a new copy OR an allele of an existing copy (MHC-class)"
    if multicopy:
        return "COPY", f"assembly resolves >=2 loci (n_loci={n_loci}, sedef_partners={n_sedef}) — structurally a copy"
    if het_range:
        return "ALLELE/HET", f"single locus (n_loci={n_loci}, no SEDEF partner) + het-range div={div}% — a het variant of a one-copy gene"
    if div is not None and div > HET_HI:
        return "DIVERGENT-NOVEL", f"single locus, no SEDEF, div={div}% too divergent for a het — a paralog the T2T lacks OR an artifact (needs DNA / 2nd haplotype)"
    return "UNRESOLVED", f"insufficient signal (n_loci={n_loci}, div={div})"


def main():
    assert SEDEF, "final.bed (SEDEF) not found"
    print(f"loading SEDEF segdup map {SEDEF} ...")
    sed = load_sedef(SEDEF)
    print(f"  {sum(len(v) for v in sed.values())} segdup intervals over {len(sed)} contigs\n")

    absent = {a["cid"]: a for a in json.load(open(f"{GW}/gw_reference_absent_copies.json"))}
    disc = {d["cid"]: d for d in json.load(open(f"{GW}/gw_discriminated.json"))}
    rows = []
    for cid in sorted(set(absent) | set(disc)):
        a, d = absent.get(cid, {}), disc.get(cid, {})
        chrom = a.get("chrom") or cid.rsplit("_", 1)[0]
        start = a.get("start", int(cid.rsplit("_", 1)[1]) if "_" in cid else 0)
        end = a.get("end", start + a.get("cons_len", 0))
        div = a.get("divergence", (d.get("div") and d["div"] * 100 if isinstance(d.get("div"), float) and d["div"] < 1 else d.get("div")))
        n_loci = d.get("n_loci")
        partners = sedef_partners(sed, chrom, start, end)
        call, reason = classify(div, n_loci, len(partners))
        rows.append(dict(cid=cid, chrom=chrom, start=start, end=end, divergence=div, n_loci=n_loci,
                         sedef_partners=len(partners), protein=a.get("protein") or d.get("protein"),
                         call=call, reason=reason))

    out = "bench/copy_vs_allele_structural.tsv"
    cols = ["cid", "chrom", "start", "end", "divergence", "n_loci", "sedef_partners", "protein", "call", "reason"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    from collections import Counter
    tally = Counter(r["call"] for r in rows)
    print(f"== STRUCTURAL copy-vs-allele on {len(rows)} reference-absent candidates ==")
    for k, v in tally.most_common():
        print(f"  {k:16} {v:4}  ({100*v/len(rows):.0f}%)")
    n_copy = tally.get("COPY", 0)
    n_allele = tally.get("ALLELE/HET", 0)
    print(f"\n  → only {n_copy} structurally support a SECOND locus (COPY); {n_allele} are single-locus het/allele;")
    print(f"     {tally.get('AMBIGUOUS',0)} MHC-class ambiguous; {tally.get('DIVERGENT-NOVEL',0)} divergent-novel (need DNA/2nd-hap).")
    print(f"  SEDEF (independent of reads) confirms a segdup partner for "
          f"{sum(1 for r in rows if r['sedef_partners']>=1)} candidates.")
    print(f"  wrote {out}")
    # the few COPY calls, named
    print("\n  COPY calls (structurally a second locus):")
    for r in rows:
        if r["call"] == "COPY":
            print(f"    {r['cid']}  div={r['divergence']}%  n_loci={r['n_loci']}  sedef_partners={r['sedef_partners']}  prot={r['protein']}")


if __name__ == "__main__":
    main()
