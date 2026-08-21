#!/usr/bin/env python3
"""Annotate a family catalog with its DUPLICATION-MECHANISM certificate.

⚠⚠ THIS IS A CERTIFICATE, NOT A MEMBERSHIP CONDITION. It changes no partition, adds no edge and
removes none. That is deliberate:

  * as an EDGE rule, SD support would make E_r depend on a THIRD input that is species-specific and
    tool-dependent — the gorilla SD catalog is a multi-day SEDEF run, and a new species needs it
    redone. The definition would stop being computable from (genome, reads);
  * it covers 21.7% of families, so as a membership condition it is inert on ~78%;
  * and no edge rule may be adopted on arm evidence alone (2026-08-20: candidates cost 0-6% on the
    164-pair arms and 3.67-12.80% genome-wide).

As a certificate it is free: it states what independent genomic evidence says about each family, and
leaves the definition alone. Same pattern as `lambda` / `cut_certified`.

MECHANISM CLASSES

  SEGMENTAL         every member pair of the family is spanned by one segmental duplication, and both
                    loci lie ENTIRELY INSIDE its two units. The family and the duplicated segment are
                    the same event.
  SD_OVERLAP        an SD spans the pair but does not contain the loci — partial or boundary support.
  RETRO_CANDIDATE   high exon-sum identity with NO genomic SD support at any floor. A retrocopy
                    duplicates the mRNA, so no genomic segment was duplicated and no SD can exist by
                    construction. ⚠ A CANDIDATE: an old segmental duplication whose signature has
                    decayed below the SD floor looks identical here. Separating them needs the
                    intron/flank divergence, which this does not compute.
  UNRESOLVED        neither — older than the SD floor, or not a real family.

⚠ DIRECTION IS NOT CLAIMED. With one genome every duplication edge is symmetric, so this reports
DNA_DUPLICATION, never DERIVED_FROM: it does not say which copy is ancestral. Nor does it exclude
that a >=90% SD signature reflects later gene conversion rather than the original event —
containment dates the SEGMENT, not the family.
"""
import bisect, collections, csv, sys

SD_FLOOR_CONTAIN = 0.90
SD_FLOOR_ANY     = 0.80
RETRO_IDENT      = 0.95     # exon-sum identity above which a no-SD family is a retrocopy candidate
MIN_UNIT         = 1000


def load_sd(path):
    by = collections.defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 34:
            continue
        try:
            ident = float(f[-1])
        except ValueError:
            continue
        a0, a1, b0, b1 = int(f[1]), int(f[2]), int(f[4]), int(f[5])
        if ident < SD_FLOOR_ANY or (a1 - a0) < MIN_UNIT or (b1 - b0) < MIN_UNIT:
            continue
        by[f[0]].append((a0, a1, f[3], b0, b1, ident))
        by[f[3]].append((b0, b1, f[0], a0, a1, ident))
    for c in by:
        by[c].sort()
    return by


def main(copies_tsv, families_tsv, sd_bed, out_tsv, ava_pafs=()):
    sd = load_sd(sd_bed)
    loc, fam = {}, collections.defaultdict(list)
    for r in csv.DictReader(open(copies_tsv), delimiter="\t"):
        k = f"{r['family_id']}~{r['copy_idx']}"
        loc[k] = (r["chrom"], int(r["start"]), int(r["end"]))
        fam[r["family_id"]].append(k)

    def link(a, b):
        """best (identity, contained) over SD pairs spanning a and b"""
        ca, sa, ea = loc[a]; cb, sb, eb = loc[b]
        best_id, best_contain = None, False
        for u0, u1, pc, p0, p1, ident in sd.get(ca, []):
            if u0 >= ea:
                break
            if u1 <= sa or pc != cb or not (p0 < eb and sb < p1):
                continue
            if best_id is None or ident > best_id:
                best_id = ident
            if sa >= u0 and ea <= u1 and sb >= p0 and eb <= p1 and ident >= SD_FLOOR_CONTAIN:
                best_contain = True
        return best_id, best_contain

    # best within-family exon-sum identity, for the retro call
    ident = collections.defaultdict(float)
    for paf, floor in ava_pafs:
        for line in open(paf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[0] == f[5]:
                continue
            de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
            i = (1.0 - de) if de is not None else int(f[9]) / max(int(f[10]), 1)
            k = (min(f[0], f[5]), max(f[0], f[5]))
            if i > ident[k]:
                ident[k] = i

    cert = {}
    for g, members in fam.items():
        best_id, contained, best_rna = None, False, 0.0
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                a, b = members[i], members[j]
                v, c = link(a, b)
                if v is not None and (best_id is None or v > best_id):
                    best_id = v
                contained |= c
                best_rna = max(best_rna, ident.get((min(a, b), max(a, b)), 0.0))
        if contained:
            mech = "SEGMENTAL"
        elif best_id is not None:
            mech = "SD_OVERLAP"
        elif best_rna >= RETRO_IDENT:
            mech = "RETRO_CANDIDATE"
        else:
            mech = "UNRESOLVED"
        cert[g] = (mech, f"{best_id:.4f}" if best_id is not None else "NA",
                   "true" if contained else "false")

    rows = list(csv.DictReader(open(families_tsv), delimiter="\t"))
    hdr = list(rows[0].keys()) + ["dup_mechanism", "sd_identity", "sd_contained"]
    with open(out_tsv, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr, delimiter="\t")
        w.writeheader()
        for r in rows:
            m, i, c = cert.get(r["family_id"], ("UNRESOLVED", "NA", "false"))
            r.update(dup_mechanism=m, sd_identity=i, sd_contained=c)
            w.writerow(r)
    n = len(rows)
    print(f"wrote {out_tsv}  ({n} families)")
    print("\n=== duplication-mechanism certificate ===")
    for m, k in collections.Counter(v[0] for v in cert.values()).most_common():
        print(f"  {m:<18} {k:>5}/{n}  = {k/n:.4f}")


if __name__ == "__main__":
    G = "/mnt/linuxdisk/home/juanfraitu/o1_gw"
    main(f"{G}/ggo_gw.copies.tsv", f"{G}/ggo_gw.families.tsv",
         "/mnt/c/Users/jfris/Desktop/final.bed", f"{G}/ggo_gw.families.sd.tsv",
         ava_pafs=[(f"{G}/gw_ava_asm20.paf", 0.80), (f"{G}/gw_ava_sens.paf", 0.60)])
