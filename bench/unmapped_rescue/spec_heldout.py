#!/usr/bin/env python3
"""Amendment 43c: the frozen species rule scored once on the held-out set (even chromosomes, X, Y, chrM, unmapped bin). Miniforge python."""
import collections
import json

O = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"
D = 0.00958


def load(p):
    best = {}
    for ln in open(p):
        x = ln.rstrip("\n").split("\t")
        if "tp:A:P" in x[12:] and (x[0] not in best or int(x[9]) > int(best[x[0]][9])):
            best[x[0]] = x
    return {k: int(x[9]) / int(x[1]) for k, x in best.items()}


def keep(k, S):
    """the frozen rule: own primary within delta of the best panel score"""
    return S["own"].get(k, 0) >= max(s.get(k, 0) for s in S.values()) - D


def main():
    lab = json.load(open(f"{O}/dna_labels.json"))
    ev = json.load(open(f"{O}/eval_binned.json"))
    S = {"own": load(f"{O}/spec_flags.cs.paf"), "human": load(f"{O}/spec_held.hsa.paf"), "chimp": load(f"{O}/spec_held.ptr.paf"),
         "orangutan": load(f"{O}/spec_held.ppy.paf"), "siamang": load(f"{O}/spec_held.ssy.paf")}
    held = {k: x for k, x in lab.items() if x["side"] == "held"}
    cls = lambda x: "real" if x["verdict"] == "TRUE" or x["lab_all"] == "DNA-SUPPORTED" else "artifact" if x["lab_all"] == "RNA-ONLY" else "undecided"
    c = collections.Counter((cls(x), keep(k, S)) for k, x in held.items())
    art = c[("artifact", True)] + c[("artifact", False)]
    real = c[("real", True)] + c[("real", False)]
    s1, s2 = c[("artifact", False)] / art, c[("real", True)] / real
    print(f"HELD-OUT flags {len(held)}: artifacts {art}, real {real}, undecided {c[('undecided', True)] + c[('undecided', False)]}")
    print(f"S1 artifacts removed {c[('artifact', False)]}/{art} = {s1:.1%} (bar >= 80%): {'MET' if s1 >= 0.8 else 'NOT met'}")
    print(f"S2 real kept {c[('real', True)]}/{real} = {s2:.1%} (bar >= 90%): {'MET' if s2 >= 0.9 else 'NOT met'}")
    for g, sel in (("assembly-TRUE", lambda x: x["verdict"] == "TRUE"), ("WRONG", lambda x: x["verdict"] == "WRONG")):
        ks = [k for k, x in held.items() if sel(x)]
        print(f"   {g}: {len(ks)} flags, kept {sum(keep(k, S) for k in ks)}")
    # recall over held-out truth loci (>= 3 mat-specific reads), chromosome from the maternal accession
    chrom = {ln.rstrip("\n").split("\t")[6]: ln.split("\t")[1] for ln in list(open("/mnt/linuxdisk/tmp/rna_allele/chrmap.tsv"))[1:] if ln.count("\t") >= 6 and ln.rstrip("\n").split("\t")[6]}
    loci = json.load(open(f"{O}/loci.json"))
    is_held = lambda l: chrom.get(l["contig"]) is not None and (not chrom[l["contig"]].isdigit() or int(chrom[l["contig"]]) % 2 == 0)
    l3 = [i for i, l in enumerate(loci) if len(l["names"]) >= 3 and is_held(l)]
    rec0 = {v["locus"] for v in ev["verdicts"] if v["verdict"] == "TRUE" and v["locus"] is not None}
    rec1 = {v["locus"] for v in ev["verdicts"] if v["verdict"] == "TRUE" and v["locus"] is not None and keep(v["k"], S)}
    r0, r1 = sum(i in rec0 for i in l3) / len(l3), sum(i in rec1 for i in l3) / len(l3)
    print(f"S3 recall over held-out truth loci with >= 3 reads ({len(l3)}): without the rule {r0:.3f}, with it {r1:.3f}; relative fall {(r0 - r1) / r0:.1%} (bar <= 10%): {'MET' if (r0 - r1) <= 0.1 * r0 else 'NOT met'}")
    wrong = [k for k, x in held.items() if x["verdict"] == "WRONG"]
    kept = [k for k in held if keep(k, S)]
    print(f"O3 convention on held-out with the rule: recall {r1:.3f} (bar 0.5); WRONG share of kept flags {sum(held[k]['verdict'] == 'WRONG' for k in kept)}/{len(kept)} = "
          f"{sum(held[k]['verdict'] == 'WRONG' for k in kept) / len(kept):.1%} (bar <= 10%); without the rule {len(wrong)}/{len(held)} = {len(wrong) / len(held):.1%}")
    ston = [i for i, l in enumerate(loci) if l["contig"] == "CM054594.2" and l["start"] < 96_100_000 and l["end"] > 95_900_000]
    print("STON1-GTF2A1L recovered with the rule:", [(i, i in rec1) for i in ston])


if __name__ == "__main__":
    main()
