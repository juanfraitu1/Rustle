#!/usr/bin/env python3
"""Amendment 44d: the frozen read-share rule (tau = 0.02) scored once on the held-out set."""
import collections
import json
import sys

sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/rustle_m2_soto/bench/unmapped_rescue")
import locus_support as L  # noqa: E402
import spec_heldout as SH  # noqa: E402

O = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"
TAU = 0.02


def main():
    lab = json.load(open(f"{O}/dna_labels.json"))
    ls = {json.loads(l)["k"]: json.loads(l) for l in open(f"{O}/locus_support.jsonl")}
    ev = json.load(open(f"{O}/eval_binned.json"))
    S = {"own": SH.load(f"{O}/spec_flags.cs.paf"), "human": SH.load(f"{O}/spec_held.hsa.paf"), "chimp": SH.load(f"{O}/spec_held.ptr.paf"),
         "orangutan": SH.load(f"{O}/spec_held.ppy.paf"), "siamang": SH.load(f"{O}/spec_held.ssy.paf")}
    held = {k: x for k, x in lab.items() if x["side"] == "held"}

    def keep(k):
        r = ls[k]
        if r["n_cons"] is None:
            return True
        s = L.share(r["n_cons"], r["n_ref"])
        return s is None or s >= TAU

    labels = {"44c (refined)": lambda x: "real" if x["verdict"] == "TRUE" or (x["lab_all"] == "DNA-SUPPORTED" and (x["med_all"] or 0) <= 25) else "artifact" if x["lab_all"] == "RNA-ONLY" else "undecided",
              "43b": lambda x: "real" if x["verdict"] == "TRUE" or x["lab_all"] == "DNA-SUPPORTED" else "artifact" if x["lab_all"] == "RNA-ONLY" else "undecided"}
    for rname, fn in (("read share >= 0.02 (frozen)", keep), ("read share AND species screen", lambda k: keep(k) and SH.keep(k, S))):
        for lname, cls in labels.items():
            c = collections.Counter((cls(x), fn(k)) for k, x in held.items())
            art = c[("artifact", True)] + c[("artifact", False)]
            real = c[("real", True)] + c[("real", False)]
            print(f"HELD-OUT [{rname}] labels {lname}: S1 artifacts removed {c[('artifact', False)]}/{art} = {c[('artifact', False)] / art:.1%}; "
                  f"S2 real kept {c[('real', True)]}/{real} = {c[('real', True)] / real:.1%}")
        tr = [k for k, x in held.items() if x["verdict"] == "TRUE"]
        print(f"   assembly-TRUE kept {sum(map(fn, tr))}/{len(tr)} = {sum(map(fn, tr)) / len(tr):.1%}")
        chrom = {ln.rstrip("\n").split("\t")[6]: ln.split("\t")[1] for ln in list(open("/mnt/linuxdisk/tmp/rna_allele/chrmap.tsv"))[1:] if ln.count("\t") >= 6 and ln.rstrip("\n").split("\t")[6]}
        loci = json.load(open(f"{O}/loci.json"))
        is_held = lambda l: chrom.get(l["contig"]) is not None and (not chrom[l["contig"]].isdigit() or int(chrom[l["contig"]]) % 2 == 0)
        l3 = [i for i, l in enumerate(loci) if len(l["names"]) >= 3 and is_held(l)]
        rec0 = {v["locus"] for v in ev["verdicts"] if v["verdict"] == "TRUE" and v["locus"] is not None}
        rec1 = {v["locus"] for v in ev["verdicts"] if v["verdict"] == "TRUE" and v["locus"] is not None and fn(v["k"])}
        r0, r1 = sum(i in rec0 for i in l3) / len(l3), sum(i in rec1 for i in l3) / len(l3)
        kept = [k for k in held if fn(k)]
        w = sum(held[k]["verdict"] == "WRONG" for k in kept)
        print(f"   S3 recall {r0:.3f} -> {r1:.3f} (relative fall {(r0 - r1) / r0:.1%}); kept flags {len(kept)}, WRONG among them {w} = {w / len(kept):.1%} (was {sum(x['verdict'] == 'WRONG' for x in held.values()) / len(held):.1%})")
        ston = [i for i, l in enumerate(loci) if l["contig"] == "CM054594.2" and l["start"] < 96_100_000 and l["end"] > 95_900_000]
        print("   STON1-GTF2A1L recovered:", [(i, i in rec1) for i in ston])


if __name__ == "__main__":
    main()
