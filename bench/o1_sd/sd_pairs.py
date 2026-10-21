#!/usr/bin/env python3
"""Are O1's families real duplications? (docs/PREREG_o1_families_sd_pairs_2026-10-21.md)

A copy is SD-PAIRED in its family when one row of an independent segmental-duplication annotation has one side overlapping one
of the copy's exons and the other side overlapping an exon of another copy of the same family (>= 1 bp each).

    sd_pairs.py run OUT.json      the prereg's three catalogs against their SD tracks (paths below), the shuffled-label null,
                                  identity bins, mainstays and failure classes
"""
import bisect
import collections
import csv
import json
import random
import re
import sys

RUNS = "/mnt/linuxdisk/tmp/rustle_figures/runs"
CHM13_SD = "/mnt/c/Users/jfris/Desktop/CHM13_SDs.SD90.bed"
CHM13_SD_NOSAT = "/mnt/c/Users/jfris/Desktop/CHM13_SDs.SD90.nosat.bed"
GGO_SD = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed"
GGO_GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff"
CHM13_GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GCF_009914755.1_T2T-CHM13v2.0_genomic.sorted.gff"
CATALOGS = [
    ("human_A119b", f"{RUNS}/human_A119b/human_A119b.fam.copies.tsv", "human"),
    ("gorilla_KB3781", f"{RUNS}/gorilla_KB3781/gorilla_KB3781.fam.copies.tsv", "gorilla"),
    ("gorilla_fibro_cat", "/mnt/linuxdisk/tmp/rna_allele/fibro_cat.copies.tsv", "gorilla"),
]
YX = {"human": ("chrY", "chrX"), "gorilla": ("NC_073248.2", "NC_073247.2")}
MAIN = re.compile(r"NPIP|TBC1D3|GOLGA6|GOLGA8|SMN[12]?\b|LRPAP1|TSPY|RBMY|DAZ|HSFY|BPY2|CDY", re.I)
MAIN_DESC = re.compile(r"nuclear pore complex|TBC1 domain family member 3|golgin A[68]|survival of motor neuron|LDL receptor related "
                       r"protein associated protein 1|testis.specific (Y.encoded )?protein|RNA.binding motif protein,? Y|deleted in "
                       r"azoospermia|heat shock transcription factor,? Y|basic charge,? Y|chromodomain protein,? Y", re.I)


def _index(copies):
    """{chrom: (sorted exon starts, [(start, end, copy index)], longest exon)}"""
    by = collections.defaultdict(list)
    for i, c in enumerate(copies):
        for s, e in c["exons"]:
            by[c["chrom"]].append((s, e, i))
    out = {}
    for ch, v in by.items():
        v.sort()
        out[ch] = ([x[0] for x in v], v, max(e - s for s, e, _ in v))
    return out


def _hits(idx, ch, s, e):
    if ch not in idx:
        return set()
    st, v, longest = idx[ch]
    lo = bisect.bisect_left(st, s - longest)
    hi = bisect.bisect_left(st, e)
    return {i for a, b, i in v[lo:hi] if a < e and s < b}


def sd_paired(copies, sd_rows):
    """Keys of the copies paired by an SD row with another copy of the same family."""
    idx = _index(copies)
    out = set()
    for c1, s1, e1, c2, s2, e2 in sd_rows:
        x, y = _hits(idx, c1, s1, e1), _hits(idx, c2, s2, e2)
        if not x or not y:
            continue
        for i in x:
            for j in y:
                if i != j and copies[i]["family"] == copies[j]["family"]:
                    out.add(copies[i]["key"])
                    out.add(copies[j]["key"])
    return out


def shuffle_families(copies, seed=1):
    labels = [c["family"] for c in copies]
    random.Random(seed).shuffle(labels)
    return [dict(c, family=f) for c, f in zip(copies, labels)]


def identity_bin(x):
    if x is None:
        return "n/a"
    return ">= 0.98" if x >= 0.98 else "0.95-0.98" if x >= 0.95 else "0.90-0.95" if x >= 0.90 else "< 0.90"


def load_copies(p):
    out = []
    for r in csv.DictReader(open(p), delimiter="\t"):
        ident = r.get("max_family_identity")
        out.append(dict(key=(r["family_id"], r["copy_idx"]), family=r["family_id"], chrom=r["chrom"],
                        exons=[tuple(int(x) for x in e.split("-")) for e in r["exons"].split(",")],
                        start=int(r["start"]), end=int(r["end"]), gene=r.get("gene_id") or r.get("tid", ""),
                        identity=float(ident) if ident not in (None, "", "NA") else None))
    return out


def load_sd(p, min_id=None):
    rows = []
    for ln in open(p):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if min_id is not None and float(f[20]) < min_id:
            continue
        rows.append((f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])))
    return rows


def gene_descriptions(gff_p):
    """[(chrom, start0, end, symbol, description)] of the annotation's genes and pseudogenes."""
    out = []
    for ln in open(gff_p):
        if ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
            continue
        name = re.search(r"Name=([^;]+)", f[8])
        desc = re.search(r"description=([^;]+)", f[8])
        out.append((f[0], int(f[3]) - 1, int(f[4]), name.group(1) if name else "", (desc.group(1) if desc else "").replace("%2C", ",")))
    return out


def mainstay_of(c, genes_idx):
    """The mainstay family name a copy belongs to (by its gene id, or by the symbol / description of an annotated gene its span
    overlaps), else ''."""
    if MAIN.search(c["gene"]):
        return MAIN.search(c["gene"]).group(0).upper()
    for (ch, s, e, sym, desc) in genes_idx.get(c["chrom"], []):
        if s < c["end"] and c["start"] < e:
            m = MAIN.search(sym) or MAIN_DESC.search(desc)
            if m:
                return (MAIN.search(sym).group(0).upper() if MAIN.search(sym) else MAIN_DESC.search(desc).group(0))
    return ""


def summarize(copies, paired, null_paired, yx, genes):
    n = collections.Counter()
    by_bin = collections.defaultdict(lambda: [0, 0, 0])            # bin -> [copies, paired, null paired]
    for c in copies:
        b = identity_bin(c["identity"])
        by_bin[b][0] += 1
        by_bin[b][1] += c["key"] in paired
        by_bin[b][2] += c["key"] in null_paired
    hi = [c for c in copies if c["identity"] is not None and c["identity"] >= 0.90]
    fam_hi = collections.defaultdict(list)
    for c in hi:
        fam_hi[c["family"]].append(c["key"] in paired)
    fails = collections.Counter()
    for c in (hi if hi else copies):
        if c["key"] in paired:
            continue
        fails["chrY" if c["chrom"] == yx[0] else "chrX" if c["chrom"] == yx[1] else
              "de novo copy" if c["gene"].startswith("DN_") else "other"] += 1
    gidx = collections.defaultdict(list)
    for g in genes:
        gidx[g[0]].append(g)
    mains = collections.defaultdict(lambda: [0, 0])
    for c in copies:
        m = mainstay_of(c, gidx)
        if m:
            mains[m.upper()][0] += 1
            mains[m.upper()][1] += c["key"] in paired
    out = dict(copies=len(copies), paired=len(paired & {c["key"] for c in copies}), families=len({c["family"] for c in copies}),
               by_identity={b: dict(copies=v[0], paired=v[1], null_paired=v[2]) for b, v in sorted(by_bin.items())},
               hi_copies=len(hi), hi_paired=sum(c["key"] in paired for c in hi), hi_null_paired=sum(c["key"] in null_paired for c in hi),
               families_all_hi_paired=sum(1 for v in fam_hi.values() if len(v) >= 2 and all(v)),
               families_with_2_hi=sum(1 for v in fam_hi.values() if len(v) >= 2),
               failures=dict(fails), mainstays={k: dict(copies=v[0], paired=v[1]) for k, v in sorted(mains.items())})
    if hi:
        out["bar_holds"] = out["hi_paired"] >= 0.9 * out["hi_copies"]
    return out


def cmd_run(out_p):
    sds = {"human": load_sd(CHM13_SD), "gorilla": load_sd(GGO_SD, min_id=0.90)}
    nosat = load_sd(CHM13_SD_NOSAT)
    genes = {"human": gene_descriptions(CHM13_GFF), "gorilla": gene_descriptions(GGO_GFF)}
    res = {}
    for name, path, sp in CATALOGS:
        copies = load_copies(path)
        paired = sd_paired(copies, sds[sp])
        null = sd_paired(shuffle_families(copies, seed=1), sds[sp])
        res[name] = summarize(copies, paired, null, YX[sp], genes[sp])
        if sp == "human":
            res[name]["nosat"] = summarize(copies, sd_paired(copies, nosat), null, YX[sp], genes[sp])
        r = res[name]
        print(f"[sd_pairs] {name}: copies {r['copies']}, SD-paired {r['paired']} ({100 * r['paired'] / r['copies']:.1f}%); "
              f">= 0.90 copies {r['hi_copies']}, paired {r['hi_paired']}"
              + (f" ({100 * r['hi_paired'] / r['hi_copies']:.1f}%), null {r['hi_null_paired']}, bar {r.get('bar_holds')}" if r['hi_copies'] else "")
              + f"; families with every >= 0.90 copy paired {r['families_all_hi_paired']} / {r['families_with_2_hi']}", flush=True)
        print(f"    by identity {r['by_identity']}", flush=True)
        print(f"    failures {r['failures']}; mainstays {r['mainstays']}", flush=True)
        if sp == "human":
            rn = r["nosat"]
            print(f"    satellite-filtered SDs: >= 0.90 paired {rn['hi_paired']} / {rn['hi_copies']} "
                  f"({100 * rn['hi_paired'] / max(1, rn['hi_copies']):.1f}%)", flush=True)
    json.dump(res, open(out_p, "w"), indent=1)
    return 0


if __name__ == "__main__":
    if sys.argv[1:2] == ["run"] and len(sys.argv) == 3:
        sys.exit(cmd_run(sys.argv[2]))
    raise SystemExit(__doc__)
