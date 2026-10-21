#!/usr/bin/env python3
"""PREREG_o1_self_hits_2026-10-21: does `mcl_families --reject-self-hits` remove register 1368's self-overlap families without
collateral? Control and fix catalogs come from the same binary, GTF and cached all-vs-all PAF.

    self_hits_eval.py run OUT.json NAME=CONTROL_PREFIX,FIX_PREFIX,SPECIES ...
        PREFIX = <dir>/<sample>.fam (its .copies.tsv, .loci.gff3, .loci.paf); SPECIES human | gorilla
"""
import collections
import csv
import gzip
import importlib.util
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SLOP = 5
spec = importlib.util.spec_from_file_location("o1_sd_pairs_eval", os.path.join(HERE, "sd_pairs.py"))
SP = importlib.util.module_from_spec(spec)
spec.loader.exec_module(SP)
CAT_GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/gencode_chm13/chm13v2.0_CAT_Liftoff.slim.gff3.gz"
MAIN = re.compile(r"^(NPIP[AB]?\d*|TBC1D3[A-Z]?|GOLGA6[A-Z]?|GOLGA8[A-Z]?|SMN[12]|LRPAP1|TSPY\d*|RBMY\d*[A-Z]?\d*|DAZ\d|HSFY\d|BPY2[A-Z]?|CDY\d[A-Z]?)$")
MAIN_DESC = re.compile(r"nuclear pore complex|TBC1 domain family member 3|golgin A[68]|survival of motor neuron|LDL receptor related "
                       r"protein associated protein 1|testis.specific Y|RNA.binding motif protein, Y|deleted in azoospermia", re.I)


def is_self_hit(qk, qs, qe, strand, tk, ts, te, slop=SLOP):
    """mcl_families' rule: one contig, + strand, aligned genomic intervals equal within slop (keys 1-based)."""
    if qk[0] != tk[0] or strand != "+":
        return False
    q0, q1 = qk[1] - 1 + qs, qk[1] - 1 + qe
    t0, t1 = tk[1] - 1 + ts, tk[1] - 1 + te
    return abs(q0 - t0) <= slop and abs(q1 - t1) <= slop


def self_only(best):
    """Copies whose best same-family evidence at >= 0.90 identity is a self-hit, with no non-self record at >= 0.90."""
    return {g for g, d in best.items() if d.get("self", 0.0) >= 0.90 and d.get("non", 0.0) < 0.90}


def b1(self_only_set, fixed_best, fixed_in_family):
    """(fixed, total): a self-only copy is fixed when it left the families or has a non-self same-family record >= 0.90."""
    n = sum(1 for g in self_only_set if g not in fixed_in_family or fixed_best.get(g, {}).get("non", 0.0) >= 0.90)
    return n, len(self_only_set)


def b2(control_paired, fixed_paired):
    """(still paired, control paired)."""
    return len(control_paired & fixed_paired), len(control_paired)


def loci_keys(gff3):
    out = {}
    for ln in open(gff3):
        f = ln.split("\t")
        if len(f) > 8 and f[2] == "gene":
            name = re.search(r"Name=([^;\n]+)", f[8]).group(1)
            out[f"{f[0]}:{f[3]}-{f[4]}"] = name
    return out


def copies(prefix):
    rows = list(csv.DictReader(open(f"{prefix}.copies.tsv"), delimiter="\t"))
    return rows, {r["gene_id"]: r["family_id"] for r in rows}


def best_by_family(paf, key2gene, fams):
    """{arm: {gene: {"self": best identity, "non": best identity}}} over same-family records, for each arm's family map."""
    best = {a: collections.defaultdict(dict) for a in fams}
    parse = lambda n: (n.rsplit(":", 1)[0], int(n.rsplit(":", 1)[1].split("-")[0]), int(n.rsplit(":", 1)[1].split("-")[1]))
    for ln in open(paf):
        f = ln.split("\t", 12)
        if f[0] == f[5]:
            continue
        gq, gt = key2gene.get(f[0]), key2gene.get(f[5])
        if gq is None or gt is None or gq == gt:
            continue
        ident = int(f[9]) / max(1, int(f[10]))
        kind = "self" if is_self_hit(parse(f[0]), int(f[2]), int(f[3]), f[4], parse(f[5]), int(f[7]), int(f[8])) else "non"
        for a, fam in fams.items():
            if gq in fam and fam.get(gq) == fam.get(gt):
                for g in (gq, gt):
                    if ident > best[a][g].get(kind, 0.0):
                        best[a][g][kind] = ident
    return best


def sd_paired_genes(rows, sd):
    cps = [dict(key=r["gene_id"], family=r["family_id"], chrom=r["chrom"],
                exons=[tuple(int(x) for x in e.split("-")) for e in r["exons"].split(",")]) for r in rows]
    return SP.sd_paired(cps, sd)


def mainstay_genes(rows, species, ggo_desc):
    """{gene_id: mainstay label} of copies overlapping a mainstay gene (human: CAT names; gorilla: RefSeq descriptions)."""
    genes = collections.defaultdict(list)
    if species == "human":
        for ln in gzip.open(CAT_GFF, "rt"):
            f = ln.split("\t")
            if len(f) > 8 and f[2] == "gene":
                m = re.search(r"gene_name=([^;\n]+)", f[8]) or re.search(r"Name=([^;\n]+)", f[8])
                if m and MAIN.match(m.group(1)):
                    genes[f[0]].append((int(f[3]) - 1, int(f[4]), m.group(1)))
    else:
        for ch, s, e, sym, desc in ggo_desc:
            if MAIN.match(sym) or MAIN_DESC.search(desc):
                genes[ch].append((s, e, sym))
    out = {}
    for r in rows:
        s, e = int(r["start"]), int(r["end"])
        hit = next((g for gs, ge, g in genes.get(r["chrom"], []) if gs < e and s < ge), None)
        if hit:
            out[r["gene_id"]] = hit
    return out


def cmd_run(out_p, specs):
    sds = {"human": SP.load_sd(SP.CHM13_SD), "gorilla": SP.load_sd(SP.GGO_SD, min_id=0.90)}
    ggo_desc = SP.gene_descriptions(SP.GGO_GFF)
    res = {}
    for spec_ in specs:
        name, rest = spec_.split("=", 1)
        ctl, fix, species = rest.split(",")
        rows_c, fam_c = copies(ctl)
        rows_f, fam_f = copies(fix)
        key2gene = loci_keys(f"{ctl}.loci.gff3")
        best = best_by_family(f"{ctl}.loci.paf", key2gene, {"control": fam_c, "fix": fam_f})
        so = self_only(best["control"])
        n1, d1 = b1(so, best["fix"], set(fam_f))
        pc, pf = sd_paired_genes(rows_c, sds[species]), sd_paired_genes(rows_f, sds[species])
        n2, d2 = b2(pc, pf)
        hi = lambda rows: [r for r in rows if r["max_family_identity"] not in ("", "NA") and float(r["max_family_identity"]) >= 0.90]
        hc, hf = hi(rows_c), hi(rows_f)
        mc, mf = mainstay_genes(rows_c, species, ggo_desc), mainstay_genes(rows_f, species, ggo_desc)
        size = lambda fam: max(collections.Counter(fam.values()).values())
        r = dict(self_only_control=d1, b1_fixed=n1, B1=d1 > 0 and n1 >= 0.9 * d1, sd_paired_control=d2, b2_kept=n2,
                 B2=d2 > 0 and n2 >= 0.98 * d2, copies=(len(rows_c), len(rows_f)), families=(len(set(fam_c.values())), len(set(fam_f.values()))),
                 largest_family=(size(fam_c), size(fam_f)),
                 hi_sd_paired=(sum(r_["gene_id"] in pc for r_ in hc), len(hc), sum(r_["gene_id"] in pf for r_ in hf), len(hf)),
                 mainstay_copies=(len(mc), len(mf)), mainstays_lost=sorted(set(mc) - set(mf)))
        r["verdict"] = "PASS" if r["B1"] and r["B2"] else "FAIL"
        res[name] = r
        print(f"[self_hits_eval] {name}: B1 {n1}/{d1} self-only copies fixed ({100 * n1 / max(1, d1):.1f}%) -> {r['B1']}; "
              f"B2 {n2}/{d2} SD-paired copies kept ({100 * n2 / max(1, d2):.1f}%) -> {r['B2']}; copies {r['copies']}, families "
              f"{r['families']}, largest {r['largest_family']}; >= 0.90 SD-paired control {r['hi_sd_paired'][0]}/{r['hi_sd_paired'][1]} "
              f"-> fix {r['hi_sd_paired'][2]}/{r['hi_sd_paired'][3]}; mainstay copies {r['mainstay_copies']}; {r['verdict']}", flush=True)
        if r["mainstays_lost"]:
            print(f"    mainstay copies no longer in a family: {r['mainstays_lost'][:20]}", flush=True)
    res["overall"] = "PASS" if all(v["verdict"] == "PASS" for k, v in res.items() if isinstance(v, dict)) else "FAIL"
    print(f"[self_hits_eval] overall {res['overall']}", flush=True)
    json.dump(res, open(out_p, "w"), indent=1)
    return 0


if __name__ == "__main__":
    if sys.argv[1:2] == ["run"] and len(sys.argv) >= 4:
        sys.exit(cmd_run(sys.argv[2], sys.argv[3:]))
    raise SystemExit(__doc__)
