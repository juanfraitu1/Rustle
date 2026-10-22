#!/usr/bin/env python3
"""PREREG_o1_self_hits_2026-10-21 (evaluation of Amendment 3): does `mcl_families --reject-self-hits` remove register 1368's self-overlap families without
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
    """(still paired, control paired), joined on gene_id."""
    return len(control_paired & fixed_paired), len(control_paired)


def copy_dicts(rows):
    return [dict(key=r["gene_id"], family=r["family_id"], chrom=r["chrom"], strand=r["strand"], start=int(r["start"]), end=int(r["end"]),
                 exons=[tuple(int(x) for x in e.split("-")) for e in r["exons"].split(",")]) for r in rows]


def exon_mates(idx, cps, c):
    """Keys of the copies in `cps` (indexed by idx) sharing exon bases with copy c, on either strand."""
    out = set()
    for s, e in c["exons"]:
        out |= {cps[i]["key"] for i in SP._hits(idx, c["chrom"], s, e)}
    return out


def b2_locus(ctl_paired, ctl_by_key, fix_paired, fix_cps):
    """(kept, control paired) by locus (Amendment 3): a control SD-paired copy is kept when it, or a fix copy sharing exon bases
    with it, is SD-paired in its fix family; a representative switch is not a loss."""
    idx = SP._index(fix_cps)
    kept = sum(1 for k in ctl_paired if k in fix_paired or exon_mates(idx, fix_cps, ctl_by_key[k]) & fix_paired)
    return kept, len(ctl_paired)


def cross_family_exon_pairs(cps):
    """(same strand, opposite strand) counts of copy pairs in different families that share exon bases."""
    idx = SP._index(cps)
    pairs = set()
    for i, c in enumerate(cps):
        for s, e in c["exons"]:
            pairs |= {(min(i, j), max(i, j)) for j in SP._hits(idx, c["chrom"], s, e) if cps[j]["family"] != c["family"]}
    same = sum(1 for i, j in pairs if cps[i]["strand"] == cps[j]["strand"])
    return same, len(pairs) - same


def g1(ctl_same, fix_same):
    """G1: same-strand cross-family exon-sharing pairs grow by at most max(5, 5% of the control's, rounded up)."""
    return fix_same <= ctl_same + max(5, (5 * ctl_same + 99) // 100)


def span_disjoint(a, b):
    """Loci keys (chrom, start, end), 1-based closed: True when they share no base."""
    return a[0] != b[0] or a[2] < b[1] or b[2] < a[1]


def g2(pairs, ctl_by_key, fix_cps):
    """(kept, total) control copies with a distant same-family partner (a non-self record >= 0.90 between span-disjoint loci):
    a pair is kept when fix copies at its two loci (the copies themselves or copies sharing exon bases with them) share a family."""
    idx = SP._index(fix_cps)
    fam = {c["key"]: c["family"] for c in fix_cps}
    mates = {}

    def m(k):
        if k not in mates:
            mates[k] = exon_mates(idx, fix_cps, ctl_by_key[k])
        return mates[k]

    kept, total = set(), set()
    for a, b in pairs:
        total |= {a, b}
        if any(x != y and fam[x] == fam[y] for x in m(a) for y in m(b)):
            kept |= {a, b}
    return len(kept), len(total)


def g3(fix_self_only, ctl_self_only):
    """G3: self-only copies left in the fix catalog, at most max(2, 10% of the control's)."""
    return fix_self_only <= max(2, ctl_self_only / 10)


def verdict(b1_counts, b2_ok, g1_ok, g2_ok, g3_ok):
    """Amendment 3: a failed bar or guard FAILS; fewer than 10 control self-only copies cannot PASS (UNDERPOWERED)."""
    n, d = b1_counts
    if not (b2_ok and g1_ok and g2_ok and g3_ok) or (d >= 10 and n < 0.9 * d):
        return "FAIL"
    return "UNDERPOWERED" if d < 10 else "PASS"


def comember_change(fam_c, fam_f):
    """Among copies in both catalogs: (new co-member pairs, lost co-member pairs, fix families holding >= 2 control families)."""
    common = sorted(g for g in fam_c if g in fam_f)

    def pairs(fam):
        by = collections.defaultdict(list)
        for g in common:
            by[fam[g]].append(g)
        return {(x, y) for v in by.values() for i, x in enumerate(v) for y in v[i + 1:]}

    pc, pf = pairs(fam_c), pairs(fam_f)
    merged = collections.defaultdict(set)
    for g in common:
        merged[fam_f[g]].add(fam_c[g])
    return len(pf - pc), len(pc - pf), sum(1 for v in merged.values() if len(v) >= 2)


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


def best_by_family(paf, key2gene, fams, pair_arm="control"):
    """({arm: {gene: {"self": best identity, "non": best identity}}} over same-family records, for each arm's family map;
    {(gene, gene)}: pair_arm's same-family copies at span-disjoint loci joined by a non-self record at >= 0.90 identity)."""
    best = {a: collections.defaultdict(dict) for a in fams}
    pairs = set()
    parse = lambda n: (n.rsplit(":", 1)[0], int(n.rsplit(":", 1)[1].split("-")[0]), int(n.rsplit(":", 1)[1].split("-")[1]))
    for ln in open(paf):
        f = ln.split("\t", 12)
        if f[0] == f[5]:
            continue
        gq, gt = key2gene.get(f[0]), key2gene.get(f[5])
        if gq is None or gt is None or gq == gt:
            continue
        ident = int(f[9]) / max(1, int(f[10]))
        q, t = parse(f[0]), parse(f[5])
        kind = "self" if is_self_hit(q, int(f[2]), int(f[3]), f[4], t, int(f[7]), int(f[8])) else "non"
        for a, fam in fams.items():
            if gq in fam and fam.get(gq) == fam.get(gt):
                for g in (gq, gt):
                    if ident > best[a][g].get(kind, 0.0):
                        best[a][g][kind] = ident
                if a == pair_arm and kind == "non" and ident >= 0.90 and span_disjoint(q, t):
                    pairs.add((min(gq, gt), max(gq, gt)))
    return best, pairs


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
        cps_c, cps_f = copy_dicts(rows_c), copy_dicts(rows_f)
        by_key_c = {c["key"]: c for c in cps_c}
        key2gene = loci_keys(f"{ctl}.loci.gff3")
        best, pairs = best_by_family(f"{ctl}.loci.paf", key2gene, {"control": fam_c, "fix": fam_f})
        so = self_only(best["control"])
        n1, d1 = b1(so, best["fix"], set(fam_f))
        left = sum(1 for g in so if g not in fam_f)
        pc, pf = SP.sd_paired(cps_c, sds[species]), SP.sd_paired(cps_f, sds[species])
        n2g = b2(pc, pf)[0]
        n2, d2 = b2_locus(pc, by_key_c, pf, cps_f)
        (cs, co), (fs, fo) = cross_family_exon_pairs(cps_c), cross_family_exon_pairs(cps_f)
        k2, t2 = g2(pairs, by_key_c, cps_f)
        fso = len(self_only(best["fix"]))
        new_p, lost_p, merged = comember_change(fam_c, fam_f)
        members = collections.defaultdict(list)
        for g, f in fam_c.items():
            members[f].append(g)
        dissolved = [sorted(v) for v in members.values() if len(v) == 2 and all(g not in fam_f for g in v)]
        dissolved_ns = [v for v in dissolved if not any(g in so for g in v)]
        hi = lambda rows: [r for r in rows if r["max_family_identity"] not in ("", "NA") and float(r["max_family_identity"]) >= 0.90]
        hc, hf = hi(rows_c), hi(rows_f)
        mc, mf = mainstay_genes(rows_c, species, ggo_desc), mainstay_genes(rows_f, species, ggo_desc)
        size = lambda fam: max(collections.Counter(fam.values()).values())
        r = dict(self_only_control=d1, b1_fixed=n1, b1_left_families=left, b1_repaired=n1 - left, B1=d1 > 0 and n1 >= 0.9 * d1,
                 sd_paired_control=d2, b2_kept_locus=n2, b2_kept_gene_id=n2g, B2=d2 > 0 and n2 >= 0.98 * d2,
                 cross_family_same_strand=(cs, fs), cross_family_opposite_strand=(co, fo), G1=g1(cs, fs),
                 distant_paired_copies=t2, distant_kept=k2, G2=k2 >= 0.98 * t2, self_only_fix=fso, G3=g3(fso, d1),
                 comember_pairs_new=new_p, comember_pairs_lost=lost_p, fix_families_merging_control_families=merged,
                 dissolved_two_copy_families=len(dissolved), dissolved_not_self_only=dissolved_ns[:20],
                 copies=(len(rows_c), len(rows_f)), families=(len(set(fam_c.values())), len(set(fam_f.values()))),
                 largest_family=(size(fam_c), size(fam_f)),
                 hi_sd_paired=(sum(r_["gene_id"] in pc for r_ in hc), len(hc), sum(r_["gene_id"] in pf for r_ in hf), len(hf)),
                 mainstay_copies=(len(mc), len(mf)), mainstays_lost=sorted(set(mc) - set(mf)))
        r["verdict"] = verdict((n1, d1), r["B2"], r["G1"], r["G2"], r["G3"])
        res[name] = r
        print(f"[self_hits_eval] {name}: B1 {n1}/{d1} self-only copies fixed ({left} left the families, {n1 - left} repaired) -> "
              f"{r['B1']}; B2 by locus {n2}/{d2} SD-paired copies kept (gene_id {n2g}) -> {r['B2']}; G1 same-strand cross-family "
              f"exon-sharing pairs {cs} -> {fs} (opposite strand {co} -> {fo}) -> {r['G1']}; G2 distant pairs {k2}/{t2} copies kept -> "
              f"{r['G2']}; G3 self-only copies in the fix {fso} -> {r['G3']}; co-member pairs +{new_p}/-{lost_p}, {merged} fix families "
              f"merge control families, {len(dissolved)} two-copy families dissolved ({len(dissolved_ns)} not self-only); copies "
              f"{r['copies']}, families {r['families']}, largest {r['largest_family']}; >= 0.90 SD-paired control {r['hi_sd_paired'][0]}/"
              f"{r['hi_sd_paired'][1]} -> fix {r['hi_sd_paired'][2]}/{r['hi_sd_paired'][3]}; mainstay copies {r['mainstay_copies']}; "
              f"{r['verdict']}", flush=True)
        if r["mainstays_lost"]:
            print(f"    mainstay copies no longer in a family: {r['mainstays_lost'][:20]}", flush=True)
    vs = [v["verdict"] for v in res.values()]
    res["overall"] = "FAIL" if "FAIL" in vs else "UNDERPOWERED" if "UNDERPOWERED" in vs else "PASS"
    print(f"[self_hits_eval] overall {res['overall']}", flush=True)
    json.dump(res, open(out_p, "w"), indent=1)
    return 0


if __name__ == "__main__":
    if sys.argv[1:2] == ["run"] and len(sys.argv) >= 4:
        sys.exit(cmd_run(sys.argv[2], sys.argv[3:]))
    raise SystemExit(__doc__)
