#!/usr/bin/env python3
"""Truth of docs/PREREG_o3_maternal_reference_2026-10-08.md section 3 (+ Amendment 1) and the read labels of section 5, for the run's reference
haplotype C.REF (env O3_REF, default mat): a locus is ABSENT iff it exists on C.OTHER and has no counterpart on C.REF.

    O3_REF=mat truth.py loci     # WR/truth/loci.tsv and lrpap1_loci.tsv
    O3_REF=mat truth.py labels   # WR/truth/labels.tsv (needs W/map/reads.<OTHER>.all.bam and WR/truth/loci.tsv)
"""
import collections
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "rna_allele"))
import common as C  # noqa: E402
import extract_reads as E  # noqa: E402


def catalog_loci(path=f"{C.TRUTH}/refabsent/bonly.tsv", other=None):
    """haplotype-only loci of the 378-family catalog lying on the truth haplotype `other` (default C.OTHER): the rows of bonly.tsv with
    hap == other (chromosomes `_pri` took from the REFERENCE haplotype); 0-based half-open on `other`"""
    other = other or C.OTHER
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["hap"] == other:
            out.append(dict(locus=r["locus"], kind="catalog", family=r["family"], chrom=r["chrom"],
                            start=int(r["start"]), end=int(r["end"]), name=r["locus"]))
    return out


def fasta_names(path):
    return [ln[1:].strip() for ln in open(path) if ln[0] == ">"]


def query_of(loci_, queries):
    """cid -> the PAF query ('NC:lo+1-hi') on the same chromosome overlapping the locus the most"""
    out = {}
    for cid, _name, chrom, lo, hi in loci_:
        best = None
        for q in queries:
            c, rng = q.rsplit(":", 1)
            s, e = (int(x) for x in rng.split("-"))
            ov = min(hi, e) - max(lo, s - 1)
            if c == chrom and ov > 0 and (best is None or ov > best[0]):
                best = (ov, q)
        if best:
            out[cid] = best[1]
    return out


def hap_intervals(loci_, chrmap, lift, hits, qof, al):
    """cid -> dict(name, mat=(acc, s, e)|None, pat=(acc, s, e)|None, sex=bool).
    On a chromosome `_pri` took from haplotype H, `_pri` coordinates ARE H coordinates (chrmap seq_check identical) and the other haplotype's
    interval is the lift (B) iff lift_frac >= .5 and a hit of the locus body (identity >= .90, coverage >= .80) overlaps it.
    hits = {'mat': {query: [(acc, s, e, ident, cov)]}, 'pat': {...}}. A chromosome with no B partner in chrmap (chrX, chrY) is `sex`; one missing from chrmap altogether is taken as chrY (pat only)."""
    out = {}
    for cid, name, chrom, lo, hi in loci_:
        d = dict(name=name, mat=None, pat=None, sex=False)
        row, L = chrmap.get(chrom), lift.get(cid)
        if row is None:
            d["pat"], d["sex"] = (al[("pat", "Y")], lo, hi), True
        else:
            h = row["same_hap"]
            o = "mat" if h == "pat" else "pat"
            d[h] = (row["same_name"], lo, hi)
            if row.get("B_name", "present") == "":      # chrX / chrY: no counterpart chromosome on the other haplotype (male: X maternal, Y paternal)
                d["sex"] = True
            elif L is not None and float(L["lift_frac"]) >= 0.5 and any(
                    x[0] == L["B_chrom"] and x[1] < int(L["B_end"]) and int(L["B_start"]) < x[2] for x in hits[o].get(qof.get(cid), [])):
                d[o] = (L["B_chrom"], int(L["B_start"]), int(L["B_end"]))
        out[cid] = d
    return out


def lrpap1_rows(loci_, iv, lift, chrmap, ref, other):
    """loci.tsv rows of the LRPAP1 loci seen from reference `ref`. kind 'lrpap1' / 'sex' = absent from `ref` (no interval there) and present on
    `other`. kind 'lrpap1_desc' (Amendment 2) = DESCRIPTIVE only: a locus on a chromosome `_pri` took from `other` whose truth_lift class is 'T?'
    (moved / structurally uncertain): it has a diverged syntenic counterpart on `ref`, so it is not an absent locus and carries no bar."""
    rows = []
    for cid, name, chrom, _lo, _hi in loci_:
        d = iv[cid]
        if d[other] is None:
            continue
        if d[ref] is None:
            kind = "sex" if d["sex"] else "lrpap1"
        elif chrmap.get(chrom, {}).get("same_hap") == other and lift.get(cid, {}).get("class") == "T?":
            kind = "lrpap1_desc"
        else:
            continue
        rows.append(dict(locus=f"LRPAP1_{cid}", kind=kind, family="LRPAP1", chrom=d[other][0], start=d[other][1], end=d[other][2], name=name))
    return rows


def label_reads(placements, loci_, fam_of, lrp_of, primaries=None):
    """placements: {read: Rec | None} the read's untied best placement on the TRUTH haplotype (accessions); loci_: loci dicts (start/end ints);
    fam_of: {read: family} of the 34-family set; lrp_of: {read: [cid]} of the LRPAP1 net; primaries: {read: Rec} the read's primary record there.
    -> {read: locus id | 'shared' (untied placement elsewhere) | 'ambiguous' (tied or unplaced)}.
    Descriptive loci (kind 'lrpap1_desc', Amendment 2) also take a TIED read of the LRPAP1 net whose primary overlaps them."""
    primaries = primaries or {}

    def on(rec, L):
        return rec is not None and rec.ref == L["chrom"] and rec.start < L["end"] and L["start"] < rec.end
    out = {}
    for n, p in placements.items():
        if p is None:
            hit = next((L["locus"] for L in loci_ if L["kind"] == "lrpap1_desc" and n in lrp_of and on(primaries.get(n), L)), None)
            out[n] = hit or "ambiguous"
            continue
        hit = None
        for L in loci_:
            if not on(p, L):
                continue
            if L["kind"] == "catalog" and fam_of.get(n) != L["family"]:
                continue
            if L["kind"] != "catalog" and n not in lrp_of:
                continue
            hit = L["locus"]
            break
        out[n] = hit or "shared"
    return out


def cmd_loci():
    import truth_lift
    os.makedirs(f"{C.WR}/truth", exist_ok=True)
    al = C.alias()
    L11 = E.loci()
    chrmap = {r["pri"]: r for r in csv.DictReader(open(f"{C.TRUTH}/chrmap.tsv"), delimiter="\t")}
    genes = f"{C.WR}/truth/lrpap1.genes.tsv"
    with open(genes, "w") as o:
        o.write("gene_id\tchrom\tstrand\texons\n")
        for cid, _n, chrom, lo, hi in L11:
            if chrom in chrmap:
                o.write(f"{cid}\t{chrom}\t+\t{lo}-{hi}\n")
    truth_lift.main(["--chrmap", f"{C.TRUTH}/chrmap.tsv", "--paf-dir", f"{C.TRUTH}/out", "--genes", genes,
                     "--out", f"{C.WR}/truth/lrpap1.lift.tsv"])
    lift = {r["gene_id"]: r for r in csv.DictReader(open(f"{C.WR}/truth/lrpap1.lift.tsv"), delimiter="\t")}
    queries = fasta_names(f"{E.LRP}/copies8.fa") + fasta_names(f"{E.LRP}/partial3.fa")
    hits = {"mat": {}, "pat": {}}
    for h in hits:
        for f in (f"copies8.{h}.paf", f"partial3.{h}.paf"):
            for q, v in C.paf_hits(f"{E.LRP}/{f}", al).items():      # asm20 -N 50 -p 0.5 (copies8) as produced on 10-04
                hits[h].setdefault(q, []).extend(v)
    iv = hap_intervals(L11, chrmap, lift, hits, query_of(L11, queries), al)
    with open(f"{C.WR}/truth/lrpap1_loci.tsv", "w") as o:
        o.write("cid\tname\tsex\tpat_acc\tpat_s\tpat_e\tmat_acc\tmat_s\tmat_e\n")
        for cid, d in iv.items():
            p, m = d["pat"] or ("", "", ""), d["mat"] or ("", "", "")
            o.write("\t".join(str(x) for x in [cid, d["name"], int(d["sex"]), *p, *m]) + "\n")
    rows = catalog_loci() + lrpap1_rows(L11, iv, lift, chrmap, C.REF, C.OTHER)
    with open(f"{C.WR}/truth/loci.tsv", "w") as o:
        o.write("locus\tkind\tfamily\tchrom\tstart\tend\tname\n")
        for r in rows:
            o.write("\t".join(str(r[k]) for k in ("locus", "kind", "family", "chrom", "start", "end", "name")) + "\n")
    print(f"reference {C.REF}, truth {C.OTHER}: LRPAP1 loci present on mat {sum(1 for d in iv.values() if d['mat'])}, on pat "
          f"{sum(1 for d in iv.values() if d['pat'])} of {len(iv)}; absent from {C.REF}: {[c for c, d in iv.items() if d[C.REF] is None]}")
    print(f"{C.REF}-absent loci: {len(rows)} ({sum(1 for r in rows if r['kind'] == 'catalog')} catalog, "
          f"{sum(1 for r in rows if r['kind'] == 'lrpap1')} LRPAP1, {sum(1 for r in rows if r['kind'] == 'sex')} sex control, "
          f"{sum(1 for r in rows if r['kind'] == 'lrpap1_desc')} descriptive)")


def cmd_labels():
    al = C.alias()
    loci_ = list(csv.DictReader(open(f"{C.WR}/truth/loci.tsv"), delimiter="\t"))
    for L in loci_:
        L["start"], L["end"] = int(L["start"]), int(L["end"])
    fam_of = {r["read"]: r["family"] for r in csv.DictReader(open(E.R34), delimiter="\t")}
    lrp_of = {r["read"]: r["cids"].split(",") for r in csv.DictReader(open(f"{C.W}/reads/R_LRP.names.tsv"), delimiter="\t")}
    recs = C.read_records(f"{C.W}/map/reads.{C.OTHER}.all.bam", al)
    labels = label_reads({n: C.place(rs) for n, rs in recs.items()}, loci_, fam_of, lrp_of,
                         {n: rs[0] for n, rs in recs.items() if rs and rs[0].primary})
    with open(f"{C.WR}/truth/labels.tsv", "w") as o:
        o.write("read\tlabel\n")
        for n, g in sorted(labels.items()):
            o.write(f"{n}\t{g}\n")
    cnt = collections.Counter(labels.values())
    old = {r["locus"]: int(r["n_reads"]) for r in csv.DictReader(open(f"{C.TRUTH}/refabsent/bonly_expressed.tsv"), delimiter="\t")}
    print(f"reads {len(labels)}: shared {cnt['shared']}, ambiguous {cnt['ambiguous']}")
    print(f"locus\tkind\treads(new, best untied placement on {C.OTHER})\treads(Amendment 10 express, both haplotypes)")
    for L in loci_:
        if cnt[L["locus"]] or old.get(L["locus"]):
            print(f"{L['locus']}\t{L['kind']}\t{cnt[L['locus']]}\t{old.get(L['locus'], '-')}")


if __name__ == "__main__":
    {"loci": cmd_loci, "labels": cmd_labels}[sys.argv[1]]()
