#!/usr/bin/env python3
"""CAT re-run step 5/6 (docs/CAT_RERUN_STEP5_PLAN_2026-10-01.md): the CAT/Liftoff v2.0 inputs of the NPIP/TBC1D3 layer-order
and nested-lattice study, written into a COPY of its results tree in the shapes the RefSeq code reads.

Subcommands (all write under --root, a copy of layer_order/npip_tbc1d3; nothing outside it is written):
  tables    light/work/refseq/{genes,exons,cds,gene_dbxref}.tsv for CAT genes (plan rows 1-4), the label table
            light/work/refseq/cat_labels.tsv, the RefSeq-name -> CAT-label map light/work/refseq/name_map.tsv (R1, Amendments 1
            and 2), the CAT expression GFF light/work/refseq/cat_expr.gff3.gz (row 22), and the heavy layer's gene/exon tables
            heavy/work/refseq_{genes,exons}.tsv + heavy/seed_members.tsv (row 20)
  members   light/members.tsv (R1 image of the 46 RefSeq corrected members, row 7), light/work/refseq/members_beside_R2.tsv,
            light/work/refseq/readthrough_over_members.txt (row 8), light/work/refseq/lit_truth.cat.tsv (row 18)
  catalogs  light/work/catalogs/{c15_17_22,c16_19_20}/nodes.tsv(.names.tsv) from the reused CAT gene-span catalogs (row 15)

Labels: a CAT gene's label is its gene_name when that name is unique among the CAT genes, else 'gene_name~<CAT gene id>';
gene_id = 'gene-' + label. Labels are display names; RefSeq -> CAT goes through R1 only (never by name).

R1 (protocol): a RefSeq record maps to the CAT gene sharing the most exonic bases on the same strand (ties: Jaccard) =
winloci_data/gencode_chm13/chm13v2.0_CAT_Liftoff.refseq_map.tsv. Amendment 1: a RefSeq record WITHOUT exon lines maps by the
largest span Jaccard on the same strand. Amendment 2: a RefSeq read-through A-B that a truth names as the B copy maps by R1
applied to its exons outside A's RefSeq record (here: PKD1P6-NPIPP1 -> CHM13_G0020725).
"""
import argparse
import bisect
import collections
import csv
import gzip
import os
import re
import subprocess

W = "/mnt/linuxdisk/home/juanfraitu"
CATD = f"{W}/winloci_data/gencode_chm13"
CAT_GENES = f"{CATD}/chm13v2.0_CAT_Liftoff.genes.tsv"
CAT_SLIM = f"{CATD}/chm13v2.0_CAT_Liftoff.slim.gff3.gz"
CAT_FULL = f"{CATD}/chm13v2.0_gencode.gff3"
REFSEQ_MAP = f"{CATD}/chm13v2.0_CAT_Liftoff.refseq_map.tsv"
ROOT_RS = f"{W}/layer_order/npip_tbc1d3"  # the frozen RefSeq results (read only)
HGNC = f"{W}/winloci_data/hgnc/hgnc_complete_set.txt"
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
STEP1_NPIP = os.path.join(REPO, "docs", "lit_subclusters_npip_dishuck_check.CAT.tsv")
STEP1_TBC = os.path.join(REPO, "docs", "tbc1d3_members.CAT.tsv")
LIT_RS = os.path.join(REPO, "docs", "lit_subclusters_npip_tbc1d3_truth.tsv")
# Amendment 2: the RefSeq read-through records a truth names as their second (B) copy
AMEND2_B_COPY = {"PKD1P6-NPIPP1": "PKD1P6"}  # read-through name -> its A part (the RefSeq record whose span is excluded)
AMEND2_EXPECT = {"PKD1P6-NPIPP1": "CHM13_G0020725"}


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def ov_bp(a, b):
    i = j = t = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo < hi:
            t += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return t


def blen(b):
    return sum(e - s for s, e in b)


def blocks(s):
    return [tuple(map(int, x.split("-"))) for x in s.split(",") if x]


def fmt_blocks(b):
    return ",".join(f"{s}-{e}" for s, e in b)


def attrs(col):
    return dict(x.split("=", 1) for x in col.rstrip().split(";") if "=" in x)


# ------------------------------------------------------------------------------------------------ CAT genes and labels
def load_cat():
    g = {}
    for r in tsv(CAT_GENES):
        g[r["gene_id"]] = dict(id=r["gene_id"], name=r["gene_name"], biotype=r["gene_biotype"], source=r["source"],
                               chrom=r["chrom"], start0=int(r["start0"]), end=int(r["end"]), strand=r["strand"],
                               blocks=blocks(r["exon_blocks"]))
    return g


def make_labels(cat):
    n = collections.Counter(g["name"] for g in cat.values())
    return {i: (g["name"] if n[g["name"]] == 1 else f"{g['name']}~{i}") for i, g in cat.items()}


class CatIndex:
    """per (chrom, strand): CAT genes sorted by start, for exon-overlap and span-overlap queries."""

    def __init__(self, cat):
        self.cat = cat
        self.by = collections.defaultdict(list)
        for i, g in cat.items():
            self.by[(g["chrom"], g["strand"])].append((g["start0"], g["end"], i))
        self.starts, self.ml = {}, {}
        for k, v in self.by.items():
            v.sort()
            self.starts[k] = [x[0] for x in v]
            self.ml[k] = max(e - s for s, e, _ in v)

    def near(self, chrom, strand, s, e):
        k = (chrom, strand)
        v = self.by.get(k, [])
        lo = bisect.bisect_left(self.starts.get(k, []), s - self.ml.get(k, 0))
        for a, b, i in v[lo:]:
            if a >= e:
                break
            if b > s:
                yield i

    def r1_exons(self, chrom, strand, ex):
        """R1 on an exon union: most shared exonic bp (ties: Jaccard). Returns (cat id, shared, jaccard) or None."""
        best = None
        for i in self.near(chrom, strand, ex[0][0], ex[-1][1]):
            sh = ov_bp(ex, self.cat[i]["blocks"])
            if sh:
                cand = (sh, sh / (blen(ex) + blen(self.cat[i]["blocks"]) - sh), i)
                if best is None or cand[:2] > best[:2]:
                    best = cand
        return None if best is None else (best[2], best[0], best[1])

    def a1_span(self, chrom, strand, s, e):
        """Amendment 1: largest span Jaccard on the same strand. Returns (cat id, overlap, jaccard) or None."""
        best = None
        for i in self.near(chrom, strand, s, e):
            g = self.cat[i]
            o = min(e, g["end"]) - max(s, g["start0"])
            if o > 0:
                cand = (o / (max(e, g["end"]) - min(s, g["start0"])), o, i)
                if best is None or cand[:2] > best[:2]:
                    best = cand
        return None if best is None else (best[2], best[1], best[0])


def refseq_tables():
    """RefSeq genes (frozen ROOT tables): id -> dict(name, chrom, start0, end, strand, description, exonless, exons)."""
    rs = {}
    for r in tsv(f"{ROOT_RS}/light/work/refseq/genes.tsv"):
        rs[r["gene_id"]] = dict(id=r["gene_id"], name=r["name"], biotype=r["biotype"], description=r["description"],
                                chrom=r["chrom"], start0=int(r["start0"]), end=int(r["end"]), strand=r["strand"])
    src = {}
    for r in tsv(f"{ROOT_RS}/heavy/work/refseq_genes.tsv"):
        src["gene-" + r["gene_id"]] = r["exon_source"]
    ex = collections.defaultdict(list)
    with open(f"{ROOT_RS}/heavy/work/refseq_exons.tsv") as fh:
        next(fh)
        for line in fh:
            g, c, s, e = line.rstrip("\n").split("\t")
            ex["gene-" + g].append((int(s), int(e)))
    for i, r in rs.items():
        r["exonless"] = src.get(i) == "gene_body"
        r["exons"] = merge(ex.get(i, [])) if not r["exonless"] else []
    return rs


def r1_image(rs, cat, cidx, rmap):
    """RefSeq gene id -> (CAT id or '', rule, shared bp, jaccard) for every RefSeq gene (R1, Amendments 1 and 2)."""
    by_name = {r["name"]: i for i, r in rs.items()}
    out = {}
    for i, r in rs.items():
        if r["name"] in AMEND2_B_COPY:
            a = rs[by_name[AMEND2_B_COPY[r["name"]]]]
            rest = []
            for s, e in r["exons"]:
                for x, y in ((s, min(e, a["start0"])), (max(s, a["end"]), e)):
                    if y > x:
                        rest.append((x, y))
            hit = cidx.r1_exons(r["chrom"], r["strand"], merge(rest)) if rest else None
            out[i] = (hit[0], "amendment2_exons_outside_" + a["name"], hit[1], hit[2]) if hit else ("", "none", 0, 0.0)
            assert out[i][0] == AMEND2_EXPECT[r["name"]], (r["name"], out[i])
            continue
        if r["exonless"]:
            hit = cidx.a1_span(r["chrom"], r["strand"], r["start0"], r["end"])
            out[i] = (hit[0], "amendment1_span_jaccard", hit[1], hit[2]) if hit else ("", "none", 0, 0.0)
            continue
        m = rmap.get(i)
        if m and m["cat_gene"]:
            out[i] = (m["cat_gene"], "R1_exons", int(m["shared_bp"]), float(m["jaccard"]))
        else:
            out[i] = ("", "none", 0, 0.0)
    return out


def load_rmap():
    return {r["refseq_gene"]: r for r in tsv(REFSEQ_MAP)}


# ------------------------------------------------------------------------------------------------ full GFF3 stream
def stream_full():
    """chm13v2.0_gencode.gff3 (streamed): gene -> source_gene; transcript -> gene; transcript -> CDS segments."""
    src, tx_gene, cds = {}, {}, collections.defaultdict(list)
    with open(CAT_FULL) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9:
                continue
            t = f[2]
            if t == "gene":
                a = attrs(f[8])
                src[a.get("gene_id", a.get("ID"))] = a.get("source_gene", "")
            elif t == "transcript":
                a = attrs(f[8])
                tx_gene[a["ID"]] = a["Parent"]
            elif t == "CDS":
                a = attrs(f[8])
                cds[a["Parent"]].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] in "012" else 0))
    return src, tx_gene, cds


def longest(cds_by_tx):
    """transcript -> [(start0, end, phase)] -> the longest CDS's segments (ties: first transcript id); = annotation_nodes /
    truth.longest, the rule refseq_tables.py used."""
    best = None
    for t in sorted(cds_by_tx):
        segs = cds_by_tx[t]
        L = sum(e - s for s, e, _ in segs)
        if best is None or L > best[0]:
            best = (L, sorted(segs))
    return best[1] if best else None


def readthrough_ids(cat, labels, rs, img):
    """R3: CAT gene_name 'A-B' with both A and B CAT gene names, or the R1 image of a RefSeq read-through."""
    names = {g["name"] for g in cat.values()}
    out = {}
    for i, g in cat.items():
        parts = g["name"].split("-")
        if len(parts) == 2 and parts[0] in names and parts[1] in names:
            out[i] = f"name {parts[0]}-{parts[1]}"
    for rid, r in rs.items():
        if "readthrough" in r["description"].lower() and img[rid][0]:
            c = img[rid][0]
            out[c] = (out[c] + "; " if c in out else "") + f"R1 image of RefSeq read-through {r['name']}"
    return out


# ------------------------------------------------------------------------------------------------ tables
def cmd_tables(a):
    root = a.root
    cat = load_cat()
    labels = make_labels(cat)
    cidx = CatIndex(cat)
    rs = refseq_tables()
    rmap = load_rmap()
    img = r1_image(rs, cat, cidx, rmap)
    rt = readthrough_ids(cat, labels, rs, img)
    print(f"[cat] genes {len(cat)}; names repeated {sum(1 for i in cat if '~' in labels[i])} genes get 'name~id' labels; "
          f"R3 read-throughs {len(rt)}")
    src, tx_gene, cds = stream_full()
    print(f"[full gff3] genes {len(src)}; transcripts {len(tx_gene)}; transcripts with CDS {len(cds)}")
    # HGNC by Ensembl id (R4)
    ens2hgnc = {}
    for r in tsv(HGNC):
        for e in r["ensembl_gene_id"].split("|"):
            if e:
                ens2hgnc.setdefault(e, r["hgnc_id"])
    d = f"{root}/light/work/refseq"
    os.makedirs(d, exist_ok=True)
    order = sorted(cat, key=lambda i: (cat[i]["chrom"], cat[i]["start0"], i))
    with open(f"{d}/genes.tsv", "w") as fg, open(f"{d}/exons.tsv", "w") as fe, open(f"{d}/cat_labels.tsv", "w") as fl, \
            open(f"{d}/gene_dbxref.tsv", "w") as fx:
        fg.write("gene_id\tname\tbiotype\tdescription\tchrom\tstart0\tend\tstrand\tcat_gene_id\tcat_gene_name\tcat_source\t"
                 "source_gene\n")
        fe.write("gene_id\tname\tchrom\tstart0\tend\tstrand\tn_blocks\texons\n")
        fl.write("label\tgene_id\tcat_gene_id\tcat_gene_name\tbiotype\tsource\tchrom\tstart0\tend\tstrand\tsource_gene\t"
                 "hgnc_id\treadthrough_R3\n")
        n_h = 0
        for i in order:
            g, lab = cat[i], labels[i]
            sg = src.get(i, "")
            ens = sg.split(".")[0] if sg.startswith("ENSG") else ""
            h = ens2hgnc.get(ens, "")
            n_h += bool(h)
            desc = f"CAT {i} {g['name']} ({g['source']})" + (f"; readthrough (R3: {rt[i]})" if i in rt else "")
            fg.write(f"gene-{lab}\t{lab}\t{g['biotype']}\t{desc}\t{g['chrom']}\t{g['start0']}\t{g['end']}\t{g['strand']}\t"
                     f"{i}\t{g['name']}\t{g['source']}\t{sg}\n")
            fe.write(f"gene-{lab}\t{lab}\t{g['chrom']}\t{g['start0']}\t{g['end']}\t{g['strand']}\t{len(g['blocks'])}\t"
                     f"{fmt_blocks(g['blocks'])}\n")
            fl.write(f"{lab}\tgene-{lab}\t{i}\t{g['name']}\t{g['biotype']}\t{g['source']}\t{g['chrom']}\t{g['start0']}\t"
                     f"{g['end']}\t{g['strand']}\t{sg}\t{h}\t{rt.get(i, '')}\n")
            fx.write(f"gene-{lab}\t{('ENSG:' + ens + ',HGNC:' + h) if h else ('ENSG:' + ens if ens else '')}\n")
    print(f"[R4] CAT genes with an HGNC id via source_gene -> ensembl_gene_id: {n_h} of {len(cat)}")
    # cds (one protein per gene: longest-CDS transcript)
    gcds = collections.defaultdict(dict)
    for t, segs in cds.items():
        gcds[tx_gene[t]][t] = segs
    n_cds = 0
    with open(f"{d}/cds.tsv", "w") as fc:
        fc.write("gene_id\tname\tbiotype\tchrom\tstrand\tcds\n")
        for i in sorted(gcds, key=lambda i: (cat[i]["chrom"], labels[i])):
            segs = longest(gcds[i])
            g = cat[i]
            fc.write(f"gene-{labels[i]}\t{labels[i]}\t{g['biotype']}\t{g['chrom']}\t{g['strand']}\t"
                     f"{','.join(f'{x}-{y}:{p}' for x, y, p in segs)}\n")
            n_cds += 1
    print(f"[cds] genes with CDS {n_cds}")
    # name map: RefSeq name -> CAT label (R1 / A1 / A2)
    with open(f"{d}/name_map.tsv", "w") as fm:
        fm.write("refseq_gene_id\trefseq_name\tcat_gene_id\tcat_label\tcat_gene_name\trule\tshared_bp\tjaccard\n")
        for rid in sorted(rs, key=lambda x: (rs[x]["chrom"], rs[x]["start0"], x)):
            c, rule, sh, jac = img[rid]
            fm.write(f"{rid}\t{rs[rid]['name']}\t{c}\t{labels.get(c, '')}\t{cat[c]['name'] if c else ''}\t{rule}\t{sh}\t"
                     f"{jac:.4f}\n")
    print(f"[name map] RefSeq genes {len(img)}; with a CAT image {sum(1 for v in img.values() if v[0])}; rules "
          f"{dict(collections.Counter(v[1].split('_exons_outside')[0] for v in img.values()))}")
    # expression GFF: gene ID=gene-<label>, transcripts and exons from the slim file
    eg = f"{d}/cat_expr.gff3"
    n = 0
    with gzip.open(CAT_SLIM, "rt") as fh, open(eg, "w") as out:
        out.write("##gff-version 3\n")
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            a = attrs(f[8])
            if f[2] == "gene":
                f[8] = f"ID=gene-{labels[a['ID']]};Name={labels[a['ID']]}"
            elif f[2] == "transcript":
                f[8] = f"ID={a['ID']};Parent=gene-{labels[a['Parent']]}"
            elif f[2] == "exon":
                f[8] = f"Parent={a['Parent']}"
            out.write("\t".join(f) + "\n")
            n += 1
    subprocess.run(["gzip", "-f", eg], check=True)
    print(f"[expr gff] {n} lines -> {eg}.gz")
    # heavy layer tables (S2 / EXPR): gene_id = label (heavy strips 'gene-'), distinct exon intervals per gene
    hx = collections.defaultdict(set)
    tx_of = {}
    with gzip.open(CAT_SLIM, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            a = attrs(f[8])
            if f[2] == "transcript":
                tx_of[a["ID"]] = a["Parent"]
            elif f[2] == "exon":
                hx[tx_of[a["Parent"]]].add((f[0], int(f[3]) - 1, int(f[4])))
    hw = f"{root}/heavy/work"
    os.makedirs(hw, exist_ok=True)
    with open(f"{hw}/refseq_genes.tsv", "w") as fg, open(f"{hw}/refseq_exons.tsv", "w") as fe:
        fg.write("gene_id\tname\tgeneid\tchrom\tstart\tend\tstrand\tbiotype\tfeature_type\tn_exons\texon_source\n")
        fe.write("gene_id\tchrom\tstart\tend\n")
        for i in order:
            g = cat[i]
            fg.write(f"{labels[i]}\t{labels[i]}\t{i}\t{g['chrom']}\t{g['start0']}\t{g['end']}\t{g['strand']}\t"
                     f"{g['biotype']}\tgene\t{len(hx[i])}\tgff_exon\n")
        for i in sorted(hx, key=lambda i: labels[i]):
            for c, s, e in sorted(hx[i]):
                fe.write(f"{labels[i]}\t{c}\t{s}\t{e}\n")
    # seeds: R1 images of the 39 RefSeq seeds (R2)
    seeds = tsv(f"{ROOT_RS}/heavy/seed_members.tsv")
    rows, drops, merged = [], [], collections.defaultdict(list)
    for r in seeds:
        c = img["gene-" + r["gene_id"]][0]
        if not c:
            drops.append(r["name"])
            continue
        merged[c].append(r["name"])
    for c, names in sorted(merged.items(), key=lambda kv: (cat[kv[0]]["chrom"], cat[kv[0]]["start0"])):
        fam = {r["family"] for r in seeds if r["name"] in names}
        assert len(fam) == 1, (c, names)
        g = cat[c]
        ex = sorted(hx[c])
        rows.append(f"{labels[c]}\t{labels[c]}\t{fam.pop()}\t{c}\t{g['chrom']}\t{g['start0']}\t{g['end']}\t{g['strand']}\t"
                    f"{g['biotype']}\t{len(ex)}\tgff_exon\t{';'.join(f'{s}-{e}' for _c, s, e in ex)}\t{','.join(names)}\n")
    with open(f"{root}/heavy/seed_members.tsv", "w") as fh:
        fh.write("gene_id\tname\tfamily\tgeneid\tchrom\tstart\tend\tstrand\tbiotype\tn_exons\texon_source\texons\t"
                 "refseq_seeds\n")
        fh.writelines(rows)
    print(f"[seeds] RefSeq seeds {len(seeds)} -> CAT seeds {len(rows)}; dropped (no CAT image) {drops}; merged "
          f"{[v for v in merged.values() if len(v) > 1]}")


# ------------------------------------------------------------------------------------------------ members, R3, literature
def cmd_members(a):
    root = a.root
    cat = load_cat()
    labels = make_labels(cat)
    cidx = CatIndex(cat)
    rs = refseq_tables()
    img = r1_image(rs, cat, cidx, load_rmap())
    d = f"{root}/light/work/refseq"
    mem = tsv(f"{ROOT_RS}/light/members.corrected.tsv")
    by_cat = collections.defaultdict(list)
    drops = []
    for r in mem:
        c, rule, sh, jac = img[r["gene_id"]]
        if not c:
            drops.append(r["name"])
            continue
        by_cat[c].append((r, rule, sh, jac))
    rows = []
    for c, lst in by_cat.items():
        fams = {x[0]["family"] for x in lst}
        assert len(fams) == 1, (c, fams)
        g = cat[c]
        rows.append({"gene_id": f"gene-{labels[c]}", "name": labels[c], "biotype": g["biotype"], "chrom": g["chrom"],
                     "start": g["start0"], "end": g["end"], "strand": g["strand"], "family": fams.pop(),
                     "member_basis": "R1 image of RefSeq member " + ",".join(x[0]["name"] for x in lst)
                                     + " (" + ",".join(x[1] for x in lst) + ")",
                     "cat_gene_id": c, "cat_gene_name": g["name"], "refseq_members": ",".join(x[0]["name"] for x in lst),
                     "refseq_member_basis": " | ".join(x[0]["member_basis"] for x in lst),
                     "rule": ",".join(x[1] for x in lst), "shared_bp": ",".join(str(x[2]) for x in lst),
                     "jaccard": ",".join(f"{x[3]:.3f}" for x in lst)})
    rows.sort(key=lambda x: (x["family"], x["chrom"], int(x["start"])))
    cols = list(rows[0].keys())
    with open(f"{root}/light/members.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for x in rows:
            fh.write("\t".join(str(x[k]) for k in cols) + "\n")
    print(f"[members] RefSeq corrected members {len(mem)} -> CAT members {len(rows)} (NPIP "
          f"{sum(x['family'] == 'NPIP' for x in rows)}, TBC1D3 {sum(x['family'] == 'TBC1D3' for x in rows)}); dropped (no CAT "
          f"gene) {drops}; merged {[x['refseq_members'] for x in rows if ',' in x['refseq_members']]}")
    for x in rows:
        print(f"   {x['family']:6s} {x['refseq_members']:28s} -> {x['cat_gene_id']} {x['name']} ({x['rule']}, shared "
              f"{x['shared_bp']}, J {x['jaccard']})")
    memb = {x["cat_gene_id"] for x in rows}
    # other CAT records overlapping a member's RefSeq record (not copies; listed, R1)
    with open(f"{d}/members_overlapping_cat_records.tsv", "w") as fh:
        fh.write("refseq_member\tcat_member\tother_cat_records_overlapping_the_refseq_record\n")
        for r in mem:
            c = img[r["gene_id"]][0]
            rr = rs[r["gene_id"]]
            other = [f"{i}:{labels[i]}" for i in sorted(set(cidx.near(rr["chrom"], "+", rr["start0"], rr["end"]))
                                                          | set(cidx.near(rr["chrom"], "-", rr["start0"], rr["end"])))
                     if i != c]
            fh.write(f"{r['name']}\t{c}:{labels.get(c, '')}\t{';'.join(other)}\n")
    # R2: CAT genes named like the family, reported beside
    with open(f"{d}/members_beside_R2.tsv", "w") as fh:
        fh.write("cat_gene_id\tlabel\tgene_name\tbiotype\tchrom\tstart0\tend\tstrand\tin_R1_member_image\n")
        n = 0
        for i in sorted(cat, key=lambda i: (cat[i]["chrom"], cat[i]["start0"])):
            nm = cat[i]["name"]
            if (nm.startswith("NPIP") or nm.startswith("TBC1D3")) and "-" not in nm:
                g = cat[i]
                fh.write(f"{i}\t{labels[i]}\t{nm}\t{g['biotype']}\t{g['chrom']}\t{g['start0']}\t{g['end']}\t{g['strand']}\t"
                         f"{'yes' if i in memb else 'no'}\n")
                n += 1
    print(f"[R2 beside] CAT genes named NPIP*/TBC1D3* without '-': {n}")
    # R3 read-throughs overlapping a member on the same strand (expr-recount --ignore)
    rt = readthrough_ids(cat, labels, rs, img)
    over = []
    for i in sorted(rt, key=lambda i: labels[i]):
        g = cat[i]
        hit = [m for m in memb if m != i and cat[m]["chrom"] == g["chrom"] and cat[m]["strand"] == g["strand"]
               and cat[m]["start0"] < g["end"] and g["start0"] < cat[m]["end"]]
        if hit:
            over.append((labels[i], rt[i], [labels[m] for m in hit]))
    with open(f"{d}/readthrough_over_members.txt", "w") as fh:
        for lab, why, hit in over:
            fh.write(f"{lab}\n")
    with open(f"{d}/readthrough_over_members.why.tsv", "w") as fh:
        fh.write("label\tR3_basis\toverlapped_members_same_strand\n")
        for lab, why, hit in over:
            fh.write(f"{lab}\t{why}\t{','.join(hit)}\n")
    print(f"[R3] read-throughs overlapping a member on the same strand: {[x[0] for x in over]}")
    # literature truth (31 RefSeq records) re-keyed through the step-1 tables (+ Amendment 2 where it applies)
    s1 = {}
    for r in tsv(STEP1_NPIP):
        s1[r["refseq_id"]] = r["cat_gene_id"]
    for r in tsv(STEP1_TBC):
        s1[r["refseq_id"]] = r["cat_gene_id"]
    lit = tsv(LIT_RS)
    out, ldrop = [], []
    for r in lit:
        rid = "gene-" + r["name"]
        c = s1.get(rid, "")
        if r["name"] in AMEND2_B_COPY:
            c = AMEND2_EXPECT[r["name"]]
        assert c == img[rid][0], (r["name"], c, img[rid])  # the step-1 tables and this module's R1 agree
        if not c:
            ldrop.append(r["name"])
            continue
        g = cat[c]
        out.append({**r, "name": labels[c], "chrom": g["chrom"], "start0": g["start0"], "end": g["end"],
                    "strand": g["strand"], "biotype": g["biotype"], "refseq_name": r["name"], "cat_gene_id": c})
    with open(f"{d}/lit_truth.cat.tsv", "w") as fh:
        cols = list(out[0].keys())
        fh.write("\t".join(cols) + "\n")
        for x in out:
            fh.write("\t".join(str(x[k]) for k in cols) + "\n")
    print(f"[literature] RefSeq records {len(lit)} -> CAT {len(out)}; dropped {ldrop}")
    # EXPR extra genes are decided after P/D/C exist (see cmd_expr_extra)


def cmd_expr_extra(a):
    """U genes (members ∪ member groups of P and D ∪ C leaves) absent from heavy/EXPR.counts.tsv -> expr_extra.txt."""
    root = a.root
    U = {r["name"] for r in tsv(f"{root}/light/members.tsv")}
    for fn, col in (("P.groups.tsv", "group_id"), ("D.groups.tsv", "group_id")):
        rows = tsv(f"{root}/light/{fn}")
        mg = {r[col] for r in rows if r["is_member"] == "yes" and "singleton" not in r[col]}
        U |= {r["name"] for r in rows if r[col] in mg}
    U |= {r["name"] for r in tsv(f"{root}/light/C.groups.tsv") if r["in_reference_trees"] == "yes"}
    expr = {r["gene_id"] for r in tsv(f"{root}/heavy/EXPR.counts.tsv")}
    extra = sorted(U - expr)
    with open(f"{root}/light/work/refseq/expr_extra.txt", "w") as fh:
        fh.writelines(f"{x}\n" for x in extra)
    print(f"[expr extra] U (pre-correction estimate) {len(U)}; absent from EXPR.counts.tsv {len(extra)}: {extra}")


# ------------------------------------------------------------------------------------------------ catalogs
CAT_CATALOGS = {"c15_17_22": f"{W}/o1_falsemerge/lit/annot_gencode", "c16_19_20": f"{W}/o1_falsemerge/lit/aj_ho/gencode"}


def cmd_catalogs(a):
    root = a.root
    cat = load_cat()
    labels = make_labels(cat)
    for tag, src in CAT_CATALOGS.items():
        out = f"{root}/light/work/catalogs/{tag}"
        os.makedirs(out, exist_ok=True)
        old_names = {r["idx"]: (r["name"], r["biotype"]) for r in tsv(f"{src}/nodes.tsv.names.tsv")}
        key = collections.defaultdict(list)
        for i, g in cat.items():
            key[(g["chrom"], g["start0"], g["end"], g["strand"], fmt_blocks(g["blocks"]), g["name"], g["biotype"])].append(i)
        for v in key.values():
            v.sort()
        used = collections.Counter()
        rows = tsv(f"{src}/nodes.tsv")
        with open(f"{out}/nodes.tsv.names.tsv", "w") as fh:
            fh.write("idx\tname\tbiotype\tcat_gene_id\n")
            for r in rows:
                nm, bt = old_names[r["idx"]]
                k = (r["chrom"], int(r["start"]), int(r["end"]), r["strand"], r["exons"], nm, bt)
                ids = key[k]
                assert ids, (tag, r["idx"], k[:5])
                i = ids[used[k] % len(ids)]
                used[k] += 1
                fh.write(f"{r['idx']}\t{labels[i]}\t{bt}\t{i}\n")
        subprocess.run(["cp", f"{src}/nodes.tsv", f"{out}/nodes.tsv"], check=True)
        print(f"[catalog {tag}] nodes {len(rows)} relabelled; coordinate+exon+name keys shared by > 1 CAT gene "
              f"{sum(1 for k, v in used.items() if v > 1)}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cmd", choices=("tables", "members", "catalogs", "expr-extra"))
    p.add_argument("--root", required=True, help="the CAT copy of layer_order/npip_tbc1d3")
    a = p.parse_args(argv)
    {"tables": cmd_tables, "members": cmd_members, "catalogs": cmd_catalogs, "expr-extra": cmd_expr_extra}[a.cmd](a)


if __name__ == "__main__":
    main()
