"""Map RefSeq CHM13 genes to the Soto et al. 2025 gene-ID convention (CAT CHM13_G* / Liftoff LOFF_G* ids).

Sources (inventory §1, §7):
  winloci_data/gencode_chm13/chm13v2.0_gencode.gff3   genome-wide CAT v2.0 annotation (source column 'CAT'); all 2,334
                                                       Soto-universe ids are present in it (checked, work/cat/)
  soto_replication/soto_gene_to_families.tsv          gene_id -> family ids (';'-separated) -> ambiguous (yes/no)
  bench/soto/soto_famCN_S1C.tsv                       Table S1C: Gene Name <-> Gene ID
Rule: exon unions (RefSeq: exon lines' gene=; CAT: exons of the gene's transcripts) on the same chrom and strand; a RefSeq
gene maps to the CAT gene with the largest shared exonic bp (ties: larger Jaccard). Reported twice: best CAT gene overall
(cat_*) and best CAT gene inside Soto's 2,334-gene universe (soto_*), with every Soto-universe gene overlapping it.
Name agreement is reported, never used to decide. soto_match_quality: strong = shared exonic bp >= 0.5 of BOTH exon
unions; partial = >= 0.5 of one; weak = neither.
"""
import bisect
import collections
import csv

SOTO = "/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication"
import os
# §6s2: repo root from THIS file, so the tool runs from any clone (register 887).
_RUSTLE_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
S1C = os.path.join(_RUSTLE_REPO, "bench", "soto", "soto_famCN_S1C.tsv")
LIGHT = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light"


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


def parse_blocks(s):
    return [tuple(map(int, x.split("-"))) for x in s.split(",")]


def load():
    cat = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{LIGHT}/work/cat/cat_genes_exons.tsv"), delimiter="\t"):
        cat[r["chrom"]].append((int(r["start0"]), int(r["end"]), r["strand"], r["gene_id"], r["name"],
                                parse_blocks(r["exons"])))
    idx = {}
    for c, v in cat.items():
        v.sort()
        idx[c] = (v, [x[0] for x in v], max(x[1] - x[0] for x in v))
    fams = {}
    for r in csv.DictReader(open(f"{SOTO}/soto_gene_to_families.tsv"), delimiter="\t"):
        fams[r["gene_id"]] = (r["family_ids_semicolon_sep"], r["ambiguous"])
    names = collections.defaultdict(set)
    name_to_ids = collections.defaultdict(set)
    for r in csv.DictReader(open(S1C), delimiter="\t"):
        names[r["Gene ID"]].add(r["Gene Name"])
        name_to_ids[r["Gene Name"]].add(r["Gene ID"])
    return idx, fams, names, name_to_ids


def map_gene(db, name, chrom, strand, exons):
    idx, fams, names, name_to_ids = db
    v, starts, ml = idx.get(chrom, ([], [], 0))
    s0, e0 = exons[0][0], exons[-1][1]
    lo, hi = bisect.bisect_left(starts, s0 - ml), bisect.bisect_left(starts, e0)
    L = sum(y - x for x, y in exons)
    hits = []
    for a0, a1, st, gid, cname, cex in v[lo:hi]:
        if st != strand or a1 <= s0:
            continue
        o = ov_bp(exons, cex)
        if o > 0:
            CL = sum(y - x for x, y in cex)
            hits.append((o, o / (L + CL - o), gid, cname, CL))
    hits.sort(key=lambda h: (h[0], h[1]), reverse=True)
    out = {k: "" for k in FIELDS}
    by_name = sorted(name_to_ids.get(name, set()))
    out["soto_ids_by_name"] = ";".join(by_name)
    out["soto_families_by_name"] = ";".join(sorted({f for i in by_name for f in fams.get(i, ("", ""))[0].split(";") if f}))
    if hits:
        o, jac, gid, cname, CL = hits[0]
        out.update(cat_gene_id=gid, cat_name=cname, cat_in_soto="yes" if gid in fams else "no", cat_ov_bp=o,
                   cat_jaccard=f"{jac:.3f}")
    soto_hits = [h for h in hits if h[2] in fams]
    if soto_hits:
        o, jac, gid, cname, CL = soto_hits[0]
        fam, amb = fams[gid]
        out.update(soto_gene_id=gid, soto_name=";".join(sorted(names.get(gid, set()))) or cname, soto_ov_bp=o,
                   soto_frac_of_refseq=f"{o / L:.3f}", soto_frac_of_soto=f"{o / CL:.3f}",
                   name_match="yes" if (name in names.get(gid, set()) or name == cname) else "no",
                   soto_families=fam, soto_ambiguous=amb,
                   soto_match_quality=("strong" if min(o / L, o / CL) >= 0.5 else "partial" if max(o / L, o / CL) >= 0.5
                                       else "weak"),
                   all_soto_overlapping=";".join(f"{h[2]}({h[3]}:{h[0]}bp:{fams[h[2]][0]})" for h in soto_hits))
    return out


FIELDS = ["cat_gene_id", "cat_name", "cat_in_soto", "cat_ov_bp", "cat_jaccard",
          "soto_gene_id", "soto_name", "soto_ov_bp", "soto_frac_of_refseq", "soto_frac_of_soto", "soto_match_quality", "name_match",
          "soto_families", "soto_ambiguous", "all_soto_overlapping", "soto_ids_by_name", "soto_families_by_name"]
