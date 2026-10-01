#!/usr/bin/env python3
"""Make the T2T-CHM13 v2.0 CAT/Liftoff annotation the default human annotation (2026-10-01).

Input: the consortium's v2.0 CAT/Liftoff GFF3 (GENCODE genes projected onto CHM13 by CAT, plus Liftoff; the annotation family Soto
et al. 2025 used: v4 on v1.0 for genes, the v2.0 CAT/Liftoff transcriptome for RNA). Outputs, next to `--out-prefix`:

  .genes.tsv           one row per gene: gene_id, gene_name, gene_biotype, source (CAT/Liftoff), chrom, start0, end, strand,
                       n_transcripts, exon-union blocks (0-based half-open, genomic order)
  .slim.gff3.gz(+.tbi) gene / transcript / exon lines only, sorted and tabix-indexed. Gene lines carry ID=Name=<gene_id> (unique;
                       CAT names repeat across paralogs) and gene_name=; exon lines carry Parent=<transcript> and gene=<gene_id>,
                       the attribute our GFF readers use to attach exons to genes.
  .vs_v4.tsv           per gene of CAT v4 (CHM13 v1.0, Soto's annotation): present in v2.0, same exon structure (exon union
                       relative to its first base, so the v1.0 -> v2.0 coordinate shift does not matter), and lengths
  .refseq_map.tsv      per RefSeq gene (v2.0): the CAT gene sharing the most exonic bases on the same strand (ties: larger
                       Jaccard), with shared bp and a quality class (strong = shared >= half of BOTH exon unions; partial = of one;
                       weak = neither; none). Names are reported, never used to match.

    python3 bench/annotation/cat_setup.py --cat chm13v2.0_gencode.gff3 --cat-v4-bed cat_v4.bed \
        --refseq chm13v2.0_RefSeq_full.gff.gz --out-prefix chm13v2.0_CAT_Liftoff
"""
import argparse
import bisect
import collections
import gzip
import subprocess


def attrs(col):
    return dict(x.split("=", 1) for x in col.rstrip().split(";") if "=" in x)


def merge(ivs):
    out = []
    for a, b in sorted(ivs):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


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


def blen(bl):
    return sum(e - s for s, e in bl)


def load_cat(path):
    genes, tx_gene, exons, tx_lines = {}, {}, collections.defaultdict(list), []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "transcript", "exon"):
                continue
            a = attrs(f[8])
            if f[2] == "gene":
                gid = a.get("gene_id", a.get("ID"))
                genes[gid] = dict(name=a.get("gene_name", a.get("Name", gid)), biotype=a.get("gene_biotype", ""), source=f[1],
                                  chrom=f[0], start0=int(f[3]) - 1, end=int(f[4]), strand=f[6], ntx=0)
            elif f[2] == "transcript":
                tid, par = a.get("ID"), a.get("Parent")
                tx_gene[tid] = par
                tx_lines.append((f[0], int(f[3]) - 1, int(f[4]), f[6], f[1], tid, par, a.get("transcript_biotype", "")))
            else:
                exons[a.get("Parent")].append((int(f[3]) - 1, int(f[4])))
    gene_ex = collections.defaultdict(list)
    for tid, ex in exons.items():
        g = tx_gene.get(tid)
        if g in genes:
            gene_ex[g].extend(ex)
    for t in tx_lines:
        if t[6] in genes:
            genes[t[6]]["ntx"] += 1
    for g in genes:
        genes[g]["blocks"] = merge(gene_ex.get(g, []))
    return genes, tx_lines, exons, tx_gene


def load_v4_bed(path):
    """CAT v4 BED (CHM13 v1.0, gene id in column 19): gene -> exon union."""
    ex = collections.defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 19:
            continue
        st = int(f[1])
        sizes = [int(x) for x in f[10].rstrip(",").split(",") if x]
        starts = [int(x) for x in f[11].rstrip(",").split(",") if x]
        ex[f[18]].extend((st + o, st + o + s) for s, o in zip(sizes, starts))
    return {g: merge(v) for g, v in ex.items()}


def load_refseq(path):
    genes, tx_gene, ex = {}, {}, collections.defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            i = a.get("ID", "")
            if i.startswith("gene-"):
                genes[i] = dict(name=a.get("Name", i[5:]), biotype=a.get("gene_biotype", ""), chrom=f[0], strand=f[6],
                                start0=int(f[3]) - 1, end=int(f[4]))
            elif f[2] == "exon":
                ex[a.get("Parent", "")].append((int(f[3]) - 1, int(f[4])))
            elif a.get("Parent", "").startswith("gene-") and i:
                tx_gene[i] = a["Parent"]
    gex = collections.defaultdict(list)
    for p, v in ex.items():
        g = p if p.startswith("gene-") else tx_gene.get(p)
        if g in genes:
            gex[g].extend(v)
    for g in genes:
        genes[g]["blocks"] = merge(gex.get(g, []))
    return genes


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("cat", "cat_v4_bed", "refseq", "out_prefix"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    P = a.out_prefix

    genes, tx_lines, exons, tx_gene = load_cat(a.cat)
    print(f"[cat v2.0] {len(genes):,} genes ({sum(g['source'] == 'Liftoff' for g in genes.values()):,} Liftoff), "
          f"{len(tx_lines):,} transcripts, {sum(len(v) for v in exons.values()):,} exon records")
    with open(P + ".genes.tsv", "w") as out:
        out.write("gene_id\tgene_name\tgene_biotype\tsource\tchrom\tstart0\tend\tstrand\tn_transcripts\texon_blocks\n")
        for gid, g in sorted(genes.items(), key=lambda kv: (kv[1]["chrom"], kv[1]["start0"], kv[0])):
            out.write(f"{gid}\t{g['name']}\t{g['biotype']}\t{g['source']}\t{g['chrom']}\t{g['start0']}\t{g['end']}\t{g['strand']}\t"
                      f"{g['ntx']}\t{','.join(f'{s}-{e}' for s, e in g['blocks'])}\n")

    rows = []
    for gid, g in genes.items():
        rows.append((g["chrom"], g["start0"] + 1, 0, f"{g['chrom']}\t{g['source']}\tgene\t{g['start0'] + 1}\t{g['end']}\t.\t{g['strand']}\t.\t"
                     f"ID={gid};Name={gid};gene_id={gid};gene_name={g['name']};gene_biotype={g['biotype']}"))
    for c, s0, e, st, src, tid, par, tb in tx_lines:
        if par not in genes:
            continue
        rows.append((c, s0 + 1, 1, f"{c}\t{src}\ttranscript\t{s0 + 1}\t{e}\t.\t{st}\t.\tID={tid};Parent={par};gene_id={par};"
                     f"transcript_biotype={tb}"))
        for k, (es, ee) in enumerate(sorted(exons.get(tid, []))):
            rows.append((c, es + 1, 2, f"{c}\t{src}\texon\t{es + 1}\t{ee}\t.\t{st}\t.\tID={tid}.exon{k + 1};Parent={tid};gene={par}"))
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    plain = P + ".slim.gff3"
    with open(plain, "w") as out:
        out.write("##gff-version 3\n")
        for r in rows:
            out.write(r[3] + "\n")
    subprocess.run(["bgzip", "-f", plain], check=True)
    subprocess.run(["tabix", "-f", "-p", "gff", plain + ".gz"], check=True)
    print(f"[slim] {len(rows):,} lines -> {plain}.gz (+ .tbi)")

    v4 = load_v4_bed(a.cat_v4_bed)
    same = diff = missing = 0
    with open(P + ".vs_v4.tsv", "w") as out:
        out.write("gene_id\tin_v2\tsame_structure\texonic_bp_v4\texonic_bp_v2\tn_blocks_v4\tn_blocks_v2\n")
        for gid, b4 in sorted(v4.items()):
            g = genes.get(gid)
            if g is None:
                missing += 1
                out.write(f"{gid}\t0\t0\t{blen(b4)}\t0\t{len(b4)}\t0\n")
                continue
            b2 = g["blocks"]
            r4 = [(s - b4[0][0], e - b4[0][0]) for s, e in b4]
            r2 = [(s - b2[0][0], e - b2[0][0]) for s, e in b2] if b2 else []
            ok = r4 == r2
            same += ok
            diff += not ok
            out.write(f"{gid}\t1\t{int(ok)}\t{blen(b4)}\t{blen(b2)}\t{len(b4)}\t{len(b2)}\n")
    print(f"[vs v4] {len(v4):,} v4 genes: same exon structure {same:,}, different {diff:,}, absent from v2.0 {missing:,}; "
          f"v2.0-only genes {len(set(genes) - set(v4)):,}")

    rs = load_refseq(a.refseq)
    idx = collections.defaultdict(list)
    for gid, g in genes.items():
        if g["blocks"]:
            idx[(g["chrom"], g["strand"])].append((g["start0"], g["end"], gid))
    for k in idx:
        idx[k].sort()
    starts = {k: [x[0] for x in v] for k, v in idx.items()}
    maxlen = {k: max(e - s for s, e, _ in v) for k, v in idx.items()}
    cls = collections.Counter()
    with open(P + ".refseq_map.tsv", "w") as out:
        out.write("refseq_gene\trefseq_name\trefseq_biotype\tchrom\tstart0\tend\tstrand\tcat_gene\tcat_name\tcat_biotype\t"
                  "shared_bp\tjaccard\tquality\n")
        for rid, r in sorted(rs.items(), key=lambda kv: (kv[1]["chrom"], kv[1]["start0"])):
            best = None
            if r["blocks"]:
                k = (r["chrom"], r["strand"])
                lst = idx.get(k, [])
                i = bisect.bisect_left(starts.get(k, []), r["start0"] - maxlen.get(k, 0))
                while i < len(lst) and lst[i][0] < r["end"]:
                    s, e, gid = lst[i]
                    if e > r["start0"]:
                        sh = ov_bp(r["blocks"], genes[gid]["blocks"])
                        if sh:
                            un = blen(r["blocks"]) + blen(genes[gid]["blocks"]) - sh
                            cand = (sh, sh / un, gid)
                            if best is None or cand[:2] > best[:2]:
                                best = cand
                    i += 1
            if best:
                sh, jac, gid = best
                q = ("strong" if sh >= 0.5 * blen(r["blocks"]) and sh >= 0.5 * blen(genes[gid]["blocks"]) else
                     "partial" if sh >= 0.5 * blen(r["blocks"]) or sh >= 0.5 * blen(genes[gid]["blocks"]) else "weak")
                g = genes[gid]
                out.write(f"{rid}\t{r['name']}\t{r['biotype']}\t{r['chrom']}\t{r['start0']}\t{r['end']}\t{r['strand']}\t{gid}\t"
                          f"{g['name']}\t{g['biotype']}\t{sh}\t{jac:.3f}\t{q}\n")
            else:
                q = "none"
                out.write(f"{rid}\t{r['name']}\t{r['biotype']}\t{r['chrom']}\t{r['start0']}\t{r['end']}\t{r['strand']}\t\t\t\t0\t0\tnone\n")
            cls[q] += 1
    print(f"[refseq map] {len(rs):,} RefSeq genes -> CAT: " + ", ".join(f"{k} {v:,}" for k, v in cls.most_common()))


if __name__ == "__main__":
    main()
