#!/usr/bin/env python3
"""Would a reference with every isoform explain O3's unexplained clusters? A proxy for a family variation graph (post hoc,
descriptive; docs/PREREG_o3_family_best_2026-10-15.md, Results).

    genomic_explain.py PREFIX COPIES.tsv COPIES.fa READS.fa BAM GENOME ALL_COPIES.tsv OUT.tsv [KEPT.tsv]

For every `unexplained` cluster of `copy_assign --psv-clusters` run PREFIX, the host is re-clustered with the prototype
(bench/o3_unexplained/psvclust.py, as copyclusters.py sites; KEPT.tsv from family_best.py restricts it to the kept pairs) and
the Rust cluster is matched to the prototype cluster sharing most reads. Register 1294's explained rule (`explains`: one locus
carries every consensus alt, covers every alt position, and has no other base where the consensus is the reference) is then
asked against three reference sets:
  G0  each family copy's spliced representative (COPIES.fa), as the method does; a check, none should explain
  G1  each family copy's GENOMIC locus (locus_start - PAD .. locus_end + PAD): every isoform, intron and PSV of the family
  G2  G1 plus every genomic locus the cluster's reads align to in BAM (primary or secondary records, merged within 5 kb, +- PAD):
      anything in the genome the aligner offers for these reads
Genomic loci are aligned to the host locus with `minimap2 -c --cs=long -x asm20` (the union of their primary records).
OUT.tsv per cluster: the explanation found at each level (`-` = none); d_host and, at G1 / G2, the NEAREST locus by
consensus divergence (g2_families = the nearest locus's catalog families, which need not be the explaining locus).
"""
import collections
import csv
import os
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "o3_unexplained"))
import psvclust as PC  # noqa: E402

MERGE_GAP, MAX_LOCI, MIN_COVER = 5_000, 30, 0.5
BASES = set("ACGT")


def explains(cons, sites, ref, cc):
    """Register 1294's rule (`host_rows`): the locus carries every alt (where it has a base), has a base at every alt, and no
    other base where the cluster's consensus is the reference."""
    alts = [p for p, a in cons.items() if a == "alt"]
    carries = all(cc.get(p, sites[p]) == sites[p] for p in alts)
    no_ref_conflict = not any(a == "ref" and p in cc and cc[p] != ref[p] for p, a in cons.items())
    return carries and no_ref_conflict and all(p in cc for p in alts)


def consensus(covs, min_reads=PC.MIN_SUPPORT):
    """The cluster's majority base at every host position that >= min_reads of its reads cover with a base."""
    cnt = collections.defaultdict(collections.Counter)
    for cov in covs:
        for p, b in cov.items():
            if b in BASES:
                cnt[p][b] += 1
    return {p: c.most_common(1)[0][0] for p, c in cnt.items() if sum(c.values()) >= min_reads}


def divergence(cons, lc):
    """(mismatches / compared, compared / consensus positions) between a cluster's consensus and a locus aligned to the host,
    over the positions where both have a base."""
    both = [p for p in cons if lc.get(p) in BASES]
    if not both:
        return 1.0, 0.0
    return sum(cons[p] != lc[p] for p in both) / len(both), len(both) / len(cons)


def nearest(cons, covs, skip=()):
    """The locus closest to the consensus among those covering >= MIN_COVER of it: (name, divergence, cover) or None."""
    best = None
    for n, lc in covs.items():
        if n in skip:
            continue
        d, cover = divergence(cons, lc)
        if cover >= MIN_COVER and (best is None or d < best[1]):
            best = (n, d, cover)
    return best


def merge_spans(spans, gap=MERGE_GAP):
    out = []
    for c, s, e in sorted(spans):
        if out and out[-1][0] == c and s <= out[-1][2] + gap:
            out[-1] = (c, out[-1][1], max(out[-1][2], e))
        else:
            out.append((c, s, e))
    return out


def align_asm(ref, queries, tmp):
    """Genomic loci against the host locus: per query, the union of its primary records' bases (best-matching record first)."""
    with open(f"{tmp}.t.fa", "w") as f:
        f.write(f">host\n{ref}\n")
    with open(f"{tmp}.q.fa", "w") as f:
        for n, s in queries.items():
            f.write(f">{n}\n{s}\n")
    out = subprocess.run(["minimap2", "-c", "--cs=long", "-x", "asm20", f"{tmp}.t.fa", f"{tmp}.q.fa"], check=True,
                         capture_output=True, text=True).stdout
    recs = collections.defaultdict(list)
    for ln in out.splitlines():
        f = ln.split("\t")
        tags = dict((x[:5], x[5:]) for x in f[12:])
        if tags.get("tp:A:") == "P":
            recs[f[0]].append((int(f[9]), int(f[7]), tags.get("cs:Z:", "")))
    cov = {}
    for q, rs in recs.items():
        bases = {}
        for _, ts, cs in sorted(rs, reverse=True):
            for p, b in PC.cs_bases(ts, cs).items():
                bases.setdefault(p, b)
        cov[q] = bases
    return cov


def main(argv):
    prefix, copies_p, copies_fa, reads_fa, bam_p, genome_p, all_p, out_p, *kept_p = argv
    kept = None
    if kept_p:
        kept = {(r["family"], r["read"]) for r in csv.DictReader(open(kept_p[0]), delimiter="\t") if r["kept"] == "1"}
    seq_of = {(n.split("|")[0], int(n.split("|")[1])): s for n, s in PC.read_fa(copies_fa).items()}
    fam = collections.defaultdict(list)
    for r in csv.DictReader(open(copies_p), delimiter="\t"):
        k = (r["family_id"], int(r["copy_idx"]))
        fam[r["family_id"]].append(dict(idx=k[1], chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]),
                                        lstart=int(r.get("locus_start") or r["start"]), lend=int(r.get("locus_end") or r["end"]),
                                        seq=seq_of[k]))
    cat = collections.defaultdict(list)
    for r in csv.DictReader(open(all_p), delimiter="\t"):
        cat[r["chrom"]].append((int(r["start"]), int(r["end"]), r["family_id"]))
    by_host = collections.defaultdict(list)
    for r in csv.DictReader(open(prefix + ".assignments.tsv"), delimiter="\t"):
        if r["catalog_copy_idx"].isdigit() and (kept is None or (r["family_id"], r["read_name"]) in kept):
            by_host[(r["family_id"], int(r["catalog_copy_idx"]))].append(r["read_name"])
    with open(prefix + ".psv_clusters.tsv") as fh:
        ucl = [c for c in csv.DictReader(fh, delimiter="\t") if c["status"] == "unexplained"]
    members = collections.defaultdict(set)
    for r in csv.DictReader(open(prefix + ".psv_cluster_reads.tsv"), delimiter="\t"):
        members[(r["family"], r["host"], r["cluster"])].add(r["read_name"])
    want = {m for c in ucl for m in members[(c["family"], c["host"], c["cluster"])]}
    spans = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.query_name in want and not rec.is_unmapped and not rec.is_supplementary:
                spans[rec.query_name].append((rec.reference_name, rec.reference_start, rec.reference_end))
    reads, genome = PC.read_fa(reads_fa), pysam.FastaFile(genome_p)
    tmp = out_p + ".tmp"
    rows = []
    for (f, host) in sorted({(c["family"], c["host"]) for c in ucl}):
        idxs = [int(x) for x in host.split(",")]
        rep = min((c for c in fam[f] if c["idx"] in idxs), key=lambda c: c["idx"])
        lo, hi = max(0, rep["start"] - PC.PAD), min(genome.get_reference_length(rep["chrom"]), rep["end"] + PC.PAD)
        ref = genome.fetch(rep["chrom"], lo, hi).upper()
        names = sorted({n for i in idxs for n in by_host.get((f, i), [])})
        cov = PC.align("host", ref, {n: reads[n] for n in names if n in reads}, tmp)
        g0 = PC.align("host", ref, {f"copy{c['idx']}": c["seq"] for c in fam[f]}, tmp)
        gen = lambda ch, s, e: genome.fetch(ch, max(0, s - PC.PAD), min(genome.get_reference_length(ch), e + PC.PAD)).upper()  # noqa: E731
        g1 = align_asm(ref, {f"copy{c['idx']}": gen(c["chrom"], c["lstart"], c["lend"]) for c in fam[f]}, tmp)
        sites, clusters, _ = PC.cluster_host(ref, cov)
        for c in (x for x in ucl if (x["family"], x["host"]) == (f, host)):
            mem = members[(f, host, c["cluster"])]
            best = max(clusters, key=lambda x: len(mem & set(x[0])), default=None)
            if best is None:
                continue
            pm, cons = best
            hit = lambda covs: next((n for n, cc in covs.items() if explains(cons, sites, ref, cc)), "-")  # noqa: E731
            loci = merge_spans([s for m in mem for s in spans.get(m, [])])
            loci = [x for x in loci if not (x[0] == rep["chrom"] and x[1] < hi and lo < x[2])]
            support = lambda x: sum(1 for m in mem if any(s[0] == x[0] and s[1] < x[2] and x[1] < s[2] for s in spans.get(m, [])))  # noqa: E731
            loci = sorted(loci, key=lambda x: (-support(x), x))[:MAX_LOCI]
            e0, e1 = hit(g0), hit(g1)
            lc = align_asm(ref, {f"{ch}:{s}-{e}": gen(ch, s, e) for ch, s, e in loci}, tmp) if loci else {}
            g2n = e1 if e1 != "-" else hit(lc)
            cons_all = consensus([cov[m] for m in pm if m in cov])
            d_host, _ = divergence(cons_all, {p: ref[p] for p in cons_all})
            n1 = nearest(cons_all, g1, skip={f"copy{i}" for i in idxs})
            n2 = nearest(cons_all, {**{k: v for k, v in g1.items() if k not in {f"copy{i}" for i in idxs}}, **lc})

            def fams_of(name):
                if name is None or name.startswith("copy"):
                    return f if name else "-"
                ch, se = name.rsplit(":", 1)
                s, e = (int(x) for x in se.split("-"))
                fs = sorted({g for a, b, g in cat.get(ch, ()) if a < e and s < b})
                return ",".join(fs) if fs else "outside_catalog"
            rows.append(dict(family=f, host=host, cluster=c["cluster"], n_reads=c["n_reads"], est_div=c["est_div"],
                             match_jaccard=round(len(mem & set(pm)) / len(mem | set(pm)), 3),
                             n_alt=sum(1 for a in cons.values() if a == "alt"), n_loci=len(loci),
                             explained_g0=e0, explained_g1=e1, explained_g2=g2n, d_host=round(d_host, 5),
                             g1_nearest=n1[0] if n1 else "-", g1_d=round(n1[1], 5) if n1 else "NA",
                             g2_nearest=n2[0] if n2 else "-", g2_d=round(n2[1], 5) if n2 else "NA",
                             g2_cover=round(n2[2], 3) if n2 else "NA", g2_families=fams_of(n2[0] if n2 else None)))
    for x in (".t.fa", ".q.fa"):
        if os.path.exists(tmp + x):
            os.remove(tmp + x)
    with open(out_p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n = len(rows)
    print(f"{n} unexplained clusters: explained by G0 {sum(r['explained_g0'] != '-' for r in rows)}, "
          f"G1 {sum(r['explained_g1'] != '-' for r in rows)}, G2 {sum(r['explained_g2'] != '-' for r in rows)} -> {out_p}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
