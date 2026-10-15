#!/usr/bin/env python3
"""What are the `copy`-labelled unexplained clusters of families whose copies are all in the genome? (descriptive, post hoc;
docs/PREREG_o3_realdata_2026-10-15.md, follow-up.)

    copyclusters.py origins LABELS PREFIX BAM COPIES.tsv OUT.tsv

origins  per `copy` cluster (status unexplained, est_div > delta) of `copy_assign --psv-clusters` run PREFIX: its reads' copy of
         origin (labels.tsv: the copy their baseline primary lay on) against the host it was clustered on; their primary alignment
         in BAM (on the origin copy's locus or elsewhere, identity = 1 - NM / aligned bases); and the share of their aligned bases
         outside the origin copy's catalog exons (exons the catalog copy lacks).
sites    the prototype's own clustering (bench/o3_unexplained/psvclust.py: align, cluster_host, linkage) re-run on each host that
         holds a `copy` cluster (the host's reads from PREFIX.assignments.tsv, every family copy, the host locus +- PAD); each copy
         cluster is matched to the prototype cluster sharing most of its reads, and its unexplained alt sites are split into
         UNCOVERED (no catalog copy has sequence there: unexplained by construction) and COVERED (some copy covers the position
         and carries the reference allele), and into inside / outside the host copy's own exons.

    copyclusters.py sites PREFIX COPIES.tsv COPIES.fa READS.fa GENOME OUT.tsv
haps     the clusters' reads aligned to KB3781's own haplotype assemblies (mat, pat; `reads` writes them, the PAFs come from
         `minimap2 -c -x splice:hq -uf --secondary=no`): per read its best hit (identity = matches / block, coverage of the read)
         and whether it lands on a lifted copy of the cluster's family, of another family, or elsewhere (the _pri chromosome's own
         haplotype has identical coordinates; the other one is reached through the frozen _pri -> B alignments, `out/chr*.paf`).

    copyclusters.py reads OUT.fa PREFIX READS.fa [PREFIX READS.fa ...]
    copyclusters.py haps PREFIX COPIES.tsv ALL_COPIES.tsv MAT.paf PAT.paf CHRMAP OUT_DIR OUT.tsv
classify each cluster from the three tables, rules applied in order (reads with a haplotype hit >= 0.995 identity and
         >= 0.9 coverage are "matched"):
           other_family_copy   most matched reads land on another family's lifted copy (reads from other loci)
           not_in_genome       most reads match neither haplotype at >= 0.995 (outside KB3781's genome: other species, artifacts)
           absent_from_ref     most matched reads land on a haplotype locus that is no lifted copy of the family (a copy the
                               reference collapses or lacks)
           catalog_lacks_exons own-family reads; >= half the unexplained sites lie where no catalog copy has sequence
           minority_segment    own-family reads; the unexplained sites are covered by a median < 0.2 of the cluster's reads
           own_copy_variant    own-family reads otherwise (allele-like differences carried by most of the cluster)

    copyclusters.py classify ORIGINS.tsv SITES.tsv HAPS.tsv OUT.tsv
"""
import collections
import csv
import json
import sys

import pysam

DELTA = 0.00958
M_OPS = (0, 7, 8)


def copy_clusters(prefix):
    with open(prefix + ".psv_clusters.tsv") as fh:
        cl = [c for c in csv.DictReader(fh, delimiter="\t") if c["status"] == "unexplained" and float(c["est_div"]) > DELTA]
    members = collections.defaultdict(list)
    with open(prefix + ".psv_cluster_reads.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            members[(r["family"], r["host"], r["cluster"])].append(r["read_name"])
    return cl, members


def aligned(rec):
    pos, out = rec.reference_start, []
    for op, n in rec.cigartuples:
        if op in M_OPS:
            out.append((pos, pos + n))
            pos += n
        elif op in (2, 3):
            pos += n
    return out


def inter(a, b):
    i = j = n = 0
    a, b = sorted(a), sorted(b)
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        n += max(0, hi - lo)
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return n


def cmd_sites(argv):
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "o3_unexplained"))
    import psvclust as PC
    prefix, copies_p, copies_fa, reads_fa, genome_p, out_p = argv
    rows_c = list(csv.DictReader(open(copies_p), delimiter="\t"))
    cfa = PC.read_fa(copies_fa)
    seq_of = {(n.split("|")[0], int(n.split("|")[1])): s for n, s in cfa.items()}
    fam_copies = collections.defaultdict(list)
    for r in rows_c:
        ex = sorted(tuple(int(x) for x in b.split("-")) for b in r["exons"].split(","))
        fam_copies[r["family_id"]].append((int(r["copy_idx"]), r["chrom"], int(r["start"]), int(r["end"]), ex,
                                           seq_of[(r["family_id"], int(r["copy_idx"]))]))
    by_host = collections.defaultdict(list)
    for r in csv.DictReader(open(prefix + ".assignments.tsv"), delimiter="\t"):
        if r["catalog_copy_idx"] not in ("", "NA", "-"):
            by_host[(r["family_id"], int(r["catalog_copy_idx"]))].append(r["read_name"])
    cl, members = copy_clusters(prefix)
    reads = PC.read_fa(reads_fa)
    genome = pysam.FastaFile(genome_p)
    tmp = out_p + ".tmp"
    out = []
    for (fam, host), group in sorted(collections.defaultdict(list, {k: [c for c in cl if (c["family"], c["host"]) == k]
                                                                     for k in {(c["family"], c["host"]) for c in cl}}).items()):
        idxs = [int(x) for x in host.split(",")]
        cps = fam_copies[fam]
        rep = min(c for c in cps if c[0] in idxs)
        _, chrom, s0, e0, hex_, _ = rep
        lo, hi = max(0, s0 - PC.PAD), min(genome.get_reference_length(chrom), e0 + PC.PAD)
        ref = genome.fetch(chrom, lo, hi).upper()
        names = sorted({n for i in idxs for n in by_host.get((fam, i), [])})
        cov = PC.align("host", ref, {n: reads[n] for n in names if n in reads}, tmp)
        ccov = PC.align("host", ref, {f"copy{c[0]}": c[5] for c in cps}, tmp)
        sites, clusters, vec = PC.cluster_host(ref, cov)
        for c in group:
            mem = set(members[(c["family"], c["host"], c["cluster"])])
            best = max(clusters, key=lambda x: len(mem & set(x[0])), default=None)
            if best is None:
                continue
            pm, cons = best
            alts = [p for p, al in cons.items() if al == "alt"]
            unexpl = [p for p in alts if not any(ccov[cn].get(p) == sites[p] for cn in ccov)]
            uncovered = [p for p in unexpl if not any(p in ccov[cn] for cn in ccov)]
            in_host = [p for p in unexpl if any(a <= p + lo < b for a, b in hex_)]
            # how many of the cluster's reads carry its unexplained signal: per unexplained site, the share of the cluster's
            # reads that cover it; and how tightly the sites bunch (distinct 200-bp windows)
            mreads = [m for m in pm if m in cov]
            share = sorted(sum(1 for m in mreads if p in cov[m]) / max(1, len(mreads)) for p in unexpl)
            out.append(dict(family=fam, host=host, cluster=c["cluster"], rust_reads=len(mem),
                            match_jaccard=round(len(mem & set(pm)) / len(mem | set(pm)), 3), rust_est_div=c["est_div"],
                            n_alt=len(alts), n_unexpl=len(unexpl), unexpl_uncovered=len(uncovered),
                            unexpl_covered=len(unexpl) - len(uncovered), unexpl_in_host_exons=len(in_host),
                            unexpl_outside_host_exons=len(unexpl) - len(in_host),
                            median_read_share=round(share[len(share) // 2], 3) if share else "NA",
                            windows_200bp=len({p // 200 for p in unexpl})))
    for x in (".t.fa", ".q.fa"):
        if os.path.exists(tmp + x):
            os.remove(tmp + x)
    with open(out_p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"{len(out)} copy clusters -> {out_p}")
    return 0


def cmd_reads(argv):
    out_p, pairs = argv[0], argv[1:]
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "o3_unexplained"))
    import psvclust as PC
    seen = set()
    with open(out_p, "w") as fh:
        for prefix, reads_fa in zip(pairs[::2], pairs[1::2]):
            cl, members = copy_clusters(prefix)
            want = {r for c in cl for r in members[(c["family"], c["host"], c["cluster"])]} - seen
            reads = PC.read_fa(reads_fa)
            for r in sorted(want):
                fh.write(f">{r}\n{reads[r]}\n")
            seen |= want
    print(f"{len(seen)} reads -> {out_p}")
    return 0


def best_hits(paf_p):
    best = {}
    with open(paf_p) as fh:
        for ln in fh:
            f = ln.split("\t")
            if "tp:A:P" not in ln:
                continue
            ident, cov = int(f[9]) / max(1, int(f[10])), (int(f[3]) - int(f[2])) / max(1, int(f[1]))
            if f[0] not in best or ident * cov > best[f[0]][0] * best[f[0]][1]:
                best[f[0]] = (ident, cov, f[5], int(f[7]), int(f[8]))
    return best


def lifter(chrmap_p, out_dir):
    """_pri (chrom, pos) -> {haplotype contig: position}: identical coordinates on the chromosome's own haplotype, and the
    frozen asm5 alignment (`out/chr*.paf`, primary records, 10 Mb query chunks) for the other one."""
    import glob
    rows = list(csv.DictReader(open(chrmap_p), delimiter="\t"))
    same = {r["pri"]: r["same_name"] for r in rows}
    recs = collections.defaultdict(list)
    for path in glob.glob(f"{out_dir}/chr*.paf"):
        with open(path) as fh:
            for ln in fh:
                f = ln.split("\t")
                if "tp:A:P" not in ln:
                    continue
                c, off = f[0].rsplit(":", 1)
                recs[c].append((int(off) + int(f[2]), int(off) + int(f[3]), f[4], f[5], int(f[7]), int(f[8])))

    def lift(chrom, pos):
        out = {}
        if chrom in same:
            out[same[chrom]] = pos
        for qs, qe, st, t, ts, te in recs.get(chrom, ()):
            if qs <= pos < qe:
                out[t] = ts + (pos - qs) if st == "+" else te - (pos - qs)
                break
        return out
    return lift


def cmd_haps(argv):
    prefix, copies_p, all_copies_p, mat_p, pat_p, chrmap_p, out_dir, out_p = argv
    lift = lifter(chrmap_p, out_dir)
    lifted = collections.defaultdict(list)          # haplotype contig -> [(lo, hi, family)]
    for r in csv.DictReader(open(all_copies_p), delimiter="\t"):
        a, b = lift(r["chrom"], int(r["start"])), lift(r["chrom"], int(r["end"]) - 1)
        for t in set(a) & set(b):
            lo, hi = sorted((a[t], b[t]))
            lifted[t].append((lo - 5000, hi + 5000, r["family_id"]))
    hits = {}
    for name, path in (("mat", mat_p), ("pat", pat_p)):
        for r, h in best_hits(path).items():
            if r not in hits or h[0] * h[1] > hits[r][1][0] * hits[r][1][1]:
                hits[r] = (name, h)
    cl, members = copy_clusters(prefix)
    rows = []
    for c in cl:
        reads = members[(c["family"], c["host"], c["cluster"])]
        where = collections.Counter()
        ids = []
        for r in reads:
            if r not in hits:
                where["no_hit"] += 1
                continue
            hap, (ident, cov, t, ts, te) = hits[r]
            ids.append(ident)
            if ident < 0.995 or cov < 0.9:
                where["no_match_0.995"] += 1
                continue
            fams = {f for lo, hi, f in lifted.get(t, ()) if lo < te and ts < hi}
            where["own_family" if c["family"] in fams else "other_family" if fams else "elsewhere"] += 1
        ids.sort()
        rows.append(dict(family=c["family"], host=c["host"], cluster=c["cluster"], n_reads=len(reads),
                         hap_median_identity=round(ids[len(ids) // 2], 4) if ids else "NA",
                         **{k: where.get(k, 0) for k in ("own_family", "other_family", "elsewhere", "no_match_0.995", "no_hit")}))
    with open(out_p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} copy clusters -> {out_p}")
    return 0


def classify(site, hap):
    n = int(hap["n_reads"])
    matched = int(hap["own_family"]) + int(hap["other_family"]) + int(hap["elsewhere"])
    if matched and int(hap["other_family"]) > max(int(hap["own_family"]), int(hap["elsewhere"])):
        return "other_family_copy"
    if int(hap["no_match_0.995"]) + int(hap["no_hit"]) > n / 2:
        return "not_in_genome"
    if int(hap["elsewhere"]) > int(hap["own_family"]):
        return "absent_from_ref"
    if int(site["unexpl_uncovered"]) >= int(site["n_unexpl"]) / 2:
        return "catalog_lacks_exons"
    if site["median_read_share"] != "NA" and float(site["median_read_share"]) < 0.2:
        return "minority_segment"
    return "own_copy_variant"


def cmd_classify(argv):
    origins_p, sites_p, haps_p, out_p = argv
    key = lambda r: (r["family"], r["host"], r["cluster"])  # noqa: E731
    o = {key(r): r for r in csv.DictReader(open(origins_p), delimiter="\t")}
    st = {key(r): r for r in csv.DictReader(open(sites_p), delimiter="\t")}
    rows = []
    for h in csv.DictReader(open(haps_p), delimiter="\t"):
        k = key(h)
        rows.append(dict(family=k[0], host=k[1], cluster=k[2], n_reads=h["n_reads"], est_div=o[k]["est_div"],
                         ref_identity=o[k]["median_identity"], hap_identity=h["hap_median_identity"],
                         category=classify(st[k], h)))
    with open(out_p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    c = collections.Counter(r["category"] for r in rows)
    r_ = collections.Counter()
    for r in rows:
        r_[r["category"]] += int(r["n_reads"])
    print(json.dumps({k: [c[k], r_[k], len({x["family"] for x in rows if x["category"] == k})] for k in sorted(c)}))
    return 0


def main(argv):
    if argv[0] == "classify":
        return cmd_classify(argv[1:])
    if argv[0] == "sites":
        return cmd_sites(argv[1:])
    if argv[0] == "reads":
        return cmd_reads(argv[1:])
    if argv[0] == "haps":
        return cmd_haps(argv[1:])
    cmd, labels_p, prefix, bam_p, copies_p, out_p = argv
    label = {}
    with open(labels_p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            label[r["read"]] = r["copy"]
    copies, idx_of = {}, {}
    with open(copies_p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ex = sorted(tuple(int(x) for x in b.split("-")) for b in r["exons"].split(","))
            key = f"{r['family_id']}:{r['copy_idx']}"
            copies[key] = (r["chrom"], int(r["start"]), int(r["end"]), ex)
            idx_of[r["tid"]] = key          # A13 labels name a copy by its panel tid (family:panel index), not copy_idx
    label = {r: idx_of.get(k, k) for r, k in label.items()}
    cl, members = copy_clusters(prefix)
    want = {r for c in cl for r in members[(c["family"], c["host"], c["cluster"])]}
    prim = {}
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.query_name in want and not rec.flag & 0x904 and rec.cigartuples:
                prim[rec.query_name] = (rec.reference_name, rec.reference_start, rec.reference_end, aligned(rec),
                                        rec.get_tag("NM") if rec.has_tag("NM") else 0)
    rows = []
    for c in cl:
        reads = members[(c["family"], c["host"], c["cluster"])]
        host_copy = f"{c['family']}:{c['host'].split(',')[0]}"
        origins = collections.Counter(label.get(r, "?") for r in reads)
        on_host = sum(n for k, n in origins.items() if k.split(":")[-1] in c["host"].split(",") and k.startswith(c["family"] + ":"))
        at_origin = outside = total = 0
        idents = []
        for r in reads:
            p, o = prim.get(r), copies.get(label.get(r, ""))
            if not p or not o:
                continue
            bl = p[3]
            al = sum(e - s for s, e in bl)
            idents.append(1 - p[4] / max(1, al))
            if p[0] == o[0] and p[1] < o[2] and o[1] < p[2]:
                at_origin += 1
                total += al
                outside += al - inter(bl, o[3])
        idents.sort()
        top = origins.most_common(1)[0]
        rows.append(dict(family=c["family"], host=c["host"], cluster=c["cluster"], n_reads=c["n_reads"], n_sites=c["n_sites"],
                         n_unexplained_alt=c["n_unexplained_alt"], est_div=c["est_div"], linkage=c["linkage"],
                         reads_from_host=round(on_host / len(reads), 3), top_origin=f"{top[0]}={top[1]}",
                         n_origins=len(origins), primary_at_origin=round(at_origin / len(reads), 3),
                         median_identity=round(idents[len(idents) // 2], 4) if idents else "NA",
                         outside_catalog_exons=round(outside / total, 3) if total else "NA"))
    with open(out_p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} copy clusters -> {out_p}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
