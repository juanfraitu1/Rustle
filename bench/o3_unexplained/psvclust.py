#!/usr/bin/env python3
"""DEVELOPMENT prototype: cluster the reads O2 places on each copy by their variant sites, and call the clusters no reference copy
explains (docs/superpowers/specs/2026-10-15-o3-unexplained-multimappers-design.md). Avoids §6ai's traps: sites come from the READS
(the absent copy's own variants included), distances use co-covered sites only, and the number of clusters comes from a support test,
not from a k+1 objective.

    psvclust.py WORLD VARIANT ARM [--score]     # reads OUT/<variant>.<arm>/A.assignments.tsv; writes .../psvclust.tsv

Per host class Y (the copies of one family with an identical spliced sequence, e.g. twins, are one class; the class is the reads'
O2 best copy, `catalog_copy_idx`):
 1. Align the class's reads AND every copy of the family (spliced sequences) to Y's padded locus (`minimap2 -c --cs=long -x splice:hq
    -uf`); per read the base at each locus position it covers.
 2. Sites: positions where >= MIN_SUPPORT reads carry the same non-reference base.
 3. Clusters: reads in order of decreasing site coverage join the compatible cluster with the most co-covered sites (conflicts at
    co-covered sites <= max(0, (co-covered - 9) // 10)), else open a new one; clusters < MIN_SUPPORT reads are dissolved into their
    best compatible cluster or kept as noise; two clusters merge when no site separates them with >= MIN_SUPPORT reads on each side.
 4. Each cluster's consensus allele at every site; its alt sites; a copy explains the cluster when it carries every alt allele the
    cluster's reads cover and no site where the cluster is reference-allele but the copy differs (co-covered by the copy).
    UNEXPLAINED = >= MIN_ALT alt sites that no copy explains AND those sites are LINKED: over every pair of them, among the host's reads
    covering both and carrying the alt at either, the fraction carrying it at both; the median over pairs (pairs with >= MIN_SUPPORT such
    reads) must be >= MIN_LINK. A missing sequence is a haplotype (its alts travel together); RNA editing and shared sequencing errors
    hit sites independently. Unlinked clusters are `unlinked` (with their A->G / T->C fraction reported: the editing signature); a
    cluster with no evaluable pair is `unresolved`. Divergence estimate = unexplained alt sites / bases the cluster covers.
"""
import argparse
import collections
import csv
import json
import os
import re
import subprocess
import sys

MIN_SUPPORT = 3
MIN_ALT = 2
MIN_LINK = 0.8
PAD = 2000


def read_fa(path):
    out, name = {}, None
    for ln in open(path):
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            out[name] = []
        else:
            out[name].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


def read_fq(path):
    out = {}
    lines = open(path).read().splitlines()
    for i in range(0, len(lines), 4):
        out[lines[i][1:]] = lines[i + 1]
    return out


def cs_bases(ts, cs):
    """`cs:Z:` long form -> {target position: read base} for the positions the alignment covers (deletions as '-'; introns skipped)"""
    out, t = {}, ts
    for op, body in re.findall(r"([=*+\-~:])([A-Za-z0-9]+)", cs):
        if op == "=":
            for i, b in enumerate(body.upper()):
                out[t + i] = b
            t += len(body)
        elif op == "*":
            out[t] = body[1].upper()
            t += 1
        elif op == "-":
            for i in range(len(body)):
                out[t + i] = "-"
            t += len(body)
        elif op == "~":
            m = re.match(r"[a-z]{2}(\d+)[a-z]{2}", body)
            t += int(m.group(1)) if m else 0
        elif op == ":":
            t += int(body)
    return out


def align(target_name, target_seq, queries, tmp):
    with open(f"{tmp}.t.fa", "w") as f:
        f.write(f">{target_name}\n{target_seq}\n")
    with open(f"{tmp}.q.fa", "w") as f:
        for n, s in queries.items():
            f.write(f">{n}\n{s}\n")
    out = subprocess.run(f"minimap2 -c --cs=long -x splice:hq -uf -N 0 {tmp}.t.fa {tmp}.q.fa 2>/dev/null", shell=True, check=True,
                         capture_output=True, text=True, executable="/bin/bash").stdout
    best = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        tags = dict((x[:5], x[5:]) for x in f[12:])
        if tags.get("tp:A:") != "P":
            continue
        m = int(f[9])
        if f[0] not in best or m > best[f[0]][0]:
            best[f[0]] = (m, cs_bases(int(f[7]), tags.get("cs:Z:", "")))
    return {k: v[1] for k, v in best.items()}


def linkage(sites_u, vec):
    """median over pairs of unexplained alt sites of P(alt at both | alt at either), among reads covering both; None if no pair has
    >= MIN_SUPPORT reads with an alt at either"""
    import statistics
    vals = []
    su = sorted(sites_u)
    for i in range(len(su)):
        for j in range(i + 1, len(su)):
            a, b = su[i], su[j]
            either = both = 0
            for v in vec.values():
                if a in v and b in v:
                    xa, xb = v[a] == "alt", v[b] == "alt"
                    if xa or xb:
                        either += 1
                        both += xa and xb
            if either >= MIN_SUPPORT:
                vals.append(both / either)
    return statistics.median(vals) if vals else None


def cluster_host(ref, reads_cov):
    """ref: the host locus; reads_cov: {read: {pos: base}} -> (sites, clusters as lists of reads)"""
    alt = collections.defaultdict(collections.Counter)
    for cov in reads_cov.values():
        for p, b in cov.items():
            if b in "ACGT" and b != ref[p]:
                alt[p][b] += 1
    sites = {p: c.most_common(1)[0][0] for p, c in alt.items() if c.most_common(1)[0][1] >= MIN_SUPPORT}
    vec = {r: {p: ("alt" if cov[p] == sites[p] else "ref" if cov[p] == ref[p] else "other") for p in sites if p in cov}
           for r, cov in reads_cov.items()}
    order = sorted(vec, key=lambda r: (-len(vec[r]), r))
    clusters = []        # [members, {site: Counter(allele)}]

    def cons(cl):
        return {p: c.most_common(1)[0][0] for p, c in cl[1].items()}

    def fit(v, cl):
        cc = cons(cl)
        co = [p for p in v if p in cc and v[p] != "other" and cc[p] != "other"]
        conf = sum(v[p] != cc[p] for p in co)
        return len(co), conf

    for r in order:
        v = vec[r]
        best = None
        for i, cl in enumerate(clusters):
            co, conf = fit(v, cl)
            if conf <= max(0, (co - 9) // 10) and (best is None or co > best[1]):
                best = (i, co)
        if best is None:
            clusters.append([[r], collections.defaultdict(collections.Counter)])
            i = len(clusters) - 1
        else:
            i = best[0]
            clusters[i][0].append(r)
        for p, a in v.items():
            clusters[i][1][p][a] += 1

    def separated(a, b):
        """a site where each cluster's consensus allele differs, supported by >= MIN_SUPPORT reads on each side"""
        ca, cb = cons(a), cons(b)
        return any(ca[p] != cb[p] and a[1][p][ca[p]] >= MIN_SUPPORT and b[1][p][cb[p]] >= MIN_SUPPORT for p in ca if p in cb)

    changed = True
    while changed:
        changed = False
        big = [c for c in clusters if len(c[0]) >= MIN_SUPPORT]
        small = [c for c in clusters if len(c[0]) < MIN_SUPPORT]
        for c in small:                                  # dissolve small clusters into their best compatible big cluster
            for r in c[0]:
                v = vec[r]
                fits = [(fit(v, b), k) for k, b in enumerate(big)]
                fits = [(co, k) for (co, conf), k in fits if conf <= max(0, (co - 9) // 10)]
                if fits:
                    k = max(fits)[1]
                    big[k][0].append(r)
                    for p, a in v.items():
                        big[k][1][p][a] += 1
        clusters = big
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                if not separated(clusters[i], clusters[j]):
                    for p, c in clusters[j][1].items():
                        clusters[i][1][p].update(c)
                    clusters[i][0] += clusters[j][0]
                    del clusters[j]
                    changed = True
                    break
            if changed:
                break
    return sites, [(c[0], cons(c)) for c in clusters], vec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("world"); ap.add_argument("variant"); ap.add_argument("arm")
    a = ap.parse_args()
    w, d = a.world, f"{a.world}/{a.variant}.{a.arm}"
    genome = read_fa(f"{w}/genome.ref.fa")
    copies_fa = read_fa(f"{w}/P.fam.copies.fa")
    reads = read_fq(f"{w}/reads.{a.variant}.fq")
    rows = list(csv.DictReader(open(f"{w}/P.fam.copies.tsv"), delimiter="\t"))
    fam_copies = collections.defaultdict(list)
    for r, (name, seq) in zip(rows, copies_fa.items()):
        fam_copies[r["family_id"]].append((int(r["copy_idx"]), r["chrom"], int(r["start"]), int(r["end"]), seq, name))
    by_host = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{d}/A.assignments.tsv"), delimiter="\t"):
        if r["catalog_copy_idx"] not in ("", "NA", "-"):
            by_host[(r["family_id"], int(r["catalog_copy_idx"]))].append(r["read_name"])
    out_rows = []
    tmp = f"{d}/psvclust.tmp"
    for fam in sorted(fam_copies):
        cps = fam_copies[fam]
        classes = collections.defaultdict(list)          # identical spliced sequences = one host class
        for c in cps:
            classes[c[4]].append(c)
        for cls in classes.values():
            rep = min(cls)                                # the class's lowest copy index is its representative locus
            _, chrom, s, e, _, _ = rep
            names = sorted({n for c in cls for n in by_host.get((fam, c[0]), [])})
            if len(names) < MIN_SUPPORT:
                continue
            lo, hi = max(0, s - PAD), min(len(genome[chrom]), e + PAD)
            ref = genome[chrom][lo:hi].upper()
            cov = align("host", ref, {n: reads[n] for n in names if n in reads}, tmp)
            ccov = align("host", ref, {f"copy{c[0]}": c[4] for c in cps}, tmp)
            sites, clusters, vec = cluster_host(ref, cov)
            for k, (members, cons) in enumerate(sorted(clusters, key=lambda x: -len(x[0]))):
                alts = [p for p, al in cons.items() if al == "alt"]
                explained_by = []
                for cn, cc in ccov.items():
                    ok = all(cc.get(p) == ref_alt for p, ref_alt in ((p, sites[p]) for p in alts) if p in cc) and \
                        not any(p in cc and cc[p] != ref[p] and al == "ref" for p, al in cons.items())
                    if ok and all(p in cc for p in alts):
                        explained_by.append(cn)
                unexpl = [p for p in alts if not any(ccov[cn].get(p) == sites[p] for cn in ccov)]
                covered = sum(1 for p in range(len(ref)) if any(p in cov[m] for m in members[:50]))
                src = collections.Counter(m.split("|")[0] for m in members)
                link = linkage(unexpl, vec) if len(unexpl) >= MIN_ALT else None
                # editing signature: ref A -> G or ref T -> C (either strand), as a fraction of the unexplained alt sites
                ag = sum((ref[p], sites[p]) in (("A", "G"), ("T", "C")) for p in unexpl) / len(unexpl) if unexpl else 0.0
                if explained_by:
                    status = "explained"
                elif len(unexpl) < MIN_ALT or link is None:
                    status = "unresolved"
                else:
                    status = "unexplained" if link >= MIN_LINK else "unlinked"
                out_rows.append(dict(family=fam, host=",".join(str(c[0]) for c in sorted(cls)), cluster=k, n_reads=len(members),
                                     n_sites=len(sites), n_alt=len(alts), n_unexplained_alt=len(unexpl), status=status,
                                     explained_by=",".join(explained_by) or "-", est_div=round(len(unexpl) / covered, 5) if covered else 0,
                                     linkage="NA" if link is None else round(link, 3), editing_frac=round(ag, 3),
                                     sources=",".join(f"{k2}:{v}" for k2, v in src.most_common())))
    for x in (".t.fa", ".q.fa"):
        if os.path.exists(tmp + x):
            os.remove(tmp + x)
    with open(f"{d}/psvclust.tsv", "w") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t", lineterminator="\n")
        wr.writeheader(); wr.writerows(out_rows)
    scen = json.load(open(f"{w}/world.json"))
    for r in out_rows:
        if r["status"] != "explained" or r["family"] in ():
            print(f"{r['family']}\t{scen[r['family']]}\thost {r['host']}\tc{r['cluster']}\t{r['status']}\tn {r['n_reads']}\talt {r['n_alt']}"
                  f"\tunexpl {r['n_unexplained_alt']}\tlink {r['linkage']}\tAG {r['editing_frac']}\tdiv {r['est_div']}\t{r['sources']}")


if __name__ == "__main__":
    main()
