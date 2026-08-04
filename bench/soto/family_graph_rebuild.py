#!/usr/bin/env python3
"""Rebuild the per-family similarity graph on an HONEST substrate.

Two defects in family_similarity_graph.py made its NPIP graph unusable as evidence:

  1. It selected exons by RefSeq `gene=` symbol equality and, when no record for that symbol
     overlapped the member window, SILENTLY fell back to the whole genomic span (introns
     included). On Soto ID_154 that fired on 4 of 14 members (NPIPP1, NPIPB7, NPIPB8,
     NPIPB14P) -- and 3 of those 4 are among the 6 isolated vertices. Because coverage is
     normalised by the SHORTER sequence, comparing an 817 bp exon-sum against a 21 kb genomic
     span depresses coverage mechanically, so part of the reported sparsity was an artifact of
     the substrate rather than a property of the loci.
  2. Its edge predicate did not match the pipeline's: identity came from matches/alignment_len
     rather than 1 - de, it kept only the max-identity record per pair, and it omitted
     --no-long-join.

This rebuild fixes both:

  * Exons are taken by COORDINATE OVERLAP with the member window, regardless of gene symbol,
    so a fusion name (PKD1P6-NPIPP1), a swapped window, or a renamed pseudogene cannot cause a
    silent substrate change. A member with zero overlapping exons is REPORTED, never patched.
  * The predicate mirrors src/rustle/vg_family nucleotide_edges: minimap2 -c -X --no-long-join
    over two tiers (asm20, then -k 11 -w 5), identity = 1 - de:f:, coverage =
    (q_end - q_start) / min(q_len, t_len), and an edge exists iff ANY single record from ANY
    tier satisfies both thresholds.
"""
import argparse, gzip, json, os, subprocess, sys, tempfile
from collections import defaultdict

REF = "/mnt/c/Users/jfris/Desktop/Reference"
FASTA = os.environ.get("SOTO_FASTA", f"{REF}/chm13v2.0.fa")
GFF = os.environ.get("SOTO_GFF", f"{REF}/chm13v2.0_RefSeq_full.gff.gz")
BED = os.environ.get("SOTO_BED", "bench/soto/80_fams.chr.bed")
SAMTOOLS = os.environ.get("SAMTOOLS", "samtools")
MINIMAP2 = os.environ.get("MINIMAP2", "minimap2")


def members(fam):
    out = []
    for ln in open(BED):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        name = f[3]
        if not name.endswith("|" + fam):
            continue
        out.append((name.split("|")[0], f[0], int(f[1]), int(f[2])))
    return sorted(out, key=lambda x: (x[1], x[2]))


def exons_in(chrom, start, end, exclude=(), gene_re=None):
    """Every annotated exon overlapping [start,end) on chrom, merged. Symbol-agnostic.

    `exclude` holds other members' windows of the same family. Symbol-agnostic selection is
    immune to fusion names and swapped symbols, but it introduces a hazard the symbol-based
    version could not have: when two truth windows OVERLAP, both members are handed the same
    genomic bases and the resulting edge is a sequence aligned against a superset of itself
    (identity 1.0000, coverage 1.0000 -- the visible signature). 6 of 83 families in
    80_fams.chr.bed have overlapping member windows; in ID_127 the only reported family was
    entirely this artifact, and in ID_400 it wrongly rescued NBPF26 from being isolated.
    Shared bases are therefore removed from BOTH members rather than credited to either.
    """
    try:
        raw = subprocess.run(["tabix", GFF, f"{chrom}:{start+1}-{end}"],
                             capture_output=True, text=True, check=True).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        raw = ""
        with gzip.open(GFF, "rt") as fh:
            for ln in fh:
                if ln.startswith("#"):
                    continue
                f = ln.split("\t")
                if f[0] != chrom or f[2] != "exon":
                    continue
                if int(f[4]) > start and int(f[3]) - 1 < end:
                    raw += ln
    iv = []
    for ln in raw.splitlines():
        f = ln.split("\t")
        if len(f) < 8 or f[2] != "exon":
            continue
        if gene_re is not None:
            g = [x[5:] for x in f[8].split(";") if x.startswith("gene=")]
            if not g or not gene_re.search(g[0]):
                continue
        a, b = max(int(f[3]) - 1, start), min(int(f[4]), end)
        if b > a:
            iv.append((a, b))
    iv.sort()
    merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    if not exclude:
        return merged
    out = []
    for a, b in merged:
        pieces = [(a, b)]
        for xs, xe in exclude:
            nxt = []
            for pa, pb in pieces:
                if xe <= pa or xs >= pb:
                    nxt.append((pa, pb))
                    continue
                if pa < xs:
                    nxt.append((pa, min(xs, pb)))
                if pb > xe:
                    nxt.append((max(xe, pa), pb))
            pieces = nxt
        out.extend([list(p) for p in pieces if p[1] > p[0]])
    return out


def fetch(chrom, iv):
    if not iv:
        return ""
    regions = [f"{chrom}:{a+1}-{b}" for a, b in iv]
    seq = []
    for i in range(0, len(regions), 200):
        r = subprocess.run([SAMTOOLS, "faidx", FASTA] + regions[i:i+200],
                           capture_output=True, text=True).stdout
        for blk in r.split(">")[1:]:
            seq.append("".join(blk.splitlines()[1:]))
    return "".join(seq).upper()


def paf_edges(fa, preset, min_id, min_cov):
    """One minimap2 tier. Returns {(a,b): (identity, coverage)} best-by-coverage per pair,
    and the set of pairs where SOME single record satisfied both thresholds."""
    cmd = [MINIMAP2, "-c", "-X", "--no-long-join", "-t", "4"] + preset + [fa, fa]
    p = subprocess.run(cmd, capture_output=True, text=True)
    best, passed = {}, set()
    for ln in p.stdout.splitlines():
        f = ln.split("\t")
        if len(f) < 12:
            continue
        q, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
        t, tlen = f[5], int(f[6])
        if q == t:
            continue
        de = None
        for tag in f[12:]:
            if tag.startswith("de:f:"):
                de = float(tag[5:])
                break
        ident = (1.0 - de) if de is not None else (int(f[9]) / max(int(f[10]), 1))
        cov = (qe - qs) / max(min(qlen, tlen), 1.0)
        key = (min(q, t), max(q, t))
        ok = ident >= min_id and cov >= min_cov
        if ok:
            passed.add(key)
        # Prefer a record that actually satisfies the predicate, so the reported
        # identity/coverage always describe the record that made the edge.
        prev = best.get(key)
        if prev is None or (ok and not prev[2]) or (ok == prev[2] and cov > prev[1]):
            best[key] = (ident, cov, ok)
    return best, passed


def components(nodes, edges):
    p = {v: v for v in nodes}

    def find(x):
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    for a, b in edges:
        p[find(a)] = find(b)
    g = defaultdict(list)
    for v in nodes:
        g[find(v)].append(v)
    return sorted((sorted(c) for c in g.values()), key=lambda c: (-len(c), c[0]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("families", nargs="+", help="Soto family IDs, e.g. ID_154 ID_354")
    ap.add_argument("--min-identity", type=float, default=0.80)
    ap.add_argument("--min-coverage", type=float, default=0.50)
    # Sensitive-tier identity floor. Default = min_identity, which is what the pipeline does
    # when --min-identity is given explicitly: homology_refine_params (denovo_pipeline.rs:2442-2449)
    # sets p.sensitive_identity = mi, commented "BOTH tiers -> effective floor = mi". Only the
    # pipeline DEFAULT differs (asm20 0.80 / sensitive 0.60), and refine_families_exon_sum
    # hardcodes 0.70 -- pass this flag to model either. Numerically inert on ID_154: the union is
    # 32 edges at 0.50, 0.60, 0.70 and 0.80 alike, because coverage is the binding constraint.
    ap.add_argument("--sensitive-identity", type=float, default=None)
    # Restrict the substrate to exons of genes matching this pattern. Symbol-agnostic selection
    # is immune to fusion names and swapped windows but admits neighbouring genes; this filter
    # removes them WITHOUT reintroducing the symbol-lookup fallback, because a member whose exons
    # all fail the filter is excluded and reported, never silently replaced by its genomic span.
    ap.add_argument("--gene-regex", default=None)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    import re as _re
    gene_re = _re.compile(a.gene_regex) if a.gene_regex else None
    result = {}
    for fam in a.families:
        mem = members(fam)
        if not mem:
            print(f"[{fam}] no members in {BED}", file=sys.stderr)
            continue
        seqs, meta, missing = {}, {}, []
        for name, c, s, e in mem:
            others = [(s2, e2) for n2, c2, s2, e2 in mem if n2 != name and c2 == c
                      and s2 < e and e2 > s]
            if others:
                print(f"[{fam}] {name}: window overlaps {len(others)} sibling window(s); "
                      f"shared bases removed from both", file=sys.stderr)
            iv = exons_in(c, s, e, exclude=others, gene_re=gene_re)
            if not iv:
                missing.append(name)
                meta[name] = {"span": e - s, "exon_sum": 0, "n_exon": 0, "substrate": "NO EXONS"}
                continue
            sq = fetch(c, iv)
            if not sq:
                raise SystemExit(f"[{fam}] {name}: samtools returned no sequence for "
                                 f"{len(iv)} interval(s) -- refusing to emit an empty node")
            seqs[name] = sq
            meta[name] = {"span": e - s, "exon_sum": len(sq), "n_exon": len(iv),
                          "substrate": "exon-sum"}
        if missing:
            print(f"[{fam}] EXCLUDED (no annotated exon in window): {missing}", file=sys.stderr)

        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
            for k, v in seqs.items():
                fh.write(f">{k}\n{v}\n")
            fa = fh.name

        allbest, passed = {}, set()
        sens = a.min_identity if a.sensitive_identity is None else a.sensitive_identity
        for preset, tier_id in ((["-x", "asm20"], a.min_identity),
                                (["-k", "11", "-w", "5"], sens)):
            best, ps = paf_edges(fa, preset, tier_id, a.min_coverage)
            passed |= ps
            for k, v in best.items():
                prev = allbest.get(k)
                if prev is None or (v[2] and not prev[2]) or (v[2] == prev[2] and v[1] > prev[1]):
                    allbest[k] = v
        os.unlink(fa)

        nodes = sorted(seqs)
        n = len(nodes)
        npairs = n * (n - 1) // 2
        cc = components(nodes, sorted(passed))
        nontrivial = [c for c in cc if len(c) >= 2]
        iso = [c[0] for c in cc if len(c) == 1]

        print(f"\n=== {fam} ===")
        print(f"  nodes {n} (excluded {len(missing)}) · pairs aligned "
              f"{len(allbest)}/{npairs} · edges passing id>={a.min_identity} AND "
              f"cov>={a.min_coverage}: {len(passed)}/{npairs}")
        idv = sorted(v[0] for v in allbest.values())
        cvv = sorted(v[1] for v in allbest.values())
        if idv:
            print(f"  median identity {idv[len(idv)//2]:.3f} · median coverage {cvv[len(cvv)//2]:.3f}")
        print(f"  connected components: {len(cc)}  sizes {[len(c) for c in cc]}")
        print(f"    with >=2 loci (become families): {len(nontrivial)} -> {nontrivial}")
        print(f"    singletons (dropped by the >=2-loci rule): {iso}")
        for c in cc:
            if len(c) >= 3:
                m = sum(1 for (x, y) in passed if x in c and y in c)
                print(f"    density of {len(c)}-part: 2*{m}/({len(c)}*{len(c)-1}) = "
                      f"{2*m/(len(c)*(len(c)-1)):.3f}")
        result[fam] = {
            "meta": meta, "excluded": missing,
            "edges": [{"a": k[0], "b": k[1], "identity": round(v[0], 4),
                       "coverage": round(v[1], 4),
                       "edge": bool(k in passed)} for k, v in sorted(allbest.items())],
            "components": cc, "families": nontrivial, "singletons": iso,
        }

    if a.out:
        json.dump(result, open(a.out, "w"), indent=1)
        print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
