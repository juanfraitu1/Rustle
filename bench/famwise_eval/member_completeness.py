#!/usr/bin/env python3
"""Are the members of known multi-copy families complete, and if not, why? (user, 2026-11-19: look at what StringTie
reports for known families; could tied secondaries seeding there, or reads claimed back from where they wrongly mapped,
restore a member where expected?)

  member_completeness.py --bam SEL.bam --gtf NAME=GTF [...] --out DIR [--truth READS.tsv --tie-split T.tsv]
                         [--species human] [--min-reads 2] [--window 10] [--make-bams]

Members: genes of the truth families of bench/denovo_eval/famtruth.py (connected components of Compara `recent` paralogue
pairs, >= 2 genes; CAT genes on CHM13). A member is FOUND by a GTF when one of its transcripts maps to it (famtruth's
map_copy: most shared exonic bp, same strand). Read evidence per member, from the BAM (every record of every read):
  primary  molecules whose primary record maps to the member
  tied     molecules whose primary is elsewhere and that have a record mapping to the member within --window ms of their best
  claim    (with --truth) molecules whose truth origin locus maps to the member while their primary maps elsewhere (or nowhere)
Class of a member no GTF finds: `reads_there` (primary >= --min-reads: the assembler or the read floor lost it), else
`seedable` (tied >= --min-reads), else `claimable` (claim >= --min-reads, needs --truth), else `nothing` (silent here, or not
in this individual). The same classes are written per GTF for the members that GTF misses.
--make-bams writes DIR/targeted.bam = every primary record + the tied secondaries that map to a member NO gtf finds (seeding
there only) and, with --truth, DIR/claimed.bam = every record, with the claimable molecules' record at their origin member
promoted to primary (the oracle ceiling of claiming back; no rule decides it).
Writes DIR/members.tsv, DIR/summary.tsv, DIR/missed_members.txt."""
import argparse, collections, csv, os, sys
import pysam
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "denovo_eval"))
import famtruth as F

def blocks(r):
    out, pos, s = [], r.reference_start, r.reference_start
    for op, n in r.cigartuples or ():
        if op in (0, 2, 7, 8):
            pos += n
        elif op == 3:
            out.append((s, pos)); pos += n; s = pos
    out.append((s, pos))
    return out

def families(G, species):
    pairs, universe = F.truth_pairs(G, F.PARALOGS, F.levels("recent"))
    fams = [c for c in F.components(pairs) if len(c) >= 2]
    fam_of = {g: i for i, c in enumerate(fams) for g in c}
    return fams, fam_of

def read_bam(path, idx, members):
    """-> {read: [(member or None, ms, primary, record index)]} for reads with >= 1 record on a member"""
    recs = collections.defaultdict(list)
    with pysam.AlignmentFile(path, check_sq=False) as fh:
        for r in fh.fetch(until_eof=True):
            if r.is_unmapped or r.is_supplementary:
                continue
            g = F.map_copy(idx, r.reference_name, "-" if r.is_reverse else "+", blocks(r))
            i = len(recs[r.query_name])
            recs[r.query_name].append((g if g in members else None, r.get_tag("ms"), not r.is_secondary, i))
    return {q: v for q, v in recs.items() if any(x[0] for x in v)}

def origin_members(truth_tsv, tie_split, idx, members):
    """-> {read: set of members its truth origin loci map to} (kinds one / several)"""
    sys.path.insert(0, HERE)
    import records as R
    T = R.load_truth(truth_tsv, tie_split)
    out = {}
    for q, t in T.items():
        if t.kind not in R.JUDGED:
            continue
        ms = set()
        for i in t.origin:
            c, s, e, st, _ = t.loci[i]
            g = F.map_copy(idx, c, st, [(s, e)])
            if g in members:
                ms.add(g)
        if ms:
            out[q] = ms
    return out

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True); ap.add_argument("--gtf", action="append", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--truth"); ap.add_argument("--tie-split"); ap.add_argument("--species", default="human")
    ap.add_argument("--min-reads", type=int, default=2); ap.add_argument("--window", type=float, default=10)
    ap.add_argument("--make-bams", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    G = F.load_genes(a.species)
    fams, fam_of = families(G, a.species)
    members = set(fam_of)
    idx = F.exon_index({g: G[g] for g in members})
    found = {n: F.expressed_from_gtf(p, idx) & members for n, p in (x.split("=", 1) for x in a.gtf)}
    recs = read_bam(a.bam, idx, members)
    prim, tied = collections.Counter(), collections.Counter()
    tied_recs = collections.defaultdict(list)                 # member -> [(read, record index)]
    prim_of = {}
    for q, v in recs.items():
        best = max(x[1] for x in v)
        p = next((x for x in v if x[2]), None)
        pm = p[0] if p else None
        prim_of[q] = pm
        if pm:
            prim[pm] += 1
        for g in {x[0] for x in v if x[0] and not x[2] and best - x[1] <= a.window and x[0] != pm}:
            tied[g] += 1
        for x in v:
            if x[0] and not x[2] and best - x[1] <= a.window and x[0] != pm:
                tied_recs[x[0]].append((q, x[3]))
    claim, claim_reads = collections.Counter(), collections.defaultdict(list)
    if a.truth:
        for q, ms in origin_members(a.truth, a.tie_split, idx, members).items():
            if q in prim_of and prim_of[q] not in ms:
                for g in ms:
                    claim[g] += 1; claim_reads[g].append(q)
    def cls(g):
        k = a.min_reads
        return ("reads_there" if prim[g] >= k else "seedable" if tied[g] >= k else
                "claimable" if claim[g] >= k else "nothing")
    any_found = set().union(*found.values())
    names = list(found)
    with open(f"{a.out}/members.tsv", "w") as fh:
        fh.write("member\tname\tfamily\tfamily_size\t" + "\t".join(f"found_{n}" for n in names) + "\tprimary\ttied\tclaim\tclass_if_missed\n")
        for g in sorted(members):
            fh.write(f"{g}\t{G[g].get('name', g)}\t{fam_of[g]}\t{len(fams[fam_of[g]])}\t" + "\t".join(str(int(g in found[n])) for n in names)
                     + f"\t{prim[g]}\t{tied[g]}\t{claim[g]}\t{cls(g)}\n")
    L = [f"families {len(fams)}, members {len(members)}; members with >= 1 primary molecule in this BAM {sum(prim[g] > 0 for g in members)}"]
    expressed = {g for g in members if prim[g] >= a.min_reads}
    for n in names:
        miss = members - found[n]
        c = collections.Counter(cls(g) for g in miss)
        L.append(f"{n}: found {len(found[n])} (of members with >= {a.min_reads} primaries: {len(found[n] & expressed)}/{len(expressed)}); "
                 f"missed {len(miss)} = " + ", ".join(f"{k} {c[k]}" for k in ("reads_there", "seedable", "claimable", "nothing")))
    miss_all = members - any_found
    c = collections.Counter(cls(g) for g in miss_all)
    L.append(f"missed by every GTF {len(miss_all)} = " + ", ".join(f"{k} {c[k]}" for k in ("reads_there", "seedable", "claimable", "nothing")))
    fams_complete = {n: sum(1 for f in fams if f <= found[n] | (members - expressed)) for n in names}
    L.append("families whose expressed members (>= min primaries) are all found: " + ", ".join(f"{n} {v}/{len(fams)}" for n, v in fams_complete.items()))
    with open(f"{a.out}/summary.tsv", "w") as fh:
        fh.write("\n".join(L) + "\n")
    with open(f"{a.out}/missed_members.txt", "w") as fh:
        for g in sorted(miss_all, key=lambda g: (cls(g), G[g].get("name", g))):
            fh.write(f"{cls(g)}\t{G[g].get('name', g)}\tprimary {prim[g]}\ttied {tied[g]}\tclaim {claim[g]}\n")
    if a.make_bams:
        target = {g for g in miss_all if cls(g) == "seedable"}
        keep_tied = {(q, i) for g in target for q, i in tied_recs[g]}
        promote = {}
        for g in miss_all:
            for q in claim_reads.get(g, ()):
                rs = recs[q]
                cand = [x for x in rs if x[0] == g]
                if cand:
                    promote[q] = max(cand, key=lambda x: x[1])[3]
        with pysam.AlignmentFile(a.bam, check_sq=False) as fh:
            wt = pysam.AlignmentFile(f"{a.out}/targeted.unsorted.bam", "wb", template=fh)
            wc = pysam.AlignmentFile(f"{a.out}/claimed.unsorted.bam", "wb", template=fh) if a.truth else None
            seen = collections.Counter()
            for r in fh.fetch(until_eof=True):
                if r.is_unmapped or r.is_supplementary:
                    continue
                i = seen[r.query_name]; seen[r.query_name] += 1
                if not r.is_secondary or (r.query_name, i) in keep_tied:
                    wt.write(r)
                if wc is not None:
                    if r.query_name in promote:
                        if i == promote[r.query_name]:
                            r.flag &= ~0x100
                        elif not r.is_secondary:
                            r.flag |= 0x100
                    wc.write(r)
            for w, name in ((wt, "targeted"), (wc, "claimed")):
                if w is None:
                    continue
                w.close()
                pysam.sort("-o", f"{a.out}/{name}.bam", f"{a.out}/{name}.unsorted.bam"); pysam.index(f"{a.out}/{name}.bam")
                os.remove(f"{a.out}/{name}.unsorted.bam")
        print(f"bams: targeted keeps {len(keep_tied)} tied records at {len(target)} seedable members; claimed promotes {len(promote)} molecules")
    print("\n".join(L))

if __name__ == "__main__":
    main()
