#!/usr/bin/env python3
"""Question 4 of docs/GOAL_multicopy_families_2026-11-19.md: borrowed structure as a LABELLED class, never expression.

  borrow.py --gtf ASSEMBLY.gtf --bam SEL.bam --fasta GENOME.fa --out PREFIX [--threads 8] [--window 10]

For every molecule with a record at another place (a paralog) within --window ms of its best, the assembly's multi-exon
transcripts at its primary's place are projected onto that paralog (minimap2 splice:hq, cov >= 0.9, id >= 0.85; the
projection of row 1519 / 1820). A projected chain is KEPT only when (a) a record of some molecule at that paralog carries a
chain that is a contiguous run of the projected chain, every junction within --tol bp (read evidence for the structure at
the copy; projected junctions can sit a few bases off), and (b) no transcript of the assembly at that paralog holds that
read chain within --tol bp (secondary records often carry junction-shifted copies of chains already assembled). The kept
chain's junctions are snapped onto the evidence read's where they matched. Kept chains are written as PREFIX.borrowed.gtf with
`fw_class "borrowed_structure"`; PREFIX.gtf = the assembly plus them. No truth is read here."""
import argparse, collections, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import records as R
from undecidable import project, read_gtf_tx, chain

def run_at(big, small, tol):
    """offset at which `small` is a contiguous run of `big` with every junction end within tol bp, else None"""
    k = len(small)
    for i in range(len(big) - k + 1):
        if all(abs(x[0] - y[0]) <= tol and abs(x[1] - y[1]) <= tol for x, y in zip(big[i:i + k], small)):
            return i
    return None

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True); ap.add_argument("--bam", required=True); ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--threads", type=int, default=8); ap.add_argument("--window", type=float, default=10)
    ap.add_argument("--tol", type=int, default=5)
    a = ap.parse_args(argv)
    tx = read_gtf_tx(a.gtf)
    multi = {t: v for t, v in tx.items() if len(v["exons"]) > 1}
    # index of assembly transcripts by contig / strand / 100 kb bin
    tb = collections.defaultdict(list)
    for t, v in multi.items():
        for k in range(v["exons"][0][0] // 100000, v["exons"][-1][1] // 100000 + 1):
            tb[(v["chrom"], v["strand"], k)].append(t)
    def at(c, st, s, e):
        return {t for k in range(s // 100000, max(s, e - 1) // 100000 + 1) for t in tb.get((c, st, k), ())
                if multi[t]["exons"][0][0] < e and s < multi[t]["exons"][-1][1]}
    asm_keys = {(v["chrom"], v["strand"], chain(v["exons"])) for v in multi.values()}
    by_cs = collections.defaultdict(list)
    for c, st, ch in asm_keys:
        by_cs[(c, st)].append(ch)
    def held(c, st, e):
        return any(run_at(ch, e, a.tol) is not None for ch in by_cs.get((c, st), ()))
    recs = R.load_records(a.bam)
    # read evidence at each place, and the (template -> target copy) pairs
    copies, cid, pairs, evidence = [], {}, collections.defaultdict(set), collections.defaultdict(set)
    for q, rs in recs.items():
        best = max(r.ms for r in rs)
        p = next((r for r in rs if r.primary), None)
        if p is None:
            continue
        temps = at(p.chrom, p.strand, p.start, p.end)
        for r in rs:
            if r.primary or best - r.ms > a.window or (r.chrom == p.chrom and r.strand == p.strand and r.start < p.end and p.start < r.end):
                continue
            if r.chain and held(r.chrom, r.strand, r.chain):
                continue                                              # (b): the copy already holds this read chain
            key = (r.chrom, r.strand, r.start // 2000)
            if key not in cid:
                cid[key] = f"t{len(copies)}"
                copies.append(dict(id=cid[key], chrom=r.chrom, strand=r.strand, start=r.start, end=r.end))
            c = copies[int(cid[key][1:])]
            c["start"], c["end"] = min(c["start"], r.start), max(c["end"], r.end)
            for t in temps:
                pairs[t].add(c["id"])
            if r.chain:
                evidence[c["id"]].add((r.chrom, r.strand, r.chain))
    proj = a.out + ".projected.gtf"
    n, nc = project(pairs, multi, copies, a.fasta, proj, a.threads)
    kept, lines, done = 0, [], set()
    pt = read_gtf_tx(proj)
    gid = {}
    for l in open(proj):
        m = re.search(r'transcript_id "([^"]+)"; template', l)
        if m:
            gid[m.group(1)] = re.search(r'gene_id "([^"]+)"', l).group(1)
    for t, v in pt.items():
        ch = chain(v["exons"]); c = gid.get(t)
        hit = next(((e, i) for e in evidence.get(c, ()) if e[0] == v["chrom"] and e[1] == v["strand"]
                    for i in [run_at(ch, e[2], a.tol)] if i is not None), None)
        if hit is None:
            continue
        e, i = hit
        ch = ch[:i] + e[2] + ch[i + len(e[2]):]                     # snap onto the read's junctions
        if (v["chrom"], v["strand"], ch) in asm_keys or (v["chrom"], v["strand"], ch) in done:
            continue
        done.add((v["chrom"], v["strand"], ch))
        ex = [(v["exons"][0][0], ch[0][0])] + [(ch[j][1], ch[j + 1][0]) for j in range(len(ch) - 1)] + [(ch[-1][1], v["exons"][-1][1])]
        if any(s >= x for s, x in ex):
            continue
        v = dict(v, exons=ex)
        kept += 1
        attr = f'gene_id "borrow.{kept}"; transcript_id "borrow.{kept}"; fw_class "borrowed_structure";'
        lines.append(f"{v['chrom']}\tborrow\ttranscript\t{v['exons'][0][0] + 1}\t{v['exons'][-1][1]}\t.\t{v['strand']}\t.\t{attr}\n")
        lines += [f"{v['chrom']}\tborrow\texon\t{s + 1}\t{e}\t.\t{v['strand']}\t.\t{attr}\n" for s, e in v["exons"]]
    open(a.out + ".borrowed.gtf", "w").writelines(lines)
    with open(a.out + ".gtf", "w") as fo:
        fo.write(open(a.gtf).read()); fo.writelines(lines)
    print(f"borrow: {len(pairs)} templates onto {nc} copies, {n} projected chains, {kept} kept as borrowed_structure "
          f"(read evidence at the copy, not already held there)")

if __name__ == "__main__":
    main()
