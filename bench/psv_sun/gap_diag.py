#!/usr/bin/env python3
"""What are the PSV columns that no local alignment supports? (docs/archive/2026-10/merged/PSV_SUN_PILEUP_2026-10-12.md, "different exon content")

For every copy X vs copy 0 of a family, in O2's frame (a flipped copy reverse-complemented), the colinear chain of
`minimap2 -cx asm20` pieces splits copy 0 into ALIGNED blocks and unaligned stretches. Each unaligned stretch of copy 0 is
  terminal_overhang   before the first / after the last block, X has <= SHORT bases there (UTR / extent difference)
  ref_only_insertion  between blocks, X has <= SHORT bases (an exon or exon part present in copy 0 only)
  both_homologous     both have sequence there and Smith-Waterman beats every one of NSHUF shuffles (divergent homolog
                      that asm20 missed)
  both_unrelated      both have sequence, SW does not beat the shuffles (alternative / non-homologous exons)
  no_alignment        no piece at all (the whole pair; same SW test, reported separately)
Each O2 decisive (copy X, column) pair with an allele is counted in the class of copy 0's offset for X.

    PSV_SUN_RUN=<dump prefix> python3 gap_diag.py OUTDIR      (needs bin/swscore)
"""
import collections, json, os, random, subprocess, sys
sys.path.insert(0, os.path.dirname(__file__))
import psv_sun as P

SHORT, NSHUF = 10, 20
SWSCORE = os.environ.get("SWSCORE", f"{P.OUT}/bin/swscore")


def chain(recs):
    """Colinear, non-overlapping '+' pieces in target order (greedy by alignment score)."""
    keep = []
    for r in sorted((r for r in recs if r["tp"] == "P" and r["strand"] == "+"), key=lambda r: -r["AS"]):
        if all(r["te"] <= k["ts"] and r["qe"] <= k["qs"] or r["ts"] >= k["te"] and r["qs"] >= k["qe"] for k in keep):
            keep.append(r)
    return sorted(keep, key=lambda r: r["ts"])


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    tab, seqs = P.copy_table(), P.copies_fa(f"{P.SUBSET}/copies.fa")
    cols, copies = P.dump()
    rng = random.Random(11)
    segs, sw_jobs = [], []  # segs: (fam, x, cls, a0, a1); sw_jobs: (seg index, kind, seqA, seqB)
    for fam, crs in sorted(copies.items()):
        t0 = tab.get(crs[0]["copy_tid"])
        if t0 is None:
            continue
        ref = seqs[(fam, t0["copy_idx"])]
        qs = {}
        for i, c in enumerate(crs[1:], 1):
            t = tab.get(c["copy_tid"])
            if t is None:
                continue
            s = seqs[(fam, t["copy_idx"])]
            qs[i] = P.rc(s) if c["strand"] != t["strand"] else s  # frame strand != annotated strand => flipped
        tq, qq = f"{outdir}/t.fa", f"{outdir}/q.fa"
        open(tq, "w").write(f">0\n{ref}\n")
        with open(qq, "w") as fh:
            for i, s in qs.items():
                fh.write(f">{i}\n{s}\n")
        byq = collections.defaultdict(list)
        for r in P.run_paf(["minimap2", "-cx", "asm20", "--cs", "-t", "1", tq, qq]):
            byq[int(r["qname"])].append(r)
        for i, s in qs.items():
            ch = chain(byq.get(i, []))
            for k in ch:
                segs.append((fam, i, "aligned", k["ts"], k["te"]))
            if not ch:
                segs.append((fam, i, "no_alignment", 0, len(ref)))
                sw_jobs.append((len(segs) - 1, ref, s))
                continue
            gaps = [(0, ch[0]["ts"], 0, ch[0]["qs"], True)]
            gaps += [(a["te"], b["ts"], a["qe"], b["qs"], False) for a, b in zip(ch, ch[1:])]
            gaps.append((ch[-1]["te"], len(ref), ch[-1]["qe"], len(s), True))
            for a0, a1, q0, q1, term in gaps:
                if a1 <= a0:
                    continue
                if q1 - q0 <= SHORT:
                    segs.append((fam, i, "terminal_overhang" if term else "ref_only_insertion", a0, a1))
                else:
                    segs.append((fam, i, "both?", a0, a1))
                    sw_jobs.append((len(segs) - 1, ref[a0:a1], s[q0:q1]))
    # Smith-Waterman: observed vs NSHUF mononucleotide shuffles of copy 0's stretch
    with open(f"{outdir}/sw_pairs.tsv", "w") as fh:
        for si, a, b in sw_jobs:
            fh.write(f"{si}:obs\t{a}\t{b}\n")
            for k in range(NSHUF):
                sh = list(a)
                rng.shuffle(sh)
                fh.write(f"{si}:s{k}\t{''.join(sh)}\t{b}\n")
    with open(f"{outdir}/sw_pairs.tsv") as fi, open(f"{outdir}/sw_scores.tsv", "w") as fo:
        subprocess.run([SWSCORE], stdin=fi, stdout=fo, check=True)
    obs, null = {}, collections.defaultdict(list)
    for line in open(f"{outdir}/sw_scores.tsv"):
        key, sc = line.split()
        si, kind = key.split(":")
        (obs.__setitem__(int(si), int(sc)) if kind == "obs" else null[int(si)].append(int(sc)))
    for si in obs:
        fam, i, cls, a0, a1 = segs[si]
        hom = obs[si] > max(null[si])
        if cls == "both?":
            cls = "both_homologous" if hom else "both_unrelated"
        else:
            cls = "no_alignment_homologous" if hom else "no_alignment_unrelated"
        segs[si] = (fam, i, cls, a0, a1)
    # count O2's decisive (copy X, column) pairs by the class of copy 0's offset
    index = collections.defaultdict(list)
    for fam, i, cls, a0, a1 in segs:
        index[(fam, i)].append((a0, a1, cls))
    count, bp = collections.Counter(), collections.Counter()
    labels = []  # per pair: (family, copy 0's genome position of the column, copy X, class, X's allele)
    for fam, i, cls, a0, a1 in segs:
        bp[cls] += a1 - a0
    for fam, crs in copies.items():
        t0 = tab.get(crs[0]["copy_tid"])
        if t0 is None:
            continue
        ex0 = P.parse_exons(t0["exons"])
        for j in range(max(len(c["alleles"]) for c in crs)):
            if not P.decisive(crs, j):
                continue
            off = P.tx_offset(ex0, t0["strand"], cols[fam].get(j, -1))
            for i, c in enumerate(crs[1:], 1):
                if j >= len(c["alleles"]) or c["alleles"][j] == ".":
                    continue
                cls = "unknown" if off is None else next((k for a0, a1, k in index.get((fam, i), []) if a0 <= off < a1), "unknown")
                count[cls] += 1
                labels.append({"family_id": fam, "gpos0": cols[fam].get(j, -1), "copy": i, "cls": cls, "allele": c["alleles"][j]})
    P.write_tsv(f"{outdir}/pair_labels.tsv", labels, ["family_id", "gpos0", "copy", "cls", "allele"])
    tot = sum(count.values())
    res = {"pairs": tot, "by_class": {k: [v, round(v / tot, 4)] for k, v in count.most_common()}, "copy0_bp_by_class": dict(bp)}
    json.dump(res, open(f"{outdir}/gap_diag.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
