#!/usr/bin/env python3
"""Yardstick of docs/archive/2026-10/merged/PREREG_psv_semiglobal_engine_2026-10-12.md: score PSV engines' DIFFERENCE pairs.

    engine_eval.py COPIES.tsv COPIES.fa OUTDIR NAME=PREFIX [NAME=PREFIX ...]
PREFIX = a `utilities psv-discover` output prefix (needs .psv_copies.tsv and .psv_pos.tsv).

Per family, copy 0 vs each copy X in copy 0's frame (X reverse-complemented when its dumped strand differs from its
catalog strand), the colinear chain of `minimap2 -cx asm20` pieces — independent of every engine. A difference pair is
(X, column) with X's allele != copy 0's. Copy 0's offset of the column:
  inside a piece  -> R if copy 0's offset lies in an internally repeated part of copy 0, or minimap2's X position in one
                     of X (off-diagonal self-alignment, `minimap2 -cx asm20 -X`; Amendment A1: positional homology is
                     ambiguous there — a tandem period or a duplicated segment aligns equally well twice); else
                     V if the engine's X genome position equals minimap2's, else M
  outside         -> the unaligned stretch around it (copy 0's and X's sequence between the flanking pieces):
                     H if Smith-Waterman beats all NSHUF shuffles (divergent homolog), else N (unrelated; also when X
                     has <= SHORT bases there or there is no piece at all and SW fails)
Writes OUTDIR/engine_eval.json and prints it.
"""
import collections, json, os, random, subprocess, sys
sys.path.insert(0, os.path.dirname(__file__))
import psv_sun as P

SHORT, NSHUF = 10, 20
SWSCORE = os.environ.get("SWSCORE", f"{P.OUT}/bin/swscore")


def chain(recs):
    keep = []
    for r in sorted((r for r in recs if r["tp"] == "P" and r["strand"] == "+"), key=lambda r: -r["AS"]):
        if all(r["te"] <= k["ts"] and r["qe"] <= k["qs"] or r["ts"] >= k["te"] and r["qs"] >= k["qe"] for k in keep):
            keep.append(r)
    return sorted(keep, key=lambda r: r["ts"])


def self_repeats(seq, tmp):
    """Offsets of `seq` covered by an off-diagonal self-alignment (an internal repeat)."""
    open(f"{tmp}/self.fa", "w").write(f">s\n{seq}\n")
    cov = set()
    for r in P.run_paf(["minimap2", "-cx", "asm20", "-X", "-N", "50", "-p", "0.1", "-t", "1", f"{tmp}/self.fa", f"{tmp}/self.fa"]):
        if (r["qs"], r["qe"]) != (r["ts"], r["te"]):
            cov.update(range(r["qs"], r["qe"]))
            cov.update(range(r["ts"], r["te"]))
    return cov


def main(copies_tsv, copies_fa, outdir, runs):
    os.makedirs(outdir, exist_ok=True)
    tab = {}
    for r in P.tsv(copies_tsv):
        tab[(r["family_id"], int(r["copy_idx"]))] = r
    fam_rows = collections.defaultdict(list)
    for r in P.tsv(copies_tsv):
        fam_rows[r["family_id"]].append(r)  # file order = the harness's copy order
    seqs = P.copies_fa(copies_fa)
    engines = {}
    for name, prefix in runs:
        alle, strand = {}, {}
        for r in P.tsv(prefix + ".psv_copies.tsv"):
            alle[(r["family_id"], int(r["copy_index"]))] = r["alleles"]
            strand[(r["family_id"], int(r["copy_index"]))] = r["strand"]
        pos = collections.defaultdict(dict)  # (fam, col) -> {copy: genome pos}
        for r in P.tsv(prefix + ".psv_pos.tsv"):
            pos[(r["family_id"], int(r["col_index"]))][int(r["copy_index"])] = int(r["genome_pos"])
        engines[name] = (alle, strand, pos)
    # families present in EVERY compared output (an engine may lack some, e.g. one it could not align)
    fam_sets = [{f for (f, _) in alle} for alle, _, _ in engines.values()]
    fams = sorted(set.intersection(*fam_sets))
    if os.environ.get("ONLY"):  # a family subset (one id per line), to run large catalogs in chunks
        keep = set(open(os.environ["ONLY"]).read().split())
        fams = [f for f in fams if f in keep]
    import time
    first_strand = next(iter(engines.values()))[1]
    # per family and copy X: minimap2 pieces, the stretch index, a copy-0-offset -> X-genome map
    stretch_jobs, stretch_of = [], {}
    pairs = {name: collections.Counter() for name in engines}
    pending = {name: [] for name in engines}  # (stretch id) per outside pair, resolved after SW
    rng = random.Random(11)
    tmp = os.path.join(outdir, "tmp")
    os.makedirs(tmp, exist_ok=True)
    for fam in fams:
        t_fam = time.time()
        rows = fam_rows[fam]
        info = []
        for k, r in enumerate(rows):
            s = seqs[(fam, r["copy_idx"])]
            flipped = first_strand.get((fam, k), r["strand"]) != r["strand"]
            info.append({"ex": P.parse_exons(r["exons"]), "strand": r["strand"], "seq": s, "flip": flipped,
                         "oseq": P.rc(s) if flipped else s})
        ref = info[0]["oseq"]
        open(f"{tmp}/t.fa", "w").write(f">0\n{ref}\n")
        with open(f"{tmp}/q.fa", "w") as fh:
            for k in range(1, len(info)):
                fh.write(f">{k}\n{info[k]['oseq']}\n")
        byq = collections.defaultdict(list)
        for r in P.run_paf(["minimap2", "-cx", "asm20", "--cs", "-t", "1", f"{tmp}/t.fa", f"{tmp}/q.fa"]):
            byq[int(r["qname"])].append(r)
        reps = [self_repeats(x["oseq"], tmp) for x in info]  # internal repeats per copy, in its oriented frame
        mm = {}  # (k) -> (pieces, map t_off -> X genome pos, gaps, map t_off -> X oriented offset)
        for k in range(1, len(info)):
            x = info[k]
            ch = chain(byq.get(k, []))
            tmap, omap = {}, {}
            for r in ch:
                for qo, t in P.anchor_pairs(r):
                    omap[t] = qo
                    own = len(x["oseq"]) - 1 - qo if x["flip"] else qo
                    g = P.genome_pos(x["ex"], x["strand"], own)
                    if g is not None:
                        tmap[t] = g
            if ch:
                gaps = [(0, ch[0]["ts"], 0, ch[0]["qs"])] + \
                       [(a["te"], b["ts"], a["qe"], b["qs"]) for a, b in zip(ch, ch[1:])] + \
                       [(ch[-1]["te"], len(ref), ch[-1]["qe"], len(x["oseq"]))]
            else:
                gaps = [(0, len(ref), 0, len(x["oseq"]))]
            mm[k] = (ch, tmap, gaps, omap)
        ex0, st0 = info[0]["ex"], info[0]["strand"]
        for name, (alle, _, pos) in engines.items():
            a0 = alle.get((fam, 0), "")
            for j in range(len(a0)):
                g0 = pos.get((fam, j), {}).get(0)
                off = None if g0 is None else P.tx_offset(ex0, st0, g0)
                for k in range(1, len(info)):
                    ak = alle.get((fam, k), "")
                    if j >= len(ak) or ak[j] == "." or a0[j] == "." or ak[j] == a0[j]:
                        continue
                    if off is None:
                        pairs[name]["unknown"] += 1
                        continue
                    ch, tmap, gaps, omap = mm[k]
                    if any(r["ts"] <= off < r["te"] for r in ch):
                        if off in reps[0] or omap.get(off) in reps[k]:
                            pairs[name]["R"] += 1
                        else:
                            pairs[name]["V" if tmap.get(off) == pos[(fam, j)].get(k) else "M"] += 1
                        continue
                    gap = next(gp for gp in gaps if gp[0] <= off < gp[1])
                    key = (fam, k, gap)
                    if key not in stretch_of:
                        a_s, b_s = ref[gap[0]:gap[1]], info[k]["oseq"][gap[2]:gap[3]]
                        stretch_of[key] = len(stretch_jobs)
                        stretch_jobs.append((a_s, b_s))
                    pending[name].append(stretch_of[key])
        if time.time() - t_fam > 5:
            print(f"[engine_eval] slow family {fam}: {time.time() - t_fam:.1f} s, {len(rows)} copies", file=sys.stderr)
    # Smith-Waterman vs shuffles for every stretch some engine put a difference pair in
    verdict = {}
    with open(f"{tmp}/sw_pairs.tsv", "w") as fh:
        for si, (a_s, b_s) in enumerate(stretch_jobs):
            if len(b_s) <= SHORT or not a_s:
                verdict[si] = "N"
                continue
            fh.write(f"{si}:obs\t{a_s}\t{b_s}\n")
            for t in range(NSHUF):
                sh = list(a_s)
                rng.shuffle(sh)
                fh.write(f"{si}:s{t}\t{''.join(sh)}\t{b_s}\n")
    with open(f"{tmp}/sw_pairs.tsv") as fi, open(f"{tmp}/sw_scores.tsv", "w") as fo:
        subprocess.run([SWSCORE], stdin=fi, stdout=fo, check=True)
    obs, null = {}, collections.defaultdict(list)
    for line in open(f"{tmp}/sw_scores.tsv"):
        key, sc = line.split()
        si, kind = key.split(":")
        if kind == "obs":
            obs[int(si)] = int(sc)
        else:
            null[int(si)].append(int(sc))
    for si in obs:
        verdict[si] = "H" if obs[si] > max(null[si]) else "N"
    res = {}
    for name in engines:
        c = pairs[name]
        for si in pending[name]:
            c[verdict[si]] += 1
        v, m = c["V"], c["M"]
        res[name] = {**c, "precision_V_over_VM": round(v / (v + m), 4) if v + m else None,
                     "difference_pairs": sum(c.values())}
    res["stretches_tested"] = len(obs)
    json.dump(res, open(os.path.join(outdir, "engine_eval.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], a[2], [tuple(x.split("=", 1)) for x in a[3:]])
