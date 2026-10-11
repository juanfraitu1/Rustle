#!/usr/bin/env python3
"""How many seeds does the DNA (guided) mode need to find every LRPAP1 copy? (docs/PREREG_lrpap1_seed_reach_2026-10-11.md)

    seed_reach.py <browser.json> <out.json>        aligns every copy as a seed with the guided pipeline's two finders, then scores seed sets

Seed i finds copy j when a passing hit of i (guided_pipeline.transcript_hits / gene_body_chains, unchanged) overlaps j's span. One pass = the union
of the seeds' finds; iteration = reachability in the finds graph. Pure functions are tested in test_seed_reach.py."""
import itertools
import json
import os
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

GENOME_FA = "/mnt/linuxdisk/tmp/rustle_figures/liftoff/gorilla/genome.fa"
GENOME_MMI = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO.splice.mmi"


def finds(hits, copies):
    """hits = [(seed, chrom, start, end)], copies = {id: (chrom, lo, hi)} -> ({seed: set of copies its hits overlap}, {seed: distinct loci outside
    every copy})"""
    f = {c: set() for c in copies}
    outside = {c: [] for c in copies}
    for q, chrom, s, e in hits:
        hit = {c for c, (ch, lo, hi) in copies.items() if ch == chrom and s < hi and lo < e}
        if hit:
            f[q] |= hit
        else:
            outside[q].append((chrom, s, e))
    extra = {}
    for q, iv in outside.items():
        n, last = 0, None
        for chrom, s, e in sorted(iv):
            if last is None or chrom != last[0] or s >= last[2]:
                n += 1
                last = (chrom, s, e)
            else:
                last = (chrom, last[1], max(last[2], e))
        extra[q] = n
    return f, extra


def one_pass(seeds, F):
    return set(seeds).union(*(F[s] for s in seeds)) if seeds else set()


def closure(seeds, F):
    seen, todo = set(seeds), list(seeds)
    while todo:
        for j in F[todo.pop()]:
            if j not in seen:
                seen.add(j)
                todo.append(j)
    return seen


def curve(F, ids, mode):
    out = []
    for k in range(1, len(ids) + 1):
        n = [len(mode(set(S), F)) for S in itertools.combinations(ids, k)]
        out.append(dict(k=k, n_sets=len(n), min=min(n), median=statistics.median(n), max=max(n),
                        all_found=sum(x == len(ids) for x in n) / len(n)))
    return out


def min_seeds(F, ids, mode):
    for k in range(1, len(ids) + 1):
        if any(len(mode(set(S), F)) == len(ids) for S in itertools.combinations(ids, k)):
            return k
    return None


def source_components(F):
    """strongly connected components of the finds graph that no other component reaches (each needs one seed under iteration)"""
    reach = {i: closure({i}, F) for i in F}
    comps, done = [], set()
    for i in sorted(F):
        if i not in done:
            c = sorted(j for j in F if j in reach[i] and i in reach[j])
            done |= set(c)
            comps.append(c)
    return [c for c in comps if not any(c[0] in reach[j] for j in F if j not in c)]


def sweep(hits, copies, ids, thresholds):
    """hits = [(seed, chrom, start, end, identity)] -> per finder identity threshold: seeds needed (one pass, iteration) and the source components"""
    rows = []
    for t in thresholds:
        F, _ = finds([h[:4] for h in hits if h[4] >= t], copies)
        rows.append(dict(threshold=t, min_seeds_one_pass=min_seeds(F, ids, one_pass), min_seeds_iterate=min_seeds(F, ids, closure),
                         components=source_components(F)))
    return rows


def _rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def main(browser, out):
    import pysam
    import guided_pipeline as G
    copies = json.load(open(browser))
    fa = pysam.FastaFile(GENOME_FA)
    work = os.path.join(os.path.dirname(os.path.abspath(out)), "seed_reach_work")
    os.makedirs(work, exist_ok=True)
    with open(f"{work}/units.fa", "w") as fu, open(f"{work}/bodies.fa", "w") as fb:
        for c in copies:
            ex = sorted(tuple(x) for x in c["model"])
            unit = "".join(fa.fetch(c["acc"], s, e) for s, e in ex).upper()
            body = fa.fetch(c["acc"], ex[0][0], ex[-1][1]).upper()
            if c["strand"] == "-":
                unit, body = _rc(unit), _rc(body)
            fu.write(f">{c['cid']}\n{unit}\n")
            fb.write(f">{c['cid']}\n{body}\n")
    for name, preset in (("units", "splice"), ("bodies", "asm20")):
        if not os.path.exists(f"{work}/{name}.paf.done"):
            with open(f"{work}/{name}.paf", "w") as fh:
                subprocess.run(["minimap2", "-c", "-x", preset, "-N", "100", "-p", "0.1", "-t", "4", GENOME_MMI, f"{work}/{name}.fa"],
                               stdout=fh, stderr=subprocess.DEVNULL, check=True)
            open(f"{work}/{name}.paf.done", "w").write("ok")
    tx = G.transcript_hits(f"{work}/units.paf")
    ch = G.gene_body_chains(f"{work}/bodies.paf")
    span = {c["cid"]: (c["acc"], c["lo"], c["hi"]) for c in copies}
    hits = [(h["q"], h["chrom"], h["s"], h["e"]) for h in tx] + [(h["q"], h["chrom"], h["xs"], h["xe"]) for h in ch]
    F, extra = finds(hits, span)
    ftx, _ = finds([(h["q"], h["chrom"], h["s"], h["e"]) for h in tx], span)
    fch, _ = finds([(h["q"], h["chrom"], h["xs"], h["xe"]) for h in ch], span)
    ident = {}
    for h in ch:                                            # best gene-body identity of seed i on copy j
        for j, (chrom, lo, hi) in span.items():
            if h["chrom"] == chrom and h["xs"] < hi and lo < h["xe"]:
                v = h["nm"] / sum(r["bl"] for r in h["recs"])
                ident[(h["q"], j)] = max(ident.get((h["q"], j), 0), v)
    ids = [c["cid"] for c in copies]
    res = dict(ids=ids, names={c["cid"]: c["name"] for c in copies}, cls={c["cid"]: c["cls"] for c in copies},
               finds={i: sorted(F[i]) for i in ids}, finds_transcript={i: sorted(ftx[i]) for i in ids},
               finds_body={i: sorted(fch[i]) for i in ids}, extra_loci=extra,
               identity={f"{i}|{j}": round(v, 4) for (i, j), v in ident.items()},
               one_pass=curve(F, ids, one_pass), iterate=curve(F, ids, closure),
               min_seeds_one_pass=min_seeds(F, ids, one_pass), min_seeds_iterate=min_seeds(F, ids, closure),
               source_components=source_components(F),
               single_seed={i: dict(one_pass=len(one_pass({i}, F)), iterate=len(closure({i}, F))) for i in ids},
               n_transcript_hits=len(tx), n_body_chains=len(ch))
    h5 = [(h["q"], h["chrom"], h["s"], h["e"], h["nm"] / h["bl"]) for h in tx] + \
         [(h["q"], h["chrom"], h["xs"], h["xe"], h["nm"] / sum(r["bl"] for r in h["recs"])) for h in ch]
    res["sweep"] = sweep(h5, span, ids, [round(0.80 + 0.0025 * i, 4) for i in range(81)])
    res["curves"] = {}
    for t in (0.80, 0.98, 0.99, 0.995, 1.0):
        Ft, _ = finds([h[:4] for h in h5 if h[4] >= t], span)
        res["curves"][f"{t:.3f}"] = dict(one_pass=curve(Ft, ids, one_pass), iterate=curve(Ft, ids, closure))
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("min_seeds_one_pass", "min_seeds_iterate", "source_components", "single_seed", "extra_loci")}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
