#!/usr/bin/env python3
"""Question 2 of docs/GOAL_multicopy_families_2026-11-19.md: score assemblies of a simulated tandem array (bench/famsim
scenario directory: truth.gtf, copies.tsv, reads.bam with the origin copy in each read name `copy|isoform|i`).

  array_score.py SCENARIO_DIR --arm NAME=GTF [...] [--out FILE]

Per arm, on the planted units (copies whose id starts with U) and the decoys (the precision control):
  expressed units with every truth chain emitted (exact intron chain, same contig and strand) / with at least one;
  silent units (expression 0) carrying any emitted multi-exon chain = false expression claims;
  multi-exon transcripts at no planted copy = off-target.
Read spread (no arm): primary alignments per unit against the reads simulated from it; the share placed on their own
unit; and the total-variation distance between the true and the primary-placed spread over units."""
import argparse, collections, csv, os, re, sys
import pysam

def gtf_chains(path):
    ex, at = collections.defaultdict(list), {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        ex[t].append((int(f[3]) - 1, int(f[4]))); at[t] = (f[0], f[6])
    out = []
    for t, e in ex.items():
        e.sort()
        out.append((at[t][0], at[t][1], e[0][0], e[-1][1], tuple((a[1], b[0]) for a, b in zip(e, e[1:])), t))
    return out

def load_copies(d):
    with open(os.path.join(d, "copies.tsv")) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scenario"); ap.add_argument("--arm", action="append", required=True); ap.add_argument("--out")
    a = ap.parse_args(argv)
    d = a.scenario
    copies = load_copies(d)
    key = lambda r: r.get("copy_id") or r.get("id") or r.get("copy")
    cid = {key(r): r for r in copies}
    truth = collections.defaultdict(set)                      # copy -> its truth chains
    span = {}
    for c, st, s, e, ch, t in gtf_chains(os.path.join(d, "truth.gtf")):
        cp = re.sub(r"[.|].*$", "", t) if t.split(".")[0] in cid or t.split("|")[0] in cid else None
        cp = next((k for k in cid if t == k or t.startswith(k + ".") or t.startswith(k + "|") or t.startswith(k + "_")), cp)
        if cp is None:
            continue
        if ch:
            truth[cp].add((c, st, ch))
        lo, hi = span.get(cp, (c, s, e))[1:], None
        span[cp] = (c, min(s, span.get(cp, (c, s, e))[1]), max(e, span.get(cp, (c, s, e))[2]))
    expr = {}
    for k, r in cid.items():
        v = r.get("expression") or r.get("reads") or "0"
        expr[k] = int(float(v))
    units = [k for k in cid if k.startswith("U")]
    decoys = [k for k in cid if not k.startswith("U") and k in truth]
    lines = [f"# {d}: units {len(units)} (expressed {sum(expr[u] > 0 for u in units)}), decoys {len(decoys)}"]
    # read spread from primaries
    true_n, placed, right = collections.Counter(), collections.Counter(), 0
    with pysam.AlignmentFile(os.path.join(d, "reads.bam")) as fh:
        for r in fh.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            src = r.query_name.split("|")[0]
            if src not in units:
                continue
            true_n[src] += 1
            at = next((u for u in units if span.get(u) and span[u][0] == r.reference_name and span[u][1] < r.reference_end and r.reference_start < span[u][2]), None)
            placed[at] += 1
            right += at == src
    tot = sum(true_n.values())
    tv = 0.5 * sum(abs(true_n[u] - placed[u]) for u in units) / tot if tot else 0
    lines.append(f"spread\tprimaries on their own unit {right}/{tot} = {right / tot if tot else 0:.3f}\ttotal variation true vs placed {tv:.3f}\t"
                 f"placed per unit {[placed[u] for u in units]}\ttrue {[true_n[u] for u in units]}")
    for spec in a.arm:
        name, g = spec.split("=", 1)
        if not os.path.exists(g):
            lines.append(f"{name}\tmissing {g}"); continue
        em = gtf_chains(g)
        keys = {(c, st, ch) for c, st, s, e, ch, t in em if ch}
        full = sum(1 for u in units if expr[u] > 0 and truth[u] and truth[u] <= keys)
        some = sum(1 for u in units if expr[u] > 0 and truth[u] & keys)
        silent = [u for u in units if expr[u] == 0]
        def at_unit(c, st, s, e):
            return next((u for u in cid if span.get(u) and span[u][0] == c and span[u][1] < e and s < span[u][2]), None)
        false_silent = len({at_unit(c, st, s, e) for c, st, s, e, ch, t in em if ch and at_unit(c, st, s, e) in silent})
        off = sum(1 for c, st, s, e, ch, t in em if ch and at_unit(c, st, s, e) is None)
        dec = sum(1 for k in decoys if truth[k] & keys)
        nexp = sum(expr[u] > 0 for u in units)
        lines.append(f"{name}\texpressed units all chains {full}/{nexp}, any chain {some}/{nexp}\tsilent units with a chain {false_silent}/{len(silent)}\t"
                     f"off-target multi-exon {off}\tdecoys found {dec}/{len(decoys)}\ttranscripts {len(em)}")
    out = "\n".join(lines) + "\n"
    if a.out:
        open(a.out, "w").write(out)
    print(out, end="")

if __name__ == "__main__":
    main()
