#!/usr/bin/env python3
"""Question 1 of docs/GOAL_multicopy_families_2026-11-19.md: a resolution call per family transcript from the reads alone
(no truth is read here).

  resolution_call.py FW.gtf SEL.bam OUT_PREFIX [--tied-share 0.5] [--window 0]

A family transcript of the `famwise` layer (transcripts sharing `fw_tx`) whose transcripts sit at >= 2 places is a
candidate copy set. Its supporters: molecules with a record carrying exactly one of its intron chains, within 10 ms of
the molecule's best (or the primary). A supporter is TIED when its records within --window ms of its best (0 = exact
ties) reach >= 2 of the family transcript's places; PLACED otherwise.
  tied share = tied / supporters;  call `copy_set` when tied share >= --tied-share, else `separable`.
Writes OUT_PREFIX.calls.tsv and OUT_PREFIX.gtf = FW.gtf with `fw_call` on the transcripts of every candidate and
`fw_tx` dropped under `separable` (its copies are then reported one by one)."""
import argparse, collections, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R
from copyset_score import load_gtf, n_places

def places(spans):
    out = []
    for c, st, s, e in sorted(spans):
        if out and out[-1][0] == c and out[-1][1] == st and s < out[-1][3]:
            out[-1][3] = max(out[-1][3], e)
        else:
            out.append([c, st, s, e])
    return out

def place_of(rec, pl):
    for i, (c, st, s, e) in enumerate(pl):
        if c == rec.chrom and st == rec.strand and s < rec.end and rec.start < e:
            return i
    return None

def call_claims(claims, recs, share, window):
    """-> [(claim, supporters, tied, tied share, call)] for claims at >= 2 places"""
    key_claim = collections.defaultdict(list)
    cand = [c for c in claims if n_places(c["spans"]) >= 2]
    for j, c in enumerate(cand):
        for k in c["keys"]:
            key_claim[k].append(j)
    sup = collections.defaultdict(set)
    for q, rs in recs.items():
        for r in R.statement(rs, "ms", 10):
            for j in key_claim.get((r.chrom, r.strand, r.chain), ()):
                sup[j].add(q)
    out = []
    for j, c in enumerate(cand):
        pl = places(c["spans"])
        tied = 0
        for q in sup[j]:
            rs = recs[q]; best = max(r.ms for r in rs)
            hit = {place_of(r, pl) for r in rs if best - r.ms <= window} - {None}
            tied += len(hit) >= 2
        n = len(sup[j])
        ts = tied / n if n else 0.0
        out.append((c, n, tied, ts, "copy_set" if n and ts >= share else "separable"))
    return out

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gtf"); ap.add_argument("bam"); ap.add_argument("out")
    ap.add_argument("--tied-share", type=float, default=0.5); ap.add_argument("--window", type=float, default=0)
    a = ap.parse_args(argv)
    claims, recs = load_gtf(a.gtf, fw=True), R.load_records(a.bam)
    res = call_claims(claims, recs, a.tied_share, a.window)
    by_key = {}
    with open(a.out + ".calls.tsv", "w") as fh:
        fh.write("places\tfirst_place\tsupporters\ttied\ttied_share\tcall\n")
        for c, n, tied, ts, cl in res:
            pl = places(c["spans"])
            fh.write(f"{len(pl)}\t{pl[0][0]}:{pl[0][2]}-{pl[0][3]}:{pl[0][1]}\t{n}\t{tied}\t{ts:.3f}\t{cl}\n")
            for k in c["keys"]:
                by_key[k] = cl
    # rewrite: the call of a transcript is the call of its family transcript, found through its fw_tx group
    grp_call, lines = {}, open(a.gtf).read().splitlines(True)
    from gtf_names import chain
    ex, fw = collections.defaultdict(list), {}
    for l in lines:
        f = l.split("\t")
        if len(f) < 9:
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        if not t:
            continue
        if f[2] == "exon":
            ex[t.group(1)].append((f[0], f[6], int(f[3]), int(f[4])))
        m = re.search(r'fw_tx "([^"]+)"', f[8])
        if m and f[2] == "transcript":
            fw[t.group(1)] = m.group(1)
    for t, e in ex.items():
        if len(e) > 1 and t in fw:
            e.sort(key=lambda v: v[2])
            k = (e[0][0], e[0][1], chain([(s, x) for _, _, s, x in e]))
            if k in by_key:
                grp_call[fw[t]] = by_key[k]
    with open(a.out + ".gtf", "w") as fo:
        for l in lines:
            m = re.search(r' fw_tx "([^"]+)";', l)
            if m and m.group(1) in grp_call:
                c = grp_call[m.group(1)]
                l = l.replace(m.group(0), ("" if c == "separable" else m.group(0)) + f' fw_call "{c}";')
            fo.write(l)
    print(f"resolution_call: {len(res)} family transcripts at >= 2 places: copy_set {sum(r[4] == 'copy_set' for r in res)}, "
          f"separable {sum(r[4] == 'separable' for r in res)} (tied share >= {a.tied_share}, window {a.window})")

if __name__ == "__main__":
    main()
