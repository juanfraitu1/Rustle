"""PREREG r1796 (docs/PREREG_famwise_origin_2026-11-18.md), row a: statements about a molecule and transcript claims
against the origin of the reads.

  origin_score.py --truth READS.tsv [--tie-split T.tsv] [--families READ_FAMILY.tsv] --bam SEL.bam
                  [--arm NAME=GTF ...] [--tx FW.tx.tsv --tx-arm NAME] --out DIR

Writes DIR/molecule.tsv, molecule_family.tsv (with --families), one_score_only.tsv, claims.tsv, claims_tx.tsv (the judged
claims of the --tx-arm arm, one row each), range.tsv, bars.tsv.
Every band is computed on ms and on AS; the bars are on ms. A statement holding a record on a place the truth never
judged has no verdict: counted, never scored. Shares are over n - no_verdict."""
import argparse, collections, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from gtf_names import chain, opn

LIMITS = ("# limits: annotation families only; clear-lead reads sampled (30,000); 12 loci per read; kept placements at AS >= 0.8 x best; "
          "HG002 is a 25% read sample on another individual's reference; `outside` reads give the upper bound on failures of assumption A\n")
WINDOWS = (("primary", None), ("0", 0), ("5", 5), ("10", 10), ("20", 20))
KINDS = ("one", "several", "one_place", "tie_unsplit", "outside", "unjudged")
VERDICTS = ("TRUE", "SHARED", "PHANTOM", "UNJUDGED")

def load_claims(gtf):
    """-> {(contig, strand, chain): [transcript ids]}, multi-exon transcripts"""
    ex, at = collections.defaultdict(list), {}
    for l in opn(gtf):
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            ex[t].append((int(f[3]), int(f[4]))); at[t] = (f[0], f[6])
    out = collections.defaultdict(list)
    for t, e in ex.items():
        if len(e) > 1:
            out[(*at[t], chain(e))].append(t)
    return out

def verdict(stmt, t):
    """-> (equals, sound), or None = no verdict (a record on a place the truth never judged)"""
    idx = {R.locus_index(r, t) for r in stmt}
    return None if None in idx else (idx == t.origin, bool(idx & t.origin))

def write(path, header, rows, notes=()):
    with open(path, "w") as fh:
        fh.write(LIMITS + "".join(f"# {n}\n" for n in notes) + "\t".join(header) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)

def share(k, n):
    lo, hi = R.wilson(k, n)
    return [k, n, f"{k / n:.4f}" if n else "NA", f"{lo:.4f}", f"{hi:.4f}"]

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True); ap.add_argument("--tie-split"); ap.add_argument("--families")
    ap.add_argument("--bam", required=True); ap.add_argument("--arm", action="append", default=[])
    ap.add_argument("--tx"); ap.add_argument("--tx-arm"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    truth, recs = R.load_truth(a.truth, a.tie_split), R.load_records(a.bam)
    both = sorted(set(truth) & set(recs))
    notes = [f"reads: truth {len(truth)}, BAM {len(recs)}, both {len(both)}"]
    fam = {}
    if a.families:
        with open(a.families) as fh:
            fam = dict(l.rstrip("\n").split("\t")[:2] for l in fh if l.strip())

    # molecules
    mol, famt, one_only = collections.Counter(), collections.Counter(), []
    b1, b1fam = collections.Counter(), collections.defaultdict(collections.Counter)
    for q in both:
        t, rs = truth[q], recs[q]
        near = {s: R.near_tied(rs, s) for s in ("ms", "AS")}
        if near["ms"] != near["AS"]:
            one_only.append((q, "ms" if near["ms"] else "AS", R.lead(rs, "ms"), R.lead(rs, "AS")))
        v = {}
        for s in ("ms", "AS"):
            for w, win in WINDOWS:
                v[s, w] = verdict(R.statement(rs, s, win), t)
                c = [1, 0, 0, 1] if v[s, w] is None else [1, *map(int, v[s, w]), 0]
                keys = [("m", s, w, "yes" if near[s] else "no", t.kind)]
                if s == "ms" and w in ("primary", "0") and q in fam:
                    keys.append(("f", fam[q], w, "yes" if near[s] else "no", t.kind))
                for k in keys:
                    for i, x in enumerate(c):
                        (mol if k[0] == "m" else famt)[(*k[1:], i)] += x
        if near["ms"] and t.kind in R.JUDGED:      # B1: the same molecules in both shares
            p, e = v["ms", "primary"], v["ms", "0"]
            if p is None or e is None:
                b1["no_verdict"] += 1
            else:
                b1["n"] += 1; b1["primary"] += p[0]; b1["exact"] += e[0]
                if q in fam:
                    f = b1fam[fam[q]]; f["n"] += 1; f["primary"] += p[0]; f["exact"] += e[0]
    cells = sorted({k[:-1] for k in mol}, key=lambda k: (k[0], [w for w, _ in WINDOWS].index(k[1]), k[2], KINDS.index(k[3])))
    write(f"{a.out}/molecule.tsv", ("score", "window", "near_tied", "kind", "n", "equals", "sound", "no_verdict"),
          [(*k, *(mol[(*k, i)] for i in range(4))) for k in cells], notes)
    if fam:
        cells = sorted({k[:-1] for k in famt})
        write(f"{a.out}/molecule_family.tsv", ("family", "window", "near_tied", "kind", "n", "equals", "sound", "no_verdict"),
              [(*k, *(famt[(*k, i)] for i in range(4))) for k in cells], notes)
    write(f"{a.out}/one_score_only.tsv", ("read", "near_tied_on", "lead_ms", "lead_AS"), one_only, notes)

    # claims: one verdict per (contig, strand, chain), whatever arm emits it
    arms = {n: load_claims(g) for n, g in (x.split("=", 1) for x in a.arm)}
    keys = set().union(*arms.values()) if arms else set()
    sup = collections.defaultdict(list)            # key -> [(origin set, locus index)] of judged supporters
    for q in both:
        t = truth[q]
        if t.kind in R.JUDGED:
            for r in R.statement(recs[q], "ms", 10):
                k, i = (r.chrom, r.strand, r.chain), R.locus_index(r, t)
                if r.chain and k in keys and i is not None:
                    sup[k].append((t.origin, i))
    def claim(k):
        s = sup.get(k, ())
        return "UNJUDGED" if not s else "TRUE" if any(o == {i} for o, i in s) else "SHARED" if any(i in o for o, i in s) else "PHANTOM"
    status = {}
    if a.tx:
        with open(a.tx) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for l in fh:
                d = dict(zip(head, l.rstrip("\n").split("\t")))
                if d["band"] == "wide":
                    status.update((t, d["status"]) for t in d["transcripts"].split(","))
    cl, detail = collections.Counter(), []
    for n, ks in arms.items():
        for k, tids in ks.items():
            st = min(status.get(t, "not_in_tx") for t in tids) if n == a.tx_arm else "."
            cl[n, st, "claims"] += 1; cl[n, st, claim(k)] += 1
            if n == a.tx_arm and k in sup:
                detail.append((",".join(tids), k[0], k[2][0][0], k[2][-1][1], k[1], st, claim(k), len(sup[k])))
    write(f"{a.out}/claims.tsv", ("arm", "status", "claims", *VERDICTS),
          [(n, st, *(cl[n, st, x] for x in ("claims", *VERDICTS))) for n, st in sorted({k[:2] for k in cl})], notes)

    if a.tx_arm:
        write(f"{a.out}/claims_tx.tsv", ("transcripts", "contig", "first_intron_start", "last_intron_end", "strand", "status", "verdict", "supporters"), sorted(detail, key=lambda d: d[1:4]), notes)

    # range: on molecules with one origin, per node of the --tx-arm arm
    nodes, rng = set(arms.get(a.tx_arm, ())), []
    for band, win in (("exact", 0), ("wide", 10)):
        lo, hi, tru, holds, n = collections.Counter(), collections.Counter(), collections.Counter(), 0, 0
        for q in both:
            t = truth[q]
            if t.kind != "one":
                continue
            n += 1
            stmt = R.statement(recs[q], "ms", win)
            at = {(r.chrom, r.strand, r.chain) for r in stmt}
            holds += any(R.locus_index(r, t) in t.origin for r in stmt)
            for k in at & nodes:
                hi[k] += 1; lo[k] += len(at) == 1
            home = [r for r in recs[q] if R.locus_index(r, t) in t.origin]
            if home:
                r = max(home, key=lambda r: r.ms)
                tru[(r.chrom, r.strand, r.chain)] += (r.chrom, r.strand, r.chain) in nodes
        seen = set(hi) | {k for k, v in tru.items() if v}
        rng.append((band, len(seen), sum(lo[k] <= tru[k] <= hi[k] for k in seen), holds, n))
    write(f"{a.out}/range.tsv", ("band", "nodes", "held", "a_holds", "a_n"), rng, notes)

    # bars
    bars = [("B1 exact band equals", *share(b1["exact"], b1["n"]), ""), ("B1 primary equals", *share(b1["primary"], b1["n"]), ""),
            ("B1 molecules without a verdict", b1["no_verdict"], "", "", "", "", "")]
    big = [f for f in b1fam.values() if f["n"] >= 20]
    won = sum(f["exact"] > f["primary"] for f in big)
    bars.append(("B1 families exact > primary", won, len(big), "", "", "", "families with >= 20 judged near-tied molecules"))
    disjoint = R.wilson(b1["exact"], b1["n"])[0] > R.wilson(b1["primary"], b1["n"])[1]
    bars.append(("B1 verdict", "", "", "", "", "", "undecidable" if not big or not b1["n"] else "PASS" if disjoint and 2 * won > len(big) else "FAIL"))
    judged = lambda st: sum(cl[a.tx_arm, st, x] for x in VERDICTS[:3])
    bad = lambda st: cl[a.tx_arm, st, "SHARED"] + cl[a.tx_arm, st, "PHANTOM"]
    bars += [(f"B2i {st} PHANTOM+SHARED", *share(bad(st), judged(st)), "") for st in ("tied_only", "certified")]
    worse = R.wilson(bad("tied_only"), judged("tied_only"))[0] > R.wilson(bad("certified"), judged("certified"))[1]
    bars.append(("B2i verdict", "", "", "", "", "", "undecidable" if judged("tied_only") < 30 else "PASS" if worse else "FAIL"))
    true_all = sum(v for (n, st, x), v in cl.items() if n == a.tx_arm and x == "TRUE")
    true_top = cl[a.tx_arm, "certified", "TRUE"] + cl[a.tx_arm, "best_only", "TRUE"]
    bars.append(("B2ii certified+best_only of TRUE", *share(true_top, true_all), ""))
    bars.append(("B2ii verdict", "", "", "", "", "", "undecidable" if not true_all else "PASS" if true_top >= 0.9 * true_all else "FAIL"))
    bars.append(("check certified PHANTOM", *share(cl[a.tx_arm, "certified", "PHANTOM"], judged("certified")), "not a bar"))
    write(f"{a.out}/bars.tsv", ("quantity", "k", "n", "share", "lo", "hi", "verdict"), bars, notes + ["bands on ms; intervals Wilson 95%"])
    print(notes[0]); print("".join(f"{b[0]}\t{b[1]}\t{b[2]}\t{b[3]}\t{b[6]}\n" for b in bars), end="")

if __name__ == "__main__":
    main()
