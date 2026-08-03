#!/usr/bin/env python3
"""Split RNA families by pairwise exon-sum COVERAGE, and score the purity/recall trade against Soto.

WHY. On the Soto DNA replication, the famCN splits that need WGS read depth turned out to be invisible to
sequence IDENTITY (median within-vs-between delta +0.0001) but obvious in alignment COVERAGE (+0.524, with
between-group coverage 0.000). That identity-vs-coverage contrast is a direct measurement and stands.

⚠ RETRACTED 2026-08-02: this docstring previously claimed "coverage alone recovered ~80% of famCN's benefit".
It did not. That figure was computed on a SHRINKING comparison set -- only genes both catalogs happened to
place in a family -- so genes the split dropped to singletons silently left the evaluation, and the anchor it
was measured against (ARI 0.720) is itself a retracted number. Re-scored on a FIXED gene set the value is
~31%, and on the DNA side coverage splitting actually LOWERS ARI (0.652 -> 0.625) while discarding 351 genes
instead of 96. Do not quote the 80%.

The pipeline's own E_r edge already requires coverage >= 0.50 to CREATE an edge. This is a different use:
a higher floor applied WITHIN an already-formed family to split it. So it can only ever split, never merge,
and recall can only fall (a member is lost solely when its family drops below 2 loci).

Scoring mirrors bench/soto/impurity_provenance.py: a "real" over-merge is one that survives assigning each
copy to its single best-overlap truth member, which strips the ~16-family artifact floor caused by Soto's
own truth intervals overlapping each other.
"""
import argparse, csv, os, subprocess, sys, tempfile
from collections import defaultdict


def load_fa(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = ln[1:].strip(), []
            else:
                buf.append(ln.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def truth(bed):
    per = defaultdict(list)
    n = 0
    with open(bed) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 4 or "|" not in f[3]:
                continue
            per[f[0]].append((f[3].split("|")[1], f[3].split("|")[0], int(f[1]), int(f[2])))
            n += 1
    return per, n


def score(assign, copies, per_truth, n_truth):
    """-> (member recall %, loose recall %, real over-merges, worst real fusion)

    ⚠ RECALL MUST ONLY COUNT MEMBERS STILL IN A FAMILY. The first version of this function credited a truth
    member as covered by ANY overlapping copy, including copies whose family the split had just dissolved
    into singletons (fid == ""). That makes "recall unchanged" true by construction and blinds the metric to
    the exact failure mode this split risks: RNA assembly incompleteness, where one copy assembles a full
    transcript and its sibling only a fragment, so the pair cannot clear the coverage floor and the family
    breaks apart. `loose` keeps the old (copy-level) number purely so the gap between the two is visible.
    """
    covered, loose, best = set(), set(), defaultdict(set)
    for (fid_orig, c, s, e), fid in assign.items():
        hits = [(min(e, te) - max(s, ts), tf, tn)
                for tf, tn, ts, te in per_truth.get(c, ()) if min(e, te) - max(s, ts) > 0]
        if not hits:
            continue
        for _, tf, tn in hits:
            loose.add((tf, tn))
            if fid:                      # only a copy still inside a >= 2-locus family counts
                covered.add((tf, tn))
        if fid:
            best[fid].add(max(hits)[1])
    impure = {k: v for k, v in best.items() if len(v) > 1}
    worst = max((len(v) for v in best.values()), default=0)
    return (100.0 * len(covered) / max(n_truth, 1), 100.0 * len(loose) / max(n_truth, 1),
            len(impure), worst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--copies", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bed", default="bench/soto/80_fams.chr.bed")
    ap.add_argument("--thresholds", default="0.0,0.5,0.7,0.9,0.95")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()

    seqs = load_fa(a.fasta)
    # header: GWFAM0|0|chr1:99596-106051|+|nexon=1  -> key on (family, chrom, start, end)
    key_of = {}
    for h in seqs:
        p = h.split("|")
        if len(p) < 3:
            continue
        c, rng = p[2].rsplit(":", 1)
        s, e = rng.split("-")
        key_of[h] = (p[0], c, int(s), int(e))

    copies = []
    with open(a.copies) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("chrom"):
                copies.append((r["family_id"], r["chrom"], int(r["start"]), int(r["end"])))
    per_fam = defaultdict(list)
    for h, k in key_of.items():
        per_fam[k[0]].append(h)

    per_truth, n_truth = truth(a.bed)

    # pairwise best coverage inside each family (one minimap2 per family)
    pair = defaultdict(dict)
    fams = [f for f, hs in per_fam.items() if len(hs) >= 2]
    for i, f in enumerate(fams, 1):
        hs = per_fam[f]
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
            for j, h in enumerate(hs):
                fh.write(f">{j}\n{seqs[h]}\n")
            p = fh.name
        out = subprocess.run([a.minimap2, "-c", "-X", "--no-long-join", "-x", "asm20",
                              "-t", str(a.threads), p, p], capture_output=True, text=True).stdout
        os.unlink(p)
        for ln in out.splitlines():
            fl = ln.split("\t")
            if len(fl) < 12:
                continue
            q, t = fl[0], fl[5]
            if q == t:
                continue
            cov = (int(fl[3]) - int(fl[2])) / min(int(fl[1]), int(fl[6]))
            k = tuple(sorted((int(q), int(t))))
            pair[f][k] = max(pair[f].get(k, 0.0), cov)
        if i % 100 == 0:
            print(f"  aligned {i}/{len(fams)} families", file=sys.stderr)

    print(f"{'coverage floor':<16}{'families':>9}{'recall %':>10}{'(loose)':>9}{'over-merges':>13}{'worst':>7}")
    print("-" * 66)
    for T in [float(x) for x in a.thresholds.split(",")]:
        assign, fi = {}, 0
        for f, hs in per_fam.items():
            n = len(hs)
            par = list(range(n))

            def find(x):
                while par[x] != x:
                    par[x] = par[par[x]]
                    x = par[x]
                return x

            for (u, v), cov in pair.get(f, {}).items():
                if cov >= T:
                    ru, rv = find(u), find(v)
                    if ru != rv:
                        par[ru] = rv
            grp = defaultdict(list)
            for j in range(n):
                grp[find(j)].append(j)
            for members in grp.values():
                # a family still needs >= 2 distinct loci; smaller groups are dropped (recall cost)
                fid = f"CS{fi}" if len(members) >= 2 else ""
                if fid:
                    fi += 1
                for j in members:
                    assign[key_of[hs[j]]] = fid
        rec, loose, over, worst = score(assign, copies, per_truth, n_truth)
        tag = "  (= no split, baseline)" if T == 0.0 else ""
        print(f"{T:<16.2f}{fi:>9}{rec:>10.1f}{loose:>9.1f}{over:>13}{worst:>7}{tag}")


if __name__ == "__main__":
    main()
