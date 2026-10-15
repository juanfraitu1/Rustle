#!/usr/bin/env python3
"""Post-MCL strongest-edge reassignment on a shipped catalog (prereg docs/PREREG_post_mcl_reassign_2026-10-15.md).

    run_reassign.py dump  SAMPLE --out DIR [--bin MCL_FAMILIES]     re-cluster from the cached all-vs-all with
                                                                      --dump-graph; the clusters must equal the shipped
                                                                      ones byte for byte (HEAVY for human_A119b)
    run_reassign.py apply SAMPLE --out DIR                           reassign; DIR/SAMPLE/SAMPLE.fam.copies.tsv is the
                                                                      shipped copy table with the new family labels
    run_reassign.py score SAMPLE --out DIR [--null 1000] [--seed 7]  F1, the pairs it joins and splits, the decoy test,
                                                                      Soto gains and losses (human), the guards
SAMPLE = human_A119b | human_testis | gorilla_KB3781.
"""
import argparse
import collections
import csv
import filecmp
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "bench", "liftoff_pairs"))
sys.path.insert(0, os.path.join(REPO, "bench", "soto_subset"))
import reassign_lib as R  # noqa: E402
import run_pairs as RP  # noqa: E402

LF, P = RP.LF, RP.P
RUNS = "/mnt/linuxdisk/tmp/rustle_figures/runs"
LIFTOFF = "/mnt/linuxdisk/tmp/rustle_figures/liftoff"
SHIPPED_HIER = "/mnt/linuxdisk/tmp/o1_close/hier_shipped"
BIN = "/mnt/linuxdisk/home/juanfraitu/rustle_target/release/mcl_families"
SAMPLES = {"human_A119b": "human", "human_testis": "human", "gorilla_KB3781": "gorilla"}
GROWTH, COVERED = 0.10, 0.02


def paths(runs, out, s):
    d = os.path.join(out, s)
    return dict(run=os.path.join(runs, s), dir=d, dump=os.path.join(d, "graph.tsv"),
                shipped_clusters=os.path.join(runs, s, f"{s}.fam.clusters.tsv"),
                re_clusters=os.path.join(d, f"{s}.fam.clusters.tsv"),
                copies=os.path.join(runs, s, f"{s}.fam.copies.tsv"), gff=os.path.join(runs, s, f"{s}.fam.loci.gff3"),
                fold=os.path.join(runs, s, f"{s}.fam.loci.tsv"), new_copies=os.path.join(d, "relabelled", s,
                                                                                         f"{s}.fam.copies.tsv"),
                moves=os.path.join(d, "moves.tsv"), apply=os.path.join(d, "apply.json"),
                score=os.path.join(d, "score.json"))


def stage_dump(a, p):
    os.makedirs(p["dir"], exist_ok=True)
    subprocess.run(["bash", os.path.join(REPO, "bench", "o1_hier", "recluster.sh"), a.bin, a.sample, p["dir"],
                    p["dump"]], check=True)
    same = filecmp.cmp(p["re_clusters"], p["shipped_clusters"], shallow=False)
    print(json.dumps(dict(sample=a.sample, reproduces_shipped_clusters=same)))
    if not same:
        raise SystemExit(f"{p['re_clusters']} differs from {p['shipped_clusters']}: the dump is not the shipped graph")


def gene_nodes(gff):
    """GFF gene Name -> its node `chrom:start-end` (1-based, as the all-vs-all names it)."""
    out = {}
    with open(gff) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                out[re.search(r"Name=([^;\n]+)", f[8]).group(1)] = f"{f[0]}:{f[3]}-{f[4]}"
    return out


def copy_nodes(copies, gff, fold, shipped):
    """Per copy-table row: its member node (gene -> span node -> folded representative). Every row must land on a node
    of its own shipped family, or the mapping is wrong and nothing downstream is valid."""
    names = gene_nodes(gff)
    with open(fold) as fh:
        rep = {r["annotation"]: r["representative"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(copies) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    nodes, bad = [], []
    for r in rows:
        n = names.get(r["gene_id"])
        n = rep.get(n, n)
        if n is None or shipped.get(n) != r["family_id"]:
            bad.append((r["family_id"], r["gene_id"], n))
        nodes.append(n)
    if bad:
        raise SystemExit(f"{len(bad)} copy rows do not map to a node of their shipped family, e.g. {bad[:3]}")
    return rows, nodes


def write_tsv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def stage_apply(a, p):
    shipped = R.read_clusters(p["shipped_clusters"])
    new, log = R.reassign(shipped, R.read_dump(p["dump"]))
    rows, nodes = copy_nodes(p["copies"], p["gff"], p["fold"], shipped)
    with open(p["copies"]) as fh:
        header = next(csv.reader(fh, delimiter="\t"))
    for r, n in zip(rows, nodes):
        r["family_id"] = new[n]
    write_tsv(p["new_copies"], header, rows)
    write_tsv(p["moves"], ["pass", "node", "from", "to", "w_in", "w_out"],
              [dict(zip(["pass", "node", "from", "to", "w_in", "w_out"], m)) for m in log["moves"]])
    size = lambda fam: sorted(collections.Counter(fam.values()).values())  # noqa: E731
    res = dict(sample=a.sample, members=len(shipped), moves=len(log["moves"]),
               moved_nodes=len({m[1] for m in log["moves"]}), passes=log["passes"], converged=log["converged"],
               lost=log["lost"], families_before=len(set(shipped.values())), families_after=len(set(new.values())),
               largest_before=size(shipped)[-1], largest_after=size(new)[-1],
               copies_relabelled=sum(1 for r, n in zip(rows, nodes) if shipped[n] != new[n]))
    write_json(p["apply"], res)
    print(json.dumps(res, sort_keys=True))


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True)
        fh.write("\n")


def soto_labels(sample, runs):
    """{member node: set(Soto family)} over the de novo universe of a human catalog (soto_subset's own mapping)."""
    import universes as U
    import cn_tables as CN
    d = os.path.join(SHIPPED_HIER, sample)
    uni = U.denovo(CN.TRUTH_V2, os.path.join(runs, sample, f"{sample}.fam.loci.gff3"),
                   os.path.join(d, "ours.hierarchy.members.tsv"), os.path.join(d, "ours.loci.tsv"))
    out = {}
    for lab, ms in uni.soto_sets.items():
        for m in ms:
            out.setdefault(m, set()).add(lab)
    return out


def soto_pairs(fam, labels):
    """Pairwise agreement with Soto over labelled members in a family: precision, recall."""
    import itertools
    ms = sorted(m for m in labels if m in fam)
    same_fam = same_lab = both = 0
    for x, y in itertools.combinations(ms, 2):
        f, t = fam[x] == fam[y], bool(labels[x] & labels[y])
        same_fam += f
        same_lab += t
        both += f and t
    return dict(members=len(ms), precision=both / same_fam if same_fam else None,
                recall=both / same_lab if same_lab else None, tp=both, pred=same_fam, truth=same_lab)


def decoy(sample, species, fam, n_null, seed):
    """Register 1284's decoy test, primary row only (genome, sequence_ID >= 0.95)."""
    loci = LF.read_loci(Path(LIFTOFF) / species / "liftoff_loci.tsv")
    sup = LF.read_support(Path(LIFTOFF) / "support" / f"{sample}.tsv")
    pairs = LF.copy_pairs(loci, 0.95, sup=sup)
    pool = P.decoy_pool(loci, sup, None)
    rows = {LF.locus_key(r): r for r in pool}
    rows.update((LF.locus_key(r), r) for pair in pairs for r in pair)
    fams = P.cover_families(list(rows.values()), fam)
    t = P.run_test(pairs, P.Pool(pool, fams), fams, n_null, seed)
    return dict(n=t["n"], k=t["k"], rate=t["rate"], null_mean=t["null_mean"], p=t["p"], passes=bool(t["p"] < 0.001))


def stage_score(a, p):
    species = SAMPLES[a.sample]
    loci = LF.read_loci(Path(LIFTOFF) / species / "liftoff_loci.tsv")
    before = LF.catalog_loci(Path(p["copies"]))
    after = LF.catalog_loci(Path(p["new_copies"]))
    f1b, f1a = RP.f1(loci, before, species), RP.f1(loci, after, species)
    prs = LF.copy_pairs(loci, 0.95)
    changes = R.pair_changes(LF.pair_families(prs, before), LF.pair_families(prs, after))
    ap = json.load(open(p["apply"]))
    nb, na = f1b["genome"]["0.95"][0], f1a["genome"]["0.95"][0]
    res = dict(sample=a.sample, f1_before=f1b, f1_after=f1a, pair_changes=changes, apply=ap,
               decoy_before=decoy(a.sample, species, before, a.null, a.seed),
               decoy_after=decoy(a.sample, species, after, a.null, a.seed),
               guards=dict(no_member_lost=not ap["lost"], covered_within=abs(na - nb) <= COVERED * nb,
                           largest_growth_ok=ap["largest_after"] <= (1 + GROWTH) * ap["largest_before"],
                           converged=ap["converged"]))
    if species == "human":
        labels = soto_labels(a.sample, a.runs)
        old = R.read_clusters(p["shipped_clusters"])
        new = dict(old)
        with open(p["moves"]) as fh:
            moved = []
            for r in csv.DictReader(fh, delimiter="\t"):
                new[r["node"]] = r["to"]
                moved.append(r["node"])
        moved = sorted(set(moved))
        res["soto_moves"] = R.label_moves(moved, old, new, labels)
        res["soto_pairs_before"], res["soto_pairs_after"] = soto_pairs(old, labels), soto_pairs(new, labels)
    write_json(p["score"], res)
    print(json.dumps({k: res[k] for k in ("f1_before", "f1_after", "pair_changes", "guards")}, sort_keys=True))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["dump", "apply", "score"])
    ap.add_argument("sample", choices=sorted(SAMPLES))
    ap.add_argument("--out", required=True)
    ap.add_argument("--bin", default=BIN)
    ap.add_argument("--runs", default=RUNS)
    ap.add_argument("--null", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args(argv)
    p = paths(a.runs, a.out, a.sample)
    {"dump": stage_dump, "apply": stage_apply, "score": stage_score}[a.stage](a, p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
