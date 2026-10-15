#!/usr/bin/env python3
"""Soto's families as our families cut by copy number: T1 and T2 on one universe.

Spec docs/superpowers/specs/2026-10-14-soto-subset-cn-design.md; prereg docs/PREREG_soto_subset_cn_2026-10-14.md.

    run_subset.py guided --out DIR [--null 1000] [--boot 1000] [--seed 7] [--inputs-only] [--no-recipe]
    run_subset.py denovo --sample human_A119b|human_testis --out DIR [same options]

DIR receives inputs.json (universe and copy-number coverage; all that --inputs-only writes), t1.json, violations.tsv,
t2.json and, for guided, recipe.json: Soto's own recipe (`soto_replication.py cluster --pair-mad`) with each
copy-number source, scored on its full gene set (where the registered rungs must reproduce) and on T2's common set.
The run stops when the gene-id join disagrees with coordinates or a registered rung does not reproduce.
"""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cn_tables as CN  # noqa: E402
import subset_lib as L  # noqa: E402
import universes as U  # noqa: E402

SOTO = os.path.join(CN.REPO, "bench", "soto")
REPLICATION = os.path.join(SOTO, "soto_replication.py")
SHARED = os.path.join(SOTO, "shared_exons_5154_exon_mapback.tsv")
F_O = "/mnt/linuxdisk/tmp/o1_layers/v2/dupmasker"
F_O_GFF = "/mnt/linuxdisk/tmp/ourhom/loci.gff3"
SHIPPED = "/mnt/linuxdisk/tmp/o1_close/hier_shipped"
RUNS = "/mnt/linuxdisk/tmp/rustle_figures/runs"
SOURCES = ("assembly", "read_depth", "s1c")
# Soto's recipe with each source; the registered rungs (ARI, exact of 491) must reproduce (registers 1162, 1169).
RECIPE = (("none", ["--no-cn-gate"], (0.7307, 345)),
          ("s1c", ["--famcn", CN.S1C], (0.9698, 479)),
          ("read_depth", ["--famcn", CN.READ_DEPTH, "--famcn-col", "famCN_sotoiv"], (0.9277, 411)),
          ("assembly", ["--famcn", "{assembly}", "--famcn-col", "famCN"], None))


def load_tables(a):
    """The three copy-number tables by gene id; stops on any id/coordinate disagreement."""
    asm, unresolved = CN.assembly_cn(a.s1e)
    tables = dict(assembly=asm, read_depth=CN.read_depth_cn(a.read_depth), s1c=CN.s1c_cn(a.s1c))
    ok, bad, absent = CN.check_coordinates(a.s1e, a.truth_v2)
    if bad:
        raise SystemExit(f"the gene-id join disagrees with coordinates on {len(bad)} S1E rows, e.g. {bad[:3]}")
    return tables, unresolved, dict(ok=ok, disagree=0, gene_absent=absent)


def build_universe(a):
    if a.universe == "guided":
        return U.guided(a.s1c, a.fo_gff, os.path.join(a.fo, "ours.hierarchy.members.tsv"),
                        os.path.join(a.fo, "ours.loci.tsv"))
    d = os.path.join(a.shipped, a.sample)
    return U.denovo(a.truth_v2, os.path.join(a.runs, a.sample, f"{a.sample}.fam.loci.gff3"),
                    os.path.join(d, "ours.hierarchy.members.tsv"), os.path.join(d, "ours.loci.tsv"))


def common_members(member_cn):
    """The members with a value in every source, sorted."""
    return sorted(set.intersection(*(set(v) for v in member_cn.values())))


def inputs_report(uni, tables, unresolved, member_cn):
    genes = sorted(set().union(*uni.genes_of.values()))
    common = set(common_members(member_cn))
    return dict(coverage=uni.coverage, members=len(uni.block_of), soto_labels=len(uni.soto_sets),
                families_ge2=sum(1 for s in uni.soto_sets.values() if len(s) >= 2),
                single_label_members=len(uni.soto_of), multi_label_members=len(uni.multi),
                cn_genes={k: CN.account(genes, v, unresolved if k == "assembly" else frozenset())
                          for k, v in tables.items()},
                cn_members={k: len(v) for k, v in member_cn.items()},
                common_members=len(common),
                families_ge2_in_common=sum(1 for s in uni.soto_sets.values() if len(s & common) >= 2),
                fam90a_in_common=sum(1 for m in common if uni.soto_of.get(m) == L.FAM90A))


def run_t2(uni, member_cn, n_null, seed):
    common = common_members(member_cn)
    out = dict(common_members=len(common), common={}, own={})
    for name in SOURCES:
        cn = {m: member_cn[name][m] for m in common}
        out["common"][name] = L.t2(uni.soto_sets, uni.soto_of, uni.block_of, cn, n_null, seed)
        own = member_cn[name]
        out["own"][name] = L.scores(uni.soto_sets, uni.soto_of, L.cn_components(uni.block_of, own), sorted(own))
    out["ours_alone"] = L.scores(uni.soto_sets, uni.soto_of, uni.block_of, common)
    single = [m for m in common if m in uni.soto_of]
    out["ceiling"] = L.scores(uni.soto_sets, uni.soto_of, L.meet(uni.soto_of, uni.block_of), single)
    out["t2_pass"] = out["common"]["assembly"]["passes"]
    out["read_depth_passes"] = out["common"]["read_depth"]["passes"]
    return out


def run_recipe(out_dir, uni, assembly, common):
    """Soto's own recipe with each copy-number source (guided only: the recipe is gene-level)."""
    sys.path.insert(0, SOTO)
    import soto_replication as R
    g1793, g2334 = os.path.join(out_dir, "g1793.tsv"), os.path.join(out_dir, "g2334.tsv")
    subprocess.run([sys.executable, REPLICATION, "genesets", "--out-eligible", g1793, "--out-full", g2334],
                   check=True)
    asm_path = os.path.join(out_dir, "cn_assembly.tsv")
    CN.write_table(asm_path, assembly)
    clean, _ = R.load_truth(R.S1C)
    universe = R.restrict_universe(clean)
    rows = {}
    for arm, args, expect in RECIPE:
        pred_path = os.path.join(out_dir, f"recipe_{arm}.tsv")
        subprocess.run([sys.executable, REPLICATION, "cluster", "--pair-mad", "--shared", SHARED, "--geneset", g1793,
                        "--full-geneset", g2334, *[x.replace("{assembly}", asm_path) for x in args],
                        "--out", pred_path], check=True)
        pred = R.load_predicted(pred_path)
        full = R.score(clean, pred, universe)
        if expect and (abs(full["ari"] - expect[0]) > 5e-5 or full["n_exact"] != expect[1]):
            raise SystemExit(f"recipe {arm}: ARI {full['ari']:.4f}, exact {full['n_exact']} does not reproduce the "
                             f"registered {expect[0]} / {expect[1]}")
        comp_of = {g: pred.get(g) or f"U:{g}" for g in uni.block_of}
        rows[arm] = dict(full_gene_set=dict(ari=full["ari"], exact=full["n_exact"], truth_families=full["n_truth_fam"]),
                         common_set=L.scores(uni.soto_sets, uni.soto_of, comp_of, common))
    return rows


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True)
        fh.write("\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("universe", choices=["guided", "denovo"])
    ap.add_argument("--sample", choices=["human_A119b", "human_testis"], help="denovo: the shipped catalog")
    ap.add_argument("--out", required=True)
    ap.add_argument("--null", type=int, default=1000)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--inputs-only", action="store_true", help="write inputs.json and stop (no statistic)")
    ap.add_argument("--no-recipe", action="store_true", help="guided: skip Soto's recipe rows")
    ap.add_argument("--s1c", default=CN.S1C)
    ap.add_argument("--s1e", default=CN.S1E)
    ap.add_argument("--read-depth", default=CN.READ_DEPTH)
    ap.add_argument("--truth-v2", default=CN.TRUTH_V2)
    ap.add_argument("--fo", default=F_O, help="guided: the F_O run directory (ours.hierarchy.members.tsv, ours.loci.tsv)")
    ap.add_argument("--fo-gff", default=F_O_GFF)
    ap.add_argument("--shipped", default=SHIPPED)
    ap.add_argument("--runs", default=RUNS)
    a = ap.parse_args(argv)
    if a.universe == "denovo" and not a.sample:
        ap.error("denovo needs --sample")
    os.makedirs(a.out, exist_ok=True)
    tables, unresolved, coord = load_tables(a)
    uni = build_universe(a)
    member_cn = {k: U.member_cn(v, uni.genes_of) for k, v in tables.items()}
    inputs = inputs_report(uni, tables, unresolved, member_cn)
    inputs["coordinate_check"] = coord
    write_json(os.path.join(a.out, "inputs.json"), dict(args=vars(a), **inputs))
    if a.inputs_only:
        print(json.dumps(inputs, indent=1, sort_keys=True))
        return 0
    t1 = L.t1(uni.soto_sets, uni.soto_of, uni.block_of, uni.multi, uni.clustered, a.null, a.boot, a.seed)
    write_json(os.path.join(a.out, "t1.json"), t1)
    with open(os.path.join(a.out, "violations.tsv"), "w") as fh:
        fh.write("label\tsize\tblocks\tclass\n")
        fh.writelines(f"{v['label']}\t{v['size']}\t{v['blocks']}\t{v['cls']}\n" for v in t1["violations"])
    t2 = run_t2(uni, member_cn, a.null, a.seed)
    write_json(os.path.join(a.out, "t2.json"), t2)
    if a.universe == "guided" and not a.no_recipe:
        write_json(os.path.join(a.out, "recipe.json"),
                   run_recipe(a.out, uni, tables["assembly"], common_members(member_cn)))
    print(json.dumps(dict(T1a=t1["t1a_pass"], T1b=t1["t1b_pass"], T2=t2["t2_pass"], contained=t1["contained"],
                          families_ge2=t1["families_ge2"], D=t1["D"], D_ci95=t1["D_ci95"],
                          t2_assembly=t2["common"]["assembly"]["observed"]), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
