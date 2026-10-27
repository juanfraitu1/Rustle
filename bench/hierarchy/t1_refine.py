#!/usr/bin/env python3
"""CLI of the T1 scorer (docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md, KEY=dupfamboundary; T1 text as committed in draft 3a, git 8ecccd72).

One subcommand per call. Each writes JSON and a short text summary under /mnt/linuxdisk/tmp/hier_run/t1/out/.
  Phase R0   inputs, hash-paf, gate0, gate1 --part fo|cnduplicon|unitcover|collect, gate2 (skipped, R3a only), gate3, gate4 --part det|calib|collect
             [--start --end --B], gate5 [--real-5a, needs GO_R1], gate6, r0 (summary)
  Phase R1   t1a, t1b, t1c [--null governing|chromosome|region] [--layer F_O|F_S] [--inflation 2.0|4.0], t1e.
             Without --synthetic they refuse real labels unless <run dir>/GO_R1 exists and the gates 0, 1, 3, 4, 5, 6 have passed.
Light calls:  RLOCK_TIMEOUT=180 bash tools/rlock.sh light python3 -B bench/hierarchy/t1_refine.py <subcommand>
Tests:        python3 -B -m unittest bench/hierarchy/test_t1_lib.py bench/hierarchy/test_t1_refine.py
"""
import argparse
import json
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
import time
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import t1_lib as T  # noqa: E402

REPO = T.REPO
RUN, INPUTS = T.RUN, T.INPUTS
OUT = RUN + "/out"
PREREG = REPO + "/docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md"
PREREG_T1_TEXT_COMMIT = "8ecccd72"        # sections 3, 4, 10 and Amendment 1 as committed in draft 3a; later commits only append amendments (2, 3) for T2
FROZEN_DUP = REPO + "/bench/soto_m2/soto_m2_duplicons.py"
FROZEN_UNIT = REPO + "/bench/soto_m2/soto_m2_unit_cover_cn.py"

# ------------------------------------------------------------------------------------------------ registered inputs (section 13)
# label, key in default_paths() or an F_O file name, sha256 (first 16 hex), registered line count (header included), registered data rows
REGISTERED_INPUTS = [
    ("colors BED", "colors", "59533d25623ebf60", 205815, 205815),
    ("duplicons BED", "duplicons_bed", "2c0d0c94890e23f1", 254199, 254198),
    ("S1C", "s1c", "776a650716fd1a1b", 2573, 2572),
    ("split", "split", "8f2cba6c34984b69", 2335, 2334),
    ("exon map-back", "mapback", "90ef554f96ec0359", 12232, None),
    ("sd98_gene_exons", "sd98_exons", "61826d144db5896b", 21769, None),
    ("sd98_v1.bed", "sd98_regions", "d3e56dfa2cf81d69", 817, 817),
    ("famcn_ours_allwssd", "famcn_ours", "9067d31ef1838764", 5155, None),
    ("gene_v2_coords.bed", "gene_v2", "8f25eb939b442432", 2332, 2332),
    ("families_cn.json", "families_cn", "625fa960d8a6ee1e", 0, None),
    ("elig.tsv", "elig", "ffd86585d04bc5e5", 1794, None),
    ("full.tsv", "full", "c2bb5a994ff9cbff", 2335, None),
    ("rt_genes.tsv", "rt_genes", "d7b107123e5ee19e", 77, 76),
    ("F_O names.tsv", "F:names.tsv", "839412bd1a5deeb0", 2334, None),
    ("F_O loci.gff3", "F:loci.gff3", "8237e7e268432583", 16482, None),
    ("F_O ours.clusters.tsv", "F:ours.clusters.tsv", "8a835d4c434493a3", 1710, None),
    ("F_O ours.loci.tsv", "F:ours.loci.tsv", "7490ec5f23f5474a", 94, None),
    ("F_O ours.pairs.tsv", "F:ours.pairs.tsv", "add12d1a70a0c101", 10046, None),
    ("F_O ours.params.tsv", "F:ours.params.tsv", "9964b9aaa6158a78", 63, None),
    ("F_O ours.graph.tsv", "F:ours.graph.tsv", "1820795ae336a003", 11031, 11031),
    ("CAT/Liftoff v2.0 genes", "cat_genes", "b365d1995d75e2ea", 64214, None),
    ("CAT/Liftoff v2.0 refseq_map", "cat_refseq", "ca87ff8aac922e47", 58564, None),
    ("CAT v4 bed", "cat_v4", "6a1fddadcc8168ab", 232792, None),
]
REGISTERED_BIG = [                      # size only (loci.paf is hashed once, in its own call: hash-paf)
    ("F_O loci.paf", "F:loci.paf", 37782223),
    ("t2t-chm13-v1.0.fa.gz", "t2t_v1", 981798956),
]
REGISTERED_BINARY_SHA1 = "e4cc13b9d7eba5bcad3cc9c8c93a3fdc8ad434d7"
REGISTERED_TOOLS = {"python3": "3.14.4", "numpy": "2.4.2", "scipy": "1.17.1", "pysam": "0.23.3", "python (not used)": "3.13.12"}
FO_NAMES_BY_ROLE = {"graph": "ours.graph.tsv", "clusters": "ours.clusters.tsv", "loci": "ours.loci.tsv", "pairs": "ours.pairs.tsv",
                    "params": "ours.params.tsv"}

# ---------------------------------------------------------------------------------------------------- registered counts (Gate 6 and others)
# name -> (value, tier, where). tier: gate6 = listed in the Gate 6 text, registered = elsewhere in the prereg text, recon = only in the input-only
# reconnaissance outputs of the revision (scratch /mnt/linuxdisk/tmp/prereg_revision/logic). Every one is compared; a difference is reported.
REGISTERED_COUNTS = {
    "genes_total": (2334, "registered", "section 4 universe V"),
    "class_multi": (149, "gate6", "Gate 6, multi-family genes"),
    "genes_with_exonic_duplicon": (2290, "registered", "Read first 1 (2,290 duplicon-bearing S1C genes)"),
    "non_multi_without_duplicon": (43, "registered", "section 4 universe V (43 more removed)"),
    "V": (2142, "registered", "section 4"),
    "V_clean": (2036, "registered", "section 4"),
    "V_famless": (106, "registered", "section 4"),
    "genes_u1_ge2_all": (1062, "registered", "Read first 1 (1,062 of 2,290)"),
    "genes_u1_ge2_V": (974, "registered", "Read first 1 (974 of 2,142)"),
    "labelled_exons_V": (12597, "gate6", "Gate 6 ties"),
    "tie_exons_V": (857, "gate6", "Gate 6 ties"),
    "tie_genes_V": (316, "gate6", "Gate 6 ties"),
    "labelled_exons_all": (13501, "gate6", "Gate 6 ties"),
    "tie_exons_all": (898, "gate6", "Gate 6 ties"),
    "U1_distinct_V": (989, "gate6", "Gate 6 (989 to 985, 989 to 989)"),
    "U1_distinct_V_union": (985, "gate6", "Gate 6"),
    "U1_distinct_V_num": (989, "gate6", "Gate 6"),
    "num_diff_exons_V": (82, "gate6", "Gate 6 (82 exons in 54 genes differ)"),
    "num_diff_genes_V": (54, "gate6", "Gate 6"),
    "num_u1_changed_V": (53, "gate6", "Gate 6 (the U1 sets of 53 genes change)"),
    "U0_distinct_V": (629, "recon", "string_only_checks.out"),
    "U1_distinct_clean": (901, "registered", "section 2 item 1 (901 on the 2,036 clean genes)"),
    "U0_distinct_clean": (564, "registered", "section 4 (901 and 564 strings)"),
    "readthrough_V": (36, "registered", "section 4 exclusion arms"),
    "readthrough_all": (76, "registered", "section 2 item 10"),
    "non_eligible_V": (388, "registered", "section 4 exclusion arms"),
    "non_eligible_all": (541, "registered", "section 2 item 10"),
    "span_pairs_V": (105, "registered", "section 4 same-locus arm"),
    "span_genes_V": (181, "registered", "section 4 same-locus arm"),
    "span_pairs_all": (233, "registered", "section 2 item 10"),
    "span_genes_all": (339, "registered", "section 2 item 10"),
    "short_V": (862, "registered", "section 4 audit"),
    "short_all": (948, "registered", "section 2 item 10"),
    "fo_clusters": (380, "registered", "section 4 F_O"),
    "fo_rows": (1709, "registered", "section 4 F_O"),
    "fo_folded_records": (93, "registered", "section 4 F_O"),
    "fo_folded_representatives": (86, "registered", "section 4 F_O"),
    "fo_size2_rows": (163, "gate6", "Gate 6 (163 of 380 clusters have two locus rows)"),
    "fo_clusters_ge2_VF": (362, "gate6", "Gate 6 (152 of the 362 clusters with at least two V_F genes)"),
    "fo_size2_VF": (152, "gate6", "Gate 6 errata (n_f = 2 in the layer universe)"),
    "fo_size2_both": (147, "gate6", "Gate 6 errata (147 of them are among the 152)"),
    "VF": (1657, "registered", "section 4 V_F"),
    "VF_rule1": (1581, "registered", "section 4 F_O rule 1"),
    "VF_rule2": (76, "registered", "section 4 F_O rule 2"),
    "V_unclustered": (485, "registered", "section 4"),
    "V_never_connected": (467, "registered", "section 4 (467 are never connected)"),
    "V_singleton_nodes": (18, "registered", "section 4 (18 are singleton nodes)"),
    "V_folded": (80, "registered", "section 4 (80 V genes are folded annotations)"),
    "V_folded_rep_not_row": (4, "registered", "section 4 (4 of them have a representative that is not a row)"),
    "graph_nodes": (1819, "gate6", "Gate 1 / section 4"),
    "graph_edges": (11031, "gate6", "Gate 1 / section 4"),
    "graph_components": (339, "gate6", "Gate 6"),
    "strata_with_VF": (325, "gate6", "Gate 6"),
    "strata_variation": (30, "gate6", "Gate 6"),
    "strata_variation_genes": (416, "gate6", "Gate 6"),
    "strata_single": (295, "gate6", "Gate 6"),
    "strata_single_genes": (1241, "gate6", "Gate 6"),
    "strata_variation_min": (4, "registered", "section 4 (4 to 42 genes, median 11)"),
    "strata_variation_max": (42, "registered", "section 4"),
    "strata_variation_median": (11, "registered", "section 4"),
    "strata_variation_ge2multi": (29, "registered", "section 4 (29 hold at least two multi-gene clusters)"),
    "strata_variation_ge2multi_genes": (408, "recon", "c2_fo_components.out"),
    "fo_clusters_touching_VF": (372, "recon", "string_only_checks.out"),
    "fo_clusters_straddling_halves": (9, "registered", "section 4 Half"),
    "fo_genes_in_straddling_clusters": (61, "registered", "section 4 Half"),
    "fo_genes_in_clusters_ge2": (1647, "registered", "section 4 Half"),
    "V_dev": (1072, "recon", "halves.out"),
    "V_heldout": (1070, "recon", "halves.out"),
    "regions_V": (445, "registered", "section 4 sensitivity nulls"),
    "fs_families_V_clean": (486, "registered", "section 4 F_S"),
    "fs_families_V_clean_ge2": (434, "registered", "section 4 F_S"),
    "fs_families_all_clean": (491, "registered", "section 4 F_S"),
    "fs_families_all_clean_ge2": (444, "registered", "section 4 F_S"),
    "V_no_famcn_either": (388, "registered", "section 4 rung U2 (388 in V; which famCN source is not stated)"),
}
REGISTERED_FQ = {"fq_groups": (405, "registered", "section 4 F_Q"), "fq_grouped_genes": (2102, "registered", "section 4 F_Q"),
                 "fq_universe": (2185, "registered", "section 4 F_Q"), "fq_ungrouped": (83, "registered", "section 4 F_Q")}
REGISTERED_5B = {"D_S_to_Q_singleton": (7, "gate5", "Gate 5b"), "S_split_singleton": (4, "gate5", "Gate 5b (4 of 444 Soto families split)"),
                 "S_scored_singleton": (444, "gate5", "Gate 5b"), "D_S_to_Q_ignored": (2, "gate5", "Gate 5b")}
CNDUP_REGISTERED = {"clusters": 33, "soto_families": 87, "genes": 642, "heldout_delta": 0.2853, "dev_delta": 0.3276}
UNITCOVER_A = dict(dev=0.439, heldout=0.473, exact=0.107, recall=0.555, precision=0.718, own=0.546, multi=0.303, empty=0.020)
UNITCOVER_UNITS = {("A", "-"): 1546, ("B", "-"): 1241, ("C", "s1c"): 2255, ("C", "ours"): 2400, ("D", "s1c"): 1715, ("D", "ours"): 1837}
GATES_FOR_R1 = ("gate0", "gate1", "gate3", "gate4", "gate5", "gate6")


# ------------------------------------------------------------------------------------------------------------------------------ helpers
def manifest_row(path, hash_it=True):
    st = os.stat(path)
    row = dict(path=path, bytes=st.st_size, mtime=time.strftime("%Y-%m-%d %H:%M", time.localtime(st.st_mtime)), sha256=None, lines=None)
    if hash_it:
        row["sha256"] = T.sha256_file(path)
        row["lines"] = T.count_lines(path)
    return row


def compare_registered(registered, obtained):
    """One row per registered name: ok iff the obtained value equals it (a missing value never does)."""
    return [dict(name=k, registered=v, obtained=obtained.get(k), ok=k in obtained and obtained[k] == v) for k, v in registered.items()]


def all_ok(rows):
    return all(r["ok"] for r in rows)


def same_bytes(pa, pb):
    """Byte-for-byte equality (filecmp caches by size and mtime, which is not safe for files rewritten in the same second)."""
    with open(pa, "rb") as fa, open(pb, "rb") as fb:
        while True:
            x, y = fa.read(1 << 20), fb.read(1 << 20)
            if x != y:
                return False
            if not x:
                return True


def compare_dirs(a, b, names):
    """The names whose file differs between two directories, or is missing in either."""
    bad = []
    for n in names:
        pa, pb = os.path.join(a, n), os.path.join(b, n)
        if not (os.path.exists(pa) and os.path.exists(pb)) or not same_bytes(pa, pb):
            bad.append(n)
    return bad


def require_go(d):
    """The decisive cells refuse real labels until the coordinator has put a GO_R1 file in `d` (after the independent code review)."""
    if not os.path.exists(os.path.join(d, "GO_R1")):
        sys.exit(f"HARD STOP: the decisive cells are not run on real labels before the independent code review ({os.path.join(d, 'GO_R1')} is absent)")


def gate_status(out_dir, name):
    p = os.path.join(out_dir, name + ".json")
    if not os.path.exists(p):
        return "MISSING"
    with open(p) as fh:
        return json.load(fh).get("status", "MISSING")


def require_gates(out_dir, gates=GATES_FOR_R1):
    """A cell is run only after its gates pass (section 3)."""
    bad = {g: gate_status(out_dir, g) for g in gates if gate_status(out_dir, g) != "PASS"}
    if bad:
        sys.exit(f"the gates that cover Phase R1 have not all passed: {bad}")


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True, default=_jsonable)
        fh.write("\n")


def _jsonable(o):
    if isinstance(o, (set, frozenset)):
        return sorted(o, key=str)
    if isinstance(o, tuple):
        return list(o)
    raise TypeError(f"not serialisable: {type(o)}")


def write_text(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text if text.endswith("\n") else text + "\n")


def say(text):
    print(text, flush=True)


def prereg_stamp():
    """The commit that holds the prereg and the sha256 of the file on disk (recorded in every result)."""
    try:
        commit = subprocess.run(["git", "-C", REPO, "log", "-1", "--format=%h", "--", "docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md"],
                                capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:
        commit = "?"
    return dict(last_commit=commit, t1_text_commit=PREREG_T1_TEXT_COMMIT, sha256=T.sha256_file(PREREG, 64) if os.path.exists(PREREG) else None)


def file_path_of(key):
    if key.startswith("F:"):
        return os.path.join(INPUTS, key[2:])
    return T.default_paths()[key]


# --------------------------------------------------------------------------------------------------------------------- inputs, Gate 0
def cmd_inputs():
    """Copy the F_O files to the hashed location (cp -p semantics), verify them against the registered sha256 before and after, make them
    read-only and write MANIFEST.tsv. loci.paf is copied and compared byte for byte; its sha256 is computed once, in its own call (hash-paf)."""
    os.makedirs(INPUTS, exist_ok=True)
    reg = {k[2:]: (sha, lines) for _, k, sha, lines, _ in REGISTERED_INPUTS if k.startswith("F:")}
    rows, ok = [], True
    for name, src in T.FO_SOURCES.items():
        dst = os.path.join(INPUTS, name)
        if name == "loci.paf":
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
            same = same_bytes(src, dst)
            row = manifest_row(dst, hash_it=False)
            row.update(name=name, source=src, verified="byte-for-byte copy" if same else "DIFFERENT", registered_bytes=37782223)
            ok &= same and row["bytes"] == 37782223
            rows.append(row)
            os.chmod(dst, 0o444)
            continue
        want = reg[name][0]
        got_src = T.sha256_file(src)
        if got_src != want:
            ok = False
            say(f"SOURCE HASH DIFFERS {name}: {got_src} != {want}")
            continue
        if not os.path.exists(dst):
            shutil.copy2(src, dst)
        row = manifest_row(dst)
        row.update(name=name, source=src, registered=want, verified="match" if row["sha256"] == want else "DIFFERENT")
        ok &= row["sha256"] == want
        rows.append(row)
        os.chmod(dst, 0o444)
    head = "name\tsha256\tbytes\tlines\tmtime\tsource\tverified"
    lines = [head] + ["\t".join(str(r.get(k)) for k in ("name", "sha256", "bytes", "lines", "mtime", "source", "verified")) for r in rows]
    mf = os.path.join(INPUTS, "MANIFEST.tsv")
    if os.path.exists(mf):
        os.chmod(mf, 0o644)
    write_text(mf, "\n".join(lines))
    os.chmod(mf, 0o444)
    say("\n".join(lines))
    say("inputs: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def cmd_hash_paf():
    """sha256 of loci.paf, once, in its own call (size, mtime and hash are recorded)."""
    p = os.path.join(INPUTS, "loci.paf")
    t0 = time.time()
    row = manifest_row(p, hash_it=False)
    row["sha256"] = T.sha256_file(p, 64)
    row["seconds"] = round(time.time() - t0, 2)
    write_json(OUT + "/hash_paf.json", row)
    say(json.dumps(row, indent=1))
    return 0


def tool_versions():
    import numpy
    import scipy
    v = {"python3": sys.version.split()[0], "numpy": numpy.__version__, "scipy": scipy.__version__}
    try:
        import pysam
        v["pysam"] = pysam.__version__
    except Exception as e:                                   # pragma: no cover
        v["pysam"] = f"missing ({e})"
    try:
        v["python (not used)"] = subprocess.run(["python", "--version"], capture_output=True, text=True, timeout=30).stdout.split()[-1]
    except Exception:                                        # pragma: no cover
        v["python (not used)"] = "?"
    return v


def sha1_file(path):
    import hashlib
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cmd_gate0():
    """Path, sha256 and row count of every registered input, the interpreter and the tools. A hash that differs from its registered value fails."""
    rows, ok = [], True
    for label, key, sha, lines, nrows in REGISTERED_INPUTS:
        p = file_path_of(key)
        r = manifest_row(p)
        r.update(label=label, registered_sha256=sha, registered_lines=lines, registered_rows=nrows)
        r["rows"] = nrows if nrows is not None and r["lines"] == lines else None      # data rows, from the registered table
        r["ok"] = r["sha256"] == sha and r["lines"] == lines
        ok &= r["ok"]
        rows.append(r)
    big = []
    for label, key, size in REGISTERED_BIG:
        p = file_path_of(key)
        r = manifest_row(p, hash_it=False)
        r.update(label=label, registered_bytes=size, ok=r["bytes"] == size)
        if key == "F:loci.paf" and os.path.exists(OUT + "/hash_paf.json"):
            with open(OUT + "/hash_paf.json") as fh:
                hp = json.load(fh)
            r["sha256"] = hp["sha256"]
            r["ok"] &= hp["bytes"] == size
        ok &= r["ok"]
        big.append(r)
    binary = dict(path=T.default_paths()["mcl_binary"], sha1=sha1_file(T.default_paths()["mcl_binary"]), registered=REGISTERED_BINARY_SHA1)
    binary["ok"] = binary["sha1"] == REGISTERED_BINARY_SHA1
    ok &= binary["ok"]
    tools = tool_versions()
    tools_ok = {k: tools.get(k) == v for k, v in REGISTERED_TOOLS.items()}
    ok &= all(tools_ok.values())
    env = dict(flag_dont_write_bytecode=sys.flags.dont_write_bytecode, PYTHONHASHSEED=os.environ.get("PYTHONHASHSEED"),
               executable=sys.executable, platform=sys.platform)
    ok &= bool(sys.flags.dont_write_bytecode)
    code = {n: T.sha256_file(os.path.join(HERE, n), 16) for n in ("t1_lib.py", "t1_refine.py", "test_t1_lib.py", "test_t1_refine.py")}
    res = dict(gate="gate0", status="PASS" if ok else "FAIL", inputs=rows, big=big, binary=binary, tools=tools, tools_ok=tools_ok, env=env,
               prereg=prereg_stamp(), code_sha256=code)
    write_json(OUT + "/gate0.json", res)
    lines = [f"{'ok' if r['ok'] else 'FAIL':4} {r['sha256']}  {r['lines']:>8} lines  {r['bytes']:>11} B  {r['label']}  {r['path']}" for r in rows]
    lines += [f"{'ok' if r['ok'] else 'FAIL':4} size {r['bytes']} B (registered {r['registered_bytes']}), mtime {r['mtime']}"
              f"{', sha256 ' + r['sha256'] if r.get('sha256') else ''}  {r['label']}" for r in big]
    lines.append(f"{'ok' if binary['ok'] else 'FAIL':4} binary sha1 {binary['sha1']}  {binary['path']}")
    lines.append("tools: " + ", ".join(f"{k} {tools[k]}{'' if tools_ok[k] else ' (REGISTERED ' + REGISTERED_TOOLS[k] + ')'}" for k in tools))
    lines.append(f"python3 -B: {bool(sys.flags.dont_write_bytecode)}; PYTHONHASHSEED={env['PYTHONHASHSEED']}")
    lines.append(f"prereg file: last commit {res['prereg']['last_commit']}, T1 text as committed in {PREREG_T1_TEXT_COMMIT}, sha256 {str(res['prereg']['sha256'])[:16]}")
    lines.append("code sha256: " + ", ".join(f"{k} {v}" for k, v in code.items()))
    lines.append("Gate 0: " + res["status"])
    write_text(OUT + "/gate0.txt", "\n".join(lines))
    say("\n".join(lines))
    return 0 if ok else 1


# ----------------------------------------------------------------------------------------------------------- universe-level tables
def load_universe(tie="lex", union=False):
    return T.build_universe(T.default_paths(), tie=tie, union=union)


def _vf_fam(U):
    return {g: U.fo_label[g] for g in U.VF}


def denominators(U):
    """Every count of Gate 6 for the loaded universe (the layer F_O and the rungs; no family label is joined to a unit string)."""
    d = {}
    G = U.genes
    cls = Counter(U.info[g]["cls"] for g in G)
    d.update(genes_total=len(G), class_clean=cls["clean"], class_multi=cls["multi"], class_famless=cls["famless"])
    d["genes_with_exonic_duplicon"] = sum(1 for g in G if U.strings[g].pooled)
    d.update(V=len(U.V), V_clean=len(U.V_clean), V_famless=len(U.V_famless))
    d["non_multi_without_duplicon"] = sum(1 for g in G if U.info[g]["cls"] != "multi" and not U.strings[g].pooled)
    d["multi_family_genes"] = cls["multi"]
    V, Vset = list(U.V), set(U.V)
    # exons and ties
    def exon_counts(genes):
        merged = labelled = ties = tie_genes = 0
        for g in genes:
            rows = T.exon_label_rows(U.colors, U.info[g]["chrom"], U.info[g]["exons"])
            merged += len(rows)
            labelled += sum(r["labelled"] for r in rows)
            t = sum(r["tie"] for r in rows)
            ties += t
            tie_genes += t > 0
        return merged, labelled, ties, tie_genes
    ma, la, ta, tga = exon_counts(G)
    mv, lv, tv, tgv = exon_counts(V)
    d.update(merged_exons_all=ma, labelled_exons_all=la, tie_exons_all=ta, tie_genes_all=tga,
             merged_exons_V=mv, labelled_exons_V=lv, tie_exons_V=tv, tie_genes_V=tgv)
    num_exons = num_genes = changed = 0
    u1_union, u1_num = set(), set()
    for g in V:
        i = U.info[g]
        rows = T.exon_label_rows(U.colors, i["chrom"], i["exons"])
        k = sum(1 for r in rows if r["labelled"] and r["lex"] != r["num"])
        num_exons += k
        num_genes += k > 0
        su = T.gene_strings(i["exons"], i["chrom"], U.colors, tie="lex", union=True).u1
        sn = T.gene_strings(i["exons"], i["chrom"], U.colors, tie="num").u1
        u1_union.add(su)
        u1_num.add(sn)
        changed += sn != U.strings[g].u1
    d.update(U1_distinct_V=len({U.strings[g].u1 for g in V}), U1_distinct_V_union=len(u1_union), U1_distinct_V_num=len(u1_num),
             num_diff_exons_V=num_exons, num_diff_genes_V=num_genes, num_u1_changed_V=changed,
             U0_distinct_V=len({U.strings[g].u0 for g in V}), U1_distinct_clean=len({U.strings[g].u1 for g in U.V_clean}),
             U0_distinct_clean=len({U.strings[g].u0 for g in U.V_clean}),
             genes_u1_ge2_all=sum(1 for g in G if len(U.strings[g].u1) >= 2), genes_u1_ge2_V=sum(1 for g in V if len(U.strings[g].u1) >= 2))
    # exclusion arms
    d.update(readthrough_V=len(U.readthrough & Vset), readthrough_all=len(U.readthrough & set(G)),
             non_eligible_V=sum(1 for g in V if not U.info[g]["eligible"]), non_eligible_all=sum(1 for g in G if not U.info[g]["eligible"]),
             short_V=sum(1 for g in V if U.info[g]["nex"] <= 2), short_all=sum(1 for g in G if U.info[g]["nex"] <= 2))
    pv, _ = T.same_locus_pairs(U.info, V)
    pa, _ = T.same_locus_pairs(U.info, G)
    d.update(span_pairs_V=len(pv), span_genes_V=len({g for p in pv for g in p}), span_pairs_all=len(pa), span_genes_all=len({g for p in pa for g in p}))
    # famCN of rung U2
    d.update(V_no_famcn_cs=sum(1 for g in V if U.info[g]["cs"] is None), V_no_famcn_co=sum(1 for g in V if U.info[g]["co"] is None),
             V_no_famcn_either=sum(1 for g in V if U.info[g]["cs"] is None or U.info[g]["co"] is None))
    # the layer F_O
    fo = U.fo
    sizes = fo.rows
    d.update(fo_clusters=len(sizes), fo_rows=sum(sizes.values()), fo_size2_rows=sum(1 for v in sizes.values() if v == 2),
             fo_folded_records=len(fo.fold), fo_folded_representatives=len(set(fo.fold.values())))
    vf_count = Counter(U.fo_label[g] for g in U.VF)
    size2_vf = {c for c, v in vf_count.items() if v == 2}
    d.update(fo_clusters_touching_VF=len(vf_count), fo_clusters_ge2_VF=sum(1 for v in vf_count.values() if v >= 2),
             fo_clusters_lt2_VF=sum(1 for v in vf_count.values() if v < 2), fo_clusters_empty_VF=len(sizes) - len(vf_count),
             fo_size2_VF=len(size2_vf), fo_size2_both=len(size2_vf & {c for c, v in sizes.items() if v == 2}))
    d.update(VF=len(U.VF), VF_rule1=sum(1 for g in U.VF if U.fo_rule[g] == 1), VF_rule2=sum(1 for g in U.VF if U.fo_rule[g] == 2),
             V_unclustered=sum(1 for g in V if U.fo_label[g] is None))
    folded = [g for g in V if fo.names.get(g) in fo.fold]
    d.update(V_folded=len(folded), V_folded_rep_not_row=sum(1 for g in folded if fo.fold[fo.names[g]] not in fo.cluster_of),
             V_never_connected=sum(1 for g in V if U.gene_node[g] is None),
             V_singleton_nodes=sum(1 for g in V if U.gene_node[g] is not None and U.fo_label[g] is None))
    d.update(graph_nodes=len(fo.nodes), graph_edges=len(fo.edges), graph_components=len(set(fo.comp.values())))
    fam = _vf_fam(U)
    var = T.variation_strata(fam, U.strata)
    var_ids = {sid for sid, _ in var}
    single = [(sid, gs) for sid, gs in U.strata if sid not in var_ids]
    ge2m = [(sid, gs) for sid, gs in var if sum(1 for v in Counter(fam[g] for g in gs).values() if v >= 2) >= 2]
    vsz = sorted(len(gs) for _, gs in var)
    d.update(strata_with_VF=len(U.strata), strata_variation=len(var), strata_variation_genes=sum(len(gs) for _, gs in var),
             strata_single=len(single), strata_single_genes=sum(len(gs) for _, gs in single),
             strata_variation_min=vsz[0] if vsz else None, strata_variation_max=vsz[-1] if vsz else None,
             strata_variation_median=statistics.median(vsz) if vsz else None,
             strata_variation_ge2multi=len(ge2m), strata_variation_ge2multi_genes=sum(len(gs) for _, gs in ge2m))
    # halves
    ge2 = [c for c, v in vf_count.items() if v >= 2]
    members = defaultdict(list)
    for g in U.VF:
        members[U.fo_label[g]].append(g)
    straddle = [c for c in ge2 if len({U.info[g]["half"] for g in members[c]}) > 1]
    d.update(V_dev=sum(1 for g in V if U.info[g]["half"] == "dev"), V_heldout=sum(1 for g in V if U.info[g]["half"] == "heldout"),
             fo_clusters_straddling_halves=len(straddle), fo_genes_in_straddling_clusters=sum(len(members[c]) for c in straddle),
             fo_genes_in_clusters_ge2=sum(len(members[c]) for c in ge2))
    d["regions_V"] = len({U.region[g] for g in V if U.region[g] is not None})
    # Soto families of the clean universe (family sizes only)
    def soto_counts(genes):
        c = Counter(next(iter(U.info[g]["fams"])) for g in genes)
        return len(c), sum(1 for v in c.values() if v >= 2)
    d["fs_families_V_clean"], d["fs_families_V_clean_ge2"] = soto_counts(U.V_clean)
    d["fs_families_all_clean"], d["fs_families_all_clean_ge2"] = soto_counts([g for g in G if U.info[g]["cls"] == "clean"])
    # the multi-family genes and the clusters that hold them
    multi = [g for g in G if U.info[g]["cls"] == "multi"]
    d.update(fo_clusters_with_multi_family_genes=len({U.fo_label[g] for g in multi if U.fo_label[g] is not None}),
             fo_members_not_in_V_multi=sum(1 for g in multi if U.fo_label[g] is not None),
             fo_members_not_in_V_no_duplicon=sum(1 for g in G if U.info[g]["cls"] != "multi" and not U.strings[g].pooled and U.fo_label[g] is not None))
    d["exons_match_families_cn"] = bool(U.exons_match_families_cn)
    return d


def strings_table(U):
    """The unit strings of the V genes (no family label of any kind), sorted by gene id: the table compared by Gate 3."""
    lines = ["gene\tu0\tu1\tordered\tn_exons"]
    for g in sorted(U.V):
        st = U.strings[g]
        lines.append("\t".join([g, st.u0 or "", ",".join(sorted(st.u1)), ">".join(st.ordered), str(st.n_exons)]))
    return "\n".join(lines) + "\n"


def multi_family_table(U):
    """The listing of the multi-family genes with their families and strings. It lies outside every statistic and is written to a file."""
    lines = ["gene\tname\tbiotype\tsoto_families\tu0\tu1\tfo_cluster\tfo_rule\treadthrough"]
    for g in sorted(U.genes):
        i = U.info[g]
        if i["cls"] != "multi":
            continue
        st = U.strings[g]
        lines.append("\t".join([g, i["name"], i["biotype"], ",".join(sorted(i["fams"])), st.u0 or "", ",".join(sorted(st.u1)),
                                U.fo_label[g] or "", str(U.fo_rule[g] or ""), str(g in U.readthrough)]))
    return "\n".join(lines) + "\n"


def fq_partition():
    """F_Q: the sequence groups collapse_cover(pair_families(..., gate=False)) exactly as fq_checks.py builds them. Frozen module, not edited."""
    sys.path.insert(0, REPO + "/bench/soto")
    import soto_replication as sr
    P = T.default_paths()
    genes_e, _ = sr.load_geneset(P["elig"])
    full, _ = sr.load_geneset(P["full"])
    edges = sr.read_edges(P["mapback"])
    _clean, ambiguous = sr.load_truth(P["s1c"])
    cover, _kept, leaf_of = sr.pair_families(edges, genes_e, full, {}, gate=False)
    pred = sr.collapse_cover(cover, leaf_of, genes_e)
    universe = (full | genes_e) - ambiguous
    Q = {g: pred[g] for g in universe if pred.get(g)}
    return Q, universe


def tripwire_rows(den, registered=REGISTERED_COUNTS):
    rows = compare_registered({k: v[0] for k, v in registered.items()}, den)
    for r in rows:
        r["tier"], r["where"] = registered[r["name"]][1], registered[r["name"]][2]
    return rows


def format_denominators(den, fq=None):
    keys = [k for k in den if not isinstance(den[k], (list, dict))]
    lines = ["DENOMINATORS (printed before any statistic)"]
    w = max(len(k) for k in keys)
    for k in keys:
        lines.append(f"  {k:<{w}}  {den[k]}")
    if fq:
        for k, v in fq.items():
            lines.append(f"  {k:<{w}}  {v}")
    return "\n".join(lines)


def format_tripwires(rows):
    lines = ["TRIPWIRES (registered value against obtained)"]
    for r in rows:
        lines.append(f"  {'ok' if r['ok'] else 'MISMATCH':8} {r['name']:<34} registered {str(r['registered']):>8}  obtained {str(r['obtained']):>8}   [{r['tier']}] {r['where']}")
    return "\n".join(lines)


def cmd_gate6():
    U = load_universe()
    den = denominators(U)
    fq = {}
    try:
        Q, universe = fq_partition()
        fq = dict(fq_groups=len(set(Q.values())), fq_grouped_genes=len(Q), fq_universe=len(universe), fq_ungrouped=len(universe - set(Q)))
    except Exception as e:                                    # pragma: no cover
        fq = dict(fq_error=str(e))
    rows = tripwire_rows(den) + tripwire_rows(fq, REGISTERED_FQ)
    mf = multi_family_table(U)
    write_text(OUT + "/multi_family_genes.tsv", mf)                     # written, not printed
    ok = all_ok(rows)
    res = dict(gate="gate6", status="PASS" if ok else "FAIL", denominators=den, fq=fq, tripwires=rows,
               multi_family_table=dict(path=OUT + "/multi_family_genes.tsv", rows=mf.count("\n") - 1), prereg=prereg_stamp())
    write_json(OUT + "/gate6.json", res)
    text = format_denominators(den, fq) + "\n\n" + format_tripwires(rows)
    text += f"\n\nmulti-family genes listed in {OUT}/multi_family_genes.tsv ({res['multi_family_table']['rows']} rows; written, not inspected)"
    bad = [r["name"] for r in rows if not r["ok"]]
    text += "\nGate 6: " + res["status"] + (f"  (mismatches: {', '.join(bad)})" if bad else "")
    write_text(OUT + "/gate6.txt", text)
    say(text)
    return 0 if ok else 1


# ------------------------------------------------------------------------------------------------------- worlds, arms and the cells
def world_of_universe(U, delete=(), half=None):
    """Everything a cell reads about V_F (or an arm / half of it), as plain dicts: labels, unit strings, strata, copy numbers, nodes, edges."""
    delete = set(delete)
    genes = [g for g in U.VF if g not in delete and (half is None or U.info[g]["half"] == half)]
    return dict(genes=genes, fam={g: U.fo_label[g] for g in genes}, u1={g: U.strings[g].u1 for g in genes},
                u0={g: U.strings[g].u0 for g in genes}, ordered={g: U.strings[g].ordered for g in genes},
                info={g: U.info[g] for g in genes}, strata=T.build_strata(genes, U.gene_node, U.fo.comp),
                half={g: U.info[g]["half"] for g in genes}, co={g: U.info[g]["co"] for g in genes}, cs={g: U.info[g]["cs"] for g in genes},
                gene_node={g: U.gene_node[g] for g in genes}, edges=U.fo.weight, region_of={g: U.region[g] for g in genes})


def arm_deletions(U):
    """The genes each exclusion arm deletes from V, defined from the readthrough list, the S1C eligibility column and gene coordinates alone."""
    V = list(U.V)
    rt = {g for g in V if g in U.readthrough}
    ne = {g for g in V if not U.info[g]["eligible"]}
    _pairs, losers = T.same_locus_pairs(U.info, V)
    return {"none": set(), "readthrough": rt, "non_eligible": ne, "same_locus": set(losers), "all": rt | ne | set(losers)}


def denominator_guard(removed, n_v):
    """An arm that removes more than 10% of V carries the flag 'denominator guard'."""
    return removed > 0.10 * n_v


def _stats(fam, strings, ordered, genes):
    out = dict(n=len(genes), m1=T.m1(fam, strings, genes), m1c=T.m1c(fam, strings, genes), m2=T.m2(fam, strings, genes),
               macro_d1=T.macro_d1(fam, strings, genes), macro_d2=T.macro_d2(fam, strings, genes))
    if ordered is not None:
        out["m1t"] = T.m1t(fam, ordered, genes)
        out["delta"], out["mult"] = T.path_defect(fam, ordered, strings, genes)
    return out


def t1a_report(w):
    """T1-a on a world: m1c (D1) and m2 (D2) printed separately, REFUTED iff either is at least 1, the full exception lists, the audit counts
    with and without 'short', and the descriptive companions (m1, m1t, macro versions, path defect, the size-2 stratum)."""
    genes, fam, u1, info = w["genes"], w["fam"], w["u1"], w["info"]
    r = _stats(fam, u1, w.get("ordered"), genes)
    cnt = Counter(fam[g] for g in genes)
    r["size2"] = _stats(fam, u1, None, [g for g in genes if cnt[fam[g]] == 2])
    exc1, exc2 = T.exceptions_d1(fam, u1, info, genes), T.exceptions_d2(fam, u1, info, genes)
    r["exceptions_d1"], r["exceptions_d2"] = exc1, exc2
    l1, l2 = Counter(e["label"] for e in exc1), Counter(e["label"] for e in exc2)
    r["labels_d1"], r["labels_d2"] = dict(l1), dict(l2)
    r["m1c_without_short"], r["m2_without_short"] = len(exc1) - l1["short"], len(exc2) - l2["short"]
    r["audit_d1"] = [dict(family=f, **rec) for f, rec in sorted(T.audit_d1(fam, u1, info, genes).items(), key=lambda kv: str(kv[0]))]
    r["audit_d2"] = [dict(string=list(k), **rec) for k, rec in sorted(T.audit_d2(fam, u1, info, genes).items())]
    r["verdict"] = "REFUTED" if r["m1c"] >= 1 or r["m2"] >= 1 else "NOT-REFUTED"
    r["verdict_scope"] = "primary"
    return r


def t1bc_report(w, stat, B=T.B_DRAWS, seed=T.SEED_T1):
    return T.run_cell(w["fam"], w["u1"], w["strata"], stat=stat, B=B, seed=seed)


def format_cell(c):
    """The text of one decisive cell. The numbers of variation-bearing and live strata come first, before any p-value."""
    cell = {"m1c": "T1-b (D1, m1c)", "m2": "T1-c (D2, m2)"}[c["stat"]]
    lines = [f"{cell}: variation-bearing strata {c['n_var']}; live strata {c['n_live']} (a stratum is live if its own statistic takes at least two values over the {c['B']} draws)",
             f"observed {c['stat']} = {c['obs']}; draws mean {c['draw_mean']:.2f}, min {c['draw_min']}, max {c['draw_max']}",
             f"p = {c['p']:.6f} (lower tail, (1 + #{{draws <= observed}}) / (B + 1)); upper tail {c['p_upper']:.6f} (no class)"]
    if c["removed_stratum"] is None:
        lines.append("robustness: no z_s is negative, the removal step is skipped; the class of the remaining strata is the class of all strata")
    else:
        lines.append(f"robustness: stratum {c['removed_stratum']} (rank {c['removed_rank']}) has the most negative z_s ({c['z'][c['removed_stratum']]:.3f}) and is removed; "
                     f"observed on the remaining genes {c['obs_removed']}, p = {c['p_removed']:.6f}")
    lines.append(f"class from all strata {c['class_all']}; class from the remaining strata {c['class_removed']}; CLASS ISSUED: {c['class']}")
    lines.append("z table: rank stratum n_genes n_clusters observed mean sd z live")
    for t in c["z_table"]:
        lines.append(f"  {t['rank']:>4} {t['stratum']:<14} {t['n_genes']:>3} {t['n_clusters']:>2} {t['obs']:>5} {t['mean']:>9.3f} {t['sd']:>8.3f} {t['z']:>8.3f} {'live' if t['live'] else 'dead'}")
    return "\n".join(lines)


def t1e_report(w, B=2000, seed=T.SEED_T1):
    """T1-e (descriptive, no p-value and no class): the pairs of the variation-bearing strata with a famCN (`co`) and a U1 label whose nodes share a
    direct edge of the dumped graph, the fixed-effects fits and the stratum bootstrap."""
    var = T.variation_strata(w["fam"], w["strata"])
    ok = [g for sid, gs in var for g in gs if w["co"].get(g) is not None and w["u1"].get(g)]
    rows = T.t1e_pairs(var, w["fam"], w["u1"], w["co"], w["gene_node"], w["edges"], w["region_of"])
    out = T.t1e_summary(rows, B=B, seed=seed, strata=[sid for sid, _ in var])      # the bootstrap resamples all variation-bearing strata, pairless ones included
    out["universe_genes"] = len(ok)
    out["variation_strata"] = len(var)
    return out


def format_t1e(s):
    lines = [f"T1-e (descriptive; no p-value, no class): universe {s['universe_genes']} genes in {s['variation_strata']} variation-bearing strata; "
             f"{s['n_pairs']} pairs in {s.get('strata', 0)} strata"]
    if s["n_pairs"]:
        lines.append(f"share with x_str = 1: {s['share_x_str']:.3f}; share with x_clu = 1: {s['share_x_clu']:.3f}; correlation of x_str and x_clu: {s['corr_x_str_x_clu']}")
        for k in sorted(k for k in s if "|" in k):
            f = s[k]
            lines.append(f"  {k}: " + "; ".join(f"{n} = {f['coef'][n]:+.4f} [{f['interval'][n][0]:+.4f}, {f['interval'][n][1]:+.4f}]" for n in f["coef"])
                         + ("  (collinear: minimum-norm solution)" if f["rank_deficient"] else ""))
    return "\n".join(lines)


def synthetic_t1_world(seed, n_strata=14):
    """A whole made-up T1 world (labels, strings, strata, copy numbers, graph, regions) for the synthetic runs of the decisive cells."""
    w = T.synthetic_world(seed, n_strata=n_strata, mode="independent", sizes=[12] * n_strata, nclu=(2, 3), vocab="ABC")
    rng = random.Random(seed + 1)
    genes = sorted(w["fam"])
    co, region, edges = {}, {}, {}
    for sid, gs in w["strata"]:
        base = rng.uniform(2, 20)
        for g in gs:
            co[g] = round(base + rng.uniform(-3, 3), 1)
            region[g] = ("chrS", int(sid[1:4]) * 1000, int(sid[1:4]) * 1000 + 500) if rng.random() < 0.5 else None
        for i, a in enumerate(gs):
            for b in gs[i + 1:]:
                if rng.random() < 0.5:
                    edges[("n_" + a, "n_" + b)] = round(rng.uniform(0.7, 1.0), 3)
    return dict(genes=genes, fam=w["fam"], u1=w["strings"], u0={g: sorted(w["strings"][g])[0] for g in genes},
                ordered={g: tuple(sorted(w["strings"][g])) for g in genes}, info=w["info"], strata=w["strata"],
                half={g: ("dev" if i % 2 else "heldout") for i, g in enumerate(genes)}, co=co, cs=dict(co), gene_node={g: "n_" + g for g in genes},
                edges=edges, region_of=region)


def run_decisive_synthetic(seed=0, B=200):
    """The four decisive cells on a made-up world (the only data they may be run on before the independent code review)."""
    w = synthetic_t1_world(seed)
    cells = T.run_cells(w["fam"], w["u1"], w["strata"], B=B)
    return {"t1a": t1a_report(w), "t1b": cells["m1c"], "t1c": cells["m2"], "t1e": t1e_report(w, B=min(B, 200))}


def guard_real_cell(run_dir=RUN, out_dir=OUT):
    """A decisive cell on real labels needs the coordinator's GO_R1 file (after the independent code review) and every covering gate to have passed."""
    require_go(run_dir)
    require_gates(out_dir)


# ------------------------------------------------------------------------------------------------------------------------------ Gate 4
def default_sizes(n=30, seed=7):
    """A made-up cluster-size shape (n strata of two to five clusters) used when no real shape is given; sizes only."""
    rng = random.Random(seed)
    out = []
    for i in range(n):
        k = rng.choice([2, 2, 2, 3, 3, 4, 5])
        out.append(sorted((rng.randint(2, 20) for _ in range(k)), reverse=True))
    out[0] = [30, 10]                                         # one clearly largest stratum
    return out


def _cells_row(name, res):
    """A p <= 0.01 clause in both cells."""
    parts, ok = [], True
    for stat in ("m1c", "m2"):
        c = res[stat]
        good = c["p"] <= 0.01
        ok &= good
        parts.append(f"{stat}: obs {c['obs']} / mean {c['draw_mean']:.1f}, p = {c['p']:.5f}, live {c['n_live']}/{c['n_var']}, class {c['class']}")
    return dict(name=name, ok=bool(ok), detail="; ".join(parts))


def gate4_det(B=T.B_DRAWS, sizes=None):
    """The deterministic synthetic controls of Gate 4 (T1). Every row is a clause of the registered text; the calibration clause (planted random
    partitions, 100 repeats) runs in blocks (gate4_calib_block) and is collected afterwards."""
    sizes = sizes or default_sizes()
    rows = []

    def add(name, ok, detail=""):
        rows.append(dict(name=name, ok=bool(ok), detail=detail))

    rng = random.Random(4)
    # 1. families equal to the string partition
    bad = 0
    for _ in range(20):
        k, n = rng.randint(2, 8), rng.randint(10, 60)
        strings = {f"g{i:03d}": frozenset({f"S{rng.randrange(k)}"}) for i in range(n)}
        fam = {g: next(iter(s)) for g, s in strings.items()}
        bad += (T.m1(fam, strings), T.m1c(fam, strings), T.m2(fam, strings)) != (0, 0, 0)
    add("families equal to the string partition give m1 = m1c = m2 = 0 (20 random universes)", bad == 0, f"failures {bad}")
    # 2. relabelling k genes
    bad = tried = 0
    for _ in range(60):
        k = rng.randint(1, 4)
        nA, nB = rng.randint(2 * k, 2 * k + 8), rng.randint(k + 1, k + 8)
        genes_a = [f"a{i:03d}" for i in range(nA)]
        genes_b = [f"b{i:03d}" for i in range(nB)]
        strings = {g: frozenset("a") for g in genes_a}
        strings.update({g: frozenset("b") for g in genes_b})
        fam = {g: "A" for g in genes_a}
        fam.update({g: "B" for g in genes_b})
        for g in genes_a[:k]:
            fam[g] = "B"                                    # k genes of A are relabelled into B; a stays the plurality of A and of its string
        tried += 1
        bad += (T.m1(fam, strings), T.m1c(fam, strings), T.m2(fam, strings)) != (k, k, k)
        # the other reading: k genes of A change their string for one of their own; A's native string stays the plurality
        fam2 = {g: "A" for g in genes_a}
        strings2 = dict.fromkeys(genes_a, frozenset("a"))
        for i, g in enumerate(genes_a[:k]):
            strings2[g] = frozenset({f"x{i}"})
        bad += (T.m1(fam2, strings2), T.m1c(fam2, strings2), T.m2(fam2, strings2)) != (k, k, 0)
    add("relabelling k genes of a family whose native string stays the plurality gives exactly k (m1, m1c, m2; both readings)", bad == 0, f"{tried} cases, failures {bad}")
    # 3. refinement
    bad = 0
    for _ in range(300):
        n, nf, ns = rng.randint(8, 40), rng.randint(1, 6), rng.randint(1, 5)
        genes = [f"g{i:03d}" for i in range(n)]
        fam = {g: f"F{rng.randrange(nf)}" for g in genes}
        s1 = {g: frozenset({f"S{rng.randrange(ns)}"}) for g in genes}
        s2 = {g: frozenset({next(iter(s1[g])) + str(rng.randrange(2))}) for g in genes}
        bad += not (T.m1(fam, s2) >= T.m1(fam, s1) and T.m2(fam, s2) <= T.m2(fam, s1))
    add("strict m1 is non-decreasing and strict m2 non-increasing when the string partition is refined (300 random partitions)", bad == 0, f"failures {bad}")
    # 4. one string, nested truncations, branches
    genes = [f"g{i}" for i in range(6)]
    fam = dict.fromkeys(genes, "F")
    one = dict.fromkeys(genes, frozenset("a"))
    trunc = dict(zip(genes, [frozenset("a"), frozenset("ab"), frozenset("abc"), frozenset("abcd"), frozenset("ab"), frozenset("a")]))
    branch = dict(zip(genes, [frozenset("ab"), frozenset("ac"), frozenset("ab"), frozenset("ac"), frozenset("a"), frozenset("ad")]))
    got = (T.m1c(fam, one), T.m1c(fam, trunc), T.m1c(fam, branch))
    add("synthetic families built as one string, nested truncations and branches give m1c 0, 0, > 0", got[0] == 0 and got[1] == 0 and got[2] > 0, f"m1c = {got}")
    # 5. path defect
    def delta(members):
        fam_ = {f"m{i}": "F" for i in range(len(members))}
        ordered_ = {f"m{i}": tuple(m) for i, m in enumerate(members)}
        strings_ = {g: frozenset(o) for g, o in ordered_.items()}
        return T.path_defect(fam_, ordered_, strings_)[0]
    got = (delta([("A", "B", "C")]), delta([("A",), ("B",)]), delta([("X", "A"), ("X", "B"), ("X", "C")]), delta([("A", "B"), ("B", "C"), ("C", "A")]))
    add("path-defect cases: path 0, two disjoint units 1, 3-leaf star 1, triangle 1", got == (0, 1, 1, 1), f"delta = {got}")
    # 6. exception lists do not depend on the input order
    genes = [f"g{i:03d}" for i in range(60)]
    base = {g: (f"F{rng.randrange(8)}", frozenset(rng.sample("abcd", rng.randint(1, 3)))) for g in genes}
    infos = {g: dict(nex=rng.randint(1, 5), interior=frozenset(rng.sample("xyz", rng.randint(0, 2)))) for g in genes}

    def lists(order):
        fam_ = {g: base[g][0] for g in order}
        str_ = {g: base[g][1] for g in order}
        inf_ = {g: infos[g] for g in order}
        return T.exceptions_d1(fam_, str_, inf_, list(order)), T.exceptions_d2(fam_, str_, inf_, list(order))
    ref = lists(genes)
    same = True
    for _ in range(10):
        order = genes[:]
        rng.shuffle(order)
        same &= lists(order) == ref
    add("the exception lists are identical under any input order of the genes (10 shuffles)", same)
    # 7, 8. UNDERPOWERED and UNINFORMATIVE
    w = T.synthetic_world(0, n_strata=12, mode="one_cluster_per_stratum")
    res = T.run_cells(w["fam"], w["strings"], w["strata"], B=min(B, 300))
    add("a family partition equal to the strata (zero variation-bearing strata) gives UNDERPOWERED in both cells",
        all(res[s]["class"] == "UNDERPOWERED" and res[s]["n_var"] == 0 for s in res), f"n_var {res['m1c']['n_var']}")
    for n in (8, 10):
        w = T.synthetic_world(n, n_strata=n, mode="one_string_per_stratum", nclu=(2, 2))
        res = T.run_cells(w["fam"], w["strings"], w["strata"], B=min(B, 300))
        add(f"{n} strata, each with two clusters and one string per stratum, give UNINFORMATIVE in both cells",
            all(res[s]["class"] == "UNINFORMATIVE" and res[s]["n_live"] == 0 and res[s]["n_var"] == n for s in res), f"n_var {res['m1c']['n_var']}, live {res['m1c']['n_live']}")
    w = T.synthetic_world(1, n_strata=7, mode="one_string_per_stratum", nclu=(2, 2))
    res = T.run_cells(w["fam"], w["strings"], w["strata"], B=min(B, 300))
    add("7 such strata give UNDERPOWERED (fewer than 8 variation-bearing strata comes first)", all(res[s]["class"] == "UNDERPOWERED" for s in res))
    # 9. a planted partition that follows the strings
    w = T.synthetic_world(0, n_strata=14, mode="follow", sizes=[10] * 14)
    rows.append(_cells_row("a planted partition that follows the strings gives p <= 0.01 in both cells", T.run_cells(w["fam"], w["strings"], w["strata"], B=B)))
    # 10-12. the documented tautology on the cluster sizes of the variation-bearing strata
    for mode, label in (("identity", "strings equal to the cluster identity"), ("pairs", "strings equal to unions of adjacent clusters (pairs, in strata of at least three clusters only; identity elsewhere)"),
                        ("split", "strings equal to a split of every cluster of at least three genes into two alternating types")):
        w = T.world_from_sizes(sizes, mode, seed=1)
        rows.append(_cells_row(f"cluster sizes of the variation-bearing strata, {label}: p <= 0.01 in both cells", T.run_cells(w["fam"], w["strings"], w["strata"], B=B)))
    big = [cl for cl in sizes if len(cl) >= 3]
    if len(big) >= 8:
        w = T.world_from_sizes(big, "pairs", seed=1)
        rows.append(_cells_row("the same with unions of adjacent clusters over the strata of at least three clusters only (the other strata left out): p <= 0.01 in both cells",
                               T.run_cells(w["fam"], w["strings"], w["strata"], B=B)))
    # 13. only the largest stratum follows the strings
    w = T.world_from_sizes(sizes, "largest_follows", seed=3)
    res = T.run_cells(w["fam"], w["strings"], w["strata"], B=B)
    largest = min(sid for sid, g in w["strata"] if len(g) == max(len(x) for _, x in w["strata"]))
    ok = True
    parts = []
    for stat in ("m1c", "m2"):
        c = res[stat]
        good = c["p"] <= 0.01 and c["class"] in ("PARTIAL", "NOT-DISTINGUISHED") and c["removed_stratum"] == largest
        ok &= good
        parts.append(f"{stat}: p = {c['p']:.5f}, removed {c['removed_stratum']} (largest {largest}), p after removal {c['p_removed']}, class {c['class']}, live {c['n_live']}")
    add("a planted set in which only the largest stratum follows the strings gives p <= 0.01 on all strata and a class below CUT-ALIGNED after the one-group removal", ok, "; ".join(parts))
    return rows


def gate4_calib_block(start, end, B=T.B_DRAWS, nrep=100, sizes=None):
    """Repeats start..end-1 of the calibration: strings independent of the cut on the cluster-size shape, B draws each, the p of both cells."""
    assert 0 <= start <= end <= nrep
    sizes = sizes or default_sizes()
    out = dict(start=start, end=end, B=B, nrep=nrep, m1c=[], m2=[])
    for r in range(start, end):
        w = T.world_from_sizes(sizes, "independent", seed=1000 + r)
        eng = T.NullEngine(w["fam"], w["strings"], w["strata"])
        obs = eng.observed()
        dr = eng.draws(B)
        out["m1c"].append(T.p_lower(obs["m1c"], dr["tot_m1c"]))
        out["m2"].append(T.p_lower(obs["m2"], dr["tot_m2"]))
    return out


def gate4_collect(blocks, nrep=100, limit=10):
    ps = {"m1c": [], "m2": []}
    seen = set()
    for b in blocks:
        for stat in ps:
            ps[stat] += b[stat]
        seen.update(range(b["start"], b["end"]))
    counts = {s: sum(1 for p in v if p <= 0.05) for s, v in ps.items()}
    ok = seen == set(range(nrep)) and all(c <= limit for c in counts.values())
    return dict(repeats=len(seen), complete=seen == set(range(nrep)), small_p_counts=counts, limit=limit, ok=bool(ok),
                mean_p={s: sum(v) / len(v) if v else None for s, v in ps.items()})


# ------------------------------------------------------------------------------------------------------------------------------ Gate 5
def gate5b_values(clean, Q):
    """D[S -> Q] for Soto families S (clean genes by family) and the sequence groups Q, with ungrouped genes as singleton cells and ignored."""
    fam = defaultdict(list)
    for g, f in clean.items():
        fam[f].append(g)
    a = T.d_partition(fam, Q.get)
    b = T.d_partition(fam, Q.get, ignore_ungrouped=True)
    return dict(D_S_to_Q_singleton=a["moved"], S_split_singleton=a["split"], S_scored_singleton=a["blocks"],
                D_S_to_Q_ignored=b["moved"], S_split_ignored=b["split"], S_scored_ignored=b["blocks"])


def soto_clean_from_families_cn(path):
    """Soto clean families as fq_checks.py reads them: the non-Unassigned labels of the `sf` field of families_cn.json."""
    with open(path) as fh:
        G = json.load(fh)["genes"]
    S = {x["id"]: {f for f in x["sf"] if not f.startswith("Unassigned")} for x in G}
    return {g: next(iter(s)) for g, s in S.items() if len(s) == 1}


def surrogate_partitions(genes, seed):
    """Made-up family partitions of the given genes (labels drawn here, no real label involved): 20, 100 and n/2 random clusters, geometric
    sizes, one big cluster with singletons, and three giant clusters."""
    rng = random.Random(seed)
    genes = list(genes)
    n = len(genes)
    out = [(f"random_{k}_clusters", {g: f"R{rng.randrange(max(1, min(k, n)))}" for g in genes}) for k in (20, 100, max(1, n // 2))]
    order = genes[:]
    rng.shuffle(order)
    lab, i, c = {}, 0, 0
    while i < n:
        size = 1 + int(rng.expovariate(1 / 6))
        for g in order[i:i + size]:
            lab[g] = f"G{c}"
        i += size
        c += 1
    out.append(("geometric_sizes", lab))
    big = {g: ("BIG" if j < n // 3 else f"S{j}") for j, g in enumerate(order)}
    out.append(("one_big_cluster_and_singletons", big))
    out.append(("three_giant_clusters", {g: f"T{rng.randrange(3)}" for g in genes}))
    return out


def gate5a_surrogate(u1, cn, seeds=(0, 1, 2), source="cs", _break=False):
    """Gate 5a on made-up family partitions of the real common universe (genes with a famCN): strict m1 and m1c do not decrease and strict m2 does
    not increase from U1 to U2. The partitions are surrogates, so no real family label meets a real unit string."""
    cu = sorted(g for g in u1 if cn.get(g) is not None and u1[g])
    s1 = {g: u1[g] for g in cu}
    s2 = T.u2_strings(s1, {g: cn[g] for g in cu})
    if _break:
        s2 = dict.fromkeys(cu, frozenset({"all"}))               # not a refinement of U1: the scorer must notice
    groups = defaultdict(set)
    for g in cu:
        groups[s2[g]].add(s1[g])
    refines = all(len(v) == 1 for v in groups.values())
    rows = []
    for seed in seeds:
        for name, fam in surrogate_partitions(cu, seed):
            a = (T.m1(fam, s1), T.m1c(fam, s1), T.m2(fam, s1))
            b = (T.m1(fam, s2), T.m1c(fam, s2), T.m2(fam, s2))
            rows.append(dict(source=source, seed=seed, partition=name, genes=len(cu), u1=dict(m1=a[0], m1c=a[1], m2=a[2]), u2=dict(m1=b[0], m1c=b[1], m2=b[2]),
                             u2_refines_u1=refines, ok=bool(a[0] <= b[0] and a[1] <= b[1] and a[2] >= b[2])))
    return rows


# ------------------------------------------------------------------------------------------------------------------------------ Gate 1
def parse_unitcover_table(text):
    rows = {}
    pat = re.compile(r"^\|\s*([ABCD])\s*\|\s*(\S+)\s*\|\s*(\d+)\s*\|\s*([\d.]+) / ([\d.]+)\s*\|\s*[\d.]+ / [\d.]+\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|"
                     r"\s*([\d.]+)\s*\|\s*([\d.]+) / ([\d.]+) / ([\d.]+)\s*\|")
    for line in text.splitlines():
        m = pat.match(line)
        if m:
            rows[(m.group(1), m.group(2))] = dict(units=int(m.group(3)), dev=float(m.group(4)), heldout=float(m.group(5)), exact=float(m.group(6)),
                                                  recall=float(m.group(7)), precision=float(m.group(8)), own=float(m.group(9)),
                                                  multi=float(m.group(10)), empty=float(m.group(11)))
    return rows


def unit_cover_row(r):
    """A result of T.unit_cover at the three decimals the frozen script prints."""
    f = lambda x: float(f"{x:.3f}")
    return dict(units=r["n_units"], dev=f(r["obs"]["dev"][0]), heldout=f(r["obs"]["heldout"][0]), exact=f(r["exact"]), recall=f(r["recall"]),
                precision=f(r["precision"]), own=f(r["clean_ok"]), multi=f(r["clean_multi"]), empty=f(r["clean_empty"]))


def unit_cover_genes(U):
    """The gene table the arms read, from the new loaders (merged sd98 exons, pooled composition, dominant labels, families_cn copy numbers)."""
    return [dict(id=g, S={f for f in U.info[g]["sf"] if not f.startswith("Unassigned")}, comp=Counter(U.strings[g].pooled), dom=set(U.strings[g].u1),
                 cn={"s1c": U.info[g]["cs"], "ours": U.info[g]["co"]}, half=U.info[g]["half"], c=U.info[g]["chrom"], s=U.info[g]["start"], e=U.info[g]["end"])
            for g in U.genes]


# ------------------------------------------------------------------------------------------------------------------------------ Gate 3
def gate3_tables(U, sizes=None, B=500):
    """The deterministic tables compared across PYTHONHASHSEED values: the unit strings, the denominators, the synthetic gate-4 clauses, a made-up
    decisive cell with its z table, the exception lists of a made-up world and a made-up T1-e. No real family label is joined to a real string."""
    sizes = sizes or default_sizes()
    dump = lambda o: json.dumps(o, sort_keys=True, indent=1, default=_jsonable)
    w = T.synthetic_world(11, n_strata=10, mode="independent", sizes=[12] * 10, nclu=(2, 3), vocab="ABC")
    cells = T.run_cells(w["fam"], w["strings"], w["strata"], B=B)
    t1e = t1e_report(synthetic_t1_world(5), B=50)
    return {"strings_table.tsv": strings_table(U), "denominators.json": dump(denominators(U)), "gate4_det.json": dump(gate4_det(B=B, sizes=sizes)),
            "synthetic_cells.json": dump(cells), "synthetic_t1e.json": dump(t1e),
            "exceptions.json": dump({"d1": T.exceptions_d1(w["fam"], w["strings"], w["info"]), "d2": T.exceptions_d2(w["fam"], w["strings"], w["info"])})}


GATE3_FILES = ["strings_table.tsv", "denominators.json", "gate4_det.json", "synthetic_cells.json", "synthetic_t1e.json", "exceptions.json"]
CANARY_WORDS = [f"canary{i:02d}" for i in range(24)]


def hash_canary():
    """The words of a set in iteration order: it differs between PYTHONHASHSEED values, which shows that the variable took effect in a child
    (it is not one of the compared tables)."""
    return " ".join(set(CANARY_WORDS))


# ------------------------------------------------------------------------------------------- layers, sensitivities, sensitivity nulls, sweep
def strata_by(genes, key):
    """Groups of genes by a key, each named by its smallest gene id and ordered by it."""
    by = defaultdict(list)
    for g in genes:
        by[key(g)].append(g)
    strata = [(min(gs), sorted(gs)) for gs in by.values()]
    strata.sort()
    return strata


def layer_world(U, layer, Q):
    """The descriptive layers F_S (Soto clean families on the clean duplicon-bearing genes of V; strata are the F_Q groups, with one stratum 'no group'
    for clean genes in none) and F_Q (the sequence groups over V, ungrouped genes as singleton families)."""
    if layer == "F_S":
        genes = sorted(U.V_clean)
        fam = {g: next(iter(U.info[g]["fams"])) for g in genes}
        strata = strata_by(genes, lambda g: Q.get(g, "no group"))
    elif layer == "F_Q":
        genes = sorted(U.V)
        fam = {g: Q.get(g, "single:" + g) for g in genes}
        strata = None
    else:
        raise ValueError(f"unknown layer {layer!r}")
    return dict(genes=genes, fam=fam, u1={g: U.strings[g].u1 for g in genes}, u0={g: U.strings[g].u0 for g in genes},
                ordered={g: U.strings[g].ordered for g in genes}, info={g: U.info[g] for g in genes}, strata=strata,
                half={g: U.info[g]["half"] for g in genes})


def null_strata(U, w, kind):
    """Strata of the governing null (components of the F_O graph) or of a sensitivity null (the chromosome, or the SD98 region of largest exonic overlap)."""
    if kind == "governing":
        return w["strata"]
    if kind == "chromosome":
        return strata_by(w["genes"], lambda g: U.info[g]["chrom"])
    if kind == "region":
        return strata_by(w["genes"], lambda g: U.region[g] if U.region[g] is not None else "no region")
    raise ValueError(f"unknown null {kind!r}")


def sensitivity_u1(U, genes):
    """The U1 strings under the printed sensitivities: union-per-label bases, numeric ID order for ties, a tie yields the SET of tied labels, the 5% near-tie
    band, and the any-overlap set (every duplicon with at least 1 bp on any exon of the gene)."""
    out = {}
    for name, kw in (("union_per_label", dict(tie="lex", union=True)), ("numeric_ties", dict(tie="num")), ("tie_set", dict(tie="set")),
                     ("near_tie_band", dict(tie="band"))):
        out[name] = {g: T.gene_strings(U.info[g]["exons"], U.info[g]["chrom"], U.colors, **kw).u1 for g in genes}
    out["any_overlap"] = {g: U.strings[g].any_overlap for g in genes}
    return out


def rung_strings(U, w):
    """The strings of rungs U0 (the principal duplicon) and U2 (our famCN `co` and the S1C famCN `cs`) for the genes of a world. The copy-number classes
    are built from every V gene with a famCN and the label in its U1 set; genes without a famCN are outside U2."""
    genes = w["genes"]
    out = dict(u0={g: frozenset({U.strings[g].u0}) for g in genes if U.strings[g].u0})
    u1V = {g: U.strings[g].u1 for g in U.V}
    for src, tag in (("co", "u2co"), ("cs", "u2cs")):
        cn = {g: U.info[g][src] for g in U.V}
        idx = T.u2_index(u1V, cn)
        s2 = T.u2_strings(u1V, cn)
        out[tag] = {g: s2[g] for g in genes if g in s2}
        out["ordered_" + tag] = T.u2_ordered({g: U.strings[g].ordered for g in genes}, cn, idx)
    return out


def sweep_dir(infl):
    return f"{OUT}/sweep_I{infl:.1f}"


def ensure_sweep(infl):
    """The F_O cut at another MCL inflation on the same graph: the frozen binary with the registered flags and --inflation. The pre-MCL graph must come
    out byte for byte as at 2.8 (it does not depend on the inflation)."""
    P = T.default_paths()
    d = sweep_dir(infl)
    os.makedirs(d, exist_ok=True)
    for n in ("ours.clusters.tsv", "ours.loci.tsv", "ours.params.tsv", "ours.graph.tsv", "names.tsv"):
        if os.path.exists(os.path.join(d, n)):
            os.remove(os.path.join(d, n))                    # files this command wrote earlier, by exact path
    cmd = [P["mcl_binary"], "--paf", INPUTS + "/loci.paf", "--gff", INPUTS + "/loci.gff3", "--min-exonic-bp", "1", "--min-shared-exon-frac", "0.60",
           "--inflation", str(infl), "--dump-graph", d + "/ours.graph.tsv", "--out", d + "/ours"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=170)
    if r.returncode != 0:
        sys.exit(f"the sweep run at inflation {infl} failed: {r.stderr[-300:]}")
    if not same_bytes(INPUTS + "/ours.graph.tsv", d + "/ours.graph.tsv"):
        sys.exit(f"the pre-MCL graph at inflation {infl} differs from the registered graph")
    shutil.copy2(INPUTS + "/names.tsv", d + "/names.tsv")
    return d


def sweep_universe(infl):
    d = ensure_sweep(infl)
    P = T.default_paths()
    P["fo_dir"] = d
    return T.build_universe(P)


def singleton_arm(U):
    """The descriptive arm that adds the unclustered V genes as singleton families: m2 over V, the number of unclustered genes and how many of them
    share a U1 string with a clustered gene (each is a D2 violation there)."""
    fam_ext = {g: U.fo_label[g] or "S:" + g for g in U.V}
    clustered = {U.strings[g].u1 for g in U.VF}
    return dict(m2=T.m2(fam_ext, {g: U.strings[g].u1 for g in U.V}, list(U.V)), unclustered=len(U.V) - len(U.VF),
                unclustered_sharing_a_string_with_a_clustered_gene=sum(1 for g in U.V if U.fo_label[g] is None and U.strings[g].u1 in clustered))


def sensitivity_null_report(U, w, kind, B=T.B_DRAWS):
    strata = null_strata(U, w, kind)
    res = T.run_general_cells(w["fam"], w["u1"], strata, B=B)
    return dict(null=kind, n_strata=len(strata), m1c=res["m1c"], m2=res["m2"])


# ------------------------------------------------------------------------------------------------------------------------- Gate 1 commands
def _read_graph_counts(path):
    edges = []
    with open(path) as fh:
        for line in fh:
            u, v, w = line.rstrip("\n").split("\t")
            edges.append((u, v, float(w)))
    comp = T.graph_components(edges)
    return dict(nodes=len({n for e in edges for n in e[:2]}), edges=len(edges), components=len(set(comp.values())))


def gate1_fo():
    """The frozen F_O binary re-run with the registered flags reproduces ours.clusters.tsv, ours.loci.tsv and ours.pairs.tsv byte for byte (and the
    graph dump); the dumped graph has 1,819 nodes, 11,031 edges and 339 components."""
    P = T.default_paths()
    outdir = OUT + "/fo_rerun"
    os.makedirs(outdir, exist_ok=True)
    names = ("ours.clusters.tsv", "ours.loci.tsv", "ours.pairs.tsv", "ours.params.tsv", "ours.graph.tsv")
    for n in names + ("run.log",):
        pth = os.path.join(outdir, n)
        if os.path.exists(pth):
            os.remove(pth)                                  # a file this command wrote in an earlier call, by exact path
    cmd = [P["mcl_binary"], "--paf", INPUTS + "/loci.paf", "--gff", INPUTS + "/loci.gff3", "--min-exonic-bp", "1", "--min-shared-exon-frac", "0.60",
           "--dump-pairs", "--dump-graph", outdir + "/ours.graph.tsv", "--out", outdir + "/ours"]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=170)
    write_text(outdir + "/run.log", f"exit {r.returncode}\n{r.stdout}{r.stderr}")
    same = {n: os.path.exists(os.path.join(outdir, n)) and same_bytes(os.path.join(INPUTS, n), os.path.join(outdir, n)) for n in names}
    counts = _read_graph_counts(INPUTS + "/ours.graph.tsv")
    rerun_counts = _read_graph_counts(outdir + "/ours.graph.tsv") if os.path.exists(outdir + "/ours.graph.tsv") else None
    want = dict(nodes=1819, edges=11031, components=339)
    ok = r.returncode == 0 and all(same[n] for n in ("ours.clusters.tsv", "ours.loci.tsv", "ours.pairs.tsv", "ours.graph.tsv")) and counts == want and rerun_counts == want
    res = dict(gate="gate1_fo", status="PASS" if ok else "FAIL", exit=r.returncode, seconds=round(time.time() - t0, 1), command=" ".join(cmd),
               byte_identical=same, graph_counts_copy=counts, graph_counts_rerun=rerun_counts, registered=want,
               note="ours.params.tsv is compared for information; the registered claim covers clusters, loci and pairs")
    write_json(OUT + "/gate1_fo.json", res)
    text = "Gate 1 (F_O binary re-run), " + res["status"] + f": exit {r.returncode}, {res['seconds']} s\n" + "\n".join(
        f"  {n}: {'identical' if same[n] else 'DIFFERENT'}" for n in names) + f"\n  graph of the copy: {counts}; graph of the re-run: {rerun_counts}; registered {want}"
    write_text(OUT + "/gate1_fo.txt", text)
    say(text)
    return 0 if ok else 1


def _cndup_clusters(P):
    """The clusters of KEY=cnduplicon built the way the frozen script builds them (frozen soto_replication calls, not edited)."""
    sys.path.insert(0, REPO + "/bench/soto")
    import soto_replication as sr
    genes, _ = sr.load_geneset(P["elig"])
    full, _ = sr.load_geneset(P["full"])
    edges = sr.read_edges(P["mapback"])
    clean, ambiguous = sr.load_truth(P["s1c"])
    half = sr.load_split(P["split"])
    members = defaultdict(set)
    for g, f in clean.items():
        if not f.startswith("Unassigned"):
            members[f].add(g)
    soto = {f: m for f, m in members.items() if len(m) >= 2}
    g2f = {g: f for f, m in soto.items() for g in m}
    cover, _kept, leaf_of = sr.pair_families(edges, genes, full, {}, gate=False)
    pred = sr.collapse_cover(cover, leaf_of, genes)
    universe = (full | genes) - ambiguous
    seqfam = defaultdict(set)
    for g in universe:
        if pred.get(g):
            seqfam[pred[g]].add(g)
    clusters = []
    for c in sorted(seqfam):
        fams = sorted({g2f[g] for g in seqfam[c] if g in g2f})
        if len(fams) >= 2:
            gs = sorted(g for f in fams for g in soto[f])
            clusters.append((c, gs, [g2f[g] for g in gs], fams))

    def assign_half(x):
        _, _gs, _labels, fams = x
        big = max(fams, key=lambda f: (len(soto[f]), f))
        hs = [half.get(g) for g in soto[big]]
        return "heldout" if hs.count("heldout") > hs.count("dev") else "dev"
    return clusters, assign_half


def gate1_cndup():
    """KEY=cnduplicon reproduced by the new loaders: the same 33 clusters, 87 Soto families and 642 genes; per-gene compositions equal to the frozen
    compositions; the frozen statistic fed the new compositions at 20 shuffles equals the frozen script's own output, and at 10,000 shuffles
    (seed 20260930) gives the registered pooled deltas."""
    import contextlib
    import io
    sys.path.insert(0, REPO + "/bench/soto_m2")
    sys.path.insert(0, REPO + "/bench/soto")
    import soto_m2_duplicons as fz
    P = T.default_paths()
    t0 = time.time()
    clusters, assign_half = _cndup_clusters(P)
    allg = {g for x in clusters for g in x[1]}
    nfam = sum(len(x[3]) for x in clusters)
    U = load_universe()
    mine = {g: dict(U.strings[g].pooled) for g in allg}
    col, maxlen = fz.load_duplicons(P["colors"])
    frozen = fz.compositions(P["sd98_exons"], col, maxlen, allg)
    comp_equal = mine == frozen
    simc = lambda x, y: fz.wjac(mine[x], mine[y])
    subsets = {"all": clusters, "dev": [x for x in clusters if assign_half(x) == "dev"], "heldout": [x for x in clusters if assign_half(x) == "heldout"]}

    def run(name, N):
        fz.N_SHUFFLE = N
        per, pooled, p, used = fz.run_test([(gs, lab) for (_, gs, lab, _) in subsets[name]], simc, random.Random(fz.SEED))
        v, frac = fz.verdict(per, p)
        return dict(pooled_delta=pooled, p=p, clusters_used=used, share_positive=frac, verdict=v)

    # the frozen script's own printout at 20 shuffles against the same path with the new compositions
    buf = io.StringIO()
    fz.N_SHUFFLE = 20
    with contextlib.redirect_stdout(buf):
        fz.main(["--geneset", P["elig"], "--full-geneset", P["full"], "--exons", P["sd98_exons"], "--duplicons", P["colors"]])
    pat = re.compile(r"### (PRIMARY[^\n]*)\n\nclusters used (\d+); pooled delta ([+-]\d+\.\d+); pooled p (\d+\.\d+); clusters with delta > 0: (\d+\.\d+)")
    printed = pat.findall(buf.getvalue())
    path20 = [run(n, 20) for n in ("all", "dev", "heldout")]
    eq20 = len(printed) == 3 and all(f"{r['pooled_delta']:+.4f}" == pr[2] and f"{r['p']:.4f}" == pr[3] and f"{r['clusters_used']}" == pr[1] and f"{r['share_positive']:.2f}" == pr[4]
                                      for r, pr in zip(path20, printed))
    big = {n: run(n, 10000) for n in ("dev", "heldout", "all")}
    reg = CNDUP_REGISTERED
    checks = dict(clusters=len(clusters) == reg["clusters"], soto_families=nfam == reg["soto_families"], genes=len(allg) == reg["genes"],
                  compositions_equal=comp_equal, frozen_printout_equal_at_20=eq20,
                  dev_delta=f"{big['dev']['pooled_delta']:.4f}" == f"{reg['dev_delta']:.4f}",
                  heldout_delta=f"{big['heldout']['pooled_delta']:.4f}" == f"{reg['heldout_delta']:.4f}")
    ok = all(checks.values())
    res = dict(gate="gate1_cnduplicon", status="PASS" if ok else "FAIL", checks=checks, clusters=len(clusters), soto_families=nfam, genes=len(allg),
               genes_without_exonic_duplicon=sum(1 for g in allg if not mine[g]), n_shuffle_10000=big, frozen_printout_at_20=printed, new_path_at_20=path20,
               seconds=round(time.time() - t0, 1), registered=reg)
    write_json(OUT + "/gate1_cnduplicon.json", res)
    text = (f"Gate 1 (KEY=cnduplicon), {res['status']}: {len(clusters)} clusters, {nfam} Soto families, {len(allg)} genes; compositions equal to the frozen ones: {comp_equal}; "
            f"frozen printout at 20 shuffles equal: {eq20}\n"
            + "\n".join(f"  {n}: pooled delta {big[n]['pooled_delta']:+.4f}, pooled p {big[n]['p']:.4f} ({big[n]['clusters_used']} clusters, share positive {big[n]['share_positive']:.2f}, {big[n]['verdict']}) at 10,000 shuffles"
                         for n in ("dev", "heldout", "all"))
            + f"\n  registered: dev {reg['dev_delta']:+.4f}, held-out {reg['heldout_delta']:+.4f}; checks {checks}")
    write_text(OUT + "/gate1_cnduplicon.txt", text)
    say(text)
    return 0 if ok else 1


def gate1_unitcover():
    """KEY=unitcovercn reproduced by the new loaders: arms A, B, C and D (both famCN sources) computed from the new gene table against the output of
    the frozen script, and arm A and the six unit counts against the registered values."""
    P = T.default_paths()
    t0 = time.time()
    U = load_universe()
    genes = unit_cover_genes(U)
    mine = {}
    for arm, src in (("A", "s1c"), ("B", "s1c"), ("C", "s1c"), ("C", "ours"), ("D", "s1c"), ("D", "ours")):
        mine[(arm, "-" if arm in ("A", "B") else src)] = unit_cover_row(T.unit_cover(genes, arm, src))
    r = subprocess.run([sys.executable, "-B", FROZEN_UNIT, "--data", P["families_cn"], "--duplicons", P["colors"]], capture_output=True, text=True, timeout=170)
    frozen = parse_unitcover_table(r.stdout)
    equal = {f"{k[0]}/{k[1]}": mine[k] == frozen.get(k) for k in mine}
    arm_a = mine[("A", "-")]
    a_ok = all(arm_a[k] == UNITCOVER_A[k] for k in UNITCOVER_A)
    units_ok = {f"{k[0]}/{k[1]}": mine[k]["units"] == v for k, v in UNITCOVER_UNITS.items()}
    ok = r.returncode == 0 and all(equal.values()) and a_ok and all(units_ok.values()) and bool(U.exons_match_families_cn)
    res = dict(gate="gate1_unitcover", status="PASS" if ok else "FAIL", new_loaders={f"{k[0]}/{k[1]}": v for k, v in mine.items()},
               frozen_script={f"{k[0]}/{k[1]}": v for k, v in frozen.items()}, equal_to_frozen_script=equal, arm_A_equals_registered=a_ok, registered_arm_A=UNITCOVER_A,
               unit_counts_equal_registered=units_ok, registered_unit_counts={f"{k[0]}/{k[1]}": v for k, v in UNITCOVER_UNITS.items()},
               exons_match_families_cn=bool(U.exons_match_families_cn), seconds=round(time.time() - t0, 1))
    write_json(OUT + "/gate1_unitcover.json", res)
    lines = [f"Gate 1 (KEY=unitcovercn), {res['status']}: arms from the new loaders against the frozen script's table and the registered values"]
    for k in sorted(mine):
        m = mine[k]
        lines.append(f"  {k[0]} {k[1]:<4} units {m['units']} (registered {UNITCOVER_UNITS[k]}), Jaccard {m['dev']:.3f} / {m['heldout']:.3f}, exact {m['exact']:.3f}, recall {m['recall']:.3f}, "
                     f"precision {m['precision']:.3f}, clean {m['own']:.3f} / {m['multi']:.3f} / {m['empty']:.3f}; equal to the frozen table: {equal[f'{k[0]}/{k[1]}']}")
    lines.append(f"  arm A equals the registered values {UNITCOVER_A}: {a_ok}; sd98 exons equal families_cn exons for every gene: {U.exons_match_families_cn}")
    write_text(OUT + "/gate1_unitcover.txt", "\n".join(lines))
    say("\n".join(lines))
    return 0 if ok else 1


def gate1_collect():
    parts = {n: gate_status(OUT, n) for n in ("gate1_fo", "gate1_cnduplicon", "gate1_unitcover")}
    ok = all(v == "PASS" for v in parts.values())
    write_json(OUT + "/gate1.json", dict(gate="gate1", status="PASS" if ok else "FAIL", parts=parts))
    text = f"Gate 1: {'PASS' if ok else 'FAIL'}  {parts}"
    write_text(OUT + "/gate1.txt", text)
    say(text)
    return 0 if ok else 1


def cmd_gate1(a):
    return {"fo": gate1_fo, "cnduplicon": gate1_cndup, "unitcover": gate1_unitcover, "collect": gate1_collect}[a.part]()


def cmd_gate2(a=None):
    """Gate 2 covers Phase R3a and the optional matched human run, not Phase R1; it needs exon sequences from two assemblies and is not run here."""
    res = dict(gate="gate2", status="SKIPPED", reason="Gate 2 covers Phase R3a and the optional matched human run only (section 3); Phase R1 does not need it. "
               "Sample rule as committed: chr1..chr22 in order, min(20, n) genes per chromosome by one random.Random(20260930).sample over the sorted ids of the genes whose v1.0 CAT v4 coordinates lie on it.")
    write_json(OUT + "/gate2.json", res)
    write_text(OUT + "/gate2.txt", "Gate 2: SKIPPED (not required for Phase R1)")
    say("Gate 2: SKIPPED (not required for Phase R1)")
    return 0


# ------------------------------------------------------------------------------------------------------------------------- Gate 3 command
def cmd_gate3_worker(a):
    os.makedirs(a.out, exist_ok=True)
    U = load_universe()
    sizes = T.cluster_sizes(_vf_fam(U), U.strata)
    for name, text in gate3_tables(U, sizes=sizes, B=500).items():
        write_text(os.path.join(a.out, name), text)
    write_text(os.path.join(a.out, "canary.txt"), hash_canary())
    write_json(os.path.join(a.out, "meta.json"), dict(PYTHONHASHSEED=os.environ.get("PYTHONHASHSEED"), python=sys.version.split()[0], bytecode_off=bool(sys.flags.dont_write_bytecode)))
    return 0


def cmd_gate3(a=None):
    """Runs under `env PYTHONHASHSEED=0 python3 -B` and `env PYTHONHASHSEED=1 python3 -B` (never -I) and compares the tables byte for byte."""
    runs = {}
    t0 = time.time()
    for hs in ("0", "1"):
        out = f"{OUT}/gate3_hs{hs}"
        env = dict(os.environ, PYTHONHASHSEED=hs)
        r = subprocess.run([sys.executable, "-B", os.path.abspath(__file__), "gate3-worker", "--out", out], env=env, capture_output=True, text=True, timeout=170)
        runs[hs] = dict(exit=r.returncode, out=out, stderr=r.stderr[-400:])
    differing = compare_dirs(runs["0"]["out"], runs["1"]["out"], GATE3_FILES)
    ran = all(r["exit"] == 0 for r in runs.values())
    metas = {hs: json.load(open(f"{runs[hs]['out']}/meta.json")) if os.path.exists(f"{runs[hs]['out']}/meta.json") else None for hs in runs}
    seeds_ok = all(metas[hs] and metas[hs]["PYTHONHASHSEED"] == hs and metas[hs]["bytecode_off"] for hs in runs)
    canary_differs = bool(compare_dirs(runs["0"]["out"], runs["1"]["out"], ["canary.txt"]))      # the seed took effect, so the comparison is not vacuous
    ok = ran and not differing and seeds_ok and canary_differs
    res = dict(gate="gate3", status="PASS" if ok else "FAIL", runs=runs, differing=differing, files=GATE3_FILES, metas=metas, canary_differs=canary_differs,
               sha256={n: T.sha256_file(f"{runs['0']['out']}/{n}") for n in GATE3_FILES if os.path.exists(f"{runs['0']['out']}/{n}")}, seconds=round(time.time() - t0, 1))
    write_json(OUT + "/gate3.json", res)
    text = (f"Gate 3 (determinism), {res['status']}: PYTHONHASHSEED=0 and =1 under python3 -B; files compared: {', '.join(GATE3_FILES)}; differing: {differing or 'none'}; "
            f"canary (a set of strings in iteration order) differs between the seeds: {canary_differs}; {res['seconds']} s")
    write_text(OUT + "/gate3.txt", text)
    say(text)
    return 0 if ok else 1


# ------------------------------------------------------------------------------------------------------------------------- Gate 4 command
def _real_sizes():
    U = load_universe()
    return T.cluster_sizes(_vf_fam(U), U.strata)


def cmd_gate4(a):
    if a.part == "det":
        t0 = time.time()
        sizes = _real_sizes()
        rows = gate4_det(B=a.B, sizes=sizes)
        ok = all(r["ok"] for r in rows)
        res = dict(gate="gate4_det", status="PASS" if ok else "FAIL", B=a.B, strata=len(sizes), genes=sum(map(sum, sizes)), rows=rows, seconds=round(time.time() - t0, 1))
        write_json(OUT + "/gate4_det.json", res)
        text = f"Gate 4, deterministic clauses (B = {a.B}, real cluster-size shape of {len(sizes)} strata, synthetic strings), {res['status']}, {res['seconds']} s\n" + "\n".join(
            f"  {'ok' if r['ok'] else 'FAIL':4} {r['name']}\n       {r['detail']}" if r["detail"] else f"  {'ok' if r['ok'] else 'FAIL':4} {r['name']}" for r in rows)
        write_text(OUT + "/gate4_det.txt", text)
        say(text)
        return 0 if ok else 1
    if a.part == "calib":
        t0 = time.time()
        sizes = _real_sizes()
        res = gate4_calib_block(a.start, a.end, B=a.B, nrep=100, sizes=sizes)
        res["seconds"] = round(time.time() - t0, 1)
        write_json(f"{OUT}/gate4_calib_{a.start:03d}_{a.end:03d}.json", res)
        say(f"Gate 4 calibration repeats {a.start}..{a.end - 1} (B = {a.B}): p <= 0.05 in m1c {sum(p <= 0.05 for p in res['m1c'])}, in m2 {sum(p <= 0.05 for p in res['m2'])} of {a.end - a.start}; {res['seconds']} s")
        return 0
    blocks = []
    for n in sorted(os.listdir(OUT)):
        if re.match(r"gate4_calib_\d{3}_\d{3}\.json$", n):
            with open(os.path.join(OUT, n)) as fh:
                blocks.append(json.load(fh))
    coll = gate4_collect(blocks)
    det = gate_status(OUT, "gate4_det")
    ok = coll["ok"] and det == "PASS"
    res = dict(gate="gate4", status="PASS" if ok else "FAIL", deterministic_clauses=det, calibration=coll)
    write_json(OUT + "/gate4.json", res)
    text = (f"Gate 4, {res['status']}: deterministic clauses {det}; calibration of planted random partitions: {coll['repeats']} repeats "
            f"(complete: {coll['complete']}), p <= 0.05 in {coll['small_p_counts']} (limit {coll['limit']} of 100 per cell), mean p {coll['mean_p']}")
    write_text(OUT + "/gate4.txt", text)
    say(text)
    return 0 if ok else 1


# ------------------------------------------------------------------------------------------------------------------------- Gate 5 command
def gate5a_real(U, Q):
    """Gate 5a on the real layers (F_O, F_Q, F_S). It joins real family labels to real unit strings, so it is not run before the independent code review."""
    rows = []
    clean = {g: next(iter(U.info[g]["fams"])) for g in U.genes if U.info[g]["cls"] == "clean"}
    layers = {"F_O": {g: U.fo_label[g] for g in U.VF}, "F_Q": {g: Q.get(g, "single:" + g) for g in U.V}, "F_S": {g: clean[g] for g in U.V_clean}}
    for src in ("cs", "co"):
        u1 = {g: U.strings[g].u1 for g in U.V if U.info[g][src] is not None}
        s2 = T.u2_strings(u1, {g: U.info[g][src] for g in u1})
        for name, fam in layers.items():
            cu = [g for g in u1 if g in fam]
            f1 = {g: u1[g] for g in cu}
            f2 = {g: s2[g] for g in cu if g in s2}
            cu = [g for g in cu if g in f2]
            a = (T.m1(fam, f1, cu), T.m1c(fam, f1, cu), T.m2(fam, f1, cu))
            b = (T.m1(fam, f2, cu), T.m1c(fam, f2, cu), T.m2(fam, f2, cu))
            rows.append(dict(source=src, layer=name, genes=len(cu), u1=a, u2=b, ok=bool(a[0] <= b[0] and a[1] <= b[1] and a[2] >= b[2])))
    return rows


def cmd_gate5(a):
    t0 = time.time()
    P = T.default_paths()
    U = load_universe()
    Q, _universe = fq_partition()
    clean = soto_clean_from_families_cn(P["families_cn"])
    v = gate5b_values(clean, Q)
    rows5b = compare_registered({k: x[0] for k, x in REGISTERED_5B.items()}, v)
    ok5b = all_ok(rows5b)
    rows5a = []
    for src in ("cs", "co"):
        rows5a += gate5a_surrogate({g: U.strings[g].u1 for g in U.V}, {g: U.info[g][src] for g in U.V}, seeds=(0, 1, 2), source=src)
    ok5a = all(r["ok"] for r in rows5a) and all(r["u2_refines_u1"] for r in rows5a)
    real = None
    if a.real_5a:
        guard_real_cell()
        real = gate5a_real(U, Q)
        ok5a = ok5a and all(r["ok"] for r in real)
    ok = ok5a and ok5b
    res = dict(gate="gate5", status="PASS" if ok else "FAIL", gate5b=dict(values=v, tripwires=rows5b, ok=ok5b), gate5a_surrogate=dict(rows=rows5a, ok=ok5a),
               gate5a_real=real if real is not None else "NOT RUN: it joins real family labels to real unit strings (hard stop)", seconds=round(time.time() - t0, 1))
    write_json(OUT + "/gate5.json", res)
    lines = [f"Gate 5, {res['status']} ({res['seconds']} s)",
             f"  5b: scorer returns D[S -> Q] = {v['D_S_to_Q_singleton']} genes ({v['S_split_singleton']} of {v['S_scored_singleton']} Soto families split) with ungrouped genes as singleton cells and "
             f"{v['D_S_to_Q_ignored']} genes ({v['S_split_ignored']} families split) when they are ignored; registered 7 (4 of 444) and 2: {'ok' if ok5b else 'MISMATCH'}",
             f"  5a (surrogate partitions on the real common universe, U1 -> U2): {len(rows5a)} checks over {sorted({r['source'] for r in rows5a})}, all ok: {ok5a}; "
             f"genes in the common universe {sorted({(r['source'], r['genes']) for r in rows5a})}; U2 refines U1: {all(r['u2_refines_u1'] for r in rows5a)}",
             "  5a on the real layers: " + ("run" if real is not None else "not run (hard stop)")]
    write_text(OUT + "/gate5.txt", "\n".join(lines))
    say("\n".join(lines))
    return 0 if ok else 1


# ------------------------------------------------------------------------------------------------------------------------- report builders
def _relabel_non_primary(rep):
    rep["verdict_scope"] = "descriptive"
    if rep["verdict"] == "NOT-REFUTED":
        rep["verdict"] = "counts are 0 (NOT-REFUTED is issued only on the primary arm)"
    return rep


def build_t1a(U, Q):
    """T1-a and its descriptive companions on the loaded universe: the primary arm with the full exception lists and the audit, the exclusion arms, the
    halves, the singleton arm, rungs U0 and U2, the string sensitivities and the descriptive layers F_Q and F_S. Returns (result, text)."""
    res = dict(arms={}, halves={}, layers={})
    for arm, dele in arm_deletions(U).items():
        rep = t1a_report(world_of_universe(U, delete=dele))
        rep.update(arm=arm, deleted=len(dele), denominator_guard=denominator_guard(len(dele), len(U.V)))
        res["arms"][arm] = rep if arm == "none" else _relabel_non_primary(rep)
    for h in ("dev", "heldout"):
        res["halves"][h] = _relabel_non_primary(t1a_report(world_of_universe(U, half=h)))
    res["singleton_arm"] = singleton_arm(U)
    w = world_of_universe(U)
    genes, fam = w["genes"], w["fam"]
    rs = rung_strings(U, w)
    rungs = {"U0": _stats(fam, rs["u0"], None, [g for g in genes if g in rs["u0"]])}
    for tag in ("u2co", "u2cs"):
        gg = [g for g in genes if g in rs[tag]]
        rungs[tag.upper()] = _stats(fam, rs[tag], rs["ordered_" + tag], gg)
        rungs[tag.upper()]["genes_without_famcn"] = len(genes) - len(gg)
    res["rungs"] = rungs
    res["sensitivities"] = {k: dict(m1c=T.m1c(fam, v, genes), m2=T.m2(fam, v, genes)) for k, v in sensitivity_u1(U, genes).items()}
    for layer in ("F_Q", "F_S"):
        lw = layer_world(U, layer, Q)
        rl = rung_strings(U, lw)
        entry = dict(genes=len(lw["genes"]), U1=_stats(lw["fam"], lw["u1"], lw["ordered"], lw["genes"]),
                     U0=_stats(lw["fam"], rl["u0"], None, [g for g in lw["genes"] if g in rl["u0"]]))
        for tag in ("u2co", "u2cs"):
            gg = [g for g in lw["genes"] if g in rl[tag]]
            entry[tag.upper()] = _stats(lw["fam"], rl[tag], rl["ordered_" + tag], gg)
        res["layers"][layer] = entry
    p = res["arms"]["none"]
    lines = [f"T1-a (F_O, rung U1, V_F = {p['n']} genes): m1c (D1) = {p['m1c']}, m2 (D2) = {p['m2']} -> {p['verdict']}",
             f"  audit labels D1 {p['labels_d1']} (m1c without short {p['m1c_without_short']}); D2 {p['labels_d2']} (m2 without short {p['m2_without_short']}); "
             f"m1 {p['m1']}, m1t {p['m1t']}, macro D1 {p['macro_d1']:.4f}, macro D2 {p['macro_d2']:.4f}, path defect {p['delta']}, multiplicity {p['mult']}",
             f"  size-2 stratum (n_f = 2): {p['size2']['n']} genes, m1c {p['size2']['m1c']}, m2 {p['size2']['m2']}, macro D1 {p['size2']['macro_d1']:.4f}, macro D2 {p['size2']['macro_d2']:.4f}"]
    lines += [f"  arm {k}: n {v['n']}, deleted {v['deleted']}{' (denominator guard)' if v['denominator_guard'] else ''}, m1c {v['m1c']}, m2 {v['m2']}" for k, v in res["arms"].items()]
    lines += [f"  half {k}: n {v['n']}, m1c {v['m1c']}, m2 {v['m2']}" for k, v in res["halves"].items()]
    s = res["singleton_arm"]
    lines.append(f"  singleton arm: m2 {s['m2']} over V; {s['unclustered']} unclustered genes, {s['unclustered_sharing_a_string_with_a_clustered_gene']} share a U1 string with a clustered gene")
    lines += [f"  rung {k}: n {v['n']}, m1 {v['m1']}, m1c {v['m1c']}, m2 {v['m2']}" + (f", without famCN {v['genes_without_famcn']}" if "genes_without_famcn" in v else "") for k, v in rungs.items()]
    lines += [f"  string sensitivity {k}: m1c {v['m1c']}, m2 {v['m2']}" for k, v in res["sensitivities"].items()]
    lines += [f"  layer {k} (descriptive): U1 m1c {v['U1']['m1c']}, m2 {v['U1']['m2']} on {v['genes']} genes; U0 m1 {v['U0']['m1']}, m2 {v['U0']['m2']}" for k, v in res["layers"].items()]
    return res, "\n".join(lines)


def level_specific_tag(main_cls, sweep_cls):
    """The headline tag of the level sweep: LEVEL-SPECIFIC when the class at another inflation differs from the class at 2.8; a sweep cell that is
    UNDERPOWERED (or UNINFORMATIVE) carries that tag and no comparison."""
    if sweep_cls in ("UNDERPOWERED", "UNINFORMATIVE"):
        return f"{sweep_cls} (no LEVEL-SPECIFIC comparison)"
    if main_cls in T.CLASS_ORDER and sweep_cls in T.CLASS_ORDER and main_cls != sweep_cls:
        return "LEVEL-SPECIFIC"
    return ""


def build_t1bc(U, stat, B=T.B_DRAWS, null="governing", layer="F_O", Q=None):
    """T1-b (stat m1c) or T1-c (stat m2). The governing null on F_O gives the cell with its class; a sensitivity null (chromosome, SD98 region) and the
    descriptive layer F_S (strata are the F_Q groups) give p and the tails and no class."""
    cell = {"m1c": "T1-b (D1, m1c)", "m2": "T1-c (D2, m2)"}[stat]
    if layer == "F_S":
        lw = layer_world(U, "F_S", Q or {})
        r = T.run_general_cells(lw["fam"], lw["u1"], lw["strata"], B=B)
        out = dict(layer="F_S", null="F_Q groups", n_strata=len(lw["strata"]), m1c=r["m1c"], m2=r["m2"])
        c = r[stat]
        return out, (f"{cell}, F_S descriptive (strata are the F_Q groups plus one stratum 'no group'; {c['n_strata']} strata, {c['n_var']} variation-bearing; no class): "
                     f"observed {c['obs']}, p = {c['p']:.6f}, upper tail {c['p_upper']:.6f}, draws mean {c['draw_mean']:.2f} (min {c['draw_min']}, max {c['draw_max']})")
    w = world_of_universe(U)
    if null != "governing":
        out = sensitivity_null_report(U, w, null, B)
        c = out[stat]
        return out, (f"{cell}, sensitivity null '{null}' ({out['n_strata']} strata, {c['n_var']} with two or more labels; no class; p <= 0.01 is expected whatever the "
                     f"organisation, prediction Q4): observed {c['obs']}, p = {c['p']:.6f}, upper tail {c['p_upper']:.6f}, draws mean {c['draw_mean']:.2f} (min {c['draw_min']}, max {c['draw_max']})")
    c = t1bc_report(w, stat, B=B)
    return c, format_cell(c)


def build_t1e(U, B=2000):
    s = t1e_report(world_of_universe(U), B=B)
    return s, format_t1e(s)


# ------------------------------------------------------------------------------------------------------------------------- the cells
def _write_cell(name, res, text):
    write_json(f"{OUT}/{name}.json", res)
    write_text(f"{OUT}/{name}.txt", text)
    say(text)


def cmd_t1a(a):
    if a.synthetic:
        r = t1a_report(synthetic_t1_world(a.seed))
        _write_cell("synthetic_t1a", r, f"T1-a on a made-up world: m1c = {r['m1c']}, m2 = {r['m2']} on {r['n']} genes -> {r['verdict']}")
        return 0
    guard_real_cell()
    U = load_universe()
    Q, _universe = fq_partition()
    res, text = build_t1a(U, Q)
    res["prereg"] = prereg_stamp()
    _write_cell("t1a", res, text)
    return 0


def cmd_t1bc(a, stat):
    name = {"m1c": "t1b", "m2": "t1c"}[stat]
    if a.synthetic:
        w = synthetic_t1_world(a.seed)
        cell = t1bc_report(w, stat, B=a.B, seed=a.seed)
        _write_cell("synthetic_" + name, cell, format_cell(cell))
        return 0
    guard_real_cell()
    infl = a.inflation
    U = sweep_universe(infl) if infl else load_universe()
    Q = fq_partition()[0] if a.layer == "F_S" else None
    res, text = build_t1bc(U, stat, B=a.B, null=a.null, layer=a.layer, Q=Q)
    suffix = ""
    if a.layer == "F_S":
        suffix = "_F_S"
    elif a.null != "governing":
        suffix = "_" + a.null
    elif infl:
        suffix = f"_I{infl:.1f}"
        main = f"{OUT}/{name}.json"
        if os.path.exists(main):
            with open(main) as fh:
                tag = level_specific_tag(json.load(fh)["class"], res["class"])
            res["level_tag"] = tag
            text += f"\nlevel sweep at inflation {infl}: class {res['class']}; {tag or 'same class as at 2.8'} (reported, not acted on)"
    res["prereg"] = prereg_stamp()
    _write_cell(name + suffix, res, text)
    return 0


def cmd_t1e(a):
    if a.synthetic:
        s = t1e_report(synthetic_t1_world(a.seed), B=min(a.B, 2000))
        _write_cell("synthetic_t1e", s, format_t1e(s))
        return 0
    guard_real_cell()
    s, text = build_t1e(load_universe(), B=2000)
    s["prereg"] = prereg_stamp()
    _write_cell("t1e", s, text)
    return 0


# ------------------------------------------------------------------------------------------------------------------------------ summary of R0
def r0_summary(out_dir=OUT):
    """One line per gate of Phase R0 and whether the gates that cover Phase R1 all passed."""
    lines, ok = [], True
    for g in ("gate0", "gate1", "gate2", "gate3", "gate4", "gate5", "gate6"):
        st = gate_status(out_dir, g)
        lines.append(f"{g}: {st}")
        if g in GATES_FOR_R1:
            ok &= st == "PASS"
    lines.append("every gate that covers Phase R1 passed: " + ("YES" if ok else "NO"))
    return "\n".join(lines), ok


def cmd_r0(a=None):
    text, ok = r0_summary()
    write_text(OUT + "/R0_SUMMARY.txt", text)
    say(text)
    return 0 if ok else 1


# ------------------------------------------------------------------------------------------------------------------------------ main
def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    for c in ("inputs", "hash-paf", "gate0", "gate2", "gate3", "gate6", "r0"):
        sub.add_parser(c)
    p = sub.add_parser("gate1")
    p.add_argument("--part", choices=["fo", "cnduplicon", "unitcover", "collect"], default="collect")
    p = sub.add_parser("gate3-worker")
    p.add_argument("--out", required=True)
    p = sub.add_parser("gate4")
    p.add_argument("--part", choices=["det", "calib", "collect"], default="det")
    p.add_argument("--start", type=int, default=0)
    p.add_argument("--end", type=int, default=20)
    p.add_argument("--B", type=int, default=T.B_DRAWS)
    p = sub.add_parser("gate5")
    p.add_argument("--real-5a", dest="real_5a", action="store_true", help="join real family labels to real strings (needs GO_R1; not run before the review)")
    for c in ("t1a", "t1b", "t1c", "t1e"):
        p = sub.add_parser(c)
        p.add_argument("--synthetic", action="store_true", help="run on a made-up world (the only mode allowed before the independent code review)")
        p.add_argument("--B", type=int, default=T.B_DRAWS)
        p.add_argument("--seed", type=int, default=0)
        if c in ("t1b", "t1c"):
            p.add_argument("--null", choices=["governing", "chromosome", "region"], default="governing", help="the governing null or a sensitivity null (no class)")
            p.add_argument("--layer", choices=["F_O", "F_S"], default="F_O", help="F_S is the descriptive layer with the F_Q groups as strata (no class)")
            p.add_argument("--inflation", type=float, default=None, help="level sweep: the F_O cut at another MCL inflation on the same graph (2.0 or 4.0)")
    return ap


COMMANDS = {"inputs": lambda a: cmd_inputs(), "hash-paf": lambda a: cmd_hash_paf(), "gate0": lambda a: cmd_gate0(), "gate1": cmd_gate1, "gate2": cmd_gate2,
            "gate3": cmd_gate3, "gate3-worker": cmd_gate3_worker, "gate4": cmd_gate4, "gate5": cmd_gate5, "gate6": lambda a: cmd_gate6(), "r0": cmd_r0,
            "t1a": cmd_t1a, "t1b": lambda a: cmd_t1bc(a, "m1c"), "t1c": lambda a: cmd_t1bc(a, "m2"), "t1e": cmd_t1e}


def main(argv=None):
    a = build_parser().parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    return COMMANDS[a.command](a)


if __name__ == "__main__":
    sys.exit(main())
