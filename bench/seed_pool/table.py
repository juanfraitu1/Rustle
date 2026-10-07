#!/usr/bin/env python3
"""The comparison tables of docs/PREREG_seed_pool_real_reads_2026-10-07.md, computed from the products of run.sh score (nothing by hand).

    table.py --work W/SUB --family NPIP --truths compara,u2,soto [--found tc|ann] [--targets "compara:CF153;u2:ID_154,ID_149;soto:ID_154"]

Prints (and writes W/SUB/table.txt, table.tsv, table.md, target.md, cost.md, matrix.tsv): (1) one row per arm: cost, the family's clusters, M1-M3 (M2n =
M2 with the exact locus itself a node), the family scores of the whole contig (bipartite and pairwise); (2) the rows of the TARGET families of each truth;
(3) the per-copy matrix (own node N, found within the own node F, locus-level L); (4) the non-dominated arms on (M1, M2, NP); (5) the registered
predictions S1-S6, S8, each HELD / FAILED / NO VERDICT (UNDERPOWERED when fewer than 8 copies are expressed). A banner says whether the gates
(gates.json) pass and are current.
M2 as the instrument computes it: the copy has an own node AND some same-strand locus overlapping it (a node or not) has the exact representative.
"""
import argparse
import csv
import json
import os
import re

ORDER = ["P", "G100", "G995", "G98", "G95", "G90", "A", "P_R1", "G100_R1", "G995_R1", "G98_R1", "G95_R1", "G90_R1", "A_R1"]
PRIMARY = ["P", "G98", "A", "G98_R1", "A_R1"]


def load(p):
    with open(p) as fh:
        return json.load(fh)


def wall(stderr_path):
    """'Elapsed (wall clock) time (h:mm:ss or m:ss): 1:05.12' of a /usr/bin/time -v log -> seconds, and the peak RSS in GB."""
    try:
        with open(stderr_path) as fh:
            txt = fh.read()
    except OSError:
        return None, None
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", txt)
    r = re.search(r"Maximum resident set size \(kbytes\): (\d+)", txt)
    secs = None
    if m:
        secs = 0.0
        for part in m.group(1).split(":"):
            secs = secs * 60 + float(part)
    return secs, (int(r.group(1)) / 1e6 if r else None)


def count_rows(path, kind):
    n = 0
    try:
        with open(path) as fh:
            for ln in fh:
                f = ln.split("\t")
                if len(f) > 2 and f[2] == kind:
                    n += 1
    except OSError:
        return None
    return n


def dominated(v, w):
    """w dominates v: >= on every coordinate and > on one."""
    return all(b >= a for a, b in zip(v, w)) and any(b > a for a, b in zip(v, w))


def dthousandths(a, b):
    """|a - b| in whole thousandths for two scores printed to three decimals. The registered S4 clause (|dF| < .01) is on these:
    abs(0.011 - 0.001) < 0.01 is True in floating point although the printed scores differ by exactly .010."""
    return abs(round(a * 1000) - round(b * 1000))


def parse_targets(text):
    """'compara:CF153;u2:ID_154,ID_149' -> {'compara': ['CF153'], 'u2': ['ID_154', 'ID_149']}."""
    out = {}
    for part in (text or "").split(";"):
        if ":" in part:
            truth, ids = part.split(":", 1)
            out[truth.strip()] = [i.strip() for i in ids.split(",") if i.strip()]
    return out


def target_rows(fs, targets):
    """The per-family rows of `fs` (a score_families.py JSON) for the requested families, in the requested order; an absent family keeps its id with None values."""
    rows = []
    for truth, ids in targets.items():
        have = {r["family_id"]: r for r in fs.get("truths", {}).get(truth, {}).get("per_family", [])}
        for fid in ids:
            r = have.get(fid)
            if r is None:
                rows.append(dict(truth=truth, family_id=fid, n_truth=None, hit=None, sens=None, prec=None, f=None))
            else:
                rows.append(dict(truth=truth, family_id=fid, n_truth=int(r["n_truth"]), hit=int(r["hit"]), sens=float(r["sens"]),
                                 prec=float(r["prec"]), f=float(r["f"])))
    return rows


def gate_banner(S):
    """One line: whether gates.json exists, passes and is newer than the last scoring (comp.json)."""
    g = os.path.join(S, "gates.json")
    if not os.path.exists(g):
        return "GATES: not run (bench/seed_pool/run.sh gates): the tables below are not certified"
    gates = load(g)
    bad = sorted(k for k, v in gates.items() if v.get("ok") is False)
    if bad:
        return "GATES: FAILED: " + ", ".join(bad) + ": the tables below are INVALID"
    # `score` writes comp.json first, then the support and family-score files: the gates are current only if newer than all of them
    products = [os.path.join(S, n) for n in os.listdir(S) if n in ("comp.json", "nodes.json") or n.startswith(("support", "fs_"))]
    newest = max((os.path.getmtime(p) for p in products if os.path.isfile(p)), default=None)
    if newest is not None and os.path.getmtime(g) < newest:
        return "GATES: older than the last scoring (comp.json, nodes.json, support*, fs_*): re-run bench/seed_pool/run.sh gates"
    unrun = sorted(k for k, v in gates.items() if v.get("ok") is None and "not run" in (v.get("note") or ""))
    na = sorted(k for k, v in gates.items() if v.get("ok") is None and k not in unrun)
    extra = [f"not applicable: {', '.join(na)}"] if na else []
    extra += [f"not run: {', '.join(unrun)}"] if unrun else []
    return "GATES: all that apply and ran pass" + (f" ({'; '.join(extra)})" if extra else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--truths", default="", help="comma list of family truths with a fs_<arm>.json; empty = none on this substrate (M4 not scored)")
    ap.add_argument("--found", default="tc", choices=["tc", "ann"],
                    help="the found criterion: tc = Amendment E (capped start + first three introns; human reads), ann = Amendment A (gorilla: no cap signal)")
    ap.add_argument("--targets", default="", help="target families per truth, e.g. 'compara:CF153;u2:ID_154,ID_149;soto:ID_154' (rows of fs_<arm>.json)")
    a = ap.parse_args()
    S, truths = a.work, [t for t in a.truths.split(",") if t]
    K = a.found                                           # column prefix in copy_support's output
    M2_KEY, M3_KEY, E_KEY = f"{K}_found_in_npip_nodes", f"locus_{K}_found", f"{K}_expressed"
    targets = parse_targets(a.targets)
    comp, sup = load(f"{S}/comp.json"), load(f"{S}/support.json")
    sup_nodes = load(f"{S}/support_nodeonly.json") if os.path.exists(f"{S}/support_nodeonly.json") else None
    arms = [x for x in ORDER if x in comp["arms"]]
    N, E = comp["copies"], sup[E_KEY]
    rows, v, tgt = [], {}, {}
    for x in arms:
        c, s = comp["arms"][x], sup["arms"][x]
        fsj = load(f"{S}/fs_{x}.json") if truths else {"truths": {}}
        fs = fsj["truths"]
        asm_s, asm_gb = wall(f"{S}/{x}/run.assemble.driver.stderr")
        fam_s, fam_gb = wall(f"{S}/{x}/run.families.driver.stderr")
        try:
            with open(f"{S}/{x}/run.paf_records") as fh:
                paf = int(fh.read())
        except OSError:
            paf = None
        m2n = sup_nodes["arms"][x][M2_KEY] if sup_nodes else None
        r = dict(arm=x, transcripts=count_rows(f"{S}/{x}/run.gtf", "transcript"), loci=c["n_loci"], paf=paf, asm_s=asm_s, asm_gb=asm_gb, fam_s=fam_s, fam_gb=fam_gb,
                 fam_clusters=len(c["family_clusters"]), nodes=c["nodes"], np=c["np"], on_copy=c["classes"]["on_copy"], in_span=c["classes"]["in_span"],
                 antisense=c["classes"]["antisense"], elsewhere=c["classes"]["elsewhere"], shared_spans=c["shared_spans"],
                 red=(c["classes"]["on_copy"] / c["copies_with_node"]) if c["copies_with_node"] else 0.0,
                 cpstar=(c["copies_with_node"] / c["nodes"]) if c["nodes"] else 0.0,
                 M1=c["copies_with_node"], M1_cs=s["old_overlap_in_npip_nodes"], M2=s[M2_KEY], M2n=m2n, M3=s[M3_KEY], M6=c.get("copies_in_kstar"),
                 ann=s.get("ann_found_in_npip_nodes"), chain=s.get("chain_found_in_npip_nodes"))
        for t in truths:
            r[f"{t}_sens"], r[f"{t}_prec"], r[f"{t}_F"], r[f"{t}_pairF"] = (fs[t]["sens"], fs[t]["prec"], fs[t]["f"], fs[t].get("pair_f"))
        tgt[x] = target_rows(fsj, targets) if truths else []
        rows.append(r)
        v[x] = r
    lines = []
    out = lines.append
    out(gate_banner(S))
    crit = "Amendment E (capped start + first three introns)" if K == "tc" else "Amendment A (the representative carries >= k annotated introns)"
    out(f"family {a.family}: N = {N} truth copies, E = {E} expressed under {crit}, arm-independent")
    hdr = (f"{'arm':9}{'tx':>7}{'loci':>6}{'PAF':>9}{'asm s':>6}{'fam s':>6}{'GB':>5} | {'fcl':>3}{'nodes':>6}{'NP':>7} ({'on/span/anti/else':>17}) {'red':>5}{'CP*':>5} | "
           f"{'M1':>3}{'M2':>4}{'M2n':>4}{'M3':>4}{'M6':>4}{'ann':>4}{'chn':>4} |")
    for t in truths:
        hdr += f" {t} sens/prec/F/pairF"
    out(hdr)
    for r in rows:
        tx = "-" if r["transcripts"] is None else r["transcripts"]
        paf = "-" if r["paf"] is None else r["paf"]
        asm = "-" if r["asm_s"] is None else f"{r['asm_s']:.0f}"
        fam = "-" if r["fam_s"] is None else f"{r['fam_s']:.0f}"
        gbs = [g for g in (r["asm_gb"], r["fam_gb"]) if g is not None]
        gb = "-" if not gbs else f"{max(gbs):.1f}"
        ann = "-" if r["ann"] is None else r["ann"]
        chn = "-" if r["chain"] is None else r["chain"]
        m2n = "-" if r["M2n"] is None else r["M2n"]
        line = (f"{r['arm']:9}{tx:>7}{r['loci']:>6}{paf:>9}{asm:>6}{fam:>6}{gb:>5} | "
                f"{r['fam_clusters']:>3}{r['nodes']:>6}{r['np']:>7.3f} ({r['on_copy']:>3}/{r['in_span']:>3}/{r['antisense']:>3}/{r['elsewhere']:>4}) {r['red']:>5.2f}{r['cpstar']:>5.2f} | "
                f"{r['M1']:>3}{r['M2']:>4}{m2n:>4}{r['M3']:>4}{('-' if r['M6'] is None else r['M6']):>4}{ann:>4}{chn:>4} |")
        for t in truths:
            pf = r[f"{t}_pairF"]
            line += f" {r[f'{t}_sens']:.3f}/{r[f'{t}_prec']:.3f}/{r[f'{t}_F']:.3f}/{(f'{pf:.3f}' if pf is not None else '-')}"
        out(line)
    shared = ", ".join("%s %d" % (r["arm"], r["shared_spans"]) for r in rows if r["shared_spans"]) or "none"
    out(f"(tx = assembled transcripts in run.gtf (the +R1 rows copy their base arm); asm s / fam s = wall seconds of the driver's assemble / families call, including any wait for the machine's heavy lock, and for an arm whose all-vs-all took several bounded calls only the LAST call (run.wrapper.log has the shard times); GB = peak RSS, kB / 10^6; "
        f"M1 own node of {N}; M2 the instrument's found count: the copy has an own node and SOME same-strand locus overlapping it, node or not, has the exact representative; "
        f"M2n = the same with the exact locus itself a node; M3 locus-level, of {E}; M6 copies in the largest family cluster K*; red = on-copy nodes per copy with a node "
        f"(several loci on one copy count as correct in NP); CP* = copies with a node / nodes; ann/chn = Amendment A/B found within own nodes; family scores are over ALL families of the "
        f"contig; shared spans per arm: {shared})")
    if targets:
        out("")
        out("target families (rows of the family scores; hit / n_truth and F per arm):")
        for t, ids in targets.items():
            for fid in ids:
                cells = []
                for x in arms:
                    rr = [q for q in tgt[x] if q["truth"] == t and q["family_id"] == fid]
                    q = rr[0] if rr else None
                    cells.append("-" if q is None or q["f"] is None else f"{x} {q['hit']}/{q['n_truth']} F{q['f']:.3f}")
                out(f"  {t} {fid}: " + "; ".join(cells))
    # per-copy matrix
    with open(f"{S}/support.copies.tsv") as fh:
        crow = list(csv.DictReader(fh, delimiter="\t"))
    out("")
    out("per copy: N = own node, F = found within the own node (M2's reading), L = locus-level found; '-' / '.' = no; e = the copy is expressed")
    out(f"{'copy':14}{'e':>2}  " + " ".join(f"{x:>8}" for x in arms))
    mat = []
    for cr in crow:
        cells = []
        for x in arms:
            own = cr.get(f"{x}_page_own_node") == "1"
            f_ = own and cr.get(f"{x}_{K}_found") == "1"
            l_ = cr.get(f"{x}_locus_{K}_found") == "1"
            cells.append(("N" if own else "-") + ("F" if f_ else ".") + ("L" if l_ else "."))
        mat.append([cr["name"], cr[E_KEY]] + cells)
        out(f"{cr['name']:14}{'e' if cr[E_KEY] == '1' else ' ':>2}  " + " ".join(f"{c:>8}" for c in cells))
    # non-dominated sets
    vec = {x: (v[x]["M1"], v[x]["M2"], v[x]["np"]) for x in arms}
    out("")
    for label, group in (("primary arms", [x for x in PRIMARY if x in arms]), ("all arms", arms)):
        nd = [x for x in group if not any(dominated(vec[x], vec[y]) for y in group if y != x)]
        out(f"non-dominated on (M1, M2, NP), {label}: {', '.join(f'{x} {vec[x][0]}/{vec[x][1]}/{vec[x][2]:.3f}' for x in nd)}")
    # predictions
    out("")
    under = E < 8
    f_keys = [t for t in ("u2", "compara") if t in truths] if a.family == "NPIP" else [t for t in ("compara", "soto") if t in truths]

    def verdict(name, text, cond, needs_m2=False):
        if needs_m2 and under:
            out(f"{name}: NO VERDICT (UNDERPOWERED, E = {E} < 8)  {text}")
        elif cond is None:
            out(f"{name}: NO VERDICT (an arm is missing)  {text}")
        else:
            out(f"{name}: {'HELD' if cond else 'FAILED'}  {text}")

    def has(*xs):
        return all(x in v for x in xs)

    g = lambda x, k: v[x][k]  # noqa: E731
    verdict("S1", "M2(P) >= M2(G98) >= M2(A), one strict" + (f"  [{g('P','M2')} >= {g('G98','M2')} >= {g('A','M2')}]" if has("P", "G98", "A") else ""),
            (g("P", "M2") >= g("G98", "M2") >= g("A", "M2") and (g("P", "M2") > g("G98", "M2") or g("G98", "M2") > g("A", "M2"))) if has("P", "G98", "A") else None, True)
    verdict("S2", "M1(G98) >= M1(P)" + (f"  [{g('G98','M1')} >= {g('P','M1')}]" if has("P", "G98") else ""), (g("G98", "M1") >= g("P", "M1")) if has("P", "G98") else None)
    verdict("S3", "NP(P) >= NP(G98) > NP(A)" + (f"  [{g('P','np'):.3f} >= {g('G98','np'):.3f} > {g('A','np'):.3f}]" if has("P", "G98", "A") else ""),
            (g("P", "np") >= g("G98", "np") > g("A", "np")) if has("P", "G98", "A") else None)

    def df(x, y):      # the largest family-score difference over the truths, in whole thousandths
        return max(dthousandths(v[x][f"{t}_F"], v[y][f"{t}_F"]) for t in f_keys)

    if f_keys:
        verdict("S4", f"M2(G98+R1) >= M2(G98)+1, M1 not lower, |dF| < .01 ({'/'.join(f_keys)})" + (
            f"  [M2 {g('G98_R1','M2')} vs {g('G98','M2')}, M1 {g('G98_R1','M1')} vs {g('G98','M1')}, max|dF| {df('G98_R1','G98') / 1000:.3f}]" if has("G98", "G98_R1") else ""),
            (g("G98_R1", "M2") >= g("G98", "M2") + 1 and g("G98_R1", "M1") >= g("G98", "M1") and df("G98_R1", "G98") < 10) if has("G98", "G98_R1") else None, True)
    else:
        verdict("S4g", "M2(G98+R1) >= M2(G98)+1, M1 and M6 not lower (no family truth: the F clause is replaced by M6)" + (
            f"  [M2 {g('G98_R1','M2')} vs {g('G98','M2')}, M1 {g('G98_R1','M1')} vs {g('G98','M1')}, M6 {g('G98_R1','M6')} vs {g('G98','M6')}]" if has("G98", "G98_R1") else ""),
            (g("G98_R1", "M2") >= g("G98", "M2") + 1 and g("G98_R1", "M1") >= g("G98", "M1") and g("G98_R1", "M6") >= g("G98", "M6")) if has("G98", "G98_R1") else None, True)
    verdict("S5", "M2(A+R1) >= M2(A)+2 and M1 not lower" + (f"  [M2 {g('A_R1','M2')} vs {g('A','M2')}, M1 {g('A_R1','M1')} vs {g('A','M1')}]" if has("A", "A_R1") else ""),
            (g("A_R1", "M2") >= g("A", "M2") + 2 and g("A_R1", "M1") >= g("A", "M1")) if has("A", "A_R1") else None, True)
    pr1, pp = f"{S}/P_R1/run.fam.clusters.tsv", f"{S}/P/run.fam.clusters.tsv"
    if os.path.exists(pr1) and os.path.exists(pp):
        with open(pr1) as f1, open(pp) as f2:
            same = f1.read() == f2.read()
    else:
        same = None
    verdict("S6", "P+R1 is identical to P", same)
    if has("G98", "G95", "G995"):
        d1 = max(abs(g("G95", "M1") - g("G98", "M1")), abs(g("G995", "M1") - g("G98", "M1")))
        d2 = max(abs(g("G95", "M2") - g("G98", "M2")), abs(g("G995", "M2") - g("G98", "M2")))
        verdict("S8", f"M1 and M2 of G95 and G995 within one copy of G98  [max |dM1| {d1}, max |dM2| {d2}]", d1 <= 1 and d2 <= 1, True)
    else:
        verdict("S8", "M1 and M2 of G95 and G995 within one copy of G98", None, True)
    text = "\n".join(lines)
    print(text)
    with open(f"{S}/table.txt", "w") as fh:
        fh.write(text + "\n")
    with open(f"{S}/table.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # markdown versions for the results document
    def fmt(val, d=3):
        return "-" if val is None else f"{val:.{d}f}"

    md = ["| arm | assembled transcripts | loci | PAF records | nodes (on-copy / in-span / antisense / elsewhere) | NP | red | CP* | M1 | M2 | M2n | M3 | M6 |"
          + "".join(f" {t} sens / prec / F (pair F) |" for t in truths),
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|" + "---|" * len(truths)]
    for r in rows:
        md.append("| %s | %s | %s | %s | %d (%d / %d / %d / %d) | %.3f | %.2f | %.2f | %d | %d | %s | %d | %s |" % (
            r["arm"].replace("_R1", "+R1"), "-" if r["transcripts"] is None else f"{r['transcripts']:,}", f"{r['loci']:,}", "-" if r["paf"] is None else f"{r['paf']:,}",
            r["nodes"], r["on_copy"], r["in_span"], r["antisense"], r["elsewhere"], r["np"], r["red"], r["cpstar"], r["M1"], r["M2"],
            "-" if r["M2n"] is None else r["M2n"], r["M3"], "-" if r["M6"] is None else r["M6"])
            + "".join(" %.3f / %.3f / %.3f (%s) |" % (r[f"{t}_sens"], r[f"{t}_prec"], r[f"{t}_F"], fmt(r[f"{t}_pairF"])) for t in truths))
    with open(f"{S}/table.md", "w") as fh:
        fh.write("\n".join(md) + "\n")
    with open(f"{S}/cost.md", "w") as fh:
        fh.write("Wall seconds of the driver call (including any wait for the heavy lock; the LAST bounded call only for a sharded arm) and peak RSS (kB / 10^6). "
                 "The first good-pool arm of a substrate also pays the one-time parse of the best-AS table (about 6 s on chr16).\n\n")
        fh.write("| arm | assemble s | assemble peak GB | families s | families peak GB |\n|---|---|---|---|---|\n")
        for r in rows:
            fh.write("| %s | %s | %s | %s | %s |\n" % (r["arm"].replace("_R1", "+R1"), fmt(r["asm_s"], 1), fmt(r["asm_gb"], 2), fmt(r["fam_s"], 1), fmt(r["fam_gb"], 2)))
    if targets:
        with open(f"{S}/target.md", "w") as fh:
            fh.write("| truth family | " + " | ".join(x.replace("_R1", "+R1") for x in arms) + " |\n|---|" + "---|" * len(arms) + "\n")
            for t, ids in targets.items():
                for fid in ids:
                    cells = []
                    for x in arms:
                        rr = [q for q in tgt[x] if q["truth"] == t and q["family_id"] == fid]
                        q = rr[0] if rr else None
                        cells.append("-" if q is None or q["f"] is None else f"{q['hit']}/{q['n_truth']}, F {q['f']:.3f}")
                    fh.write(f"| {t} {fid} | " + " | ".join(cells) + " |\n")
    with open(f"{S}/matrix.tsv", "w") as fh:
        fh.write("copy\tE\t" + "\t".join(arms) + "\n")
        for m in mat:
            fh.write("\t".join(m) + "\n")


if __name__ == "__main__":
    main()
