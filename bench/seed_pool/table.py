#!/usr/bin/env python3
"""The comparison tables of docs/PREREG_seed_pool_real_reads_2026-10-07.md, computed from the products of run.sh score (nothing by hand).

    table.py --work W/SUB --family NPIP --truths compara,u2,soto

Prints (and writes W/SUB/table.txt, table.tsv, matrix.tsv): (1) one row per arm: cost, the family's clusters, M1-M3, the family scores;
(2) the per-copy matrix (own node N, E-found within the own node F, locus-level E-found L); (3) the non-dominated arms on (M1, M2, NP);
(4) the registered predictions S1-S6, S8, each HELD / FAILED / NO VERDICT (UNDERPOWERED when fewer than 8 copies are E-expressed).
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
    """'Elapsed (wall clock) time (h:mm:ss or m:ss): 1:05.12' -> seconds, and the peak RSS in GB."""
    try:
        txt = open(stderr_path).read()
    except OSError:
        return None, None
    m = re.search(r"Elapsed \(wall clock\) time[^:]*: ([\d:.]+)", txt)
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--truths", default="", help="comma list of family truths with a fs_<arm>.json; empty = none on this substrate (M4 not scored)")
    ap.add_argument("--found", default="tc", choices=["tc", "ann"],
                    help="the found criterion: tc = Amendment E (capped start + first three introns; human reads), ann = Amendment A (gorilla: no cap signal)")
    a = ap.parse_args()
    S, truths = a.work, [t for t in a.truths.split(",") if t]
    K = a.found                                           # column prefix in copy_support's output
    M2_KEY, M3_KEY, E_KEY = f"{K}_found_in_npip_nodes", f"locus_{K}_found", f"{K}_expressed"
    comp, sup = load(f"{S}/comp.json"), load(f"{S}/support.json")
    arms = [x for x in ORDER if x in comp["arms"]]
    N, E = comp["copies"], sup[E_KEY]
    rows, v = [], {}
    for x in arms:
        c, s = comp["arms"][x], sup["arms"][x]
        fs = load(f"{S}/fs_{x}.json")["truths"] if truths else {}
        asm_s, asm_gb = wall(f"{S}/{x}/run.assemble.driver.stderr")
        fam_s, fam_gb = wall(f"{S}/{x}/run.families.driver.stderr")
        try:
            paf = int(open(f"{S}/{x}/run.paf_records").read())
        except OSError:
            paf = None
        r = dict(arm=x, transcripts=count_rows(f"{S}/{x}/run.gtf", "transcript"), loci=c["n_loci"], paf=paf, asm_s=asm_s, fam_s=fam_s,
                 fam_clusters=len(c["family_clusters"]), nodes=c["nodes"], np=c["np"], on_copy=c["classes"]["on_copy"], in_span=c["classes"]["in_span"],
                 antisense=c["classes"]["antisense"], elsewhere=c["classes"]["elsewhere"], shared_spans=c["shared_spans"],
                 M1=c["copies_with_node"], M1_cs=s["old_overlap_in_npip_nodes"], M2=s[M2_KEY], M3=s[M3_KEY], M6=c.get("copies_in_kstar"),
                 ann=s.get("ann_found_in_npip_nodes"), chain=s.get("chain_found_in_npip_nodes"))
        for t in truths:
            r[f"{t}_sens"], r[f"{t}_prec"], r[f"{t}_F"], r[f"{t}_pairF"] = (fs[t]["sens"], fs[t]["prec"], fs[t]["f"], fs[t].get("pair_f"))
        rows.append(r)
        v[x] = r
    lines = []
    out = lines.append
    crit = "Amendment E (capped start + first three introns)" if K == "tc" else "Amendment A (the representative carries >= k annotated introns)"
    out(f"family {a.family}: N = {N} truth copies, E = {E} expressed under {crit}, arm-independent")
    hdr = f"{'arm':9}{'tx':>7}{'loci':>6}{'PAF':>9}{'asm s':>6}{'fam s':>7} | {'fcl':>3}{'nodes':>6}{'NP':>6} ({'on/span/anti/else':>17}) | {'M1':>3}{'M2':>4}{'M3':>4}{'M6':>4}{'ann':>4}{'chn':>4} |"
    for t in truths:
        hdr += f" {t} sens/prec/F/pairF"
    out(hdr)
    for r in rows:
        tx = "-" if r["transcripts"] is None else r["transcripts"]
        paf = "-" if r["paf"] is None else r["paf"]
        asm = "-" if r["asm_s"] is None else f"{r['asm_s']:.0f}"
        fam = "-" if r["fam_s"] is None else f"{r['fam_s']:.0f}"
        ann = "-" if r["ann"] is None else r["ann"]
        chn = "-" if r["chain"] is None else r["chain"]
        line = (f"{r['arm']:9}{tx:>7}{r['loci']:>6}{paf:>9}{asm:>6}{fam:>7} | "
                f"{r['fam_clusters']:>3}{r['nodes']:>6}{r['np']:>6.2f} ({r['on_copy']:>3}/{r['in_span']:>3}/{r['antisense']:>3}/{r['elsewhere']:>4}) | "
                f"{r['M1']:>3}{r['M2']:>4}{r['M3']:>4}{('-' if r['M6'] is None else r['M6']):>4}{ann:>4}{chn:>4} |")
        for t in truths:
            pf = r[f"{t}_pairF"]
            line += f" {r[f'{t}_sens']:.3f}/{r[f'{t}_prec']:.3f}/{r[f'{t}_F']:.3f}/{(f'{pf:.3f}' if pf is not None else '-')}"
        out(line)
    shared = ", ".join("%s %d" % (r["arm"], r["shared_spans"]) for r in rows if r["shared_spans"]) or "none"
    out(f"(M1 own node of {N}; M2 E-found within own nodes and M3 locus-level E-found, of {E}; M6 copies in the largest family cluster K*; ann/chn = Amendment A/B found within own nodes; "
        f"shared spans per arm: {shared})")
    # (2) per-copy matrix
    with open(f"{S}/support.copies.tsv") as fh:
        crow = list(csv.DictReader(fh, delimiter="\t"))
    out("")
    out("per copy: N = own node, F = E-found within the own node, L = locus-level E-found; '-' / '.' = no; e = the copy is E-expressed")
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
    # (3) non-dominated sets
    vec = {x: (v[x]["M1"], v[x]["M2"], v[x]["np"]) for x in arms}
    out("")
    for label, group in (("primary arms", [x for x in PRIMARY if x in arms]), ("all arms", arms)):
        nd = [x for x in group if not any(dominated(vec[x], vec[y]) for y in group if y != x)]
        out(f"non-dominated on (M1, M2, NP), {label}: {', '.join(f'{x} {vec[x][0]}/{vec[x][1]}/{vec[x][2]:.2f}' for x in nd)}")
    # (4) predictions
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
    verdict("S3", "NP(P) >= NP(G98) > NP(A)" + (f"  [{g('P','np'):.2f} >= {g('G98','np'):.2f} > {g('A','np'):.2f}]" if has("P", "G98", "A") else ""),
            (g("P", "np") >= g("G98", "np") > g("A", "np")) if has("P", "G98", "A") else None)

    def df(x, y):
        return max(abs(v[x][f"{t}_F"] - v[y][f"{t}_F"]) for t in f_keys)

    if f_keys:
        verdict("S4", f"M2(G98+R1) >= M2(G98)+1, M1 not lower, |dF| < .01 ({'/'.join(f_keys)})" + (
            f"  [M2 {g('G98_R1','M2')} vs {g('G98','M2')}, M1 {g('G98_R1','M1')} vs {g('G98','M1')}, max|dF| {df('G98_R1','G98'):.3f}]" if has("G98", "G98_R1") else ""),
            (g("G98_R1", "M2") >= g("G98", "M2") + 1 and g("G98_R1", "M1") >= g("G98", "M1") and df("G98_R1", "G98") < 0.01) if has("G98", "G98_R1") else None, True)
    else:
        verdict("S4g", "M2(G98+R1) >= M2(G98)+1, M1 and M6 not lower (no family truth: the F clause is replaced by M6)" + (
            f"  [M2 {g('G98_R1','M2')} vs {g('G98','M2')}, M1 {g('G98_R1','M1')} vs {g('G98','M1')}, M6 {g('G98_R1','M6')} vs {g('G98','M6')}]" if has("G98", "G98_R1") else ""),
            (g("G98_R1", "M2") >= g("G98", "M2") + 1 and g("G98_R1", "M1") >= g("G98", "M1") and g("G98_R1", "M6") >= g("G98", "M6")) if has("G98", "G98_R1") else None, True)
    verdict("S5", "M2(A+R1) >= M2(A)+2 and M1 not lower" + (f"  [M2 {g('A_R1','M2')} vs {g('A','M2')}, M1 {g('A_R1','M1')} vs {g('A','M1')}]" if has("A", "A_R1") else ""),
            (g("A_R1", "M2") >= g("A", "M2") + 2 and g("A_R1", "M1") >= g("A", "M1")) if has("A", "A_R1") else None, True)
    p_same = os.path.exists(f"{S}/P_R1/run.fam.clusters.tsv") and open(f"{S}/P_R1/run.fam.clusters.tsv").read() == open(f"{S}/P/run.fam.clusters.tsv").read()
    verdict("S6", "P+R1 is identical to P", p_same if os.path.exists(f"{S}/P_R1/run.fam.clusters.tsv") else None)
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
    with open(f"{S}/matrix.tsv", "w") as fh:
        fh.write("copy\tE\t" + "\t".join(arms) + "\n")
        for m in mat:
            fh.write("\t".join(m) + "\n")


if __name__ == "__main__":
    main()
