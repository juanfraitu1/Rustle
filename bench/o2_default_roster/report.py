#!/usr/bin/env python3
"""O2 read-truth report on one simulated roster (docs/PREREG_o2_default_roster_2026-10-06.md): the registered scorer's per-read verdicts
(`bench/score.py reads --per-read`, OWN / PRIMARY / ANY, and the union-certificate arm), joined with the simulated BAM by figures/_o2.py, tabulated
for the MAPQ-0 reads whose source copy lies in the target set and for every MAPQ-0 read of the roster, by identity band, with the registered bars.

    report.py --work DIR --contig chr16 --target NPIP|TBC1D3|Y|none --sim PREFIX --catalog CAT.tsv --o2 PREFIX --u2 PREFIX [--logs DIR]

Target set T = the roster copies (same chromosome and strand) whose exons overlap >= 1 bp
  NPIP / TBC1D3  the exon union of a CAT/Liftoff copy of that family (copy_recovery_tools_cat/ann: copies.hsa.tsv + truth.hsa.gtf)
  Y              the body of a RefSeq chrY gene named DAZ*, RBMY*, TSPY*, BPY2*, HSFY*, VCY*, CDY*, PRY* or XKRY* (families_gw/.../genes_only.gff)
Writes DIR/report.json and prints the tables. Everything is computed here from the scorer's own outputs; the cross-check of figures/_o2.py
(its tallies must equal the ALL rows of `score.py reads`) runs first and stops the report when it fails.
"""
import argparse
import collections
import csv
import json
import os
import re
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "figures"))
sys.path.insert(0, os.path.join(REPO, "bench"))
import _o2  # noqa: E402
import copy_support as cs  # noqa: E402

ANN = "/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools_cat/ann"
GENES = "/mnt/linuxdisk/tmp/rustle_figures/families_gw/species/human/genes_only.gff"
YAG = re.compile(r"^(DAZ|RBMY|TSPY|BPY2|HSFY|VCY|CDY|PRY|XKRY)")
BANDS = ["identical", "99.5-100%", "99-99.5%", "98-99%", "<98%", "unknown"]
LOW_BANDS = ("98-99%", "<98%")          # divergence >= 1%: the bins the registered per-bin bar applies to
MIN_ASSIGNED, MIN_BAND, B2_ACC, B2_BAND_ACC = 30, 10, 0.95, 0.90
VERDICT = {"lost": "not_scored", "no_own_row": "not_scored", "no_primary_row": "not_scored"}


def overlap(a, b):
    i = j = 0
    while i < len(a) and j < len(b):
        if min(a[i][1], b[j][1]) > max(a[i][0], b[j][0]):
            return True
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return False


def target_regions(target, contig):
    """[(strand, [[s0, e], ...])] on `contig`: the target's truth intervals."""
    out = []
    if target in ("NPIP", "TBC1D3"):
        copies = [r for r in csv.DictReader(open(f"{ANN}/copies.hsa.tsv"), delimiter="\t") if r["family"] == target and r["chrom"] == contig]
        genes = {r["isoform_gene"] for r in copies} | {r["cid"] for r in copies}
        tx = cs.gtf_transcripts(f"{ANN}/truth.hsa.gtf", genes)
        for c in copies:
            t = tx.get(c["cid"]) or tx.get(c["isoform_gene"], {})
            u = cs.merge([e for ex in t.values() for e in ex]) or [[int(c["terr_lo0"]), int(c["terr_hi"])]]
            out.append((c["strand"], u))
    elif target == "Y":
        for ln in open(GENES):
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != contig or f[2] not in ("gene", "pseudogene"):
                continue
            nm = re.search(r"(?:^|;)Name=([^;]+)", f[8])
            if nm and YAG.match(nm.group(1)):
                out.append((f[6], [[int(f[3]) - 1, int(f[4])]]))
    return out


def roster(catalog):
    rows = {}
    for r in csv.DictReader(open(catalog), delimiter="\t"):
        ex = [[int(a), int(b)] for a, b in (x.split("-") for x in r["exons"].split(",") if x)]
        rows[(r["family_id"], r["copy_idx"])] = dict(chrom=r["chrom"], strand=r["strand"], exons=cs.merge(sorted(ex)))
    return rows


def tally(recs, key):
    c = collections.Counter(r[key] for r in recs)
    n = len(recs)
    asg = c["correct"] + c["wrong"] + c["conflict"]
    return dict(n=n, correct=c["correct"], wrong=c["wrong"], conflict=c["conflict"], abstain=c["abstain"], not_scored=c["not_scored"],
                assigned=asg, accuracy=(c["correct"] / asg if asg else None), coverage=(asg / n if n else None))


def table(recs, keys):
    out = {}
    for k in keys:
        out[k] = dict(all=tally(recs, k), bands={b: tally([r for r in recs if r["band"] == b], k) for b in BANDS if any(r["band"] == b for r in recs)})
    return out


def bar_b1(t):
    a = t["assigned"]
    if a < MIN_ASSIGNED:
        return dict(verdict="UNDERPOWERED", assigned=a, wrong=t["wrong"] + t["conflict"], upper95=(3.0 / a if a and not (t["wrong"] + t["conflict"]) else None))
    return dict(verdict="PASS" if t["wrong"] + t["conflict"] == 0 else "FAIL", assigned=a, wrong=t["wrong"] + t["conflict"])


def bar_b2(tb):
    t = tb["all"]
    if t["assigned"] < MIN_ASSIGNED:
        return dict(verdict="UNDERPOWERED", assigned=t["assigned"], accuracy=t["accuracy"])
    bad = []
    if t["accuracy"] < B2_ACC:
        bad.append(f"overall {t['accuracy']:.3f} < {B2_ACC}")
    for b in LOW_BANDS:
        x = tb["bands"].get(b)
        if x and x["assigned"] >= MIN_BAND and x["accuracy"] < B2_BAND_ACC:
            bad.append(f"{b} {x['accuracy']:.3f} < {B2_BAND_ACC} (n={x['assigned']})")
    return dict(verdict="PASS" if not bad else "FAIL", assigned=t["assigned"], accuracy=t["accuracy"], failing=bad)


def show(title, tb):
    print(f"\n{title}")
    print(f"  {'reading':8s} {'MAPQ-0':>7s} {'correct':>8s} {'wrong':>6s} {'confl':>6s} {'abstain':>8s} {'other':>6s} {'assigned':>9s} {'accuracy':>9s} {'coverage':>9s}")
    for k, v in tb.items():
        t = v["all"]
        print(f"  {k:8s} {t['n']:7d} {t['correct']:8d} {t['wrong']:6d} {t['conflict']:6d} {t['abstain']:8d} {t['not_scored']:6d} {t['assigned']:9d} "
              f"{('%.3f' % t['accuracy']) if t['accuracy'] is not None else '-':>9s} {('%.3f' % t['coverage']) if t['coverage'] is not None else '-':>9s}")
    for k, v in tb.items():
        print(f"  by identity band, {k}: " + "; ".join(
            f"{b} n={x['n']} asg={x['assigned']} ok={x['correct']} bad={x['wrong'] + x['conflict']} cov={('%.2f' % x['coverage']) if x['coverage'] is not None else '-'}"
            for b, x in v["bands"].items()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--target", required=True, choices=("NPIP", "TBC1D3", "Y", "none"))
    ap.add_argument("--sim", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--o2", required=True)
    ap.add_argument("--u2", required=True)
    ap.add_argument("--logs", default=None)
    ap.add_argument("--exclude-family", action="append", default=[], help="source family left out of the assignment run (Amendment 2); its reads are dropped from every table")
    a = ap.parse_args()
    runs = dict(species="human", sample=f"human_{a.contig}", sim=a.sim, catalog=a.catalog, o2=a.o2, u2=a.u2, catalog_scope="dev")
    logs = _o2.score_logs(runs, a.logs or os.path.join(a.work, "score"))
    recs = _o2.per_read(runs, logs)
    notes = _o2.crosscheck(recs, logs)
    prim = {}
    for r in csv.DictReader(open(logs["o2_per_read"]), delimiter="\t"):
        prim[r["read_name"]] = VERDICT.get(r["PRIMARY"], r["PRIMARY"])
    for r in recs:
        r["primary_outcome"] = prim.get(r["read"], "not_scored")
    z = [r for r in recs if r["mapq"] == 0 and r["family"] not in a.exclude_family]
    n_excluded = sum(1 for r in recs if r["mapq"] == 0 and r["family"] in a.exclude_family)
    ros = roster(a.catalog)
    regions = target_regions(a.target, a.contig) if a.target != "none" else []
    inT = {k for k, v in ros.items() if any(v["chrom"] == a.contig and v["strand"] == s and overlap(v["exons"], iv) for s, iv in regions)}
    zt = [r for r in z if (r["family"], r["copy"]) in inT]
    keys = ["own_outcome", "primary_outcome", "any_outcome", "union_outcome"]
    names = {"own_outcome": "OWN", "primary_outcome": "PRIMARY", "any_outcome": "ANY", "union_outcome": "UNION"}
    out = dict(contig=a.contig, target=a.target, n_sim_reads=len(recs), n_mapq0=len(z), excluded_families=a.exclude_family, n_mapq0_excluded=n_excluded, roster_copies=len(ros), target_copies=len(inT),
               target_copy_ids=sorted("|".join(k) for k in inT), crosscheck=notes)
    out["aligner_primary_true_copy_among_mapq0"] = dict(
        all=(sum(r["aligner"] == "true_copy" for r in z) / len(z) if z else None),
        target=(sum(r["aligner"] == "true_copy" for r in zt) / len(zt) if zt else None))
    out["all"] = {names[k]: v for k, v in table(z, keys).items()}
    out["target_table"] = {names[k]: v for k, v in table(zt, keys).items()} if a.target != "none" else {}
    print(f"contig {a.contig}: {len(recs)} simulated reads, {len(z)} at MAPQ 0 (+{n_excluded} from the excluded families {a.exclude_family}), roster {len(ros)} copies; target {a.target}: {len(inT)} copies, {len(zt)} MAPQ-0 reads")
    for n in notes:
        print("  " + n)
    show(f"ALL MAPQ-0 reads of the roster ({a.contig})", {names[k]: v for k, v in table(z, keys).items()})
    if a.target != "none":
        tt = table(zt, keys)
        show(f"TARGET {a.target}: MAPQ-0 reads whose source copy is one of its {len(inT)} roster copies", {names[k]: v for k, v in tt.items()})
        b1 = bar_b1(tt["own_outcome"]["all"])
        b2 = {names[k]: bar_b2(tt[k]) for k in ("primary_outcome", "any_outcome")}
        b3 = bar_b1(tt["union_outcome"]["all"])
        out["bars"] = dict(B1_OWN=b1, B2=b2, B3_UNION_beside=b3)
        print(f"\nBARS {a.target}: B1 OWN (0 wrong of >= {MIN_ASSIGNED} assigned) {b1}")
        for k, v in b2.items():
            print(f"      B2 {k} (accuracy >= {B2_ACC}, bands {LOW_BANDS} >= {B2_BAND_ACC} with >= {MIN_BAND}; >= {MIN_ASSIGNED} assigned) {v}")
        print(f"      B3 UNION, beside, no bar: {b3}")
        print(f"      aligner-primary on the target's MAPQ-0 reads: {out['aligner_primary_true_copy_among_mapq0']['target']}")
    json.dump(out, open(os.path.join(a.work, "report.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
