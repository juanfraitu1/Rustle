#!/usr/bin/env python3
"""Amendment 9 on a bed: IsoCon-style correction of the cluster consensus sequences, then the registered fidelity metric before and after.

    run_polish.py <bed> polish    # W/<bed>/registered/polish/{cons.polished.fa, polish.json}; resumable (exit 75 = run again, one time-boxed batch per call)
    run_polish.py <bed> score     # aligns original and polished consensus to the unmasked genome (--cs) -> polish/score.json
bed = bedA (dev) | bedH (held-out). Heavy: the score step loads the 15 GB genome index; run both steps under tools/rlock.sh heavy."""
import csv
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fidelity as F  # noqa: E402
import polish as P  # noqa: E402
import score as S  # noqa: E402
import seeds as SD  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
GENOME_IDX = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO.splice.mmi"
PANELS = {"bedA": "/home/juanfra/winloci_scratch/o3_excise/panel.json", "bedH": "/mnt/linuxdisk/tmp/rna_allele/linktest/panel.json"}
BATCH_SECONDS = 420


import re as _re


def R_CIG(cigar):
    return _re.findall(r"(\d+)([MIDNSHP=X])", cigar)


def read_cons(path):
    out, name = {}, None
    for ln in open(path):
        if ln[0] == ">":
            name = ln[1:].strip().split()[0]
            out[name] = []
        else:
            out[name].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


def registered(d):
    """the folder holding clusters.tsv and cons.fa: <bed>/registered, or the folder itself (the LRPAP1 example)"""
    return f"{d}/registered" if os.path.isdir(f"{d}/registered") else d


def cluster_reads(d):
    cl = {}
    for r in csv.DictReader(open(f"{registered(d)}/clusters.tsv"), delimiter="\t"):
        cl.setdefault("cl" + r["cluster"], []).append(r["read"])
    return cl


def polish_all(bed):
    d = f"{W}/{bed}"
    out = f"{registered(d)}/polish"
    os.makedirs(f"{out}/tmp", exist_ok=True)
    seqs = SD.read_fa(f"{d}/pool.fa")
    cons = read_cons(f"{registered(d)}/cons.fa")
    cl = cluster_reads(d)
    done = json.load(open(f"{out}/polish.json")) if os.path.exists(f"{out}/polish.json") else {}
    t0 = time.time()
    for name, c0 in cons.items():
        if name in done:
            continue
        if time.time() - t0 > BATCH_SECONDS:
            json.dump(done, open(f"{out}/polish.json", "w"), indent=1)
            print(f"{len(done)} of {len(cons)} clusters done; run again")
            sys.exit(75)
        key = name.split("|")[0]
        rf = f"{out}/tmp/reads.fa"
        SD.write_fa(rf, seqs, sorted(cl[key]))

        def align(c):
            cf = f"{out}/tmp/cons.fa"
            open(cf, "w").write(f">c\n{c}\n")
            return subprocess.run(f"minimap2 -ax map-hifi --eqx -t 2 {cf} {rf}", shell=True, capture_output=True, text=True, check=True).stdout.splitlines()

        c1, rounds, applied, minority = P.polish(c0, align)
        done[name] = dict(seq=c1, rounds=rounds, applied=len(applied), kinds={k: sum(1 for a in applied if a[0] == k) for k in ("sub", "del", "ins")},
                          minority=[list(m) for m in minority], len_before=len(c0), len_after=len(c1), reads=len(cl[key]))
    json.dump(done, open(f"{out}/polish.json", "w"), indent=1)
    with open(f"{out}/cons.polished.fa", "w") as o:
        for name, v in done.items():
            o.write(f">{name}\n{v['seq']}\n")
    print(f"polished {len(done)} clusters; corrections {sum(v['applied'] for v in done.values())}; clusters changed {sum(1 for v in done.values() if v['applied'])}")


def best_records(paf_lines, prefix):
    """{cluster key: record} best PAF record per query by matches (same rule as fidelity.best_hits), with the NM tag and cs string"""
    b = {}
    for ln in paf_lines:
        f = ln.rstrip("\n").split("\t")
        if not f[0].startswith(prefix):
            continue
        q = f[0][len(prefix):].split("|")[0]
        m = int(f[9])
        if q not in b or m > b[q]["m"]:
            tags = {t[:2]: t[5:] for t in f[12:]}
            b[q] = dict(m=m, ident=m / max(1, int(f[10])), cov=(int(f[3]) - int(f[2])) / int(f[1]), ref=f[5], start=int(f[7]), end=int(f[8]),
                        nm=int(tags.get("NM", 0)), cs=tags.get("cs", ""), qstart=int(f[2]), qend=int(f[3]), qlen=int(f[1]), strand=f[4])
    return b


def score_bed(bed):
    d = f"{W}/{bed}"
    out = f"{d}/registered/polish"
    pol = json.load(open(f"{out}/polish.json"))
    orig = read_cons(f"{d}/registered/cons.fa")
    with open(f"{out}/both.fa", "w") as o:
        for name, s in orig.items():
            o.write(f">orig.{name}\n{s}\n>pol.{name}\n{pol[name]['seq']}\n")
    paf = f"{out}/both.paf"
    if not os.path.exists(paf + ".done"):
        subprocess.run(f"minimap2 -c --cs -x splice:hq -uf -N 5 -t 4 {GENOME_IDX} {out}/both.fa > {paf}", shell=True, check=True)
        open(paf + ".done", "w").write("ok")
    lines = open(paf).read().splitlines()
    rec = {"orig": best_records(lines, "orig."), "pol": best_records(lines, "pol.")}
    lab = {r["read"]: r for r in csv.DictReader(open(f"{d}/labels.tsv"), delimiter="\t")}
    truth = {n: ("bg" if x["role"] == "bg" else (x["family"] or None)) for n, x in lab.items()}
    erased = {p["fam"]: tuple(p["mask"][:3]) for p in json.load(open(PANELS[bed]))}
    cl = cluster_reads(d)
    rows, edit_cols = {}, {}
    for name in orig:
        key = name.split("|")[0]
        fam = S.majority(cl[key], truth)
        if fam is None or fam not in erased:
            continue
        r = {}
        for which in ("orig", "pol"):
            h = rec[which].get(key)
            e = erased[fam]
            ok = h is not None and h["ref"] == e[0] and h["start"] < e[2] and e[1] < h["end"]
            r[which] = dict(on=ok, idcov=(h["ident"] * h["cov"]) if ok else None, nm=h["nm"] if ok else None)
        rows[key] = dict(fam=fam, **r)
        h = rec["pol"].get(key)
        if h and r["pol"]["on"]:
            edit_cols[key] = P.cs_edit_columns(h["cs"], h["qstart"], h["qend"], h["qlen"], h["strand"])
    bar = 0.999
    on = [k for k, v in rows.items() if v["orig"]["on"] and v["pol"]["on"]]
    res = dict(bed=bed, clusters_with_family=len(rows), on_erased_copy=dict(orig=sum(v["orig"]["on"] for v in rows.values()), pol=sum(v["pol"]["on"] for v in rows.values())),
               idcov_ge_0_999=dict(orig=sum(1 for v in rows.values() if v["orig"]["on"] and v["orig"]["idcov"] >= bar), pol=sum(1 for v in rows.values() if v["pol"]["on"] and v["pol"]["idcov"] >= bar)),
               lost=sum(1 for k in on if rows[k]["pol"]["idcov"] < rows[k]["orig"]["idcov"] - 1e-12), gained=sum(1 for k in on if rows[k]["pol"]["idcov"] > rows[k]["orig"]["idcov"] + 1e-12),
               crossed_up=sum(1 for k in on if rows[k]["orig"]["idcov"] < bar <= rows[k]["pol"]["idcov"]), crossed_down=sum(1 for k in on if rows[k]["pol"]["idcov"] < bar <= rows[k]["orig"]["idcov"]),
               nm=dict(orig=sum(rows[k]["orig"]["nm"] for k in on), pol=sum(rows[k]["pol"]["nm"] for k in on)))
    # regression check: the original numbers must reproduce the stored fidelity result
    stored = json.load(open(f"{d}/registered/fidelity.json"))["all_clusters"]
    res["stored_fidelity"] = stored
    below = [k for k in on if rows[k]["pol"]["idcov"] < bar]
    tot = hit = 0
    for k in below:
        cols = edit_cols.get(k, set())
        mino = [m for m in next(v for n, v in pol.items() if n.split("|")[0] == k)["minority"]]
        pos = [m[1] for m in mino]
        tot += len(cols)
        hit += sum(any(abs(c - p) <= 1 for p in pos) for c in cols)
    res["d0"] = dict(clusters_below_after=len(below), edit_columns=tot, allele_like_columns=hit, share=(hit / tot if tot else None))
    res["descriptive"] = dict(clusters=len(pol), changed=sum(1 for v in pol.values() if v["applied"]), corrections=sum(v["applied"] for v in pol.values()),
                              kinds={k: sum(v["kinds"][k] for v in pol.values()) for k in ("sub", "del", "ins")},
                              length_cut_by_more_than_10_percent=sum(1 for v in pol.values() if v["len_after"] < 0.9 * v["len_before"]),
                              clusters_with_a_minority_variant=sum(1 for v in pol.values() if v["minority"]))
    json.dump(dict(summary=res, clusters=rows), open(f"{out}/score.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    {"polish": polish_all, "score": score_bed}[sys.argv[2]](sys.argv[1])
