#!/usr/bin/env python3
"""W/artifact/data.json and index.html for the artifact. Every number comes from the two runs' fate/side/score/compare files and the excision
table; nothing is typed in. Needs both W/mat and W/pat to be complete (Tasks 4-9 for each reference, then side.py and compare.py)."""
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

EXC = "/home/juanfra/winloci_scratch/o3_excise/per_family2.json"
LABEL = {"LRPAP1_p12": "LRPAP1 5′ fragment p12 (descriptive)", "LRPAP1_c07": "LRPAP1 chrY copy (sex control)"}
TITLE = {"mat": "The mother's genome is the reference: copies only the father has are missing",
         "pat": "The father's genome is the reference: copies only the mother has are missing"}


def compact_reads(rows):
    """[read, fate, ref, start, de] -> [0, fate, None, None, de]: the page needs fate and divergence only"""
    return [[0, r[1], None, None, r[4]] for r in rows]


def compact_pairs(rows):
    """[read, fate, de ref, MAPQ ref, de other, MAPQ other] -> the same without the read name"""
    return [[0] + list(r[1:]) for r in rows]


def parse_unm(text):
    """W/unm.txt -> {'total': N, 'mat': K, 'pat': K, 'mat_region': {acc, start, end, n}} (the largest region of the reads mapping on mat) or None"""
    out, hap = {}, None
    for ln in text.splitlines():
        m = re.match(r"R_unm on (mat|pat): (\d+) reads, (\d+) with", ln)
        if m:
            hap = m.group(1)
            out["total"], out[hap] = int(m.group(2)), int(m.group(3))
            continue
        r = re.match(r"\s+region (\S+):(\d+)-(\d+): (\d+) reads", ln)
        if r and hap == "mat" and "mat_region" not in out:
            out["mat_region"] = dict(acc=r.group(1), start=int(r.group(2)), end=int(r.group(3)), n=int(r.group(4)))
    return out or None


def light(row):
    """the shared-control row without its (29k) read list"""
    return {k: v for k, v in row.items() if k != "reads"}


def load(path, default=None):
    return json.load(open(path)) if os.path.exists(path) else default


def direction(a):
    """one direction: reference haplotype `a`; `b` = the haplotype that has the copies"""
    b = "pat" if a == "mat" else "mat"
    WA = f"{C.W}/{a}"
    fate = load(f"{WA}/fate/fate.json")
    side = load(f"{WA}/fate/side.json", {})
    kinds = {r["locus"]: r for r in csv.DictReader(open(f"{WA}/truth/loci.tsv"), delimiter="\t")}
    loci = []
    for k, v in fate["loci"].items():
        if v["n"] == 0:
            continue
        s = side.get(k, {})
        loci.append(dict(id=k, label=LABEL.get(k, k), kind=v["kind"], n=v["n"], fates=v["fates"], fractions=v["fractions"],
                         verdict=v["verdict"], de_median=v["de_median"], paralog_identity=v["paralog_identity"], reads=compact_reads(v["reads"]),
                         pairs=compact_pairs(s.get("pairs", [])), de_other_median=s.get("de_other_median"), mapq_ref_median=s.get("mapq_ref_median"),
                         mapq_other_median=s.get("mapq_other_median"), ref_track=s.get("ref_track"), other_track=s.get("other_track")))
    loci.sort(key=lambda x: (x["kind"] == "sex", -x["n"]))
    arms = {arm: load(f"{WA}/score/{arm}.json") for arm in ("isocon", "inhouse")}
    cmp_ = load(f"{C.W}/compare_{a}.json", {}).get("loci", {})
    return dict(ref=a, other=b, title=TITLE[a], loci=loci, shared=light(fate["shared"]),
                arms={k: {x: v[x] for x in ("R1", "R2", "R3", "R4", "candidates")} for k, v in arms.items() if v}, p12=(arms.get("isocon") or {}).get("p12"),
                compare=cmp_)


def main():
    unm = parse_unm(open(f"{C.W}/unm.txt").read()) if os.path.exists(f"{C.W}/unm.txt") else None
    if unm and "mat_region" in unm:
        al = C.alias()
        unm["mat_region"]["chrom"] = next((num for (h, num), acc in al.items() if h == "mat" and acc == unm["mat_region"]["acc"]), None)
    out = dict(directions={a: direction(a) for a in ("mat", "pat") if os.path.exists(f"{C.W}/{a}/fate/fate.json")},
               context=[dict(fam=r["fam"], unaln=r["unaln"], conc=r["conc"], mig_de=r["mig_de"]) for r in json.load(open(EXC))],
               unm=unm)
    os.makedirs(f"{C.W}/artifact", exist_ok=True)
    json.dump(out, open(f"{C.W}/artifact/data.json", "w"))
    html = open(f"{HERE}/template.html").read().replace("/*DATA*/null", json.dumps(out))
    open(f"{C.W}/artifact/index.html", "w").write(html)
    print(f"data.json: directions {list(out['directions'])}, loci {[len(d['loci']) for d in out['directions'].values()]}, "
          f"{len(out['context'])} excision families; index.html {len(html) // 1024} KB")


if __name__ == "__main__":
    main()
