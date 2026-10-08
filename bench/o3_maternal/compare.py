#!/usr/bin/env python3
"""Chain side-by-side of the two runs (docs/PREREG_o3_maternal_reference_2026-10-08.md, Amendment 1).

For the loci absent from the reference of run A (A = mat: copies only the father has; A = pat: the reverse) show the SAME copies as seen by
run A (the reference lacks them: are their transcripts flagged as new, are they recovered?) and by run B (the reference has them: are their
transcripts simply found in it?). Locus coordinates are on the truth haplotype of run A, which is the reference of run B.

    compare.py mat | pat      # writes W/compare_<A>.json; needs W/<A>/{truth,isoc,inhouse,score} and W/<B>/{isoc,inhouse}
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

MIN_SCORE = 0.999


def outputs_on(best, loci_, min_score=MIN_SCORE):
    """{locus: [output names]}: outputs whose best hit (identity x coverage >= min_score) overlaps the locus; best = common.best_hits() on the
    haplotype the locus lies on"""
    out = {L["locus"]: [] for L in loci_}
    for q, h in sorted(best.items()):
        if h[0] < min_score:
            continue
        for L in loci_:
            if h[1] == L["chrom"] and h[2] < L["end"] and L["start"] < h[3]:
                out[L["locus"]].append(q)
    return out


def new_outputs(contig_rows):
    """names of the IsoCon outputs the link step kept as new copies (contigs.tsv rows with linked == '0')"""
    return {r["output"].split()[0] for r in contig_rows if r["linked"] == "0"}


def recovering(score_rows, locus):
    """[(candidate, n_transcripts)] of the flagged candidates of a score JSON that recover the locus"""
    return [dict(candidate=r["candidate"], n_transcripts=r["n_transcripts"]) for r in score_rows if locus in r["recovers"]]


def inhouse_near(rows, loci_, al):
    """{locus: [dict(candidate, n_clusters, d, flagged)]}: stage candidates (cands_all.tsv rows) whose nearest locus overlaps the locus"""
    out = {L["locus"]: [] for L in loci_}
    for r in rows:
        chrom, rng = r["nearest_locus"].rsplit(":", 1)
        s, e = (int(x) for x in rng.split("-"))
        acc = C.accession(chrom, al)
        for L in loci_:
            if acc == L["chrom"] and s < L["end"] and L["start"] < e:
                out[L["locus"]].append(dict(candidate=r["candidate"], n_clusters=int(r["n_clusters"]), d=float(r["d"]), flagged=int(r["flagged"])))
    return out


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def main(a):
    b = "pat" if a == "mat" else "mat"
    al = C.alias()
    WA, WB = f"{C.W}/{a}", f"{C.W}/{b}"
    loci_ = tsv(f"{WA}/truth/loci.tsv")
    for L in loci_:
        L["start"], L["end"] = int(L["start"]), int(L["end"])
    best = lambda p: C.best_hits(p, al) if os.path.exists(p) and os.path.getsize(p) else {}
    on_a = outputs_on(best(f"{WA}/isoc/outputs.other.paf"), loci_)      # run A: outputs against the truth haplotype
    on_b = outputs_on(best(f"{WB}/isoc/outputs.pri.paf"), loci_)        # run B: outputs against ITS reference = the same haplotype
    new_a, new_b = new_outputs(tsv(f"{WA}/isoc/contigs.tsv")), new_outputs(tsv(f"{WB}/isoc/contigs.tsv"))
    sc = {arm: json.load(open(f"{WA}/score/{arm}.json")) for arm in ("isocon", "inhouse") if os.path.exists(f"{WA}/score/{arm}.json")}
    near = inhouse_near(tsv(f"{WB}/inhouse/cands_all.tsv"), loci_, al) if os.path.exists(f"{WB}/inhouse/cands_all.tsv") else {}
    out = {}
    for L in loci_:
        k = L["locus"]
        out[k] = dict(
            kind=L["kind"],
            isocon=dict(run_a=dict(outputs=len(on_a[k]), flagged_new=sum(1 for o in on_a[k] if o in new_a),
                                   recovering=recovering(sc.get("isocon", {}).get("rows", []), k)),
                        run_b=dict(outputs=len(on_b[k]), flagged_new=sum(1 for o in on_b[k] if o in new_b))),
            inhouse=dict(run_a=dict(recovering=recovering(sc.get("inhouse", {}).get("rows", []), k)),
                         run_b=dict(near=near.get(k, []))))
    json.dump(dict(a=a, b=b, loci=out), open(f"{C.W}/compare_{a}.json", "w"), indent=1)
    print(f"run A = {a} (reference lacks the copy), run B = {b} (reference has it)")
    print("locus\tkind\tIsoCon outputs on the locus: A (flagged new) | B (flagged new)\trecovered by IsoCon / in-house candidate (A)\tin-house near the locus (B)")
    for k, v in out.items():
        i, h = v["isocon"], v["inhouse"]
        print(f"{k}\t{v['kind']}\t{i['run_a']['outputs']} ({i['run_a']['flagged_new']}) | {i['run_b']['outputs']} ({i['run_b']['flagged_new']})\t"
              f"{len(i['run_a']['recovering'])} / {len(h['run_a']['recovering'])}\t{len(h['run_b']['near'])}")


if __name__ == "__main__":
    main(sys.argv[1])
