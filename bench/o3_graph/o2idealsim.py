#!/usr/bin/env python3
"""The frozen corrector on ideal-expression reads simulated from both of KB3781's haplotypes
(docs/archive/2026-10/merged/PREREG_o2_corrector_idealsim_2026-10-21.md; user: "ok lets do that").

    o2idealsim.py prepare GWFAM.tsv FAMILIES OUT_DIR
        the catalog's copies of FAMILIES (comma-separated) -> OUT_DIR/copies.tsv (o2psv.load_copies), families.json
        (donorsim's, set "ideal", one coding model per copy) and sel.json
    o2idealsim.py instances COPIES.tsv ORTH.tsv CONTESTED.tsv OUT.tsv
        donorsim's ORTHOLOGS rows: every copy's reference instance (hap "ref") and its counterpart on his other haplotype
        (o2personal.py orthologs --repaired), flagged contested when it overlapped another copy's before resolution
    o2idealsim.py truth TRUTH.tsv INSTANCES.tsv OUT_PREFIX
        donorsim's truth rows by stratum -> OUT_PREFIX.{all,ref,alt}.tsv; reads of contested instances left out
    o2idealsim.py bars OUT.json NAME=ALL_SCORE.json,REF_SCORE.json ...
        the prereg's bars per replicate (o2correct.py score JSONs) and the overall verdict
"""
import collections
import csv
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SET = "ideal"


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_is", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def family_set(rows, set_name=SET):
    """donorsim's FAMILIES.json for catalog copy rows (o2personal.gwfam_rows): one coding model per copy ("main"), the
    copies of each family listed by catalog index so donorsim names them FAMILY.INDEX as the copies table does."""
    fams = collections.defaultdict(dict)
    for r in rows:
        f, k = r["copy"].rsplit(".", 1)
        fams[f][int(k)] = dict(chrom=r["chrom"], start=r["start"], end=r["end"], strand=r["strand"],
                               transcripts={"main": [list(e) for e in r["exons"]]}, main="main", biotype="protein_coding")
    out = []
    for f in sorted(fams):
        ks = sorted(fams[f])
        assert ks == list(range(len(ks))), f"{f}: catalog indexes are not 0..n-1"
        out.append(dict(family=f, copies=[fams[f][k] for k in ks]))
    return {"sets": {set_name: out}}


def instances(rows, orth, contested):
    """Every copy's reference instance (its own exons on the reference, hap "ref") and its counterparts on his other
    haplotype; contested = {(copy, hap)} whose counterpart overlapped another copy's before resolution."""
    out = [dict(copy=r["copy"], hap="ref", contig=r["chrom"], start=r["start"], end=r["end"], strand="+", mapq=60,
                isoforms={"main": [list(e) for e in r["exons"]]}, contested=0) for r in rows]
    for o in orth:
        out.append(dict(copy=o["copy"], hap=o["hap"], contig=o["contig"], start=o["start"], end=o["end"], strand=o["strand"],
                        mapq=o["mapq"], isoforms={"main": [list(e) for e in o["exons"]]},
                        contested=int((o["copy"], o["hap"]) in contested)))
    return out


def strata(truth_rows, contested):
    """The scored truth rows: ref (reads of reference instances), alt (his other haplotype's, contested instances left
    out) and all (both)."""
    ref = [t for t in truth_rows if t["hap"] == "ref"]
    alt = [t for t in truth_rows if t["hap"] != "ref" and (t["truth"], t["hap"]) not in contested]
    return dict(all=ref + alt, ref=ref, alt=alt)


def replicate_bars(s_all, s_ref):
    """The prereg's bars for one replicate from o2correct.py score JSONs (all strata; the ref stratum)."""
    p, r = s_all["pooled"], s_ref["pooled"]
    reads = p.get("reads", 0)
    wrong = reads - p.get("minimap2 right", 0)
    fx, br = p.get("fixes", 0), p.get("breaks", 0)
    b = dict(P=reads >= 50 and wrong >= 10, R=fx > br and fx - br >= 0.02 * wrong, S=fx >= 5 * br,
             H=r.get("breaks", 0) <= r.get("fixes", 0) + 2,
             F=len(s_all["families_better"]) >= 2 and not s_all["families_worse"])
    b["verdict"] = "UNDERPOWERED" if not b["P"] else "PASS" if all(b[k] for k in "RSHF") else "FAIL"
    return b


def overall(verdicts):
    if verdicts and all(v == "PASS" for v in verdicts):
        return "PASS"
    return "FAIL" if "FAIL" in verdicts else "UNDERPOWERED"


def cmd_prepare(gwfam_p, families, out_dir):
    PC = _sibling("o2personal")
    keep = families.split(",")
    rows = [r for r in PC.gwfam_rows(gwfam_p) if r["family"] in keep]
    os.makedirs(out_dir, exist_ok=True)
    PC._write_copies(rows, f"{out_dir}/copies.tsv")
    json.dump(family_set(rows), open(f"{out_dir}/families.json", "w"), indent=1)
    json.dump(sorted(keep), open(f"{out_dir}/sel.json", "w"))
    n = collections.Counter(r["family"] for r in rows)
    print(f"[o2idealsim] prepare -> {out_dir}: " + ", ".join(f"{f} {k}" for f, k in sorted(n.items())), flush=True)
    return 0


def cmd_instances(copies_p, orth_p, contested_p, out_p):
    PC = _sibling("o2personal")
    con = {(r["copy"], r["hap"]) for r in csv.DictReader(open(contested_p), delimiter="\t")}
    inst = instances(PC.read_copies(copies_p), PC.read_orthologs(orth_p), con)
    with open(out_p, "w") as fh:
        fh.write("copy\thap\tcontig\tstart\tend\tstrand\tmapq\tisoforms\tcontested\n")
        for i in inst:
            fh.write(f"{i['copy']}\t{i['hap']}\t{i['contig']}\t{i['start']}\t{i['end']}\t{i['strand']}\t{i['mapq']}\t"
                     f"{json.dumps(i['isoforms'], separators=(',', ':'))}\t{i['contested']}\n")
    n = collections.Counter((i["hap"] if i["hap"] == "ref" else "alt") + (" contested" if i["contested"] else "") for i in inst)
    print(f"[o2idealsim] instances -> {out_p}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)
    return 0


def cmd_truth(truth_p, inst_p, out_prefix):
    con = {(r["copy"], r["hap"]) for r in csv.DictReader(open(inst_p), delimiter="\t") if r["contested"] == "1"}
    rd = csv.DictReader(open(truth_p), delimiter="\t")
    rows = list(rd)
    for name, sub in strata(rows, con).items():
        with open(f"{out_prefix}.{name}.tsv", "w") as fh:
            w = csv.DictWriter(fh, fieldnames=rd.fieldnames, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(sub)
    s = strata(rows, con)
    print(f"[o2idealsim] truth -> {out_prefix}.{{all,ref,alt}}.tsv: reads {len(rows)}; ref {len(s['ref'])}, alt {len(s['alt'])}, "
          f"contested left out {len(rows) - len(s['all'])}", flush=True)
    return 0


def cmd_bars(out_p, reps):
    res = {}
    for spec in reps:
        name, paths = spec.split("=", 1)
        a, r = paths.split(",")
        res[name] = replicate_bars(json.load(open(a)), json.load(open(r)))
        print(f"[o2idealsim] {name}: " + ", ".join(f"{k} {v}" for k, v in res[name].items()), flush=True)
    res["overall"] = overall([v["verdict"] for v in res.values()])
    print(f"[o2idealsim] overall: {res['overall']}", flush=True)
    json.dump(res, open(out_p, "w"), indent=1)
    return 0


def main(argv):
    if argv[:1] == ["prepare"] and len(argv) == 4:
        return cmd_prepare(*argv[1:])
    if argv[:1] == ["instances"] and len(argv) == 5:
        return cmd_instances(*argv[1:])
    if argv[:1] == ["truth"] and len(argv) == 4:
        return cmd_truth(*argv[1:])
    if argv[:1] == ["bars"] and len(argv) >= 3:
        return cmd_bars(argv[1], argv[2:])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
