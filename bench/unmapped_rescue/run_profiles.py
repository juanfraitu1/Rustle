#!/usr/bin/env python3
"""HMM profile arm on bed H (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 5).   run_profiles.py   -> W/bedH/profiles/{families.hmm,hmm.tbl,result.json}"""
import collections
import csv
import json
import os
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import attribute as A  # noqa: E402
import profiles as P  # noqa: E402
import score as S  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
LT = "/mnt/linuxdisk/tmp/rna_allele/linktest"
PY = "/home/juanfra/miniforge3/bin/python"
HM = "/home/juanfra/miniforge3/envs/hmmer/bin"
MAFFT = "/home/juanfra/miniforge3/bin/mafft"


def sh(c):
    subprocess.run(c, shell=True, check=True)


def main():
    d = f"{W}/bedH"
    out = f"{d}/profiles"
    os.makedirs(out, exist_ok=True)
    panel = json.load(open(f"{LT}/panel.json"))
    # 1. survivor reads (primary placement overlapping the surviving copy), <= 100 longest
    if not os.path.exists(f"{out}/surv_reads.fa"):
        bam = pysam.AlignmentFile(f"{LT}/R.bam")
        with open(f"{out}/surv_reads.fa", "w") as o:
            for p in panel:
                for c, s, e, name in p["keep"]:
                    reads = {}
                    for rd in bam.fetch(c, s, e):
                        if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                            continue
                        reads[rd.query_name] = rd.get_forward_sequence()
                    for n, sq in sorted(reads.items(), key=lambda kv: -len(kv[1]))[:100]:
                        o.write(f">{name}|{n}\n{sq}\n")
    if not os.path.exists(f"{out}/surv_cons.fa"):
        sh(f"{PY} {HERE}/survivor_consensus.py {out}/surv_reads.fa {out}/surv_cons.fa")
    seqs, cur = collections.defaultdict(str), None
    for ln in open(f"{out}/surv_cons.fa"):
        if ln[0] == ">":
            cur = ln[1:].strip()
        else:
            seqs[cur] += ln.strip()
    fam_members = collections.defaultdict(list)
    for name in seqs:
        fam_members[A.family_of(name)].append(name)
    # 2. per-family MSA + hmmbuild
    if not os.path.exists(f"{out}/families.hmm"):
        open(f"{out}/families.hmm", "w").close()
        for fam, names in sorted(fam_members.items()):
            fa, afa, hmm = f"{out}/{fam}.fa", f"{out}/{fam}.afa", f"{out}/{fam}.hmm"
            with open(fa, "w") as o:
                for n in names:
                    o.write(f">{n}\n{seqs[n]}\n")
            if len(names) >= 2:
                sh(f"{MAFFT} --auto --quiet {fa} > {afa}")
            else:
                sh(f"cp {fa} {afa}")
            sh(f"{HM}/hmmbuild --dna --informat afa -n {fam} --cpu 4 {hmm} {afa} > /dev/null")
            sh(f"cat {hmm} >> {out}/families.hmm")
    # 3. nhmmer of every family HMM against the cluster consensus sequences
    if not os.path.exists(f"{out}/hmm.tbl"):
        sh(f"{HM}/nhmmer --tblout {out}/hmm.tbl -E 1e-3 --cpu 4 {out}/families.hmm {d}/registered/cons.fa > {out}/nhmmer.out")
    rows = [(t.split("|")[0],) + r[1:] for r in P.parse_tblout(open(f"{out}/hmm.tbl")) for t in (r[0],)]
    # 4. evaluate against the restricted baseline in the same universe
    lab = {x["read"]: x for x in csv.DictReader(open(f"{d}/labels.tsv"), delimiter="\t")}
    truth = {n: ("bg" if x["role"] == "bg" else (x["family"] or None)) for n, x in lab.items()}
    dread = {n for n, x in lab.items() if x["role"] == "D"}
    cl = {}
    for x in csv.DictReader(open(f"{d}/registered/clusters.tsv"), delimiter="\t"):
        cl.setdefault("cl" + x["cluster"], []).append(x["read"])
    universe = set(fam_members)
    hs = [h for h in A.read_blastn(f"{d}/registered/cons.blastn.tsv") if A.family_of(h[1]) in universe]
    surv = {p["fam"]: len(p["keep"]) for p in panel}
    arms = {"restricted_dc_megablast_cover": A.cover_scores(hs), "hmm_cover": P.cover_scores_from_rows(rows)}
    bits = P.bit_scores_from_rows(rows)

    def ev(sc, margin, subset=None):
        att = A.attribute_cover(sc, margin)
        sel = {c: rs for c, rs in cl.items() if subset is None or S.majority(rs, truth) in subset}
        m = S.rescue_metrics(sel, {c: (att[c][0] if c in att else None) for c in sel}, truth, {r for rs in sel.values() for r in rs if r in dread})
        m["wrong_fraction_of_joined"] = m["wrong_joins"] / max(1, m["rescued_correct"] + m["wrong_joins"])
        return m
    res = {"families_with_profile": len(universe), "survivor_consensus": len(seqs)}
    big = {f for f in universe if surv.get(f, 0) >= 3}
    for name, sc in arms.items():
        res[name] = ev(sc, 1.10)
        res[name + "_families_ge3_survivors"] = ev(sc, 1.10, big)
    res["hmm_bits_margin_1.10"] = ev(bits, 1.10)
    res["n_families_ge3_survivors"] = len(big)
    json.dump(res, open(f"{out}/result.json", "w"), indent=1)
    for k in ("restricted_dc_megablast_cover", "hmm_cover", "hmm_bits_margin_1.10", "restricted_dc_megablast_cover_families_ge3_survivors", "hmm_cover_families_ge3_survivors"):
        m = res[k]
        print(f"{k:55s} rescued {m['rescued_correct']:5d} wrong {m['wrong_joins']:4d} ({m['wrong_fraction_of_joined']:.4f}) attributed {m['clusters_attributed']:2d} abstained {m['clusters_abstained']:2d} copies {m['copies_reached']}/{m['copies_with_unmapped_reads']}")
    print("families with profile:", res["families_with_profile"], "survivor consensus sequences:", res["survivor_consensus"], "families >=3 survivors:", len(big))


if __name__ == "__main__":
    main()
