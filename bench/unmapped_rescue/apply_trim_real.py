#!/usr/bin/env python3
"""Rule T (Amendment 13) on a real bed, descriptive: decisions, agreement with the genome-based finding, identity x coverage after the trim.
    apply_trim_real.py <bedA|bedH>      (run with /home/juanfra/miniforge3/bin/python or python3; needs polish/both.paf and polish/score.json from run_polish.py)"""
import collections
import csv
import json
import sys

import attribute as A
import run_polish as RP
import trim5g as T

W = "/mnt/linuxdisk/tmp/o3_rescue"


def library_gate():
    """the 5' clip signature of the real library's cleanly aligned survivor reads (linktest/R.bam, role S) and the gate (Amendment 15)"""
    import os
    import subprocess
    import libsig
    cache = f"{W}/library_signature.json"
    if os.path.exists(cache):
        j = json.load(open(cache))
        return j["gate"], j["p"], j["signature"]
    lab = {r["read"]: r for r in csv.DictReader(open("/mnt/linuxdisk/tmp/rna_allele/linktest/labels.tsv"), delimiter="\t")}
    surv = {n for n, r in lab.items() if r["role"] == "S"}
    p = subprocess.Popen("samtools view -F 2308 /mnt/linuxdisk/tmp/rna_allele/linktest/R.bam", shell=True, stdout=subprocess.PIPE, text=True)
    sig = libsig.signature(p.stdout, keep=lambda n: n in surv)
    ok, pv = libsig.gate(sig)
    json.dump(dict(gate=ok, p=pv, signature=sig), open(cache, "w"), indent=1)
    return ok, pv, sig


def main(bed, mode="t"):
    d = f"{W}/{bed}/registered"
    cons = RP.read_cons(f"{d}/cons.fa")
    hs = collections.defaultdict(list)
    allh = list(A.read_blastn(f"{d}/cons.blastn.tsv"))
    for q, t, a, b in allh:
        hs[q].append((q, t, a, b))
    att = A.attribute_cover(A.cover_scores(allh), 1.10)
    rec = RP.best_records(open(f"{d}/polish/both.paf").read().splitlines(), "orig.")
    sc = json.load(open(f"{d}/polish/score.json"))["clusters"]
    if mode == "t2":
        import subprocess
        paf = f"{d}/polish/cons.targets.paf"
        subprocess.run(f"minimap2 -c -x splice:hq -uf -N 5 -t 4 {W}/{bed}/targets.fa {d}/cons.fa > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
        fams = {n.split("|")[0]: att.get(n.split("|")[0], (None,))[0] for n in cons}
        pref2 = T.prefixes_from_paf(open(paf).read().splitlines(), fams)
    if mode in ("t3", "t4"):
        gate_ok, gate_p, sig = library_gate()
        print(f"library gate: {'OPEN' if gate_ok else 'closed'} p={gate_p:.3g}; clean reads {sig['reads']}, pure clips {sig['pure']}, G lengths {sig['g_lengths']}")
    if mode == "t4":
        import libsig
        import seeds as SD
        art = libsig.artifact_distribution(sig)
        pool = SD.read_fa(f"{W}/{bed}/pool.fa")
        clr = RP.cluster_reads(f"{W}/{bed}")
        jh, amounts, vs_t3, under, withclip = collections.Counter(), collections.Counter(), collections.Counter(), 0, 0
    over = collections.Counter()
    why = collections.Counter()
    agree = collections.Counter()
    ok_before = ok_after = n = 0
    for name, seq in cons.items():
        k = name.split("|")[0]
        if k not in sc or not sc[k]["orig"]["on"]:
            continue
        n += 1
        h = rec[k]
        qs, qe, qlen = h["qstart"], h["qend"], h["qlen"]
        genome_clip = qs if (0 < qs <= 3 and set(seq[:qs]) == {"G"}) else 0
        fam = att.get(k, (None,))[0]
        if mode == "t4":
            new, tl, jj = T.correct_cluster(seq, [pool[r] for r in clr[k]], art) if gate_ok else (seq, 0, None)
            reason = "trimmed" if tl else ("gate_closed" if not gate_ok else "none")
            jh[jj] += 1
            amounts[tl] += 1
            vs_t3[T.trim_leading_g(seq)[1] - tl] += 1
            if genome_clip:
                withclip += 1
                under += tl < genome_clip
        elif mode == "t3":
            new, tl, reason = T.trim_leading_g(seq) if gate_ok else (seq, 0, "gate_closed")
            if tl:
                over[tl - genome_clip] += 1
        elif fam is None:
            new, tl, reason = seq, 0, "abstained"
        else:
            new, tl, reason = T.decide(seq, pref2.get(k) if mode == "t2" else T.best_copy_prefix(hs.get(k, []), fam))
        why[reason] += 1
        agree[("genome G clip" if genome_clip else "no genome G clip", "T trims" if tl else "T does not trim")] += 1
        if genome_clip and fam is not None:
            agree[("attributed with a genome G clip", "T trims" if tl else "T does not trim")] += 1
        if genome_clip and not tl:
            agree[("missed", reason)] += 1
        ident = h["ident"]
        before = ident * (qe - qs) / qlen
        aligned_after = (qe - qs) - max(0, tl - qs)          # trimmed bases that the aligner had covered are no longer counted
        after = ident * aligned_after / (qlen - tl)
        ok_before += before >= 0.999
        ok_after += after >= 0.999
    print(f"{bed} [{mode}]: {n} consensus sequences on the erased copy; T decisions {dict(why)}")
    for kk, v in sorted(agree.items(), key=str):
        print("  ", kk, v)
    if mode == "t4":
        print(f"  j_hat (templated G's per cluster): {dict(sorted(jh.items(), key=str))}; bases removed: {dict(sorted(amounts.items()))}; T3 removes this many more than T4: {dict(sorted(vs_t3.items()))}")
        print(f"  consensus sequences with a genome G clip: {withclip}; T4 removes fewer bases than the unaligned clip in {under} ({under / max(1, withclip):.1%}; bar: <= 10%)")
    if mode == "t3":
        trims = sum(over.values())
        print("  trimmed run minus the unaligned genome G clip (bases of templated sequence lost), over the trims:", dict(sorted(over.items())),
              f"-> <= 1 base in {sum(v for k, v in over.items() if k <= 1)} of {trims}")
    print(f"  identity x coverage >= 0.999: before {ok_before}, after the trim {ok_after}")


if __name__ == "__main__":
    main(*sys.argv[1:])
