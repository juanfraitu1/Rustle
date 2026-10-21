#!/usr/bin/env python3
"""docs/PREREG_o3_npip_maternal_copy_2026-10-21.md: the frozen reference-free chain (Amendments 7-9 of
PREREG_rna_allele_haplotype_count, control_test.py) on KB3781's NPIP reads with nothing masked, and its answer key from his
assemblies. Steps, in order (the IsoCon run and the minimap2 alignments are the prereg's commands, run between them):

    npip_o3.py net --w W                 W/net.fa: every read with an alignment record on one of the 25 NPIP copies, from its primary
                                         record in sequencing orientation (no cap, no unmapped reads); W/panel.json (control_test's)
    npip_o3.py frozen STEP --w W         control_test.STEP (outputs, contigs, lift, classify) on W
    npip_o3.py truth --w W               W/net.{mat,pat}.bam (his assemblies) -> W/target_reads.txt (R_T)
    npip_o3.py score --w W               candidates, their IsoCon reads and R_T -> bars E1-E5 -> W/score.json
"""
import argparse
import collections
import csv
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BAM = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
GWFAM = "/mnt/linuxdisk/tmp/rna_allele/fibro_cat.copies.tsv"
TRUTH = "/mnt/linuxdisk/tmp/rna_allele"
R1360 = "/mnt/linuxdisk/tmp/o3_graph/mainstay/positional/truth_ms.tsv"   # the 29 NPIP off-copy reads of register 1360
TARGET = ("CM054600.2", 6720000, 6820000)
MARGIN, MAX_INTRON = 5, 100000
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def sequencing_orientation(seq, is_reverse):
    return seq.translate(COMP)[::-1] if is_reverse else seq


def _inside(a, target):
    return a[0] == target[0] and a[1] < target[2] and target[1] < a[2]


def in_target(alns, target, margin=MARGIN, max_intron=MAX_INTRON):
    """Whether a read's best alignment in his genome lies in the target by >= margin ms over its best outside it. alns:
    [(accession, start, end, ms, longest intron)]; an alignment with an intron >= max_intron is no placement."""
    kept = [a for a in alns if a[4] < max_intron]
    inside = [a[3] for a in kept if _inside(a, target)]
    outside = [a[3] for a in kept if not _inside(a, target)]
    if not inside:
        return False
    return not outside or max(inside) - max(outside) >= margin


def hit_in_target(hap_hit, target):
    """classify's hap_hit ("hap:accession:start-end score") overlapping the target."""
    if not hap_hit:
        return False
    _, acc, span = hap_hit.split()[0].split(":")
    s, e = (int(x) for x in span.split("-"))
    return _inside((acc, s, e), target)


def bars(cands, rt):
    """E1-E5 of the prereg. cands: [{candidate, n_contigs, cls, in_T, reads}]; rt: the target reads."""
    counted = [c for c in cands if c["n_contigs"] >= 2]
    hits = [c for c in counted if c["cls"] == "a_haplotype_only" and c["in_T"]]
    false = [c["candidate"] for c in counted if c["cls"] != "a_haplotype_only"]
    other = [c["candidate"] for c in counted if c["cls"] == "a_haplotype_only" and not c["in_T"]]
    b = dict(E1=len(hits) >= 1, E2=len(hits) == 1, E3=False, E4=False, E5=not false, held=0, target_reads=len(rt),
             candidate_reads=0, false=false, other_true=other)
    if len(hits) == 1:
        held = len(hits[0]["reads"] & rt)
        b.update(held=held, candidate_reads=len(hits[0]["reads"]), E3=bool(rt) and held >= 0.5 * len(rt),
                 E4=bool(hits[0]["reads"]) and held >= 0.9 * len(hits[0]["reads"]))
    b["verdict"] = "PASS" if all(b[k] for k in ("E1", "E2", "E3", "E4", "E5")) else "FAIL"
    return b


def _copies():
    return [(r["chrom"], int(r["start"]), int(r["end"]), f"NPIP:{r['copy_idx']}")
            for r in csv.DictReader(open(GWFAM), delimiter="\t") if r["family_id"] == "NPIP"]


def net(a):
    import pysam
    cps = _copies()
    os.makedirs(a.w, exist_ok=True)
    json.dump([dict(fam="NPIP", mask=list(cps[0]), keep=[list(c) for c in cps[1:]])], open(f"{a.w}/panel.json", "w"))
    names = set()
    with pysam.AlignmentFile(BAM) as bam:
        for c, s, e, _ in cps:
            for rec in bam.fetch(c, s, e):
                names.add(rec.query_name)
    with open(f"{a.w}/net.names", "w") as o:
        o.write("".join(f"{n}\n" for n in sorted(names)))
    subprocess.run(f"samtools view -@ 4 -F 2308 -N {a.w}/net.names {BAM} -o {a.w}/net.prim.sam", shell=True, check=True)
    n = 0
    with open(f"{a.w}/net.fa", "w") as o:
        for ln in open(f"{a.w}/net.prim.sam"):
            f = ln.split("\t")
            o.write(f">{f[0]}\n{sequencing_orientation(f[9], int(f[1]) & 16 != 0)}\n")
            n += 1
    print(f"[npip_o3] net: {len(cps)} copies, {len(names)} reads with a record on a copy, {n} primaries -> {a.w}/net.fa", flush=True)


def _control():
    sys.path.insert(0, HERE)
    import control_test
    return control_test


def frozen(a):
    getattr(_control(), a.step)(argparse.Namespace(w=a.w, l=a.w))


def _alias():
    """{(hap, chromosome number): accession} for the haplotype splice indexes' names chrN_<hap>_hsa*."""
    out = {}
    for h in ("mat", "pat"):
        for acc, num, _ in csv.reader(open(f"{TRUTH}/{h}.len.tsv"), delimiter="\t"):
            out[(h, num)] = acc
    return out


def truth(a):
    import pysam
    alias, alns = _alias(), collections.defaultdict(list)
    for h in ("mat", "pat"):
        with pysam.AlignmentFile(f"{a.w}/net.{h}.bam") as bam:
            for rec in bam.fetch(until_eof=True):
                if rec.is_unmapped or rec.is_supplementary or not rec.has_tag("ms"):
                    continue
                m = re.fullmatch(r"chr(\w+?)_(mat|pat)_hsa[^_]*", rec.reference_name)
                acc = alias[(m.group(2), m.group(1))] if m else rec.reference_name
                intron = max((n for op, n in rec.cigartuples if op == 3), default=0)
                alns[rec.query_name].append((acc, rec.reference_start, rec.reference_end, rec.get_tag("ms"), intron))
    rt = sorted(q for q, v in alns.items() if in_target(v, TARGET))
    with open(f"{a.w}/target_reads.txt", "w") as o:
        o.write("".join(f"{q}\n" for q in rt))
    r29 = {r["read"] for r in csv.DictReader(open(R1360), delimiter="\t") if r["family"] == "NPIP" and r["status"] == "off-copy"}
    print(f"[npip_o3] truth: {len(alns)} net reads aligned to his assemblies; target reads {len(rt)}; register 1360's {len(r29)} "
          f"off-copy NPIP reads among them {len(r29 & set(rt))}", flush=True)


def score(a):
    C = _control()
    rt = {ln.strip() for ln in open(f"{a.w}/target_reads.txt") if ln.strip()}
    comp, _ = C.comps_at(argparse.Namespace(w=a.w, l=a.w), C.DELTA)
    out_of = {r["contig"]: r["output"] for r in csv.DictReader(open(f"{a.w}/contigs.tsv"), delimiter="\t")}
    tid = lambda s: re.search(r"transcript_(\d+)_", s).group(1)
    reads_of = collections.defaultdict(set)
    for ln in open(f"{a.w}/iso/NPIP/cluster_info.tsv"):
        f = ln.rstrip("\n").split("\t")
        reads_of[tid(f[1])].add(f[0])
    members = collections.defaultdict(list)
    for c, cid in comp.items():
        members[cid].append(c)
    cands = []
    for r in csv.DictReader(open(f"{a.w}/candidates.tsv"), delimiter="\t"):
        rd = set().union(*(reads_of[tid(out_of[c])] for c in members[r["candidate"]]))
        cands.append(dict(candidate=r["candidate"], n_contigs=int(r["n_contigs"]), cls=r["class"], hap_hit=r["hap_hit"],
                          in_T=hit_in_target(r["hap_hit"], TARGET), reads=rd))
    b = bars(cands, rt)
    r29 = {r["read"] for r in csv.DictReader(open(R1360), delimiter="\t") if r["family"] == "NPIP" and r["status"] == "off-copy"}
    hit = [c for c in cands if c["n_contigs"] >= 2 and c["cls"] == "a_haplotype_only" and c["in_T"]]
    b["r1360_reads"] = len(r29)
    b["r1360_reads_held"] = len(r29 & hit[0]["reads"]) if len(hit) == 1 else None
    b["candidates"] = [dict({k: v for k, v in c.items() if k != "reads"}, reads=len(c["reads"]), target_reads=len(c["reads"] & rt))
                       for c in cands]
    json.dump(b, open(f"{a.w}/score.json", "w"), indent=1)
    for c in b["candidates"]:
        print(f"    {c['candidate']}: contigs {c['n_contigs']}, {c['cls']}, {c['hap_hit']}, in T {c['in_T']}, reads {c['reads']} "
              f"(target {c['target_reads']})", flush=True)
    print("[npip_o3] " + ", ".join(f"{k} {b[k]}" for k in ("E1", "E2", "E3", "E4", "E5", "held", "target_reads", "candidate_reads",
                                                           "false", "other_true", "r1360_reads", "r1360_reads_held", "verdict")),
          flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("net", "frozen", "truth", "score"))
    ap.add_argument("step", nargs="?")
    ap.add_argument("--w", required=True)
    a = ap.parse_args(argv)
    {"net": net, "frozen": frozen, "truth": truth, "score": score}[a.cmd](a)


if __name__ == "__main__":
    main()
