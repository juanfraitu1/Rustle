#!/usr/bin/env python3
"""docs/archive/2026-10/merged/PREREG_o3_npip_maternal_copy_2026-10-21.md: the frozen reference-free chain (Amendments 7-9 of
PREREG_rna_allele_haplotype_count, control_test.py) on KB3781's NPIP reads with nothing masked, and its answer key from his
assemblies. Steps, in order (the IsoCon run and the minimap2 alignments are the prereg's commands, run between them):

    npip_o3.py net --w W                 W/net.fa: every read with an alignment record on one of the 25 NPIP copies, from its primary
                                         record in sequencing orientation (no cap, no unmapped reads); W/panel.json (control_test's)
    npip_o3.py frozen STEP --w W         control_test.STEP (outputs, contigs, lift, classify) on W
    npip_o3.py truth --w W               W/net.{mat,pat}.bam (his assemblies) -> W/target_reads.txt (R_T)
    npip_o3.py score --w W               candidates, their IsoCon reads and R_T -> bars E1-E5 -> W/score.json

docs/archive/2026-10/merged/PREREG_o3_npip_candidates_realign_2026-10-21.md (the counted candidates join the reference; the net is realigned):
    npip_o3.py realign-prep --w W        W/realign/C.fa (GGO + the candidates' contigs), W/realign/units.tsv (target reads' units)
    npip_o3.py score-realign --w W       W/realign/{R0,C}.bam -> bars H1-H3 and the tied-read counts -> W/realign/score.json
    npip_o3.py rel-prep --w W            Amendment 1: W/realign/rel.{a,b}.fa (the relatives' reads, from rel.names)
    npip_o3.py rel-score --w W           W/realign/{R0,C}.rel.bam -> bar H3' per set -> W/realign/rel_score.json
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
GGO = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
CAND = {"iso_NPIP_14": "B", "iso_NPIP_15": "B", "iso_NPIP_16": "C", "iso_NPIP_17": "C"}   # register 1363's counted candidates
UNITS = dict(A=(6720000, 6744300), B=(6744300, 6756100), C=(6766000, 6822000))           # his maternal chr18 (the audit's)
SLACK = 9                                                                                  # tied: within 10 AS (Eichler)
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


def unit_of(blocks):
    """The units (A, B, C of the maternal-only locus) a read's aligned blocks cover by >= 50 bases."""
    def cov(iv):
        return sum(min(e, iv[1]) - max(s, iv[0]) for s, e in blocks if s < iv[1] and iv[0] < e)
    return "".join(u for u in ("A", "B", "C") if cov(UNITS[u]) >= 50)


def realign_bars(units, calls, other):
    """H1-H3: units {target read: unit}, calls {read: "B" / "C" (a candidate) / "ref" / "unmapped"} in arm C, other {non-target
    read: call}."""
    c = [q for q, u in units.items() if u == "C"]
    b = [q for q, u in units.items() if u == "B"]
    hc, hb = sum(calls.get(q) == "C" for q in c), sum(calls.get(q) == "B" for q in b)
    fm = sum(v in ("B", "C") for v in other.values())
    out = dict(C_reads=len(c), C_home=hc, B_reads=len(b), B_home=hb, other=len(other), false_moves=fm,
               H1=bool(c) and hc >= 0.5 * len(c), H2=bool(b) and hb >= 0.5 * len(b), H3=fm <= 0.01 * len(other))
    out["verdict"] = "PASS" if out["H1"] and out["H2"] and out["H3"] else "FAIL"
    return out


def tied(best_as):
    """Whether a read's second-best locus ({locus: best AS}) lies within SLACK AS of its best."""
    v = sorted(best_as.values(), reverse=True)
    return len(v) >= 2 and v[0] - v[1] <= SLACK


def realign_prep(a):
    import pysam
    d = f"{a.w}/realign"
    os.makedirs(d, exist_ok=True)
    seq, cur = {}, None
    for ln in open(f"{a.w}/contigs_L.fa"):
        if ln[0] == ">":
            cur = ln[1:].strip(); seq[cur] = []
        else:
            seq[cur].append(ln.strip())
    with open(f"{d}/candidates.fa", "w") as o:
        for c in sorted(CAND):
            o.write(f">{c}\n{''.join(seq[c])}\n")
    subprocess.run(f"cat {GGO} {d}/candidates.fa > {d}/C.fa", shell=True, check=True)
    rt = {ln.strip() for ln in open(f"{a.w}/target_reads.txt") if ln.strip()}
    best = {}
    with pysam.AlignmentFile(f"{a.w}/net.mat.bam", check_sq=False) as bam:
        for r in bam.fetch(until_eof=True):
            if (r.query_name not in rt or r.is_unmapped or r.is_supplementary or r.reference_name != "chr18_mat_hsa16"
                    or not (r.reference_start < TARGET[2] and TARGET[1] < r.reference_end)):
                continue
            if r.query_name not in best or r.get_tag("ms") > best[r.query_name][0]:
                best[r.query_name] = (r.get_tag("ms"), r.get_blocks())
    with open(f"{d}/units.tsv", "w") as o:
        o.write("read\tunit\n")
        for q in sorted(rt):
            o.write(f"{q}\t{unit_of(best[q][1]) if q in best else ''}\n")
    n = collections.Counter(unit_of(v[1]) for v in best.values())
    print(f"[npip_o3] realign-prep: {d}/C.fa = GGO + {len(CAND)} contigs; target reads {len(rt)} by unit {dict(n)}", flush=True)


def _arm(bam_p, copies):
    """{read: call of its primary} and {read: {locus: best AS}} of one arm."""
    import pysam
    calls, loci = {}, collections.defaultdict(dict)
    with pysam.AlignmentFile(bam_p) as bam:
        for r in bam.fetch(until_eof=True):
            q = r.query_name
            if r.is_unmapped:
                calls.setdefault(q, "unmapped")
                continue
            if r.is_supplementary:
                continue
            ch = r.reference_name
            if ch in CAND:
                locus = CAND[ch]
            else:
                locus = next((g for c, s, e, g in copies if c == ch and s < r.reference_end and r.reference_start < e),
                             f"{ch}:{r.reference_start // 100000}")
            if not r.is_secondary:
                calls[q] = CAND.get(ch, "ref")
            if r.has_tag("AS"):
                loci[q][locus] = max(loci[q].get(locus, -10 ** 9), r.get_tag("AS"))
    return calls, loci


def score_realign(a):
    d = f"{a.w}/realign"
    units = {r["read"]: r["unit"] for r in csv.DictReader(open(f"{d}/units.tsv"), delimiter="\t")}
    copies = _copies()
    arms = {k: _arm(f"{d}/{k}.bam", copies) for k in ("R0", "C")}
    calls_c = arms["C"][0]
    other = {q: c for q, c in calls_c.items() if q not in units}
    b = realign_bars(units, {q: calls_c.get(q) for q in units}, other)
    r29 = {r["read"] for r in csv.DictReader(open(R1360), delimiter="\t") if r["family"] == "NPIP" and r["status"] == "off-copy"}
    by_unit = collections.defaultdict(collections.Counter)
    for q, u in units.items():
        by_unit[u or "none"][calls_c.get(q, "missing")] += 1
    b["calls_by_unit"] = {u: dict(v) for u, v in sorted(by_unit.items())}
    b["r1360_29_calls"] = dict(collections.Counter(calls_c.get(q, "missing") for q in r29))
    b["target_tied"] = {k: sum(tied(arms[k][1].get(q, {})) for q in units) for k in arms}
    b["other_calls_R0"] = dict(collections.Counter(c for q, c in arms["R0"][0].items() if q not in units))
    b["false_movers"] = sorted(q for q, c in other.items() if c in ("B", "C"))
    json.dump(b, open(f"{d}/score.json", "w"), indent=1)
    print("[npip_o3] " + ", ".join(f"{k} {b[k]}" for k in ("C_home", "C_reads", "B_home", "B_reads", "false_moves", "other", "H1",
                                                           "H2", "H3", "verdict")), flush=True)
    print(f"    by unit {b['calls_by_unit']}; register 1360's 29: {b['r1360_29_calls']}; target reads tied (within 10 AS) "
          f"R0 {b['target_tied']['R0']} -> C {b['target_tied']['C']}", flush=True)


def rel_tally(rows):
    """rows {read: (primary call, best AS, best AS on a candidate or None)} -> reads, moved (primary on a candidate), reachable
    (any candidate record), near-tied (a candidate alignment within SLACK AS of the read's best)."""
    n = collections.Counter(reads=len(rows))
    for call, best, cand in rows.values():
        n["moved"] += call in ("B", "C")
        n["reachable"] += cand is not None
        n["near_tied"] += cand is not None and best is not None and cand >= best - SLACK
    return dict(reads=n["reads"], moved=n["moved"], reachable=n["reachable"], near_tied=n["near_tied"])


def _rel_rows(bam_p):
    import pysam
    rows = {}
    with pysam.AlignmentFile(bam_p) as bam:
        for r in bam.fetch(until_eof=True):
            q = r.query_name
            call, best, cand = rows.get(q, ("unmapped", None, None))
            if r.is_unmapped or r.is_supplementary:
                rows.setdefault(q, (call, best, cand))
                continue
            a = r.get_tag("AS") if r.has_tag("AS") else None
            if not r.is_secondary:
                call = CAND.get(r.reference_name, "ref")
            if a is not None:
                best = a if best is None else max(best, a)
                if r.reference_name in CAND:
                    cand = a if cand is None else max(cand, a)
            rows[q] = (call, best, cand)
    return rows


def rel_prep(a):
    d = f"{a.w}/realign"
    subprocess.run(f"samtools view -@ 4 -F 2308 -N {d}/rel.names {BAM} -o {d}/rel.prim.sam", shell=True, check=True)
    recs = [ln.split("\t") for ln in open(f"{d}/rel.prim.sam")]
    half = (len(recs) + 1) // 2
    for k, part in (("a", recs[:half]), ("b", recs[half:])):
        with open(f"{d}/rel.{k}.fa", "w") as o:
            for f in part:
                o.write(f">{f[0]}\n{sequencing_orientation(f[9], int(f[1]) & 16 != 0)}\n")
    print(f"[npip_o3] rel-prep: {len(recs)} primaries -> {d}/rel.a.fa ({half}), rel.b.fa ({len(recs) - half})", flush=True)


def rel_score(a):
    d = f"{a.w}/realign"
    sets = {r["read"]: r["set"] for r in csv.DictReader(open(f"{d}/rel.sets.tsv"), delimiter="\t")}
    out = {}
    for arm in ("R0", "C"):
        rows = _rel_rows(f"{d}/{arm}.rel.bam")
        out[arm] = dict(all=rel_tally(rows))
        for name in ("rrn3", "c_rel"):
            out[arm][name] = rel_tally({q: v for q, v in rows.items() if name in sets.get(q, "")})
        if arm == "C":
            out["moved_reads"] = sorted(q for q, v in rows.items() if v[0] in ("B", "C"))
    tot = out["C"]["all"]["reads"]
    out["H3_prime"] = out["C"]["all"]["moved"] <= 0.01 * tot
    json.dump(out, open(f"{d}/rel_score.json", "w"), indent=1)
    for arm in ("R0", "C"):
        print(f"[npip_o3] {arm}: " + "; ".join(f"{k} {v}" for k, v in out[arm].items()), flush=True)
    print(f"[npip_o3] H3' (moved <= 1% of {tot}): {out['H3_prime']}", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("net", "frozen", "truth", "score", "realign-prep", "score-realign", "rel-prep", "rel-score"))
    ap.add_argument("step", nargs="?")
    ap.add_argument("--w", required=True)
    a = ap.parse_args(argv)
    {"net": net, "frozen": frozen, "truth": truth, "score": score, "realign-prep": realign_prep,
     "score-realign": score_realign, "rel-prep": rel_prep, "rel-score": rel_score}[a.cmd](a)


if __name__ == "__main__":
    main()
