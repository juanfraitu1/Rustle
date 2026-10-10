#!/usr/bin/env python3
"""O3 deliverable: the flagged reference-absent candidates of a run as unplaced transcripts, to load beside the assembler's GTF (docs/PREREG_unmapped_rescue_2026-10-08.md
Amendment 48: O3 = detection; placement is DNA's job).

    o3_export.py <run_dir> <out_prefix>      writes <out_prefix>.fa, .gtf, .nearest.bed, .tsv (resumable read-share step; exit 75 = run again)

Each candidate is a single-exon transcript on its own sequence (the consensus, in transcript orientation) with attributes: class, supporting reads, identity x coverage
on the reference, nearest reference locus (where its reads align on the reference), share of that locus' reads it carries (Amendment 44), placement "requires_WGS".
The .nearest.bed puts each candidate at its nearest reference locus for a genome browser. Pure functions are tested in test_o3_export.py."""
import json
import os
import subprocess
import sys
import time

FLAG = ("DIVERGED", "NOVEL", "ELSEWHERE")
SOURCE = "rustle_o3"


def candidate_id(i):
    return f"O3cand_{i:06d}"


def flagged(rows):
    """the rows a run flags (supported, not held by the reference within the allele cutoff), in run order"""
    return [r for r in rows if r["cls"] in FLAG]


def _f(v, nd=4):
    return "NA" if v is None else f"{v:.{nd}f}"


def _locus(n):
    return "NA" if not n else f"{n[0]}:{n[1] + 1}-{n[2]}({n[3]})"


def gtf_lines(c):
    """transcript + exon lines of one candidate on its own sequence (1-based, inclusive)"""
    L = len(c["seq"])
    attrs = (f'gene_id "{c["id"]}"; transcript_id "{c["id"]}.t1"; o3_class "{c["cls"]}"; supporting_reads "{c["reads"]}"; '
             f'ref_identity_x_coverage "{_f(c.get("R"))}"; median_read_divergence "{_f(c.get("median_read_divergence"))}"; '
             f'nearest_ref_locus "{_locus(c.get("nearest"))}"; locus_read_share "{_f(c.get("read_share"))}"; placement "requires_WGS";')
    score = str(min(1000, int(c["reads"])))
    t = "\t".join([c["id"], SOURCE, "transcript", "1", str(L), score, "+", ".", attrs])
    e = "\t".join([c["id"], SOURCE, "exon", "1", str(L), score, "+", ".", attrs + ' exon_number "1";'])
    return [t, e]


def bed_line(c):
    n = c.get("nearest")
    if not n:
        return None
    return "\t".join([n[0], str(n[1]), str(n[2]), c["id"], str(min(1000, int(c["reads"]))), n[3]])


def fasta_record(c):
    return (f">{c['id']} class={c['cls']} reads={c['reads']} ref_idcov={_f(c.get('R'))} nearest={_locus(c.get('nearest'))} placement=requires_WGS\n"
            + "\n".join(c["seq"][i:i + 80] for i in range(0, len(c["seq"]), 80)) + "\n")


def main(run_dir, out, budget=520, cap=300):
    import pysam
    HERE = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, HERE)
    import discover as D
    import locus_support as L
    import mat_truth as T
    import run_polish as RP
    import seeds as SD
    t0 = time.time()
    rows = json.load(open(f"{run_dir}/classes.json"))
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{run_dir}/cons.fa").items()}
    best = {}
    for ln in open(f"{run_dir}/cons.R.paf"):
        f = ln.rstrip("\n").split("\t")
        if "tp:A:P" in f[12:] and (f[0] not in best or int(f[9]) > int(best[f[0]][9])):
            best[f[0]] = f
    fl = flagged(rows)
    # read share at the nearest locus (resumable)
    sp = f"{out}.share.jsonl"
    share = {}
    if os.path.exists(sp):
        for ln in open(sp):
            k, v = json.loads(ln)
            share[k] = v
    if len(share) < len(fl):
        bam = pysam.AlignmentFile(L.BAM, "rb")
        pri = pysam.FastaFile(D.PRIMARY_FA)
        tmp = f"{out}.tmp"
        os.makedirs(tmp, exist_ok=True)
        RCT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        with open(sp, "a") as o:
            for r in fl:
                k = r["k"]
                if k in share:
                    continue
                if time.time() - t0 > budget:
                    print(f"read share: paused at {len(share)} of {len(fl)}; run again")
                    sys.exit(75)
                v = None
                f = best.get(k)
                if f:
                    cg = next(t[5:] for t in f[12:] if t.startswith("cg:Z:"))
                    pv = T.assembly_version(lambda a, b, c=f[5]: pri.fetch(c, a, b).upper(), int(f[7]), cg)
                    pv = pv.translate(RCT)[::-1] if f[4] == "-" else pv
                    names, seqs = [], {}
                    for a in bam.fetch(f[5], int(f[7]), int(f[8])):
                        if a.is_secondary or a.is_supplementary or a.query_sequence is None:
                            continue
                        q = a.query_sequence
                        seqs[a.query_name] = q.translate(RCT)[::-1] if a.is_reverse else q
                        names.append(a.query_name)
                    SD.write_fa(f"{tmp}/r.fa", seqs, L.sample_every(names, cap))
                    SD.write_fa(f"{tmp}/t.fa", {"cons": cons[k], "ref": pv}, ["cons", "ref"])
                    lines = subprocess.run(f"minimap2 -c -x splice:hq -uf --secondary=no -t 4 {tmp}/t.fa {tmp}/r.fa", shell=True, stdout=subprocess.PIPE,
                                           stderr=subprocess.DEVNULL, text=True).stdout.splitlines()
                    v = L.share(*L.assign(lines, "cons", "ref"))
                o.write(json.dumps([k, v]) + "\n")
                o.flush()
                share[k] = v
    cands = []
    for i, r in enumerate(fl, 1):
        f = best.get(r["k"])
        cands.append(dict(id=candidate_id(i), seq=cons[r["k"]], cls=r["cls"], reads=r["reads"], R=r["R"], median_read_divergence=r["median_read_divergence"],
                          nearest=(f[5], int(f[7]), int(f[8]), f[4]) if f else None, read_share=share.get(r["k"]), run_cluster=r["k"]))
    with open(f"{out}.fa", "w") as fa, open(f"{out}.gtf", "w") as gtf, open(f"{out}.nearest.bed", "w") as bed, open(f"{out}.tsv", "w") as tsv:
        gtf.write(f"##description: O3 reference-absent candidates (unplaced; placement requires WGS); source {run_dir}\n")
        tsv.write("id\tclass\treads\tref_idcov\tmedian_read_divergence\tnearest_ref_locus\tlocus_read_share\tlength\trun_cluster\n")
        for c in cands:
            fa.write(fasta_record(c))
            gtf.write("\n".join(gtf_lines(c)) + "\n")
            b = bed_line(c)
            if b:
                bed.write(b + "\n")
            tsv.write(f"{c['id']}\t{c['cls']}\t{c['reads']}\t{_f(c['R'])}\t{_f(c['median_read_divergence'])}\t{_locus(c['nearest'])}\t{_f(c['read_share'])}\t"
                      f"{len(c['seq'])}\t{c['run_cluster']}\n")
    print(f"{len(cands)} candidates written to {out}.{{fa,gtf,nearest.bed,tsv}}; with a nearest locus {sum(c['nearest'] is not None for c in cands)}; "
          f"read share >= 0.02: {sum((c['read_share'] or 0) >= 0.02 for c in cands)}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
