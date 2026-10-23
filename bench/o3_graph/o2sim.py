#!/usr/bin/env python3
"""Milestone 3: O2 (forced haplotype-to-copy matches) on simulated gorilla families with exact read origins.
Per family: every copy's dominant fibroblast isoform on both haplotypes (o3sim_one), N reads per allele, the family's mini reference;
case `control` (nothing erased, every copy in the table) and case `missing` (the copy with the most tied reads masked and left out).
The port (default options) runs; its AS-tied reads (>= 2 piles) are scored against their origin (o2eval_sim.py's rule), with
minimap2's primary as the baseline.
    o2sim.py OUT_DIR FAM[,FAM...] [N]"""
import collections, csv, json, os, random, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import o3sim_one as S
import pysam
out, fams = sys.argv[1], sys.argv[2].split(","); N = int(sys.argv[3]) if len(sys.argv) > 3 else 100
rows = S.tsv(S.SRC); by = collections.defaultdict(list)
for r in rows: by[r["fam"]].append(r)
cm = {r["pri"]: r for r in S.tsv(f"{S.RA}/chrmap.tsv")}
pri, bam = pysam.FastaFile(S.PRI), pysam.AlignmentFile(S.FIBRO)
tot = collections.Counter()
for f in fams:
    copies, tx = {}, {}
    for r in sorted(by[f], key=lambda r: (r["chrom"], int(r["clean_start"]))):
        iso = S.dominant_isoform(bam, r["chrom"], int(r["clean_start"]), int(r["clean_end"]))
        if not iso: continue
        bl, st, k, n = iso
        t = S.transcripts(pri, cm, r["chrom"], bl, st)
        copies[r["gene"]] = dict(chrom=r["chrom"], start=int(r["clean_start"]), end=int(r["clean_end"]))
        tx[f"{r['gene']}|A"] = t["A"]
        if t["B"]: tx[f"{r['gene']}|B"] = t["B"]
    for case in ("control", "missing"):
        D = f"{out}/{f}_{case}"; os.makedirs(D, exist_ok=True)
        E = None
        if case == "missing":
            c0 = f"{out}/{f}_control"
            o2 = list(csv.DictReader(open(f"{c0}/port.o2.tsv"), delimiter="\t"))
            cnt = collections.Counter(r["read"].split("|")[0] for r in o2)
            E = cnt.most_common(1)[0][0] if cnt else None
            if E is None: continue
        where = S.mini_reference(copies, E, case == "missing", D)
        S.simulate(tx, [(t, N) for t in tx], random.Random(1), f"{D}/reads.fa")
        S.align(D, f"{D}/reads.fa", f"{D}/reads.bam")
        table = [c for c in copies if c != E]
        with open(f"{D}/copies.tsv", "w") as fh:
            fh.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
            for k, c in enumerate(table):
                seg, s0, e0 = where[c]; fh.write(f"{f}\t{k}\t{c}\t{seg}\t{s0}\t{e0}\t1\t+\t0\t{s0}-{e0}\n")
        names = [l[1:].strip() for l in open(f"{D}/reads.fa") if l.startswith(">")]
        with open(f"{D}/net.tsv", "w") as fh:
            fh.write("family\tread\tused\trole\tvia\tclass\tin_reads_tsv\n")
            for n in names: fh.write(f"{f}\t{n}\t1\tS\tA\tNA\t0\n")
        subprocess.run([sys.executable, f"{HERE}/o3port.py", "run", "--bam", f"{D}/reads.bam", "--copies", f"{D}/copies.tsv",
                        "--genome", f"{D}/ref.fa", "--net", f"{D}/net.tsv", "--out", f"{D}/port"], check=True,
                       stderr=subprocess.DEVNULL)
        r = subprocess.run([sys.executable, f"{HERE}/o2eval_sim.py", D], capture_output=True, text=True).stdout.strip()
        print(r, flush=True)
