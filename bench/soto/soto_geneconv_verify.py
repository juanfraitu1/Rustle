#!/usr/bin/env python3
"""(a) Verify the gene-conversion signature on the 31 "exon-homogenized" floor members.

Claim being tested: these members' reads multi-map (exon sequence indistinguishable from a sibling) yet the
full genomic sequence is <99.9% identical -> the divergence is in INTRONS, exons were homogenized. If exons
are conserved while introns diverge, and the member is a PSEUDOGENE (no purifying selection possible), the
only explanation is gene conversion. Purifying selection is the coding-gene alternative.

Decisive metric per member B vs its nearest sibling A (both from soto_members.fa; exons from soto_reads.bam):
  identity_exon   = align exon-sum(B) vs exon-sum(A)   (read-defined exon blocks on the reference sequence)
  identity_genome = align B-genomic   vs A-genomic
  homogenized iff identity_exon >> identity_genome (exons far more conserved than the surrounding genome)

Verdict: gene-conversion (homogenized + pseudogene), homogenized-coding (could be conversion OR purifying
selection), or uniform (exon ~ genome divergence -> the 'exon-homogenized' label was wrong).

Run: /home/juanfra/miniforge3/bin/python bench/soto/soto_geneconv_verify.py
"""
import csv
import os
import re
import subprocess
from collections import defaultdict

import pysam

D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
MEMFA = f"{D}/soto_members.fa"
BAM = f"{D}/soto_reads.bam"
PAF = "bench/soto/a119b_member_pairs.paf"
DECOMP = "bench/soto/soto_floor_decomposition.tsv"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
SCRATCH = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/c40f2f02-e69d-4a8f-82ea-9ed7974a0cb9/scratchpad"
PSEUDO = re.compile(r"P\d*$")  # trailing pseudogene marker (…P, …P1, …P2)


def load_fa(path):
    seqs = {}; hdr = None; buf = []
    for line in open(path):
        if line.startswith(">"):
            if hdr:
                seqs[hdr] = "".join(buf)
            hdr = line[1:].strip(); buf = []
        else:
            buf.append(line.strip())
    if hdr:
        seqs[hdr] = "".join(buf)
    return seqs


def best_siblings():
    best = defaultdict(lambda: (0.0, None))
    for line in open(PAF):
        f = line.split("\t")
        if len(f) < 11:
            continue
        q, t, m, al = f[0], f[5], int(f[9]), int(f[10])
        idp = m / al if al else 0
        if idp > best[q][0]:
            best[q] = (idp, t)
    return best


def exon_sum(seqs, header, chrom, start, end):
    """concat of read-defined exon blocks of the member's reference sequence (spliced, transcription-agnostic)."""
    if header not in seqs:
        return None
    gseq = seqs[header]
    bam = pysam.AlignmentFile(BAM, "rb")
    ivs = []
    for a in bam.fetch(chrom, max(0, start), end):
        if a.is_unmapped or a.is_secondary or a.is_supplementary:
            continue
        for (bs, be) in a.get_blocks():
            s, e = max(bs, start) - start, min(be, end) - start
            if e > s:
                ivs.append((s, e))
    bam.close()
    if not ivs:
        return None
    ivs.sort(); merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return "".join(gseq[s:e] for s, e in merged if e <= len(gseq))


def mm2_id(qseq, tseq, tag):
    if not qseq or not tseq or len(qseq) < 50 or len(tseq) < 50:
        return None
    qf, tf = f"{SCRATCH}/gc_q_{tag}.fa", f"{SCRATCH}/gc_t_{tag}.fa"
    open(qf, "w").write(f">q\n{qseq}\n"); open(tf, "w").write(f">t\n{tseq}\n")
    out = subprocess.run([MM2, "-cx", "asm20", "-t", "2", tf, qf],
                         capture_output=True, text=True).stdout
    best = None
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        matches, alnlen = int(f[9]), int(f[10])
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        idp = (1 - de) if de is not None else matches / alnlen
        if best is None or alnlen > best[1]:
            best = (idp, alnlen)
    os.remove(qf); os.remove(tf)
    return best


def main():
    seqs = load_fa(MEMFA)
    sib = best_siblings()
    members = [r for r in csv.DictReader(open(DECOMP), delimiter="\t")
               if r["cause"].startswith("expressed-K=0: exon-homogenized")]
    rows = []
    tax = defaultdict(int)
    for i, r in enumerate(members):
        gene, chrom = r["gene"], r["chrom"]
        s, e = int(r["start"]), int(r["end"])
        hdr = f"{chrom}:{s+1}-{e}"
        sib_id, sib_hdr = sib.get(hdr, (0.0, None))
        if sib_hdr is None:
            tax["no-sibling-alignment"] += 1
            rows.append([r["family_id"], gene, chrom, s, e, "NA", "NA", "NA", "no-sibling"])
            continue
        # sibling genomic coords from its header chrom:start1-end
        m = re.match(r"(.+):(\d+)-(\d+)", sib_hdr)
        sc, ss, se = m.group(1), int(m.group(2)) - 1, int(m.group(3))
        ex_b = exon_sum(seqs, hdr, chrom, s, e)
        ex_a = exon_sum(seqs, sib_hdr, sc, ss, se)
        id_exon = mm2_id(ex_b, ex_a, f"e{i}")
        id_gen = mm2_id(seqs.get(hdr), seqs.get(sib_hdr), f"g{i}")
        ie = id_exon[0] if id_exon else None
        ig = id_gen[0] if id_gen else None
        pseudo = bool(PSEUDO.search(gene)) or gene.startswith("AC") or gene.startswith("AL")
        if ie is None or ig is None:
            verdict = "unmeasurable (exon-sum too short)"
        elif ie >= 0.990 and (ie - ig) >= 0.010:
            verdict = "GENE-CONVERSION (pseudogene: exons homogenized, introns diverged)" if pseudo \
                      else "homogenized-coding (conversion OR purifying selection)"
        elif ie >= 0.990 and ig >= 0.990:
            verdict = "young-identical (exon AND genome ~identical)"
        else:
            verdict = "uniform-divergence (exon ~ genome; 'homogenized' label suspect)"
        tax[verdict.split(" (")[0]] += 1
        rows.append([r["family_id"], gene, chrom, s, e,
                     f"{ie:.4f}" if ie is not None else "NA",
                     f"{ig:.4f}" if ig is not None else "NA",
                     f"{ie-ig:+.4f}" if (ie is not None and ig is not None) else "NA",
                     "pseudo" if pseudo else "coding", verdict])

    rows.sort(key=lambda x: x[-1])
    with open("bench/soto/soto_geneconv_verify.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family_id", "gene", "chrom", "start", "end", "identity_exon", "identity_genome",
                    "delta_exon_minus_genome", "kind", "verdict"])
        w.writerows(rows)

    print("\n=== (a) gene-conversion verification of 31 exon-homogenized floor members ===")
    for v, n in sorted(tax.items(), key=lambda kv: -kv[1]):
        print(f"  {n:3d}  {v}")
    print(f"  wrote bench/soto/soto_geneconv_verify.tsv")


if __name__ == "__main__":
    main()
