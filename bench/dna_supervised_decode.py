#!/usr/bin/env python
"""F3 — DNA-supervised copy decoding: assign RNA reads using the DNA-derived PSV signatures (SDA/Vollger idea:
DNA pre-phases PSVs->copies, turning the NP-hard unsupervised phasing into a supervised nearest-signature decode),
and measure NON-CIRCULAR accuracy.

Why: the O2 'silver standard' is circular (agreement with minimap2's own primary) and uses RNA-assembled copy
profiles. The DNA-derived PSV catalog (T2T reference) defines each copy's distinguishing alleles INDEPENDENTLY of
the RNA reads. So we can (1) build per-copy signatures from DNA, (2) decode RNA reads against them, and (3)
validate with a HELD-OUT-DNA-COLUMN cross-validation: split each read's DNA-defined PSV columns into a TRAIN half
(drives the call) and a TEST half (confirms it). High held-out agreement, with neither the copies nor the columns
nor the labels coming from the RNA, is non-circular evidence the reads genuinely carry the DNA copy signatures.

Method (per co-located family with a DNA catalog JSON):
  - load {fid}.json -> ref0 (chrom, iv=[rs,re]) + matrix[label][ref0_localpos]=base (DNA substitutions vs ref0).
  - exonic PSV columns = exonic ref0-local positions where >=2 copies differ (genomic pos = rs + localpos).
  - per-copy signature over those columns: ref0 carries the reference base; others carry matrix base or ref.
  - REALIGN the family's RNA reads to the ref0 sequence (minimap2) -> every read is in ref0's frame; read its
    allele at each PSV genomic position.
  - decode each read to the DNA copy whose signature it matches (copy_assign.assign_read LLR over observed cols).
  - HELD-OUT CV: split observed PSV cols even/odd; assign on TRAIN; check TEST cols rank the same DNA copy first.

Run: /home/juanfra/miniforge3/bin/python bench/dna_supervised_decode.py
"""
import os
import subprocess
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.dirname(__file__))
import dna_psv_catalog as dc          # copy_interval, extract_seq, exon_intervals, _exonic_localpos, parse_member
import copy_assign as ca              # assign_read (LLR significance gate)

DNADIR = "/home/juanfra/winloci_scratch/dna_catalog"
FAM_TSV = dc.FAM_TSV
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
MM2 = os.environ.get("RUSTLE_MINIMAP2", "/home/juanfra/miniforge3/bin/minimap2")
SAM = "/home/juanfra/miniforge3/bin/samtools"
COLOC_SPAN = 3_000_000      # treat a family as co-located if all members are within this on one chrom
MIN_PSV = 2                 # need >=2 exonic PSV cols to split train/test

# Cache all GFF exon intervals once (per chrom), replacing dna_psv_catalog's per-family awk scan (the bottleneck).
_EXON_CACHE = None
def _load_exons():
    global _EXON_CACHE
    if _EXON_CACHE is not None:
        return _EXON_CACHE
    import collections
    cache = collections.defaultdict(list)
    with open(dc.GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) > 4 and f[2] == "exon":
                cache[f[0]].append((int(f[3]) - 1, int(f[4])))
    for c in cache:
        cache[c].sort()
    _EXON_CACHE = cache
    return cache

def cached_exon_intervals(chrom, s, e):
    return [(a, b) for (a, b) in _load_exons().get(chrom, []) if a < e and b > s]


def colocated_families():
    """families from denovo_families.tsv whose members are all on one chrom within COLOC_SPAN (clean tandems)."""
    fams = []
    for line in open(FAM_TSV):
        if line.startswith("family_id"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 3:
            continue
        fid = p[0]
        members = [x for x in (dc.parse_member(m) for m in p[2].split(",")) if x]
        if not (2 <= len(members) <= 8):
            continue
        chroms = {c for c, _ in members}
        starts = [s for _, s in members]
        if len(chroms) == 1 and (max(starts) - min(starts)) <= COLOC_SPAN:
            fams.append((fid, members))
    return fams


def build_signatures(fid, members):
    """Return (chrom, ref_genome, psv_gpos[list], copy_vecs{label:{col:allele}}) or None if unusable.
    psv_gpos[i] = genomic position of PSV column i; copy_vecs over column index i."""
    jp = os.path.join(DNADIR, f"{fid}.json")
    if not os.path.exists(jp):
        return None
    import json
    J = json.load(open(jp))
    chrom, (rs, re_) = J["chrom"], J["iv"]
    matrix = J["matrix"]                       # {label: {localpos(str): base}}
    exons = cached_exon_intervals(chrom, rs, re_)
    exonic = dc._exonic_localpos(chrom, (rs, re_), exons)
    if not exonic:
        return None
    fa = pysam.FastaFile(GENOME)
    ref_seq = fa.fetch(chrom, rs, re_).upper()
    labels = ["ref0"] + list(matrix.keys())    # ref0 = implicit reference everywhere
    # candidate PSV cols = exonic local positions where some copy has a substitution
    cand = sorted({int(p) for lab in matrix for p in matrix[lab]} & exonic)
    psv_gpos, copy_vecs = [], {lab: {} for lab in labels}
    col = 0
    for p in cand:
        refb = ref_seq[p] if 0 <= p < len(ref_seq) else None
        if refb is None:
            continue
        alleles = {}
        for lab in labels:
            alleles[lab] = refb if lab == "ref0" else matrix[lab].get(str(p), refb)
        if len({a for a in alleles.values()}) < 2:    # not actually distinguishing
            continue
        for lab in labels:
            copy_vecs[lab][col] = alleles[lab]
        psv_gpos.append(rs + p)
        col += 1
    if col < MIN_PSV:
        return None
    return chrom, ref_seq, psv_gpos, copy_vecs


def read_alleles_via_ref0(fid, chrom, members, ref_seq, psv_gpos, rs):
    """Realign all family reads to ref0; return list of {col_idx: allele} (one per read, observed PSV cols)."""
    obs_list = []
    with tempfile.TemporaryDirectory() as td:
        ref0_fa = os.path.join(td, "ref0.fa")
        with open(ref0_fa, "w") as fh:
            fh.write(f">ref0\n{ref_seq}\n")
        reads_fa = os.path.join(td, "reads.fa")
        # union of member loci +-2kb
        regions = " ".join(f"{chrom}:{max(0,s-2000)}-{s+40000}" for _, s in members)
        seen = set()
        with open(reads_fa, "w") as out:
            r = subprocess.run(f"{SAM} view -F 0x900 {BAM} {regions}", shell=True, capture_output=True, text=True)
            for ln in r.stdout.splitlines():
                f = ln.split("\t")
                if len(f) < 10 or f[0] in seen:
                    continue
                seen.add(f[0])
                out.write(f">{f[0]}\n{f[9]}\n")
        if not seen:
            return obs_list
        bam = os.path.join(td, "r2.bam")
        subprocess.run(f"{MM2} -ax map-hifi --secondary=no -t 4 {ref0_fa} {reads_fa} 2>/dev/null "
                       f"| {SAM} sort -o {bam} - && {SAM} index {bam}", shell=True)
        gpos_to_col = {g: i for i, g in enumerate(psv_gpos)}
        local_psv = {g - rs: i for i, g in enumerate(psv_gpos)}   # ref0-local pos -> col idx
        af = pysam.AlignmentFile(bam, "rb")
        for aln in af.fetch():
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            obs = {}
            for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
                if rpos in local_psv and qpos is not None:
                    obs[local_psv[rpos]] = aln.query_sequence[qpos].upper()
            if obs:
                obs_list.append(obs)
    return obs_list


def argmax_set(obs, copy_vecs):
    _, _, _, _, logL = ca.assign_read(obs, copy_vecs)
    top = max(logL.values())
    return {n for n, v in logL.items() if abs(v - top) < 1e-9}


def main():
    fams = colocated_families()
    print(f"{len(fams)} co-located DNA-catalog families (<= {COLOC_SPAN//10**6}Mb span, 2-8 members)")
    tot_reads = tot_decoded = 0
    tot_cv = tot_confirm = 0
    chance_acc = 0.0
    nfam_used = 0
    per = []
    for fid, members in fams:
        sig = build_signatures(fid, members)
        if not sig:
            continue
        chrom, ref_seq, psv_gpos, copy_vecs = sig
        import json
        rs = json.load(open(os.path.join(DNADIR, f"{fid}.json")))["iv"][0]   # ref0 genomic start (iv[0])
        K = len(copy_vecs)
        obs_list = read_alleles_via_ref0(fid, chrom, members, ref_seq, psv_gpos, rs)
        if not obs_list:
            continue
        nfam_used += 1
        fam_cv = fam_confirm = fam_dec = 0
        for obs in obs_list:
            tot_reads += 1
            # decisive cols = cols where >=2 alleles among copies (all our cols qualify by construction)
            cols = sorted(obs)
            if len(cols) >= 1:
                best, margin, _, n_dec, _ = ca.assign_read(obs, copy_vecs)
                if n_dec >= 1 and margin >= ca.MARGIN:
                    tot_decoded += 1
                    fam_dec += 1
            # HELD-OUT CV: need >=2 observed cols
            if len(cols) >= 2:
                train = {c: obs[c] for i, c in enumerate(cols) if i % 2 == 0}
                test = {c: obs[c] for i, c in enumerate(cols) if i % 2 == 1}
                if train and test:
                    tot_cv += 1
                    fam_cv += 1
                    btr, _, _, _, _ = ca.assign_read(train, copy_vecs)
                    if btr in argmax_set(test, copy_vecs):
                        tot_confirm += 1
                        fam_confirm += 1
        chance_acc += (1.0 / K) * fam_cv
        if fam_cv:
            per.append((fid, K, fam_cv, fam_confirm, round(100 * fam_confirm / fam_cv, 1)))

    print(f"\nfamilies used: {nfam_used} | reads decoded (margin pass): {tot_decoded}/{tot_reads}")
    if tot_cv:
        chance = 100 * chance_acc / tot_cv
        conf = 100 * tot_confirm / tot_cv
        print(f"\nHELD-OUT-DNA-COLUMN cross-validation (DNA defines copies AND columns; no RNA self-ref, no silver):")
        print(f"  reads with >=2 DNA PSV cols: {tot_cv}")
        print(f"  held-out confirmation: {tot_confirm}/{tot_cv} = {conf:.1f}%")
        print(f"  weighted 1/K chance baseline: {chance:.1f}%   ->   enrichment {conf/chance:.1f}x")
        print(f"\n  => DNA-supervised decoding is corroborated by held-out DNA columns at {conf/chance:.1f}x chance,")
        print(f"     using a reference (copies + distinguishing columns) built entirely from DNA, independent of")
        print(f"     the RNA reads and of minimap2's primary placement (the circular silver standard).")
    print(f"\nper-family (top by reads): fid K cv confirm conf%")
    for r in sorted(per, key=lambda x: -x[2])[:20]:
        print(f"  {r[0]:>10} K={r[1]} cv={r[2]} confirm={r[3]} {r[4]}%")


if __name__ == "__main__":
    main()
