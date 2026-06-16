#!/usr/bin/env python3
"""
Read-coherence + PSV copy-split on REAL CO-LOCATED / TANDEM paralog families in the
GGO long-read data.  Adapted from bench/copy_split_realdata.py (the RABL2 driver):
same algorithm port of `split_readchain_by_psv` (src/rustle/vg_family/copy_split.rs),
same CIGAR/allele extraction, same PairwiseAligner PSV discovery.

Difference vs the RABL2 driver: RABL2's two copies sit on DIFFERENT contigs (separate
loci, each with their own primary reads).  Here the copies are TANDEM on the SAME contig
(NC_073235.2), so we

  - generalise COPIES to N copies (not 2),
  - discover PSVs from an all-pairs alignment of the N copy sequences (union of columns),
  - pull reads from the WHOLE tandem span once (so a read is counted once, not per copy),
    assigning each alignment to the copy whose genomic span it overlaps most,
  - flag the CO-LOCATED signature: do the copies share an intron chain? (the thing RABL2
    does NOT do, since RABL2A/B chains live on different contigs).

Two real families (autosomal NC_073235.2, ~100x library baseline but these families are
lowly expressed; DAZ/RBMY are unusable here per the task note):
  - LOC101144552 cluster : multi-exon paralog copies LOC129532018 / LOC129532017 / LOC101144552
  - LOC101132628 cluster : LOC101132628 / LOC115933736 / LOC115933737

Deterministic.  Usage:  python bench/copy_split_colocated.py
"""

import subprocess
import statistics
from collections import defaultdict, Counter

import pysam
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

# ----------------------------------------------------------------------------- inputs
SAMTOOLS = "/home/juanfra/miniforge3/bin/samtools"
BAM = "/home/juanfra/winloci_scratch/GGO.bam"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
CONTIG = "NC_073235.2"

# Families.  Each copy = a distinct annotated paralog gene (genomic span, strand) inside the
# tandem cluster.  Spans taken from ref.gtf (exon-derived, 1-based inclusive).
FAMILIES = {
    "LOC101144552": {
        "region": (144691596, 144756221),  # task-given tandem region
        "copies": [
            # name, start1, end1, strand   (the multi-exon protein-coding paralogs in the cluster)
            ("LOC129532018", 144622575, 144683955, "+"),
            ("LOC129532017", 144691596, 144721773, "+"),
            ("LOC101144552", 144779454, 144823106, "+"),
        ],
    },
    "LOC101132628": {
        "region": (135066522, 135115177),
        "copies": [
            ("LOC101132628", 135066522, 135069510, "+"),
            ("LOC115933736", 135094425, 135098162, "+"),
            ("LOC115933737", 135111301, 135115177, "+"),
        ],
    },
}

K_SWEEP = [2, 3]
MIN_READS_SWEEP = [2]


# ----------------------------------------------------------------------------- step 1
def faidx_seq(contig, start1, end1):
    region = f"{contig}:{start1}-{end1}"
    out = subprocess.run(
        [SAMTOOLS, "faidx", FASTA, region],
        capture_output=True, text=True, check=True,
    ).stdout
    seq = "".join(l.strip() for l in out.splitlines() if not l.startswith(">"))
    return seq.upper()


def extract_copy_seqs(copies):
    """Return {name: forward-oriented (transcription strand) sequence}."""
    seqs = {}
    for (name, s1, e1, strand) in copies:
        s = faidx_seq(CONTIG, s1, e1)
        if strand == "-":
            s = str(Seq(s).reverse_complement())
        seqs[name] = s
    return seqs


def pairwise_align(seq_a, seq_b):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    aln = aligner.align(seq_a, seq_b)[0]
    return aln


def copy_identity(aln, seq_a, seq_b):
    """% identity over aligned (ungapped) columns."""
    blocks_a, blocks_b = aln.aligned
    same = tot = 0
    for (a0, a1), (b0, b1) in zip(blocks_a, blocks_b):
        for i in range(a1 - a0):
            ba, bb = seq_a[a0 + i], seq_b[b0 + i]
            if ba in "ACGT" and bb in "ACGT":
                tot += 1
                if ba == bb:
                    same += 1
    return (100.0 * same / tot if tot else 0.0), tot


def offset_to_genomic(copy, off):
    """Map a forward-oriented-sequence offset to a 0-based genomic coord + the base AS IT
    APPEARS ON THE FORWARD GENOME STRAND (since we read alleles off the genomic pileup)."""
    name, s1, e1, strand = copy
    length = e1 - s1 + 1
    if strand == "+":
        return (s1 - 1) + off
    else:
        return (s1 - 1) + (length - 1 - off)


def discover_psvs(copies, seqs):
    """All-pairs PSV discovery.  A PSV column = a genomic position on the REFERENCE copy
    (copy[0]) where at least one other copy differs at the aligned base.  We anchor PSV
    columns on the reference copy's coordinate and record, per copy, the genomic position +
    forward-genome-strand base.  Returns:
       psv_cols: list of dict { copyname -> {gpos0, genome_base} }  (None if a copy has a gap)
    plus per-copy %identity to the reference.
    """
    ref = copies[0]
    ref_seq = seqs[ref[0]]
    # ref_off -> {copyname: (off_in_copy)}; only positions aligned in >=1 pair are kept
    ref_off_map = defaultdict(dict)  # ref_off -> {copyname: copy_off}
    idents = {}
    for other in copies[1:]:
        oseq = seqs[other[0]]
        aln = pairwise_align(ref_seq, oseq)
        idents[other[0]], _ = copy_identity(aln, ref_seq, oseq)
        ba, bb = aln.aligned
        for (a0, a1), (b0, b1) in zip(ba, bb):
            for i in range(a1 - a0):
                ai, bi = a0 + i, b0 + i
                if ref_seq[ai] in "ACGT" and oseq[bi] in "ACGT" and ref_seq[ai] != oseq[bi]:
                    ref_off_map[ai][other[0]] = bi

    psv_cols = []
    for ref_off in sorted(ref_off_map.keys()):
        rec = {}
        # reference copy
        rec[ref[0]] = dict(
            gpos0=offset_to_genomic(ref, ref_off),
            genome_base=(ref_seq[ref_off] if ref[3] == "+"
                         else str(Seq(ref_seq[ref_off]).reverse_complement())),
        )
        for other in copies[1:]:
            if other[0] in ref_off_map[ref_off]:
                ooff = ref_off_map[ref_off][other[0]]
                ob = seqs[other[0]][ooff]
                rec[other[0]] = dict(
                    gpos0=offset_to_genomic(other, ooff),
                    genome_base=(ob if other[3] == "+"
                                 else str(Seq(ob).reverse_complement())),
                )
            else:
                rec[other[0]] = None
        psv_cols.append(rec)
    return psv_cols, idents


# ----------------------------------------------------------------------------- step 2
def cigar_introns_and_bases(aln, want_positions):
    """Walk one alignment; return (introns[list of (d0,a0)], base_at{gpos0->base})."""
    ref_pos = aln.reference_start
    introns = []
    base_at = {}
    seq = aln.query_sequence or ""
    qpos = 0
    for op, length in aln.cigartuples or []:
        if op in (0, 7, 8):
            for k in range(length):
                g = ref_pos + k
                if g in want_positions:
                    base_at[g] = seq[qpos + k] if (qpos + k) < len(seq) else None
            ref_pos += length
            qpos += length
        elif op == 1:
            qpos += length
        elif op == 2:
            ref_pos += length
        elif op == 3:
            introns.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op == 4:
            qpos += length
        elif op == 5:
            pass
        elif op == 6:
            pass
    return introns, base_at


def assign_copy(aln, copies):
    """Assign an alignment to the copy whose genomic span it overlaps most (by bp)."""
    rs, re_ = aln.reference_start, aln.reference_end or aln.reference_start
    best, bestov = None, 0
    for (name, s1, e1, strand) in copies:
        ov = max(0, min(re_, e1) - max(rs, s1 - 1))
        if ov > bestov:
            bestov, best = ov, name
    return best


def pull_reads(family, psv_cols):
    """Pull every alignment across the WHOLE tandem cluster span (incl secondary), once.
    Build a global allele vector (one column per PSV) by reading the base at whichever
    copy-position falls under this alignment."""
    copies = family["copies"]
    span_lo = min(c[1] for c in copies)
    span_hi = max(c[2] for c in copies)
    n_cols = len(psv_cols)

    # genomic position -> column index (a column may have a position under several copies)
    want = {}
    for col, rec in enumerate(psv_cols):
        for name, d in rec.items():
            if d is not None:
                want.setdefault(d["gpos0"], []).append(col)

    bam = pysam.AlignmentFile(BAM, "rb")
    reads = []
    seen = set()
    for aln in bam.fetch(CONTIG, span_lo - 1, span_hi):
        if aln.is_unmapped:
            continue
        key = (aln.query_name, aln.reference_start, aln.flag)
        if key in seen:
            continue
        seen.add(key)
        introns, base_at = cigar_introns_and_bases(aln, set(want.keys()))
        psv_alleles = [None] * n_cols
        for gpos0, b in base_at.items():
            if b is None:
                continue
            for col in want[gpos0]:
                # only set if this column actually expects this copy's position; but since
                # different copies' PSV positions are distinct genomic coords, a read over
                # one copy contributes only to that copy's coordinate's column.
                psv_alleles[col] = ord(b.upper())
        reads.append(dict(
            name=aln.query_name,
            assigned_copy=assign_copy(aln, copies),
            is_secondary=aln.is_secondary,
            is_supplementary=aln.is_supplementary,
            strand="-" if aln.is_reverse else "+",
            mapq=aln.mapping_quality,
            ref_start=aln.reference_start,
            intron_chain=tuple(introns),
            psv_alleles=psv_alleles,
        ))
    bam.close()
    return reads


# ----------------------------------------------------------------------------- step 3 (port)
def read_consistent_with(read, copy):
    for r, cc in zip(read, copy):
        if r is not None and cc is not None:
            if r != cc:
                return False
        elif r is not None and cc is None:
            return False
    return True


def distinct_columns(a, b):
    return sum(1 for x, y in zip(a, b) if x is not None and y is not None and x != y)


def select_identifiable_copies(candidates, min_psv_k):
    chosen = []
    for cand in candidates:
        if all(distinct_columns(cand, c) >= min_psv_k for c in chosen):
            chosen.append(cand)
    return chosen


def consensus_haplotype(group, n_cols):
    out = []
    for col in range(n_cols):
        tally = Counter()
        for r in group:
            b = r["psv_alleles"][col]
            if b is not None:
                tally[b] += 1
        out.append(sorted(tally.items(), key=lambda kv: (-kv[1], kv[0]))[0][0] if tally else None)
    return out


def _vec_key(v):
    return tuple(-1 if b is None else b for b in v)


def split_readchain_by_psv(reads, n_psv_columns, min_psv_k, min_reads_per_copy):
    by_chain = defaultdict(list)
    for r in reads:
        by_chain[r["intron_chain"]].append(r)
    out = []
    for chain in sorted(by_chain.keys()):
        group = by_chain[chain]
        vector_counts = Counter()
        for r in group:
            if all(a is None for a in r["psv_alleles"]):
                continue
            vector_counts[tuple(r["psv_alleles"])] += 1
        candidates = [list(v) for v, cnt in sorted(vector_counts.items(), key=lambda kv: _vec_key(kv[0]))
                      if cnt >= min_reads_per_copy]
        split_copies = (select_identifiable_copies(candidates, min_psv_k)
                        if n_psv_columns >= min_psv_k else [])
        if len(split_copies) >= 2:
            counts = [0] * len(split_copies)
            assigned_reads = [[] for _ in split_copies]
            for r in group:
                consistent, unique = None, True
                for i, copy in enumerate(split_copies):
                    if read_consistent_with(r["psv_alleles"], copy):
                        if consistent is not None:
                            unique = False
                            break
                        consistent = i
                if unique and consistent is not None:
                    counts[consistent] += 1
                    assigned_reads[consistent].append(r)
            for copy, count, ar in zip(split_copies, counts, assigned_reads):
                out.append(dict(intron_chain=chain, allele_vector=copy,
                                read_count=count, identifiable=True, reads=ar))
        else:
            out.append(dict(intron_chain=chain,
                            allele_vector=consensus_haplotype(group, n_psv_columns),
                            read_count=len(group), identifiable=False, reads=list(group)))
    out.sort(key=lambda d: (d["intron_chain"], _vec_key(d["allele_vector"])))
    return out


# ----------------------------------------------------------------------------- reporting
def vec_str(v):
    return "".join("-" if b is None else chr(b) for b in v)


def n_observed(r):
    return sum(1 for a in r["psv_alleles"] if a is not None)


def run_family(fam_name, family):
    print("=" * 78)
    print(f"FAMILY {fam_name}   tandem region {CONTIG}:{family['region'][0]}-{family['region'][1]}")
    print("=" * 78)
    copies = family["copies"]
    print(f"\n[STEP 1] copies ({len(copies)}):")
    for (n, s, e, st) in copies:
        print(f"  {n:16s} {s}-{e} ({e-s+1}bp) strand {st}")

    seqs = extract_copy_seqs(copies)
    psv_cols, idents = discover_psvs(copies, seqs)
    print(f"\n  pairwise identity to reference copy {copies[0][0]}:")
    for n, pid in idents.items():
        print(f"    {n:16s} {pid:.2f}%")
    n_cols = len(psv_cols)
    print(f"  PSV columns (positions where >=1 copy differs from ref): {n_cols}")

    # PSV spacing on the reference copy
    refpos = sorted(rec[copies[0][0]]["gpos0"] for rec in psv_cols)
    if len(refpos) >= 2:
        sp = [b - a for a, b in zip(refpos, refpos[1:])]
        print(f"  PSV spacing on ref copy (bp): min={min(sp)} median={int(statistics.median(sp))} "
              f"max={max(sp)} mean={sum(sp)/len(sp):.0f}")

    print(f"\n[STEP 2] reads over whole cluster span (incl secondary):")
    reads = pull_reads(family, psv_cols)
    by_copy = Counter(r["assigned_copy"] for r in reads)
    cls = Counter(("sec" if r["is_secondary"] else ("supp" if r["is_supplementary"] else "pri")) for r in reads)
    print(f"  total alignments={len(reads)}   class={dict(cls)}")
    print(f"  per-copy assignment: {dict(by_copy)}")
    # distinct molecules with a primary in this cluster
    prim_names = set(r["name"] for r in reads if not r["is_secondary"] and not r["is_supplementary"])
    print(f"  distinct molecules with a primary in cluster: {len(prim_names)}")

    # identifiability budget
    span_hist = Counter(n_observed(r) for r in reads)
    print(f"\n[STEP 2b] identifiability budget (PSV columns observed per alignment):")
    for k in sorted(span_hist):
        print(f"    spans {k:2d} PSV cols : {span_hist[k]} alignments")
    for K in K_SWEEP:
        nspan = sum(1 for r in reads if n_observed(r) >= K)
        print(f"    alignments spanning >= K={K} PSV cols: {nspan} / {len(reads)}")

    # CO-LOCATED signal: do copies share an intron chain?
    print(f"\n[STEP 2c] co-located signature (shared intron chain across copies?):")
    chain_to_copies = defaultdict(set)
    chain_to_reads = defaultdict(list)
    for r in reads:
        if r["intron_chain"]:
            chain_to_copies[r["intron_chain"]].add(r["assigned_copy"])
            chain_to_reads[r["intron_chain"]].append(r)
    shared = {c: cp for c, cp in chain_to_copies.items() if len(cp) >= 2}
    print(f"    distinct spliced intron chains: {len(chain_to_copies)}")
    print(f"    chains assigned across >=2 copies (SHARED): {len(shared)}")
    # also: a multimapper signal — same read NAME on two copies (pri+sec)
    name_to_copies = defaultdict(set)
    for r in reads:
        name_to_copies[r["name"]].add(r["assigned_copy"])
    multi = {n: c for n, c in name_to_copies.items() if len(c) >= 2}
    print(f"    read MOLECULES placed on >=2 copies (multimappers): {len(multi)}")
    for n, c in sorted(multi.items())[:12]:
        ivs = sorted([f"{r['assigned_copy']}({'sec' if r['is_secondary'] else 'pri'},mq{r['mapq']},{len(r['intron_chain'])}i)"
                      for r in reads if r["name"] == n])
        print(f"      {n[-10:]:10s}: {ivs}")

    print(f"\n[STEP 3+4] copy-split sweep:")
    for K in K_SWEEP:
        for MR in MIN_READS_SWEEP:
            out = split_readchain_by_psv(reads, n_cols, K, MR)
            n_ident = sum(1 for d in out if d["identifiable"])
            n_split_chains = len(set(d["intron_chain"] for d in out if d["identifiable"]))
            print(f"\n  --- K={K} min_reads={MR} ---")
            print(f"  distinct chains={len(set(d['intron_chain'] for d in out))}  emitted units={len(out)}  "
                  f"identifiable(split) copies={n_ident} across {n_split_chains} chain(s)")
            for d in out:
                tag = "SPLIT " if d["identifiable"] else "merged"
                ic = d["intron_chain"]
                chain_desc = f"{len(ic)}i" if ic else "mono"
                loci = Counter(r["assigned_copy"] for r in d["reads"])
                print(f"    [{tag}] {chain_desc:5s} hap={vec_str(d['allele_vector'])[:60]} "
                      f"reads={d['read_count']:3d} copies={dict(loci)}")
    return dict(n_copies=len(copies), n_psv=n_cols, n_align=len(reads),
                idents=idents, span_hist=dict(span_hist),
                shared_chains=len(shared), multimappers=len(multi))


def main():
    summary = {}
    for fam_name, family in FAMILIES.items():
        summary[fam_name] = run_family(fam_name, family)
        print()
    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    for fam, s in summary.items():
        ge2 = sum(v for k, v in s["span_hist"].items() if k >= 2)
        ge3 = sum(v for k, v in s["span_hist"].items() if k >= 3)
        print(f"{fam}: copies={s['n_copies']} PSVs={s['n_psv']} aligns={s['n_align']} "
              f"span>=2PSV={ge2} span>=3PSV={ge3} shared_chains={s['shared_chains']} "
              f"multimappers={s['multimappers']}")


if __name__ == "__main__":
    main()
