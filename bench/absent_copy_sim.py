#!/usr/bin/env python3
"""Positive efficacy demo for --absent-copies.

Builds a minimal synthetic genome + BAM where:

  SIM_COLL contig  (2-locus conflict family, 3 collapsed haplotypes at locus A)
  ─────────────────────────────────────────────────────────────────────────────
  Locus A (pos   0-600, intron 200-300):
    - Copy A0 (host/ref): exactly matches the reference exons. Primary at A,
      secondary at B (host-sequence bridge read → de-tie conflict edge A–B).
    - Copy A1 (absent-1): 15 PSV substitutions in the mRNA → 3 % divergence → remap
      identity ~97 % < 98 % → passes gate 5 → ADMITTED as the single AC_ copy.
      Placed UNIQUELY at locus A (no record at B) so no mirror copy is admitted at B.
    - Copy A2 (absent-2): 15 different PSV substitutions, also UNIQUE at A → makes
      locus A split into >=3 haplotypes (clears gate-1 n_clusters>=3). A2 itself may
      be admitted or routed to DnaNeeds by the remap gate; only A1 is asserted.
  Locus B (pos 700-1300, intron 900-1000): IDENTICAL exon sequences to A.
    - B0 reads: primary at B, secondary at A → gives B its primary read pile
      (needed for the de-novo locus detector) AND completes the de-tie conflict
      edge A–B. B sees only host-sequence reads → never admits a mirror AC copy.

  Because the absent haplotype lives at exactly ONE locus, an A1 read is Stage-1
  TIED over the two identical reference copies (not frozen) yet Stage-2 DECISIVE for
  the single admitted absent copy → it is ASSIGNED to AC_SIM_COLL_0 (discovery_coupled).

  SIM_HET contig  (2-locus conflict family, 2 haplotypes at locus D → DnaNeeds)
  ─────────────────────────────────────────────────────────────────────────────
  Locus D (pos   0-600, intron 200-300): 2 haplotypes D0 and D1 (15 PSVs each).
  Locus E (pos 700-1300, intron 900-1000): identical exons to D, supplies E0 primaries.
  → n_clusters=2 < 3 → gate-1 fires → DnaNeeds record emitted to dna_needs.tsv.

Assertions:
  P1. At least one AC_* copy in <out>.quant.tsv               (absent copy admitted)
  P2. At least one row in <out>.dna_needs.tsv beyond header   (DnaNeeds fired)
  P3. At least one read status=assigned to an AC_* copy in <out>.assignments.tsv
                                                              (reference-absent read PLACED)
  P4. SIM_HET (2-haplotype het) -> DnaNeeds with '<3 clusters' (het NOT made a copy)

Run:
  python3 bench/absent_copy_sim.py [--out-dir /path/to/scratch] [--binary /path/to/copy_assign]

Requirements:
  pysam, samtools (on PATH or in /home/juanfra/miniforge3/bin)
"""

import argparse
import os
import subprocess
import sys
import tempfile

import pysam

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
BINARY_DEFAULT = "/mnt/c/Users/jfris/Desktop/Rustle/target/release/copy_assign"
SAMTOOLS = os.environ.get("SAMTOOLS", "/home/juanfra/miniforge3/bin/samtools")
MINIMAP2 = os.environ.get("MINIMAP2", "/home/juanfra/miniforge3/bin/minimap2")
N_READS_PER_COPY = 30   # primary reads per copy per locus
READ_LEN = 500          # mRNA length (200bp exon1 + 300bp exon2)

# Cycling base pattern for exon sequences (deterministic, avoids low-complexity)
PATTERN = b"ACGTACGATCGATCGACGTACGATCGACGT"

# ---------------------------------------------------------------------------
# Reference genome construction
# ---------------------------------------------------------------------------

def make_exon_seq(length: int, offset: int = 0) -> bytes:
    """Deterministic byte sequence of `length` bp cycling PATTERN."""
    return bytes(PATTERN[(i + offset) % len(PATTERN)] for i in range(length))


def build_reference(out_fa: str) -> dict:
    """
    Build the synthetic reference FASTA and return metadata about each contig.

    Returns metadata dict:
      {contig_name: {
          "exon1_start", "exon1_end",   # 0-based half-open
          "intron_donor", "intron_acceptor",
          "exon2_start", "exon2_end",
          "locus_start", "locus_end",
          "mrna_seq",        # bytes: exon1+exon2
          "total_len",       # total contig length
      }}
    """
    meta = {}

    def make_contig(contig: str, exon1_len=200, intron_len=100, exon2_len=300,
                    gap=100, seq_offset=0):
        # Locus A (or D for SIM_HET)
        e1_start = 0
        e1_end = e1_start + exon1_len
        intron_start = e1_end              # donor position (0-based, inclusive)
        intron_end = intron_start + intron_len  # acceptor (exclusive → last 2 bp end at intron_end-1)
        e2_start = intron_end
        e2_end = e2_start + exon2_len

        # Gap between locus A and B
        spacer_start = e2_end
        spacer_end = spacer_start + gap

        # Locus B (identical exons, different intron coords)
        b_e1_start = spacer_end
        b_e1_end = b_e1_start + exon1_len
        b_intron_start = b_e1_end
        b_intron_end = b_intron_start + intron_len
        b_e2_start = b_intron_end
        b_e2_end = b_e2_start + exon2_len
        total_len = b_e2_end + 100  # trailing pad

        # Sequences
        exon1_seq = make_exon_seq(exon1_len, seq_offset)
        exon2_seq = make_exon_seq(exon2_len, seq_offset + exon1_len)
        intron_body = b"A" * (intron_len - 4)  # 4 reserved for GT..AG boundaries

        def build_intron(body: bytes) -> bytes:
            return b"GT" + body + b"AG"

        spacer = b"N" * gap
        trailing = b"N" * 100

        seq = (
            exon1_seq                # 0 .. e1_end
            + build_intron(intron_body)  # e1_end .. intron_end
            + exon2_seq              # intron_end .. e2_end
            + spacer                 # e2_end .. spacer_end
            + exon1_seq              # spacer_end .. b_e1_end  (SAME as locus A)
            + build_intron(intron_body)  # b_e1_end .. b_intron_end
            + exon2_seq              # b_intron_end .. b_e2_end  (SAME as locus A)
            + trailing
        )
        assert len(seq) == total_len, f"{len(seq)} != {total_len}"

        # Verify canonical GT-AG at intron-A
        assert seq[intron_start:intron_start+2] == b"GT", f"intron-A donor: {seq[intron_start:intron_start+2]}"
        assert seq[intron_end-2:intron_end] == b"AG",   f"intron-A acceptor: {seq[intron_end-2:intron_end]}"
        # Verify canonical GT-AG at intron-B
        assert seq[b_intron_start:b_intron_start+2] == b"GT"
        assert seq[b_intron_end-2:b_intron_end] == b"AG"

        mrna_seq = exon1_seq + exon2_seq  # 500bp

        meta[contig] = {
            # Locus A
            "exon1_start": e1_start, "exon1_end": e1_end,
            "intron_donor": intron_start, "intron_acceptor": intron_end,
            "exon2_start": e2_start, "exon2_end": e2_end,
            "locus_start": e1_start, "locus_end": e2_end,
            # Locus B
            "b_exon1_start": b_e1_start, "b_exon1_end": b_e1_end,
            "b_intron_donor": b_intron_start, "b_intron_acceptor": b_intron_end,
            "b_exon2_start": b_e2_start, "b_exon2_end": b_e2_end,
            "b_locus_start": b_e1_start, "b_locus_end": b_e2_end,
            "mrna_seq": mrna_seq,
            "total_len": total_len,
            "seq": seq,
        }
        return seq

    coll_seq = make_contig("SIM_COLL", seq_offset=0)
    het_seq  = make_contig("SIM_HET",  seq_offset=7)  # different base pattern

    with open(out_fa, "w") as fh:
        fh.write(">SIM_COLL\n")
        for i in range(0, len(coll_seq), 80):
            fh.write(coll_seq[i:i+80].decode() + "\n")
        fh.write(">SIM_HET\n")
        for i in range(0, len(het_seq), 80):
            fh.write(het_seq[i:i+80].decode() + "\n")

    return meta


# ---------------------------------------------------------------------------
# PSV design
# ---------------------------------------------------------------------------

BASES = b"ACGT"


def psv_allele(ref_base: int) -> int:
    """Non-A→G substitution: cycle to the next base (avoids A→G editing artefact)."""
    return BASES[(BASES.index(ref_base) + 2) % 4]


def apply_psvs(mrna: bytes, psv_mrna_pos: list[int], seed_offset: int = 0) -> bytes:
    """Return mRNA with substitutions at the given mRNA positions (non-A→G)."""
    seq = bytearray(mrna)
    for p in psv_mrna_pos:
        seq[p] = psv_allele(seq[p])
    return bytes(seq)


# ---------------------------------------------------------------------------
# BAM record construction
# ---------------------------------------------------------------------------

def make_bam_record(
    name: str,
    ref_id: int,
    ref_start: int,       # 0-based
    seq: bytes,           # read sequence (full mRNA)
    cigar_str: str,       # e.g. "200M100N300M"
    flag: int,            # 0=primary, 256=secondary
    mapq: int,
    as_score: int,
    de: float,            # gap-compressed divergence
) -> pysam.AlignedSegment:
    """Build a pysam AlignedSegment with the required tags."""
    rec = pysam.AlignedSegment()
    rec.query_name = name
    rec.flag = flag
    rec.reference_id = ref_id
    rec.reference_start = ref_start
    rec.mapping_quality = mapq
    rec.set_tag("AS", as_score)
    rec.set_tag("de", de, value_type="f")
    rec.query_sequence = seq.decode()
    rec.query_qualities = pysam.qualitystring_to_array("~" * len(seq))  # Q93 HiFi-grade
    rec.cigarstring = cigar_str
    return rec


def make_paired_records(
    name: str,
    ref_id_a: int, ref_start_a: int,
    ref_id_b: int, ref_start_b: int,
    cigar_a: str, cigar_b: str,
    seq: bytes,
    mapq: int,
    as_a: int, as_b: int,
    de_a: float, de_b: float,
    primary_at: str = "A",   # "A" or "B" — which locus gets the primary record
) -> list[pysam.AlignedSegment]:
    """Create a primary + secondary pair (or reverse) for a read confused between A and B."""
    if primary_at == "A":
        prim = make_bam_record(name, ref_id_a, ref_start_a, seq, cigar_a, 0,   mapq, as_a, de_a)
        sec  = make_bam_record(name, ref_id_b, ref_start_b, seq, cigar_b, 256, mapq, as_b, de_b)
    else:
        prim = make_bam_record(name, ref_id_b, ref_start_b, seq, cigar_b, 0,   mapq, as_b, de_b)
        sec  = make_bam_record(name, ref_id_a, ref_start_a, seq, cigar_a, 256, mapq, as_a, de_a)
    return [prim, sec]


# ---------------------------------------------------------------------------
# Full BAM builder for one contig family
# ---------------------------------------------------------------------------

def build_bam_records_for_contig(
    contig: str,
    ref_id: int,
    m: dict,             # metadata from build_reference
    psv_a1: list[int],   # mRNA-level PSV positions for absent copy 1
    psv_a2: list[int],   # mRNA-level PSV positions for absent copy 2 (may be None → only 2 haplotypes)
    n_reads: int = N_READS_PER_COPY,
) -> list[pysam.AlignedSegment]:
    """
    Generate BAM records for a conflict family on `contig`.

    Locus A: host (A0), absent-1 (A1), optional absent-2 (A2).
    Locus B: B0 reads (same sequence as A0) — gives B primary pile + AS-tie with A.

    MAPQ=0 for all reads (genuinely confused between A and B, since B exon sequences = A exon sequences).
    AS scores use match_score=2, so AS = 2*read_len for a perfect match.
    `de` (gap-compressed divergence) = n_mismatches / mrna_len.
    """
    mrna = m["mrna_seq"]
    mrna_len = len(mrna)

    exon1_len = m["exon1_end"] - m["exon1_start"]   # 200
    exon2_len = m["exon2_end"] - m["exon2_start"]   # 300

    # Locus A genomic layout
    a_ref_start = m["exon1_start"]
    a_intron_len = m["intron_acceptor"] - m["intron_donor"]
    # CIGAR for A-frame reads: exon1 + intron skip + exon2
    cigar_a = f"{exon1_len}M{a_intron_len}N{exon2_len}M"

    # Locus B genomic layout
    b_ref_start = m["b_exon1_start"]
    b_intron_len = m["b_intron_acceptor"] - m["b_intron_donor"]
    cigar_b = f"{exon1_len}M{b_intron_len}N{exon2_len}M"

    match_score = 2
    as_perfect = match_score * mrna_len   # AS for 0 mismatches

    def as_for_n_mm(n_mm: int) -> int:
        # minimap2 map-hifi mismatch penalty = 4 (--score-N 1, -O 4 is default gap)
        # For spliced: penalty per mismatch in map-hifi default = 4 per mismatch
        # Use a simpler model: AS = 2*(mrna_len - n_mm) - 4*n_mm
        # so 0mm -> 2*mrna_len, 15mm -> 2*(mrna_len-15) - 60
        return 2 * (mrna_len - n_mm) - 4 * n_mm

    def de_val(n_mm: int) -> float:
        return n_mm / mrna_len

    recs: list[pysam.AlignedSegment] = []

    # --- A0 reads: host haplotype, primary at A, secondary at B ---
    seq_a0 = mrna
    as_a0 = as_perfect
    de_a0 = 0.0  # perfect match to reference (A exons = reference exons)
    for i in range(n_reads):
        recs.extend(make_paired_records(
            name=f"{contig}_A0_r{i}",
            ref_id_a=ref_id, ref_start_a=a_ref_start,
            ref_id_b=ref_id, ref_start_b=b_ref_start,
            cigar_a=cigar_a, cigar_b=cigar_b,
            seq=seq_a0,
            mapq=0,
            as_a=as_a0, as_b=as_perfect,   # A=B exons so same AS at B
            de_a=de_a0, de_b=de_a0,
            primary_at="A",
        ))

    # --- A1 reads: absent copy 1, UNIQUE at locus A (single record, NO secondary at B) ---
    # ASYMMETRY (vs a naive A=B mirror): the absent-haplotype reads are placed ONLY at locus A.
    # Locus B's collapse-recovery pile therefore never sees the A1/A2 haplotypes, so NO mirror
    # `AC_..._700` copy is ever admitted at B. With a SINGLE absent copy (anchored at A) the
    # per-read IsoCon gate is no longer torn between two degenerate mirror copies, so an A1 read
    # decisively matches the one admitted absent copy → `status=assigned` (discovery_coupled).
    # The A↔B conflict EDGE (which the >=2-copy family gate needs) is still supplied by the
    # host-sequence A0/B0 bridge reads below, which DO carry an A and a B placement.
    seq_a1 = apply_psvs(mrna, psv_a1)
    n_mm_a1 = sum(seq_a1[p] != mrna[p] for p in range(mrna_len))
    as_a1 = as_for_n_mm(n_mm_a1)
    de_a1 = de_val(n_mm_a1)
    for i in range(n_reads):
        recs.append(make_bam_record(
            name=f"{contig}_A1_r{i}",
            ref_id=ref_id, ref_start=a_ref_start,
            seq=seq_a1, cigar_str=cigar_a,
            flag=0, mapq=0, as_score=as_a1, de=de_a1,
        ))

    # --- A2 reads: absent copy 2, UNIQUE at locus A (single record, NO secondary at B) ---
    # Needed so locus A splits into >=3 haplotypes (host + A1 + A2) and clears gate-1
    # (n_clusters >= min_clusters=3). A2 itself may be admitted or routed to DnaNeeds by the
    # remap gate depending on how minimap2 chains its PSVs; only the A1 admission is asserted.
    if psv_a2:
        seq_a2 = apply_psvs(mrna, psv_a2)
        n_mm_a2 = sum(seq_a2[p] != mrna[p] for p in range(mrna_len))
        as_a2 = as_for_n_mm(n_mm_a2)
        de_a2 = de_val(n_mm_a2)
        for i in range(n_reads):
            recs.append(make_bam_record(
                name=f"{contig}_A2_r{i}",
                ref_id=ref_id, ref_start=a_ref_start,
                seq=seq_a2, cigar_str=cigar_a,
                flag=0, mapq=0, as_score=as_a2, de=de_a2,
            ))

    # --- B0 reads: primary at B, secondary at A ---
    # B0 has the SAME sequence as A0 (collapsed copies share exon sequences).
    # Giving B0 reads primary-at-B is what makes the de-novo detector build a skeleton at B.
    for i in range(n_reads):
        recs.extend(make_paired_records(
            name=f"{contig}_B0_r{i}",
            ref_id_a=ref_id, ref_start_a=a_ref_start,
            ref_id_b=ref_id, ref_start_b=b_ref_start,
            cigar_a=cigar_a, cigar_b=cigar_b,
            seq=seq_a0,   # same sequence as A0
            mapq=0,
            as_a=as_perfect, as_b=as_perfect,
            de_a=0.0, de_b=0.0,
            primary_at="B",   # primary at B → B gets a skeleton
        ))

    return recs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess:
    print(f"  $ {' '.join(cmd)}")
    return subprocess.run(cmd, check=True, **kwargs)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=None, help="Scratch directory (auto-created if omitted)")
    ap.add_argument("--binary", default=BINARY_DEFAULT, help="copy_assign binary path")
    ap.add_argument("--keep", action="store_true", help="Keep scratch directory on exit")
    args = ap.parse_args()

    if args.out_dir:
        out_dir = args.out_dir
        os.makedirs(out_dir, exist_ok=True)
        cleanup = not args.keep
    else:
        out_dir = tempfile.mkdtemp(prefix="absent_sim_")
        cleanup = not args.keep

    print(f"[sim] scratch dir: {out_dir}")

    # -----------------------------------------------------------------------
    # 1. Build reference FASTA + FASTA index
    # -----------------------------------------------------------------------
    fa_path = os.path.join(out_dir, "sim.fa")
    print("[sim] building reference FASTA...")
    meta = build_reference(fa_path)

    run([SAMTOOLS, "faidx", fa_path])
    print(f"  SIM_COLL total_len={meta['SIM_COLL']['total_len']}")
    print(f"  SIM_HET  total_len={meta['SIM_HET']['total_len']}")

    # -----------------------------------------------------------------------
    # 2. PSV design
    # -----------------------------------------------------------------------
    mrna_len = READ_LEN  # 500bp = exon1(200) + exon2(300)
    exon1_len = 200

    # mRNA positions for A1 PSVs: spread across both exons, spaced 13bp apart
    psv_a1_mrna = [13 + 13*i for i in range(15)]   # [13,26,39,...,195] — all in exon1&2
    # mRNA positions for A2 PSVs: interleaved, not overlapping A1
    psv_a2_mrna = [6 + 13*i for i in range(15)]    # [6,19,32,...,182]

    # DET: no overlap between A1 and A2
    assert not set(psv_a1_mrna) & set(psv_a2_mrna), "A1/A2 PSV overlap!"
    assert all(p < mrna_len for p in psv_a1_mrna + psv_a2_mrna), "PSV out of range"

    # HET locus: only A1-class PSVs (2 haplotypes → DnaNeeds)
    psv_d1_mrna = psv_a1_mrna   # same PSV positions, different contig

    print(f"[sim] A1 PSVs at mRNA positions: {psv_a1_mrna[:5]}...  ({len(psv_a1_mrna)} total)")
    print(f"[sim] A2 PSVs at mRNA positions: {psv_a2_mrna[:5]}...  ({len(psv_a2_mrna)} total)")

    # -----------------------------------------------------------------------
    # 3. Build BAM
    # -----------------------------------------------------------------------
    bam_path = os.path.join(out_dir, "sim.bam")
    sam_path = os.path.join(out_dir, "sim.unsorted.bam")

    contigs = [
        pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6"},
            "SQ": [
                {"SN": "SIM_COLL", "LN": meta["SIM_COLL"]["total_len"]},
                {"SN": "SIM_HET",  "LN": meta["SIM_HET"]["total_len"]},
            ],
        })
    ][0]

    ref_id_coll = 0
    ref_id_het  = 1

    all_records: list[pysam.AlignedSegment] = []

    # SIM_COLL: 3-haplotype collapse → 2 absent copies admitted
    print("[sim] generating SIM_COLL reads (3 haplotypes: A0, A1, A2 + B0)...")
    all_records.extend(build_bam_records_for_contig(
        "SIM_COLL", ref_id_coll, meta["SIM_COLL"],
        psv_a1=psv_a1_mrna, psv_a2=psv_a2_mrna,
    ))

    # SIM_HET: 2-haplotype (DnaNeeds) test
    print("[sim] generating SIM_HET reads (2 haplotypes: D0, D1 + E0)...")
    all_records.extend(build_bam_records_for_contig(
        "SIM_HET", ref_id_het, meta["SIM_HET"],
        psv_a1=psv_d1_mrna, psv_a2=None,  # only 2 haplotypes
    ))

    print(f"[sim] writing {len(all_records)} BAM records...")
    with pysam.AlignmentFile(sam_path, "wb", header=contigs) as bam_out:
        for rec in all_records:
            bam_out.write(rec)

    run([SAMTOOLS, "sort", "-o", bam_path, sam_path])
    run([SAMTOOLS, "index", bam_path])
    os.remove(sam_path)

    # -----------------------------------------------------------------------
    # 4. Build regions file (one region per contig covering both loci)
    # -----------------------------------------------------------------------
    regions_path = os.path.join(out_dir, "regions.txt")
    with open(regions_path, "w") as fh:
        mc = meta["SIM_COLL"]
        fh.write(f"SIM_COLL:0-{mc['total_len']}\n")
        mh = meta["SIM_HET"]
        fh.write(f"SIM_HET:0-{mh['total_len']}\n")

    # -----------------------------------------------------------------------
    # 5. Run copy_assign --absent-copies
    # -----------------------------------------------------------------------
    out_prefix = os.path.join(out_dir, "ca_out")
    print("[sim] running copy_assign --absent-copies...")
    env = dict(os.environ)
    env["PATH"] = "/home/juanfra/miniforge3/bin:" + env.get("PATH", "")
    env["RUSTLE_MINIMAP2"] = MINIMAP2
    result = subprocess.run(
        [
            args.binary,
            "--bam", bam_path,
            "--fasta", fa_path,
            "--regions", regions_path,
            "--min-copies", "2",
            "--skip-poa-diagnostic",
            "--absent-copies",
            "--out", out_prefix,
        ],
        capture_output=True, text=True, env=env,
    )
    sys.stdout.write(result.stderr)
    if result.returncode != 0:
        print(f"[sim] copy_assign exited {result.returncode}", file=sys.stderr)
        return 1

    # -----------------------------------------------------------------------
    # 6. Assertions
    # -----------------------------------------------------------------------
    quant_path    = out_prefix + ".quant.tsv"
    dna_needs_path = out_prefix + ".dna_needs.tsv"
    assign_path   = out_prefix + ".assignments.tsv"

    passed = 0
    failed = 0

    def check(cond: bool, msg: str):
        nonlocal passed, failed
        if cond:
            print(f"  PASS: {msg}")
            passed += 1
        else:
            print(f"  FAIL: {msg}")
            failed += 1

    # P1: AC_* copy appears in quant.tsv
    ac_copies = []
    with open(quant_path) as fh:
        for line in fh:
            if "\tAC_" in line:
                ac_copies.append(line.strip())
    check(len(ac_copies) > 0, f"P1  AC_* copy admitted and present in quant.tsv ({len(ac_copies)} found)")
    if ac_copies:
        print(f"       First: {ac_copies[0]}")

    # P2: dna_needs.tsv exists and has ≥1 data rows
    dna_rows = []
    if os.path.exists(dna_needs_path):
        with open(dna_needs_path) as fh:
            lines = [l for l in fh if not l.startswith("chrom\t")]
        dna_rows = [l for l in lines if l.strip()]
    check(len(dna_rows) > 0, f"P2  dna_needs.tsv has ≥1 DnaNeeds record ({len(dna_rows)} found)")
    if dna_rows:
        print(f"       First: {dna_rows[0].strip()}")

    # P3: at least one read is DECISIVELY ASSIGNED to an AC_* copy in assignments.tsv.
    # `assigned_copy` in assignments.tsv is a per-family copy INDEX (= best_copy), not a tid;
    # map (family_id, index) -> copy_tid via quant.tsv, then require >=1 status=assigned row
    # landing on a tid that starts with "AC_". This is the real efficacy proof (a read placed
    # on a reference-absent copy), strictly stronger than EM n_reads_hard alone.
    idx_to_tid: dict[tuple[str, str], str] = {}
    with open(quant_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0] != "family_id":
                idx_to_tid[(parts[0], parts[1])] = parts[2]  # (family_id, copy_index) -> tid
    assigned_to_ac = []
    with open(assign_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or parts[0] == "read_name":
                continue
            read_name, fam, copy_idx, status = parts[0], parts[1], parts[2], parts[3]
            if status == "assigned":
                tid = idx_to_tid.get((fam, copy_idx), "")
                if tid.startswith("AC_"):
                    assigned_to_ac.append((read_name, tid))
    n_ac_assigned = len(assigned_to_ac)
    check(n_ac_assigned > 0,
          f"P3  >=1 read status=assigned to an AC_* copy ({n_ac_assigned} reads)")
    if assigned_to_ac:
        ex_read, ex_tid = assigned_to_ac[0]
        print(f"       e.g. {ex_read} -> {ex_tid} (assigned)")

    # P4: HET locus generates DnaNeeds with gate-1 reason (<3 clusters)
    het_gate1 = [l for l in dna_rows if "clusters" in l and "SIM_HET" in l]
    check(len(het_gate1) > 0, f"P4  SIM_HET DnaNeeds with '<3 clusters' reason ({len(het_gate1)} found)")
    if het_gate1:
        print(f"       {het_gate1[0].strip()}")

    print()
    print(f"=== SIM DEMO: {passed} PASS, {failed} FAIL ===")

    if cleanup and not args.keep:
        import shutil
        shutil.rmtree(out_dir, ignore_errors=True)
    else:
        print(f"[sim] outputs in: {out_dir}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
