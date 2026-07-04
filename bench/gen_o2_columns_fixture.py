#!/usr/bin/env python3
"""Golden byte-parity fixture generator for the O2 VG-materialization step-2 Rust port
(src/rustle/vg_family/o2_columns.rs, porting bench/psv_graph_genomewide.py::column_alleles,
the pysam-pileup BAM parser).

Run deterministically:
    PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_o2_columns_fixture.py

Ships small BAMs (+ .bai) into src/rustle/vg_family/testdata/ and emits
o2_columns_fixture.json = { bams: [{file, ref, kind, golden}], excluded_read_names }.
`golden` is the EXACT pysam `column_alleles` output: {ref_pos(str) -> {read_name -> base}}.
The Rust test opens each shipped BAM, runs its own column_alleles, and asserts map-equality.

Two fixture sources (matching the FULL O2 port pattern):
  (A) REAL BAMs mined from the persisted o2vg family scratch dirs -- genuine minimap2
      asm20 copy alignments (copies.bam, with real D/I/large-I CIGARs) and splice:hq read
      alignments (reads.bam, WITH N introns, soft-clips and reverse-strand reads).
  (B) One ADVERSARIAL synthetic BAM built from hand-written SAM exercising every CIGAR op
      (M/=/X, D, N-intron, I, S, H, P) and every excluded flag (SECONDARY/SUPPLEMENTARY/
      DUP/QCFAIL/UNMAP), plus reverse strand (no revcomp), multi-read columns, full-vs-
      partial overlap, an off-contig read, and an in-header contig with no reads.

Internal safety: for EVERY shipped BAM the REAL pysam column_alleles is cross-checked
cell-for-cell against a plain per-read CIGAR walk (the exact algorithm the Rust port
implements). If they ever diverge the generator aborts -- so a passing generator proves
the Rust port (which is that walk) matches pysam.
"""
import collections
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pysam  # noqa: E402
import psv_graph_genomewide as pg  # noqa: E402  the REAL column_alleles under test

TESTDATA = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata")
OUT = os.path.join(TESTDATA, "o2_columns_fixture.json")
O2VG = "/home/juanfra/winloci_scratch/o2vg"
SAMTOOLS = "/home/juanfra/miniforge3/bin/samtools"


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


# ---- plain per-read CIGAR walk == the exact algorithm the Rust port implements ----
# BAM op codes: M=0 I=1 D=2 N=3 S=4 H=5 P=6 EQ=7 X=8
def plain_walk(bampath, ref):
    out = {}
    bam = pysam.AlignmentFile(bampath, "rb")
    for r in bam.fetch(until_eof=True):
        if r.flag & 0xF04:                      # UNMAP|SECONDARY|QCFAIL|DUP|SUPPLEMENTARY
            continue
        if r.reference_name != ref:
            continue
        seq = r.query_sequence
        if seq is None:
            continue
        refpos = r.reference_start
        qpos = 0
        for (op, length) in (r.cigartuples or []):
            if op in (0, 7, 8):                 # M, =, X  -> aligned base per ref pos
                for k in range(length):
                    out.setdefault(refpos + k, {})[r.query_name] = seq[qpos + k]
                refpos += length
                qpos += length
            elif op in (2, 3):                  # D, N    -> reference-only
                refpos += length
            elif op in (1, 4):                  # I, S    -> query-only
                qpos += length
            elif op in (5, 6):                  # H, P    -> neither
                pass
            else:
                raise SystemExit(f"unknown CIGAR op {op}")
    bam.close()
    return out


def golden_for(bampath, ref):
    """Run the REAL pysam column_alleles, cross-check the plain walk, return JSON golden."""
    gold = {p: dict(d) for p, d in pg.column_alleles(bampath, ref).items()}
    walk = {p: dict(d) for p, d in plain_walk(bampath, ref).items()}
    if gold != walk:
        # localize
        for p in sorted(set(gold) | set(walk)):
            if gold.get(p) != walk.get(p):
                raise SystemExit(f"[FATAL] walk != pysam in {bampath} ref={ref} at pos {p}: "
                                 f"pysam={gold.get(p)} walk={walk.get(p)}")
    ncells = sum(len(d) for d in gold.values())
    print(f"    {os.path.basename(bampath)} ref={ref}: {len(gold)} cols, {ncells} cells  (pysam==walk OK)")
    return {str(p): dict(sorted(d.items())) for p, d in sorted(gold.items())}, ncells


# ------------------------------------------------------------------ (A) real BAMs
REAL = [
    # (src family dir, src filename, shipped filename, ref, note)
    ("tmp_fam93", "copies.bam", "o2cols_fam93_copies.bam", "copy0",
     "minimap2 asm20 copy alignment (revcomp copy1 onto copy0)"),
    ("tmp_fam93", "reads.bam", "o2cols_fam93_reads.bam", "copy0",
     "splice:hq reads with N introns, insertions, soft-clips, reverse strand"),
    ("tmp_fam77", "copies.bam", "o2cols_fam77_copies.bam", "copy0",
     "asm20 copy alignment with rich D / large 120I / multi-D CIGAR"),
]


def ship_real():
    entries = []
    for famdir, srcname, dst, ref, note in REAL:
        src = os.path.join(O2VG, famdir, srcname)
        if not os.path.exists(src):
            raise SystemExit(f"missing real BAM {src}")
        dstpath = os.path.join(TESTDATA, dst)
        shutil.copy(src, dstpath)
        # ensure an index exists (pysam column_alleles pileup needs it); ship it too
        srcbai = src + ".bai"
        dstbai = dstpath + ".bai"
        if os.path.exists(srcbai):
            shutil.copy(srcbai, dstbai)
        else:
            sh(f"{SAMTOOLS} index {dstpath}")
        golden, ncells = golden_for(dstpath, ref)
        entries.append(dict(file=dst, ref=ref, kind="real", note=note, golden=golden))
    return entries


# ------------------------------------------------------------- (B) adversarial BAM
# One BAM, header contigs: adv (60), empty (30), other (30). QUAL provided (len==SEQ)
# so no base-quality subtlety; min_base_quality=0 keeps every base anyway.
ADV_HEADER = (
    "@HD\tVN:1.6\tSO:coordinate\n"
    "@SQ\tSN:adv\tLN:60\n"
    "@SQ\tSN:empty\tLN:30\n"
    "@SQ\tSN:other\tLN:30\n"
)

# (qname, flag, rname, pos(1-based), cigar, seq)  -- QUAL auto = 'I'*len(seq) (or * if unmapped)
ADV_READS = [
    # plain full-span forward match: ref 0..19
    ("r_full",         0,     "adv", 1,  "20M",           "ACGTACGTACGTACGTACGT"),
    # partial overlap (ref 10..19) -> columns 10..19 carry TWO reads
    ("r_partial",      0,     "adv", 11, "10M",           "GGGGGCCCCC"),
    # DELETION: ref 0..4 (M) , ref 5..6 (D, no base) , ref 7..11 (M)
    ("r_del",          0,     "adv", 1,  "5M2D5M",        "AAAAATTTTT"),
    # INTRON (N splice): ref 0..4 (M) , skip 5..14 , ref 15..19 (M)
    ("r_intron",       0,     "adv", 1,  "5M10N5M",       "CCCCCGGGGG"),
    # INSERTION: ref 0..4 (M) , 3 inserted bases (no ref) , ref 5..9 (M)
    ("r_ins",          0,     "adv", 1,  "5M3I5M",        "TTTTTNNNAAAAA"),
    # SOFT-CLIP both ends: 3S skip, ref 0..9 (M), 4S skip
    ("r_softclip",     0,     "adv", 1,  "3S10M4S",       "GGGACGTACGTACTTTT"),
    # HARD-CLIP both ends (clipped bases absent from SEQ): ref 0..9 (M)
    ("r_hardclip",     0,     "adv", 1,  "3H10M2H",       "CACACACACA"),
    # EQ/X ops: ref 0..3 (=), ref 4..5 (X), ref 6..9 (M)
    ("r_eqx",          0,     "adv", 1,  "4=2X4M",        "AAAACCGGGG"),
    # REVERSE strand (SEQ already forward-in-ref; NO revcomp): ref0 must be 'T' (=seq[0])
    ("r_reverse",      16,    "adv", 1,  "10M",           "TGTGTGTGTG"),
    # reverse + intron + soft-clip combo: 2S, ref4..9 (M), skip 10..13, ref14..19 (M), 3S
    ("r_combo",        16,    "adv", 5,  "2S6M4N6M3S",    "AAACGTACGTACGTCCC"),
    # ---- excluded by flag: must NOT contribute a single base ----
    ("r_secondary",    0x100, "adv", 1,  "10M",           "AAAAAAAAAA"),
    ("r_supplementary",0x800, "adv", 1,  "10M",           "CCCCCCCCCC"),
    ("r_dup",          0x400, "adv", 1,  "10M",           "GGGGGGGGGG"),
    ("r_qcfail",       0x200, "adv", 1,  "10M",           "TTTTTTTTTT"),
    ("r_unmap",        0x4,   "adv", 1,  "*",             "AAAAAAAAAA"),
    # off-contig read: on "other", must NOT appear when querying "adv"
    ("r_other",        0,     "other", 1, "10M",          "ACGTACGTAC"),
    # (contig "empty" intentionally has NO reads)
]

# names that must never appear in the adv column map (flag- or contig-excluded)
EXCLUDED_NAMES = ["r_secondary", "r_supplementary", "r_dup", "r_qcfail", "r_unmap", "r_other"]


def build_adv():
    sampath = os.path.join(TESTDATA, "o2cols_adv.sam")
    with open(sampath, "w") as fh:
        fh.write(ADV_HEADER)
        for (qname, flag, rname, pos, cigar, seq) in ADV_READS:
            unmapped = bool(flag & 0x4)
            qual = "*" if (unmapped or cigar == "*") else ("I" * len(seq))
            mapq = 0 if unmapped else 60
            fh.write(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n")
    bam = os.path.join(TESTDATA, "o2cols_adv.bam")
    sh(f"{SAMTOOLS} view -bS {sampath} | {SAMTOOLS} sort -o {bam} -")
    sh(f"{SAMTOOLS} index {bam}")
    os.remove(sampath)

    entries = []
    golden_adv, _ = golden_for(bam, "adv")
    entries.append(dict(file="o2cols_adv.bam", ref="adv", kind="adversarial",
                        note="all CIGAR ops + all excluded flags + reverse/no-revcomp + multi-read cols",
                        golden=golden_adv))
    # in-header contig with no reads -> empty map
    golden_empty, _ = golden_for(bam, "empty")
    assert golden_empty == {}, f"empty contig not empty: {golden_empty}"
    entries.append(dict(file="o2cols_adv.bam", ref="empty", kind="adversarial",
                        note="in-header contig with zero alignments (empty map)", golden=golden_empty))

    # ---- explicit adversarial assertions on the golden itself ----
    g = {int(p): d for p, d in golden_adv.items()}
    # excluded reads leaked nowhere
    for p, d in g.items():
        for nm in d:
            assert nm not in EXCLUDED_NAMES, f"excluded read {nm} leaked at pos {p}"
    # reverse-strand: NO revcomp -> r_reverse base at ref0 is seq[0]='T' (revcomp would be 'C')
    assert g[0]["r_reverse"] == "T", f"reverse read revcomp'd? got {g[0].get('r_reverse')}"
    assert g[9]["r_reverse"] == "G", g[9].get("r_reverse")   # seq[9]='G'
    # deletion: no base at ref5/ref6 from r_del
    assert "r_del" not in g.get(5, {}) and "r_del" not in g.get(6, {}), "r_del leaked into deleted ref"
    assert g[7]["r_del"] == "T", "r_del base after deletion wrong"
    # intron: no base at ref 5..14 from r_intron; base resumes at ref15
    for p in range(5, 15):
        assert "r_intron" not in g.get(p, {}), f"r_intron leaked into intron at {p}"
    assert g[15]["r_intron"] == "G", "r_intron base after intron wrong"
    # insertion: inserted bases skipped -> ref5 from r_ins is 'A' (seq[8]), not 'N'
    assert g[5]["r_ins"] == "A", f"r_ins insertion mis-walked, ref5={g[5].get('r_ins')}"
    # soft-clip: ref0 from r_softclip is seq[3]='A'
    assert g[0]["r_softclip"] == "A", "r_softclip left-clip offset wrong"
    # hard-clip: ref0 from r_hardclip is seq[0]='C' (H not in SEQ)
    assert g[0]["r_hardclip"] == "C", "r_hardclip offset wrong"
    # multi-read column: ref10 carries both r_full and r_partial
    assert "r_full" in g[10] and "r_partial" in g[10], "multi-read column missing a read"
    # EQ/X: ref4 (X) from r_eqx is seq[4]='C'
    assert g[4]["r_eqx"] == "C", "r_eqx X-op base wrong"
    # combo (reverse+intron+softclip): ref4 = seq[2]='A', ref14 = seq[8]='G', intron 10..13 empty
    assert g[4]["r_combo"] == "A" and g[14]["r_combo"] == "G", "r_combo walk wrong"
    for p in range(10, 14):
        assert "r_combo" not in g.get(p, {}), f"r_combo leaked into intron at {p}"
    print("    adversarial golden assertions: all pass")
    return entries


def main():
    os.makedirs(TESTDATA, exist_ok=True)
    print("[A] real BAMs:")
    real_entries = ship_real()
    print("[B] adversarial BAM:")
    adv_entries = build_adv()

    fixture = dict(
        note=("Golden byte-parity fixture for src/rustle/vg_family/o2_columns.rs "
              "(port of bench/psv_graph_genomewide.py::column_alleles). golden = exact pysam "
              "column_alleles output {ref_pos -> {read_name -> base}}. Generated by "
              "bench/gen_o2_columns_fixture.py (PYTHONHASHSEED=0)."),
        excluded_read_names=EXCLUDED_NAMES,
        bams=real_entries + adv_entries,
    )
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, separators=(",", ":"))
    total = sum(sum(len(d) for d in e["golden"].values()) for e in fixture["bams"])
    print(f"\nwrote {OUT}")
    print(f"  bams: {len(fixture['bams'])}  (real={len(real_entries)}, adversarial={len(adv_entries)})")
    print(f"  total column-cells across golden: {total}")
    print(f"  fixture json size: {os.path.getsize(OUT)} bytes")


if __name__ == "__main__":
    main()
