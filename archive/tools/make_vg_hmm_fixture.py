#!/usr/bin/env python3
"""Generate synthetic VG-HMM test fixtures (deterministic, no external simulators).

Outputs:
  test_data/vg_hmm/synthetic_family.fa   — two "gene copies" (reference + perturbed paralog)
  test_data/vg_hmm/synthetic_reads.bam   — 50 unmapped reads spanning the novel paralog
  test_data/vg_hmm/synthetic_reads.bam.bai (not produced — unmapped BAMs don't need an index)

The synthetic paralog is built by:
  1. Taking a deterministic 3000-bp "exonic" template (ACGT repeat pattern).
  2. Inserting a 50-bp copy-specific exon at position 1000.
  3. Perturbing 5% of bases at random (seed=42) — PacBio-HiFi-like divergence.

Reads are 1000-bp substrings of the perturbed sequence with 1% per-base noise (seed=43).
All reads are written as **unmapped** BAM records (FLAG=4) with no reference coordinates.
The BAM SQ header lists the reference copy as the "genome" so rustle's reader is happy.

Usage:
  python3 tools/make_vg_hmm_fixture.py [--out-dir test_data/vg_hmm]
"""

import argparse
import random
import struct
import gzip
import io
import sys
import os
import subprocess
from pathlib import Path

# ── helpers ───────────────────────────────────────────────────────────────────

BASES = b"ACGT"
COMP  = {ord('A'): 'T', ord('C'): 'G', ord('G'): 'C', ord('T'): 'A', ord('N'): 'N'}


def make_reference_exons(seed: int = 0) -> bytes:
    """3000-bp deterministic reference sequence (two 1500-bp exons)."""
    rng = random.Random(seed)
    bases = [rng.choice(b"ACGT") for _ in range(3000)]
    return bytes(bases)


def perturb(seq: bytes, frac: float, seed: int) -> bytes:
    """Randomly substitute `frac` fraction of bases (in-place copy)."""
    rng = random.Random(seed)
    arr = bytearray(seq)
    n_mut = int(len(seq) * frac)
    positions = rng.sample(range(len(seq)), k=n_mut)
    for pos in positions:
        current = arr[pos]
        choices = [b for b in b"ACGT" if b != current]
        arr[pos] = rng.choice(choices)
    return bytes(arr)


def add_noise(seq: bytes, frac: float, seed: int) -> bytes:
    """Add per-base substitution noise (simulates HiFi ~1% error rate)."""
    rng = random.Random(seed)
    arr = bytearray(seq)
    for i in range(len(arr)):
        if rng.random() < frac:
            choices = [b for b in b"ACGT" if b != arr[i]]
            arr[i] = rng.choice(choices)
    return bytes(arr)


# ── BAM writer (minimal, no external deps beyond stdlib) ─────────────────────
# We write a valid BGZF/BAM file using only stdlib + bgzf compression.
# Each BGZF block is independently gzip-compressed with a special extra field.

def _bgzf_block(data: bytes) -> bytes:
    """Wrap `data` in a single BGZF block."""
    # BGZF extra field: SI1=66 SI2=67 SLEN=2 BSIZE=<total-block-size - 1>
    # We'll compute BSIZE after compressing.
    compressed = zlib_compress(data)
    # BGZF header: ID1 ID2 CM FLG MTIME XFL OS XLEN=6 SI1 SI2 SLEN=2 BSIZE(placeholder)
    xlen = 6
    header_before_bsize = bytes([
        0x1f, 0x8b,             # ID1, ID2
        8,                      # CM = deflate
        4,                      # FLG = FEXTRA
        0, 0, 0, 0,             # MTIME
        0,                      # XFL
        0xff,                   # OS = unknown
        xlen, 0,                # XLEN
        66, 67,                 # SI1, SI2
        2, 0,                   # SLEN = 2
    ])  # BSIZE will follow
    # total block size = len(header_before_bsize) + 2 (for BSIZE) + len(compressed) + 8 (CRC+ISIZE)
    total = len(header_before_bsize) + 2 + len(compressed) + 8
    bsize = total - 1  # BSIZE field = block_size - 1
    block = header_before_bsize + struct.pack("<H", bsize) + compressed

    # Append CRC32 and ISIZE
    import zlib as _zlib
    crc = _zlib.crc32(data) & 0xFFFFFFFF
    isize = len(data) & 0xFFFFFFFF
    block += struct.pack("<II", crc, isize)
    return block


def zlib_compress(data: bytes) -> bytes:
    """Deflate-compress data (raw, no zlib header)."""
    import zlib as _zlib
    obj = _zlib.compressobj(_zlib.Z_DEFAULT_COMPRESSION, _zlib.DEFLATED, -15)
    return obj.compress(data) + obj.flush()


def bgzf_eof() -> bytes:
    """Standard BGZF EOF block (28 bytes)."""
    return bytes([
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
        0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
        0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00,
    ])


def encode_seq_qual(seq: bytes) -> tuple[bytes, bytes]:
    """Encode sequence in 4-bit BAM format and return (seq_bytes, qual_bytes)."""
    nib_map = {ord('A'): 1, ord('C'): 2, ord('G'): 4, ord('T'): 8, ord('N'): 15}
    n = len(seq)
    out = bytearray((n + 1) // 2)
    for i, b in enumerate(seq):
        nib = nib_map.get(b, 15)
        if i % 2 == 0:
            out[i // 2] = nib << 4
        else:
            out[i // 2] |= nib
    qual = bytes([0xff] * n)  # 0xff = no quality stored
    return bytes(out), qual


def bam_record(qname: bytes, flag: int, seq: bytes) -> bytes:
    """Build a minimal BAM record for an unmapped read."""
    seq_encoded, qual = encode_seq_qual(seq)
    l_read_name = len(qname) + 1  # null-terminated
    n_cigar_op = 0
    l_seq = len(seq)

    # BAM core is 9 int32 fields (36 bytes) — standard layout
    # ref: https://samtools.github.io/hts-specs/SAMv1.pdf §4.2
    refID    = struct.pack("<i", -1)
    pos      = struct.pack("<i", -1)  # 0-based; -1 for unmapped
    l_read_name_b = struct.pack("<B", l_read_name)
    mapq     = struct.pack("<B", 255)
    bin_val  = struct.pack("<H", 4680)  # unmapped bin
    n_cigar  = struct.pack("<H", 0)
    flag_b   = struct.pack("<H", flag)
    l_seq_b  = struct.pack("<i", l_seq)
    next_ref = struct.pack("<i", -1)
    next_pos = struct.pack("<i", -1)
    tlen     = struct.pack("<i", 0)

    core = (refID + pos + l_read_name_b + mapq + bin_val +
            n_cigar + flag_b + l_seq_b + next_ref + next_pos + tlen)
    data = qname + b'\x00' + seq_encoded + qual

    rec_size = len(core) + len(data)
    return struct.pack("<i", rec_size) + core + data


def write_bam(path: Path, sq_name: str, sq_len: int, records: list[bytes]) -> None:
    """Write a minimal valid BAM file with one SQ header line."""
    # BAM text header
    text_header = f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{sq_name}\tLN:{sq_len}\n"
    text_bytes = text_header.encode()

    # BAM binary header: "BAM\1" + int32 l_text + text + int32 n_ref + (ref blocks)
    bam_header = b"BAM\x01"
    bam_header += struct.pack("<i", len(text_bytes))
    bam_header += text_bytes
    bam_header += struct.pack("<i", 1)  # n_ref = 1
    sq_name_b = sq_name.encode() + b'\x00'
    bam_header += struct.pack("<i", len(sq_name_b))
    bam_header += sq_name_b
    bam_header += struct.pack("<i", sq_len)

    # Concatenate all record bytes
    all_records = b"".join(records)

    # Write as BGZF: one block for header, blocks for records
    payload = bam_header + all_records
    # Split into 64 KB blocks
    MAX_BLOCK = 65536
    out = b""
    offset = 0
    while offset < len(payload):
        chunk = payload[offset : offset + MAX_BLOCK]
        out += _bgzf_block(chunk)
        offset += len(chunk)
    out += bgzf_eof()

    path.write_bytes(out)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="test_data/vg_hmm")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Build reference copy (copy 0) — 3000 bp.
    ref_seq = make_reference_exons(seed=0)
    ref_name = "copy0"

    # 2. Build novel paralog (copy 1):
    #    - Insert 50-bp copy-specific exon at position 1000.
    COPY_SPECIFIC_EXON = (b"ACGTACGTACGTACGTACGT" * 3)[:50]  # deterministic 50-bp insert
    novel_base = ref_seq[:1000] + COPY_SPECIFIC_EXON + ref_seq[1000:]
    # Perturb 5% of bases.
    novel_seq = perturb(novel_base, frac=0.05, seed=42)
    novel_name = "copy1_novel"

    # 3. Write FASTA.
    fa_path = out_dir / "synthetic_family.fa"
    with fa_path.open("w") as fh:
        # Write 80-char-wide FASTA.
        fh.write(f">{ref_name}\n")
        for i in range(0, len(ref_seq), 80):
            fh.write(ref_seq[i:i+80].decode() + "\n")
        fh.write(f">{novel_name}\n")
        for i in range(0, len(novel_seq), 80):
            fh.write(novel_seq[i:i+80].decode() + "\n")

    print(f"Wrote {fa_path} ({fa_path.stat().st_size} bytes)")

    # 4. Simulate 50 reads ~1000 bp each from the novel paralog sequence.
    READ_LEN = 1000
    N_READS = 50
    rng = random.Random(43)
    read_records = []
    for i in range(N_READS):
        start = rng.randint(0, len(novel_seq) - READ_LEN)
        raw_read = novel_seq[start : start + READ_LEN]
        noisy_read = add_noise(raw_read, frac=0.01, seed=43 + i)
        qname = f"synth_read_{i:04d}".encode()
        flag = 4  # unmapped
        read_records.append(bam_record(qname, flag, noisy_read))

    bam_path = out_dir / "synthetic_reads.bam"
    write_bam(bam_path, sq_name=ref_name, sq_len=len(ref_seq), records=read_records)
    print(f"Wrote {bam_path} ({bam_path.stat().st_size} bytes, {N_READS} reads)")
    print("Done.")


if __name__ == "__main__":
    main()
