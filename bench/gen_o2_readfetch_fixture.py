#!/usr/bin/env python3
"""Golden-fixture generator for the O2 VG-materialization READ-FETCH port (step 4).

Runs the EXACT read-fetch block of bench/o2_vg_visualization.py::materialize_family
(~lines 139-161) against the real pg.BAM for a handful of real families, and emits
src/rustle/vg_family/testdata/o2_readfetch_fixture.json so the Rust
o2_materialize::fetch_reads can be checked byte-for-byte:
  * regions  (dedup_copies output, sorted by -length -- the EXACT input to fetch_reads)
  * best     qname -> [region_idx, de_as_raw_f64_bits]  (map; order-independent)
  * seqs     insertion-ordered qname list, PRE-cap and POST-cap(=6000)
  * lowered-cap survivor lists to exercise the READ_CAP truncation path (no family
    hits 6000 naturally on this substrate).

The seq STRINGS are intentionally omitted from the golden (long); parity is asserted
on the qname ORDER + which survive the cap, exactly as the task specifies.

Run:  PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_o2_readfetch_fixture.py
"""
import collections
import json
import os
import struct
import sys

os.environ.setdefault("PYTHONHASHSEED", "0")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pysam                       # noqa: E402
import psv_graph_genomewide as pg  # noqa: E402

OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata",
                   "o2_readfetch_fixture.json")

# Families spanning single-region / small / medium / large read counts, incl. the
# three flagship families. (probe counts: 100->77 seqs, 39->195, 2->332[1 region],
# 22->519, 1->1553.)
FAM_IDS = [100, 39, 2, 22, 1]
# Lowered-cap cases (fam_id, cap) to drive truncation; each cap < that family's n_seqs.
LOWERED = [(39, 100), (1, 1000), (100, 50)]


def load_families():
    fams = collections.defaultdict(list)
    with open(pg.FAM_TSV) as f:
        next(f)
        for line in f:
            fi, lid, c, s, e, nr = line.rstrip("\n").split("\t")
            fams[int(fi)].append((lid, c, int(s), int(e), int(nr)))
    return fams


def de_bits(de):
    """Raw IEEE-754 f64 bits of pysam's de (a Python double = widened BAM f32)."""
    # Sanity: de must be an exact single-precision value (BAM 'f' tag), i.e. widening
    # a f32 to f64 is lossless -> Rust `get_tag_float(f32) as f64` reproduces it bit-for-bit.
    assert struct.unpack("<f", struct.pack("<f", de))[0] == de, f"de {de!r} not exact f32"
    return struct.unpack("<Q", struct.pack("<d", de))[0]


def fetch(bam, contig_len, regions, read_cap):
    """VERBATIM copy of the materialize_family read-fetch block (lines 139-161)."""
    best, seqs = {}, {}
    for i, (c, s, e, _) in enumerate(regions):
        s2, e2 = max(0, s - pg.PAD), e + pg.PAD
        for r in bam.fetch(c, s2, min(contig_len.get(c, e2), e2)):
            if r.is_unmapped or r.is_supplementary:
                continue
            if r.reference_end is None or r.reference_end < s or r.reference_start > e:
                continue
            de = dict(r.get_tags()).get("de")
            if de is None:
                continue
            if r.query_name not in best or de < best[r.query_name][1]:
                best[r.query_name] = (i, de)
            if not r.is_secondary and r.query_sequence and r.query_name not in seqs:
                seqs[r.query_name] = r.query_sequence
    precap = list(seqs.keys())
    if len(seqs) > read_cap:
        seqs = dict(list(seqs.items())[:read_cap])
    postcap = list(seqs.keys())
    return best, precap, postcap


def main():
    genome = pysam.FastaFile(pg.GENOME)
    contig_len = dict(zip(genome.references, genome.lengths))
    fams = load_families()
    bam = pysam.AlignmentFile(pg.BAM, "rb")

    out = {
        "bam_path": pg.BAM,
        "genome": pg.GENOME,
        "pad": pg.PAD,
        "read_cap": pg.READ_CAP,
        "contig_len": {c: int(l) for c, l in contig_len.items()},
        "families": [],
    }
    for fid in FAM_IDS:
        regions = pg.dedup_copies(fams[fid])
        regions = sorted(regions, key=lambda r: -(r[2] - r[1]))
        best, precap, postcap = fetch(bam, contig_len, regions, pg.READ_CAP)
        rec = {
            "fam_id": fid,
            "regions": [[c, int(s), int(e), int(nr)] for (c, s, e, nr) in regions],
            "best": {qn: [int(i), str(de_bits(de))] for qn, (i, de) in best.items()},
            "seqs_order_precap": precap,
            "seqs_order_postcap": postcap,
            "n_best": len(best),
            "n_seqs_precap": len(precap),
            "hit_read_cap": len(precap) > pg.READ_CAP,
        }
        out["families"].append(rec)
        print(f"fam {fid}: regions={len(regions)} best={len(best)} "
              f"seqs={len(precap)} (postcap {len(postcap)}, hit_cap={rec['hit_read_cap']})")

    out["lowered_cap"] = []
    for fid, cap in LOWERED:
        regions = pg.dedup_copies(fams[fid])
        regions = sorted(regions, key=lambda r: -(r[2] - r[1]))
        best, precap, postcap = fetch(bam, contig_len, regions, cap)
        out["lowered_cap"].append({
            "fam_id": fid, "read_cap": cap,
            "seqs_order_postcap": postcap, "n_seqs_precap": len(precap),
            "truncated": len(precap) > cap,
        })
        print(f"lowered fam {fid} cap {cap}: precap={len(precap)} "
              f"postcap={len(postcap)} truncated={len(precap) > cap}")

    bam.close()
    genome.close()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=0, sort_keys=True)
    sz = os.path.getsize(OUT)
    print(f"wrote {OUT} ({sz} bytes)")


if __name__ == "__main__":
    main()
