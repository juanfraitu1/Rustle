#!/usr/bin/env python3
"""Real reads in the integrated pipeline's world format (2026-10-18; user: "can we prove that O2 and O3 work? Maybe using the
YAGs or the mainstays?"). A real library's reads of a family set, with a truth from their unique placement, become a world
like the simulator's -- realworld.py build --real-reads -- so fastset runs and scores it unchanged.

Truth: a read whose primary alignment (in the library's own CHM13 alignments) is placed uniquely (MAPQ >= MIN_MAPQ) and whose
aligned blocks overlap a copy's exons belongs to the copy it overlaps most; every other read is in the world but unscored.
Masking test (REMOVED with real reads): per family, one coding copy with >= MIN_READS truth reads, chosen with a seed, is masked
from the reference; its truth reads must be flagged by O3 or abstained by O2, not placed on a present copy.

    realreads.py pool BAM FAMILIES.json SET OUT_DIR      -> OUT_DIR/reads.fq, truth.tsv (read, copy, family), counts.json
    realreads.py choose OUT_DIR FAMILIES.json SET SEED   -> OUT_DIR/masked_<SEED>.json ({family: gene_id})
"""
import collections
import json
import os
import random
import sys

MIN_MAPQ, MIN_READS, PAD = 60, 20, 1000


def copy_exons(fams):
    """{chrom: [(start, end, copy id)]} -- every copy's exons (all its isoforms), copy ids as the simulator numbers them."""
    out = collections.defaultdict(list)
    for fam in fams:
        for k, c in enumerate(fam["copies"]):
            for s, e in sorted({tuple(x) for v in c["transcripts"].values() for x in v}):
                out[c["chrom"]].append((s, e, f"{fam['family']}.{k}"))
    return out


def truth_copy(chrom, blocks, mapq, idx, min_mapq=MIN_MAPQ):
    """The copy a uniquely placed read belongs to: the one whose exons its aligned blocks overlap most; None when the
    placement is not unique or lies on no exon."""
    if mapq < min_mapq:
        return None
    ov = collections.Counter()
    for s, e in blocks:
        for s2, e2, cid in idx.get(chrom, ()):
            if s < e2 and s2 < e:
                ov[cid] += min(e, e2) - max(s, s2)
    return max(sorted(ov), key=lambda c: ov[c]) if ov else None


def choose_masked(fams, counts, seed, min_reads=MIN_READS):
    """{family: gene_id of the copy to mask}: a seeded choice among the family's CODING copies with >= min_reads truth reads
    (the simulator's REMOVED masks a coding copy too); a family with none is left out of the test."""
    out = {}
    for fam in fams:
        ok = [c for k, c in enumerate(fam["copies"]) if c.get("biotype", "protein_coding") == "protein_coding"
              and counts.get(f"{fam['family']}.{k}", 0) >= min_reads]
        if ok:
            out[fam["family"]] = random.Random(f"{seed}:{fam['family']}:masked").choice(ok)["gene_id"]
    return out


def cmd_pool(bam_p, fam_p, set_name, out_dir):
    """Every primary alignment overlapping a copy's span (+- PAD) -> its read (FASTQ) and, if uniquely placed on an exon,
    its truth copy."""
    import pysam
    fams = json.load(open(fam_p))["sets"][set_name]
    idx = copy_exons(fams)
    fam_of = {f"{f['family']}.{k}": f["family"] for f in fams for k in range(len(f["copies"]))}
    os.makedirs(out_dir, exist_ok=True)
    seen, counts = set(), collections.Counter()
    with pysam.AlignmentFile(bam_p) as bam, open(f"{out_dir}/reads.fq", "w") as fq, open(f"{out_dir}/truth.tsv", "w") as tr:
        tr.write("read\tcopy\tfamily\n")
        for fam in fams:
            for c in fam["copies"]:
                for rec in bam.fetch(c["chrom"], max(0, c["start"] - PAD), c["end"] + PAD):
                    if rec.is_secondary or rec.is_supplementary or rec.is_unmapped or rec.query_name in seen:
                        continue
                    if rec.query_sequence is None:
                        continue
                    seen.add(rec.query_name)
                    seq, qual = rec.query_sequence, rec.query_qualities
                    if rec.is_reverse:                       # the read as sequenced
                        comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
                        seq = seq[::-1].translate(comp)
                        qual = qual[::-1] if qual is not None else None
                    q = "".join(chr(x + 33) for x in qual) if qual is not None else "~" * len(seq)
                    fq.write(f"@{rec.query_name}\n{seq}\n+\n{q}\n")
                    cid = truth_copy(rec.reference_name, rec.get_blocks(), rec.mapping_quality, idx)
                    if cid:
                        tr.write(f"{rec.query_name}\t{cid}\t{fam_of[cid]}\n")
                        counts[cid] += 1
    json.dump(dict(reads=len(seen), truth=sum(counts.values()), per_copy=dict(sorted(counts.items()))),
              open(f"{out_dir}/counts.json", "w"), indent=1)
    print(f"[realreads] {out_dir}: {len(seen)} reads, {sum(counts.values())} with a truth copy", flush=True)
    return 0


def cmd_choose(out_dir, fam_p, set_name, seed):
    fams = json.load(open(fam_p))["sets"][set_name]
    counts = json.load(open(f"{out_dir}/counts.json"))["per_copy"]
    got = choose_masked(fams, counts, seed)
    json.dump(got, open(f"{out_dir}/masked_{seed}.json", "w"), indent=1)
    print(f"[realreads] {len(got)} of {len(fams)} families have a copy to mask -> {out_dir}/masked_{seed}.json", flush=True)
    return 0


def main(argv):
    if len(argv) == 5 and argv[0] == "pool":
        return cmd_pool(*argv[1:])
    if len(argv) == 5 and argv[0] == "choose":
        return cmd_choose(argv[1], argv[2], argv[3], int(argv[4]))
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
