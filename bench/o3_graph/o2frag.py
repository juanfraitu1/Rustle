#!/usr/bin/env python3
"""Real-read fragments with a truth (2026-10-19, the O2 regimen's real-read test). A real library's reads that minimap2
places clearly on one copy (MAPQ 60) are cut into seeded windows; a window keeps the read's real sequence -- its sequencing
errors and the donor's own alleles -- but spans fewer PSVs, so realigned it is often almost tied. Its truth is its full
read's copy, kept only when the full read, realigned with the windows, still lands on that copy clear by MARGIN AS of any
other locus (on HG002, reads clear by >= 30 agree with the donor's assembly 98.4-99.8%).

    o2frag.py sources BAM COPIES.tsv OUT_DIR [--per-family 60] [--seed 1] [--min-len 1000]
                                         -> OUT_DIR/sources.fa (original orientation), sources.tsv (read, copy, family)
    o2frag.py fragments OUT_DIR [--per-read 3] [--min 300] [--max 1200] [--seed 1]  -> OUT_DIR/fragments.fa (read|f<i>|start|len)
    o2frag.py truth OUT_DIR ALN.bam COPIES.tsv [--margin 30]  -> OUT_DIR/truth.tsv (read, copy): every fragment of a kept source
"""
import collections
import csv
import hashlib
import importlib.util
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_frag", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def windows(name, length, k, lo, hi, seed):
    """k seeded windows (start, length) of a read: lengths uniform in [lo, min(hi, length)], starts uniform; none when the read
    is shorter than lo."""
    if length < lo:
        return []
    rng = random.Random(f"{seed}:{name}")
    out = []
    for _ in range(k):
        n = rng.randint(lo, min(hi, length))
        out.append((rng.randint(0, length - n), n))
    return out


def source_truth(alns, copy, margin):
    """The source read's copy when its primary lies on `copy` and every alignment on another copy or off every copy scores at
    least `margin` AS below it; else None. alns: [{AS, copy, prim}]."""
    prim = next((a for a in alns if a["prim"]), None)
    if prim is None or prim["copy"] != copy:
        return None
    rivals = [a["AS"] for a in alns if a["copy"] != copy]
    return copy if not rivals or max(rivals) <= prim["AS"] - margin else None


def source_of(name):
    return name.split("|")[0] if "|" in name else None


def _rank(seed, name):
    return hashlib.sha1(f"{seed}:{name}".encode()).hexdigest()


def _neg(rank):
    """A key that orders ranks in reverse, and its own inverse (heapq is a min-heap)."""
    return "".join(chr(255 - ord(ch)) for ch in rank)


def keep_lowest(heap, n, rank, item):
    """Keep in heap the n items of lowest rank seen so far (a deterministic sample by a seeded hash)."""
    import heapq
    if len(heap) < n:
        heapq.heappush(heap, (_neg(rank), item))
    elif rank < _neg(heap[0][0]):
        heapq.heapreplace(heap, (_neg(rank), item))


def cmd_sources(bam_p, copies_p, out_dir, per_family=60, seed=1, min_len=1000):
    import pysam
    P, S = _sibling("o2psv"), _sibling("summary")
    cps = P.load_copies(copies_p)
    fam_of = {c: c.rsplit(".", 1)[0] for v in cps.values() for c, *_ in v}
    got = collections.defaultdict(list)                  # family -> keep_lowest heap of (read, copy, sequence)
    seen = set()
    with pysam.AlignmentFile(bam_p) as bam:
        for chrom, v in cps.items():
            for c, strand, ex, (s, e) in v:
                for rec in bam.fetch(chrom, s, e):
                    if (rec.is_secondary or rec.is_supplementary or rec.is_unmapped or rec.mapping_quality < 60
                            or rec.query_sequence is None or rec.infer_read_length() < min_len or rec.query_name in seen):
                        continue
                    if S.attribute(rec.get_blocks(), P._strand(rec), v) != c:
                        continue
                    seen.add(rec.query_name)
                    q = rec.query_sequence
                    keep_lowest(got[fam_of[c]], per_family, _rank(seed, rec.query_name),
                                (rec.query_name, c, q[::-1].translate(COMP) if rec.is_reverse else q))
    os.makedirs(out_dir, exist_ok=True)
    n = 0
    with open(f"{out_dir}/sources.fa", "w") as fa, open(f"{out_dir}/sources.tsv", "w") as tsv:
        tsv.write("read\tcopy\tfamily\n")
        for f in sorted(got):
            for _, (r, c, seq) in sorted((_neg(k), x) for k, x in got[f]):
                fa.write(f">{r}\n{seq}\n")
                tsv.write(f"{r}\t{c}\t{f}\n")
                n += 1
    print(f"[o2frag] {n} source reads from {len(got)} families -> {out_dir}/sources.fa", flush=True)
    return 0


def _fasta(p):
    name, seq = None, []
    for ln in open(p):
        if ln.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = ln[1:].strip().split()[0], []
        else:
            seq.append(ln.strip())
    if name:
        yield name, "".join(seq)


def cmd_fragments(out_dir, per_read=3, lo=300, hi=1200, seed=1):
    n = 0
    with open(f"{out_dir}/fragments.fa", "w") as fh:
        for name, seq in _fasta(f"{out_dir}/sources.fa"):
            for i, (s, k) in enumerate(windows(name, len(seq), per_read, lo, hi, seed)):
                fh.write(f">{name}|f{i}|{s}|{k}\n{seq[s:s + k]}\n")
                n += 1
    print(f"[o2frag] {n} fragments -> {out_dir}/fragments.fa", flush=True)
    return 0


def cmd_truth(out_dir, aln_p, copies_p, margin=30):
    import pysam
    P, S = _sibling("o2psv"), _sibling("summary")
    cps = P.load_copies(copies_p)
    src = {r["read"]: r["copy"] for r in csv.DictReader(open(f"{out_dir}/sources.tsv"), delimiter="\t")}
    alns = collections.defaultdict(list)
    frags = collections.defaultdict(list)
    with pysam.AlignmentFile(aln_p) as bam:
        for rec in bam.fetch(until_eof=True):
            q = rec.query_name
            if q in src and not rec.is_unmapped and not rec.is_supplementary and rec.has_tag("AS"):
                alns[q].append(dict(AS=rec.get_tag("AS"), prim=not rec.is_secondary,
                                    copy=S.attribute(rec.get_blocks(), P._strand(rec), cps.get(rec.reference_name, ()))))
            elif source_of(q) in src and not rec.is_secondary and not rec.is_supplementary:
                frags[source_of(q)].append(q)
    kept = {q: source_truth(a, src[q], margin) for q, a in alns.items()}
    n = collections.Counter("kept" if v else "dropped" for v in kept.values())
    with open(f"{out_dir}/truth.tsv", "w") as fh:
        fh.write("read\tcopy\n")
        for q, c in sorted(kept.items()):
            if c:
                for f in sorted(set(frags.get(q, ()))):
                    fh.write(f"{f}\t{c}\n")
                    n["fragments"] += 1
    print(f"[o2frag] sources {dict(n)} -> {out_dir}/truth.tsv", flush=True)
    return 0


def _opt(argv, name, default, cast=int):
    return cast(argv[argv.index(name) + 1]) if name in argv else default


def main(argv):
    if argv[:1] == ["sources"] and len(argv) >= 4:
        return cmd_sources(argv[1], argv[2], argv[3], _opt(argv, "--per-family", 60), _opt(argv, "--seed", 1),
                           _opt(argv, "--min-len", 1000))
    if argv[:1] == ["fragments"] and len(argv) >= 2:
        return cmd_fragments(argv[1], _opt(argv, "--per-read", 3), _opt(argv, "--min", 300), _opt(argv, "--max", 1200),
                             _opt(argv, "--seed", 1))
    if argv[:1] == ["truth"] and len(argv) >= 4:
        return cmd_truth(argv[1], argv[2], argv[3], _opt(argv, "--margin", 30))
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
