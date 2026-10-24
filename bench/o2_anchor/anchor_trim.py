#!/usr/bin/env python3
"""O2 anchor-trim truth (docs/archive/2026-10/merged/PREREG_o2_anchor_2026-10-13.md, arm 1): real reads placed uniquely (MAPQ 60) on a catalog
copy lose the part of the read that is not shared with a sibling copy; O2 must recover the copy from the shared part.

    anchor_trim.py CATALOG_PREFIX BAM GENOME_FA OUTDIR [PER_COPY] [MAX_FAMILIES] [MAX_LINES]

Per catalog family with >= 2 copies: primary MAPQ-60 reads (-F 2308 -q 60) whose reference span overlaps exactly one family
(truth copy = the family copy with the largest overlap), at most PER_COPY per copy (lowest md5 of the read name). Each read
(original orientation) is aligned to the genomic spans of its family's copies (minimap2 -c -x splice:hq -uf, secondaries
kept); the trimmed read = the read interval covered by BOTH the truth copy's best hit and the best sibling hit, kept iff
>= MIN_SHARED bp long and >= MIN_REMOVED bp were removed. Writes OUTDIR/trimmed.fa (names `<read>|trim`) and
OUTDIR/truth.tsv. MAX_FAMILIES (0 = all) keeps the families with the lowest md5(family_id); MAX_LINES (0 = all) reads at most
that many records per copy before sampling (Amendment 1b, human held-out only).
"""
import collections, hashlib, os, re, subprocess, sys

MIN_SHARED = 300
MIN_REMOVED = 100
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def original_seq(seq, flag):
    """BAM SEQ is reference-forward; a reverse-strand record's original read is its reverse complement."""
    return seq.translate(COMP)[::-1] if int(flag) & 16 else seq


def ref_end(pos0, cigar):
    """0-based exclusive reference end of an alignment starting at pos0."""
    return pos0 + sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MDN=X")


def truth_copy(copies, chrom, s, e):
    """copies = [(copy_id, chrom, start, end)]: the copy overlapping [s, e) the most, None if none overlaps."""
    best, bo = None, 0
    for cid, c, cs, ce in copies:
        if c == chrom:
            o = min(e, ce) - max(s, cs)
            if o > bo:
                best, bo = cid, o
    return best


def shared_interval(hits, truth, min_shared=MIN_SHARED):
    """hits = PAF records of ONE read against its family's copies. -> (qs, qe, sibling) = the read interval covered by the
    truth copy's best hit and the best sibling hit (most matching bases), or None."""
    th = [h for h in hits if h["tname"] == truth]
    sh = [h for h in hits if h["tname"] != truth]
    if not th or not sh:
        return None
    t = max(th, key=lambda h: h["nmatch"])
    s = max(sh, key=lambda h: h["nmatch"])
    qs, qe = max(t["qs"], s["qs"]), min(t["qe"], s["qe"])
    if qe - qs < min_shared:
        return None
    return qs, qe, s["tname"]


def tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in fh]


def paf(text):
    out = []
    for l in text.splitlines():
        f = l.split("\t")
        if len(f) >= 12:
            out.append({"qname": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "tname": f[5],
                        "nmatch": int(f[9])})
    return out


def region_lines(bam, region, max_lines=0):
    """Primary MAPQ-60 records of a region, streamed; at most max_lines (0 = all)."""
    p = subprocess.Popen(["samtools", "view", "-F", "2308", "-q", "60", bam, region], stdout=subprocess.PIPE, text=True)
    out = []
    for line in p.stdout:
        out.append(line.rstrip("\n"))
        if max_lines and len(out) >= max_lines:
            p.kill()
            break
    p.stdout.close()
    rc = p.wait()
    if rc not in (0, -9) or (rc == -9 and not max_lines):
        raise RuntimeError(f"samtools view {region} exited {rc}")
    return out


def pick_families(fams, max_families):
    """{family: copies} -> the families with >= 2 copies, the max_families lowest md5(family_id) if max_families > 0."""
    keep = sorted((f for f, cs in fams.items() if len(cs) >= 2), key=lambda f: hashlib.md5(f.encode()).hexdigest())
    return set(keep[:max_families] if max_families > 0 else keep)


def main(prefix, bam, genome, out, per_copy=20, max_families=0, max_lines=0):
    per_copy, max_families, max_lines = int(per_copy), int(max_families), int(max_lines)
    os.makedirs(out, exist_ok=True)
    rows = tsv(prefix + ".fam.copies.tsv")
    fams = collections.defaultdict(list)
    for r in rows:
        fams[r["family_id"]].append((r["copy_idx"], r["chrom"], int(r["start"]), int(r["end"])))
    by_chrom = collections.defaultdict(list)
    for f, cs in fams.items():
        for cid, c, s, e in cs:
            by_chrom[c].append((f, s, e))
    n = collections.Counter()
    tf, tt = open(os.path.join(out, "trimmed.fa"), "w"), open(os.path.join(out, "truth.tsv"), "w")
    tt.write("read_name\torig_read\tfamily_id\ttruth_copy\tsibling_copy\tread_len\tshared_len\tremoved\tchrom\tstart\tend\n")
    chosen = pick_families(fams, max_families)
    for fam, cs in sorted(fams.items()):
        if fam not in chosen:
            continue
        picked = {}
        for cid, c, s, e in cs:
            cand = []
            for line in region_lines(bam, f"{c}:{s + 1}-{e}", max_lines):
                f = line.split("\t")
                q, flag, chrom, p0, cig, seq = f[0], f[1], f[2], int(f[3]) - 1, f[5], f[9]
                re_ = ref_end(p0, cig)
                if truth_copy(cs, chrom, p0, re_) != cid:
                    continue  # collected under its own truth copy
                if any(g != fam and gs < re_ and p0 < ge for g, gs, ge in by_chrom[chrom]):
                    n["skip_two_families"] += 1
                    continue
                cand.append((hashlib.md5(q.encode()).hexdigest(), q, original_seq(seq, flag), chrom, p0, re_))
            for x in sorted(cand)[:per_copy]:
                picked[x[1]] = (cid,) + x[2:]
        if not picked:
            continue
        n["reads_picked"] += len(picked)
        rfa, tfa = os.path.join(out, "_reads.fa"), os.path.join(out, "_targets.fa")
        with open(rfa, "w") as fh:
            fh.writelines(f">{q}\n{v[1]}\n" for q, v in picked.items())
        with open(tfa, "w") as fh:
            for cid, c, s, e in cs:
                seq = "".join(subprocess.run(["samtools", "faidx", genome, f"{c}:{s + 1}-{e}"], capture_output=True,
                                             text=True, check=True).stdout.splitlines()[1:])
                fh.write(f">{cid}\n{seq}\n")
        hits = collections.defaultdict(list)
        for h in paf(subprocess.run(["minimap2", "-c", "-x", "splice:hq", "-uf", "--secondary=yes", "-N", "50", "-p",
                                     "0.1", "-t", "2", tfa, rfa], capture_output=True, text=True, check=True).stdout):
            hits[h["qname"]].append(h)
        for q, (cid, seq, chrom, p0, re_) in sorted(picked.items()):
            si = shared_interval(hits.get(q, []), cid)
            if si is None:
                n["no_shared_interval"] += 1
                continue
            qs, qe, sib = si
            if len(seq) - (qe - qs) < MIN_REMOVED:
                n["nothing_to_remove"] += 1
                continue
            n["trimmed"] += 1
            tf.write(f">{q}|trim\n{seq[qs:qe]}\n")
            tt.write(f"{q}|trim\t{q}\t{fam}\t{cid}\t{sib}\t{len(seq)}\t{qe - qs}\t{len(seq) - (qe - qs)}\t{chrom}\t{p0}\t{re_}\n")
        n["families_with_reads"] += 1
    tf.close()
    tt.close()
    for f in ("_reads.fa", "_targets.fa"):
        p = os.path.join(out, f)
        if os.path.exists(p):
            os.remove(p)
    print(dict(n))


if __name__ == "__main__":
    main(*sys.argv[1:8])
