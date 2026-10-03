#!/usr/bin/env python3
"""o3_augment.py — the augmentation step of the pipeline's `candidates` stage (spec 2026-10-02-o3-candidates-design.md §4, §7).

usage: o3_augment.py --fasta G.fa --copies P.fam.copies.tsv --copies-fa P.fam.copies.fa --regions P.fam.copies.regions
                     --cand P.cand --out P.aug

Reads the FLAGGED rows of P.cand.candidates.tsv (o3_candidates) and their exon-union sequences in P.cand.contigs.fa, and
writes, each candidate `cand_<family>_<k>` becoming a contig of its own:
  P.aug.fa            G.fa followed by one record per flagged candidate (its union sequence, one line)
  P.aug.copies.tsv    P.fam.copies.tsv followed by one row per flagged candidate, in the table's own column order:
                      copy_idx after the family's last, tid = chrom = cand_<f>_<k>, start 0, end = locus_end = <len>, one exon
                      0-<len> on +, n_reads 0, max_family_identity 0, source o3_candidate, gene_id ., core_hull NA (the
                      copy_assign parser reads `.` as a malformed hull), sd_depth/core_bp/rep_frac 0, member_status
                      candidate, locus_start 0; any other column NA
  P.aug.copies.fa     P.fam.copies.fa followed by `>{family}|{idx}|cand_<f>_<k>:0-<len>|+|nexon=1` + the union sequence
  P.aug.regions.txt   copy_assign --regions for the candidate families: their P.fam.copies.regions intervals (second
                      column) merged per contig into disjoint intervals (copy_assign refuses overlapping regions, and
                      families' hulls overlap and nest), then `cand_<f>_<k>:0-<len>` per candidate (0-based half-open,
                      as copy_assign reads a region: the BAM query starts at 0 + 1)
  P.aug.families.txt  the candidate families, one id per line, in candidates-table order

Refuses (exit 2, nothing written) when a candidate name already names a sequence of G.fa (G.fa.fai, else G.fa's headers),
when a flagged candidate has no contig or a contig whose length is not its union_len, when a candidate's family has no
row in the copies table or the regions file, and on malformed inputs. Every product is written to `<path>.tmp` and renamed
only once all of them are complete; an I/O error (an unreadable input, a full disk) removes every `.tmp` file and every
product already renamed, and exits 2. Python 3 standard library only.
"""
import argparse
import os
import shutil
import sys

COPY_REQUIRED = ["family_id", "copy_idx", "tid", "chrom", "start", "end", "n_exon", "strand", "n_reads", "exons"]
CAND_REQUIRED = ["family", "candidate", "flagged", "union_len"]


def die(msg):
    print(f"o3_augment: {msg}", file=sys.stderr)
    sys.exit(2)


def read_table(path, required, what):
    """(header, rows as dicts) of a tab-separated table with a header line; `required` columns must be present."""
    try:
        with open(path) as f:
            lines = [l.rstrip("\n") for l in f]
    except OSError as e:
        die(f"cannot read {what} {path}: {e}")
    if not lines or not lines[0]:
        die(f"{what} {path} is empty (expected a header line)")
    head = lines[0].split("\t")
    missing = [c for c in required if c not in head]
    if missing:
        die(f"{what} {path} has no column {', '.join(missing)} (header: {lines[0]!r})")
    rows = []
    for n, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < len(head):
            die(f"{what} {path} line {n} has {len(f)} fields, the header {len(head)}")
        rows.append(dict(zip(head, f)))
    return head, rows


def read_fasta(path, what):
    """{name: sequence} of a FASTA (name = the header's first word)."""
    seqs, name, parts = {}, None, []
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if name is not None:
                        seqs[name] = "".join(parts)
                    name, parts = (line[1:].split() or [""])[0], []
                    if name in seqs:
                        die(f"{what} {path} holds {name} twice")
                elif line:
                    if name is None:
                        die(f"{what} {path} has sequence before its first header")
                    parts.append(line)
    except OSError as e:
        die(f"cannot read {what} {path}: {e}")
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def genome_names(fasta):
    """The sequence names of the genome FASTA, and where they were read from (its .fai when present)."""
    try:
        with open(fasta, "rb") as f:
            if f.read(2) == b"\x1f\x8b":
                die(f"{fasta} is compressed: the augmentation appends contigs to a plain FASTA")
    except OSError as e:
        die(f"cannot read --fasta {fasta}: {e}")
    fai = fasta + ".fai"
    try:
        if os.path.exists(fai):
            with open(fai) as f:
                return {l.split("\t", 1)[0] for l in f if l.strip()}, fai
        names = set()
        with open(fasta) as f:
            for line in f:
                if line.startswith(">"):
                    names.add((line[1:].split() or [""])[0])
    except OSError as e:
        die(f"cannot read the sequence names of --fasta {fasta}: {e}")
    return names, f"the headers of {fasta}"


def parse_region(text, where):
    """chrom, lo, hi of `chrom:lo-hi` (the chromosome is everything before the LAST colon)."""
    chrom, sep, span = text.rpartition(":")
    lo, dash, hi = span.partition("-")
    if not sep or not dash or not chrom or not lo.isdigit() or not hi.isdigit():
        die(f"bad region {text!r} in {where} (expected chrom:start-end)")
    return chrom, int(lo), int(hi)


def merge_regions(intervals):
    """Disjoint `chrom:lo-hi` strings covering `intervals` [(chrom, lo, hi)]: per contig, sorted, an interval that
    overlaps or touches the previous one is merged into it (the driver's fam_regions applies the same rule)."""
    out = []
    for chrom, lo, hi in sorted(intervals):
        if out and out[-1][0] == chrom and lo <= out[-1][2]:
            out[-1][2] = max(out[-1][2], hi)
        else:
            out.append([chrom, lo, hi])
    return [f"{c}:{lo}-{hi}" for c, lo, hi in out]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for flag in ("--fasta", "--copies", "--copies-fa", "--regions", "--cand", "--out"):
        ap.add_argument(flag, required=True)
    a = ap.parse_args()

    # the flagged candidates and their union sequences
    _, cand_rows = read_table(f"{a.cand}.candidates.tsv", CAND_REQUIRED, "candidates table")
    flagged = [r for r in cand_rows if r["flagged"] == "1"]
    if not flagged:
        die(f"{a.cand}.candidates.tsv has no flagged candidate: nothing to augment")
    contigs = read_fasta(f"{a.cand}.contigs.fa", "contigs")
    names, names_src = genome_names(a.fasta)
    seen = set()
    for r in flagged:
        cid = r["candidate"]
        if cid in names:
            die(f"candidate {cid} already names a sequence of {a.fasta} ({names_src}); refusing to write {a.out}.*")
        if cid in seen:
            die(f"candidate {cid} is flagged twice in {a.cand}.candidates.tsv")
        seen.add(cid)
        if cid not in contigs:
            die(f"flagged candidate {cid} has no record in {a.cand}.contigs.fa")
        if not r["union_len"].isdigit() or len(contigs[cid]) != int(r["union_len"]) or not contigs[cid]:
            die(f"{cid}: contig of {len(contigs[cid])} bp but union_len {r['union_len']!r} in {a.cand}.candidates.tsv")

    # the copy table: candidate rows continue each family's copy numbering
    head, copy_rows = read_table(a.copies, COPY_REQUIRED, "copies table")
    next_idx = {}
    for r in copy_rows:
        if not r["copy_idx"].isdigit():
            die(f"{a.copies}: bad copy_idx {r['copy_idx']!r} for {r['family_id']}")
        next_idx[r["family_id"]] = max(next_idx.get(r["family_id"], 0), int(r["copy_idx"]) + 1)
    new_rows, new_fa, cand_regions, families = [], [], [], []
    for r in flagged:
        fid, cid, n = r["family"], r["candidate"], len(contigs[r["candidate"]])
        if fid not in next_idx:
            die(f"candidate {cid}: family {fid} has no row in {a.copies}")
        idx = next_idx[fid]
        next_idx[fid] += 1
        val = {"family_id": fid, "copy_idx": str(idx), "tid": cid, "chrom": cid, "start": "0", "end": str(n), "n_exon": "1",
               "strand": "+", "n_reads": "0", "exons": f"0-{n}", "max_family_identity": "0", "source": "o3_candidate",
               "gene_id": ".", "core_hull": "NA", "sd_depth": "0", "core_bp": "0", "rep_frac": "0",
               "member_status": "candidate", "locus_start": "0", "locus_end": str(n)}
        new_rows.append("\t".join(val.get(c, "NA") for c in head))
        new_fa.append(f">{fid}|{idx}|{cid}:0-{n}|+|nexon=1\n{contigs[cid]}\n")
        cand_regions.append(f"{cid}:0-{n}")
        if fid not in families:
            families.append(fid)

    # the regions of the candidate families, merged, then the candidate contigs
    fam_set, intervals, with_region = set(families), [], set()
    try:
        with open(a.regions) as f:
            for n, line in enumerate(f, start=1):
                if not line.strip():
                    continue
                fs = line.rstrip("\n").split("\t")
                if len(fs) < 2:
                    die(f"{a.regions} line {n}: expected `family<TAB>chrom:start-end`")
                if fs[0] in fam_set:
                    intervals.append(parse_region(fs[1], f"{a.regions} line {n}"))
                    with_region.add(fs[0])
    except OSError as e:
        die(f"cannot read --regions {a.regions}: {e}")
    lost = [fid for fid in families if fid not in with_region]
    if lost:
        die(f"{a.regions} has no region for the candidate famil{'y' if len(lost) == 1 else 'ies'} {', '.join(lost)}")

    # write every product to .tmp, then rename them all; an I/O error (an unreadable --fasta / --copies-fa, a full disk, an
    # unwritable --out directory) removes every .tmp file and every product this run already renamed, then exits 2
    tmp, done = {}, []

    def out(suffix):
        tmp[suffix] = f"{a.out}.{suffix}.tmp"
        return tmp[suffix]

    try:
        with open(out("fa"), "wb") as w, open(a.fasta, "rb") as g:
            shutil.copyfileobj(g, w, 1 << 22)
            size = g.tell()
            if size:
                g.seek(size - 1)
                if g.read(1) != b"\n":
                    w.write(b"\n")
            for r in flagged:
                w.write(f">{r['candidate']}\n{contigs[r['candidate']]}\n".encode())
        for suffix, src, extra in (("copies.tsv", a.copies, [l + "\n" for l in new_rows]), ("copies.fa", a.copies_fa, new_fa)):
            with open(src) as s, open(out(suffix), "w") as w:
                text = s.read()
                w.write(text if not text or text.endswith("\n") else text + "\n")
                w.writelines(extra)
        with open(out("regions.txt"), "w") as w:
            w.writelines(l + "\n" for l in merge_regions(intervals) + cand_regions)
        with open(out("families.txt"), "w") as w:
            w.writelines(fid + "\n" for fid in families)
        for suffix, path in tmp.items():
            os.replace(path, f"{a.out}.{suffix}")
            done.append(f"{a.out}.{suffix}")
    except OSError as e:
        for path in list(tmp.values()) + done:
            try:
                os.remove(path)
            except OSError:
                pass
        die(f"cannot make {a.out}.*: {e} (nothing written)")
    print(f"o3_augment: {len(flagged)} candidate contig(s) of {len(families)} famil{'y' if len(families) == 1 else 'ies'} "
          f"-> {a.out}.{{fa,copies.tsv,copies.fa,regions.txt,families.txt}}", file=sys.stderr)


if __name__ == "__main__":
    main()
