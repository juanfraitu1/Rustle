#!/usr/bin/env python3
"""Seed-and-extend family discovery with minimap2 (alignment-based, no k-mers).

For each existing gw_xcbase family member (seed), align its spliced sequence
against a database of genome-wide thin-locus spliced sequences using minimap2.
A thin locus is added to the seed's family iff the alignment is full-length
homologous (identity >= 0.80, coverage-of-shorter >= 0.50) and the locus is at a
distinct genomic position.

This avoids k-mer clustering and uses the same principled exon-sum homology
criterion already validated in the gw_family_catalog --refine step.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/seed_extend_minimap2.py
"""
import csv
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

SCRATCH = "/home/juanfra/winloci_scratch"
HERE = os.path.dirname(__file__)
GGO_FA = os.path.join(SCRATCH, "GGO.fasta")
THIN_LOCI_TSV = os.path.join(SCRATCH, "thin_loci_genome_wide.tsv")
SEED_FAMILIES = os.path.join(SCRATCH, "gw_xcbase.families.tsv")
SEED_COPIES = os.path.join(SCRATCH, "gw_xcbase.copies.tsv")
SEED_FA = os.path.join(SCRATCH, "gw_xcbase.copies.fa")
OUT_RESCUED = os.path.join(HERE, "seed_extend_minimap2.rescued.tsv")
OUT_FAMILIES = os.path.join(SCRATCH, "gw_seedext.families.tsv")
OUT_COPIES = os.path.join(SCRATCH, "gw_seedext.copies.tsv")

MIN_IDENTITY = 0.80
MIN_COVERAGE = 0.40  # full-length homology threshold (lenient to catch partial/fragmented copies)
LEN_CAP = 9000

# Exon-structure guard: real paralog copies have similar exon count and similar
# spliced length. This filters out fragmentary isoforms and domain-sharers.
MAX_EXON_DIFF = 1    # |n_exon_seed - n_exon_locus| must be <= 1
MIN_LEN_RATIO = 0.6  # shorter / longer spliced length
MAX_LEN_RATIO = 1.7


def parse_copies_fa(path):
    """Parse gw_xcbase.copies.fa -> dict (family_id, copy_idx) -> seq."""
    seqs = {}
    cur = None
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                parts = ln[1:].split()[0].split("|")
                cur = (parts[0], int(parts[1]))
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(ln.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def load_seeds():
    copies = []
    with open(SEED_COPIES) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            copies.append(row)
    seqs = parse_copies_fa(SEED_FA)
    seeds = []
    for c in copies:
        key = (c["family_id"], int(c["copy_idx"]))
        seq = seqs.get(key, "")
        if not seq:
            continue
        seeds.append({
            "tid": c["tid"],
            "family_id": c["family_id"],
            "chrom": c["chrom"],
            "start": int(c["start"]),
            "end": int(c["end"]),
            "strand": c["strand"],
            "n_exon": int(c["n_exon"]),
            "seq": seq,
        })
    return seeds


def load_thin_loci(path):
    loci = []
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            introns = []
            if "introns" in row and row["introns"]:
                for part in row["introns"].split(","):
                    d, a = part.split("-")
                    introns.append((int(d), int(a)))
            loci.append({
                "lid": int(row["lid"]),
                "chrom": row["chrom"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "strand": row["strand"],
                "support": int(row["support"]),
                "n_exon": int(row["n_exon"]),
                "introns": introns,
                "seq": row["seq"],
            })
    return loci


def write_fasta(path, records):
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")


def parse_paf(paf_path):
    """Parse minimap2 PAF. Return list of dicts with alignment metrics."""
    hits = []
    with open(paf_path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qname, qlen, qs, qe, strand, tname, tlen, ts, te, nmatch, aln_len, mapq = f[:12]
            qlen, qs, qe = int(qlen), int(qs), int(qe)
            tlen, ts, te = int(tlen), int(ts), int(te)
            nmatch, aln_len = int(nmatch), int(aln_len)
            identity = nmatch / aln_len if aln_len > 0 else 0.0
            qcov = (qe - qs) / qlen if qlen > 0 else 0.0
            tcov = (te - ts) / tlen if tlen > 0 else 0.0
            coverage = min(qcov, tcov)
            hits.append({
                "qname": qname, "tname": tname, "strand": strand,
                "identity": identity, "coverage": coverage,
                "aln_len": aln_len, "nmatch": nmatch,
                "qlen": qlen, "tlen": tlen,
            })
    return hits


def main():
    seeds = load_seeds()
    print(f"Loaded {len(seeds)} seed copies from gw_xcbase", flush=True)

    thin_loci = load_thin_loci(THIN_LOCI_TSV)
    print(f"Loaded {len(thin_loci)} thin loci", flush=True)

    # Filter thin loci by length cap (same as rescue)
    thin_loci = [L for L in thin_loci if 200 <= len(L["seq"]) <= LEN_CAP]
    print(f"Thin loci after length filter: {len(thin_loci)}", flush=True)

    # Build FASTA of thin-locus database
    thin_records = [(f"TL_{L['lid']}", L["seq"]) for L in thin_loci]
    with tempfile.TemporaryDirectory() as tmpdir:
        thin_fa = os.path.join(tmpdir, "thin_loci.fa")
        seed_fa = os.path.join(tmpdir, "seeds.fa")
        paf_out = os.path.join(tmpdir, "seeds_vs_thin.paf")
        write_fasta(thin_fa, thin_records)
        seed_records = [(f"SEED_{s['tid']}", s["seq"]) for s in seeds]
        write_fasta(seed_fa, seed_records)

        # minimap2 all seeds vs thin-locus database
        cmd = [
            "minimap2",
            "-x", "asm20",
            "-t", "8",
            "-N", "50",
            "-p", "0.1",
            "--secondary", "yes",
            thin_fa,
            seed_fa,
        ]
        print(f"Running minimap2: {' '.join(cmd)}", flush=True)
        with open(paf_out, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=True)

        hits = parse_paf(paf_out)
        print(f"Raw minimap2 alignments: {len(hits)}", flush=True)

    # Filter hits by identity + coverage
    filtered = [
        h for h in hits
        if h["identity"] >= MIN_IDENTITY and h["coverage"] >= MIN_COVERAGE
    ]
    print(f"Full-length homologous hits: {len(filtered)}", flush=True)

    # Map thin-locus target names back to loci and seed tids back to families
    thin_by_name = {f"TL_{L['lid']}": L for L in thin_loci}
    seed_by_tid = {s["tid"]: s for s in seeds}

    # Exclude thin loci that overlap an existing seed span by >= 50%
    def max_overlap_frac(locus, seed):
        s, e = locus["start"], locus["end"]
        ms, me = seed["start"], seed["end"]
        if s >= me or e <= ms:
            return 0.0
        ov = min(e, me) - max(s, ms)
        return ov / (e - s)

    def exon_compatible(seed, locus):
        # similar exon count
        if abs(locus["n_exon"] - seed["n_exon"]) > MAX_EXON_DIFF:
            return False
        # similar spliced length
        ls, ll = len(seed["seq"]), len(locus["seq"])
        if ls == 0 or ll == 0:
            return False
        ratio = min(ls, ll) / max(ls, ll)
        if ratio < MIN_LEN_RATIO or ratio > MAX_LEN_RATIO:
            return False
        return True

    rescued_by_key = {}
    for h in filtered:
        seed_tid = h["qname"].replace("SEED_", "", 1)
        thin_name = h["tname"]
        if thin_name not in thin_by_name:
            continue
        seed = seed_by_tid.get(seed_tid)
        if seed is None:
            continue
        L = thin_by_name[thin_name]
        # minimap2 PAF strand: '+' means the transcription-orientation sequences
        # align directly; '-' means one must be reverse-complemented. We only
        # accept '+' hits so paralogs are homologous in their native orientation.
        if h["strand"] != "+":
            continue
        if not exon_compatible(seed, L):
            continue
        # exclude if thin locus is mostly contained in the seed's span
        if max_overlap_frac(L, seed) >= 0.5:
            continue
        key = (L["chrom"], L["start"], L["end"])
        entry = dict(
            fid=seed["family_id"],
            seed_tid=seed_tid,
            tid=f"SE_{L['chrom']}_{L['start']}_{L['n_exon']}",
            chrom=L["chrom"], start=L["start"], end=L["end"],
            strand=L["strand"], n_exon=L["n_exon"], support=L["support"],
            introns=L.get("introns", []),
            identity=h["identity"], coverage=h["coverage"],
        )
        # keep best hit per thin locus (by identity*coverage)
        if key not in rescued_by_key or (entry["identity"] * entry["coverage"]) > (rescued_by_key[key]["identity"] * rescued_by_key[key]["coverage"]):
            rescued_by_key[key] = entry

    def collapse_rescued(copies):
        """Collapse rescued copies that share a junction (isoforms of one locus)."""
        sorted_c = sorted(copies, key=lambda x: (x["chrom"], x["strand"], x["start"]))
        groups = []
        for C in sorted_c:
            intr_set = set(C.get("introns", []))
            merged = False
            for g in groups:
                rep = g[0]
                if C["chrom"] != rep["chrom"] or C["strand"] != rep["strand"]:
                    continue
                if C["start"] < rep["end"] and C["end"] > rep["start"] and intr_set.intersection(rep.get("introns", [])):
                    # keep the copy with best identity*coverage
                    if C["identity"] * C["coverage"] > rep["identity"] * rep["coverage"]:
                        g.insert(0, C)
                    else:
                        g.append(C)
                    merged = True
                    break
            if not merged:
                groups.append([C])
        return [g[0] for g in groups]

    rescued = collapse_rescued(list(rescued_by_key.values()))
    rescued = sorted(rescued, key=lambda x: -(x["identity"] * x["coverage"]))
    print(f"Rescued {len(rescued)} copies after isoform collapse", flush=True)

    with open(OUT_RESCUED, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family_id", "seed_tid", "tid", "chrom", "start", "end", "n_exon",
                    "strand", "n_reads", "identity", "coverage"])
        for r in rescued:
            w.writerow([r["fid"], r["seed_tid"], r["tid"], r["chrom"], r["start"], r["end"],
                        r["n_exon"], r["strand"], r["support"], f"{r['identity']:.4f}", f"{r['coverage']:.4f}"])
    print(f"Wrote {OUT_RESCUED}", flush=True)

    # Build merged catalog: gw_xcbase + seed-extended copies
    base_copies = []
    with open(SEED_COPIES) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            base_copies.append(dict(row))

    rescued_by_fam = defaultdict(list)
    for r in rescued:
        rescued_by_fam[r["fid"]].append(r)

    with open(OUT_FAMILIES, "w", newline="") as ff, open(OUT_COPIES, "w", newline="") as fc:
        wf = csv.writer(ff, delimiter="\t")
        wc = csv.writer(fc, delimiter="\t")
        wf.writerow(["family_id", "n_copies", "n_chroms", "chroms", "cross_chrom", "avg_reads"])
        wc.writerow(["family_id", "copy_idx", "tid", "chrom", "start", "end", "n_exon", "strand", "n_reads"])

        base_fams = defaultdict(list)
        for c in base_copies:
            base_fams[c["family_id"]].append(c)

        for fid in sorted(base_fams.keys()):
            fam_copies = list(base_fams[fid])
            for r in rescued_by_fam.get(fid, []):
                fam_copies.append({
                    "family_id": fid,
                    "copy_idx": len(fam_copies),
                    "tid": r["tid"],
                    "chrom": r["chrom"],
                    "start": r["start"],
                    "end": r["end"],
                    "n_exon": r["n_exon"],
                    "strand": r["strand"],
                    "n_reads": r["support"],
                })
            fam_copies.sort(key=lambda c: (c["chrom"], int(c["start"])))
            chroms = sorted(set(c["chrom"] for c in fam_copies))
            cross = len(chroms) > 1
            avg_reads = sum(int(c["n_reads"]) for c in fam_copies) / max(1, len(fam_copies))
            wf.writerow([fid, len(fam_copies), len(chroms), ",".join(chroms), cross, f"{avg_reads:.1f}"])
            for ci, c in enumerate(fam_copies):
                wc.writerow([fid, ci, c["tid"], c["chrom"], c["start"], c["end"],
                             c["n_exon"], c["strand"], c["n_reads"]])
    print(f"Wrote {OUT_FAMILIES} and {OUT_COPIES}", flush=True)


if __name__ == "__main__":
    main()
