#!/usr/bin/env python3
"""Non-coding-aware promotion of reference-absent copies: promote credible collapses on a
COLLAPSE-QUALITY bar (divergence band + genome coverage + own-locus + balanced co-segregation)
instead of a protein-BLASTx hit. Protein/ORF are LABELS, not gates. Purely additive: re-scores the
734 pre-ORF-gate consensuses already in gw_promoted/cons.fa and writes gw_noncoding_copies.json.
Every record is a FLAGGED reference-divergent candidate (copy_vs_allele='candidate-DNA-needed');
copy-vs-allele is not resolvable from RNA and needs DNA parCN. See
docs/superpowers/specs/2026-07-18-noncoding-promotion-design.md."""
import os
import sys
import json
import subprocess
from collections import Counter, defaultdict

sys.path.insert(0, "bench")
from copy_vs_allele_structural import load_sedef, sedef_partners  # stdlib-only module

# ---- paths (overridable in main via argparse) ----
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
CAT = "/home/juanfra/winloci_scratch/refabsent"
OUT = f"{CAT}/gw_promoted"
GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff"
SEDEF = next((p for p in ["/mnt/c/Users/jfris/Desktop/final.bed",
                          "/mnt/c/Users/jfris/Desktop/Rustle/final.bed"] if os.path.exists(p)), None)

# ---- the collapse-quality bar (single-edit recall/precision knob) ----
GID_LO, GID_HI = 0.60, 0.97   # genome_id band: <LO = repeat/chimera artifact, >=HI = ~reference allele
COV_MIN = 0.90                # consensus must map near full-length (not a chimeric fragment)
MIN_COLS = 8                  # co-segregation breadth (alt columns)
AF_LO, AF_HI = 0.15, 0.60     # overall alt-read fraction band (balanced collapse)
MIN_ALT = 5                   # minimum alt reads
MIN_DEPTH = 8                 # minimum primary reads (matches hidden_copy_scan.MIN_DEPTH)

# ---- honesty rail (set unconditionally on every promoted record) ----
CANDIDATE = "candidate-DNA-needed"
STATUS = "flagged-reference-divergent-candidate"


def classify_noncoding(rec):
    """Pure collapse-quality gate. rec keys: genome_id, genome_cov, own_locus (bool),
    alt_cols, alt_read_fraction, alt_reads, n_primary.
    Returns (promote: bool, call: str, reason: str). Only call=='noncoding-candidate' promotes."""
    gid = rec["genome_id"]
    cov = rec["genome_cov"]
    if gid >= GID_HI:
        return (False, "~REF", f"genome_id {gid:.3f} >= {GID_HI} (essentially the reference allele)")
    if gid < GID_LO:
        return (False, "artifact", f"genome_id {gid:.3f} < {GID_LO} (repeat/chimera, not a divergent copy)")
    if cov < COV_MIN:
        return (False, "chimera", f"genome_cov {cov:.2f} < {COV_MIN} (partial/chimeric consensus)")
    if not rec["own_locus"]:
        return (False, "not-own-locus", "best genome hit is not the candidate's own locus")
    if rec["alt_cols"] < MIN_COLS:
        return (False, "thin", f"alt_cols {rec['alt_cols']} < {MIN_COLS} (too few co-segregating columns)")
    if not (AF_LO <= rec["alt_read_fraction"] <= AF_HI):
        return (False, "thin", f"alt_read_fraction {rec['alt_read_fraction']:.3f} outside "
                               f"[{AF_LO},{AF_HI}] (unbalanced)")
    if rec["alt_reads"] < MIN_ALT or rec["n_primary"] < MIN_DEPTH:
        return (False, "thin", f"depth alt_reads={rec['alt_reads']} (min {MIN_ALT}) / "
                               f"n_primary={rec['n_primary']} (min {MIN_DEPTH})")
    return (True, "noncoding-candidate",
            f"divergent ({100*(1-gid):.1f}%) full-length (cov {cov:.2f}) balanced collapse "
            f"({rec['alt_cols']} cols, af {rec['alt_read_fraction']:.2f}) at own locus")


# ---- ORF label (copied from promote_genomewide to keep this module import-light) ----
_COMP = str.maketrans("ACGTN", "TGCAN")
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
from itertools import product as _product
_CODON = {"".join(c): _AAS[i] for i, c in enumerate(_product("TCAG", repeat=3))}


def orf_aa(seq):
    """Longest ATG..stop ORF length in aa, scanning all 3 frames on both strands."""
    best = 0
    for strand in (seq, seq.translate(_COMP)[::-1]):
        for f in range(3):
            prot = "".join(_CODON.get(strand[i:i+3], "X") for i in range(f, len(strand)-2, 3))
            for p in prot.split("*"):
                k = p.find("M")
                if k >= 0:
                    best = max(best, len(p) - k)
    return best


def coding_potential(orf, cons_len):
    """Label (never a gate): 'noncoding' if the ORF spans <30% of the transcript."""
    frac = (orf * 3) / cons_len if cons_len else 0.0
    return "noncoding" if frac < 0.30 else "coding-capable/novel-protein"


def load_gff_from_lines(lines):
    """Build the biotype index from an iterable of GFF lines (test seam + file loader core)."""
    idx = defaultdict(list)
    for ln in lines:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        bt = "unknown"
        for kv in f[8].split(";"):
            if kv.startswith("gene_biotype="):
                bt = kv.split("=", 1)[1]
        try:
            s, e = int(f[3]), int(f[4])
        except ValueError:
            continue
        idx[f[0]].append((s, e, bt))
    for c in idx:
        idx[c].sort()
    return idx


def load_gff_genes(path):
    """chrom -> sorted [(start, end, biotype)] for gene features; {} if path missing/unreadable."""
    if not path or not os.path.exists(path):
        return {}
    with open(path) as fh:
        return load_gff_from_lines(fh)


def biotype_of(idx, chrom, start, end):
    """Annotation LABEL of a gene overlapping [start,end]; never a gate.
    'unknown' if no annotation loaded, 'intergenic' if loaded but no overlap."""
    if not idx:
        return "unknown"
    ivs = idx.get(chrom, [])
    hits = [bt for (s, e, bt) in ivs if s <= end and e >= start]
    if not hits:
        return "intergenic"
    if any(b == "lncRNA" for b in hits):
        return "lncRNA"
    if any(b == "protein_coding" for b in hits):
        return "protein_coding-gene-body"
    if any(b == "pseudogene" for b in hits):
        return "pseudogene"
    return hits[0]


def read_fasta(path):
    """header (first whitespace-delimited token after '>') -> concatenated sequence."""
    seqs, cur, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(buf)
            cur = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


_TSV_COLS = ["cid", "chrom", "start", "end", "track", "genome_hit", "genome_id", "genome_cov",
             "divergence", "alt_cols", "alt_read_fraction", "alt_reads", "n_primary", "orf_aa",
             "coding_potential", "biotype", "sedef_partner", "protein", "copy_vs_allele",
             "status", "reason"]


def write_tsv(records, path):
    with open(path, "w") as fh:
        fh.write("\t".join(_TSV_COLS) + "\n")
        for r in records:
            fh.write("\t".join(str(r.get(c, "")) for c in _TSV_COLS) + "\n")


def build_record(cid, flag, hit, cons_len, orf, biotype, sedef_partner, protein, reason):
    """Assemble a promoted non-coding record. copy_vs_allele and status are set UNCONDITIONALLY."""
    tgt, pos, gid, gcov = hit
    return dict(
        cid=cid, chrom=flag["chrom"], start=flag["start"], end=flag["end"], track="noncoding",
        genome_hit=f"{tgt}:{pos}", genome_id=round(gid, 3), genome_cov=round(gcov, 2),
        divergence=round(100 * (1 - gid), 1),
        alt_cols=flag["n_alt_positions"], alt_read_fraction=round(flag["alt_read_fraction"], 3),
        alt_reads=flag["n_alt_reads"], n_primary=flag["n_primary_reads"],
        orf_aa=orf, coding_potential=coding_potential(orf, cons_len),
        biotype=biotype, sedef_partner=sedef_partner, protein=protein,
        copy_vs_allele=CANDIDATE, status=STATUS, reason=reason)


def best_hits_from_paf(paf_text):
    """cid -> (target, pos, genome_id, genome_cov, mapq) for the best hit by identity*coverage."""
    best = {}
    for line in paf_text.splitlines():
        c = line.split("\t")
        if len(c) < 12:
            continue
        try:
            cid, qlen, tname, tstart = c[0], int(c[1]), c[5], int(c[7])
            nm, al, mapq = int(c[9]), int(c[10]), int(c[11])
        except (ValueError, IndexError):
            continue
        gid = nm / al if al else 0.0
        gcov = al / qlen if qlen else 0.0
        sc = gid * gcov
        if cid not in best or sc > best[cid][5]:
            best[cid] = (tname, tstart, gid, gcov, mapq, sc)
    return {cid: (t, p, gid, gcov, mapq) for cid, (t, p, gid, gcov, mapq, _sc) in best.items()}


def promote_all(cons, flags, hits, gff, sed, carried, exclude=frozenset()):
    """Pure orchestration: classify every consensus, assemble promoted records, tally calls.
    `exclude` = cids already promoted to the coding catalog; they are skipped (tallied
    'already-coding'), so the non-coding track surfaces only what the protein gate MISSED
    rather than duplicating the coding catalog."""
    promoted, tally = [], Counter()
    for cid, seq in cons.items():
        if cid in exclude:
            tally["already-coding"] += 1
            continue
        flag = flags.get(cid)
        if flag is None:
            tally["no-flag"] += 1
            continue
        h = hits.get(cid)
        if h is None:
            tally["no-hit"] += 1
            continue
        tgt, pos, gid, gcov, _mapq = h
        own = (tgt == flag["chrom"]) and abs(pos - flag["start"]) < 200000
        rec = dict(genome_id=gid, genome_cov=gcov, own_locus=own,
                   alt_cols=flag["n_alt_positions"], alt_read_fraction=flag["alt_read_fraction"],
                   alt_reads=flag["n_alt_reads"], n_primary=flag["n_primary_reads"])
        promote, call, reason = classify_noncoding(rec)
        tally[call] += 1
        if not promote:
            continue
        parts = sedef_partners(sed, flag["chrom"], flag["start"], flag["end"]) if sed else []
        sp = f"{parts[0][0]}:{parts[0][1]}" if parts else None
        promoted.append(build_record(
            cid, flag, (tgt, pos, gid, gcov), len(seq), orf_aa(seq),
            biotype_of(gff, flag["chrom"], flag["start"], flag["end"]),
            sp, carried.get(cid, "not-tested"), reason))
    promoted.sort(key=lambda r: r["genome_id"])
    return promoted, tally


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Non-coding-aware reference-absent promotion.")
    ap.add_argument("--cons", default=f"{OUT}/cons.fa")
    ap.add_argument("--out", default=f"{OUT}/gw_noncoding_copies.json")
    ap.add_argument("--threads", default="6")
    a = ap.parse_args()

    # union-load both flag files so every cons.fa cid finds its flag record
    flags = {}
    for fn in [f"{CAT}/genomewide_flags.json", f"{CAT}/genomewide_flags_new.json"]:
        if os.path.exists(fn):
            for f in json.load(open(fn)):
                flags[f"{f['chrom']}_{f['start']}"] = f
    cons = read_fasta(a.cons)
    gff = load_gff_genes(GFF)
    sed = load_sedef(SEDEF) if SEDEF and os.path.exists(SEDEF) else {}
    carried = {}
    rap = f"{OUT}/gw_reference_absent_copies.json"
    if os.path.exists(rap):
        carried = {r["cid"]: r.get("protein") for r in json.load(open(rap))}
    exclude = set()
    disc = f"{OUT}/gw_discriminated.json"
    if os.path.exists(disc):
        exclude = {r["cid"] for r in json.load(open(disc))}   # already-promoted coding copies

    paf = subprocess.run([MM2, "-cx", "splice:hq", "--eqx", "-N1", "-t", a.threads, FASTA, a.cons],
                         capture_output=True, text=True, timeout=3600).stdout
    hits = best_hits_from_paf(paf)

    promoted, tally = promote_all(cons, flags, hits, gff, sed, carried, exclude)
    json.dump(promoted, open(a.out, "w"), indent=1)
    write_tsv(promoted, a.out[:-5] + ".tsv" if a.out.endswith(".json") else a.out + ".tsv")

    bt = Counter(r["biotype"] for r in promoted)
    print(f"non-coding promotion: {len(promoted)} flagged reference-divergent candidates "
          f"(track=noncoding; all copy_vs_allele={CANDIDATE})")
    print(f"  by biotype: {dict(bt)}")
    print("  call tally (incl. exclusions, logged not dropped):")
    for k, v in tally.most_common():
        print(f"    {k:18s} {v}")
    print(f"  wrote {a.out} (+ .tsv)")


if __name__ == "__main__":
    main()
