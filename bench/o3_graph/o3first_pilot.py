#!/usr/bin/env python3
"""EXPLORATORY pilot of the order O1 -> O3 -> O2 on HG002 (2026-10-22; doc: docs/O3_FIRST_PILOT_2026-10-22.md).

O3 runs first, on each family's whole read net (`candidate_copies`, the opt-in pipeline stage), before any tied read is placed.
This file is the glue around it:
- `copies POOL GENOME OUT`: the O2-VC copy table (copy family chrom start end strand exons) -> the catalog format the stage reads
  (OUT.fam.copies.tsv; OUT.fam.copies.fa = each copy's spliced exon sum in transcription orientation, as gw_family_catalog writes it).
  Copies on chrM are dropped, as `o2vc.py run --exclude-chrM` drops them.
"""
import csv
import sys

import pysam

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def union_blocks(exons):
    """The O2-VC table lists every exon variant of a copy (793 of its 1,686 copies overlap); the catalog wants ascending, disjoint
    blocks, so overlapping or touching blocks are merged into the copy's exon union (the stage's own representative, R13)."""
    out = []
    for a, b in sorted(tuple(map(int, x.split("-"))) for x in exons.split(",") if x):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def cmd_copies(pool_p, genome_p, out):
    fa = pysam.FastaFile(genome_p)
    rows = [r for r in csv.DictReader(open(pool_p), delimiter="\t") if r["chrom"] != "chrM"]
    seen = set()
    with open(f"{out}.fam.copies.tsv", "w") as t, open(f"{out}.fam.copies.fa", "w") as s:
        t.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
        for r in rows:
            idx = int(r["copy"].rsplit(".", 1)[1])
            key = (r["family"], idx)
            if key in seen:
                sys.exit(f"duplicate copy index {key}")
            seen.add(key)
            blocks = union_blocks(r["exons"])
            seq = "".join(fa.fetch(r["chrom"], a, b) for a, b in blocks).upper()
            if r["strand"] == "-":
                seq = revcomp(seq)
            ex = ",".join(f"{a}-{b}" for a, b in blocks)
            t.write(f"{r['family']}\t{idx}\t{r['copy']}\t{r['chrom']}\t{r['start']}\t{r['end']}\t{len(blocks)}\t{r['strand']}\t0\t{ex}\n")
            s.write(f">{r['family']}|{idx}|{r['chrom']}:{r['start']}-{r['end']}|{r['strand']}|nexon={len(blocks)}\n{seq}\n")
    print(f"{len(rows)} copies, {len({r['family'] for r in rows})} families -> {out}.fam.copies.*")


MIN_OVERLAP = 100      # aligned bases a record must put on a copy's exon union to count as hitting it


def load_copies(tsv_p):
    copies = {}
    for r in csv.DictReader(open(tsv_p), delimiter="\t"):
        blocks = [tuple(map(int, b.split("-"))) for b in r["exons"].split(",") if b]
        copies[r["tid"]] = (r["family_id"], r["chrom"], blocks)
    return copies


def overlap(blocks, exons):
    return sum(max(0, min(b, d) - max(a, c)) for a, b in blocks for c, d in exons)


def cmd_silent(bam_p, tsv_p, fa_p, out_p):
    """The advisor's claim (2026-10-23): a copy with no reads of its own is still found through the multimapping reads of its
    expressed relatives. Per family: expressed copies (>= 3 reads with their primary on the copy) and silent copies (0); for each
    silent copy, the reads of the family's expressed copies with a secondary record on it, and its identity / coverage against
    the nearest expressed copy (minimap2 -c -x asm20 on the copies' spliced sequences; no hit = below asm20's reach)."""
    import collections
    import os
    import subprocess
    import tempfile
    copies = load_copies(tsv_p)
    bam = pysam.AlignmentFile(bam_p)
    prim = collections.defaultdict(set)          # copy -> reads with their primary record on it
    sec = collections.defaultdict(set)           # copy -> reads with a secondary record on it
    for cid, (fam, chrom, exons) in copies.items():
        for rec in bam.fetch(chrom, exons[0][0], exons[-1][1]):
            if rec.is_unmapped or rec.is_supplementary:
                continue
            if overlap(rec.get_blocks(), exons) < MIN_OVERLAP:
                continue
            (sec if rec.is_secondary else prim)[cid].add(rec.query_name)
    read_fam = collections.defaultdict(set)      # read -> families of its primary copy
    for cid, names in prim.items():
        for n in names:
            read_fam[n].add(copies[cid][0])
    by_fam = collections.defaultdict(list)
    for cid, (fam, _, _) in copies.items():
        by_fam[fam].append(cid)
    seqs = {}
    for h, s in ((rec.name, rec.sequence) for rec in pysam.FastxFile(fa_p)):
        f, i = h.split("|")[:2]
        seqs[(f, int(i))] = s
    seq_of = {cid: seqs[(fam, int(cid.rsplit(".", 1)[1]))] for cid, (fam, _, _) in copies.items()}
    rows = []
    tmp = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(out_p)))
    for fam, cids in sorted(by_fam.items()):
        expressed = [c for c in cids if len(prim[c]) >= 3]
        silent = [c for c in cids if not prim[c]]
        if not expressed or not silent:
            continue
        t, q = os.path.join(tmp, "t.fa"), os.path.join(tmp, "q.fa")
        with open(t, "w") as fh:
            fh.writelines(f">{c}\n{seq_of[c]}\n" for c in expressed)
        with open(q, "w") as fh:
            fh.writelines(f">{c}\n{seq_of[c]}\n" for c in silent)
        paf = subprocess.run(["minimap2", "-c", "-x", "asm20", "--secondary=no", t, q], capture_output=True, text=True).stdout
        best = {}
        for line in paf.splitlines():
            p = line.split("\t")
            qn, qlen, qs, qe, tn, tlen, ts, te, nm, bl = p[0], int(p[1]), int(p[2]), int(p[3]), p[5], int(p[6]), int(p[7]), int(p[8]), int(p[9]), int(p[10])
            if qn not in best or nm > best[qn][3]:
                best[qn] = (tn, nm / bl, (qe - qs) / qlen, nm, (te - ts) / tlen)
        fam_reads = set().union(*(prim[c] for c in expressed))
        for c in silent:
            hit = sec[c] & fam_reads
            tn, ident, qcov, _, tcov = best.get(c, ("", 0.0, 0.0, 0, 0.0))
            rows.append((fam, c, len(seq_of[c]), len(expressed), len(hit), tn, ident, qcov, tcov))
    with open(out_p, "w") as fh:
        fh.write("family\tsilent_copy\tlen\texpressed_copies\tfamily_reads_with_secondary\tnearest_expressed\tidentity\tcov_silent\tcov_expressed\n")
        for r in rows:
            fh.write("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
    print(f"{len(rows)} silent copies in {len({r[0] for r in rows})} families with an expressed copy -> {out_p}")


def cmd_strainy_inputs(bam_p, tsv_p, fa_p, families, out_dir):
    """Strainy arm (2026-10-23, advisor): per family, OUT/F.fa = its copies' spliced exon unions named by copy, the collapsed
    sequences Strainy phases; OUT/F.fq = every read with a record (primary or secondary) of >= MIN_OVERLAP bases on a copy,
    its sequence taken from its primary record."""
    import collections
    import os
    os.makedirs(out_dir, exist_ok=True)
    fams = set(families.split(","))
    copies = {c: v for c, v in load_copies(tsv_p).items() if v[0] in fams}
    bam = pysam.AlignmentFile(bam_p)
    fam_of = collections.defaultdict(set)
    for cid, (fam, chrom, exons) in copies.items():
        for rec in bam.fetch(chrom, exons[0][0], exons[-1][1]):
            if not rec.is_unmapped and not rec.is_supplementary and overlap(rec.get_blocks(), exons) >= MIN_OVERLAP:
                fam_of[rec.query_name].add(fam)
    for h, s in ((rec.name, rec.sequence) for rec in pysam.FastxFile(fa_p)):
        f, i = h.split("|")[:2]
        if f in fams:
            with open(os.path.join(out_dir, f"{f}.fa"), "a") as fh:
                fh.write(f">{f}.{i}\n{s}\n")
    out = {f: open(os.path.join(out_dir, f"{f}.fq"), "w") for f in fams}
    n = collections.Counter()
    for rec in pysam.AlignmentFile(bam_p).fetch(until_eof=True):
        if rec.is_secondary or rec.is_supplementary or rec.query_name not in fam_of:
            continue
        q = "".join(chr(x + 33) for x in rec.query_qualities) if rec.query_qualities is not None else "I" * rec.query_length
        for f in fam_of[rec.query_name]:
            out[f].write(f"@{rec.query_name}\n{rec.query_sequence}\n+\n{q}\n")
            n[f] += 1
    for fh in out.values():
        fh.close()
    print(f"{len(fams)} families, reads per family: min {min(n.values(), default=0)} median {sorted(n.values())[len(n) // 2] if n else 0} max {max(n.values(), default=0)}")


UNCLUSTERED = "1000000"   # Strainy's label for reads in no cluster (params.UNCLUSTERED_GROUP_N)
PAR = {"chrX": ((0, 2394410), (153925834, 154259566)), "chrY": ((0, 2458320), (62122809, 62460029))}   # hg002truth.PAR


def strainy_multiplicity(fam_dir, copy):
    """Haplotypes Strainy finds on one copy: the most clusters (each >= 2 reads, Strainy's unclustered group excluded) that
    cover one of Strainy's own SNP positions together. Counting at SNP positions keeps one haplotype that Strainy split at a
    low-heterozygosity stretch from counting twice (its pieces only meet where no SNP is). None = Strainy did not phase it."""
    import glob
    import os
    cl_p = glob.glob(os.path.join(fam_dir, "intermediate", "clusters", f"clusters_{copy}_*.csv"))
    if not cl_p or os.path.getsize(cl_p[0]) == 0:
        return None
    cluster = {r["ReadName"]: r["Cluster"] for r in csv.DictReader(open(cl_p[0])) if r["Cluster"] != UNCLUSTERED}
    vcf = os.path.join(fam_dir, "intermediate", "vcf", f"vcf_{copy}_filtered.vcf")
    snps = sorted({int(l.split("\t")[1]) - 1 for l in open(vcf) if not l.startswith("#")}) if os.path.exists(vcf) else []
    spans = []
    for rec in pysam.AlignmentFile(os.path.join(fam_dir, "alignment_phased.bam")).fetch(copy):
        if not rec.is_secondary and not rec.is_supplementary and rec.query_name in cluster:
            spans.append((rec.reference_start, rec.reference_end, cluster[rec.query_name]))
    if not snps:
        return min(1, len(set(c for *_, c in spans)))
    best = 0
    for x in snps:
        n = {}
        for s, e, c in spans:
            if s <= x < e:
                n[c] = n.get(c, 0) + 1
        best = max(best, sum(1 for v in n.values() if v >= 2))
    return best


def cmd_score(work):
    """Family-level table: Strainy (a copy with more concurrent haplotypes than its allele slots), candidate_copies (flagged / any
    candidate), and the 2026-10-22 spike's tier 1 and today's flag, on the spike's positives and its 73 others (PAR excluded)."""
    import collections
    import glob
    import json
    import os
    par = set(json.load(open("/mnt/linuxdisk/tmp/o3_graph/famident/par_families.json")))
    pos_a = ["S19SPDYE1", "S19RASA4B", "S19USP17L11", "S19POTEM", "S19BMS1P1"]
    pos_b = pos_a + ["S19NPIPB2"]
    fams = [f for b in ["batch_pos"] + [f"batch_o{i}" for i in range(6)] for f in open(f"{work}/{b}.txt").read().split(",")]
    others = [f for f in fams if f not in pos_b and f not in par]
    copies = load_copies(f"{work}/hg002.fam.copies.tsv")
    by_fam = collections.defaultdict(list)
    for cid, (fam, chrom, exons) in copies.items():
        by_fam[fam].append((cid, chrom, exons[0][0], exons[-1][1]))
    cand = collections.defaultdict(lambda: [0, 0])
    for p in glob.glob(f"{work}/*.cand.candidates.tsv"):
        for r in csv.DictReader(open(p), delimiter="\t"):
            cand[r["family"]][0] += 1
            cand[r["family"]][1] += int(r["flagged"])
    spike = {r["family"]: r for r in csv.DictReader(open(
        "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/e542c6ba-28cf-42e7-becf-fd99d54f830f/scratchpad/hg002.tsv"), delimiter="\t")}
    rows = {}
    with open(f"{work}/strainy_family.tsv", "w") as out:
        out.write("family\tstrainy_status\tcopies_phased\tflagged_copies\tmax_haplotypes\tcc_candidates\tcc_flagged\n")
        for f in fams:
            done = open(f"{work}/strainy_out/{f}.done").read().split()[-1] if os.path.exists(f"{work}/strainy_out/{f}.done") else "?"
            status = "ok" if done in ("0", "") else f"exit{done}"
            phased, flagged, mx = 0, [], 0
            if status == "ok":
                for cid, chrom, s, e in by_fam[f]:
                    m = strainy_multiplicity(f"{work}/strainy_out/{f}", cid)
                    if m is None:
                        continue
                    phased += 1
                    mx = max(mx, m)
                    cap = 1 if chrom in PAR and not any(s < pe and e > ps for ps, pe in PAR[chrom]) else 2
                    if m > cap:
                        flagged.append(f"{cid}:{m}/{cap}")
            rows[f] = dict(status=status, phased=phased, flagged=flagged, mx=mx)
            out.write(f"{f}\t{status}\t{phased}\t{','.join(flagged)}\t{mx}\t{cand[f][0]}\t{cand[f][1]}\n")
    arms = [("Strainy: haplotypes > allele slots", lambda f: bool(rows[f]["flagged"])),
            ("candidate_copies: flagged candidate", lambda f: cand[f][1] > 0),
            ("candidate_copies: any candidate", lambda f: cand[f][0] > 0),
            ("spike tier 1 (10-22)", lambda f: f in spike and int(spike[f]["tier1"]) > 0),
            ("today's flag (O2-VC)", lambda f: f in spike and int(spike[f]["flags"]) > 0)]
    print(f"{len(fams)} families: positives (a) {len(pos_a)}, (b) {len(pos_b)}; others {len(others)}")
    print(f"Strainy phased >= 1 copy in {sum(1 for f in fams if rows[f]['phased'])} families; not run to the end: "
          f"{[f for f in fams if rows[f]['status'] != 'ok']}")
    for name, fn in arms:
        a = [f for f in pos_a if fn(f)]
        b = [f for f in pos_b if fn(f)]
        o = [f for f in others if fn(f)]
        print(f"{name}: (a) {len(a)}/{len(pos_a)} {a} | (b) {len(b)}/{len(pos_b)} | others {len(o)}/{len(others)} {o}")
    print("Strainy per positive:", {f: (rows[f]["phased"], rows[f]["mx"], rows[f]["flagged"]) for f in pos_b})


def main(argv):
    if len(argv) == 4 and argv[0] == "copies":
        cmd_copies(*argv[1:])
    elif len(argv) == 2 and argv[0] == "score":
        cmd_score(argv[1])
    elif len(argv) == 6 and argv[0] == "strainy-inputs":
        cmd_strainy_inputs(*argv[1:])
    elif len(argv) == 5 and argv[0] == "silent":
        cmd_silent(*argv[1:])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
