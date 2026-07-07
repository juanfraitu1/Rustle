#!/usr/bin/env python3
"""psv_split_stringtie.py — PSV-aware splitting of StringTie transcripts in multi-copy families.

For each StringTie transcript that overlaps a known multi-copy family, fetch its
supporting reads and assign them to copies using PSVs in genomic coordinates. If a
transcript's reads belong to multiple copies, split it into copy-specific isoforms.

Run example (RABL2, catalog fam 70):
  PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/psv_split_stringtie.py --family 70
"""
import argparse
import collections
import csv
import os
import subprocess
import sys
import tempfile

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = "/home/juanfra/winloci_scratch"

REFINE_TSV = os.path.join(HERE, "family_rna_refine.tsv")
GENOME = os.path.join(SCRATCH, "GGO.fasta")
BAM = os.path.join(SCRATCH, "GGO.bam")
STRINGTIE_GTF = os.path.join(SCRATCH, "genome_st.gtf")
REFSEQ_GFF = os.path.join(SCRATCH, "GGO_tx.gff")
OUT_DIR = os.path.join(HERE, "psv_split_stringtie")

MIN_READS_PER_COPY = 3


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def load_family_members(family_id):
    members = []
    with open(REFINE_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if int(r["family_id"]) == family_id:
                members.append((r["member_dn"], r["chrom"], int(r["start"]), int(r["end"]), r["member_gene"]))
    return members


def parse_gff_attrs(field):
    attrs = {}
    for tok in field.split(";"):
        tok = tok.strip()
        if not tok:
            continue
        if "=" in tok:
            k, v = tok.split("=", 1)
            attrs[k] = v
    return attrs


def load_refseq_transcripts(chrom, start, end, gene_name=None):
    """Return list of (transcript_id, strand, exons, spliced_len, biotype) from RefSeq GFF.
       If gene_name is given, prefer transcripts whose gene attribute matches it."""
    if not os.path.exists(REFSEQ_GFF):
        return []
    genes = {}
    tx_parent = {}
    txs = {}
    with open(REFSEQ_GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] != chrom:
                continue
            s, e = int(p[3]), int(p[4])
            typ = p[2]
            strand = p[6]
            attrs = parse_gff_attrs(p[8])
            gid = attrs.get("ID", "")
            gene = attrs.get("gene", "")
            parent = attrs.get("Parent", "")
            if typ == "gene":
                genes[gid] = {"start": s, "end": e, "gene": gene}
            elif typ in ("mRNA", "transcript", "lnc_RNA", "misc_RNA"):
                tid = gid
                # only keep transcripts overlapping the region
                if e < start or s > end:
                    continue
                tx_parent[tid] = parent
                txs[tid] = {"strand": strand, "exons": [], "type": typ,
                            "start": s, "end": e, "gene": gene}
            elif typ == "exon" and parent in txs:
                txs[parent]["exons"].append((s, e))

    # compute gene for each transcript from parent gene
    for tid, t in txs.items():
        if not t["gene"] and tx_parent[tid] in genes:
            t["gene"] = genes[tx_parent[tid]]["gene"]

    out = []
    for tid, t in txs.items():
        if not t["exons"]:
            continue
        exons = sorted(t["exons"])
        spliced_len = sum(e - s + 1 for s, e in exons)
        # prefer mRNA, then transcript, then others
        biotype = t["type"]
        match = 1 if gene_name and t["gene"] == gene_name else 0
        out.append((tid, t["strand"], exons, spliced_len, biotype, match, t["start"], t["end"]))
    # sort: gene match first, then mRNA, then spliced length
    biotype_rank = {"mRNA": 0, "transcript": 1, "lnc_RNA": 2, "misc_RNA": 3}
    out.sort(key=lambda x: (-x[5], biotype_rank.get(x[4], 9), -x[3]))
    return [(tid, strand, exons, spliced_len, 0.0) for tid, strand, exons, spliced_len, *_ in out]


def load_stringtie_transcripts(chrom, start, end):
    """Return list of (transcript_id, strand, exons, spliced_len, cov) from StringTie GTF."""
    txs = {}
    with open(STRINGTIE_GTF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] != chrom:
                continue
            s, e = int(p[3]), int(p[4])
            if e < start or s > end:
                continue
            typ = p[2]
            strand = p[6]
            attrs = p[8]
            tid = None
            cov = 0.0
            for tok in attrs.split(";"):
                tok = tok.strip()
                if tok.startswith("transcript_id"):
                    tid = tok.split('"')[1]
                if tok.startswith("cov"):
                    try:
                        cov = float(tok.split('"')[1])
                    except (IndexError, ValueError):
                        pass
            if tid is None:
                continue
            if typ == "transcript":
                txs[tid] = {"strand": strand, "exons": [], "cov": cov}
            elif typ == "exon" and tid in txs:
                txs[tid]["exons"].append((s, e))
    out = []
    for tid, t in txs.items():
        if not t["exons"]:
            continue
        exons = sorted(t["exons"])
        spliced_len = sum(e - s + 1 for s, e in exons)
        out.append((tid, t["strand"], exons, spliced_len, t["cov"]))
    return out


def revcomp(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
            "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}
    return "".join(comp.get(b, b) for b in reversed(seq))


def fetch_spliced_sequence(genome, chrom, exons, strand="+"):
    """Return concatenated exon sequence in 5'->3' transcript orientation."""
    parts = []
    for s, e in sorted(exons):
        s2 = max(0, s - 1)
        parts.append(genome.fetch(chrom, s2, e).upper())
    seq = "".join(parts)
    if strand == "-":
        seq = revcomp(seq)
    return seq


def extract_copy_sequences(members, tmp_dir, use_spliced=True):
    """Write ref.fa (best spliced transcript) and copies.fa.
       Return (names, members, refseq, seq_info).
       seq_info[nm] = {"type": "spliced"|"genomic", "chrom", "exons" or (start,end)}"""
    genome = pysam.FastaFile(GENOME)
    contig_len = dict(zip(genome.references, genome.lengths))

    seq_info = {}
    spliced_seqs = {}

    for i, (dn, chrom, s, e, gene) in enumerate(members):
        nm = f"copy{i}"
        if use_spliced:
            # prefer RefSeq annotated transcript matching the gene symbol
            txs = load_refseq_transcripts(chrom, s, e, gene_name=gene if gene != "NA" else None)
            if not txs:
                txs = load_stringtie_transcripts(chrom, s, e)
            if txs:
                tid, strand, exons, spliced_len, cov = txs[0]
                # keep only exons that overlap the member region (with small margin)
                margin = 100
                region_s, region_e = s, e
                exons = [(es, ee) for es, ee in exons if ee >= region_s - margin and es <= region_e + margin]
                if not exons:
                    # fall through to genomic
                    pass
                else:
                    seq = fetch_spliced_sequence(genome, chrom, exons, strand)
                    spliced_seqs[nm] = seq
                    seq_info[nm] = {"type": "spliced", "chrom": chrom, "exons": exons,
                                    "strand": strand, "tid": tid, "cov": cov}
                    print(f"    {nm}: spliced {tid} {strand} {len(exons)} exons, {len(seq)} bp, cov={cov:.2f}")
                    continue
        # fallback to genomic interval
        s2, e2 = max(0, s), min(contig_len.get(chrom, e), e)
        seq = genome.fetch(chrom, s2, e2).upper()
        spliced_seqs[nm] = seq
        seq_info[nm] = {"type": "genomic", "chrom": chrom, "start": s2, "end": e2}
        print(f"    {nm}: genomic {chrom}:{s2}-{e2} ({len(seq)} bp)")

    genome.close()

    # reference = longest spliced sequence with an annotated gene symbol
    # (avoids choosing a spurious NA outlier that dominates by length)
    annotated = [n for n in spliced_seqs if members[int(n[4:])][4] != "NA"]
    if not annotated:
        annotated = list(spliced_seqs.keys())
    names_sorted = sorted(annotated, key=lambda n: -len(spliced_seqs[n]))
    names = names_sorted + [n for n in spliced_seqs if n not in names_sorted]
    refseq = spliced_seqs[names[0]]
    with open(os.path.join(tmp_dir, "ref.fa"), "w") as o:
        o.write(f">{names[0]}\n{refseq}\n")
    with open(os.path.join(tmp_dir, "copies.fa"), "w") as o:
        for nm in names[1:]:
            o.write(f">{nm}\n{spliced_seqs[nm]}\n")

    return names, members, refseq, seq_info


def align_copies_minimap(tmp_dir):
    sh(f"samtools faidx {tmp_dir}/ref.fa")
    sh(f"minimap2 -ax asm20 {tmp_dir}/ref.fa {tmp_dir}/copies.fa 2>/dev/null "
       f"| samtools sort -o {tmp_dir}/copies.bam - && samtools index {tmp_dir}/copies.bam")


def align_copies_mafft(tmp_dir):
    """Run mafft on ref.fa + copies.fa, produce aligned.fa."""
    sh(f"cat {tmp_dir}/ref.fa {tmp_dir}/copies.fa > {tmp_dir}/all.fa")
    sh(f"mafft --auto --quiet {tmp_dir}/all.fa > {tmp_dir}/aligned.fa 2>/dev/null")


def read_fasta(path):
    """Return {name: seq} from FASTA."""
    out = {}
    name = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            out[name] = "".join(seq)
    return out


def call_psvs_mafft(tmp_dir, names, refseq):
    """Call PSV columns from mafft MSA.
       Returns (psvs, n_aligned) where psvs = [(ref_pos, {name: allele})]."""
    aln = read_fasta(f"{tmp_dir}/aligned.fa")
    # ensure all names present
    for nm in names:
        if nm not in aln:
            return [], 0

    # build ungapped position maps: for each column, what is the 0-based ungapped position in each seq
    ungapped = {nm: -1 for nm in names}
    maps = {nm: [] for nm in names}  # maps[nm][col] = ungapped pos or None
    seqs = [aln[nm] for nm in names]
    n_cols = len(seqs[0])
    for col in range(n_cols):
        for nm in names:
            base = aln[nm][col]
            if base != "-":
                ungapped[nm] += 1
                maps[nm].append(ungapped[nm])
            else:
                maps[nm].append(None)

    psvs = []
    aligned_count = 0
    for col in range(n_cols):
        bases = {}
        all_defined = True
        for nm in names:
            pos = maps[nm][col]
            if pos is None:
                all_defined = False
                break
            bases[nm] = aln[nm][col].upper()
        if not all_defined:
            continue
        aligned_count += 1
        if len(set(bases.values())) >= 2:
            # ref_pos is the ungapped position in reference (names[0])
            ref_pos = maps[names[0]][col]
            psvs.append((ref_pos, bases))
    return psvs, aligned_count


def build_ref_to_query_map_mafft(tmp_dir, names):
    """Build {copy_name: {ref_pos: query_pos}} from mafft MSA."""
    aln = read_fasta(f"{tmp_dir}/aligned.fa")
    ref_nm = names[0]
    ungapped = {nm: -1 for nm in names}
    ref_map = []
    query_maps = {nm: {} for nm in names[1:]}
    n_cols = len(aln[ref_nm])
    for col in range(n_cols):
        for nm in names:
            if aln[nm][col] != "-":
                ungapped[nm] += 1
        ref_pos = ungapped[ref_nm]
        for nm in names[1:]:
            if aln[nm][col] != "-":
                query_maps[nm][ref_pos] = ungapped[nm]
    return query_maps


def build_ref_to_query_map(tmp_dir, names):
    """For each copy j>0, build dict {ref_pos: query_pos} from aligned pairs."""
    maps = {nm: {} for nm in names[1:]}
    bam = pysam.AlignmentFile(f"{tmp_dir}/copies.bam", "rb")
    for aln in bam.fetch(names[0]):
        if aln.is_unmapped or aln.is_supplementary or aln.is_secondary:
            continue
        nm = aln.query_name
        for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
            maps[nm][rpos] = qpos
    bam.close()
    return maps


def copy_alleles_from_bam(bampath, ref_name):
    """{ref_pos: {copy_name: base}} from aligned pairs."""
    out = collections.defaultdict(dict)
    bam = pysam.AlignmentFile(bampath, "rb")
    for aln in bam.fetch(ref_name):
        if aln.is_unmapped or aln.is_supplementary or aln.is_secondary:
            continue
        qs = aln.query_sequence
        if qs is None:
            continue
        for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
            out[rpos][aln.query_name] = qs[qpos]
    bam.close()
    return out


def call_psvs(tmp_dir, names, refseq):
    """Return list of (ref_pos, {copy_name: allele}) PSV columns."""
    copy_cols = copy_alleles_from_bam(f"{tmp_dir}/copies.bam", names[0])
    aligned = {names[0]}
    for d in copy_cols.values():
        aligned.update(d)

    psvs = []
    for p in sorted(copy_cols):
        cc = copy_cols[p]
        if not all(nm in cc for nm in names[1:]):
            continue
        hap = {names[0]: refseq[p]}
        hap.update(cc)
        if len(set(hap.values())) >= 2:
            psvs.append((p, hap))
    return psvs, len(aligned)


def spliced_to_genomic(qpos, exons, strand="+"):
    """Map a 0-based position in the spliced sequence back to genomic coordinate."""
    if strand == "-":
        # position 0 is the last base of the highest-coordinate exon
        total = sum(e - s + 1 for s, e in exons)
        rev_pos = total - 1 - qpos
        qpos = rev_pos
    offset = 0
    for s, e in sorted(exons):
        exon_len = e - s + 1
        if qpos < offset + exon_len:
            return s + (qpos - offset)
        offset += exon_len
    return None


def psvs_to_genomic(psvs, names, seq_info, ref_to_query):
    """Convert PSVs from reference-copy coordinates to genomic coordinates per copy.
       Returns {copy_name: [(chrom, gpos, allele), ...]}."""
    out = collections.defaultdict(list)
    ref_info = seq_info[names[0]]
    for p, hap in psvs:
        # reference copy
        if ref_info["type"] == "spliced":
            gpos = spliced_to_genomic(p, ref_info["exons"], ref_info.get("strand", "+"))
            chrom = ref_info["chrom"]
        else:
            gpos = ref_info["start"] + p
            chrom = ref_info["chrom"]
        if gpos is not None:
            out[names[0]].append((chrom, gpos, hap[names[0]]))
        # other copies
        for nm in names[1:]:
            qpos = ref_to_query[nm].get(p)
            if qpos is None:
                continue
            info = seq_info[nm]
            if info["type"] == "spliced":
                gpos = spliced_to_genomic(qpos, info["exons"], info.get("strand", "+"))
                chrom = info["chrom"]
            else:
                gpos = info["start"] + qpos
                chrom = info["chrom"]
            if gpos is not None:
                out[nm].append((chrom, gpos, hap[nm]))
    return out


def parse_stringtie_region(chrom, region_start, region_end):
    """Return transcript_id -> {chrom,start,end,strand,exons}."""
    transcripts = {}
    with open(STRINGTIE_GTF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] != chrom:
                continue
            start, end = int(p[3]), int(p[4])
            if end < region_start or start > region_end:
                continue
            typ = p[2]
            strand = p[6]
            tid = None
            for tok in p[8].split(";"):
                tok = tok.strip()
                if tok.startswith("transcript_id"):
                    tid = tok.split('"')[1]
                    break
            if tid is None:
                continue
            if typ == "transcript":
                transcripts[tid] = {"chrom": chrom, "start": start, "end": end,
                                    "strand": strand, "exons": []}
            elif typ == "exon" and tid in transcripts:
                transcripts[tid]["exons"].append((start, end))
    return transcripts


def fetch_reads_over_locus(bam, chrom, start, end):
    """Yield primary non-supplementary alignments overlapping [start,end]."""
    for r in bam.fetch(chrom, start, end):
        if r.is_unmapped or r.is_supplementary or r.is_secondary:
            continue
        if r.reference_end is None or r.reference_end < start or r.reference_start > end:
            continue
        yield r


def read_allele_at(read, chrom, gpos):
    """Return base at genomic position (chrom,gpos) or None."""
    if read.reference_name != chrom:
        return None
    qs = read.query_sequence
    if qs is None:
        return None
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if rpos == gpos:
            return qs[qpos]
    return None


def filter_psv_sites(copy_psv_sites, chrom, start, end):
    """Return copy_psv_sites restricted to positions within [start,end] on chrom.
       All copy names present in the input are preserved (empty list if no sites)."""
    out = {}
    for nm, sites in copy_psv_sites.items():
        out[nm] = [(c, g, a) for c, g, a in sites if c == chrom and start <= g <= end]
    return out


def assign_read_to_copy(read, copy_psv_sites, names):
    """Assign a read to the copy whose PSV alleles it matches best.
       copy_psv_sites = {copy_name: [(chrom, gpos, allele), ...]}
       Returns (copy_name, n_covered, n_mismatches) or None.
    """
    if read.reference_name is None:
        return None
    qs = read.query_sequence
    if qs is None:
        return None
    read_alleles = {rpos: qs[qpos] for qpos, rpos in read.get_aligned_pairs(matches_only=True)}

    best = []
    best_mm = None
    for nm in names:
        sites = copy_psv_sites.get(nm, [])
        covered = 0
        mm = 0
        for chrom, gpos, allele in sites:
            if chrom != read.reference_name:
                continue
            obs = read_alleles.get(gpos)
            if obs is None:
                continue
            covered += 1
            if obs != allele:
                mm += 1
        if covered == 0:
            continue
        if best_mm is None or mm < best_mm:
            best_mm = mm
            best = [(nm, covered, mm)]
        elif mm == best_mm:
            best.append((nm, covered, mm))

    if len(best) != 1:
        return None
    return best[0]


def split_transcripts(family_id, members, names, copy_psv_sites, out_gtf, out_tsv):
    # Per-chromosome aggregated family span
    by_chrom = collections.defaultdict(lambda: [None, None])
    for _, chrom, s, e, _ in members:
        if by_chrom[chrom][0] is None or s < by_chrom[chrom][0]:
            by_chrom[chrom][0] = s
        if by_chrom[chrom][1] is None or e > by_chrom[chrom][1]:
            by_chrom[chrom][1] = e

    bam = pysam.AlignmentFile(BAM, "rb")
    total_tx = 0
    split_tx = 0
    rows = []

    with open(out_gtf, "w") as o_gtf:
        o_gtf.write("##gff-version 2\n")
        o_gtf.write(f"# PSV-split StringTie transcripts for family {family_id}\n")

        for chrom, (cs, ce) in by_chrom.items():
            transcripts = parse_stringtie_region(chrom, cs - 5000, ce + 5000)
            tx_psv_sites = filter_psv_sites(copy_psv_sites, chrom, cs - 5000, ce + 5000)
            for tid, t in transcripts.items():
                total_tx += 1
                reads = list(fetch_reads_over_locus(bam, chrom, t["start"], t["end"]))
                # only consider PSV sites that fall inside this transcript
                local_psv = filter_psv_sites(tx_psv_sites, chrom, t["start"], t["end"])
                assignments = collections.Counter()
                covered = 0
                for r in reads:
                    res = assign_read_to_copy(r, local_psv, names)
                    if res is not None:
                        copy, n_psv, _ = res
                        assignments[copy] += 1
                        covered += 1

                copies_supported = [c for c, n in assignments.items() if n >= MIN_READS_PER_COPY]
                if len(copies_supported) >= 2:
                    split_tx += 1
                    gene_base = f"PSVSPLIT.{family_id}.{tid}"
                    for copy in copies_supported:
                        new_tid = f"{gene_base}.{copy}"
                        o_gtf.write(f"{chrom}\tPSVSPLIT\ttranscript\t{t['start']}\t{t['end']}\t1000\t{t['strand']}\t.\t"
                                    f"gene_id \"{gene_base}\"; transcript_id \"{new_tid}\"; copy \"{copy}\"; "
                                    f"reads \"{assignments[copy]}\";\n")
                        for j, (s, e) in enumerate(sorted(t["exons"]), 1):
                            o_gtf.write(f"{chrom}\tPSVSPLIT\texon\t{s}\t{e}\t1000\t{t['strand']}\t.\t"
                                        f"gene_id \"{gene_base}\"; transcript_id \"{new_tid}\"; exon_number \"{j}\";\n")

                rows.append((tid, chrom, t["start"], t["end"], len(reads), covered,
                             len(copies_supported), ",".join(copies_supported)))

    bam.close()

    with open(out_tsv, "w") as o:
        o.write("transcript_id\tchrom\tstart\tend\tn_reads\tn_psv_assigned\tn_copies_supported\tcopies\n")
        for row in rows:
            o.write("\t".join(str(x) for x in row) + "\n")

    return total_tx, split_tx


def main():
    parser = argparse.ArgumentParser(description="PSV-aware split StringTie transcripts in multi-copy families")
    parser.add_argument("--family", type=int, required=True, help="catalog family_id")
    parser.add_argument("--keep-tmp", action="store_true", help="keep temporary directory for debugging")
    parser.add_argument("--aligner", choices=["minimap2", "mafft"], default="mafft",
                        help="alignment engine for copy comparison (default: mafft)")
    args = parser.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    members = load_family_members(args.family)
    if not members:
        print(f"[!] family {args.family} not found in {REFINE_TSV}", file=sys.stderr)
        sys.exit(1)

    print(f"[+] family {args.family}: {len(members)} members")
    for i, (dn, c, s, e, g) in enumerate(members):
        print(f"    member {i}: {dn} {c}:{s}-{e} ({e-s} bp) gene={g}")

    tmp_factory = tempfile.TemporaryDirectory(dir=SCRATCH, prefix="psv_split_")
    if args.keep_tmp:
        tmp = os.path.join(SCRATCH, f"psv_split_fam{args.family}_{os.urandom(4).hex()}")
        os.makedirs(tmp, exist_ok=True)
    else:
        tmp = tmp_factory.__enter__()

    try:
        names, members, refseq, seq_info = extract_copy_sequences(members, tmp, use_spliced=True)
        print(f"[+] extracted {len(names)} copies; reference = {names[0]} ({len(refseq)} bp)")

        if args.aligner == "minimap2":
            align_copies_minimap(tmp)
            psvs, n_aligned = call_psvs(tmp, names, refseq)
            ref_to_query = build_ref_to_query_map(tmp, names)
        else:
            align_copies_mafft(tmp)
            psvs, n_aligned = call_psvs_mafft(tmp, names, refseq)
            ref_to_query = build_ref_to_query_map_mafft(tmp, names)
        print(f"[+] {len(psvs)} PSV columns from {n_aligned} aligned positions/columns")

        copy_psv_sites = psvs_to_genomic(psvs, names, seq_info, ref_to_query)
        print(f"[+] mapped PSVs to genomic coordinates for {len(copy_psv_sites)} copies")

        out_gtf = os.path.join(OUT_DIR, f"fam{args.family}_psv_split.gtf")
        out_tsv = os.path.join(OUT_DIR, f"fam{args.family}_psv_split.tsv")
        total, split = split_transcripts(args.family, members, names, copy_psv_sites, out_gtf, out_tsv)
    finally:
        if not args.keep_tmp:
            tmp_factory.__exit__(None, None, None)

    print(f"[+] StringTie transcripts in family region: {total}")
    print(f"[+] transcripts split into copy-specific isoforms: {split}")
    print(f"[+] wrote {out_gtf}")
    print(f"[+] wrote {out_tsv}")


if __name__ == "__main__":
    main()
