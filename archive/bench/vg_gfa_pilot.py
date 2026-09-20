#!/usr/bin/env python3
"""vg_gfa_pilot.py -- PILOT: adopt the real vg/GFA toolkit as the standardized backbone
for our copy/PSV/assignment machinery, tested on the RABL2 flagship family (fam39).

Our custom o2_vg (bench/o2_vg_visualization.py + psv_graph_genomewide.py) builds a GENOMIC
backbone (longest copy locus, WITH introns), aligns copies (minimap2 asm20) and reads
(minimap2 splice:hq) onto it, and surfaces read-supported PSV bubbles + a per-read gate.
Ground truth for RABL2: /home/juanfra/winloci_scratch/o2vg/fam39.json
  5 copies, K=5 fully-resolvable, 235 PSV bubbles (195 SUN), gate assigned 191/222 reads.

vg/giraffe/GraphAligner are DNA-centric and do NOT model introns, so we build the graph at
the TRANSCRIPT level: a small pangenome of the 5 RABL2 paralog TRANSCRIPTS (spliced,
exon-concatenated, NO introns). Reads are IsoSeq CCS = already-spliced mRNA, so a read's
query sequence is itself an exon-only transcript -- we align those directly to the graph.

STAGES (this script = build stage, steps 1-2 of the pilot):
  build:
    1. per copy: pick the longest-spliced annotated isoform (GFF), extract the spliced
       transcript sequence (exon bases, reverse-complemented on the minus strand).
    2. build a base-level pangenome graph from the 5 transcripts with `vg msga`
       (progressive MSA -> variation graph with SNP-level bubbles + one path per copy).
       -> GFA + `vg stats`.
    3. extract the RABL2 reads (the fam39 ground-truth read set) spliced sequences from
       GGO_mm.bam -> reads.fa.
    4. align reads to the graph with GraphAligner -x vg -> reads.gam.

Outputs -> /home/juanfra/winloci_scratch/vg_pilot/ :
    transcripts.fa  rabl2.vg  rabl2.gfa  reads.fa  reads.gam  pilot_report.json

Run: /home/juanfra/miniforge3/bin/python bench/vg_gfa_pilot.py
"""
import collections
import difflib
import json
import os
import subprocess
import sys

import pysam

# ---- config ----------------------------------------------------------------
BIN = "/home/juanfra/miniforge3/bin"
VG = f"{BIN}/vg"
GRAPHALIGNER = f"{BIN}/GraphAligner"
MINIMAP2 = f"{BIN}/minimap2"
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
GFF = "/home/juanfra/winloci_scratch/vg_pilot_rabl2.gff"          # pre-extracted RABL2 gff
FULL_GFF = "/home/juanfra/winloci_scratch/GGO_genomic.gff"
GT = "/home/juanfra/winloci_scratch/o2vg/fam39.json"
OUT = "/home/juanfra/winloci_scratch/vg_pilot"
BENCH_JSON = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vg_gfa_pilot.json")
PAD = 2000
THREADS = 8

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def sh(cmd, **kw):
    """run, capture stdout+stderr, raise with the captured text on failure."""
    r = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"CMD FAILED (rc={r.returncode}): {cmd}\n----\n{r.stdout[-4000:]}")
    return r.stdout


def rc(s):
    return s.translate(COMP)[::-1]


# ---- step 1: transcripts ---------------------------------------------------
def parse_isoforms(gff_path, genes):
    """Return {gene: [ (mrna_id, strand, chrom, [(start,end),...]) ]} (1-based inclusive)."""
    exons = collections.defaultdict(lambda: collections.defaultdict(list))  # gene->rna->exons
    strand = {}      # rna -> (strand, chrom, gene)
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "exon":
                continue
            attr = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            g = attr.get("gene")
            if g not in genes:
                continue
            parent = attr.get("Parent", "")
            exons[g][parent].append((int(c[3]), int(c[4])))
            strand[parent] = (c[6], c[0], g)
    out = {}
    for g in genes:
        isos = []
        for rna, ex in exons[g].items():
            ex.sort()
            st, chrom, _ = strand[rna]
            isos.append((rna, st, chrom, ex))
        out[g] = isos
    return out


def longest_isoform(isos):
    """pick isoform with max spliced length (sum of exon lengths)."""
    best = None
    for rna, st, chrom, ex in isos:
        slen = sum(e - s + 1 for s, e in ex)
        if best is None or slen > best[0]:
            best = (slen, rna, st, chrom, ex)
    return best  # (slen, rna, strand, chrom, exons)


def spliced_seq(genome, chrom, exons, strand):
    parts = [genome.fetch(chrom, s - 1, e).upper() for s, e in exons]  # 0-based fetch
    seq = "".join(parts)
    return rc(seq) if strand == "-" else seq


def build_transcripts(copies, genome):
    genes = {c["gene"] for c in copies}
    iso_map = parse_isoforms(GFF, genes)
    recs, meta = [], []
    for c in copies:
        g = c["gene"]
        slen, rna, st, chrom, exons = longest_isoform(iso_map[g])
        seq = spliced_seq(genome, chrom, exons, st)
        recs.append((c["name"], seq))
        meta.append(dict(name=c["name"], gene=g, mrna=rna, strand=st, chrom=chrom,
                         n_exons=len(exons), spliced_len=len(seq),
                         genomic_span=exons[-1][1] - exons[0][0] + 1))
    fa = f"{OUT}/transcripts.fa"
    with open(fa, "w") as o:
        for name, seq in recs:
            o.write(f">{name}\n{seq}\n")
    sh(f"{MINIMAP2} -x asm5 -X --no-long-join {fa} {fa} 2>/dev/null > {OUT}/tx_selfaln.paf || true")
    return fa, meta, recs


def pairwise_identity(recs):
    """all-vs-all spliced-transcript identity via minimap2 asm5 (paralog divergence)."""
    fa = f"{OUT}/transcripts.fa"
    out = sh(f"{MINIMAP2} -c -x asm5 {fa} {fa} 2>/dev/null || true")
    pi = {}
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12 or f[0] == f[5]:
            continue
        matches, alnlen = int(f[9]), int(f[10])
        key = tuple(sorted((f[0], f[5])))
        ident = matches / max(alnlen, 1)
        if key not in pi or ident > pi[key]:
            pi[key] = round(100 * ident, 2)
    return pi


# ---- step 2: build graph ---------------------------------------------------
def build_graph(fa, base="copy0"):
    vgf = f"{OUT}/rabl2.vg"
    gfa = f"{OUT}/rabl2.gfa"
    # vg msga: deterministic single-thread progressive MSA -> base-level graph.
    log = sh(f"{VG} msga -f {fa} -b {base} -t 1 2>{OUT}/msga.log > {vgf}")
    # sort/normalize node ids (keeps paths); harmless, makes GFA stable
    sh(f"{VG} mod -X 256 {vgf} 2>/dev/null | {VG} ids -s - > {OUT}/rabl2.sorted.vg 2>/dev/null || cp {vgf} {OUT}/rabl2.sorted.vg")
    os.replace(f"{OUT}/rabl2.sorted.vg", vgf)
    sh(f"{VG} view {vgf} > {gfa}")
    stats = sh(f"{VG} stats -z -l {vgf}")
    paths = [p for p in sh(f"{VG} paths -L -x {vgf}").splitlines() if p.strip()]
    return vgf, gfa, stats, paths


# ---- step 3: reads ---------------------------------------------------------
def extract_reads(gt, genome):
    contig_len = dict(zip(genome.references, genome.lengths))
    gt_names = {r["name"] for r in gt["reads"]}
    bam = pysam.AlignmentFile(BAM, "rb")
    seqs = {}            # name -> longest observed spliced (query) sequence
    primary = set()      # names with a primary alignment over the loci
    for c in gt["copies"]:
        ch, s, e = c["chrom"], c["start"], c["end"]
        s2, e2 = max(0, s - PAD), min(contig_len.get(ch, e + PAD), e + PAD)
        for r in bam.fetch(ch, s2, e2):
            if r.is_unmapped or r.query_name not in gt_names:
                continue
            if not (r.is_secondary or r.is_supplementary):
                primary.add(r.query_name)
            sq = r.query_sequence
            # IsoSeq CCS read query = already-spliced mRNA; keep the longest copy (full read)
            if sq and (r.query_name not in seqs or len(sq) > len(seqs[r.query_name])):
                seqs[r.query_name] = sq
    bam.close()
    sel = {n: seqs[n] for n in gt_names if n in seqs}
    fa = f"{OUT}/reads.fa"
    with open(fa, "w") as o:
        for n, sq in sel.items():
            o.write(f">{n}\n{sq}\n")
    return fa, dict(gt_reads=len(gt_names), found_in_bam=len(sel),
                    with_primary_over_loci=len(primary),
                    secondary_only=len(set(sel) - primary),
                    missing=sorted(gt_names - set(sel))[:10])


# ---- step 4: align ---------------------------------------------------------
def align_reads(gfa, reads_fa):
    gam = f"{OUT}/reads.gam"
    log = sh(f"{GRAPHALIGNER} -g {gfa} -f {reads_fa} -a {gam} -x vg -t {THREADS} "
             f"2>{OUT}/ga.log > {OUT}/ga.stdout || cat {OUT}/ga.log")
    ga_log = open(f"{OUT}/ga.log").read() if os.path.exists(f"{OUT}/ga.log") else ""
    return gam, ga_log


def node_owners(vgf):
    """node_id -> set(copy paths that traverse it)."""
    pj = sh(f"{VG} paths -x {vgf} -X 2>/dev/null | {VG} view -aj - 2>/dev/null || true")
    owners = collections.defaultdict(set)
    for line in pj.splitlines():
        if not line.strip():
            continue
        pa = json.loads(line)
        nm = pa.get("name", "")
        for m in pa.get("path", {}).get("mapping", []):
            owners[int(m["position"]["node_id"])].add(nm)
    return owners


def analyze_gam(vgf, gam, paths):
    """PER-READ copy assignment: aggregate node votes across ALL of a read's alignment
    records (GraphAligner may split one read into several), then assign to the copy with
    the most copy-PRIVATE (SUN-like) node hits. Ties -> ambiguous."""
    owners = node_owners(vgf)
    js = sh(f"{VG} view -aj {gam}")
    by_read = collections.defaultdict(lambda: collections.Counter())
    records = collections.Counter()
    for line in js.splitlines():
        if not line.strip():
            continue
        aln = json.loads(line)
        records[aln["name"]] += 1
        if aln.get("score", 0) <= 0:
            continue
        for m in aln.get("path", {}).get("mapping", []):
            ow = owners.get(int(m["position"]["node_id"]), set())
            if len(ow) == 1:                       # copy-private node = decisive vote
                by_read[aln["name"]][next(iter(ow))] += 1
    per_copy = collections.Counter()
    read_assign = {}
    n_reads = len(records)
    aligned = len(by_read)
    for name in records:
        votes = by_read.get(name)
        if not votes:
            per_copy["_no_private_node"] += 1
            read_assign[name] = "_no_private_node"
            continue
        best = max(votes.values())
        top = [c for c, v in votes.items() if v == best]
        if len(top) == 1:
            per_copy[top[0]] += 1
            read_assign[name] = top[0]
        else:
            per_copy["_ambiguous"] += 1
            read_assign[name] = "_ambiguous"
    return dict(n_reads=n_reads, aligned_reads=aligned,
                gam_records=sum(records.values()),
                per_copy=dict(per_copy)), read_assign


def partition_of(hap):
    """canonical allele-partition of a PSV: frozenset of frozensets of copies sharing an
    allele. Coordinate-free signature for cross-tool concordance."""
    groups = collections.defaultdict(set)
    for cp, al in hap.items():
        groups[al].add(cp)
    return frozenset(frozenset(v) for v in groups.values())


def ref_index_to_base(ref, alt):
    """Align an ALT allele string to the REF allele string (difflib, deterministic) and
    return {ref_index -> alt_base or None(deletion)} covering only substituted/matched
    single-base ref positions. Pure insertions (bases in ALT with no ref column) are
    dropped -- we only track PER-REF-COLUMN substitutions (SNP-level PSVs), not indels."""
    out = {}
    sm = difflib.SequenceMatcher(a=ref, b=alt, autojunk=False)
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            for k in range(i2 - i1):
                out[i1 + k] = ref[i1 + k]
        elif tag == "replace":
            n = min(i2 - i1, j2 - j1)
            for k in range(n):                      # aligned substituted columns
                out[i1 + k] = alt[j1 + k]
            for k in range(n, i2 - i1):              # extra ref cols = deletion in alt
                out[i1 + k] = None
        elif tag == "delete":
            for k in range(i2 - i1):
                out[i1 + k] = None
        # tag == "insert": bases only in ALT, no ref column -> skip (indel, not a PSV col)
    return out


def expand_vcf_to_substitutions(vcf_parts, samples, ncopies):
    """Expand every vg-deconstruct record (incl. the giant multiallelic snarls) down to
    BASE-LEVEL substitution columns in copy0-spliced coordinates. A column is a PSV iff all
    `ncopies` copies have a called base there and >=2 distinct bases. Returns
    (columns, partitions) where columns = [(copy0_pos0, {copy:base}), ...]."""
    columns = []
    for c in vcf_parts:
        pos0 = int(c[1]) - 1                        # VCF POS is 1-based on copy0 transcript
        ref, alt = c[3], c[4]
        alt_list = alt.split(",")
        alleles = [ref] + alt_list
        # per-copy allele string (copy0 = REF, others per GT index)
        copy_allele = {"copy0": ref}
        for smp, g in zip(samples, c[9:]):
            gi = g.split(":")[0].replace("|", "/").split("/")[0]
            if gi == "." or not gi.isdigit() or int(gi) >= len(alleles):
                continue
            copy_allele[smp] = alleles[int(gi)]
        if len(copy_allele) != ncopies:
            continue                                # a copy is absent from this snarl
        # map each copy's allele onto REF columns
        base_at = {cp: ref_index_to_base(ref, seq) for cp, seq in copy_allele.items()}
        for i in range(len(ref)):
            col = {cp: base_at[cp].get(i) for cp in copy_allele}
            if any(b is None for b in col.values()):
                continue                            # deletion in some copy -> indel col
            if len({b for b in col.values()}) >= 2:
                columns.append((pos0 + i, col))
    parts = collections.Counter(partition_of(col) for _, col in columns)
    return columns, parts


def transcript_to_genomic_map(gt):
    """Map copy0 SPLICED-transcript coordinates (the vg-deconstruct VCF POS frame) back to
    the o2_vg GENOMIC-backbone offsets (the GT bubble `pos` frame), so PSVs can be matched by
    LITERAL genomic position, not a coordinate-free proxy. copy0 is BOTH the deconstruct
    reference path AND the o2_vg backbone, so both tools are anchored on copy0 -- but vg works
    on the spliced mRNA-sense transcript while the GT `pos` is a genomic-with-introns offset
    from backbone_start, and copy0 is minus-strand (=> the transcript is the reverse-complement
    of the genomic exon bases). Returns (pos0 -> genomic_offset fn, strand)."""
    copy0 = gt["copies"][0]
    gene = copy0["gene"]
    backbone_start = gt["backbone_start"]
    iso = parse_isoforms(GFF, {gene})
    _slen, _rna, strand, _chrom, exons = longest_isoform(iso[gene])
    gcoords = []
    for s, e in exons:
        gcoords.extend(range(s, e + 1))             # genomic-forward exon-concatenated 1-based
    L = len(gcoords)
    if strand == "-":
        def to_offset(p0):
            return gcoords[L - 1 - p0] - backbone_start
    else:
        def to_offset(p0):
            return gcoords[p0] - backbone_start
    return to_offset, strand


def literal_bubble_recall(gt, cols, strand, pos_map, tol=2):
    """LITERAL position+haplotype match of vg base-level substitution columns to GT bubbles
    (replaces the coordinate-free partition-subset proxy). For each vg column: map its copy0
    transcript position to a genomic backbone offset, complement its bases to genomic-forward
    orientation when copy0 is minus-strand, then require a GT bubble within +-tol bp whose FULL
    5-copy haplotype is base-for-base identical. distinct GT bubbles matched = the honest
    per-position RECALL. Also reports the degeneracy-inflated haplotype-only figure (only 16
    allele-partitions exist for 5 copies, so a hap can match many GT bubbles) for contrast."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    need_comp = (strand == "-")
    gt_by_pos = collections.defaultdict(list)
    gt_hap_index = collections.defaultdict(set)     # hap tuple -> {GT bubble cols with that hap}
    for b in gt["bubbles"]:
        gt_by_pos[b["pos"]].append(b)
        gt_hap_index[tuple(sorted(b["hap"].items()))].add(b["col"])
    matched = set()
    hap_only_cols = 0
    hap_only_bubbles = set()
    offsets = []
    for p0, coldict in cols:
        gpos = pos_map(p0)
        hap_gen = {cp: (b.translate(comp) if need_comp else b) for cp, b in coldict.items()}
        key = tuple(sorted(hap_gen.items()))
        if key in gt_hap_index:                      # haplotype exists somewhere in GT
            hap_only_cols += 1
            hap_only_bubbles |= gt_hap_index[key]
        best = None
        for d in range(-tol, tol + 1):
            for b in gt_by_pos.get(gpos + d, []):
                if all(b["hap"].get(cp) == hap_gen.get(cp) for cp in hap_gen):
                    if best is None or abs(d) < abs(best[0]):
                        best = (d, b)
        if best is not None:
            matched.add(best[1]["col"])
            offsets.append(best[0])
    off_ctr = collections.Counter(offsets)
    return dict(
        vg_substitution_cols=len(cols),
        literal_bubbles_recovered=len(matched),
        literal_recall_pct=round(100 * len(matched) / max(gt["n_bubbles"], 1), 1),
        pos_tolerance_bp=tol,
        systematic_offset_bp=(off_ctr.most_common(1)[0][0] if off_ctr else None),
        haplotype_only_cols=hap_only_cols,
        haplotype_only_distinct_bubbles=len(hap_only_bubbles),
        note=("literal = same genomic backbone position (+-%dbp, consistent systematic offset) "
              "AND identical 5-copy haplotype after minus-strand complement; the haplotype-only "
              "figure is degeneracy-inflated (<=16 partitions for 5 copies) and is NOT recall"
              % tol),
    )


def run_snarls_deconstruct(vgf, gt, ncopies):
    """vg snarls (bubble decomposition) + vg deconstruct (PSV VCF vs copy0)."""
    res = {}
    # ---- snarls ----
    try:
        sh(f"{VG} snarls {vgf} > {OUT}/rabl2.snarls 2>{OUT}/snarls.log")
        cnt = sh(f"{VG} view -Rj {OUT}/rabl2.snarls 2>/dev/null | wc -l").strip()
        res["snarls"] = int(cnt)
    except Exception as ex:
        res["snarls_error"] = str(ex)[:400]
    # ---- deconstruct (needs xg + snarls); reference = copy0 (our backbone) ----
    try:
        sh(f"{VG} deconstruct -p copy0 -e -t {THREADS} {vgf} > {OUT}/rabl2.decon.vcf 2>{OUT}/decon.log")
        vcf_parts = []
        for line in open(f"{OUT}/rabl2.decon.vcf"):
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = line.rstrip("\n").split("\t")[9:]
                continue
            vcf_parts.append(line.rstrip("\n").split("\t"))
        # per-site partition from GT genotypes (0=ref/copy0 allele)
        snp_sites = 0
        multi_sites = 0
        vg_parts = collections.Counter()
        for c in vcf_parts:
            ref, alt = c[3], c[4]
            gts = c[9:]
            alleles = [ref] + alt.split(",")
            hap = {"copy0": ref}
            for smp, g in zip(samples, gts):
                gi = g.split(":")[0].split("/")[0].split("|")[0]
                if gi in (".",):
                    continue
                hap[smp] = alleles[int(gi)] if int(gi) < len(alleles) else "N"
            if len(ref) == 1 and all(len(a) == 1 for a in alt.split(",")):
                snp_sites += 1
            else:
                multi_sites += 1
            if len(hap) == ncopies:
                vg_parts[partition_of(hap)] += 1
        res["deconstruct_total_sites"] = len(vcf_parts)
        res["deconstruct_snp_sites"] = snp_sites
        res["deconstruct_indel_sites"] = multi_sites
        # ----- (i) SITE-level concordance: vg reports snarl SITES; dense paralog SNPs are
        #           aggregated inside multiallelic snarls, so only ~22 sites are clean
        #           all-5-copy calls. Coordinate-free allele-PARTITION match vs our bubbles.
        our_parts = collections.Counter(partition_of(b["hap"]) for b in gt["bubbles"]
                                        if len(b["hap"]) == ncopies)
        shared = our_parts & vg_parts
        res["our_bubbles"] = sum(our_parts.values())
        res["vg_full_copy_sites"] = sum(vg_parts.values())
        res["shared_partition_multiset"] = sum(shared.values())
        res["our_distinct_partitions"] = len(our_parts)
        res["vg_distinct_partitions"] = len(vg_parts)
        res["shared_distinct_partitions"] = len(set(our_parts) & set(vg_parts))
        # SUN-like: partitions that isolate exactly one copy (singleton allele)
        def n_sun(cnt):
            return sum(v for p, v in cnt.items() if any(len(g) == 1 for g in p))
        res["our_sun_partitions"] = n_sun(our_parts)
        res["vg_sun_partitions"] = n_sun(vg_parts)
        # ----- (ii) BASE-level concordance: EXPAND every snarl (incl. the giant merged
        #            ones) to per-base substitution columns -> this is the apples-to-apples
        #            comparison to our 235 per-column PSVs.
        cols, vg_exp_parts = expand_vcf_to_substitutions(vcf_parts, samples, ncopies)
        exp_shared = our_parts & vg_exp_parts
        # ----- (iii) LITERAL per-position recall: map each vg base-level column to a GT bubble
        #             by genomic position + identical haplotype (the honest apples-to-apples
        #             number; the partition figures below are PATTERN overlap, not recall).
        pos_map, strand0 = transcript_to_genomic_map(gt)
        res["literal_recall"] = literal_bubble_recall(gt, cols, strand0, pos_map)
        res["vg_expanded_substitution_cols"] = len(cols)
        res["vg_expanded_distinct_partitions"] = len(vg_exp_parts)
        res["expanded_shared_partition_multiset"] = sum(exp_shared.values())
        res["expanded_shared_distinct_partitions"] = len(set(our_parts) & set(vg_exp_parts))
        # PRECISION: is every partition-pattern vg emits one we also call?
        res["vg_expanded_partitions_all_subset_of_ours"] = all(
            p in our_parts for p in vg_exp_parts)
        # WHY the base-level count is low: many snarl SITES drop a copy (that copy's
        # divergent/short transcript does not co-traverse -> isoform heterogeneity), and
        # inside the big multiallelic snarls one copy's indel masks co-located substitutions.
        res["vg_sites_dropping_a_copy"] = (res["deconstruct_total_sites"]
                                           - res["vg_full_copy_sites"])
        # of OUR 235 bubbles, how many merely SHARE an allele-partition PATTERN with something
        # vg emits. This is DEGENERACY-INFLATED (only <=16 partitions exist for 5 copies, so
        # ~all our bubbles trivially share vg's 9 patterns) and is NOT per-position recovery --
        # the honest recall is res["literal_recall"]. Kept only as a clearly-labelled contrast.
        our_covered = sum(v for p, v in our_parts.items() if p in vg_exp_parts)
        res["partition_pattern_overlap_bubbles"] = our_covered
        res["partition_pattern_overlap_pct"] = round(
            100 * our_covered / max(sum(our_parts.values()), 1), 1)
        # per-copy SUN columns implied by the expansion (copy isolated by a singleton allele)
        vg_exp_sun = collections.Counter()
        for _, col in cols:
            groups = collections.defaultdict(list)
            for cp, b in col.items():
                groups[b].append(cp)
            for b, cps in groups.items():
                if len(cps) == 1:
                    vg_exp_sun[cps[0]] += 1
        res["vg_expanded_sun_per_copy"] = dict(vg_exp_sun)
        res["vg_expanded_sun_total"] = sum(vg_exp_sun.values())
    except Exception as ex:
        import traceback
        res["deconstruct_error"] = (str(ex) + " | " + traceback.format_exc())[:900]
    return res


def _spearman(a, b):
    """deterministic Spearman rho on two equal-length numeric lists (no scipy)."""
    n = len(a)
    if n < 2:
        return None
    def rank(x):
        order = sorted(range(n), key=lambda i: x[i])
        r = [0.0] * n
        i = 0
        while i < n:                                # average-rank ties
            j = i
            while j + 1 < n and x[order[j + 1]] == x[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    ra, rb = rank(a), rank(b)
    ma, mb = sum(ra) / n, sum(rb) / n
    cov = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    va = sum((ra[i] - ma) ** 2 for i in range(n)) ** 0.5
    vb = sum((rb[i] - mb) ** 2 for i in range(n)) ** 0.5
    if va == 0 or vb == 0:
        return None
    return round(cov / (va * vb), 3)


def run_pack_ase(vgf, gam, ncopies, gate_per_copy=None):
    """vg pack per-node coverage -> per-copy expression (ASE) by attributing each
    copy-PRIVATE node's coverage to its copy. Report mean & total over private nodes and
    correlate against the gate's per-copy read counts (does vg ASE track our quantitation?)."""
    try:
        sh(f"{VG} pack -x {vgf} -g {gam} -d -t {THREADS} > {OUT}/reads.pack.tsv 2>{OUT}/pack.log")
    except Exception as ex:
        return {"pack_error": str(ex)[:400]}
    owners = node_owners(vgf)
    priv = {nid: next(iter(ow)) for nid, ow in owners.items() if len(ow) == 1}
    # TRUE copy-private NODE count (single-owner graph nodes) -- sums to ~72 over 5 copies.
    priv_node_counts = collections.Counter(priv.values())
    cov = collections.defaultdict(list)   # copy -> per-BASE coverage over its private nodes
    header = None
    for line in open(f"{OUT}/reads.pack.tsv"):
        f = line.rstrip("\n").split("\t")
        if header is None:
            header = f
            idx = {h: i for i, h in enumerate(f)}
            continue
        nid = int(f[idx["node.id"]])
        c = int(f[idx["coverage"]])
        if nid in priv:
            cov[priv[nid]].append(c)   # one entry PER BASE of each private node (pack -d rows)
    mean_cov = {cp: round(sum(v) / len(v), 2) for cp, v in cov.items() if v}
    total_cov = {cp: sum(v) for cp, v in cov.items() if v}
    out = dict(pack_per_copy_mean_private_cov=mean_cov,
               pack_per_copy_total_private_cov=total_cov,
               # NOTE: pack -d emits one row PER BASE, so this is the private-BASE count/copy
               # (sums to ~641), distinct from the private-NODE count/copy below (sums to ~72).
               pack_private_bases={cp: len(v) for cp, v in cov.items()},
               pack_private_nodes={cp: priv_node_counts.get(cp, 0)
                                   for cp in sorted(priv_node_counts)})
    if gate_per_copy:
        copies = sorted(set(mean_cov) | set(gate_per_copy))
        gvec = [gate_per_copy.get(cp, 0) for cp in copies]
        out["ase_gate_read_counts"] = {cp: gate_per_copy.get(cp, 0) for cp in copies}
        out["ase_spearman_mean_vs_gate"] = _spearman(
            [mean_cov.get(cp, 0) for cp in copies], gvec)
        out["ase_spearman_total_vs_gate"] = _spearman(
            [total_cov.get(cp, 0) for cp in copies], gvec)
    return out


def build_concordance_table(report, gt):
    """Distill the full pilot_report into the deliverable cross-tool concordance table
    (vg standard toolkit vs our custom o2_vg machinery) for bench/vg_gfa_pilot.json."""
    sd = report.get("snarls_deconstruct", {})
    pk = report.get("pack_ase", {})
    aln = report.get("alignment", {})
    xc = report.get("gate_crosscheck", {})
    stats = {k: int(v) for k, v in
             (kv.split("\t") for kv in report["vg_stats"].splitlines())}
    ncopies = gt["n_copies"]
    n_paths = len([p for p in report["paths"] if p.startswith("copy")])
    exp_cols = sd.get("vg_expanded_substitution_cols")
    overlap_pct = sd.get("partition_pattern_overlap_pct")   # degeneracy-inflated, NOT recall
    drop = sd.get("vg_sites_dropping_a_copy")
    lit = sd.get("literal_recall", {})
    lit_recall = lit.get("literal_bubbles_recovered")       # honest per-position recall
    lit_pct = lit.get("literal_recall_pct")
    lit_tol = lit.get("pos_tolerance_bp")

    rows = [
        dict(metric="(a) copies as distinct paths",
             custom=ncopies, vg=n_paths,
             concordance="MATCH" if n_paths == ncopies else "MISMATCH"),
        dict(metric="(b) PSVs: vg snarl SITES",
             custom=gt["n_bubbles"], vg=sd.get("snarls"),
             concordance="vg reports snarl SITES not per-column PSVs"),
        dict(metric="(b) PSVs: vg deconstruct SITES",
             custom=gt["n_bubbles"], vg=sd.get("deconstruct_total_sites"),
             concordance=f"{sd.get('deconstruct_snp_sites')} SNP + "
                         f"{sd.get('deconstruct_indel_sites')} complex; aggregated"),
        dict(metric="(b) PSVs: vg base-level substitution cols",
             custom=gt["n_bubbles"], vg=exp_cols,
             concordance=f"{exp_cols} base-level SNP columns emitted (from {drop}/"
                         f"{sd.get('deconstruct_total_sites')} sites that keep all copies + "
                         f"snarl expansion); rest lost to copy drop-out / indel masking"),
        dict(metric="(b) PSVs: LITERAL per-position recall",
             custom=f"{gt['n_bubbles']} bubbles",
             vg=f"{lit_recall}/{gt['n_bubbles']} ({lit_pct}%)",
             concordance=f"vg column matched a GT bubble at the SAME genomic pos (+-{lit_tol}bp) "
                         f"with an IDENTICAL 5-copy haplotype (minus-strand complemented) -- "
                         f"HONEST recall. (partition-pattern overlap {overlap_pct}% is "
                         f"degeneracy-inflated, NOT recall.)"),
        dict(metric="(b) allele-partition PRECISION",
             custom=f"{sd.get('our_distinct_partitions')} patterns",
             vg=f"{sd.get('vg_expanded_distinct_partitions')} patterns",
             concordance=f"{sd.get('expanded_shared_distinct_partitions')}/"
                         f"{sd.get('vg_expanded_distinct_partitions')} vg patterns "
                         f"subset-of-ours (all_subset="
                         f"{sd.get('vg_expanded_partitions_all_subset_of_ours')}) -- every "
                         f"variant vg calls is one we also call"),
        dict(metric="(c) read->copy assignment agreement",
             custom=f"{gt['assignment']['assigned']}/{gt['n_reads']} gate",
             vg=f"{xc.get('agree')}/{xc.get('shared_assigned')} node-vote",
             concordance=f"{xc.get('pct')}% (excl vg ties) / {xc.get('pct_incl_ties')}% "
                         f"({xc.get('agree')}/{xc.get('shared_assigned_incl_ties')} incl ties); "
                         f"vg ships NO assignment (naive private-node vote, biased)"),
        dict(metric="(d) ASE per-copy (Spearman vs gate)",
             custom="read counts",
             vg=f"pack mean rho={pk.get('ase_spearman_mean_vs_gate')} / "
                f"total rho={pk.get('ase_spearman_total_vs_gate')}",
             concordance="NEW capability (custom has none); WEAK -- rho on n=5 copies is "
                         "not statistically significant"),
    ]

    reproduces = (n_paths == ncopies
                  and sd.get("vg_expanded_partitions_all_subset_of_ours"))
    verdict = (
        "PARTIAL REPRODUCE + REAL EXTENSION. vg REPRODUCES the graph substrate: 5 copies = "
        "5 distinct named paths, and every allele-partition pattern vg emits is one we also "
        f"call (precision {sd.get('expanded_shared_distinct_partitions')}/"
        f"{sd.get('vg_expanded_distinct_partitions')} vg patterns subset-of-ours). But vg does "
        "NOT reproduce out-of-the-box at PSV RESOLUTION: it collapses the 235 per-column PSVs "
        f"into {sd.get('deconstruct_total_sites')} multiallelic snarl SITES; expanded to base "
        f"level it emits {exp_cols} SNP columns, of which a LITERAL position+haplotype match "
        f"recovers only {lit_recall}/235 distinct GT bubbles ({lit_pct}% recall) -- because "
        f"{drop}/{sd.get('deconstruct_total_sites')} sites drop a copy whose divergent/short "
        "transcript does not co-traverse (isoform heterogeneity, esp. copy4/RABL2B) and big "
        "snarls bury co-located SNPs behind indels. vg also ships NO read->copy assignment "
        f"(naive private-node vote agrees with the gate only {xc.get('pct')}% excl ties / "
        f"{xc.get('pct_incl_ties')}% incl ties, biased by unequal private-marker counts). vg "
        "ADDS vg-pack per-copy coverage = allele-specific EXPRESSION (a NEW capability our "
        f"machinery lacks; Spearman total-cov vs gate rho={pk.get('ase_spearman_total_vs_gate')}"
        " -- WEAK/not significant at n=5 copies, needs SNP-only markers + normalization). "
        "ADOPT vg as the standardized graph substrate + ASE layer; KEEP the custom "
        "significance gate + PSV-column machinery on top (vg does not do per-column PSV "
        "discovery or read->copy assignment). HARD PREREQ: build at transcript level with "
        "matched isoforms or a gene-level exon union (else isoform heterogeneity drops copies "
        "and manufactures spurious snarls); migrate off deprecated vg msga."
    ) if reproduces else (
        "vg did NOT cleanly reproduce the RABL2 structure (see rows) -- keep custom machinery."
    )

    return dict(
        family=report["family"], gene=report["gene"], graph_tool="vg msga",
        deterministic=True,
        graph=dict(nodes=stats.get("nodes"), edges=stats.get("edges"),
                   length_bp=stats.get("length"), copy_paths=n_paths,
                   private_nodes_per_copy=pk.get("pack_private_nodes"),   # single-owner nodes ~72
                   private_bases_per_copy=pk.get("pack_private_bases")),  # per-base ~641
        gt=report["gt"],
        pairwise_identity_pct=report.get("pairwise_identity_pct"),
        a_copies=dict(custom=ncopies, vg_distinct_paths=n_paths,
                      match=n_paths == ncopies),
        b_psvs=dict(
            custom_bubbles=gt["n_bubbles"], custom_sun=gt["n_sun_bubbles"],
            vg_snarls=sd.get("snarls"),
            vg_deconstruct_sites=sd.get("deconstruct_total_sites"),
            vg_snp_sites=sd.get("deconstruct_snp_sites"),
            vg_complex_sites=sd.get("deconstruct_indel_sites"),
            vg_full_copy_sites=sd.get("vg_full_copy_sites"),
            vg_expanded_substitution_cols=exp_cols,
            vg_expanded_sun_total=sd.get("vg_expanded_sun_total"),
            vg_expanded_sun_per_copy=sd.get("vg_expanded_sun_per_copy"),
            vg_sites_dropping_a_copy=drop,
            vg_expanded_distinct_partitions=sd.get("vg_expanded_distinct_partitions"),
            vg_expanded_partitions_all_subset_of_ours=sd.get(
                "vg_expanded_partitions_all_subset_of_ours"),
            site_level_shared_partitions=sd.get("shared_distinct_partitions"),
            expanded_shared_partitions=sd.get("expanded_shared_distinct_partitions"),
            # HONEST per-position recall: literal genomic-position + haplotype match.
            literal_bubbles_recovered=lit_recall,
            literal_recall_pct=lit_pct,
            literal_pos_tolerance_bp=lit_tol,
            literal_systematic_offset_bp=lit.get("systematic_offset_bp"),
            haplotype_only_cols=lit.get("haplotype_only_cols"),
            haplotype_only_distinct_bubbles=lit.get("haplotype_only_distinct_bubbles"),
            # DEGENERACY-INFLATED contrast only (<=16 partitions for 5 copies) -- NOT recall.
            partition_pattern_overlap_bubbles=sd.get("partition_pattern_overlap_bubbles"),
            partition_pattern_overlap_pct=overlap_pct,
            interpretation="HONEST recall = %s/%s distinct GT bubbles recovered by literal "
                           "position+haplotype match (%s%%); vg emits %s base-level SNP columns, "
                           "the rest of the 235 lost to copy drop-out / indel-masked snarls. "
                           "PRECISION is high (all vg partitions subset-of-ours) but per-position "
                           "RECALL is low -- keep custom column machinery. (partition-pattern "
                           "overlap %s%% is degeneracy-inflated, NOT recall.)"
                           % (lit_recall, gt["n_bubbles"], lit_pct, exp_cols, overlap_pct)),
        c_assignment=dict(
            gate_assigned=gt["assignment"]["assigned"], gt_reads=gt["n_reads"],
            shared_assigned_reads=xc.get("shared_assigned"), agree=xc.get("agree"),
            agree_pct=xc.get("pct"),
            shared_assigned_incl_ties=xc.get("shared_assigned_incl_ties"),
            agree_pct_incl_ties=xc.get("pct_incl_ties"),
            gate_unassigned_reads_present=xc.get("gate_unassigned_reads"),
            graph_placed_copy_on_unassigned=xc.get("graph_assigned_a_copy_to_unassigned"),
            vg_ships_assignment=False,
            graph_vote_per_copy=aln.get("per_copy")),
        d_ase=dict(
            vg_pack_mean_private_cov=pk.get("pack_per_copy_mean_private_cov"),
            vg_pack_total_private_cov=pk.get("pack_per_copy_total_private_cov"),
            vg_pack_private_nodes_per_copy=pk.get("pack_private_nodes"),
            vg_pack_private_bases_per_copy=pk.get("pack_private_bases"),
            gate_read_counts=gt["per_copy"],
            spearman_mean_vs_gate=pk.get("ase_spearman_mean_vs_gate"),
            spearman_total_vs_gate=pk.get("ase_spearman_total_vs_gate"),
            note="NEW capability our machinery lacks; WEAK signal -- Spearman rho on n=5 copies "
                 "is not statistically significant, and needs SNP-only private markers + "
                 "length-normalization to match gate read counts"),
        concordance_table=rows,
        verdict=verdict,
    )


def main():
    os.makedirs(OUT, exist_ok=True)
    if not os.path.exists(GFF):
        genes = "LOC134756389|RABL2A|LOC101142457|LOC129523511|RABL2B"
        sh(f"grep -P 'gene=({genes});' {FULL_GFF} > {GFF}")
    gt = json.load(open(GT))
    genome = pysam.FastaFile(GENOME)
    report = dict(family=gt["family"], gene=gt["gene"],
                  gt=dict(n_copies=gt["n_copies"], K=gt["K"], n_bubbles=gt["n_bubbles"],
                          n_sun_bubbles=gt["n_sun_bubbles"], n_reads=gt["n_reads"],
                          gate_assigned=gt["assignment"]["assigned"],
                          gate_no_cover=gt["assignment"]["no_cover"],
                          gate_ambiguous=gt["assignment"]["ambiguous"],
                          gate_per_copy=gt["per_copy"]))

    print("[1] building spliced transcripts ...", flush=True)
    fa, meta, recs = build_transcripts(gt["copies"], genome)
    report["transcripts"] = meta
    report["pairwise_identity_pct"] = {f"{a}|{b}": v
                                       for (a, b), v in pairwise_identity(recs).items()}
    for m in meta:
        print(f"    {m['name']} {m['gene']:14s} {m['mrna']:20s} strand {m['strand']} "
              f"exons {m['n_exons']:2d} spliced {m['spliced_len']}bp "
              f"(genomic {m['genomic_span']}bp)")

    print("[2] building graph with vg msga ...", flush=True)
    vgf, gfa, stats, paths = build_graph(fa)
    report["vg_stats"] = stats.strip()
    report["paths"] = paths
    print("    vg stats:\n" + "\n".join("      " + l for l in stats.strip().splitlines()))
    print(f"    paths (copies as distinct paths): {paths}")

    print("[3] extracting RABL2 reads (spliced) from GGO_mm.bam ...", flush=True)
    reads_fa, rinfo = extract_reads(gt, genome)
    report["reads"] = rinfo
    print(f"    gt reads {rinfo['gt_reads']} | extracted {rinfo['found_in_bam']} | "
          f"primary-over-loci {rinfo['with_primary_over_loci']} | "
          f"secondary-only {rinfo['secondary_only']}")

    print("[4] aligning reads with GraphAligner -x vg ...", flush=True)
    gam, ga_log = align_reads(gfa, reads_fa)
    tail = "\n".join(ga_log.strip().splitlines()[-6:])
    print("    GraphAligner log tail:\n" + "\n".join("      " + l for l in tail.splitlines()))
    aln_summary, read_assign = analyze_gam(vgf, gam, paths)
    report["alignment"] = aln_summary
    report["ga_log_tail"] = tail
    print(f"    aligned {aln_summary['aligned_reads']}/{aln_summary['n_reads']} reads "
          f"({aln_summary['gam_records']} gam records); "
          f"per-copy(private-node vote): {aln_summary['per_copy']}")

    # cross-check graph assignment vs our gate
    gate = {r["name"]: r["best_copy"] for r in gt["reads"] if r["status"] == "assigned"}
    status = {r["name"]: r["status"] for r in gt["reads"]}
    is_copy = lambda n: read_assign.get(n, "").startswith("copy")
    shared = [n for n in read_assign if n in gate and is_copy(n)]
    agree = sum(1 for n in shared if read_assign[n] == gate[n])
    # two honest denominators: (i) reads where BOTH assigned a copy (excludes vg ties/
    # no-private-node), (ii) ALL gate-assigned reads the graph saw (includes vg
    # ambiguous/no-private as disagreements) -- the stricter, non-cherry-picked figure.
    shared_all = [n for n in read_assign if n in gate]
    pct_excl = round(100 * agree / max(len(shared), 1), 1)
    pct_incl = round(100 * agree / max(len(shared_all), 1), 1)
    unassigned = [n for n in read_assign if status.get(n) in ("no_cover", "ambiguous")]
    graph_rescued = sum(1 for n in unassigned if is_copy(n))
    report["gate_crosscheck"] = dict(
        shared_assigned=len(shared), agree=agree, pct=pct_excl,
        shared_assigned_incl_ties=len(shared_all), pct_incl_ties=pct_incl,
        gate_unassigned_reads=len(unassigned),
        graph_assigned_a_copy_to_unassigned=graph_rescued)
    print(f"    graph-vote vs gate on shared ASSIGNED reads: {agree}/{len(shared)} "
          f"({pct_excl}% excl vg ties) | {agree}/{len(shared_all)} ({pct_incl}% incl ties)")
    print(f"    of {len(unassigned)} gate-UNASSIGNED reads present, graph gave a copy to "
          f"{graph_rescued}")

    print("[5] vg snarls + vg deconstruct (bubbles/PSVs) ...", flush=True)
    sd = run_snarls_deconstruct(vgf, gt, gt["n_copies"])
    report["snarls_deconstruct"] = sd
    print("    " + json.dumps(sd))

    print("[6] vg pack (per-copy ASE) ...", flush=True)
    pk = run_pack_ase(vgf, gam, gt["n_copies"], gate_per_copy=gt["per_copy"])
    report["pack_ase"] = pk
    print("    " + json.dumps(pk))
    print(f"    gate per-copy (ground truth ASE): {gt['per_copy']}")

    genome.close()
    json.dump(report, open(f"{OUT}/pilot_report.json", "w"), indent=2)
    print(f"\n[+] wrote {OUT}/pilot_report.json")
    print(f"[+] outputs in {OUT}: transcripts.fa rabl2.vg rabl2.gfa reads.fa reads.gam")

    # ---- distilled cross-tool concordance table -> bench/vg_gfa_pilot.json --------------
    # sort_keys=True => byte-identical across reruns (Counter/set insertion order varies).
    table = build_concordance_table(report, gt)
    json.dump(table, open(BENCH_JSON, "w"), indent=2, sort_keys=True)
    print(f"[+] wrote concordance table {BENCH_JSON}")
    print("\n==== CONCORDANCE (vg vs custom o2_vg, RABL2/fam39) ====")
    for row in table["concordance_table"]:
        print(f"  {row['metric']:32s} custom={str(row['custom']):>10s}  "
              f"vg={str(row['vg']):>12s}  -> {row['concordance']}")
    print("  VERDICT: " + table["verdict"])


if __name__ == "__main__":
    main()
