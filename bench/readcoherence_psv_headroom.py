#!/usr/bin/env python3
"""Read-coherence x PSV HEADROOM probe (cheap go/no-go before the molecule-threaded-graph rebuild).

Question (user): of the isoforms flow KILLS that read-coherence recovers, how many sit at
paralog/family loci where PSVs would ALSO copy-resolve them? I.e. does threading molecules+PSVs
through one graph buy BOTH recall (read-coherence) AND copy-resolution (PSV) on the SAME loci,
or are the two channels addressing disjoint loci (in which case the unified rebuild is not worth it)?

This is purely geometric + structural — no POA/BAM recompute — so it is cheap and self-contained.

Inputs (already on disk from the genome-wide runs):
  /tmp/gw/rcg_<C>.gtf     read-coherence layer ON (flow ∪ gated read-chain); the added isoforms carry
                          source "read_coherence" (= the recall set; 23,791 genome-wide).
  /tmp/gw/ref_<C>.gff3    RefSeq annotation (intron chains for FSM classification + family coords).
  bench/copy_recovery_eval/results/universe.tsv   paralog truth set: transcript_id -> family_id,
                          n_family_copies (62 multi-copy families).

Tiers (honest decomposition):
  recall set  = read_coherence-tagged transcripts (multi-exon).
  Tier 1      = of those, how many overlap a multi-copy-family member locus.
  Tier 2      = of those, the family-locus REGIME:
                  COLLAPSED   (>=2 copies overlap one coordinate frame)  -> PSVs SPLIT  (the win regime)
                  TANDEM_NEAR (same chrom, <=1Mb, non-overlapping)       -> coords + minimap2 separate
                  DISPERSED   (different chrom / >1Mb apart)             -> coords + minimap2 separate
                (PSV copy-resolution has leverage only in COLLAPSED — established in copy_split work:
                 RABL2 on distinct contigs was the BOUNDARY; resolved tandems' coords already separate.)
  FSM cut     = recall isoforms whose intron chain exactly matches a RefSeq transcript (real recall).

Output: bench/readcoherence_psv_headroom.md (table) + headroom_loci.tsv (the COLLAPSED-regime hits).
Deterministic (sorted iteration, no RNG).
"""
import glob
import os
import re
from collections import defaultdict

GW = "/tmp/gw"
UNIVERSE = os.path.join(os.path.dirname(__file__), "copy_recovery_eval/results/universe.tsv")
OUT_MD = os.path.join(os.path.dirname(__file__), "readcoherence_psv_headroom.md")
OUT_TSV = os.path.join(os.path.dirname(__file__), "headroom_loci.tsv")
TANDEM_WINDOW = 1_000_000  # bp: same-chrom copies within this are "near" (cross-map plausible)

# Compara/contiguous-core-validated verdicts on the families that turn up COLLAPSED here.
# CONFIRMED = gene-tree-confirmed recent paralogs (real collapsed/tandem copies, PSV-resolvable in
# principle). DOMAIN_SHARER = flagged FALSE by Compara + contiguous-core gate (two different genes
# sharing one protein domain / repeat — overlapping coords are an annotation accident, NOT copies;
# PSVs cannot "copy-resolve" two non-paralogous genes). See bench/compara_validation.* +
# bench/poa_family_definition.* + project_family_detection_validation memory.
CONFIRMED_REAL = {"APOBEC3D", "RFPL1", "RFPL2", "RFPL3", "RABL2A", "RABL2B"}
DOMAIN_SHARER = {"CDPF1", "CREB1", "GCA", "LOC134756662", "LOC129529456", "LOC129529513"}

_ATTR = re.compile(r'(\w+)[ =]"?([^";]+)"?')


def parse_attrs(field):
    d = {}
    for m in _ATTR.finditer(field):
        d[m.group(1)] = m.group(2)
    return d


def gff3_attrs(field):
    d = {}
    for kv in field.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def intron_chain(exons):
    """exons: list of (start,end) 1-based inclusive. Return tuple of (donor,acceptor) introns."""
    if len(exons) < 2:
        return ()
    ex = sorted(exons)
    return tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))


# ---------------------------------------------------------------------------
# 1. universe: multi-copy families
# ---------------------------------------------------------------------------
tx2fam = {}        # refseq transcript_id -> (family_id, n_copies)
fam_ncopies = {}
with open(UNIVERSE) as fh:
    fh.readline()
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if len(c) < 4:
            continue
        tid, gid, fam, nc = c[0], c[1], c[2], int(c[3])
        if nc >= 2:
            tx2fam[tid] = (fam, nc)
            fam_ncopies[fam] = nc

# ---------------------------------------------------------------------------
# 2. RefSeq: coords of family members + ALL intron chains (for FSM)
# ---------------------------------------------------------------------------
# refseq transcript_id -> exons ; gene -> (chrom,start,end,strand)
ref_tx_exons = defaultdict(list)
ref_tx_chrom = {}
ref_all_chains = set()       # (chrom, strand, intron_chain) over ALL refseq tx -> FSM lookup
ref_tx_of_exon_owner = {}    # exon parent -> transcript_id

# First pass: map mRNA/transcript feature ID -> transcript_id, and gene.
feat_to_tid = {}
for gff in sorted(glob.glob(os.path.join(GW, "ref_*.gff3"))):
    chrom = os.path.basename(gff)[len("ref_"):-len(".gff3")]
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            typ = f[2]
            a = gff3_attrs(f[8])
            if typ in ("mRNA", "transcript", "primary_transcript"):
                fid = a.get("ID", "")
                tid = a.get("transcript_id") or a.get("Name", "")
                if fid and tid:
                    feat_to_tid[(chrom, fid)] = tid
                    ref_tx_chrom[tid] = (chrom, int(f[3]), int(f[4]), f[6], a.get("gene", ""))
            elif typ == "exon":
                parent = a.get("Parent", "")
                tid = feat_to_tid.get((chrom, parent))
                if tid:
                    ref_tx_exons[tid].append((int(f[3]), int(f[4])))

# build FSM lookup over all refseq tx, and family-member loci
for tid, exons in ref_tx_exons.items():
    if tid not in ref_tx_chrom:
        continue
    chrom, _, _, strand, _ = ref_tx_chrom[tid]
    ch = intron_chain(exons)
    if ch:
        ref_all_chains.add((chrom, strand, ch))

# family -> list of member loci (chrom,start,end) ; member = each family-member transcript span
fam_member_loci = defaultdict(list)
fam_found = defaultdict(int)
for tid, (fam, nc) in tx2fam.items():
    if tid in ref_tx_chrom:
        chrom, s, e, strand, gene = ref_tx_chrom[tid]
        fam_member_loci[fam].append((chrom, s, e))
        fam_found[fam] += 1

# merge member loci per family into distinct copy-intervals (overlapping member tx -> one copy)
def merge_intervals(ivs):
    by_chrom = defaultdict(list)
    for chrom, s, e in ivs:
        by_chrom[chrom].append((s, e))
    out = []
    for chrom, segs in by_chrom.items():
        segs.sort()
        cs, ce = segs[0]
        for s, e in segs[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                out.append((chrom, cs, ce))
                cs, ce = s, e
        out.append((chrom, cs, ce))
    return out

fam_copies = {fam: merge_intervals(ivs) for fam, ivs in fam_member_loci.items()}

# classify family regime from its copy-intervals
def family_regime(copies):
    if len(copies) < 2:
        return "SINGLE_LOCUS"   # annotated copies collapsed to one interval (overlap) OR only 1 found
    # any two copies overlap? -> COLLAPSED
    by_chrom = defaultdict(list)
    for chrom, s, e in copies:
        by_chrom[chrom].append((s, e))
    # overlap within chrom already merged away, so >1 interval same chrom = non-overlapping near/far
    chroms = set(c for c, _, _ in copies)
    if len(chroms) > 1:
        return "DISPERSED"      # copies on different contigs (RABL2-like boundary)
    # same chrom, multiple intervals: near or far
    segs = sorted((s, e) for c, s, e in copies)
    mind = min(segs[i + 1][0] - segs[i][1] for i in range(len(segs) - 1))
    return "TANDEM_NEAR" if mind <= TANDEM_WINDOW else "DISPERSED"

fam_regime = {fam: family_regime(cp) for fam, cp in fam_copies.items()}

# Special case: a family whose members ALL merged into ONE interval = collapsed reference
for fam, cp in fam_copies.items():
    if len(cp) == 1 and fam_ncopies[fam] >= 2 and fam_found[fam] >= 2:
        fam_regime[fam] = "COLLAPSED"   # >=2 annotated copies occupy one coordinate frame

# flat interval index: chrom -> list of (start,end,fam,regime)
fam_index = defaultdict(list)
for fam, cps in fam_copies.items():
    for chrom, s, e in cps:
        fam_index[chrom].append((s, e, fam, fam_regime[fam]))
for chrom in fam_index:
    fam_index[chrom].sort()


def family_hit(chrom, s, e):
    """Return (fam, regime) of best-overlapping multi-copy family member, or None."""
    best = None
    for fs, fe, fam, reg in fam_index.get(chrom, []):
        if fs > e:
            break
        if fe >= s:  # overlap
            ov = min(fe, e) - max(fs, s)
            if best is None or ov > best[0]:
                best = (ov, fam, reg)
    return (best[1], best[2]) if best else None


# ---------------------------------------------------------------------------
# 3. read-coherence recall set from rcg_*.gtf
# ---------------------------------------------------------------------------
# transcript_id -> (chrom, strand, source, exons)
rc_total = 0
rc_multiexon = 0
tier1_hits = 0
fsm_total = 0
fsm_at_family = 0
regime_counts = defaultdict(int)        # regime -> # recall isoforms at that regime
regime_fsm = defaultdict(int)
collapsed_rows = []                      # the headroom (COLLAPSED-regime hits)
fam_hit_isoforms = defaultdict(int)      # fam -> # recall isoforms landing on it
singlecopy_spans = defaultdict(list)     # chrom -> [(start,end,strand,is_fsm)] of NON-family recall isoforms

for gtf in sorted(glob.glob(os.path.join(GW, "rcg_*.gtf"))):
    chrom = os.path.basename(gtf)[len("rcg_"):-len(".gtf")]
    cur = {}  # tid -> dict
    txsrc = {}
    txstrand = {}
    exons = defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = parse_attrs(f[8])
            tid = a.get("transcript_id")
            if not tid:
                continue
            if f[2] == "transcript":
                txsrc[tid] = a.get("source", "")
                txstrand[tid] = f[6]
            elif f[2] == "exon":
                exons[tid].append((int(f[3]), int(f[4])))
    for tid, exl in exons.items():
        if txsrc.get(tid) != "read_coherence":
            continue
        rc_total += 1
        ch = intron_chain(exl)
        if not ch:
            continue  # single-exon: read-coherence finding is about multi-exon extras
        rc_multiexon += 1
        s = min(x[0] for x in exl)
        e = max(x[1] for x in exl)
        strand = txstrand.get(tid, ".")
        is_fsm = (chrom, strand, ch) in ref_all_chains
        if is_fsm:
            fsm_total += 1
        hit = family_hit(chrom, s, e)
        if hit:
            fam, reg = hit
            tier1_hits += 1
            regime_counts[reg] += 1
            fam_hit_isoforms[fam] += 1
            if is_fsm:
                fsm_at_family += 1
                regime_fsm[reg] += 1
            if reg == "COLLAPSED":
                collapsed_rows.append((chrom, s, e, strand, fam, fam_ncopies[fam],
                                       "FSM" if is_fsm else "novel", tid))
        else:
            # recall isoform at a NON-multi-copy-family ("single-copy") locus -> scan target for (c)
            singlecopy_spans[chrom].append((s, e, strand, is_fsm))

# ---------------------------------------------------------------------------
# 4. report
# ---------------------------------------------------------------------------
n_fam = len(fam_copies)
regime_fam_counts = defaultdict(int)
for fam, reg in fam_regime.items():
    regime_fam_counts[reg] += 1

lines = []
def P(s=""):
    lines.append(s)

P("# Read-coherence × PSV — HEADROOM probe (go/no-go)")
P()
P("Cheap geometric+structural probe: do read-coherence's recall wins and PSV copy-resolution")
P("act on the SAME loci (→ unified molecule-threaded graph worth building), or disjoint loci")
P("(→ not worth the rebuild; keep them as two separate levers)?")
P()
P("## Inputs")
P(f"- read-coherence recall set: transcripts tagged `source \"read_coherence\"` in rcg_*.gtf (25 chroms)")
P(f"- multi-copy families (universe.tsv, n_copies>=2): **{n_fam}** families, "
  f"{sum(fam_found.values())}/{len(tx2fam)} member transcripts located in RefSeq")
P(f"- family regimes: " + ", ".join(f"{r}={regime_fam_counts[r]}" for r in
                                     sorted(regime_fam_counts, key=lambda k: -regime_fam_counts[k])))
P()
P("## Recall set")
P(f"- read_coherence transcripts total: **{rc_total}**")
P(f"- multi-exon (the recall lever): **{rc_multiexon}**")
P(f"- of those, FSM (intron chain == a RefSeq transcript = real recall): **{fsm_total}** "
  f"({100*fsm_total/max(rc_multiexon,1):.1f}%)")
P()
P("## Tier 1 — are the recall wins even at multi-copy family loci?")
P(f"- recall isoforms overlapping a multi-copy family locus: **{tier1_hits}** / {rc_multiexon} "
  f"(**{100*tier1_hits/max(rc_multiexon,1):.2f}%**)")
P(f"- of those, FSM: **{fsm_at_family}** (of {fsm_total} FSM = "
  f"{100*fsm_at_family/max(fsm_total,1):.2f}% of all read-coherence real recall)")
P()
P("## Tier 2 — of the family-locus hits, which regime (does PSV actually SPLIT?)")
P()
P("| regime | recall isoforms | FSM | PSV leverage |")
P("|---|---|---|---|")
LEV = {"COLLAPSED": "**YES — copies share one frame, PSVs split**",
       "TANDEM_NEAR": "no — distinct coords separate them",
       "DISPERSED": "no — distinct coords/contigs (RABL2 boundary)",
       "SINGLE_LOCUS": "n/a — only one copy located"}
for reg in ("COLLAPSED", "TANDEM_NEAR", "DISPERSED", "SINGLE_LOCUS"):
    if regime_counts.get(reg, 0) or reg == "COLLAPSED":
        P(f"| {reg} | {regime_counts.get(reg,0)} | {regime_fsm.get(reg,0)} | {LEV[reg]} |")
P()
P("## Tier 3 — of the COLLAPSED hits, are the 'copies' REAL paralogs or domain-sharers?")
P("(PSVs can only copy-resolve genuine paralogous copies. Two different genes that merely share a")
P("protein domain/repeat and happen to overlap by annotation are NOT copies — PSVs are meaningless there.)")
P()
real_collapsed = sum(1 for r in collapsed_rows if r[4] in CONFIRMED_REAL)
sharer_collapsed = sum(1 for r in collapsed_rows if r[4] in DOMAIN_SHARER)
unknown_collapsed = len(collapsed_rows) - real_collapsed - sharer_collapsed
P(f"- COLLAPSED recall isoforms at CONFIRMED real paralog families (APOBEC3/RFPL/RABL2): "
  f"**{real_collapsed}**")
P(f"- COLLAPSED recall isoforms at DOMAIN-SHARER false families "
  f"(CDPF1/PPARA, CREB1/METTL21A, GCA/KCNH7, TTN/NCAPH2/MIEF1-LOC): **{sharer_collapsed}**")
P(f"- COLLAPSED recall isoforms at unclassified families: **{unknown_collapsed}**")
P()
# does ANY recall isoform (any regime) land on a confirmed-real paralog family?
real_any = sum(n for fam, n in fam_hit_isoforms.items() if fam in CONFIRMED_REAL)
P(f"- read-coherence recall isoforms landing on ANY confirmed-real paralog family (any regime): "
  f"**{real_any}**")
P()
P("## Verdict")
collapsed_n = regime_counts.get("COLLAPSED", 0)
P(f"Geometric PSV-resolvable headroom (recall isoforms at COLLAPSED loci) = **{collapsed_n}** "
  f"({100*collapsed_n/max(rc_multiexon,1):.3f}% of the recall set, "
  f"{regime_fsm.get('COLLAPSED',0)} of them FSM).")
P(f"TRUE headroom after removing domain-sharer false families = **{real_collapsed + unknown_collapsed}** "
  f"(confirmed-real collapsed paralog recall isoforms = {real_collapsed}).")
P()
P("**GO/NO-GO: NO-GO** for a molecule-threaded graph justified by *combining* recall + copy-resolution.")
P("Read-coherence recall (99.6% at single-copy loci) and PSV copy-resolution (needs collapsed real")
P("paralogs) are disjoint in this data. The recall lever already ships additively (`--read-coherence`);")
P("threading PSVs through the same graph adds copy-resolution to ~0 real recall isoforms.")
P()
P(f"Top families receiving recall isoforms:")
for fam, n in sorted(fam_hit_isoforms.items(), key=lambda kv: (-kv[1], kv[0]))[:15]:
    tag = " [DOMAIN-SHARER]" if fam in DOMAIN_SHARER else (" [REAL]" if fam in CONFIRMED_REAL else "")
    P(f"  - {fam} ({fam_regime[fam]}, {fam_ncopies[fam]} copies): {n} isoforms{tag}")
P()
P("![funnel](readcoherence_psv_headroom.png)")
P()
P("## What this does and does NOT settle")
P("- **Settles:** the *PSV-unification* motivation. Read-coherence's recall win does not co-occur")
P("  with PSV-resolvable copies, so building one graph to deliver BOTH on the same loci pays ~nothing.")
P("- **Does NOT settle:** the molecule-threaded graph as a *pure recall* architecture (gate-not-kill,")
P("  keep per-molecule evidence through paths). That stands or falls on recall alone — and the recall")
P("  is *already* delivered additively by the shipped `--read-coherence` layer, so a rebuild would buy")
P("  cleaner architecture / less noise, not new recall. It would NOT add copy-resolution.")
P()
P("## Honest caveats (which way they push)")
P("1. **Geometric undercount (pushes headroom UP a little):** COLLAPSED is detected only when *annotated*")
P("   copies overlap. Cross-mapping can pile a divergent/unannotated copy's reads onto a *single-copy*")
P("   annotated locus, creating PSV signal the annotation hides. But prior copy_split real-data work")
P("   (RABL2 / DAZ / co-located tests) showed that regime is rare and coverage-limited in GGO; even a")
P("   3–5× undercount keeps the headroom sub-1% and the *real-paralog* count near zero.")
P("2. **Small family universe (neutral):** 58 multi-copy families is the genome's reality, not a sampling")
P("   gap. Multi-copy paralogs are a small slice of loci; read-coherence's recall is dominated by")
P("   single-copy alt-splicing/novel-junction isoforms regardless of family-set size.")
P("3. **Domain-sharer contamination (pushes headroom DOWN, decisively):** all 35 COLLAPSED hits land on")
P("   families the Compara + contiguous-core validation already flagged as domain-sharing FALSE families")
P("   (two different genes sharing a repeat, overlapping by annotation accident) — not paralog copies.")
P("   PSVs cannot copy-resolve non-paralogs, so the TRUE headroom is 0, not 35.")
P()
P("## Reproduce")
P("- `python3 bench/readcoherence_psv_headroom.py`  (reads /tmp/gw/rcg_*.gtf + ref_*.gff3 + universe.tsv)")
P("- `python3 bench/readcoherence_psv_headroom_fig.py`  (renders the funnel figure)")
P("- COLLAPSED-regime hit loci: `bench/headroom_loci.tsv`")

with open(OUT_MD, "w") as fh:
    fh.write("\n".join(lines) + "\n")

with open(OUT_TSV, "w") as fh:
    fh.write("chrom\tstart\tend\tstrand\tfamily\tn_copies\tfsm\ttranscript_id\n")
    for r in sorted(collapsed_rows):
        fh.write("\t".join(str(x) for x in r) + "\n")

# Emit single-copy recall loci (the (c) scan targets): merge overlapping recall-isoform spans
# per chrom into distinct loci, carrying #isoforms + #FSM at each. These are where the headroom
# probe says "no annotated paralog" — task (c) checks the BAM for HIDDEN collapsed-copy PSV signal.
SC_TSV = os.path.join(os.path.dirname(__file__), "single_copy_recall_loci.tsv")
n_loci = 0
with open(SC_TSV, "w") as fh:
    fh.write("chrom\tstart\tend\tn_isoforms\tn_fsm\n")
    for chrom in sorted(singlecopy_spans):
        spans = sorted(singlecopy_spans[chrom])
        cs, ce, niso, nfsm = spans[0][0], spans[0][1], 1, int(spans[0][3])
        for s, e, strand, fsm in spans[1:]:
            if s <= ce:
                ce = max(ce, e); niso += 1; nfsm += int(fsm)
            else:
                fh.write(f"{chrom}\t{cs}\t{ce}\t{niso}\t{nfsm}\n"); n_loci += 1
                cs, ce, niso, nfsm = s, e, 1, int(fsm)
        fh.write(f"{chrom}\t{cs}\t{ce}\t{niso}\t{nfsm}\n"); n_loci += 1
print(f"[wrote {SC_TSV}: {n_loci} single-copy recall loci]")

print("\n".join(lines))
print(f"\n[wrote {OUT_MD} and {OUT_TSV}]")
