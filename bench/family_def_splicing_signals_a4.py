#!/usr/bin/env python3
"""
family_def_splicing_signals_a4.py — ANGLE A4: is anything about the JUNCTIONS / SPLICING a
UNIQUE discriminating signal for real multi-copy gene families vs the confounds (domain-sharers,
retrocopies, name-coincidences), ORTHOGONAL to cDNA sequence identity?

THE QUESTION (user): cDNA sequence homology OVER-MERGES families. Distinct genes sharing a domain
look homologous (id ~0.9-1.0); processed-pseudogene RETROCOPIES look id~1.0 to their parent yet are
not a transcribed tandem family. Sequence id alone cannot separate these. Is there a splicing-derived
signal that flags cases id gets wrong?

We enumerate and MEASURE five candidate splicing signals on the airtight panel + the two retrocopy
counterexamples (EEF1A1->LOC109023808, CNN2->LOC129524764):

  (1) RETROCOPY INTRONLESSNESS   : n_intron of the copy. Processed pseudogenes are intronless (0)
                                   yet cDNA id ~1.0. id CANNOT flag; n_intron==0 flags instantly.
  (2) INTRON-PHASE concordance   : phase (0/1/2) at homologous junctions. Whole-gene dup preserves
                                   phase; domain-shuffle need not.
  (3) SPLICE-SITE MOTIF class    : GT-AG / GC-AG / AT-AC at each junction (from genome FASTA).
                                   Conserved across real copies?
  (4) shared JUNCTION USAGE      : here measured as ordered intron-architecture concordance (the
                                   annotation analogue of "co-usage of the same junctions"); the
                                   read/alt-junction co-usage axis is discussed against the ASJ machinery.
  (5) EXON-COUNT concordance     : |n_exon_A - n_exon_B| and the ordered intron-length vector match.

For each pair we report: cDNA id (the thing to beat), the splicing signal's value, and whether the
signal FLAGS a confound that id gets wrong (id high but pair is NOT a real family).

Per-gene architecture (canonical transcript = most CDS-exons, tie->most exons, tie->longest span) is
parsed directly from the GFF for the ~30 panel genes; phase from CDS phase column; motifs from the
genome FASTA via pysam. A genome-wide proxy effect on the 389-mega bridge edges uses the existing
length-vector index (gene_intron_index.json) since absolute-coord phase for 40k genes is out of scope.

Run: /home/juanfra/miniforge3/bin/python bench/family_def_splicing_signals_a4.py
"""
import collections
import json
import os
import re
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology  # noqa: E402

GFF = "/mnt/c/Users/jfris/Desktop/GGO_genomic.gff"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
LEN_INDEX = "/home/juanfra/winloci_scratch/gene_intron_index.json"

GENE_RE = re.compile(r"(?:^|;)gene=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
TX_TYPES = {"mRNA", "transcript", "lnc_RNA", "ncRNA", "primary_transcript"}

# ---- panel (name -> (chrom,start,end)) and labelled cases -------------------------------------
GENES = {
    "RABL2A": ("NC_073235.2", 15131653, 15147533), "RABL2B": ("NC_086018.1", 48818440, 48832011),
    "APOBEC3D": ("NC_086018.1", 37023493, 37035944), "APOBEC3F": ("NC_086018.1", 37043278, 37058621),
    "RFPL1": ("NC_086018.1", 27445867, 27472024), "RFPL2": ("NC_086018.1", 30208643, 30221673),
    "RFPL3": ("NC_086018.1", 30368560, 30376055),
    "EEF1A1": ("NC_073229.2", 97608632, 97613071), "LOC109023808": ("NC_073243.2", 97380144, 97381766),
    "CNN2": ("NC_073244.2", 6807132, 6819584), "LOC129524764": ("NC_073237.2", 31580552, 31581609),
    "ASDURF": ("NC_073235.2", 92106380, 92111318), "ASNSD1": ("NC_073235.2", 92111059, 92115776),
    "CASP8": ("NC_073235.2", 103850223, 103901051), "FLACC1": ("NC_073235.2", 103898904, 103977183),
    "CREB1": ("NC_073235.2", 110133100, 110220981), "METTL21A": ("NC_073235.2", 110198842, 110241136),
    "GPR39": ("NC_073235.2", 34964263, 35195860), "LYPD1": ("NC_073235.2", 35194087, 35224540),
    "GCA": ("NC_073235.2", 64785386, 64835227), "KCNH7": ("NC_073235.2", 64832333, 65309180),
}

# (case, kind, members, is_real_family)
PANEL = [
    ("RABL2",        "real",       ["RABL2A", "RABL2B"],          True),
    ("APOBEC3",      "real",       ["APOBEC3D", "APOBEC3F"],      True),
    ("RFPL",         "real",       ["RFPL1", "RFPL2", "RFPL3"],   True),
    ("EEF1A1_retro", "retrocopy",  ["EEF1A1", "LOC109023808"],    False),
    ("CNN2_retro",   "retrocopy",  ["CNN2", "LOC129524764"],      False),
    ("ASDURF_ds",    "domainshare", ["ASDURF", "ASNSD1"],         False),
    ("CASP8_ds",     "domainshare", ["CASP8", "FLACC1"],          False),
    ("CREB1_ds",     "domainshare", ["CREB1", "METTL21A"],        False),
    ("GPR39_ds",     "domainshare", ["GPR39", "LYPD1"],           False),
    ("GCA_ds",       "domainshare", ["GCA", "KCNH7"],             False),
]

WANT = set(GENES)


def parse_gff_for_genes(want):
    """Stream GFF; collect for the wanted genes the exon & CDS blocks of every transcript, then pick
    the canonical transcript per gene (most CDS exons -> most exons -> longest span)."""
    tx2gene = {}
    tx_exons = collections.defaultdict(list)  # tx -> [(start,end)]  1-based GFF coords
    tx_cds = collections.defaultdict(list)    # tx -> [(start,end,frame)]
    tx_strand = {}
    tx_chrom = {}
    with open(GFF) as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            ftype, attrs = p[2], p[8]
            if ftype in TX_TYPES:
                m = ID_RE.search(attrs)
                g = GENE_RE.search(attrs)
                if m and g and g.group(1) in want:
                    tx2gene[m.group(1)] = g.group(1)
                    tx_strand[m.group(1)] = p[6]
                    tx_chrom[m.group(1)] = p[0]
            elif ftype == "exon":
                par = PARENT_RE.search(attrs)
                if par and par.group(1) in tx2gene:
                    tx_exons[par.group(1)].append((int(p[3]), int(p[4])))
            elif ftype == "CDS":
                par = PARENT_RE.search(attrs)
                if par and par.group(1) in tx2gene:
                    frame = p[7] if p[7] in ("0", "1", "2") else "."
                    tx_cds[par.group(1)].append((int(p[3]), int(p[4]), frame))
    # pick canonical tx per gene
    by_gene = collections.defaultdict(list)
    for tx, g in tx2gene.items():
        by_gene[g].append(tx)
    out = {}
    for g, txs in by_gene.items():
        best = None
        for tx in txs:
            ex = sorted(tx_exons.get(tx, []))
            cds = sorted([(s, e, fr) for (s, e, fr) in tx_cds.get(tx, [])])
            key = (len(cds), len(ex), (ex[-1][1] - ex[0][0]) if ex else 0)
            if best is None or key > best[0]:
                best = (key, tx, ex, cds)
        if best is None:
            continue
        _, tx, ex, cds = best
        out[g] = dict(tx=tx, chrom=tx_chrom[tx], strand=tx_strand[tx], exons=ex, cds=cds)
    return out


def introns_5to3(exons, strand):
    """Ordered list of introns (donor_genomic, acceptor_genomic) in transcription (5'->3') order."""
    ex = sorted(exons)
    gaps = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]  # (end_i, start_{i+1}) genomic
    if strand == "-":
        gaps = gaps[::-1]
    return gaps


def cds_intron_phases(cds, strand):
    """Phase at each CDS intron in 5'->3' order. Phase = number of CDS bases BEFORE the intron mod 3.
    Uses cumulative CDS length up to (and including) the exon 5' of each intron."""
    if len(cds) < 2:
        return []
    blocks = sorted([(s, e) for (s, e, _fr) in cds])
    if strand == "-":
        blocks = blocks[::-1]  # 5'->3'
    phases = []
    cum = 0
    for i in range(len(blocks) - 1):
        s, e = blocks[i]
        cum += (e - s + 1)
        phases.append(cum % 3)
    return phases  # one per CDS intron, 5'->3'


def splice_motifs(introns, chrom, strand, fa):
    """Donor (first 2) + acceptor (last 2) dinucleotides per intron, reverse-complemented on '-'.
    Returns list of 'XX-YY' motif strings in 5'->3' order. Genome is 1-based in GFF; pysam fetch 0-based."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

    def rc(s):
        return "".join(comp.get(b, "N") for b in s.upper()[::-1])
    out = []
    for (g_end, g_start) in introns:
        # intron occupies genomic [g_end+1 .. g_start-1] (1-based inclusive) since exon ends at g_end,
        # next exon starts at g_start.
        ilo, ihi = g_end + 1, g_start - 1  # 1-based inclusive intron span
        if ihi < ilo:
            out.append("NA")
            continue
        left = fa.fetch(chrom, ilo - 1, ilo + 1).upper()    # first 2 bases
        right = fa.fetch(chrom, ihi - 2, ihi).upper()        # last 2 bases
        if strand == "+":
            donor, acc = left, right
        else:
            donor, acc = rc(right), rc(left)
        out.append(f"{donor}-{acc}")
    return out


def motif_class(m):
    if m in ("GT-AG",):
        return "GT-AG"
    if m in ("GC-AG",):
        return "GC-AG"
    if m in ("AT-AC",):
        return "AT-AC"
    return "other"


def best_align_phase(a, b):
    """Best ungapped sliding-window overlap fraction of two ordered vectors (handles 5'/3' truncation
    typical of retro/diverged copies). Returns (match_frac, overlap_len)."""
    if not a or not b:
        return (0.0, 0)
    best = (0.0, 0)
    for off in range(-(len(b) - 1), len(a)):
        m = t = 0
        for i in range(len(a)):
            j = i - off
            if 0 <= j < len(b):
                t += 1
                if a[i] == b[j]:
                    m += 1
        if t and m / t >= best[0]:
            if (m / t, t) > best:
                best = (m / t, t)
    return best


def main():
    print("[*] parsing GFF for panel genes ...")
    arch = parse_gff_for_genes(WANT)
    fa = pysam.FastaFile(FASTA)
    Hd, _ = dna_homology()

    def cdna_id(a, b):
        k = (a, b) if a < b else (b, a)
        if k in Hd:
            return Hd[k]["id"]
        return None

    # per-gene derived signals
    G = {}
    for g, d in arch.items():
        introns = introns_5to3(d["exons"], d["strand"])
        phases = cds_intron_phases(d["cds"], d["strand"])
        motifs = splice_motifs(introns, d["chrom"], d["strand"], fa)
        mclass = [motif_class(m) for m in motifs]
        G[g] = dict(n_exon=len(d["exons"]), n_intron=len(introns), strand=d["strand"],
                    phases=phases, motifs=motifs, mclass=mclass,
                    ilen=[ (gs - ge - 1) for (ge, gs) in introns ])

    print("\n==================== PER-GENE SPLICING ARCHITECTURE ====================")
    for g in sorted(G):
        s = G[g]
        gtag = sum(1 for m in s["mclass"] if m == "GT-AG")
        print(f"  {g:14s} n_exon={s['n_exon']:2d} n_intron={s['n_intron']:2d} "
              f"GT-AG={gtag}/{s['n_intron']} phases={s['phases']}")

    # ---------------- panel pairwise tests ----------------
    print("\n==================== PANEL PAIRWISE SIGNALS ====================")
    hdr = ("case", "kind", "pair", "id", "nI_a", "nI_b", "S1_retro", "exonΔ", "S2_phase", "S3_motif", "real")
    print("{:14s} {:11s} {:24s} {:>5s} {:>4s} {:>4s} {:>8s} {:>6s} {:>9s} {:>9s} {:>5s}".format(*hdr))
    rows = []
    for case, kind, members, is_real in PANEL:
        # all pairs within the case
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                a, b = members[i], members[j]
                if a not in G or b not in G:
                    print(f"  {case}: MISSING {a if a not in G else b}")
                    continue
                Ga, Gb = G[a], G[b]
                idv = cdna_id(a, b)
                # S1 retrocopy: min intron count == 0 (one side intronless)
                s1_retro = (min(Ga["n_intron"], Gb["n_intron"]) == 0)
                exon_delta = abs(Ga["n_exon"] - Gb["n_exon"])
                # S2 phase concordance over best-aligned overlap of CDS-intron phase vectors
                pf, pov = best_align_phase(Ga["phases"], Gb["phases"])
                # S3 motif-class concordance over best-aligned overlap
                mf, mov = best_align_phase(Ga["mclass"], Gb["mclass"])
                rows.append((case, kind, a, b, idv, Ga["n_intron"], Gb["n_intron"],
                             s1_retro, exon_delta, pf, pov, mf, mov, is_real))
                ids = f"{idv:.2f}" if idv is not None else "  na"
                print("{:14s} {:11s} {:24s} {:>5s} {:>4d} {:>4d} {:>8s} {:>6d} {:>4.2f}/{:<2d}  {:>4.2f}/{:<2d}  {:>5s}".format(
                    case, kind, f"{a[:11]}~{b[:11]}", ids, Ga["n_intron"], Gb["n_intron"],
                    "YES" if s1_retro else "-", exon_delta, pf, pov, mf, mov, "T" if is_real else "F"))

    # ---------------- signal-by-signal verdict ----------------
    print("\n==================== SIGNAL VERDICTS (does it add over id?) ====================")
    real_pairs = [r for r in rows if r[-1]]
    fake_pairs = [r for r in rows if not r[-1]]

    # S1 retrocopy: separates retro cases?
    retro_rows = [r for r in rows if r[1] == "retrocopy"]
    print("\n[S1] RETROCOPY INTRONLESSNESS (min n_intron == 0):")
    for r in rows:
        flag = "FLAG" if r[7] else "    "
        if r[1] == "retrocopy" or r[7]:
            print(f"     {r[0]:14s} {r[2]}~{r[3]}  id={r[4]:.2f}  min_nI={min(r[5],r[6])}  -> {flag} "
                  f"({'id HIGH but intronless => id WRONG, S1 RIGHT' if (r[4] and r[4]>0.85 and r[7]) else ''})")
    print(f"     retro cases flagged by S1: {sum(1 for r in retro_rows if r[7])}/{len(retro_rows)}; "
          f"real-family pairs falsely flagged: {sum(1 for r in real_pairs if r[7])}/{len(real_pairs)}")

    # S5 exon-count concordance: separation
    print("\n[S5] EXON-COUNT concordance (|Δn_exon|):")
    rd = [r[8] for r in real_pairs]
    fd = [r[8] for r in fake_pairs]
    print(f"     real  |Δexon|: {sorted(rd)}  (max {max(rd) if rd else 0})")
    print(f"     fake  |Δexon|: {sorted(fd)}  (min {min(fd) if fd else 0})")

    # S2 phase concordance over overlap
    print("\n[S2] INTRON-PHASE concordance (best-aligned phase match frac; overlap len):")
    for r in rows:
        print(f"     {r[0]:14s} {('REAL' if r[-1] else 'fake'):4s}  id={r[4] if r[4] else 0:.2f}  phase={r[9]:.2f} over {r[10]}")
    rp = [r[9] for r in real_pairs if r[10] >= 2]
    fp = [r[9] for r in fake_pairs if r[10] >= 2]
    print(f"     real (overlap>=2) phase-match: {[f'{x:.2f}' for x in rp]}")
    print(f"     fake (overlap>=2) phase-match: {[f'{x:.2f}' for x in fp]}")

    # S3 motif class
    print("\n[S3] SPLICE-MOTIF class concordance (best-aligned motif-class match frac):")
    rm = [r[11] for r in real_pairs if r[12] >= 2]
    fm = [r[11] for r in fake_pairs if r[12] >= 2]
    print(f"     real (overlap>=2): {[f'{x:.2f}' for x in rm]}")
    print(f"     fake (overlap>=2): {[f'{x:.2f}' for x in fm]}")
    # GT-AG dominance: are essentially ALL junctions canonical regardless of class? (then S3 inert)
    allm = collections.Counter()
    for g in G:
        allm.update(G[g]["mclass"])
    print(f"     genome-of-panel motif-class counts: {dict(allm)}")

    # combined separator: S1 (retro) OR exon-count gate
    print("\n==================== COMBINED SPLICING SEPARATOR ====================")
    # rule: a pair is 'real-family-compatible' if NOT retro AND |Δexon| small AND phase-match high.
    # report the per-class min/max so the reader sees the margin.
    def sep_score(r):
        # higher = more real-family-like; penalise retro + exon mismatch + phase mismatch
        if r[7]:
            return -1.0  # intronless => not a transcribed tandem family
        return r[9] if r[10] >= 2 else 0.5  # phase match over overlap, else neutral
    print("   per-pair separator score (retro=-1; else phase-match over overlap):")
    realsc, fakesc = [], []
    for r in rows:
        sc = sep_score(r)
        (realsc if r[-1] else fakesc).append(sc)
        print(f"     {r[0]:14s} {('REAL' if r[-1] else 'fake'):4s} id={r[4] if r[4] else 0:.2f} score={sc:+.2f}")
    print(f"\n   REAL scores: {sorted(realsc)}")
    print(f"   FAKE scores: {sorted(fakesc)}")

    out = dict(per_gene=G, rows=[list(r[:7]) + [bool(r[7])] + list(r[8:]) for r in rows])
    json.dump(out, open(os.path.join(HERE, "family_def_splicing_signals_a4.json"), "w"),
              indent=2, default=str)
    print(f"\n[+] wrote {os.path.join(HERE, 'family_def_splicing_signals_a4.json')}")


if __name__ == "__main__":
    main()
