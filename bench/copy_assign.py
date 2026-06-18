#!/usr/bin/env python3
"""Copy ASSIGNMENT (resolution): assign each read -- especially hard multimappers that minimap2 leaves
at MAPQ 0 -- to a specific paralog COPY, using paralog-sequence variants (PSVs) and, where copies differ
in splicing, copy-specific JUNCTIONS (the copy-level analogue of allele-specific junctions).

This is the Canzar "resolve the ambiguity, pick one mapping" step, on top of the de-novo family layer
(denovo_families_split.tsv): the family tells us the copy SET; this assigns the reads to copies.

METHOD (realistic, BAM-driven -- the version that generalises to real data, unlike the oracle in
build_sim5x.py which knows the designed PSV indices/alleles):
  1. PSV discovery from the copy SEQUENCES (all-pairs alignment anchored on copy[0]); each copy gets an
     allele vector A_c over the PSV columns (in transcription/seq space, strand-correct).
  2. copy-specific JUNCTIONS: introns observed (from reads) in some copies but not others -> extra
     discriminating columns (a junction is "present"/"absent" per copy).
  3. for each read: walk its CIGAR in the frame of the copy it aligned to, read its base at each PSV
     position + its intron set -> an observed feature vector.
  4. assign by likelihood: logL(read|copy c) = sum over spanned PSV columns of log P(obs base | A_c, e)
     + sum over junctions of log P(obs junction | copy c junction set). Assign = argmax_c logL.
  5. IDENTIFIABILITY gate (the theorem): a read is RESOLVABLE iff it spans >= 1 feature where the
     candidate copies differ. Reads spanning 0 such features are genuinely TIED ("5 equally good places").
     Among resolvable reads, require a log-likelihood-ratio margin over the runner-up to call ASSIGNED.

Deterministic. Run with /home/juanfra/miniforge3/bin/python.  Usage: python copy_assign.py sim5x
"""
import math
import os
import sys
from collections import defaultdict

import pysam
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

ERR = 0.003          # per-base error assumed in the likelihood (HiFi)
MARGIN = 2.0         # min log-likelihood-ratio over runner-up to call ASSIGNED (else AMBIGUOUS)
SIM5X = "/home/juanfra/winloci_scratch/sim5x"

_COMP = str.maketrans("ACGTN", "TGCAN")


def rc_base(b):
    return b.translate(_COMP)


def pairwise_align(a, b):
    al = PairwiseAligner()
    al.mode = "global"
    al.match_score = 2
    al.mismatch_score = -1
    al.open_gap_score = -5
    al.extend_gap_score = -1
    return al.align(a, b)[0]


def read_copy_seqs(fa, copies):
    """copies: list of (name, contig, start1, end1, strand). Return {name: transcription-strand seq}."""
    seqs = {}
    for name, contig, s1, e1, strand in copies:
        s = fa.fetch(contig, s1 - 1, e1).upper()
        if strand == "-":
            s = str(Seq(s).reverse_complement())
        seqs[name] = s
    return seqs


def seqoff_to_genomic(copy, off):
    """transcription-strand seq offset -> 0-based genomic coord on the forward genome."""
    name, contig, s1, e1, strand = copy
    if strand == "+":
        return (s1 - 1) + off
    return (s1 - 1) + ((e1 - s1) - off)


def discover_psvs(copies, seqs, off2gen):
    """All-pairs vs copies[0]. Returns columns: list of {copyname: (gpos0, seq_base) or None}.
    off2gen[name](off) maps a copy's seq offset -> genomic 0-based coord (genomic span for sim5x; spliced
    exon map for real transcripts). A column = a ref offset where >=1 copy DIFFERS. At each column EVERY
    copy gets its aligned base (which may EQUAL the ref allele -- crucial: a copy that matches the ref must
    inherit that allele, not be left missing, else the likelihood would favour copies it shares no columns
    with). None = a genuine alignment gap only. seq_base is transcription/seq space (comparable across copies)."""
    ref = copies[0]
    rseq = seqs[ref[0]]
    aligned = {}        # other -> {ref_off: copy_off} over ALL ungapped aligned positions
    diff_off = set()    # ref offsets where >=1 copy differs (the PSV columns)
    for other in copies[1:]:
        oseq = seqs[other[0]]
        aln = pairwise_align(rseq, oseq)
        amap = {}
        for (a0, a1), (b0, b1) in zip(*aln.aligned):
            for i in range(a1 - a0):
                ai, bi = a0 + i, b0 + i
                amap[ai] = bi
                if rseq[ai] in "ACGT" and oseq[bi] in "ACGT" and rseq[ai] != oseq[bi]:
                    diff_off.add(ai)
        aligned[other[0]] = amap
    cols = []
    for ro in sorted(diff_off):
        rec = {ref[0]: (off2gen[ref[0]](ro), rseq[ro])}
        for other in copies[1:]:
            oo = aligned[other[0]].get(ro)
            if oo is None:
                rec[other[0]] = None                      # genuine gap: copy not aligned here
            else:
                rec[other[0]] = (off2gen[other[0]](oo), seqs[other[0]][oo])
        cols.append(rec)
    return cols


def read_bases_and_introns(aln, want_gpos):
    """Walk one alignment; return (base_at{gpos0->fwd-genome base}, introns set[(d0,a0)])."""
    ref_pos = aln.reference_start
    qpos = 0
    seq = aln.query_sequence or ""
    base_at = {}
    introns = set()
    for op, length in aln.cigartuples or []:
        if op in (0, 7, 8):
            for k in range(length):
                g = ref_pos + k
                if g in want_gpos:
                    base_at[g] = seq[qpos + k] if (qpos + k) < len(seq) else None
            ref_pos += length
            qpos += length
        elif op == 1:
            qpos += length
        elif op == 2:
            ref_pos += length
        elif op == 3:
            introns.add((ref_pos, ref_pos + length))
            ref_pos += length
        elif op == 4:
            qpos += length
        # 5 (hard clip), 6 (pad): no-op
    return base_at, introns


def assign_read(obs, copy_vecs, read_bounds=None, copy_bounds=None, e=ERR, jw=5.0, btol=4):
    """obs: {col_idx -> observed seq-space base} for spanned PSV columns; copy_vecs: {copy -> {col -> base}}.
    Optional JUNCTION evidence: read_bounds = set of the read's intron-boundary spliced offsets;
    copy_bounds = {copy -> set of that copy's intron-boundary offsets}. A read boundary matching a copy's
    boundary set (within btol) is +jw for that copy, a non-match is -jw (a read using a junction the copy
    lacks argues against it). Returns (best, margin, n_span, n_decisive, logL). n_decisive counts BOTH
    PSV columns and junction boundaries where the candidate copies differ (the identifiability gate)."""
    names = list(copy_vecs)
    lp_match, lp_mis = math.log(1 - e), math.log(e / 3)
    spanned = [c for c in obs if obs[c] is not None]
    n_decisive = 0
    logL = {n: 0.0 for n in names}
    # --- PSV term ---
    for c in spanned:
        alleles = {copy_vecs[n].get(c) for n in names if copy_vecs[n].get(c) is not None}
        if len(alleles) > 1:
            n_decisive += 1
        for n in names:
            a = copy_vecs[n].get(c)
            if a is not None:
                logL[n] += lp_match if obs[c] == a else lp_mis
    # --- junction (copy-specific splicing) term ---
    if read_bounds and copy_bounds:
        def has(b, s):
            return any(abs(b - x) <= btol for x in s)
        for b in read_bounds:
            present = [n for n in names if has(b, copy_bounds.get(n, ()))]
            if 0 < len(present) < len(names):     # a boundary that some copies have and others lack
                n_decisive += 1
            for n in names:
                logL[n] += jw if has(b, copy_bounds.get(n, ())) else -jw
    order = sorted(names, key=lambda n: logL[n], reverse=True)
    best = order[0]
    margin = logL[order[0]] - logL[order[1]] if len(order) > 1 else float("inf")
    return best, margin, len(spanned), n_decisive, logL


# --------------------------------------------------------------------------- sim5x validation
def sim5x_copies(fa, K):
    contig = f"SIM5X_K{K}"
    clen = fa.get_reference_length(contig)
    len_g = (clen - 4 * 2000) // 5     # contig = 5*len_g + 4*SPACER(2000)
    unit = len_g + 2000
    copies = [(f"copy{c}", contig, c * unit + 1, c * unit + len_g, "+") for c in range(5)]
    return copies, unit


def validate_sim5x():
    print(f"{'K':>2} {'reads':>6} {'PSVcols':>7} {'resolvable%':>11} {'acc|assigned':>12} "
          f"{'acc|argmax':>10} {'tied%':>6}")
    fa_paths = {}
    rows = []
    for K in (0, 1, 2, 4, 8):
        ref = f"{SIM5X}/sim5x_K{K}.ref.fa"
        bam = f"{SIM5X}/sim5x_K{K}.bam"
        if not os.path.exists(bam):
            continue
        fa = pysam.FastaFile(ref)
        copies, unit = sim5x_copies(fa, K)
        seqs = read_copy_seqs(fa, copies)
        off2gen = {c[0]: (lambda cp: (lambda off: seqoff_to_genomic(cp, off)))(c) for c in copies}
        cols = discover_psvs(copies, seqs, off2gen)
        # per-copy allele vector in column space
        copy_vecs = {c[0]: {} for c in copies}
        # per-copy genomic PSV positions (to read a read aligned to that copy)
        copy_gpos = {c[0]: {} for c in copies}
        for j, rec in enumerate(cols):
            for name, v in rec.items():
                if v is not None:
                    gpos, base = v
                    copy_vecs[name][j] = base
                    copy_gpos[name][j] = gpos
        b = pysam.AlignmentFile(bam, "rb")
        n = correct_assigned = n_assigned = correct_argmax = n_argmax = n_resolvable = n_tied = 0
        for aln in b.fetch():
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                continue
            n += 1
            true_c = int(aln.query_name.split("_c")[1].split("_")[0])
            mapped_c = aln.reference_start // unit
            mc = f"copy{mapped_c}"
            want = set(copy_gpos[mc].values())
            base_at, _ = read_bases_and_introns(aln, want)
            # observed allele vector in column space (read aligned to mapped copy mc)
            gpos2col = {g: j for j, g in copy_gpos[mc].items()}
            obs = {gpos2col[g]: base_at[g] for g in base_at if g in gpos2col}
            best, margin, n_span, n_dec, _ = assign_read(obs, copy_vecs)
            best_c = int(best.replace("copy", ""))
            if n_dec >= 1:
                n_resolvable += 1
                n_argmax += 1
                correct_argmax += (best_c == true_c)
                if margin >= MARGIN:
                    n_assigned += 1
                    correct_assigned += (best_c == true_c)
            else:
                n_tied += 1
        ra = (100 * n_resolvable / n) if n else 0
        aa = (correct_assigned / n_assigned) if n_assigned else float("nan")
        am = (correct_argmax / n_argmax) if n_argmax else float("nan")
        tp = (100 * n_tied / n) if n else 0
        print(f"{K:>2} {n:>6} {len(cols):>7} {ra:>10.1f}% {aa:>12.3f} {am:>10.3f} {tp:>5.1f}%")
        rows.append((K, n, len(cols), ra, aa, am, tp))
    print("\nreading: resolvable% = reads spanning >=1 PSV where copies differ (identifiability gate);")
    print("acc|assigned = accuracy among reads passing the likelihood-margin; K=0 has no PSVs -> 0% resolvable (all tied).")
    return rows


# --------------------------------------------------------------------------- real GGO families
GGO_FA = "/home/juanfra/winloci_scratch/GGO.fasta"
GGO_BAM = "/home/juanfra/winloci_scratch/GGO.bam"
SPLICED_FA = "/home/juanfra/winloci_scratch/denovo_transcripts.fa"
SKEL_F = "/home/juanfra/winloci_scratch/denovo_skeletons.tsv"
SPLIT = os.path.join(os.path.dirname(__file__), "denovo_families_split.tsv")
META_F = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"


def load_meta():
    meta = {}
    for line in open(META_F):
        if line.startswith("id"):
            continue
        tid, c, s, e, st, ne, nr = line.rstrip("\n").split("\t")
        meta[tid] = (c, int(s), int(e), st, int(nr))
    return meta


def load_skel():
    skel = {}
    for line in open(SKEL_F):
        if line.startswith("chrom"):
            continue
        c, s, e, ne, nr, intr = line.rstrip("\n").split("\t")
        skel[(c, int(s), int(e))] = ([tuple(map(int, x.split("-"))) for x in intr.split(";")] if intr else [])
    return skel


def copy_idx2gen(tid, meta, skel):
    """Map a spliced-seq offset -> genomic 0-based coord, in transcription order (exon-only)."""
    c, s, e, strand, nr = meta[tid]
    introns = skel.get((c, s, e), [])
    exons = []
    prev = s
    for d, a in introns:
        exons.append((prev, d))
        prev = a
    exons.append((prev, e))
    gpos = []
    if strand == "-":
        for xs, xe in reversed(exons):
            gpos.extend(range(xe - 1, xs - 1, -1))
    else:
        for xs, xe in exons:
            gpos.extend(range(xs, xe))
    return gpos


RESCUED = os.path.join(os.path.dirname(__file__), "denovo_rescued_copies.tsv")
RESCUED_FA = os.path.join(os.path.dirname(__file__), "denovo_rescued_copies.fa")


def load_rescued(meta, skel, spliced):
    """Merge family-aware rescued copies (denovo_rescued_copies.*) into meta/skel/spliced; return
    {family_id: [rescued tids]} so they can be added to their family rosters."""
    if not os.path.exists(RESCUED):
        return {}
    seqs = {}
    fa = pysam.FastaFile(RESCUED_FA)
    for r in fa.references:
        seqs[r] = fa.fetch(r)
    by_fam = defaultdict(list)
    for line in open(RESCUED):
        if line.startswith("family_id"):
            continue
        fid, tid, c, s, e, strand, nex, sup, cr, intr = line.rstrip("\n").split("\t")
        s, e = int(s), int(e)
        meta[tid] = (c, s, e, strand, int(sup))
        skel[(c, s, e)] = [tuple(map(int, x.split("-"))) for x in intr.split(";")] if intr else []
        spliced[tid] = seqs.get(tid, "")
        by_fam[fid].append(tid)
    return by_fam


def colocated_families(meta, win=5_000_000, min_copies=3, rescued_fam=None):
    """Return [(family_id, dom_gene, [copies...])] for same-chrom tandem clusters (the hard regime)."""
    rescued_fam = rescued_fam or {}
    out = []
    for line in open(SPLIT):
        if line.startswith("family_id"):
            continue
        f = line.rstrip("\n").split("\t")
        fid, n, nchr, dom, conc, dens, core, arts, klass, members = f
        if klass != "family":
            continue
        mem = members.split(",") + rescued_fam.get(fid, [])
        bychr = defaultdict(list)
        for t in mem:
            if t in meta:
                bychr[meta[t][0]].append(t)
        for c, ts in bychr.items():
            if len(ts) < min_copies:
                continue
            ts.sort(key=lambda t: meta[t][1])
            if meta[ts[-1]][1] - meta[ts[0]][1] <= win:
                copies = [(t, meta[t][0], meta[t][1] + 1, meta[t][2], meta[t][3]) for t in ts]
                out.append((fid, dom, copies))
    return out


def copy_boundaries(tid, meta, skel):
    """Set of intron-boundary offsets in transcription-spliced space (cumulative exon lengths)."""
    c, s, e, strand, nr = meta[tid]
    introns = skel.get((c, s, e), [])
    exons = []
    prev = s
    for d, a in introns:
        exons.append((prev, d))
        prev = a
    exons.append((prev, e))
    lens = [xe - xs for xs, xe in (reversed(exons) if strand == "-" else exons)]
    bnd = set()
    cum = 0
    for L in lens[:-1]:
        cum += L
        bnd.add(cum)
    return bnd


def assign_family(spliced, bam, copies, meta, skel):
    """PSV + copy-specific-junction copy-assignment on one co-located family, using SPLICED (exon-only)
    copy sequences + a spliced<->genomic exon map. Returns a stats dict (PSV-only and PSV+junction)."""
    seqs, idx2gen, gen2off, strand, cbounds, ok = {}, {}, {}, {}, {}, []
    for c in copies:
        tid = c[0]
        s = spliced.get(tid)
        g = copy_idx2gen(tid, meta, skel)
        if s and len(s) == len(g):
            seqs[tid] = s
            idx2gen[tid] = g
            gen2off[tid] = {gp: i for i, gp in enumerate(g)}
            strand[tid] = c[4]
            cbounds[tid] = copy_boundaries(tid, meta, skel)
            ok.append(c)
    if len(ok) < 2:
        return dict(n=0, mq0=0, psv_cols=0, resolvable=0, assigned=0, uniq=0, uniq_agree=0,
                    resolvable_j=0, assigned_j=0, junction_only=0, uniq_j=0, uniq_agree_j=0)
    copies = ok
    off2gen = {c[0]: (lambda nm: (lambda off: idx2gen[nm][off]))(c[0]) for c in copies}
    cols = discover_psvs(copies, seqs, off2gen)
    copy_vecs = {c[0]: {} for c in copies}
    copy_gpos = {c[0]: {} for c in copies}
    for j, rec in enumerate(cols):
        for name, v in rec.items():
            if v is not None:
                copy_vecs[name][j] = v[1]
                copy_gpos[name][j] = v[0]
    spans = {c[0]: (c[1], c[2] - 1, c[3]) for c in copies}
    contig = copies[0][1]
    lo = min(c[2] - 1 for c in copies)
    hi = max(c[3] for c in copies)
    n = uniq = uniq_agree = resolvable = assigned = mq0 = 0
    uniq_j = uniq_agree_j = resolvable_j = assigned_j = junction_only = 0
    for aln in bam.fetch(contig, lo, hi):
        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
            continue
        mc, bovl = None, 0
        for name, (ct, s0, e0) in spans.items():
            ov = min(aln.reference_end, e0) - max(aln.reference_start, s0)
            if ov > bovl:
                bovl, mc = ov, name
        if mc is None:
            continue
        n += 1
        if aln.mapping_quality == 0:
            mq0 += 1
        want = set(copy_gpos[mc].values())
        base_at, introns = read_bases_and_introns(aln, want)
        g2c = {g: j for j, g in copy_gpos[mc].items()}
        obs = {}
        for g, b in base_at.items():
            if g in g2c and b is not None:
                obs[g2c[g]] = rc_base(b) if strand[mc] == "-" else b
        # read intron boundaries in transcription-spliced space (strand-agnostic via gen2off)
        g2o = gen2off[mc]
        rbounds = set()
        for d0, a0 in introns:
            o1, o2 = g2o.get(d0 - 1), g2o.get(a0)
            if o1 is not None and o2 is not None:
                rbounds.add(max(o1, o2))
        # PSV-only
        _, m_p, _, nd_p, _ = assign_read(obs, copy_vecs)
        # PSV + junctions
        best, m_j, _, nd_j, _ = assign_read(obs, copy_vecs, read_bounds=rbounds, copy_bounds=cbounds)
        if nd_p >= 1:
            resolvable += 1
            if m_p >= MARGIN:
                assigned += 1
        if nd_j >= 1:
            resolvable_j += 1
            if nd_p == 0:
                junction_only += 1
            if m_j >= MARGIN:
                assigned_j += 1
                if aln.mapping_quality > 0:
                    uniq_j += 1
                    uniq_agree_j += (best == mc)
    return dict(n=n, mq0=mq0, psv_cols=len(cols), resolvable=resolvable, assigned=assigned,
                uniq=uniq, uniq_agree=uniq_agree, resolvable_j=resolvable_j, assigned_j=assigned_j,
                junction_only=junction_only, uniq_j=uniq_j, uniq_agree_j=uniq_agree_j)


def load_spliced():
    fa = pysam.FastaFile(SPLICED_FA)
    return {r: fa.fetch(r) for r in fa.references}


def validate_real(limit=25):
    meta = load_meta()
    skel = load_skel()
    spliced = load_spliced()
    rescued_fam = load_rescued(meta, skel, spliced)
    fams = colocated_families(meta, rescued_fam=rescued_fam)
    bam = pysam.AlignmentFile(GGO_BAM, "rb")
    n_resc = sum(len(v) for v in rescued_fam.values())
    print(f"merged {n_resc} family-aware rescued copies into the rosters")
    print(f"co-located families (>=3 same-chrom tandem copies): {len(fams)}; running top {limit} by size\n")
    print(f"{'family':9s} {'gene':9s} {'cop':>3} {'reads':>6} {'MAPQ0%':>6} {'resolv_PSV':>10} "
          f"{'resolv_+J':>9} {'J_only':>6} {'assign_+J':>9} {'uniqAgree':>10}")
    agg = defaultdict(int)
    for fid, dom, copies in sorted(fams, key=lambda x: -len(x[2]))[:limit]:
        try:
            st = assign_family(spliced, bam, copies, meta, skel)
        except Exception as ex:
            print(f"{fid:9s} {dom:9s} ERROR {type(ex).__name__}: {ex}")
            continue
        if st["n"] == 0:
            continue
        f = lambda k: 100 * st[k] / st["n"]
        ua = f"{st['uniq_agree_j']}/{st['uniq_j']}" if st["uniq_j"] else "-"
        print(f"{fid:9s} {dom:9s} {len(copies):>3} {st['n']:>6} {f('mq0'):>5.0f}% {f('resolvable'):>9.0f}% "
              f"{f('resolvable_j'):>8.0f}% {st['junction_only']:>6} {f('assigned_j'):>8.0f}% {ua:>10}")
        for k in ("n", "mq0", "resolvable", "assigned", "resolvable_j", "assigned_j", "junction_only",
                  "uniq_j", "uniq_agree_j"):
            agg[k] += st[k]
    N = max(agg["n"], 1)
    print(f"\nAGGREGATE over {limit} families: reads={agg['n']:,} MAPQ0={100*agg['mq0']/N:.0f}%")
    print(f"  PSV-only:      resolvable={100*agg['resolvable']/N:.0f}%  assigned={100*agg['assigned']/N:.0f}%")
    print(f"  PSV+junction:  resolvable={100*agg['resolvable_j']/N:.0f}%  assigned={100*agg['assigned_j']/N:.0f}%"
          f"  (+{agg['junction_only']:,} reads resolved by a copy-specific junction where PSVs alone could not)")
    print(f"  unique-mapper agreement (PSV+junction): {agg['uniq_agree_j']:,}/{agg['uniq_j']:,} "
          f"({100*agg['uniq_agree_j']/max(agg['uniq_j'],1):.1f}%)  <- silver-standard accuracy")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "sim5x"
    if mode == "sim5x":
        validate_sim5x()
    elif mode == "real":
        validate_real()
    else:
        sys.exit(f"unknown mode {mode}")
