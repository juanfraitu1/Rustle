#!/usr/bin/env python3
"""O2's PSV certificate on the reads minimap2 leaves almost tied (2026-10-19, /goal: "establish a working O2 regimen from
both simulations and real reads, O2 must get almost tied multi-mapping reads and help assign via PSVs from which copy did they
actually come from"). It reads any BAM: a simulated world's E0.bam (=/X CIGARs) or a real library's CHM13 alignments (M + NM).

A read is O2's when its alignments within 10 AS of its best (Eichler's rule) lie on two or more copies of one family and on
no locus off it. Each alignment belongs to one copy (summary.attribute: strand, then exon overlap). Per tied copy, the PSV
evidence is the substitutions of its best tied alignment. Clips and indels are left out: on real reads a donor's copy can
hold sequence CHM13's lacks (HG002 development: counting clipped bases picked the wrong copy for 225 of 228 reads at a
3+ edit margin), and HiFi errors are mostly homopolymer indels. The read is placed on the copy with the fewest substitutions
when every other tied copy has at least MIN_PSV = 2 more -- two agreeing PSVs, since one can be the donor's own allele
(HG002: one-PSV placements 85% right, two-PSV 98-99%) -- and that copy explains the read within the HiFi substitution floor
(FIT = 0.1% of its aligned bases). Otherwise the read keeps the set of copies no two PSVs exclude; a read its best copy does
not explain is O3's (psv_unexplained).

    o2psv.py copies-world WORLD OUT.tsv               the copy table of a simulated world (copies.tsv + isoforms.json)
    o2psv.py copies-pool FAMILIES.json SET OUT.tsv    the copy table of a family pool
    o2psv.py run BAM COPIES.tsv OUT.tsv               O2's reads: read, family, call, reason, candidates, set, primary
    o2psv.py svp READS.bam TEACH.bam COPIES.tsv GENOME.fa OUT.tsv [--k 3] [--exclude TSV] [--names TSV]
                 [--fit RATE | --fit-quantile Q]   the explained floor: a rate, or the sample's own Q-quantile (error_floor)
                 [--teach-margin M]                the AS margin a read placed on a copy needs to validate it (10)
                                                      the same reads decided by sample-validated PSVs (the regimen)
    o2psv.py score OUT.tsv TRUTH.tsv ... [--donor]    truth: read + copy (worlds, fragments), or hg002truth's classes

Sample-validated PSVs (svp): one PSV is enough in one genome (round 23) but not on a real donor, whose copy may carry another
copy's allele (HG002: one-PSV placements 85-94% right). So a PSV counts against copy Y only where the sample's own confident
reads on Y (Eichler's well placed: every rival >= 10 AS lower) cover it at least k times and none shows the read's base: the
sample's Y lacks that allele. The copy no such PSV excludes takes the read; several stay a set.
"""
import collections
import csv
import importlib.util
import json
import os
import sys

SLACK, MIN_PSV, FIT = 9, 2, 0.001
HERE = os.path.dirname(os.path.abspath(__file__))
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_psv", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def substitutions(cigartuples, nm):
    """(substitutions, aligned bases) of one alignment: its X operations under =/X CIGARs, else NM less its inserted and
    deleted bases."""
    aligned = sum(n for op, n in cigartuples if op in (0, 7, 8))
    if any(op in (7, 8) for op, _ in cigartuples):
        return sum(n for op, n in cigartuples if op == 8), aligned
    indel = sum(n for op, n in cigartuples if op in (1, 2))
    return nm - indel, aligned


def tied_copies(alns):
    """{copy: (substitutions, aligned bases) of its best tied alignment} over a read's alignments [{AS, copy, mism, alen}]
    within SLACK of its best AS; None when one of them lies off every copy or the copies span two families."""
    best = max(a["AS"] for a in alns)
    out = {}
    for a in alns:
        if a["AS"] < best - SLACK:
            continue
        if a["copy"] is None:
            return None
        o = out.get(a["copy"])
        if o is None or (a["mism"], -a["alen"]) < (o[0], -o[1]):
            out[a["copy"]] = (a["mism"], a["alen"])
    if len({c.rsplit(".", 1)[0] for c in out}) > 1:
        return None
    return out


def quantile_rate(rates, q):
    """The q-quantile of substitution rates (FIT when there are none)."""
    r = sorted(rates)
    return r[int(q * (len(r) - 1))] if r else FIT


def error_floor(bam_p, q, n=20000, min_aligned=500):
    """The substitution rate q of the sample's uniquely placed reads (MAPQ 60 primaries, >= min_aligned aligned bases, the
    first n in the file) stay within: a read whose best copy explains it worse is explained by no copy (O3's). It holds the
    sample's sequencing errors and the donor's alleles, so it is the sample's own floor, not a constant."""
    import pysam
    rates = []
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary or rec.is_unmapped or rec.mapping_quality < 60 or not rec.has_tag("NM"):
                continue
            m, a = substitutions(rec.cigartuples, rec.get_tag("NM"))
            if a >= min_aligned:
                rates.append(m / a)
                if len(rates) >= n:
                    break
    return quantile_rate(rates, q)


def psv_decide(cands, fit=FIT):
    """(call or "-", reason, set) of a read tied on cands {copy: (substitutions, aligned bases)}."""
    order = sorted(cands, key=lambda c: (cands[c][0], c))
    best = order[0]
    keep = sorted(c for c in cands if cands[c][0] - cands[best][0] < MIN_PSV)
    if cands[best][0] > fit * cands[best][1]:
        return "-", "psv_unexplained", sorted(cands)
    if len(keep) == 1:
        return best, "psv_certificate", keep
    return "-", "psv_set", keep


def aligned_pairs(start, cigar):
    """[(read position, reference position)] of an alignment's aligned (M/=/X) bases. Read positions count soft clips and a
    leading hard clip, so every alignment of a read shares its coordinates."""
    out, q, r = [], 0, start
    for i, (op, n) in enumerate(cigar):
        if op == 5 and i == 0:
            q += n
        elif op in (0, 7, 8):
            out.extend((q + j, r + j) for j in range(n))
            q += n
            r += n
        elif op in (1, 4):
            q += n
        elif op in (2, 3):
            r += n
    return out


def psv_columns(length, prim_rev, alns, fetch):
    """{read position (the primary's orientation): {copy: (chrom, reference position, the copy's base read in the read's
    orientation, flipped)}} at the positions every tied copy aligns and their bases differ. alns: {copy: (chrom, start,
    cigar, reverse)}; fetch(chrom, start, end) -> reference sequence."""
    cols = collections.defaultdict(dict)
    for c, (chrom, start, cigar, rev) in alns.items():
        pairs = aligned_pairs(start, cigar)
        if not pairs:
            continue
        lo, hi = min(p for _, p in pairs), max(p for _, p in pairs) + 1
        ref = fetch(chrom, lo, hi).upper()
        flip = rev != prim_rev
        for qp, rp in pairs:
            b = ref[rp - lo]
            cols[(length - 1 - qp) if flip else qp][c] = (chrom, rp, b.translate(COMP) if flip else b, flip)
    return {k: d for k, d in cols.items() if len(d) == len(alns) and len({v[2] for v in d.values()}) > 1}


def excluded(cols, read, teach, k, own=None):
    """The tied copies the read's own PSVs exclude in this sample: at a column where the read carries another tied copy's base,
    copy Y is excluded when the sample's reads placed on Y (teach: {(copy, chrom, position): {read: genome-strand base}}),
    the read itself (own) left out, cover Y's position at least k times and none shows the read's base -- the sample's Y lacks
    the read's allele there. A read misplaced on Y can only void a PSV, never validate one."""
    out = set()
    for q, d in cols.items():
        rb = read[q].upper()
        if not any(v[2] == rb for v in d.values()):
            continue
        for y, (chrom, pos, b, flip) in d.items():
            if b == rb:
                continue
            obs = [b2 for r, b2 in teach.get((y, chrom, pos), {}).items() if r != own]
            if len(obs) >= k and (rb.translate(COMP) if flip else rb) not in obs:
                out.add(y)
    return out


def svp_decide(cands, excl, fit=FIT):
    """(call or "-", reason, set) once the sample-validated PSVs exclude excl among the tied copies cands {copy:
    (substitutions, aligned bases)}: the one copy left takes the read; several stay a set; none left means the read carries
    alleles of several copies that their own reads lack (a mosaic: O3's); a read its best copy does not explain within FIT is
    O3's too."""
    best = min(cands, key=lambda c: (cands[c][0], c))
    if cands[best][0] > fit * cands[best][1]:
        return "-", "psv_unexplained", sorted(cands)
    left = sorted(c for c in cands if c not in excl)
    if not left:
        return "-", "svp_contradiction", sorted(cands)
    if len(left) == 1:
        return left[0], "svp_certificate", left
    return "-", "svp_set", left


def contradicted(cols, read):
    """The tied copies the read's own bases speak against: at some PSV column the read carries another copy's base, not
    theirs (a sequencing error or a donor allele on one PSV makes a read look like a mosaic of two copies)."""
    out = set()
    for q, d in cols.items():
        rb = read[q].upper()
        if any(v[2] == rb for v in d.values()):
            out |= {c for c, v in d.items() if v[2] != rb}
    return out


def regimen_decide(cands, excl, contra, fit=FIT):
    """O2's regimen (2026-10-19 /goal) on a read tied within 10 AS: (call or "-", reason, set). A read its best copy does not
    explain within fit is O3's. It is placed on the copy with at least MIN_PSV fewer substitutions than every other
    (psv_two: two PSVs, no sample evidence needed), or on the one copy the sample-validated PSVs leave (excl: excluded())
    when the read's own bases contradict it nowhere (psv_validated). Disagreeing rules abstain (psv_conflict). Otherwise the
    read keeps the copies neither rule removes (psv_set)."""
    best = min(cands, key=lambda c: (cands[c][0], c))
    if cands[best][0] > fit * cands[best][1]:
        return "-", "psv_unexplained", sorted(cands)
    keep = sorted(c for c in cands if cands[c][0] - cands[best][0] < MIN_PSV)
    two = best if len(keep) == 1 else None
    left = sorted(c for c in cands if c not in excl)
    val = left[0] if len(left) == 1 and left[0] not in contra else None
    if two and val and two != val:
        return "-", "psv_conflict", sorted(cands)
    if two:
        return two, "psv_two", [two]
    if val:
        return val, "psv_validated", [val]
    held = [c for c in keep if c not in excl]
    return "-", "psv_set", held if held and not set(held) <= contra else keep


def fate(call, held, truth, cls=None):
    """(fate, set holds the truth) of one read. truth: its copy (None for a donor class without one); cls: hg002truth's
    class -- unique / unique_elsewhere are judged on the donor's copy, own_tie placements are unsupported, elsewhere and
    no_counterpart placements wrong, and poor_match / unmapped reads (the donor's genome fits them badly) unscored."""
    if cls in (None, "unique", "unique_elsewhere"):
        f = "abstain" if call == "-" else "right" if call == truth else "wrong"
        return f, truth in held
    if cls == "own_tie":
        return ("abstain" if call == "-" else "unsupported"), None
    if cls in ("poor_match", "unmapped"):
        return "unscored", None
    return ("abstain" if call == "-" else "wrong"), False


def copies_world(w, out_p):
    exons = collections.defaultdict(set)
    for cps in json.load(open(f"{w}/isoforms.json")).values():
        for c in cps:
            for ex in c["isoforms"].values():
                exons[c["copy"]] |= {tuple(x) for x in ex}
    with open(out_p, "w") as fh:
        fh.write("copy\tfamily\tchrom\tstart\tend\tstrand\texons\n")
        for c in csv.DictReader(open(f"{w}/copies.tsv"), delimiter="\t"):
            if c["in_reference"] == "1":
                ex = sorted(exons.get(c["copy"], ())) or [(int(c["start"]), int(c["end"]))]
                fh.write(f"{c['copy']}\t{c['family']}\t{c['chrom']}\t{c['start']}\t{c['end']}\t{c['strand']}\t"
                         + ",".join(f"{s}-{e}" for s, e in ex) + "\n")
    return 0


def copies_pool(fam_p, set_name, out_p):
    with open(out_p, "w") as fh:
        fh.write("copy\tfamily\tchrom\tstart\tend\tstrand\texons\n")
        for f in json.load(open(fam_p))["sets"][set_name]:
            for k, c in enumerate(f["copies"]):
                ex = sorted({tuple(x) for v in c["transcripts"].values() for x in v})
                fh.write(f"{f['family']}.{k}\t{f['family']}\t{c['chrom']}\t{c['start']}\t{c['end']}\t{c['strand']}\t"
                         + ",".join(f"{s}-{e}" for s, e in ex) + "\n")
    return 0


def load_copies(p):
    """{chrom: [(copy, strand, exons, span)]} for summary.attribute."""
    out = collections.defaultdict(list)
    for r in csv.DictReader(open(p), delimiter="\t"):
        ex = [tuple(int(x) for x in e.split("-")) for e in r["exons"].split(",")]
        out[r["chrom"]].append((r["copy"], r["strand"], ex, (int(r["start"]), int(r["end"]))))
    return out


def _strand(rec):
    """The transcript's genomic strand: the alignment's, flipped when minimap2's ts tag says the read is antisense."""
    s = "-" if rec.is_reverse else "+"
    if rec.has_tag("ts") and rec.get_tag("ts") == "-":
        s = "+" if s == "-" else "-"
    return s


def cmd_run(bam_p, copies_p, out_p):
    import pysam
    att = _sibling("summary").attribute
    cps = load_copies(copies_p)
    alns = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_supplementary or not rec.has_tag("AS"):
                continue
            mism, alen = substitutions(rec.cigartuples, rec.get_tag("NM") if rec.has_tag("NM") else 0)
            alns[rec.query_name].append(dict(AS=rec.get_tag("AS"), mism=mism, alen=alen, prim=not rec.is_secondary,
                                             copy=att(rec.get_blocks(), _strand(rec), cps.get(rec.reference_name, ()))))
    n = collections.Counter()
    with open(out_p, "w") as fh:
        fh.write("read\tfamily\tcall\treason\tcandidates\tset\tprimary\n")
        for q, al in sorted(alns.items()):
            prim = next((a["copy"] for a in al if a["prim"]), None)
            if prim is None:
                continue
            tc = tied_copies(al)
            if tc is None:
                n["tied off the family"] += 1
                continue
            if len(tc) < 2:
                n["one copy (minimap2's)"] += 1
                continue
            call, reason, held = psv_decide(tc)
            n[reason] += 1
            fh.write(f"{q}\t{prim.rsplit('.', 1)[0]}\t{call}\t{reason}\t{','.join(sorted(tc))}\t{','.join(held)}\t{prim}\n")
    print(f"[o2psv] {bam_p} -> {out_p}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)
    return 0


def _domain(bam_p, cps, att, names=None):
    """{read: {copy: (AS, substitutions, aligned, chrom, start, cigar, reverse)} of its best tied alignment per copy} for the
    reads O2 decides (tied within SLACK on >= 2 copies of one family, none off them), their primaries' (sequence, reverse,
    copy); two passes so only those reads' CIGARs are kept. names: only these reads."""
    import pysam
    light = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_supplementary or not rec.has_tag("AS") or names is not None and rec.query_name not in names:
                continue
            mism, alen = substitutions(rec.cigartuples, rec.get_tag("NM") if rec.has_tag("NM") else 0)
            light[rec.query_name].append(dict(AS=rec.get_tag("AS"), mism=mism, alen=alen, prim=not rec.is_secondary,
                                              copy=att(rec.get_blocks(), _strand(rec), cps.get(rec.reference_name, ()))))
    keep = {q for q, al in light.items() if any(a["prim"] and a["copy"] for a in al)
            and (lambda tc: tc is not None and len(tc) >= 2)(tied_copies(al))}
    del light
    alns, prim = collections.defaultdict(list), {}
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            q = rec.query_name
            if q not in keep or rec.is_unmapped or rec.is_supplementary or not rec.has_tag("AS"):
                continue
            mism, alen = substitutions(rec.cigartuples, rec.get_tag("NM") if rec.has_tag("NM") else 0)
            cp = att(rec.get_blocks(), _strand(rec), cps.get(rec.reference_name, ()))
            alns[q].append((rec.get_tag("AS"), mism, alen, cp, rec.reference_name, rec.reference_start,
                            tuple(rec.cigartuples), rec.is_reverse))
            if not rec.is_secondary:
                prim[q] = (rec.query_sequence, rec.is_reverse, cp)
    out = {}
    for q, al in alns.items():
        if q not in prim or prim[q][0] is None:
            continue
        best = max(a[0] for a in al)
        per = {}
        for a in al:
            if a[0] >= best - SLACK and (a[3] not in per or (a[1], -a[2]) < (per[a[3]][1], -per[a[3]][2])):
                per[a[3]] = a
        out[q] = per
    return out, prim


def teacher_bases(bam_p, cps, att, need, skip=frozenset(), min_margin=SLACK + 1):
    """{(copy, chrom, position): {read: genome-strand base}} of each copy's placed reads -- primary on the copy, every
    alignment on another copy of its family at least min_margin AS lower (SLACK + 1: Eichler's well placed; 1: any read
    minimap2 resolved) -- at the positions in need ({copy: {(chrom, position)}}); skip: reads that may not teach (a
    fragment's own source read)."""
    import pysam
    span = {c: (ch, s, e) for ch, v in cps.items() for c, _, _, (s, e) in v}
    fam = collections.defaultdict(list)
    for c in span:
        fam[c.rsplit(".", 1)[0]].append(c)
    out = collections.defaultdict(dict)
    with pysam.AlignmentFile(bam_p) as bam:
        for f in sorted({c.rsplit(".", 1)[0] for c in need}):
            best = collections.defaultdict(dict)                 # read -> {copy: best AS}; primary copy apart
            prim = {}
            for c in fam[f]:
                ch, s, e = span[c]
                for rec in bam.fetch(ch, s, e):
                    if rec.is_unmapped or rec.is_supplementary or not rec.has_tag("AS") or rec.query_name in skip:
                        continue
                    cp = att(rec.get_blocks(), _strand(rec), cps.get(rec.reference_name, ()))
                    if cp is None or cp.rsplit(".", 1)[0] != f:
                        continue
                    q = rec.query_name
                    best[q][cp] = max(best[q].get(cp, -10 ** 9), rec.get_tag("AS"))
                    if not rec.is_secondary:
                        prim[q] = (cp, rec.get_tag("AS"))
            teach = collections.defaultdict(set)
            for q, (cp, a) in prim.items():
                if all(v <= a - min_margin for c2, v in best[q].items() if c2 != cp):
                    teach[cp].add(q)
            for c in fam[f]:
                if c not in need or not teach.get(c):
                    continue
                by = collections.defaultdict(set)
                for ch, p in need[c]:
                    by[ch].add(p)
                for ch, want in by.items():
                    ps = sorted(want)
                    for rec in bam.fetch(ch, ps[0], ps[-1] + 1):
                        if rec.is_secondary or rec.is_supplementary or rec.query_name not in teach[c] or rec.query_sequence is None:
                            continue
                        sq = rec.query_sequence
                        for qp, rp in rec.get_aligned_pairs(matches_only=True):
                            if rp in want:
                                out[(c, ch, rp)][rec.query_name] = sq[qp].upper()
    return out


def cmd_svp(reads_p, teach_p, copies_p, genome_p, out_p, k=3, exclude_p=None, names_p=None, fit=None, fit_q=None,
            min_margin=SLACK + 1):
    """O2 with sample-validated PSVs: o2psv.py svp READS.bam TEACH.bam COPIES.tsv GENOME.fa OUT.tsv [--k K]
    [--exclude TSV (first column: reads that may not teach)] [--names TSV (first column: the reads to decide)]."""
    import pysam
    att = _sibling("summary").attribute
    cps = load_copies(copies_p)
    skip = {ln.split("\t")[0] for ln in open(exclude_p)} if exclude_p else frozenset()
    names = {ln.split("\t")[0].strip() for ln in open(names_p)} if names_p else None
    dom, prim = _domain(reads_p, cps, att, names)
    fa = pysam.FastaFile(genome_p)
    cols, need = {}, collections.defaultdict(set)
    for q, per in dom.items():
        seq, rev, _ = prim[q]
        cols[q] = psv_columns(len(seq), rev, {c: (a[4], a[5], a[6], a[7]) for c, a in per.items()}, fa.fetch)
        for d in cols[q].values():
            for c, (ch, p, _, _) in d.items():
                need[c].add((ch, p))
    tb = teacher_bases(teach_p, cps, att, need, skip, min_margin)
    if fit is None:
        fit = error_floor(teach_p, fit_q) if fit_q is not None else FIT
    n = collections.Counter()
    with open(out_p, "w") as fh:
        fh.write("read\tfamily\tcall\treason\tcandidates\tset\tprimary\n")
        for q in sorted(dom):
            per = dom[q]
            cands = {c: (a[1], a[2]) for c, a in per.items()}
            call, reason, held = regimen_decide(cands, excluded(cols[q], prim[q][0], tb, k, own=q),
                                                contradicted(cols[q], prim[q][0]), fit)
            n[reason] += 1
            p = prim[q][2]
            fh.write(f"{q}\t{p.rsplit('.', 1)[0]}\t{call}\t{reason}\t{','.join(sorted(per))}\t{','.join(held)}\t{p}\n")
    print(f"[o2psv] svp {reads_p} (teachers {teach_p}, margin >= {min_margin}, k {k}, fit {fit:.5f}) -> {out_p}: " + ", ".join(f"{a} {b}" for a, b in sorted(n.items())),
          flush=True)
    return 0


SCORE_KEYS = ("reads", "right", "wrong", "unsupported", "abstain", "set holds truth", "judged sets", "candidates", "set size",
              "minimap2 right", "minimap2 wrong", "minimap2 unsupported")


def score(rows, truth, donor=False):
    """Counts of O2's calls and sets, and of minimap2's primary, against truth {read: copy} or, donor, {read: (class,
    copy)}; reads without a truth row are left out."""
    n = dict.fromkeys(SCORE_KEYS, 0)
    for r in rows:
        t = truth.get(r["read"])
        if t is None:
            continue
        cls, tc = t if donor else (None, t)
        held = r["set"].split(",")
        f, h = fate(r["call"], held, tc, cls)
        if f == "unscored":
            continue
        n["reads"] += 1
        n[f] += 1
        if h is not None:
            n["judged sets"] += 1
            n["set holds truth"] += h
        n["candidates"] += len(r["candidates"].split(","))
        n["set size"] += len(held)
        m, _ = fate(r["primary"], [r["primary"]], tc, cls)
        n[f"minimap2 {m}"] += 1
    return n


def _truth(paths, donor):
    out = {}
    for p in paths:
        for r in csv.DictReader(open(p), delimiter="\t"):
            if donor:
                out[r["read"]] = (r["class"], r["truth"] or None)
            else:
                out[r["read"]] = r["copy"]
    return out


def cmd_score(out_p, truth_ps, donor=False):
    n = score(csv.DictReader(open(out_p), delimiter="\t"), _truth(truth_ps, donor), donor)
    placed = n["right"] + n["wrong"]
    mm = n["minimap2 right"] + n["minimap2 wrong"]
    print(", ".join(f"{k} {v}" for k, v in n.items()))
    print(f"O2: placed {placed}, precision {n['right'] / max(1, placed):.4f}; sets hold the truth "
          f"{n['set holds truth']} / {n['judged sets']}; copies {n['candidates'] / max(1, n['reads']):.2f} -> "
          f"{n['set size'] / max(1, n['reads']):.2f} | minimap2: precision {n['minimap2 right'] / max(1, mm):.4f}", flush=True)
    return 0


def main(argv):
    donor, argv = "--donor" in argv, [a for a in argv if a != "--donor"]
    if argv[:1] == ["copies-world"] and len(argv) == 3:
        return copies_world(argv[1], argv[2])
    if argv[:1] == ["copies-pool"] and len(argv) == 4:
        return copies_pool(argv[1], argv[2], argv[3])
    if argv[:1] == ["run"] and len(argv) == 4:
        return cmd_run(argv[1], argv[2], argv[3])
    if argv[:1] == ["svp"] and len(argv) >= 6:
        opt = lambda k, d=None: argv[argv.index(k) + 1] if k in argv else d  # noqa: E731
        return cmd_svp(argv[1], argv[2], argv[3], argv[4], argv[5], int(opt("--k", 3)), opt("--exclude"), opt("--names"),
                       float(opt("--fit")) if opt("--fit") else None,
                       float(opt("--fit-quantile")) if opt("--fit-quantile") else None, int(opt("--teach-margin", SLACK + 1)))
    if argv[:1] == ["score"] and len(argv) >= 3:
        return cmd_score(argv[1], argv[2:], donor)
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
