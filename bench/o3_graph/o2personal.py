#!/usr/bin/env python3
"""Personal copies against the corrector's pileup on gorilla KB3781's own reads with exact truth (PREREG_o2_personal_copies_
2026-10-20; user 2026-10-20: "yes please, also maybe we can use gorilla for this since we have the fibroblasts and T2T genome?").

A copy's PERSONAL COPY is its sequence in the sample: the reference span with the events (substitutions, left-shifted indels)
its confident reads carry -- the reads o2correct.py's corrector learns from -- split into at most two haplotypes by the reads
that carry them. A tied read is realigned to the personal copies of its tied copies and goes to the copy with the unique best
AS; a shared best or a tied copy it does not align is left to minimap2. The ORACLE does the same with the copies' true
sequences in the animal's own two assemblies; the TRUTH is the read's best placement in them.

    o2personal.py copies GWFAM.tsv GENOME.fa OUT_PREFIX         OUT.copies.tsv (o2psv.load_copies) + OUT.tx.fa (main transcripts)
    o2personal.py ava TX.fa COPIES.tsv OUT_DIR [--budget S]     each family all-vs-all (minimap2 -X -c -x asm20) -> FAMILY.paf
    o2personal.py qualify AVA.paf COPIES.tsv OUT_PREFIX         families whose every member pair aligns at [0.90, 1.00):
                                                               OUT.json {family: copies, identities, half}, OUT.copies.tsv, OUT.bed
    o2personal.py reads ARM2.tsv SUB.bam OUT.fa                 the domain reads (o2correct.py donor's rows) as sequenced
    o2personal.py queries COPIES.tsv GENOME.fa CHRMAP.tsv OUT    each copy's span +- 20 kb, by ALT haplotype -> OUT.{mat,pat}.fa
    o2personal.py orthologs COPIES.tsv CHRMAP.tsv OUT.tsv --haps mat=PAF:FAI pat=PAF:FAI       the copies' ALT orthologs
    o2personal.py truth ARM2.tsv COPIES.tsv ORTH.tsv CHRMAP.tsv QUAL.json OUT.tsv --bams mat=BAM:FAI pat=BAM:FAI
    o2personal.py personal CACHE_DIR COPIES.tsv GENOME.fa OUT    personal copies -> OUT.phased.fa, OUT.unphased.fa, OUT.haps.tsv
    o2personal.py oracle COPIES.tsv ORTH.tsv GENOME.fa OUT.fa --haps mat=FASTA pat=FASTA      the copies' true sequences
    o2personal.py decide ARM2.tsv REALIGNED.bam OUT.tsv          per read: the tied copy with the unique best AS, or none
    o2personal.py arms ARM2.tsv SUB.bam COPIES.tsv GENOME.fa OUT.tsv    O2 (fewest edits), highest AS, the tie's structure
    o2personal.py score ARM2.tsv ARMS.tsv RUN.tsv TRUTH.tsv QUAL.json HALF OUT.json     Amendment 1's tallies and verdict
    (personal, oracle, decide: the arms Amendment 1 dropped.)
"""
import collections
import csv
import hashlib
import importlib.util
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
FLANK_P = 5000                                                  # bp of reference kept either side of a copy's span


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_pc", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


C = _sibling("o2correct")                                       # K, SHARED, HET, HOM: the corrector's thresholds


def rc(s):
    return s.translate(COMP)[::-1]


def gwfam_rows(p):
    """The KB3781 catalog's copies (rna_allele/fibro_cat.copies.tsv) as {copy: FAMILY.INDEX, family, chrom, start, end,
    strand, exons}."""
    out = []
    for r in csv.DictReader(open(p), delimiter="\t"):
        out.append(dict(copy=f"{r['family_id']}.{r['copy_idx']}", family=r["family_id"], chrom=r["chrom"],
                        start=int(r["start"]), end=int(r["end"]), strand=r["strand"],
                        exons=[tuple(int(x) for x in e.split("-")) for e in r["exons"].split(",")]))
    return out


def spliced(fetch, exons, strand):
    s = "".join(fetch(a, b) for a, b in exons).upper()
    return rc(s) if strand == "-" else s


def half(family):
    """'dev' or 'heldout', fixed by the family's name (sha1 parity) before any result is seen."""
    return "dev" if int(hashlib.sha1(family.encode()).hexdigest(), 16) % 2 == 0 else "heldout"


def _covers(blocks, p):
    return any(s <= p < e for s, e in blocks)


def _depths(reads, positions):
    """{position: reads whose blocks cover it} (blocks of one read do not overlap: alignment_events' segments)."""
    import bisect
    st = sorted(s for _, bl in reads for s, _ in bl)
    en = sorted(e for _, bl in reads for _, e in bl)
    return {p: bisect.bisect_right(st, p) - bisect.bisect_right(en, p) for p in positions}


def variant_events(reads):
    """(carried, heterozygous) events of a copy's confident reads [(events, blocks)]: those the corrector calls explained (>=
    SHARED reads, >= HET of a depth >= K), split at HOM of the depth."""
    n = collections.Counter(e for ev, _ in reads for e in ev)
    d = _depths(reads, {e[1] for e in n})
    ok = {e for e, k in n.items() if k >= C.SHARED and d[e[1]] >= C.K and k >= C.HET * d[e[1]]}
    hom = {e for e in ok if n[e] >= C.HOM * d[e[1]]}
    return hom, ok - hom


def consensus(reads):
    """The events more than half the reads covering them carry (>= SHARED of them, depth >= K): one unphased personal copy."""
    n = collections.Counter(e for ev, _ in reads for e in ev)
    d = _depths(reads, {e[1] for e in n})
    return {e for e, k in n.items() if k >= C.SHARED and d[e[1]] >= C.K and k > 0.5 * d[e[1]]}


def _foot(e):
    return (e[1], e[1] + (1 if e[0] == "S" else e[2] if e[0] == "D" else 0))


def _clash(e, f):
    (a, b), (c, d) = _foot(e), _foot(f)
    if e[0] == "I" and f[0] == "I":
        return a == c
    if e[0] == "I":
        return c < a < d
    if f[0] == "I":
        return a < c < b
    return a < d and c < b


def compatible(counts):
    """The events {event: reads carrying it} one sequence can carry together: most carried first, each kept unless it
    overlaps one kept (an insertion goes before its base, so it does not clash with a substitution of that base)."""
    kept = []
    for e in sorted(counts, key=lambda e: (-counts[e], e[1], repr(e))):
        if not any(_clash(e, f) for f in kept):
            kept.append(e)
    return set(kept)


def phase(reads, het):
    """Two groups of reads by the heterozygous events they carry (2-means over the events each read covers, seeded by the
    most evenly split event) -> (events of group 0, events of group 1), each event in the group where most of its covering
    reads carry it."""
    het = sorted(het, key=lambda e: (e[1], repr(e)))
    states = [{e: int(e in ev) for e in het if _covers(bl, e[1])} for ev, bl in reads]

    def split(e):
        return min(sum(1 for s in states if s.get(e) == 1), sum(1 for s in states if s.get(e) == 0))

    seed = max(het, key=lambda e: (split(e), -e[1]))
    group = [s.get(seed) for s in states]
    prof = None
    for _ in range(50):
        prof = []
        for g in (0, 1):
            p = {}
            for e in het:
                k1 = sum(1 for s, x in zip(states, group) if x == g and s.get(e) == 1)
                k0 = sum(1 for s, x in zip(states, group) if x == g and s.get(e) == 0)
                if k1 + k0:
                    p[e] = int(k1 > k0)
            prof.append(p)
        new = []
        for s, x in zip(states, group):
            dist = [sum(1 for e, v in s.items() if e in prof[g] and prof[g][e] != v) for g in (0, 1)]
            new.append(0 if dist[0] < dist[1] else 1 if dist[1] < dist[0] else x)
        if new == group:
            break
        group = new
    return ({e for e, v in prof[0].items() if v}, {e for e, v in prof[1].items() if v})


def haplotypes(reads):
    """A copy's personal haplotypes: its carried events plus each phased group's heterozygous ones, made compatible -> one
    or two event sets."""
    hom, het = variant_events(reads)
    n = collections.Counter(e for ev, _ in reads for e in ev)
    if not het:
        return [compatible({e: n[e] for e in hom})]
    a, b = phase(reads, het)
    haps = [compatible({e: n[e] for e in hom | g}) for g in (a, b)]
    return [haps[0]] if haps[0] == haps[1] else haps


def apply_events(ref, lo, events):
    """ref = the reference from position lo, with the events (reference coordinates) written in; events outside it, or
    overlapping one already written, are skipped."""
    out, cur, last_ins = [], lo, None
    hi = lo + len(ref)
    for e in sorted(events, key=lambda e: (e[1], e[0] != "I", repr(e))):
        p = e[1]
        if p < cur or p >= hi or (e[0] != "I" and p < lo):
            continue
        if e[0] == "I":
            if p == last_ins:
                continue
            out.append(ref[cur - lo:p - lo] + e[2])
            cur, last_ins = p, p
        elif e[0] == "S":
            out.append(ref[cur - lo:p - lo] + e[2])
            cur = p + 1
        else:
            if p + e[2] > hi:
                continue
            out.append(ref[cur - lo:p - lo])
            cur = p + e[2]
    out.append(ref[cur - lo:])
    return "".join(out)


def edits_of(cigar, nm):
    """An alignment's edits: NM (substituted, inserted and deleted bases; introns are free) plus the read bases it clips."""
    return nm + sum(n for op, n in cigar if op in (4, 5))


def truth_of(alns):
    """A read's origin in the animal's own diploid genome from its alignments [(edits, reference copy or None, 'ref'|'alt')]
    (Amendment 1): the one copy whose alignments hold the fewest edits, strictly fewer than any other copy's or any off-copy
    locus's -> ('determined', copy, 'refcons' when a reference-haplotype alignment is among them, else 'altonly'); two copies
    (or a copy and an off-copy locus) at the fewest -> 'undetermined'; the fewest off every copy -> 'off-copy'."""
    if not alns:
        return "unmapped", None, None
    best = min(e for e, _, _ in alns)
    w = [(e, c, r) for e, c, r in alns if e == best]
    cs = {c for _, c, _ in w}
    if len(cs) > 1:
        return "undetermined", None, None
    c = cs.pop()
    if c is None:
        return "off-copy", None, None
    return "determined", c, "refcons" if any(r == "ref" for _, _, r in w) else "altonly"


def window_edits(events):
    """Bases substituted, inserted and deleted among alignment_events' events."""
    return sum(1 if e[0] == "S" else len(e[2]) if e[0] == "I" else e[2] for e in events)


def fewest(scores, prim):
    """The copy with the unique lowest score, else minimap2's pick."""
    m = min(scores.values())
    top = [c for c, v in scores.items() if v == m]
    return top[0] if len(top) == 1 else prim


def highest(scores, prim):
    """The copy with the unique highest score, else minimap2's pick."""
    m = max(scores.values())
    top = [c for c, v in scores.items() if v == m]
    return top[0] if len(top) == 1 else prim


def match_columns(cigar, flip, length):
    """(read positions an =, read positions an X) of an --eqx alignment, in the primary's orientation (flip: the alignment
    runs the other way)."""
    eq, x, q = set(), set(), 0
    for op, n in cigar:
        if op == 7:
            eq.update(range(q, q + n))
        elif op == 8:
            x.update(range(q, q + n))
        if op in (0, 1, 4, 5, 7, 8):
            q += n
    if flip:
        eq, x = {length - 1 - i for i in eq}, {length - 1 - i for i in x}
    return eq, x


def contradicted(w, cols, window):
    """Whether the read's bases contradict copy w anywhere in the window: a read base w mismatches (X) that a rival's
    alignment matches (=), i.e. a PSV column where the read carries the rival's base. cols: {copy: (=, X) positions}."""
    lo, hi = window
    return any(lo <= i < hi and any(i in eq for c, (eq, _) in cols.items() if c != w) for i in cols[w][1])


def window_primary(spans, length):
    """The read bases, in the primary's orientation, every tied alignment aligns: spans {copy: ((first, last + 1) in its own
    orientation, whether it runs opposite the primary)}."""
    p = [(length - e, length - s) if rev else (s, e) for (s, e), rev in spans.values()]
    return max(s for s, _ in p), min(e for _, e in p)


def structure(introns):
    """'intron-differ' when the tied copies' alignments {copy: introns} differ in intron count, else 'same'."""
    return "intron-differ" if len(set(introns.values())) > 1 else "same"


def decide(best_as, tie):
    """(copy, reason) of a realignment {copy: best AS}: the tied copy with the unique best AS; a shared best ('tie') or a
    tied copy without an alignment ('incomplete') leaves the read to minimap2 (None)."""
    if any(c not in best_as for c in tie):
        return None, "incomplete"
    s = {c: best_as[c] for c in tie}
    m = max(s.values())
    top = [c for c, v in s.items() if v == m]
    return (top[0], "decided") if len(top) == 1 else (None, "tie")


# ---------------------------------------------------------------------------------------------------------------- commands

def _write_copies(rows, p):
    with open(p, "w") as fh:
        fh.write("copy\tfamily\tchrom\tstart\tend\tstrand\texons\n")
        for r in rows:
            fh.write(f"{r['copy']}\t{r['family']}\t{r['chrom']}\t{r['start']}\t{r['end']}\t{r['strand']}\t"
                     + ",".join(f"{a}-{b}" for a, b in r["exons"]) + "\n")


def read_copies(p):
    return [dict(r, start=int(r["start"]), end=int(r["end"]),
                 exons=[tuple(int(x) for x in e.split("-")) for e in r["exons"].split(",")])
            for r in csv.DictReader(open(p), delimiter="\t")]


def cmd_copies(gwfam_p, genome_p, out):
    import pysam
    rows = gwfam_rows(gwfam_p)
    fa = pysam.FastaFile(genome_p)
    _write_copies(rows, f"{out}.copies.tsv")
    with open(f"{out}.tx.fa", "w") as fh:
        for r in rows:
            fh.write(f">{r['copy']}\n{spliced(lambda a, b: fa.fetch(r['chrom'], a, b), r['exons'], r['strand'])}\n")
    print(f"[o2personal] {len(rows)} copies of {len({r['family'] for r in rows})} families -> {out}.copies.tsv, {out}.tx.fa",
          flush=True)
    return 0


PRETEST = 20                                                    # a large family is first tested on this many members


def _ava(names, seqs, threads=2):
    import subprocess
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        with open(f"{d}/f.fa", "w") as fh:
            for n in names:
                fh.write(f">{n}\n{seqs[n]}\n")
        return subprocess.run(["minimap2", "-X", "-c", "-x", "asm20", "-t", str(threads), f"{d}/f.fa", f"{d}/f.fa"],
                              capture_output=True, text=True, check=True).stdout


def cmd_ava(tx_p, copies_p, out_dir, budget=None):
    """Each family's members aligned all-vs-all (minimap2 -X -c -x asm20, the human list's command, one family at a time:
    only same-family pairs count) -> OUT_DIR/FAMILY.paf. A family of more than PRETEST members is first aligned on PRETEST
    of them (sha1 order); a pair failing there fails the family, and only that subset is kept. Resumes; exit 3 at budget."""
    import time
    import pysam
    PF = _sibling("o2perfamily")
    t0 = time.time()
    seqs = {r.name: r.sequence for r in pysam.FastxFile(tx_p)}
    fam = collections.defaultdict(list)
    for r in read_copies(copies_p):
        fam[r["family"]].append(r["copy"])
    os.makedirs(out_dir, exist_ok=True)
    n = collections.Counter()
    for f in sorted(fam):
        cps, done = fam[f], os.path.join(out_dir, f"{f}.paf")
        if len(cps) < 2 or os.path.exists(done):
            continue
        if budget is not None and time.time() - t0 > budget:
            print(f"[o2personal] ava: stopped at the {budget} s budget before {f}; rerun to resume", flush=True)
            return 3
        paf = None
        if len(cps) > PRETEST:
            sub = sorted(cps, key=lambda c: hashlib.sha1(c.encode()).hexdigest())[:PRETEST]
            paf = _ava(sub, seqs)
            if PF.qualifies(len(sub), list(PF.pair_identities(paf.splitlines()).values())):
                paf = None
            else:
                n["failed the pretest"] += 1
        if paf is None:
            paf = _ava(cps, seqs)
            n["aligned whole"] += 1
        open(done + ".tmp", "w").write(paf)
        os.replace(done + ".tmp", done)
    print(f"[o2personal] ava -> {out_dir}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)
    return 0


def cmd_qualify(paf_p, copies_p, out):
    PF = _sibling("o2perfamily")
    rows = read_copies(copies_p)
    fam = collections.defaultdict(list)
    for r in rows:
        fam[r["family"]].append(r)
    ids = collections.defaultdict(list)
    for (a, b), v in PF.pair_identities(open(paf_p)).items():
        ids[a.rsplit(".", 1)[0]].append(v)
    q = {f: dict(copies=len(v), min_identity=min(ids[f]), max_identity=max(ids[f]), half=half(f))
         for f, v in sorted(fam.items()) if PF.qualifies(len(v), ids[f])}
    json.dump(q, open(f"{out}.json", "w"), indent=1)
    keep = [r for r in rows if r["family"] in q]
    _write_copies(keep, f"{out}.copies.tsv")
    with open(f"{out}.bed", "w") as fh:
        for r in sorted(keep, key=lambda r: (r["chrom"], r["start"])):
            fh.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['copy']}\n")
    n = collections.Counter(v["half"] for v in q.values())
    print(f"[o2personal] {len(q)} of {len(fam)} families qualify ({len(keep)} copies; dev {n['dev']}, held-out "
          f"{n['heldout']}) -> {out}.json, .copies.tsv, .bed", flush=True)
    return 0


def names_by_length(sq, fai):
    """{name in an index or BAM header: the assembly's contig} -- the same name, else the one contig of that length."""
    by = collections.defaultdict(list)
    for n, ln in fai.items():
        by[ln].append(n)
    return {n: n if n in fai else by[ln][0] for n, ln in sq.items() if n in fai or len(by.get(ln, ())) == 1}


def chrmap_rows(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def reference_contigs(chrmap_p):
    """{haplotype contig the reference is identical to: (reference chromosome, its haplotype)}."""
    return {r["same_name"]: (r["pri"], r["same_hap"]) for r in chrmap_rows(chrmap_p) if r["seq_check"] == "identical"}


def alt_of(chrmap_p):
    """{reference chromosome: (ALT haplotype, the homologous ALT contig)}; chrX and chrY have none."""
    return {r["pri"]: (r["B_hap"], r["B_name"]) for r in chrmap_rows(chrmap_p) if r["B_name"]}


def alt_rows(orth, strand_of):
    """The ALT orthologs [{copy, contig, start, end, strand, exons}] as copy rows {contig: [(copy, strand, exons, span)]}
    for summary.attribute: the copy's strand, turned when its ortholog lies reversed."""
    out = collections.defaultdict(list)
    for o in orth:
        s = strand_of[o["copy"]]
        if o["strand"] == "-":
            s = "+" if s == "-" else "-"
        out[o["contig"]].append((o["copy"], s, [tuple(e) for e in o["exons"]], (o["start"], o["end"])))
    return dict(out)


def resolve_by_de(orth):
    """{(copy, hap): ortholog} where two copies' counterparts overlap on one contig: the copy whose query aligned with the
    lower divergence (de) keeps it, the other gets none (Amendment 3 of PREREG_o2_corrector_families; donorsim.drop_ambiguous
    drops both)."""
    kept = {}
    for k, o in sorted(orth.items(), key=lambda kv: (kv[1].get("de") if kv[1].get("de") is not None else 1.0, kv[0])):
        if not any(k[1] == j[1] and o["tname"] == p["tname"] and o["start"] < p["end"] and p["start"] < o["end"]
                   for j, p in kept.items()):
            kept[k] = o
    return kept


def read_orthologs(p):
    return [dict(r, start=int(r["start"]), end=int(r["end"]), mapq=int(r["mapq"]), exons=json.loads(r["exons"]))
            for r in csv.DictReader(open(p), delimiter="\t")]


def tally(rows, arm):
    """Counts of an arm's calls against the truth and against minimap2's on the same reads."""
    n = collections.Counter()
    for r in rows:
        c, m, t = r[arm], r["minimap2"], r["truth"]
        n["reads"] += 1
        n["right"] += c == t
        if c != m:
            n["fixes" if c == t else "breaks" if m == t else "both wrong"] += 1
    return n


def families_vs(rows, a, b):
    """(families where arm a is right more often than arm b, families where less often)."""
    d = collections.Counter()
    for r in rows:
        d[r["family"]] += (r[a] == r["truth"]) - (r[b] == r["truth"])
    return sum(1 for v in d.values() if v > 0), sum(1 for v in d.values() if v < 0)


def families_fix_break(rows, arm):
    """(families where the arm fixes more of minimap2's calls than it breaks, families where it breaks more)."""
    d = collections.Counter()
    for r in rows:
        if r[arm] != r["minimap2"]:
            d[r["family"]] += (r[arm] == r["truth"]) - (r["minimap2"] == r["truth"])
    return sum(1 for v in d.values() if v > 0), sum(1 for v in d.values() if v < 0)


def verdict(b, min_gain):
    """PASS when every bar holds: gain >= min_gain points; fixes >= 5 x breaks; >= 11 families better and no more than a
    fifth as many (rounded up) worse; on refcons reads no more than 0.2 points below minimap2."""
    ok = (b["gain"] >= min_gain and b["fixes"] >= 5 * b["breaks"] and b["better"] >= 11
          and b["worse"] <= math.ceil(b["better"] / 5) and b["refcons_drop"] <= 0.2)
    return "PASS" if ok else "FAIL"


def read_domain(arm2_p):
    """{read: (family, tied copies, minimap2's copy, corrector's copy, AS margin)} from o2correct.py donor's output."""
    out = {}
    for r in csv.DictReader(open(arm2_p), delimiter="\t"):
        out[r["read"]] = (r["family"], set(r["candidates"].split(",")), r["minimap2"], r["corrector"], int(r["margin"]))
    return out


def cmd_reads(arm2_p, bam_p, out_p):
    """The domain reads' sequences, as sequenced (the primary record, reverse-complemented back when reverse)."""
    import pysam
    dom = read_domain(arm2_p)
    seen = set()
    with pysam.AlignmentFile(bam_p) as bam, open(out_p, "w") as fh:
        for rec in bam.fetch(until_eof=True):
            q = rec.query_name
            if q in dom and q not in seen and not (rec.is_secondary or rec.is_supplementary or rec.is_unmapped):
                s = rec.query_sequence
                fh.write(f">{q}\n{rc(s) if rec.is_reverse else s}\n")
                seen.add(q)
    print(f"[o2personal] {len(seen)} of {len(dom)} domain reads -> {out_p}", flush=True)
    return 0 if len(seen) == len(dom) else 1


def _fai(p):
    return {x[0]: int(x[1]) for x in (ln.split("\t") for ln in open(p))}


def _region(r, length):
    return max(0, r["start"] - FLANK_P), min(length, r["end"] + FLANK_P)


def query_lo(start, flank):
    return max(0, start - flank)


def query_name(copy, chrom, start, flank):
    """A copy's query window name, copy|chrom|lo (donorsim's), shared by queries and orthologs."""
    return f"{copy}|{chrom}|{query_lo(start, flank)}"


def cmd_queries(copies_p, genome_p, chrmap_p, out, flank=None):
    """Each copy's reference span +- donorsim.FLANK, to be aligned to its ALT haplotype's whole assembly ->
    OUT.mat.fa / OUT.pat.fa (named copy|chrom|lo, donorsim's names)."""
    import pysam
    DS = _sibling("donorsim")
    alt = alt_of(chrmap_p)
    fa = pysam.FastaFile(genome_p)
    fh = {h: open(f"{out}.{h}.fa", "w") for h in ("mat", "pat")}
    n = collections.Counter()
    for r in read_copies(copies_p):
        if r["chrom"] not in alt:
            n["no ALT"] += 1
            continue
        h = alt[r["chrom"]][0]
        f = DS.FLANK if flank is None else flank
        lo, hi = query_lo(r["start"], f), min(fa.get_reference_length(r["chrom"]), r["end"] + f)
        fh[h].write(f">{query_name(r['copy'], r['chrom'], r['start'], f)}\n{fa.fetch(r['chrom'], lo, hi).upper()}\n")
        n[h] += 1
    for f in fh.values():
        f.close()
    print(f"[o2personal] queries -> {out}.{{mat,pat}}.fa: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)
    return 0


def cmd_orthologs(copies_p, chrmap_p, out_p, haps, repaired=False, flank=None):
    """Per copy, its ALT ortholog (donorsim.pick_ortholog on the copy's own exons, then donorsim.drop_ambiguous), kept only
    on the homologous ALT contig or an unplaced one. haps: {hap: (PAF, the hap's FASTA index)}."""
    DS = _sibling("donorsim")
    alt = alt_of(chrmap_p)
    placed = {r["same_name"] for r in chrmap_rows(chrmap_p)} | {r["B_name"] for r in chrmap_rows(chrmap_p) if r["B_name"]}
    refset = set(reference_contigs(chrmap_p))                   # repaired: any contig but the reference's own
    orth, n = {}, collections.Counter()
    rows = read_copies(copies_p)
    for h, (paf_p, fai_p) in haps.items():
        paf = DS.read_paf(paf_p)
        lens = {}
        for ln in open(paf_p):
            x = ln.split("\t")
            lens[x[5]] = int(x[6])
        names = names_by_length(lens, _fai(fai_p))
        for r in rows:
            if r["chrom"] not in alt or alt[r["chrom"]][0] != h:
                continue
            f = DS.FLANK if flank is None else flank
            lo = query_lo(r["start"], f)
            o = DS.pick_ortholog(paf.get(query_name(r["copy"], r["chrom"], r["start"], f), []),
                                 {"main": [(a - lo, b - lo) for a, b in r["exons"]]},
                                 None if repaired else DS.MAX_EXON_GROWTH)
            if o is None:
                n[f"{h}: none"] += 1
                continue
            cm = names.get(o["tname"])
            refused = cm is None or (cm in refset if repaired else (cm != alt[r["chrom"]][1] and cm in placed))
            if refused:
                n[f"{h}: on another chromosome (refused)"] += 1
                continue
            o["tname"] = cm
            orth[(r["copy"], h)] = o
            n[f"{h}: ortholog"] += 1
    kept = resolve_by_de(orth) if repaired else DS.drop_ambiguous(orth)
    n["ambiguous (dropped)"] = len(orth) - len(kept)
    with open(out_p, "w") as fh:
        fh.write("copy\thap\tcontig\tstart\tend\tstrand\tmapq\texons\n")
        for (c, h), o in sorted(kept.items()):
            fh.write(f"{c}\t{h}\t{o['tname']}\t{o['start']}\t{o['end']}\t{o['strand']}\t{o['mapq']}\t"
                     f"{json.dumps(o['isoforms']['main'], separators=(',', ':'))}\n")
    print(f"[o2personal] orthologs -> {out_p}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items()))
          + f"; {len(kept)} of {len(rows)} copies", flush=True)
    return 0


def cmd_truth(arm2_p, copies_p, orth_p, chrmap_p, qual_p, out_p, bams):
    """Each domain read's truth (truth_of, on edits: Amendment 1) from its alignments to the animal's two assemblies. bams:
    {hap: (BAM, FAI)}. Prints the development half's statuses only."""
    import pysam
    P, S = _sibling("o2psv"), _sibling("summary")
    dom = read_domain(arm2_p)
    qual = json.load(open(qual_p))
    refc = reference_contigs(chrmap_p)
    rows = read_copies(copies_p)
    idx_ref = C.span_index(P.load_copies(copies_p))
    idx_alt = C.span_index(alt_rows(read_orthologs(orth_p), {r["copy"]: r["strand"] for r in rows}))
    alns = collections.defaultdict(list)
    for h, (bam_p, fai_p) in bams.items():
        with pysam.AlignmentFile(bam_p) as bam:
            names = names_by_length({s["SN"]: s["LN"] for s in bam.header.to_dict()["SQ"]}, _fai(fai_p))
            for rec in bam.fetch(until_eof=True):
                q = rec.query_name
                if q not in dom or rec.is_unmapped or rec.is_supplementary or not rec.has_tag("NM"):
                    continue
                cm = names.get(rec.reference_name)
                if cm in refc:
                    role, cands = "ref", C.overlapping(idx_ref, refc[cm][0], rec.reference_start, rec.reference_end)
                else:
                    role = "alt"
                    cands = C.overlapping(idx_alt, cm, rec.reference_start, rec.reference_end) if cm else []
                c = S.attribute(rec.get_blocks(), P._strand(rec), cands) if cands else None
                alns[q].append((edits_of(rec.cigartuples, rec.get_tag("NM")), c, role))
    n = collections.Counter()
    with open(out_p, "w") as fh:
        fh.write("read\tfamily\thalf\tstatus\ttruth\torigin\tin_tie\tbest_edits\trunner_up\talignments\n")
        for q, (f, tie, mm, _, _) in sorted(dom.items()):
            st, t, o = truth_of(alns.get(q, []))
            inn = int(t in tie) if t else ""
            best = min((e for e, _, _ in alns[q]), default="")
            ru = min((e for e, c, _ in alns[q] if c != t), default="") if st == "determined" else ""
            fh.write(f"{q}\t{f}\t{qual[f]['half']}\t{st}\t{t or ''}\t{o or ''}\t{inn}\t{best}\t{ru}\t{len(alns.get(q, []))}\n")
            n[(qual[f]["half"], st if st != "determined" else f"determined {o}{'' if inn else ', outside the tie'}")] += 1
    print(f"[o2personal] truth -> {out_p}: " + "; ".join(f"{h} {s} {v}" for (h, s), v in sorted(n.items()) if h == "dev"), flush=True)
    return 0


def cmd_personal(cache_dir, copies_p, genome_p, out):
    """Each copy's personal haplotypes (haplotypes()) and unphased consensus from the corrector's confident reads (its
    --cache pickles) -> OUT.phased.fa (copy|P0, copy|P1), OUT.unphased.fa (copy|U0), OUT.haps.tsv. A copy without
    confident reads keeps the reference."""
    import glob
    import pickle
    import pysam
    fa = pysam.FastaFile(genome_p)
    reads = collections.defaultdict(list)
    for p in sorted(glob.glob(os.path.join(cache_dir, "*.pkl"))):
        for cp, ev, bl, _ in pickle.load(open(p, "rb")):
            reads[cp].append((ev, bl))
    n = collections.Counter()
    with open(f"{out}.phased.fa", "w") as ph, open(f"{out}.unphased.fa", "w") as up, open(f"{out}.haps.tsv", "w") as ht:
        ht.write("copy\treads\tcarried\theterozygous\thaplotypes\tevents\n")
        for r in read_copies(copies_p):
            c, rd = r["copy"], reads.get(r["copy"], [])
            lo, hi = _region(r, fa.get_reference_length(r["chrom"]))
            ref = fa.fetch(r["chrom"], lo, hi).upper()
            hom, het = variant_events(rd) if rd else (set(), set())
            haps = haplotypes(rd) if rd else [set()]
            cnt = collections.Counter(e for ev, _ in rd for e in ev)
            uni = compatible({e: cnt[e] for e in consensus(rd)}) if rd else set()
            for i, hv in enumerate(haps):
                ph.write(f">{c}|P{i}\n{apply_events(ref, lo, hv)}\n")
            up.write(f">{c}|U0\n{apply_events(ref, lo, uni)}\n")
            ht.write(f"{c}\t{len(rd)}\t{len(hom)}\t{len(het)}\t{len(haps)}\t"
                     + "|".join(",".join(f"{e[0]}{e[1]}{e[2]}" for e in sorted(hv, key=lambda e: (e[1], repr(e)))) for hv in haps)
                     + "\n")
            n[f"{len(haps)} haplotype{'s' if len(haps) > 1 else ''}"] += 1
            n["with reads"] += bool(rd)
    print(f"[o2personal] personal copies -> {out}.phased.fa, .unphased.fa: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())),
          flush=True)
    return 0


def cmd_oracle(copies_p, orth_p, genome_p, out_p, hap_fa):
    """The copies' true sequences: copy|R the reference span +- FLANK_P, copy|A its ALT ortholog +- FLANK_P. hap_fa:
    {hap: FASTA}."""
    import pysam
    fa = pysam.FastaFile(genome_p)
    hf = {h: pysam.FastaFile(p) for h, p in hap_fa.items()}
    orth = {o["copy"]: o for o in read_orthologs(orth_p)}
    n = collections.Counter()
    with open(out_p, "w") as fh:
        for r in read_copies(copies_p):
            lo, hi = _region(r, fa.get_reference_length(r["chrom"]))
            fh.write(f">{r['copy']}|R\n{fa.fetch(r['chrom'], lo, hi).upper()}\n")
            o = orth.get(r["copy"])
            if o:
                g = hf[o["hap"]]
                a, b = _region(o, g.get_reference_length(o["contig"]))
                fh.write(f">{r['copy']}|A\n{g.fetch(o['contig'], a, b).upper()}\n")
            n["with ALT" if o else "reference only"] += 1
    print(f"[o2personal] oracle copies -> {out_p}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)
    return 0


def cmd_decide(arm2_p, bam_p, out_p):
    """Each domain read's call from its realignment to per-copy targets named copy|tag (decide())."""
    import pysam
    dom = read_domain(arm2_p)
    best = collections.defaultdict(dict)
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            q = rec.query_name
            if q not in dom or rec.is_unmapped or rec.is_supplementary or not rec.has_tag("AS"):
                continue
            c, a = rec.reference_name.split("|")[0], rec.get_tag("AS")
            best[q][c] = max(best[q].get(c, -10 ** 9), a)
    n = collections.Counter()
    with open(out_p, "w") as fh:
        fh.write("read\tcall\treason\tscores\n")
        for q, (f, tie, mm, _, _) in sorted(dom.items()):
            c, why = decide(best.get(q, {}), tie)
            n[why] += 1
            n["overrules"] += c is not None and c != mm
            sc = ",".join(f"{x}:{best[q][x]}" for x in sorted(best.get(q, {})) if x in tie)
            fh.write(f"{q}\t{c or ''}\t{why}\t{sc}\n")
    print(f"[o2personal] decisions -> {out_p}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)
    return 0


def _bin(x):
    return "0.90-0.95" if x < 0.95 else "0.95-0.98" if x < 0.98 else "0.98-0.99" if x < 0.99 else "0.99-1.00"


def cmd_arms(arm2_p, bam_p, copies_p, genome_p, out_p):
    """Per domain read, the allele-blind arms from the library's own alignments (Amendment 1). Each tied copy is its
    highest-AS alignment (the corrector's choice, o2correct.collect's pass 2). O2: the unique fewest edits (window_edits)
    over the read bases every tied alignment aligns (o2correct.common_window); max_as: the unique highest AS; structure:
    whether those alignments differ in intron count; consistent (Amendment 2): O2's pick only where the read's bases
    contradict it nowhere in the common window (contradicted()), else minimap2's. A read missing a tied alignment keeps
    minimap2's pick."""
    import pysam
    P, S = _sibling("o2psv"), _sibling("summary")
    dom = read_domain(arm2_p)
    idx = C.span_index(P.load_copies(copies_p))
    fa = pysam.FastaFile(genome_p)
    recs, seqs = collections.defaultdict(dict), {}
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            q = rec.query_name
            if q not in dom or rec.is_unmapped or rec.is_supplementary or not rec.has_tag("AS"):
                continue
            cands = C.overlapping(idx, rec.reference_name, rec.reference_start, rec.reference_end)
            c = S.attribute(rec.get_blocks(), P._strand(rec), cands) if cands else None
            a = rec.get_tag("AS")
            if c in dom[q][1] and (c not in recs[q] or a > recs[q][c][0]):
                recs[q][c] = (a, rec.reference_start, tuple(rec.cigartuples), rec.is_reverse, rec.reference_name)
            if not rec.is_secondary:
                seqs[q] = (rec.query_sequence, rec.is_reverse)

    def fetcher(chrom):
        return lambda lo, hi: fa.fetch(chrom, lo, hi)

    n = collections.Counter()
    with open(out_p, "w") as fh:
        fh.write("read\tfewest_edits\tconsistent\tmax_as\tstructure\tedits\tas\tintrons\n")
        for q, (f, tie, mm, _, _) in sorted(dom.items()):
            if q not in seqs or set(recs[q]) != tie:
                n["incomplete"] += 1
                fh.write(f"{q}\t{mm}\t{mm}\t{mm}\tincomplete\t\t\t\n")
                continue
            seq, prev = seqs[q]
            win = C.common_window({c: (C.read_span(cig), rev != prev) for c, (_, _, cig, rev, _) in recs[q].items()}, len(seq))
            ed, sc, intr, cols = {}, {}, {}, {}
            for c, (a, st, cig, rev, ch) in recs[q].items():
                ev, _ = C.alignment_events(st, cig, seq if rev == prev else rc(seq), fetcher(ch), win[c])
                ed[c], sc[c], intr[c] = window_edits(ev), a, sum(1 for op, _ in cig if op == 3)
                cols[c] = match_columns(cig, rev != prev, len(seq))
            fe, ma, st_ = fewest(ed, mm), highest(sc, mm), structure(intr)
            wp = window_primary({c: (C.read_span(cig), rev != prev) for c, (_, _, cig, rev, _) in recs[q].items()}, len(seq))
            co = fe if fe == mm or not contradicted(fe, cols, wp) else mm
            n["consistent overrules"] += co != mm
            n["O2 overrules"] += fe != mm
            n["highest AS overrules"] += ma != mm
            n[st_] += 1

            def col(d):
                return ",".join(f"{c}:{d[c]}" for c in sorted(d))
            fh.write(f"{q}\t{fe}\t{co}\t{ma}\t{st_}\t{col(ed)}\t{col(sc)}\t{col(intr)}\n")
    print(f"[o2personal] arms -> {out_p}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)
    return 0


def _bin(x):
    return "0.90-0.95" if x < 0.95 else "0.95-0.98" if x < 0.98 else "0.98-0.99" if x < 0.99 else "0.99-1.00"


ARMS = ("minimap2", "o2", "o2_consistent", "fewest_subs", "max_as", "corrector")


def score_rows(dom, truth, qual, arms, run, which):
    """(status counts of the half's domain reads, the scored reads S: determined, the truth among the tied copies)."""
    status, rows = collections.Counter(), []
    for q, (f, tie, mm, corr, margin) in sorted(dom.items()):
        if qual[f]["half"] != which:
            continue
        t = truth[q]
        status[t["status"] if t["status"] != "determined" else
               f"determined {t['origin']}" + ("" if t["in_tie"] == "1" else ", truth outside the tie")] += 1
        if t["status"] != "determined" or t["in_tie"] != "1":
            continue
        a = arms[q]
        rows.append(dict(read=q, family=f, truth=t["truth"], minimap2=mm, corrector=corr, o2=a["fewest_edits"],
                         o2_consistent=a.get("consistent") or a["fewest_edits"],
                         max_as=a["max_as"], fewest_subs=run[q]["corrector"] if q in run else mm, origin=t["origin"],
                         margin=margin, structure=a["structure"], bin=_bin(qual[f]["min_identity"])))
    return status, rows


def bars_for(rows, arm, min_gain=2.0):
    """The pre-registered bars (Amendment 1) of one arm against minimap2 on the scored reads."""
    acc = {a: sum(r[a] == r["truth"] for r in rows) / len(rows) if rows else 0.0 for a in (arm, "minimap2")}
    rc_rows = [r for r in rows if r["origin"] == "refcons"]
    rc = {a: sum(r[a] == r["truth"] for r in rc_rows) / len(rc_rows) if rc_rows else 0.0 for a in (arm, "minimap2")}
    n = tally(rows, arm)
    better, worse = families_fix_break(rows, arm)
    b = dict(gain=round(100 * (acc[arm] - acc["minimap2"]), 2), fixes=n["fixes"], breaks=n["breaks"], better=better,
             worse=worse, refcons_drop=round(100 * (rc["minimap2"] - rc[arm]), 2), refcons_reads=len(rc_rows))
    return dict(b, verdict=verdict(b, min_gain))


def cmd_score(arm2_p, arms_p, run_p, truth_p, qual_p, which, out_p):
    """The amended scoring (Amendment 1) on one half: every arm against the edit truth on S, by stratum and family; the
    verdict for O2 (fewest edits) against minimap2, the other arms' numbers without verdicts."""
    dom = read_domain(arm2_p)
    qual = json.load(open(qual_p))
    truth = {r["read"]: r for r in csv.DictReader(open(truth_p), delimiter="\t")}
    arms = {r["read"]: r for r in csv.DictReader(open(arms_p), delimiter="\t")}
    run = {r["read"]: r for r in csv.DictReader(open(run_p), delimiter="\t")}
    status, rows = score_rows(dom, truth, qual, arms, run, which)

    def block(sub):
        out = {"reads": len(sub)}
        for a in ARMS:
            n = tally(sub, a)
            out[a] = dict(right=n["right"], accuracy=round(n["right"] / len(sub), 4) if sub else None, fixes=n["fixes"],
                          breaks=n["breaks"], both_wrong=n["both wrong"])
        return out

    strata = {"all": rows}
    for k, f in (("refcons", lambda r: r["origin"] == "refcons"), ("altonly", lambda r: r["origin"] == "altonly"),
                 ("exact", lambda r: r["margin"] == 0), ("near", lambda r: r["margin"] > 0),
                 ("intron-differ", lambda r: r["structure"] == "intron-differ"),
                 ("same structure", lambda r: r["structure"] == "same"),
                 ("incomplete", lambda r: r["structure"] == "incomplete")):
        strata[k] = [r for r in rows if f(r)]
    for b_ in ("0.90-0.95", "0.95-0.98", "0.98-0.99", "0.99-1.00"):
        strata["identity " + b_] = [r for r in rows if r["bin"] == b_]
    fams = sorted({r["family"] for r in rows})
    mm_wrong = sum(r["minimap2"] != r["truth"] for r in rows)
    power = len(rows) >= 300 and len(fams) >= 11 and mm_wrong >= 50
    bars = dict(P=dict(reads=len(rows), families=len(fams), minimap2_wrong=mm_wrong, holds=power))
    for a in ARMS[1:]:
        bars[a] = bars_for(rows, a)
        if not power:
            bars[a]["verdict"] = "UNDERPOWERED"
    res = dict(half=which, status=dict(status), strata={k: block(v) for k, v in strata.items()},
               families={f: block([r for r in rows if r["family"] == f]) for f in fams}, bars=bars,
               verdict=bars["o2"]["verdict"])
    json.dump(res, open(out_p, "w"), indent=1)
    print(f"[o2personal] {which}: domain status " + ", ".join(f"{k} {v}" for k, v in sorted(status.items())), flush=True)
    for k in ("all", "refcons", "altonly", "exact", "near", "intron-differ", "same structure"):
        s_ = res["strata"][k]
        print(f"  {k:15s} {s_['reads']:5d} reads | " + " | ".join(
            f"{a} {s_[a]['accuracy']}" + (f" +{s_[a]['fixes']}/-{s_[a]['breaks']}" if a != "minimap2" else "")
            for a in ARMS) if s_["reads"] else f"  {k:15s}     0 reads", flush=True)
    print(f"  P {bars['P']}", flush=True)
    for a in ARMS[1:]:
        print(f"  {a}: {bars[a]}", flush=True)
    print(f"  VERDICT (O2 = fewest edits vs minimap2): {res['verdict']}", flush=True)
    return 0


def _pairs(a, flag):
    """{key: (left, right)} from KEY=LEFT:RIGHT arguments after flag."""
    return {x.split("=", 1)[0]: tuple(x.split("=", 1)[1].split(":")) for x in a[a.index(flag) + 1:] if "=" in x}


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    cmd, a = argv[0], argv[1:]
    if cmd == "copies":
        return cmd_copies(*a[:3])
    if cmd == "ava":
        return cmd_ava(*a[:3], budget=float(a[a.index("--budget") + 1]) if "--budget" in a else None)
    if cmd == "qualify":
        return cmd_qualify(*a[:3])
    if cmd == "reads":
        return cmd_reads(*a[:3])
    fl = int(a[a.index("--flank") + 1]) if "--flank" in a else None
    if cmd == "queries":
        return cmd_queries(*a[:4], flank=fl)
    if cmd == "orthologs":
        return cmd_orthologs(*a[:3], _pairs(a, "--haps"), repaired="--repaired" in a, flank=fl)
    if cmd == "truth":
        return cmd_truth(*a[:6], _pairs(a, "--bams"))
    if cmd == "personal":
        return cmd_personal(*a[:4])
    if cmd == "oracle":
        return cmd_oracle(*a[:4], {k: v[0] for k, v in _pairs(a, "--haps").items()})
    if cmd == "decide":
        return cmd_decide(*a[:3])
    if cmd == "arms":
        return cmd_arms(*a[:5])
    if cmd == "score":
        return cmd_score(*a[:7])
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
