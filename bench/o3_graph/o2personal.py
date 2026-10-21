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
MARGIN_T = 5                                                    # one substitution under splice:hq (A 1 + B 4)


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
    return {p: sum(1 for _, bl in reads if _covers(bl, p)) for p in positions}


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


def truth_of(alns, margin=MARGIN_T):
    """A read's origin in the animal's own diploid genome from its alignments [(AS, reference copy or None, 'ref'|'alt')]:
    the one copy every alignment within `margin` of the best lies on -> ('determined', copy, 'refcons' when the reference's
    own copy is among them, else 'altonly'); two copies (or a copy and an off-copy locus) -> 'undetermined'; the best off
    every copy -> 'off-copy'."""
    if not alns:
        return "unmapped", None, None
    best = max(a for a, _, _ in alns)
    w = [(a, c, r) for a, c, r in alns if a > best - margin]
    cs = {c for _, c, _ in w}
    if len(cs) > 1:
        return "undetermined", None, None
    c = cs.pop()
    if c is None:
        return "off-copy", None, None
    return "determined", c, "refcons" if any(r == "ref" for _, _, r in w) else "altonly"


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
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
