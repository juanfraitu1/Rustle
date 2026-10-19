#!/usr/bin/env python3
"""A truth for AS-tied reads from the donor's own genome (2026-10-18; user: "we need to prove it in that way ... it just needs to
assign reads to copies based on PSVs for tied AS reads"). HG002: PacBio Kinnex reads and the donor's diploid T2T assembly
(hg002v1.1, maternal and paternal haplotypes), with chains to CHM13 v2.0.

A read is TIED in CHM13 when its best alignment score is shared by alignments on two or more copies (copy_assign's exact AS
tie); copies by exon overlap (realreads.copy_exons). The user's split (2026-10-18: "O2 can only work for tied AS reads,
everything else is already resolved by minimap, then O3 is just an extension in which if we notice the PSVs nothing in the
genome explains them"): a read whose best-AS alignments all lie on one copy is RESOLVED by minimap2 (the premise, tested on
the untied reads with the same truth); a tied read is O2's, which may only choose among the copies it ties on (scoped). Its truth comes from its alignments to the donor's own haplotypes: the
best-scoring ones, lifted to CHM13 through the chain, name the copy --
  unique          every best placement (either haplotype: the copy's two alleles) lifts onto one CHM13 copy -> the truth
  own_tie         they lift onto two or more copies: the donor's copies tie too -> undecidable from the read
  elsewhere       they lift onto CHM13, but onto no copy of the family: the read comes from another CHM13 locus (whose CHM13
                  allele, unlike the donor's, loses to the tied copies) -> never a copy to assign
  no_counterpart  none lifts onto CHM13 at all: a copy of the donor's that CHM13 lacks -> O3's case, never a copy to assign
  unmapped        no alignment to the donor's genome
  poor_match      the donor's best placement is no near-perfect match (NM > max(2, 0.3% of its aligned bases), arm 4's
                  floor, bench/o2_truth/haplotype_truth.py): the read's sequence is not in the donor's assembly as such
The donor margin of a unique read (reported, no class) is how far the best placement off its truth copy falls behind.

    hg002truth.py tied BAM FAMILIES.json SET OUT.tsv           -> the reads tied in BAM (CHM13 coordinates): read, AS, copies,
                                                                  and minimap2's primary (its coin toss among them) + MAPQ
    hg002truth.py untied BAM FAMILIES.json SET OUT.tsv         -> the reads minimap2 resolves, in the same format (copies = the
                                                                  one copy its best-AS alignments lie on)
    hg002truth.py truth TIED.tsv MAT.bam PAT.bam CHAIN FAMILIES.json SET OUT.tsv -> each tied read's class and truth copy
    hg002truth.py world TRUTH.tsv FAMILIES.json SET MIN OUT_DIR [--par keep|exclude|only]
                                                               -> OUT_DIR/families.json (set `hg002`: the families with >= MIN
                                                                  scored tied reads; pseudoautosomal ones kept, excluded or
                                                                  alone) and OUT_DIR/world_truth.tsv (the reads with a donor
                                                                  copy, realreads' truth format)
    hg002truth.py bars O2PROOF.json SCORE.json OUT.json       -> the pre-registered bars T1-T5 and each system's O2 verdict
    hg002truth.py score TRUTH.tsv TIED.tsv WORLD OUT.json [--o2 NAME]  -> per method (minimap2's primary, the primary at MAPQ
                                                                  >= 1, current, integrated: summary.world_calls with the O2
                                                                  table NAME, o2.integv.m2.tsv by default) and truth class,
                                                                  the fates (tied_fate) of the world's tied reads
"""
import collections
import csv
import importlib.util
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("o3_graph_realreads_h", os.path.join(HERE, "realreads.py"))
R = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(R)


METHODS = ("mm2", "mm2_q1", "current", "integrated", "current_tied", "integrated_tied")
SYSTEMS = ("integrated", "current", "integrated_tied", "current_tied")
FLOOR, T1_MAX, T3_MAX, T4_MIN, T5_MAX, VALID_MIN = 100, 0.02, 0.05, 0.90, 0.05, 0.25   # PREREG_o2_hg002_tied_2026-10-18
PREMISE_MIN = 0.99
SCORED = ("unique", "unique_elsewhere", "own_tie", "no_counterpart", "elsewhere")


def _summary():
    spec = importlib.util.spec_from_file_location("o3_graph_summary_h", os.path.join(HERE, "summary.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _copy_of(chrom, blocks, idx):
    return R.truth_copy(chrom, blocks, 0, idx, min_mapq=0)


def tied_copies(alns, idx):
    """(best AS, the copies its best-scoring alignments lie on) of one read's alignments [(chrom, blocks, AS, ...)]; tied
    when two or more."""
    if not alns:
        return None, set()
    best = max(a[2] for a in alns)
    return best, {c for ch, bl, sc, *_ in alns if sc == best for c in [_copy_of(ch, bl, idx)] if c}


def resolved_copy(alns, idx):
    """minimap2's resolution of a read: the copy its best-AS alignments all lie on; None when they tie across copies or
    with a locus off the copies (or lie off them)."""
    if not alns:
        return None
    best = max(a[2] for a in alns)
    cps = {_copy_of(ch, bl, idx) for ch, bl, sc, *_ in alns if sc == best}
    return next(iter(cps)) if len(cps) == 1 and None not in cps else None


def scoped(call, tied):
    """The user's O2 scope: a tied read's call stands only on one of the copies it ties on; otherwise it is an abstention."""
    return call if call in tied else None


def premise_bar(rows):
    """The premise that minimap2 resolves the untied reads: its copy is the donor's (`unique`) rather than another
    (`unique_elsewhere`) in >= PREMISE_MIN of the untied reads with a donor copy."""
    n = collections.Counter(r["class"] for r in rows)
    k = n["unique"] + n["unique_elsewhere"]
    return dict(value=n["unique"] / k if k else None, n=k, powered=k >= FLOOR,
                holds=(n["unique"] / k >= PREMISE_MIN) if k else None, classes=dict(n))


def near_perfect(nm, aligned):
    """Arm 4's floor for a donor placement to count as the read's own locus: NM <= max(2, 0.003 x aligned bases)."""
    return nm <= max(2, 0.003 * aligned)


OFF = "off"                                                   # a best placement lifting onto CHM13 outside every copy


def _placed_copy(contig, bl, lift, idx):
    """The CHM13 copy a donor placement lifts onto (the one its lifted blocks overlap most), OFF when it lifts onto no copy,
    None when it does not lift."""
    pieces = [(ch, s2, e2) for s, e in bl for ch, s2, e2 in lift(contig, s, e)]
    cps = collections.Counter()
    for ch, s2, e2 in pieces:
        c = _copy_of(ch, [(s2, e2)], idx)
        if c:
            cps[c] += e2 - s2
    return max(sorted(cps), key=lambda c: cps[c]) if cps else (OFF if pieces else None)


def donor_margin(own, lift, idx, truth):
    """The best AS minus the best AS of a placement that does not lift onto the truth copy (another copy, CHM13 outside
    the copies, or nowhere); None when every placement lifts onto it."""
    rivals = [sc for contig, bl, sc, *_ in own if _placed_copy(contig, bl, lift, idx) != truth]
    return max(a[2] for a in own) - max(rivals) if rivals else None


def own_truth(own, lift, idx):
    """(class, truth copy) of a read from its alignments to the donor's haplotypes [(contig, blocks, AS)]; lift(contig, s,
    e) -> [(chrom, s, e)] on CHM13."""
    if not own:
        return "unmapped", None
    best = max(a[2] for a in own)
    found = {_placed_copy(contig, bl, lift, idx) for contig, bl, sc, *_ in own if sc == best}
    real = found - {None, OFF}
    if not real:
        return ("elsewhere" if OFF in found else "no_counterpart"), None
    if len(found) > 1:
        return "own_tie", None
    return "unique", next(iter(real))


def tied_fate(cls, truth, assigned, in_o3, tied):
    """A tied read's fate under one method (an absent-copy call takes precedence over a placement, as summary.read_fate):
    unique / unique_elsewhere: right only on the donor's copy, and an O3 call is an error (its copy is in CHM13); own_tie:
    placing it on one of the copies it ties on is a guess, anywhere else or an O3 call an error, abstaining right;
    no_counterpart (a donor copy CHM13 lacks): any placement is an error, an O3 call right; elsewhere (another CHM13 locus):
    a placement on the family is an error, an O3 call is reported apart."""
    if in_o3:
        return {"no_counterpart": "correct", "elsewhere": "o3"}.get(cls, "error")
    if assigned is None:
        return "abstain"
    if cls in ("unique", "unique_elsewhere"):
        return "correct" if assigned == truth else "error"
    return "guess" if cls == "own_tie" and assigned in tied else "error"


def _alns(bam_path, primaries=None):
    """{read: [(contig, blocks, AS, NM, aligned bases)]} of every primary and secondary alignment; primaries, a dict, gets
    {read: (contig, blocks, MAPQ)} of each primary."""
    import pysam
    out = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_supplementary:
                continue
            aligned = sum(n for op, n in r.cigartuples if op in (0, 7, 8))
            out[r.query_name].append((r.reference_name, r.get_blocks(), r.get_tag("AS"),
                                     r.get_tag("NM") if r.has_tag("NM") else 10 ** 6, aligned))   # no NM: never near-perfect
            if primaries is not None and not r.is_secondary:
                primaries[r.query_name] = (r.reference_name, r.get_blocks(), r.mapping_quality)
    return out


def cmd_untied(bam_p, fam_p, set_name, out_p):
    (fams,) = [json.load(open(fam_p))["sets"][set_name]]
    idx = R.copy_exons(fams)
    fam_of = {f"{f['family']}.{k}": f["family"] for f in fams for k in range(len(f["copies"]))}
    n, prim = collections.Counter(), {}
    with open(out_p, "w") as fh:
        fh.write("read\tas\tfamily\tcopies\tprimary\tmapq\n")
        for q, al in _alns(bam_p, prim).items():
            n["reads"] += 1
            cp = resolved_copy(al, idx)
            if cp:
                ch, bl, mq = prim[q]
                fh.write(f"{q}\t{max(a[2] for a in al)}\t{fam_of[cp]}\t{cp}\t{_copy_of(ch, bl, idx) or ''}\t{mq}\n")
                n["resolved"] += 1
    print(f"[hg002truth] {dict(n)} -> {out_p}", flush=True)
    return 0


def cmd_tied(bam_p, fam_p, set_name, out_p):
    (fams,) = [json.load(open(fam_p))["sets"][set_name]]
    idx = R.copy_exons(fams)
    fam_of = {f"{f['family']}.{k}": f["family"] for f in fams for k in range(len(f["copies"]))}
    n, prim = collections.Counter(), {}
    with open(out_p, "w") as fh:
        fh.write("read\tas\tfamily\tcopies\tprimary\tmapq\n")
        for q, al in _alns(bam_p, prim).items():
            best, cps = tied_copies(al, idx)
            n["reads"] += 1
            fs = {fam_of[c] for c in cps}
            if len(cps) >= 2 and len(fs) == 1:
                ch, bl, mq = prim[q]
                fh.write(f"{q}\t{best}\t{next(iter(fs))}\t{','.join(sorted(cps))}\t{_copy_of(ch, bl, idx) or ''}\t{mq}\n")
                n["tied"] += 1
            elif len(cps) >= 2:
                n["tied_across_families"] += 1
    print(f"[hg002truth] {dict(n)} -> {out_p}", flush=True)
    return 0


def cmd_truth(tied_p, mat_p, pat_p, chain_p, fam_p, set_name, out_p):
    spec = importlib.util.spec_from_file_location("o3_graph_chainlift_h", os.path.join(HERE, "chainlift.py"))
    cl = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cl)
    ch = cl.Chain(chain_p)
    fams = json.load(open(fam_p))["sets"][set_name]
    idx = R.copy_exons(fams)
    lift = lambda contig, s, e: [(q, a, b) for q, a, b, _ in ch.lift(contig, s, e)]  # noqa: E731
    own = _alns(mat_p)
    for q, v in _alns(pat_p).items():
        own[q] += v
    n = collections.Counter()
    with open(out_p, "w") as fh:
        fh.write("read\tfamily\ttied_copies\tclass\ttruth\tdonor_nm\tdonor_margin\n")
        for r in csv.DictReader(open(tied_p), delimiter="\t"):
            al = own.get(r["read"], [])
            cls, cp = own_truth(al, lift, idx)
            nm, margin = "", ""
            if al:
                top = max(a[2] for a in al)
                b = min((a for a in al if a[2] == top), key=lambda a: a[3])
                nm = b[3]
                if not near_perfect(b[3], b[4]):
                    cls, cp = "poor_match", None
            if cls == "unique":
                m = donor_margin(al, lift, idx, cp)
                margin = "" if m is None else m
                if cp not in r["copies"].split(","):
                    cls = "unique_elsewhere"                # the donor places it on a copy CHM13 did not tie it with
            n[cls] += 1
            fh.write(f"{r['read']}\t{r['family']}\t{r['copies']}\t{cls}\t{cp or ''}\t{nm}\t{margin}\n")
    print(f"[hg002truth] {dict(n)} -> {out_p}", flush=True)
    return 0


PAR = {"chrX": ((0, 2394410), (153925834, 154259566)), "chrY": ((0, 2458320), (62122809, 62460029))}   # CHM13 v2.0 (T2T-Y)


def pseudoautosomal(fam):
    """Every copy inside a pseudoautosomal region: in a male donor the chrX and chrY copies are the maternal and paternal
    alleles of one recombining locus, so telling them apart is phasing, not paralog assignment."""
    return all(any(s <= c["start"] and c["end"] <= e for s, e in PAR.get(c["chrom"], ())) for c in fam["copies"])


def select_families(fams, truth_rows, min_scored, par="keep"):
    """The families (in pool order) with >= min_scored tied reads of a scored class; par: keep, exclude or only the
    pseudoautosomal ones."""
    n = collections.Counter(r["family"] for r in truth_rows if r["class"] in SCORED)
    return [f for f in fams if n[f["family"]] >= min_scored
            and (par == "keep" or pseudoautosomal(f) == (par == "only"))]


def world_truth(truth_rows):
    """[(read, copy, family)] of the tied reads the donor places on one copy (unique, unique_elsewhere)."""
    return [(r["read"], r["truth"], r["family"]) for r in truth_rows if r["class"] in ("unique", "unique_elsewhere")]


def cmd_world(truth_p, fam_p, set_name, min_scored, out_dir, par="keep"):
    d = json.load(open(fam_p))
    rows = list(csv.DictReader(open(truth_p), delimiter="\t"))
    keep = select_families(d["sets"][set_name], rows, min_scored, par)
    names = {f["family"] for f in keep}
    os.makedirs(out_dir, exist_ok=True)
    json.dump(dict({k: v for k, v in d.items() if k != "sets"}, sets={"hg002": keep}), open(f"{out_dir}/families.json", "w"))
    wt = [t for t in world_truth(rows) if t[2] in names]
    with open(f"{out_dir}/world_truth.tsv", "w") as fh:
        fh.write("read\tcopy\tfamily\n" + "".join("\t".join(t) + "\n" for t in wt))
    print(f"[hg002truth] {len(keep)} of {len(d['sets'][set_name])} families, {sum(len(f['copies']) for f in keep)} copies; "
          f"{len(wt)} reads with a donor copy -> {out_dir}", flush=True)
    return 0


def _fam_call(calls, family):
    """A system's placement of a read of `family`: None, its copy, `no-truth-copy` or `wrong-family` (as o2proof)."""
    if not calls:
        return None
    if family in calls:
        return calls[family] or "no-truth-copy"
    return "wrong-family"


def _world_reads(fq):
    with open(fq) as fh:
        return {h[1:].split()[0] for i, h in enumerate(fh) if i % 4 == 0}


def cmd_score(truth_p, tied_p, world, out_p, o2_name="o2.integv.m2.tsv"):
    wc = _summary().world_calls(world, o2_name)
    tied = {r["read"]: r for r in csv.DictReader(open(tied_p), delimiter="\t")}
    reads = _world_reads(f"{world}/reads.E0.fq")
    e2e, o2 = ({m: collections.defaultdict(collections.Counter) for m in METHODS} for _ in range(2))
    o3 = {m: collections.Counter() for m in METHODS}
    skipped = collections.Counter()
    for r in csv.DictReader(open(truth_p), delimiter="\t"):
        q, f, cls = r["read"], r["family"], r["class"]
        if cls in ("unmapped", "poor_match") or q not in reads:
            skipped[cls if cls in ("unmapped", "poor_match") else "not_in_world"] += 1
            continue
        t, ts = tied[q], set(r["tied_copies"].split(","))
        prim = t["primary"] or "wrong-locus"
        calls = dict(mm2=(prim, False), mm2_q1=(prim if int(t["mapq"]) >= 1 else None, False),
                     current=(_fam_call(wc["current"].get(q), f), q in wc["current_o3"]),
                     integrated=(_fam_call(wc["integrated"].get(q), f), q in wc["integrated_o3"]))
        for m in ("current", "integrated"):
            calls[f"{m}_tied"] = (scoped(calls[m][0], ts), calls[m][1])
        for m, (a, in_o3) in calls.items():
            e2e[m][cls][tied_fate(cls, r["truth"] or None, a, in_o3, ts)] += 1
            o2[m][cls][tied_fate(cls, r["truth"] or None, a, False, ts)] += 1
            if in_o3:
                o3[m][cls] += 1
    res = dict(e2e={m: {c: dict(v) for c, v in e2e[m].items()} for m in METHODS},
               o2={m: {c: dict(v) for c, v in o2[m].items()} for m in METHODS},
               o3={m: dict(o3[m]) for m in METHODS}, skipped=dict(skipped))
    json.dump(res, open(out_p, "w"), indent=1)
    for m in METHODS:
        print(f"  {m:10s} " + " | ".join(f"{c}: " + " ".join(f"{k} {v}" for k, v in sorted(res['e2e'][m][c].items()))
                                        for c in sorted(res["e2e"][m])), flush=True)
    print(f"  skipped {dict(skipped)}", flush=True)
    return 0


def wilson(k, n, z=1.96):
    """The Wilson score interval of k successes in n."""
    if n == 0:
        return 0.0, 1.0
    p, d = k / n, 1 + z * z / n
    c, h = (p + z * z / (2 * n)) / d, z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (0.0 if k == 0 else max(0.0, c - h)), (1.0 if k == n else min(1.0, c + h))




def _bar(value, n, ok):
    return dict(value=value, n=n, powered=n >= FLOOR, holds=ok(value) if n else None)


def bars(o2p, score, systems=SYSTEMS):
    """The pre-registered bars from o2proof's pooled card (the world's truth reads) and hg002truth's score: per system T1
    (wrong / placed), T2 (its Wilson upper bound below the coin toss's lower one), T3 (placed own_tie reads), T4 (correct
    decidable reads), T5 (O3 calls on reads whose copy is present) and the O2 verdict; `valid`: mm2 wrong on >= 25%."""
    def tot(card, field):
        return sum(v["all"][field] for v in card.values())
    mm2 = o2p["pooled"]["card"]["mm2"]
    mm2_n = tot(mm2, "correct") + tot(mm2, "wrong") + tot(mm2, "abstain")
    valid = bool(mm2_n) and tot(mm2, "wrong") / mm2_n >= VALID_MIN
    out = dict(valid=valid, mm2_wrong=tot(mm2, "wrong"), mm2_n=mm2_n)
    for m in systems:
        card = o2p["pooled"]["card"][m]
        wrong, placed = tot(card, "wrong"), tot(card, "correct") + tot(card, "wrong")
        b = dict(T1=_bar(wrong / placed if placed else None, placed, lambda v: v <= T1_MAX))
        up, lo = wilson(wrong, placed)[1], wilson(out["mm2_wrong"], mm2_n)[0]
        b["T2"] = dict(value=(up, lo), n=placed, powered=placed >= FLOOR and valid,
                       holds=(up < lo) if placed and valid else None)
        tie = score["o2"][m].get("own_tie", {})
        n_tie = sum(tie.values())
        b["T3"] = _bar((n_tie - tie.get("abstain", 0)) / n_tie if n_tie else None, n_tie, lambda v: v <= T3_MAX)
        dec = card.get("decidable", {}).get("all", dict(correct=0, wrong=0, abstain=0))
        n_dec = sum(dec.values())
        b["T4"] = _bar(dec["correct"] / n_dec if n_dec else None, n_dec, lambda v: v >= T4_MIN)
        present = ("unique", "unique_elsewhere", "own_tie")
        n_p = sum(sum(score["e2e"][m].get(c, {}).values()) for c in present)
        calls = sum(score["o3"][m].get(c, 0) for c in present)
        b["T5"] = _bar(calls / n_p if n_p else None, n_p, lambda v: v <= T5_MAX)
        o2_bars = [b[t] for t in ("T1", "T2", "T3", "T4")]
        if any(x["powered"] and x["holds"] is False for x in o2_bars):
            b["verdict"] = "FAIL"
        elif all(b[t]["powered"] and b[t]["holds"] for t in ("T1", "T2")):
            b["verdict"] = "PASS"
        else:
            b["verdict"] = "UNDERPOWERED"
        out[m] = b
    return out


def cmd_bars(o2p_p, score_p, out_p):
    res = bars(json.load(open(o2p_p)), json.load(open(score_p)))
    json.dump(res, open(out_p, "w"), indent=1)
    print(f"valid (mm2 wrong >= {VALID_MIN}): {res['valid']} ({res['mm2_wrong']} / {res['mm2_n']})")
    for m in SYSTEMS:
        print(f"{m}: O2 {res[m]['verdict']}; " + "; ".join(
            f"{t} {res[m][t]['value']} n {res[m][t]['n']} {'powered' if res[m][t]['powered'] else 'UNDERPOWERED'} "
            f"{'holds' if res[m][t]['holds'] else 'FAILS' if res[m][t]['holds'] is False else '-'}"
            for t in ("T1", "T2", "T3", "T4", "T5")), flush=True)
    return 0


def main(argv):
    if len(argv) == 4 and argv[0] == "bars":
        return cmd_bars(*argv[1:])
    if argv[:1] == ["score"] and (len(argv) == 5 or len(argv) == 7 and argv[5] == "--o2"):
        return cmd_score(*argv[1:5], o2_name=argv[6] if len(argv) == 7 else "o2.integv.m2.tsv")
    if argv[:1] == ["world"] and (len(argv) == 6 or len(argv) == 8 and argv[6] == "--par"):
        return cmd_world(argv[1], argv[2], argv[3], int(argv[4]), argv[5], argv[7] if len(argv) == 8 else "keep")
    if len(argv) == 5 and argv[0] == "tied":
        return cmd_tied(*argv[1:])
    if len(argv) == 5 and argv[0] == "untied":
        return cmd_untied(*argv[1:])
    if len(argv) == 8 and argv[0] == "truth":
        return cmd_truth(*argv[1:])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
