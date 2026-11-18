#!/usr/bin/env python3
"""Test C of the cross-copy splicing design (docs/superpowers/specs/2026-11-18-cross-copy-splicing-design.md, section 4):
do the reads behind a call come from the copy they are counted at?

    origin.py loci FAM.psi.tsv OUT.bed          the member loci of LARGE and SAME calls (family in column 4): the reads to give a truth
    origin.py score FAM.psi.tsv COLUMNS.tsv FAM.loci.gff3 PSI.introns.tsv SEL.bam TRUTH.tsv TIE_SPLIT.tsv OUT_PREFIX [SITES.tsv [DEV.sites.tsv]]

SITES = site_cause.py's table of the sample: its MISKEYED calls are left out. DEV.sites.tsv = that table of the development
sample, when it sits on the same genome: calls that are units there are left out.

TRUTH / TIE_SPLIT = bench/o2_truth/neutral_truth.py `score` and tie_split.py on the reads of SEL.bam (every record of each).
A call's bunch at copy X = the introns of >= 3 reads sharing the anchor site of X's key junction, plus that junction.
A read counts at X when every record of its wide band (primary + records within 10 `ms` of its best) carries one intron
of the bunch: the `placed_w` of tools/psi_bunches.py. Its verdict, from the truth:
    SUPPORTED  one origin copy, and it is X (same contig and strand, an exon base shared with X's column exons)
    AGAINST    one origin copy elsewhere, or several copies fit equally (tie, class SEVERAL COPIES)
    OUTSIDE    its home is none of its reference places (a copy the reference lacks): printed beside, in no denominator
    none       every other truth kind, or a read the truth run did not score
C1 (reads): SUPPORTED over SUPPORTED + AGAINST on the bunches of LARGE calls, each read once per copy and site.
C2 (calls): a LARGE call is judged with >= 20 SUPPORTED reads in the bunch at U and at S; origin PSI = those on the copy's
introns under the key (the ones the stated range sums) over all of them; HOLDS when origin PSI(U) - origin PSI(S) >= half the call's gap.
Stdout: counts first (step 0), then C1, the OUTSIDE share, C2; SAME calls' reads beside."""
import collections, os, sys, types
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "..", "tools"), os.path.join(HERE, "..", "famwise_eval")]
import calls as K
from psi_family import load_introns, site_index, site_of, load_columns
from records import load_truth, load_records, statement, wilson

MIN_READS, MIN_SUPPORTED, C1_SHARE, C1_LOW, C2_SHARE, C2_CALLS = 3, 20, 0.95, 0.90, 0.90, 20

def placed_introns(recs):
    """the introns (chrom, strand, first base, last base; 1-based) every record of the read's wide band carries"""
    return set(frozenset.intersection(*(frozenset((r.chrom, r.strand, a + 1, b - 1) for a, b in r.chain) for r in statement(recs, "ms", 10))))

def read_verdict(truth, member_exons):
    """truth: records.Truth or None; member_exons: [(chrom, strand, start, end)] of copy X, 1-based closed"""
    if truth is None or truth.kind not in ("one", "several", "outside"):
        return None
    if truth.kind != "one":
        return "OUTSIDE" if truth.kind == "outside" else "AGAINST"
    at = lambda l: any(c == l[0] and st == l[3] and l[1] < e and s - 1 < l[2] for c, st, s, e in member_exons)
    return "SUPPORTED" if any(at(truth.loci[i]) for i in truth.origin) else "AGAINST"

def table(path):
    lines = open(path).read().splitlines(); h = lines[0].split("\t")
    return [dict(zip(h, l.split("\t"))) for l in lines[1:]]

def score(psi, columns, gff, dump_p, bam, truth_p, tie_split, out, sites=None, dev=None):
    C = [c for c in K.load_calls(psi) if c.cls in ("LARGE", "SAME")]
    drop = {(r["family"], r["anchor"], r["col_anchor"], r["col_other"]) for r in (table(sites) if sites else ()) if r["verdict_1"] == "MISKEYED"}
    seen = {r["unit"] for r in (table(dev) if dev else ())}
    left = len(C); C = [c for c in C if c[:4] not in drop and K.unit_name(c) not in seen]; left -= len(C)
    dump = load_introns(dump_p); idx = site_index(dump)
    exons = collections.defaultdict(list)
    for (chrom, strand), ivs in load_columns(columns, gff).items():
        for s, e, _, mem, _, _ in ivs:
            exons[mem].append((chrom, strand, s, e))
    where = collections.defaultdict(list)       # intron -> [(call index, role)]
    for n, c in enumerate(C):
        for role, key in (("U", c.u_intron), ("S", c.s_intron)):
            for i in {key, *(i for i in idx[site_of(key, c.anchor)] if dump[i][3] >= MIN_READS)}:
                where[i].append((n, role))
    truth, rows = load_truth(truth_p, tie_split), []
    for read, recs in load_records(bam).items():
        for i in placed_introns(recs):
            for n, role in where.get(i, ()):
                c, t = C[n], truth.get(read)
                rows.append((read, n, role, i, t.kind if t else "none", read_verdict(t, exons[c.u if role == "U" else c.s]) or "none"))
    with open(out + ".reads.tsv.tmp", "w") as fh:
        fh.write("read\tfamily\tcall\tcls\trole\tintron\tkind\tverdict\n")
        fh.writelines(f"{r}\t{C[n].family}\t{C[n].anchor}:{C[n].col_anchor}:{C[n].col_other}\t{C[n].cls}\t{role}\t{i[0]}:{i[1]}:{i[2]}-{i[3]}\t{k}\t{v}\n" for r, n, role, i, k, v in rows)
    os.replace(out + ".reads.tsv.tmp", out + ".reads.tsv")
    # C1: each read once per copy and anchor site
    once = {}
    for r, n, role, i, k, v in rows:
        c = C[n]; key = c.u_intron if role == "U" else c.s_intron
        once[(r, site_of(key, c.anchor), c.cls)] = (c.family, role, k, v)
    tab = lambda cls, f=lambda x: True: collections.Counter(v for (_, _, cl), (fam, role, k, v) in once.items() if cl == cls and f((fam, role, k)))
    L = tab("LARGE")
    calls = []
    for n, c in enumerate(C):
        if c.cls != "LARGE":
            continue
        sup = {role: [i for r, m, ro, i, k, v in rows if m == n and ro == role and v == "SUPPORTED"] for role in "US"}
        nu, ns, ku, ks = len(sup["U"]), len(sup["S"]), sum(i in c.u_set for i in sup["U"]), sum(i in c.s_set for i in sup["S"])
        judged = nu >= MIN_SUPPORTED and ns >= MIN_SUPPORTED
        og = round(ku / nu - ks / ns, 4) if judged else "NA"
        calls.append(dict(family=c.family, anchor=c.anchor, col_anchor=c.col_anchor, col_other=c.col_other, u=c.u, s=c.s, gap=c.gap, copy_specific=int(c.copy_specific),
                          sup_u=nu, key_u=ku, sup_s=ns, key_s=ks, origin_gap=og, verdict="UNJUDGED" if not judged else "HOLDS" if og >= c.gap / 2 else "FAILS"))
    with open(out + ".calls.tsv.tmp", "w") as fh:
        fh.write("\t".join(calls[0] if calls else ["family"]) + "\n"); fh.writelines("\t".join(map(str, d.values())) + "\n" for d in calls)
    os.replace(out + ".calls.tsv.tmp", out + ".calls.tsv")
    res = types.SimpleNamespace(calls=calls, c1=(L["SUPPORTED"], L["SUPPORTED"] + L["AGAINST"]), outside=L["OUTSIDE"], none=L["none"])
    sh = lambda k, n: f"{k} / {n}" + (f" = {k / n:.3f} [{wilson(k, n)[0]:.3f}, {wilson(k, n)[1]:.3f}]" if n else "")
    v = collections.Counter(d["verdict"] for d in calls); j = v["HOLDS"] + v["FAILS"]
    print(f"calls left out (miskeyed, or units of the development sample): {left}")
    print(f"step 0  LARGE calls {len(calls)}; reads on their bunches: judged {res.c1[1]}, outside {res.outside}, no verdict {res.none}; judged calls {j}")
    k, n = res.c1; low = wilson(k, n)[0]
    print(f"C1  SUPPORTED {sh(k, n)}  ->  " + ("undecidable (no judged read)" if not n else "PASS" if k / n >= C1_SHARE and low >= C1_LOW else "FAIL"))
    for name, f in (("kind one", lambda x: x[2] == "one"), ("kind several", lambda x: x[2] == "several"), ("at U", lambda x: x[1] == "U"), ("at S", lambda x: x[1] == "S")):
        t = tab("LARGE", f); print(f"    {name}: {sh(t['SUPPORTED'], t['SUPPORTED'] + t['AGAINST'])}")
    for fam in sorted({x[0] for (_, _, cl), x in once.items() if cl == "LARGE"}):
        t = tab("LARGE", lambda x: x[0] == fam); print(f"    family {fam}: {sh(t['SUPPORTED'], t['SUPPORTED'] + t['AGAINST'])}; outside {t['OUTSIDE']}")
    print(f"    beside: OUTSIDE {sh(res.outside, n + res.outside)} of judged + outside reads")
    t = tab("SAME"); print(f"    beside: SAME calls' bunches SUPPORTED {sh(t['SUPPORTED'], t['SUPPORTED'] + t['AGAINST'])}; outside {t['OUTSIDE']}")
    print(f"C2  HOLDS {sh(v['HOLDS'], j)} of judged LARGE calls (unjudged {v['UNJUDGED']})  ->  "
          + (f"undecidable (under {C2_CALLS} judged calls)" if j < C2_CALLS else "PASS" if v["HOLDS"] / j >= C2_SHARE else "FAIL"))
    cs = [d for d in calls if d["copy_specific"] and d["verdict"] != "UNJUDGED"]
    print(f"    copy-specific calls: HOLDS {sh(sum(d['verdict'] == 'HOLDS' for d in cs), len(cs))}")
    return res

def main(argv):
    if argv[0] == "loci":
        C = [c for c in K.load_calls(argv[1]) if c.cls in ("LARGE", "SAME")]
        span = lambda m: (m.rsplit(":", 1)[0], *map(int, m.rsplit(":", 1)[1].split("-")))
        rows = sorted({(*span(m), c.family) for c in C for m in (c.u, c.s)})
        with open(argv[2] + ".tmp", "w") as fh:
            fh.writelines(f"{chrom}\t{s - 1}\t{e}\t{fam}\n" for chrom, s, e, fam in rows)
        os.replace(argv[2] + ".tmp", argv[2])
    else:
        score(*argv[1:9], *argv[9:11])

if __name__ == "__main__":
    main(sys.argv[1:])
