#!/usr/bin/env python3
"""Test A of the cross-copy splicing design (docs/superpowers/specs/2026-11-18-cross-copy-splicing-design.md, section 3):
where the table says copy U uses a junction more than copy S, do the splice sites say the same?

    site_cause.py score FAM.psi.tsv PSI.introns.tsv FAM.loci.gated.paf GTF FASTA OUT.tsv
    site_cause.py pool --dev GENOME=DEV.tsv GENOME=OUT.tsv [GENOME=OUT.tsv ...]

Advantage of the key junction at a copy = strength of its site at the end it does not share with its bunch, minus the
strongest such site among its competitors there (the introns of >= 3 reads sharing its anchor site, of either copy, each
carried to the other copy through the pair's alignment). dA = advantage at U - advantage at S; a cause in the sites
predicts dA > 0. A competitor without an image is left out when it holds under 10% of its site's certain reads, else the
call is UNLOCATED; so is a call with a site too near a contig end or holding N.
A junction no read carries at a copy exists there only as an image through the alignment, which can slip by a few bases
beside a gap (development sample: 2 of 32 agreeing calls were a 2-base slip onto an intact site). When either copy's
side of a junction is such an image, the site is read at BOTH copies as the strongest within 6 bases, so identical
sequences give identical strengths. A site lost with another within 6 bases reads as no difference.
MISKEYED (calls.flag_miskeyed): S uses the image of U's junction under another key; counted, no verdict.
`pool` counts a unit (genome, U's junction, S) once (the sample with more certain reads at U), leaves out every unit of
the development sample, and prints bar A1 and reading A2."""
import collections, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "..", "tools"), os.path.join(HERE, "..", "famwise_eval")]
import calls as K, sites as S
from psi_family import load_introns, site_index, site_of, load_pair_records, project_intron

MIN_READS, MIN_SHARE, TAUS, SLIP = 3, 0.10, (0.5, 1.0, 2.0), 6
COLS = ("family anchor col_anchor col_other cls gap u s u_chrom u_strand u_start u_end reads_u copy_specific silent a_u a_s dA "
        "verdict_0.5 verdict_1 verdict_2 canonical_lost n_competitors dropped_competitors pair_identity s_chrom s_strand s_start s_end unit").split()

def advantage(key_strength, competitor_strengths):
    """NaN when the key site or any competitor has no window: the unit is then not located"""
    if any(math.isnan(x) for x in (key_strength, *competitor_strengths)):
        return float("nan")
    return key_strength - max(competitor_strengths) if competitor_strengths else key_strength

def verdict(dA, tau=1.0):
    return "UNLOCATED" if math.isnan(dA) else "AGREES" if dA >= tau else "OPPOSES" if dA <= -tau else "NO SITE DIFFERENCE"

def pair_identity(paf, members):
    """{frozenset(a, b): matching / aligned bases of the pair's record with the most matching bases}"""
    best = {}
    for l in open(paf):
        f = l.split("\t", 11)
        if f[0] != f[5] and f[0] in members and f[5] in members:
            k = frozenset((f[0], f[5]))
            if k not in best or int(f[9]) > best[k][0]:
                best[k] = (int(f[9]), int(f[10]))
    return {k: m / n for k, (m, n) in best.items()}

def compared(c, dump, idx, rec):
    """([(intron at U, intron at S)], dropped, unlocated): the competitors of the key junction at either copy, each at both"""
    pairs, dropped, unlocated = [], 0, False
    for src, dst, key, other, flip in ((c.u, c.s, c.u_intron, c.s_intron, False), (c.s, c.u, c.s_intron, c.u_intron, True)):
        at = idx[site_of(key, c.anchor)]
        total = sum(dump[i][2] for i in at) or 1
        for i in sorted(at):
            if i == key or dump[i][3] < MIN_READS or any(p[flip] == i for p in pairs):
                continue
            j = project_intron(rec, src, dst, i[2], i[3])
            if j is None:
                if dump[i][2] / total < MIN_SHARE: dropped += 1
                else: unlocated = True
                continue
            j = (other[0], other[1], *j)
            if j != other and not any(p[not flip] == j for p in pairs):
                pairs.append((j, i) if flip else (i, j))
    return pairs, dropped, unlocated

def score(psi, dump_p, paf, gtf, fasta, out):
    import pysam
    C = K.load_calls(psi)
    dump = load_introns(dump_p); idx = site_index(dump)
    members = {m: i[1] for c in C for m, i in ((c.u, c.u_intron), (c.s, c.s_intron))}
    rec, ident = load_pair_records(paf, members), pair_identity(paf, members)
    C = K.flag_miskeyed(C, dump, rec)
    pw, fa = S.sample_pwms(gtf, fasta), pysam.FastaFile(fasta)
    unread = lambda i: not dump.get(i, (0, 0, 0, 0))[3]
    def site(i, anchor, search):
        """(strength, canonical) of the intron's end away from its anchor; with `search`, of the strongest site within SLIP bases"""
        first = (anchor == "donor") == (i[1] == "-")      # that end is the intron's first base on the reference
        out = []
        for d in (sorted(range(-SLIP, SLIP + 1), key=abs) if search else [0]):
            j = (i[0], i[1], i[2] + d, i[3]) if first else (i[0], i[1], i[2], i[3] + d)
            x = S.other_end(fa, j, anchor)
            if d == 0 and x is None:
                return float("nan"), False                # the position itself has no window: not located
            if x is not None:
                out.append((S.score(pw[1] if anchor == "donor" else pw[0], x), S.canonical(x, anchor)))
        return max(out)                                   # canonical = of the strongest site found
    def both(pair, anchor):
        """the two copies' sites of one junction, read by one rule: the search runs on both when either is only an image"""
        search = unread(pair[0]) or unread(pair[1])
        return site(pair[0], anchor, search), site(pair[1], anchor, search)
    rows = []
    for c in C:
        pairs, dropped, unlocated = compared(c, dump, idx, rec)
        key = both((c.u_intron, c.s_intron), c.anchor); comp = [both(p, c.anchor) for p in pairs]
        a_u = advantage(key[0][0], [x[0][0] for x in comp]); a_s = advantage(key[1][0], [x[1][0] for x in comp])
        dA = float("nan") if unlocated else a_u - a_s
        vs = ["MISKEYED"] * len(TAUS) if c.miskeyed else [verdict(dA, t) for t in TAUS]
        rows.append([*c[:8], *c.u_intron, c.reads_u, int(c.copy_specific), int(c.silent), f"{a_u:.3f}", f"{a_s:.3f}", f"{dA:.3f}",
                     *vs, int(key[0][1] and not key[1][1]), len(pairs), dropped,
                     f"{ident.get(frozenset((c.u, c.s)), float('nan')):.4f}", *c.s_intron, K.unit_name(c)])
    with open(out + ".tmp", "w") as fh:
        fh.write("\t".join(COLS) + "\n"); fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    os.replace(out + ".tmp", out)
    n = collections.Counter((r[4], r[19]) for r in rows)
    for cls in ("LARGE", "SMALL", "SAME"):
        print(cls, *(f"{v} {n[(cls, v)]}" for v in ("AGREES", "OPPOSES", "NO SITE DIFFERENCE", "UNLOCATED", "MISKEYED")), sep="\t")

def read(path):
    lines = open(path).read().splitlines(); h = lines[0].split("\t")
    return [dict(zip(h, l.split("\t"))) for l in lines[1:]]

def events(rows):
    """LARGE units that share U's junction or S's anchor site on a genome are one event (a key gives one call, between its two
    extreme copies, so samples with other copies state one event as different units); the unit with the largest gap speaks.
    -> (events, AGREES, OPPOSES at tau 1)"""
    root = {}
    def find(x):
        while root.setdefault(x, x) != x:
            x = root[x]
        return x
    for r in rows:
        u, s = r["unit"].split("|")
        root[find((r["genome"], "S", s))] = find((r["genome"], "U", u))
    best = {}
    for r in rows:
        e = find((r["genome"], "U", r["unit"].split("|")[0]))
        if e not in best or float(r["gap"]) > float(best[e]["gap"]):
            best[e] = r
    v = collections.Counter(r["verdict_1"] for r in best.values())
    return len(best), v["AGREES"], v["OPPOSES"]

def pool(dev, samples):
    """dev = (genome, DEV.tsv); samples = [(genome, OUT.tsv)]. -> dict of counts; units of the development sample left out."""
    uid = lambda g, r: (g, r["unit"])
    seen = {uid(dev[0], r) for r in read(dev[1])}
    units = {}
    for g, path in samples:
        for r in read(path):
            k = uid(g, r)
            if k not in seen and (k not in units or int(r["reads_u"]) > int(units[k]["reads_u"])):
                units[k] = dict(r, genome=g)
    res = {"units": units}
    L = [r for r in units.values() if r["cls"] == "LARGE"]; same = [r for r in units.values() if r["cls"] == "SAME"]
    diff = lambda rs, col="verdict_1": sum(r[col] in ("AGREES", "OPPOSES") for r in rs)
    loc = lambda rs: [r for r in rs if r["verdict_1"] not in ("UNLOCATED", "MISKEYED")]
    res.update(units_large=len(L), agrees=sum(r["verdict_1"] == "AGREES" for r in L), opposes=sum(r["verdict_1"] == "OPPOSES" for r in L),
               site_diff_large=(diff(L), len(loc(L))), site_diff_same=(diff(same), len(loc(same))), miskeyed=sum(r["verdict_1"] == "MISKEYED" for r in L), events=events(L))
    return res

def report(res):
    from records import wilson
    sh = lambda k, n: f"{k} / {n}" + (f" = {k / n:.3f} [{wilson(k, n)[0]:.3f}, {wilson(k, n)[1]:.3f}]" if n else "")
    U = list(res["units"].values()); L = [r for r in U if r["cls"] == "LARGE"]
    a, o = res["agrees"], res["opposes"]
    print(f"units {len(U)} (LARGE {len(L)}, SAME {sum(r['cls'] == 'SAME' for r in U)}, SMALL {sum(r['cls'] == 'SMALL' for r in U)}); development units left out")
    low = wilson(a, a + o)[0] if a + o else float("nan")
    print(f"A1  LARGE units with a site difference at tau 1: AGREES {sh(a, a + o)}  ->  "
          + ("undecidable (under 30 units)" if a + o < 30 else "PASS" if low > 0.5 else "FAIL"))
    n, ea, eo = res["events"]
    print(f"    beside, linked units counted once ({n} events; units sharing U's junction or S's anchor site, the largest gap speaks): AGREES {sh(ea, ea + eo)}")
    for t in ("0.5", "2"):
        v = collections.Counter(r[f"verdict_{t}"] for r in L)
        print(f"    tau {t}: AGREES {sh(v['AGREES'], v['AGREES'] + v['OPPOSES'])}")
    print(f"    canonical dinucleotide at U and not at S: {sum(r.get('canonical_lost') == '1' for r in L)} of {len(L)} LARGE units; UNLOCATED {sum(r['verdict_1'] == 'UNLOCATED' for r in L)}; MISKEYED {res['miskeyed']}")
    print(f"A2  share with a site difference (located units): LARGE {sh(*res['site_diff_large'])}; SAME {sh(*res['site_diff_same'])}")
    for name, f in (("pair identity >= 0.98", lambda r: float(r["pair_identity"]) >= 0.98), ("pair identity < 0.98", lambda r: float(r["pair_identity"]) < 0.98)):
        parts = []
        for cls in ("LARGE", "SAME"):
            rs = [r for r in U if r["cls"] == cls and r["verdict_1"] not in ("UNLOCATED", "MISKEYED") and f(r)]
            parts.append(f"{cls} {sh(sum(r['verdict_1'] in ('AGREES', 'OPPOSES') for r in rs), len(rs))}")
        print(f"    {name}: " + "; ".join(parts))
    print("per genome (LARGE): genome, units, AGREES, OPPOSES, NO SITE DIFFERENCE, UNLOCATED, MISKEYED")
    for g in sorted({r["genome"] for r in L}):
        v = collections.Counter(r["verdict_1"] for r in L if r["genome"] == g)
        print(f"    {g}\t{sum(v.values())}\t{v['AGREES']}\t{v['OPPOSES']}\t{v['NO SITE DIFFERENCE']}\t{v['UNLOCATED']}\t{v['MISKEYED']}")
    fam = collections.Counter((r["genome"], r["family"]) for r in L)
    top = fam.most_common(3)
    print(f"LARGE units in the three largest families: {sum(n for _, n in top)} of {len(L)} ({', '.join(f'{g}:{f} {n}' for (g, f), n in top)})")

def main(argv):
    if argv[0] == "score":
        score(*argv[1:7])
    else:
        assert argv[1] == "--dev", __doc__
        report(pool(tuple(argv[2].split("=", 1)), [tuple(a.split("=", 1)) for a in argv[3:]]))

if __name__ == "__main__":
    main(sys.argv[1:])
