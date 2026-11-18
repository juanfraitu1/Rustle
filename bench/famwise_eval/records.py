"""PREREG r1796 (docs/PREREG_famwise_origin_2026-11-18.md): the shared reader of bench/famwise_eval.
A truth table of bench/o2_truth/neutral_truth.py (optionally with the tie split) and a selected-reads BAM become
per-molecule records carrying both scores (ms and AS). No path is fixed here; the scripts beside this one take them."""
import collections, csv, os, sys
import pysam
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from famwise import record_chain

Truth = collections.namedtuple("Truth", "stratum status kind origin loci")
Rec = collections.namedtuple("Rec", "chrom start end strand chain ms AS primary")
JUDGED = ("one", "several")            # kinds with a verdict on the origin
KEEP = 0.8                             # neutral_truth.py judges placements at AS >= 0.8 x the best
SPLIT = {"SEVERAL COPIES": "several", "ONE PLACE": "one_place"}

def load_truth(reads_tsv, tie_split=None):
    """-> {read: Truth}. loci = [(chrom, start, end, strand, AS)], 0-based half-open, best AS first; origin = indexes."""
    cls = {}
    if tie_split:
        with open(tie_split) as fh:
            cls = {r["read"]: r["class"] for r in csv.DictReader(fh, delimiter="\t")}
    out = {}
    with open(reads_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            loci = []
            for item in filter(None, r["loci"].split(";")):
                _, chrom, span, strand, a = item.split(":")
                s, e = span.split("-")
                loci.append((chrom, int(s), int(e), strand, int(a[2:])))
            origin = frozenset(int(x) for x in r["origin_loci"].split(",") if x)
            kind = {"decidable": "one", "outside": "outside"}.get(r["status"], "unjudged")
            if r["status"] == "tie":
                kind = SPLIT.get(cls.get(r["read"]), "tie_unsplit")
            if kind not in ("outside", "unjudged") and (not origin or max(origin) >= len(loci)):
                kind = "unjudged"      # an origin beyond the 12 written loci
            out[r["read"]] = Truth(r["stratum"], r["status"], kind, origin, loci)
    return out

def load_records(bam):
    """-> {read: [Rec]}: every mapped, non-supplementary record (the selected-reads BAMs are unsorted, no index)."""
    out = collections.defaultdict(list)
    with pysam.AlignmentFile(bam, check_sq=False) as fh:
        for r in fh.fetch(until_eof=True):
            if r.is_unmapped or r.is_supplementary:
                continue
            out[r.query_name].append(Rec(r.reference_name, r.reference_start, r.reference_end, "-" if r.is_reverse else "+",
                                         record_chain(r.reference_start, r.cigartuples or ()), r.get_tag("ms"), r.get_tag("AS"), not r.is_secondary))
    return dict(out)

def locus_index(rec, t):
    """index of the truth locus the record lies on (contig, strand, overlap), or None = the truth never judged it:
    no such locus, or one the truth lists below its kept placements (AS < 0.8 x best) that is not an origin."""
    for i, (chrom, s, e, strand, a) in enumerate(t.loci):
        if chrom == rec.chrom and strand == rec.strand and s < rec.end and rec.start < e:
            return i if a >= KEEP * t.loci[0][4] or i in t.origin else None
    return None

def statement(recs, score, window):
    """window None = the primary record alone; else the primary plus every record within `window` of the best `score`."""
    v = lambda r: r.ms if score == "ms" else r.AS
    best = max(map(v, recs))
    return [r for r in recs if r.primary or (window is not None and best - v(r) <= window)]

def places(recs, score):
    """best `score` per place, best first; records on one contig and strand with overlapping spans are one place
    (the truth's own merge, neutral_truth.loci_of)"""
    out = []
    for r in sorted(recs, key=lambda r: -(r.ms if score == "ms" else r.AS)):
        for p in out:
            if p[0] == r.chrom and p[1] == r.strand and p[2] < r.end and r.start < p[3]:
                p[2], p[3] = min(p[2], r.start), max(p[3], r.end)
                break
        else:
            out.append([r.chrom, r.strand, r.start, r.end, r.ms if score == "ms" else r.AS])
    return [p[4] for p in out]

def lead(recs, score):
    """best minus second-best place on `score`; None with one place"""
    v = places(recs, score)
    return v[0] - v[1] if len(v) > 1 else None

def near_tied(recs, score, within=10):
    """the molecule's two best places are within `within` on `score`"""
    d = lead(recs, score)
    return d is not None and d <= within

def wilson(k, n, z=1.96):
    """Wilson 95% interval of k / n -> (lo, hi); (0.0, 1.0) at n = 0"""
    if not n:
        return 0.0, 1.0
    p, d = k / n, 1 + z * z / n
    c, h = (p + z * z / (2 * n)) / d, z * (p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5 / d
    return max(0.0, c - h), min(1.0, c + h)
