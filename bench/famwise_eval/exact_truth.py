"""PREREG r1797 / r1798 (docs/PREREG_famwise_sim_2026-11-18.md): the exact read-origin truth of simulated reads, in the
reads.tsv format of bench/o2_truth/neutral_truth.py, so records.load_truth and origin_bam.py read it unchanged.

  exact_truth.py BAM OUT.tsv --origins ORIGIN.tsv [--null] [--cf CF_READS.tsv]
  exact_truth.py BAM OUT.tsv --gtf TRUTH.gtf [--name-re REGEX] [--null] [--cf CF_READS.tsv]

BAM        every record of the reads (any order, AS tag). Loci, their order (best AS first), strata (on AS: S0 one locus,
           S1 lead > 10, S2 lead <= 10) and the primary pick are the truth scorer's own (neutral_truth.records / loci_of).
--origins  header line (first field "read") with the columns read, chrom, strand, blocks ("s-e,s-e": 0-based half-open
           exon blocks of the read's source gene on the reference); other columns ignored; one row per read.
--gtf      the source transcript is named in the read name (--name-re, group 1), its exons are the GTF's (row 1799).
Origin     the read's locus on the source's contig and strand whose own blocks share a base with a source exon; of several
           the best-AS one. status decidable; no such locus: outside (no home record).
--null     the origin is the locus of the primary record, whatever the source (input N of row 1798).
OUT.tsv    the truth scorer's columns (no cap on reads or loci), then source (chrom:strand, or the transcript id) and
           home (loci that qualified before the best was taken).
Stdout     reads per stratum and status; the placement table of row 1798: primary at the origin over all reads, by lead
           class (clear = one locus or lead > 10) and, with --cf, by the status the catalogue-free truth gave the read.
Known limit: loci_of lets a locus grow over a later one; records.locus_index then puts the home record on the earlier locus
and origin_bam.py calls such a decidable read no_home (primary only): 93 of 42,137 on the fly near-tied reads."""
import argparse, collections, csv, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "o2_truth")]
import neutral_truth as NT
import records as R

HEADER = "read\tstratum\tlead\tstatus\tpicks\torigin_loci\tloci\tsource\thome"
LEADS = ("clear", "1-10", "0"); ABSENT = "not in the table"

def table_sources(path):
    """-> {read: (label, chrom, strand, blocks)}"""
    out = {}
    with open(path) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        if (rows.fieldnames or [""])[0] != "read" or not {"chrom", "strand", "blocks"} <= set(rows.fieldnames):
            sys.exit(f"{path}: no header line 'read ... chrom strand blocks'")
        for r in rows:
            if r["read"] in out:
                sys.exit(f"{path}: read {r['read']!r} has two rows")
            out[r["read"]] = (f"{r['chrom']}:{r['strand']}", r["chrom"], r["strand"], [tuple(map(int, b.split("-"))) for b in r["blocks"].split(",")])
    return out

def gtf_sources(reads, gtf, name_re):
    """-> {read: (transcript id, chrom, strand, exon blocks)}, the transcript from the read name"""
    want = {}
    for q in reads:
        m = re.match(name_re, q)
        if not m:
            sys.exit(f"read name {q!r} does not match {name_re!r}")
        want[q] = m.group(1)
    need, exons, at = set(want.values()), collections.defaultdict(list), {}
    with open(gtf) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            m = re.search(r'transcript_id "([^"]+)"', f[8]) if len(f) > 8 and f[2] == "exon" else None
            if m and m.group(1) in need:
                exons[m.group(1)].append((int(f[3]) - 1, int(f[4]))); at[m.group(1)] = (f[0], f[6])
    if need - set(exons):
        sys.exit(f"{gtf}: {len(need - set(exons))} source transcripts have no exon line, e.g. {min(need - set(exons))!r}")
    return {q: (t, *at[t], sorted(exons[t])) for q, t in want.items()}

def home(loci, chrom, strand, blocks):
    """indexes of the loci on the source's contig and strand whose own blocks share a base with a source exon"""
    return [i for i, L in enumerate(loci) if L["chrom"] == chrom and L["strand"] == strand
            and any(x < v and u < y for x, y in L["blocks"] for u, v in blocks)]

def judge(recs, src, null=False):
    """-> per read (read, stratum, lead, status, primary index, origin indexes, loci, source label, home count)"""
    for q, rs in recs.items():
        if q not in src:
            sys.exit(f"read {q!r} of the BAM has no source")
        label, *where = src[q]
        loci = NT.loci_of(rs)
        lead = loci[0]["AS"] - loci[1]["AS"] if len(loci) > 1 else None
        pi = next((i for i, L in enumerate(loci) if L["primary"]), None)
        if pi is None:
            sys.exit(f"read {q!r} has no primary record")
        h = [pi] if null else home(loci, *where)      # loci are best AS first
        yield q, "S0" if lead is None else "S1" if lead > 10 else "S2", lead, "decidable" if h else "outside", pi, h[:1], loci, label, len(h)

def placement(key, n, k):
    return "\t".join((*key, str(n), str(k), f"{k / n:.4f}" if n else "-", "-".join(f"{x:.4f}" for x in R.wilson(k, n))))

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam"); ap.add_argument("out")
    g = ap.add_mutually_exclusive_group(required=True); g.add_argument("--origins"); g.add_argument("--gtf")
    ap.add_argument("--name-re", default=r"^PacBio_simulated_read_\d+_(.+)$"); ap.add_argument("--null", action="store_true"); ap.add_argument("--cf")
    a = ap.parse_args(argv)
    recs, cf = NT.records(a.bam), None
    src = gtf_sources(recs, a.gtf, a.name_re) if a.gtf else table_sources(a.origins)
    if a.cf:
        with open(a.cf) as fh:
            cfo = {r["read"]: (r["status"], r["origin_loci"].split(",")) for r in csv.DictReader(fh, delimiter="\t")}
        cf = {q: v[0] for q, v in cfo.items()}
    lines, T, P, C = [HEADER], collections.Counter(), collections.defaultdict(lambda: [0, 0]), collections.defaultdict(lambda: [0, 0, 0])
    for q, st, lead, status, pi, origin, loci, label, n in judge(recs, src, a.null):
        right = pi in origin
        T[st, "reads"] += 1; T[st, status] += 1
        if status == "decidable":
            T[st, "primary right" if right else "primary wrong"] += 1
        if cf is not None and q in cfo:      # the catalogue-free call against the exact origin (loci are the same list in both tables)
            c = C[cfo[q][0]]; c[0] += 1; c[1] += bool(origin) and str(origin[0]) in cfo[q][1]; c[2] += not origin
        lc = "clear" if lead is None or lead > 10 else "1-10" if lead else "0"
        for key in [(lc,)] + ([(lc, cf.get(q, ABSENT))] if cf is not None else []):
            P[key][0] += 1; P[key][1] += right
        lines.append("\t".join(map(str, (q, st, lead, status, f"primary={pi}", ",".join(map(str, origin)),
                     ";".join(f"{i}:{L['chrom']}:{L['s']}-{L['e']}:{L['strand']}:AS{L['AS']}" for i, L in enumerate(loci)), label, n))))
    with open(a.out, "w") as fh:      # written whole: a read that fails above leaves no half table
        fh.write("\n".join(lines) + "\n")
    cols = ("reads", "outside", "decidable", "primary right", "primary wrong")
    print("stratum\t" + "\t".join(cols))
    for st in ("S0", "S1", "S2"):
        print(st + "\t" + "\t".join(str(T[st, c]) for c in cols))
    print("\nlead\treads\tprimary at origin\tshare\twilson95")
    for lc in LEADS:
        print(placement((lc,), *P[lc,]))
    if cf is not None:
        print("\nlead\tcatalogue-free\treads\tprimary at origin\tshare\twilson95")
        for key in sorted((k for k in P if len(k) == 2), key=lambda k: (LEADS.index(k[0]), k[1] == ABSENT, k[1])):
            print(placement(key, *P[key]))
        print("\ncatalogue-free\treads\texact origin is the call / in the tie\tshare\twilson95\texact outside")
        for st in sorted(C, key=lambda s: (("decidable", "tie", "outside") + (s,)).index(s)):
            print(placement((st,), *C[st][:2]) + f"\t{C[st][2]}")

if __name__ == "__main__":
    main()
