#!/usr/bin/env python3
"""A famsim scenario's LITERAL origin as a truth table in neutral_truth.py's reads.tsv format (for copyset_score.py):
each read's loci = its records' places (merged on contig + strand + overlap, best AS first), origin = the place overlapping
the copy named in the read (`copy|isoform|i`); status `decidable` (one origin) or `outside` (no record at its copy).
Under this truth a claim over identical copies is WIDE whenever some of them are silent: the source, not what the reads fit.

  sim_truth.py SCENARIO_DIR OUT.tsv"""
import collections, csv, os, sys
import pysam
d, out = sys.argv[1:3]
span = {r["copy_id"]: (r["contig"], int(r["start"]) - 1, int(r["end"]), r["strand"]) for r in csv.DictReader(open(os.path.join(d, "copies.tsv")), delimiter="\t")}
recs = collections.defaultdict(list)
with pysam.AlignmentFile(os.path.join(d, "reads.bam")) as fh:
    for r in fh.fetch(until_eof=True):
        if not r.is_unmapped and not r.is_supplementary:
            recs[r.query_name].append((r.reference_name, r.reference_start, r.reference_end, "-" if r.is_reverse else "+", r.get_tag("AS")))
with open(out, "w") as w:
    w.write("read\tstratum\tlead\tstatus\tpicks\torigin_loci\tloci\n")
    for q, rs in recs.items():
        pl = []
        for c, s, e, st, a in sorted(rs, key=lambda x: -x[4]):
            for p in pl:
                if p[0] == c and p[3] == st and p[1] < e and s < p[2]:
                    p[1], p[2] = min(p[1], s), max(p[2], e); break
            else:
                pl.append([c, s, e, st, a])
        pl = pl[:12]
        src = span.get(q.split("|")[0])
        o = [i for i, p in enumerate(pl) if src and p[0] == src[0] and p[1] < src[2] and src[1] < p[2]]
        loci = ";".join(f"{i}:{p[0]}:{p[1]}-{p[2]}:{p[3]}:AS{p[4]}" for i, p in enumerate(pl))
        w.write(f"{q}\tS2\t0\t{'decidable' if o else 'outside'}\t-\t{','.join(map(str, o[:1]))}\t{loci}\n")
