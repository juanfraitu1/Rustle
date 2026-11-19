#!/usr/bin/env python3
"""Borrowing as SHARED ISOFORMS (user, 2026-11-19: "isoforms that more than one family [member] produces"): an isoform is
shared when its spliced sequence matches an isoform assembled at another place end to end (both covered >= --cov, no
internal gap >= --max-gap bp in the cDNA-to-cDNA alignment: a gap that size is an exon one has and the other lacks).

  family_isoform.py --base P.gtf --candidates P1.gtf --fasta GENOME.fa --out PREFIX [--threads 8] [--cov 0.9] [--max-gap 20]

Writes PREFIX.shared.tsv (base isoforms with a match at another place: the descriptive count), PREFIX.promoted.gtf (candidate
multi-exon chains absent from the base whose isoform is assembled at another place by the base, labelled
fw_class "family_isoform") and PREFIX.gtf = base + promoted. No truth is read here."""
import argparse, collections, os, re, subprocess, sys
import pysam
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from undecidable import read_gtf_tx, chain, revcomp

def cdna(fa, v):
    s = "".join(fa.fetch(v["chrom"], a, b) for a, b in v["exons"])
    return revcomp(s) if v["strand"] == "-" else s

def matches(paf, cov, max_gap):
    """-> {query: set of targets} for alignments covering both >= cov with no I/D >= max_gap"""
    out = collections.defaultdict(set)
    for l in paf:
        f = l.rstrip("\n").split("\t")
        if len(f) < 12 or f[0] == f[5]:
            continue
        ql, qs, qe, tl, ts, te = int(f[1]), int(f[2]), int(f[3]), int(f[6]), int(f[7]), int(f[8])
        if (qe - qs) < cov * ql or (te - ts) < cov * tl:
            continue
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
        if any(int(n) >= max_gap for n, o in re.findall(r"(\d+)([ID])", cg)):
            continue
        out[f[0]].add(f[5])
    return out

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True); ap.add_argument("--candidates", required=True); ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--cov", type=float, default=0.9); ap.add_argument("--max-gap", type=int, default=20)
    a = ap.parse_args(argv)
    fa = pysam.FastaFile(a.fasta)
    base = {t: v for t, v in read_gtf_tx(a.base).items() if len(v["exons"]) > 1}
    bkeys = {(v["chrom"], v["strand"], chain(v["exons"])) for v in base.values()}
    cand = {t: v for t, v in read_gtf_tx(a.candidates).items() if len(v["exons"]) > 1 and (v["chrom"], v["strand"], chain(v["exons"])) not in bkeys}
    tmp = os.environ.get("TMPDIR", os.path.dirname(os.path.abspath(a.out)))
    bf, cf = f"{tmp}/famiso.{os.getpid()}.b.fa", f"{tmp}/famiso.{os.getpid()}.c.fa"
    with open(bf, "w") as fh:
        fh.writelines(f">B|{t}\n{cdna(fa, v)}\n" for t, v in base.items())
    with open(cf, "w") as fh:
        fh.writelines(f">C|{t}\n{cdna(fa, v)}\n" for t, v in cand.items())
    run = lambda q: subprocess.run(["minimap2", "-x", "asm20", "-c", "-N", "50", "-p", "0.5", "-t", str(a.threads), bf, q],
                                   check=True, capture_output=True, text=True).stdout.splitlines()
    def elsewhere(v, w):
        return not (v["chrom"] == w["chrom"] and v["exons"][0][0] < w["exons"][-1][1] and w["exons"][0][0] < v["exons"][-1][1])
    mb = matches(run(bf), a.cov, a.max_gap)
    shared = {q for q, ts in mb.items() if any(elsewhere(base[q[2:]], base[t[2:]]) for t in ts)}
    with open(a.out + ".shared.tsv", "w") as fh:
        fh.write("transcript\tshared_with\n")
        for q in sorted(shared):
            fh.write(q[2:] + "\t" + ",".join(sorted(t[2:] for t in mb[q] if elsewhere(base[q[2:]], base[t[2:]]))) + "\n")
    mc = matches(run(cf), a.cov, a.max_gap)
    prom = [q[2:] for q, ts in mc.items() if any(elsewhere(cand[q[2:]], base[t[2:]]) for t in ts)]
    lines = []
    for i, t in enumerate(sorted(prom), 1):
        v = cand[t]
        attr = f'gene_id "famiso.{i}"; transcript_id "famiso.{i}"; fw_class "family_isoform"; from "{t}";'
        lines.append(f"{v['chrom']}\tfamiso\ttranscript\t{v['exons'][0][0] + 1}\t{v['exons'][-1][1]}\t.\t{v['strand']}\t.\t{attr}\n")
        lines += [f"{v['chrom']}\tfamiso\texon\t{s + 1}\t{e}\t.\t{v['strand']}\t.\t{attr}\n" for s, e in v["exons"]]
    open(a.out + ".promoted.gtf", "w").writelines(lines)
    with open(a.out + ".gtf", "w") as fo:
        fo.write(open(a.base).read()); fo.writelines(lines)
    os.remove(bf); os.remove(cf)
    print(f"family_isoform: base multi-exon isoforms {len(base)}, shared with another place {len(shared)}; "
          f"candidates absent from the base {len(cand)}, promoted {len(prom)}")

if __name__ == "__main__":
    main()
