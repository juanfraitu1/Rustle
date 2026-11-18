#!/usr/bin/env python3
"""PREREG r1799 (docs/PREREG_famwise_sim_2026-11-18.md): the simulated transcripts of a truth GTF and their reads, in the
form sim_chains.py takes (one family per gene; every read counts as carrying its transcript in full, the simulator's
truncation is not known per read).

  gtf_origins.py TRUTH.gtf EXACT_READS.tsv OUT_PREFIX [--tid-prefix P]

EXACT_READS.tsv: exact_truth.py --gtf (column source = the transcript id). Writes OUT_PREFIX.fam.json and
OUT_PREFIX.origin.tsv; --tid-prefix keeps the reads and transcripts whose id starts with P (fly: FBtr genes, FBti transposons)."""
import argparse, collections, csv, json, re

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gtf"); ap.add_argument("exact"); ap.add_argument("out"); ap.add_argument("--tid-prefix", default="")
    a = ap.parse_args(argv)
    with open(a.exact) as fh:
        reads = [(r["read"], r["source"]) for r in csv.DictReader(fh, delimiter="\t") if r["source"].startswith(a.tid_prefix)]
    want, exons, at = {t for _, t in reads}, collections.defaultdict(list), {}
    with open(a.gtf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            if t in want:
                exons[t].append([int(f[3]) - 1, int(f[4])])
                at[t] = (f[0], f[6], re.search(r'gene_id "([^"]+)"', f[8]).group(1))
    missing = want - set(at)
    assert not missing, f"{len(missing)} source transcripts are not in {a.gtf}, e.g. {sorted(missing)[:3]}"
    genes = collections.defaultdict(dict)
    for t in sorted(at):
        genes[at[t][2]][t] = sorted(exons[t])
    fam = [dict(family=g, copies=[dict(chrom=at[next(iter(tx))][0], strand=at[next(iter(tx))][1], start=min(e[0][0] for e in tx.values()),
                                       end=max(e[-1][1] for e in tx.values()), main=next(iter(tx)), transcripts=tx, gene_id=g, gene_name=g)]) for g, tx in sorted(genes.items())]
    with open(a.out + ".fam.json", "w") as fh:
        json.dump({"sets": {"sim": fam}}, fh)
    with open(a.out + ".origin.tsv", "w") as fh:
        fh.write("read\tchrom\tstrand\tblocks\tmember\tfamily\thap\tisoform\tfull\n")
        for q, t in reads:
            c, s, g = at[t]
            fh.write(f"{q}\t{c}\t{s}\t\t{g}\t{g}\tref\t{t}\t1\n")
    print(f"{len(reads)} reads of {len(at)} transcripts in {len(genes)} genes -> {a.out}.fam.json, {a.out}.origin.tsv")

if __name__ == "__main__":
    main()
