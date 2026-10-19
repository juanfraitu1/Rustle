#!/usr/bin/env python3
"""The sample's own copies for O2's ambiguous reads (2026-10-19; user: "lets ... focus on building on those families that are more
likely to work"; design approved the same day). O2 applies only to AS-tied multimapping reads, and on real data most of them
carry the donor's own alleles, which CHM13's copies lack equally (registers 1182, 1329): CHM13's PSVs cannot place them. So the
family's reads are clustered by the alleles they share (read-backed phasing over N copies: SDA / Chaisson 2017, IsoCon), each
cluster is one of the sample's copy-haplotypes, and a cluster takes the copy minimap2 resolves its untied reads to (98% right
against the donor's genome, register 1330). A tied read takes its cluster's copy, or abstains.

    samplecopy.py run FAMILIES.json SET FAMILY CHM13.bam TIED.tsv UNTIED.tsv GENOME.fa OUT_DIR
                                         -> OUT_DIR/<family>.calls.tsv (read, call, posterior, cluster, label) + .clusters.json

Steps (all reads of the family, tied or not):
1. every read aligned to one copy's locus (the anchor: the copy with most untied reads), minimap2 splice:hq;
2. informative columns: substitutions where two alleles each have >= max(min_reads, min_frac x coverage) reads;
3. a mixture of haplotypes over those columns (EM, the error rate estimated), started at kmax = 2 x copies and reduced by
   merging clusters that differ at fewer than two columns (one difference can be one error) and dissolving clusters under
   min_size reads;
4. labels: a cluster's copy is the one its untied members resolve to by a clear majority;
5. a tied read takes its cluster's label when its posterior is >= 0.99 and the label is a copy it ties on.
"""
import collections
import csv
import json
import math
import os
import random
import subprocess
import sys
import tempfile

BASES = "ACGT"


def informative_columns(reads, min_reads=3, min_frac=0.1):
    """Columns where two or more alleles each have >= max(min_reads, min_frac x the column's coverage) reads; reads are
    {column: base}."""
    cnt = collections.defaultdict(collections.Counter)
    for r in reads:
        for c, b in r.items():
            cnt[c][b] += 1
    out = []
    for c, v in cnt.items():
        floor = max(min_reads, min_frac * sum(v.values()))
        if sum(1 for n in v.values() if n >= floor) >= 2:
            out.append(c)
    return sorted(out)


def _loglik(read, hap, le, lm):
    return sum(le if hap.get(c) == b else lm for c, b in read.items())


def _em(reads, haps, iters=50, fixed_e=None):
    """EM of a haplotype mixture from initial haplotypes [{col: base}] -> (haps, posteriors, weights, error, loglik); fixed_e:
    the error rate held at the reads' background (background_error) rather than fitted."""
    k = len(haps)
    w = [1.0 / k] * k
    e = fixed_e if fixed_e is not None else 0.01
    post = []
    ll = -math.inf
    for _ in range(iters):
        le, lm = math.log(1 - e), math.log(e / 3)
        post, ll_new = [], 0.0
        for r in reads:
            lp = [math.log(max(w[j], 1e-300)) + _loglik(r, haps[j], le, lm) for j in range(k)]
            m = max(lp)
            z = sum(math.exp(x - m) for x in lp)
            post.append([math.exp(x - m) / z for x in lp])
            ll_new += m + math.log(z)
        w = [sum(p[j] for p in post) / len(reads) for j in range(k)]
        cnt = [collections.defaultdict(collections.Counter) for _ in range(k)]
        mism, obs = 0.0, 0.0
        for r, p in zip(reads, post):
            for j in range(k):
                if p[j] > 1e-9:
                    for c, b in r.items():
                        cnt[j][c][b] += p[j]
        haps = [{c: max(sorted(v), key=v.__getitem__) for c, v in cnt[j].items()} or haps[j] for j in range(k)]
        for r, p in zip(reads, post):
            for j in range(k):
                mism += p[j] * sum(1 for c, b in r.items() if haps[j].get(c) != b)
            obs += len(r)
        if fixed_e is None:
            e = min(max(mism / obs if obs else 0.01, 1e-4), 0.25)
        if abs(ll_new - ll) < 1e-6:
            ll = ll_new
            break
        ll = ll_new
    return haps, post, w, e, ll


def _seeds(reads, k, seed):
    """k starting haplotypes: a seeded first read, then each next the read farthest from those chosen (Hamming over shared
    columns), its missing columns filled with the column's majority."""
    maj = {}
    cnt = collections.defaultdict(collections.Counter)
    for r in reads:
        for c, b in r.items():
            cnt[c][b] += 1
    for c, v in cnt.items():
        maj[c] = max(sorted(v), key=v.__getitem__)
    rng = random.Random(seed)
    chosen = [rng.randrange(len(reads))]

    def dist(a, b):
        return sum(1 for c in a if c in b and a[c] != b[c])
    while len(chosen) < min(k, len(reads)):
        best = max(range(len(reads)), key=lambda i: (min(dist(reads[i], reads[j]) for j in chosen), -i))
        if min(dist(reads[best], reads[j]) for j in chosen) == 0:
            break
        chosen.append(best)
    return [{**maj, **reads[i]} for i in chosen]


def _support(post, reads, j):
    """{column: posterior weight of cluster j's reads covering it}."""
    s = collections.Counter()
    for r, p in zip(reads, post):
        for c in r:
            s[c] += p[j]
    return s


def fit_clusters(reads, cols, kmax, seed=0, min_size=3, min_diff=2, iters=50, error=None):
    """The sample's copy-haplotypes over `cols`: EM from kmax farthest-point seeds, then clusters closer than min_diff
    columns (compared where both have >= 1 read of support) merged and clusters under min_size reads dissolved, EM rerun
    until nothing changes. -> {haplotypes, posteriors, weights, error, loglik}."""
    reads = [{c: b for c, b in r.items() if c in set(cols)} for r in reads]
    haps = _seeds(reads, kmax, seed)
    while True:
        haps, post, w, e, ll = _em(reads, haps, iters, error)
        k = len(haps)
        size = [sum(p[j] for p in post) for j in range(k)]
        small = [j for j in range(k) if size[j] < min_size]
        if small and k > 1:
            drop = min(small, key=lambda j: size[j])
            haps = [h for j, h in enumerate(haps) if j != drop]
            continue
        sup = [_support(post, reads, j) for j in range(k)]
        pair = None
        for a in range(k):
            for b in range(a + 1, k):
                shared = [c for c in cols if sup[a][c] >= 1 and sup[b][c] >= 1]
                if sum(1 for c in shared if haps[a].get(c) != haps[b].get(c)) < min_diff:
                    pair = (a, b) if pair is None or size[a] + size[b] < sum(size[x] for x in pair) else pair
        if pair:
            a, b = pair
            keep = a if size[a] >= size[b] else b
            haps = [h for j, h in enumerate(haps) if j not in pair] + [haps[keep]]
            continue
        return dict(haplotypes=haps, posteriors=post, weights=w, error=e, loglik=ll)


def _wilson_lower(k, n, z=1.96):
    if n == 0:
        return 0.0
    p, d = k / n, 1 + z * z / n
    return ((p + z * z / (2 * n)) - z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / d


def label_clusters(member, resolved, n_clusters):
    """{cluster: the copy its members' minimap2 resolutions (untied reads; None for tied ones) name by a CONFIDENT majority --
    the 95% Wilson lower bound of its share above one half -- or None}; member: each read's cluster (or None). (Development,
    TBC1D3: a label from 2 untied reads of 2 was wrong for every tied read it placed.)"""
    cnt = [collections.Counter() for _ in range(n_clusters)]
    for m, r in zip(member, resolved):
        if m is not None and r is not None:
            cnt[m][r] += 1
    out = {}
    for j, v in enumerate(cnt):
        top = max(sorted(v), key=v.__getitem__) if v else None
        out[j] = top if top is not None and _wilson_lower(v[top], sum(v.values())) > 0.5 else None
    return out


MIN_SITES = 2                                                       # O2's certificate: one site can be one error


def decisive_sites(read, j, haps, labels):
    """How many of the read's own sites back cluster j against its closest differently-labelled rival (an unlabelled cluster
    counts as different): columns where the read carries j's allele and the rival's differs. 0 with no such rival."""
    rivals = [k for k in range(len(haps)) if k != j and labels.get(k) != labels.get(j)]
    if not rivals:
        return 0
    return min(sum(1 for c, b in read.items() if haps[j].get(c) == b and haps[k].get(c) != b) for k in rivals)


def assign_tied(post, labels, tied, min_post=0.99, sites=None):
    """A tied read's call: its most likely cluster's label when the posterior is >= min_post, the label is one of the copies
    the read ties on and (when given) the read's decisive sites number >= MIN_SITES; else None (abstain)."""
    j = max(range(len(post)), key=post.__getitem__)
    lab = labels.get(j)
    if sites is not None and sites < MIN_SITES:
        return None
    return lab if post[j] >= min_post and lab is not None and lab in tied else None


def read_pileup(rec):
    """{reference position: read base} over a pysam alignment's aligned (=, X, M) positions; insertions, deletions and clips
    are left out (substitutions only: HiFi's errors are mostly homopolymer indels)."""
    q = rec.query_sequence
    return {rp: q[qp] for qp, rp in rec.get_aligned_pairs(matches_only=True)}


def choose_anchor(resolved_counts, copies):
    """The anchor copy: the one most untied reads resolve to (the best-supported frame); the first copy when none do."""
    return max(copies, key=lambda c: (resolved_counts.get(c, 0), -copies.index(c)))


def background_error(reads, ref, informative):
    """The reads' mismatch rate to the anchor's bases ({position: base}) off the informative columns: the sequencing error
    the mixture's likelihood should assume (fitting it instead lets an under-split model call its own mixing error)."""
    inf = set(informative)
    mism = obs = 0
    for r in reads:
        for c, b in r.items():
            if c in inf or c not in ref:
                continue
            obs += 1
            mism += b != ref[c]
    return mism / obs if obs else 0.0


def evidence_posteriors(reads, haps, e):
    """Each read's posterior over the clusters from its own alleles alone (a uniform prior): the clusters' sizes must not
    place a read that carries no informative site (one large cluster would otherwise take every such read)."""
    le, lm = math.log(1 - e), math.log(e / 3)
    out = []
    for r in reads:
        lp = [_loglik(r, h, le, lm) for h in haps]
        m = max(lp)
        z = sum(math.exp(x - m) for x in lp)
        out.append([math.exp(x - m) / z for x in lp])
    return out


PAD, KMAX_CAP, MIN_POST = 1000, 24, 0.99
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def cmd_run(fam_p, set_name, family, bam_p, tied_p, untied_p, genome_p, out_dir, threads=2):
    import pysam
    fam = {f["family"]: f for f in json.load(open(fam_p))["sets"][set_name]}[family]
    copies = [f"{family}.{k}" for k in range(len(fam["copies"]))]
    seqs = {}
    with pysam.AlignmentFile(bam_p) as bam:
        for c in fam["copies"]:
            for rec in bam.fetch(c["chrom"], max(0, c["start"] - PAD), c["end"] + PAD):
                if rec.is_secondary or rec.is_supplementary or rec.is_unmapped or rec.query_sequence is None:
                    continue
                if rec.query_name not in seqs:
                    q = rec.query_sequence
                    seqs[rec.query_name] = q[::-1].translate(COMP) if rec.is_reverse else q
    tied = {r["read"]: set(r["copies"].split(",")) for r in csv.DictReader(open(tied_p), delimiter="\t")
            if r["family"] == family}
    resolved = {r["read"]: r["copies"] for r in csv.DictReader(open(untied_p), delimiter="\t") if r["family"] == family}
    anchor = choose_anchor(collections.Counter(resolved[q] for q in seqs if q in resolved), copies)
    a = fam["copies"][copies.index(anchor)]
    with tempfile.TemporaryDirectory(dir=out_dir) as d:
        s0 = max(0, a["start"] - PAD)
        with pysam.FastaFile(genome_p) as fa:
            locus = fa.fetch(a["chrom"], s0, a["end"] + PAD)
        open(f"{d}/anchor.fa", "w").write(f">anchor\n{locus}\n")
        open(f"{d}/reads.fa", "w").write("".join(f">{q}\n{x}\n" for q, x in seqs.items()))
        subprocess.run(f"minimap2 -a -x splice:hq -uf --eqx --secondary=no -t {threads} {d}/anchor.fa {d}/reads.fa > {d}/aln.sam "
                       f"2> {d}/mm2.log", shell=True, check=True)
        pile = {}
        with pysam.AlignmentFile(f"{d}/aln.sam") as al:
            for rec in al:
                if not (rec.is_unmapped or rec.is_secondary or rec.is_supplementary):
                    pile[rec.query_name] = read_pileup(rec)
    names = sorted(pile)
    reads = [pile[q] for q in names]
    kmax = min(2 * len(copies), KMAX_CAP)
    cols = informative_columns(reads, min_reads=3, min_frac=0.5 / kmax)    # half of one copy-haplotype's expected share
    ref = dict(enumerate(locus.upper()))                                   # the anchor FASTA is the locus: positions from 0
    err = max(background_error(reads, ref, cols), 1e-4)
    fit = fit_clusters(reads, cols, kmax=kmax, error=err)
    sub = [{c: b for c, b in r.items() if c in set(cols)} for r in reads]
    ev = evidence_posteriors(sub, fit["haplotypes"], fit["error"])
    member = [max(range(len(p)), key=p.__getitem__) if max(p) >= MIN_POST else None for p in ev]
    labels = label_clusters(member, [resolved.get(q) for q in names], len(fit["haplotypes"]))
    with open(f"{out_dir}/{family}.calls.tsv", "w") as fh:
        fh.write("read\tcall\tposterior\tcluster\tlabel\n")
        for q, p in zip(names, ev):
            if q in tied:
                j = max(range(len(p)), key=p.__getitem__)
                ns = decisive_sites(sub[names.index(q)], j, fit["haplotypes"], labels)
                fh.write(f"{q}\t{assign_tied(p, labels, tied[q], MIN_POST, ns) or '-'}\t{p[j]:.4f}\t{j}\t{labels.get(j) or '-'}\n")
    lab_counts = [collections.Counter(resolved[q] for q, m in zip(names, member) if m == j and q in resolved)
                  for j in range(len(fit["haplotypes"]))]
    json.dump(dict(family=family, anchor=anchor, reads=len(names), tied=sum(q in tied for q in names), columns=len(cols),
                   error=fit["error"], clusters=[dict(size=sum(m == j for m in member), label=labels[j], untied=dict(lab_counts[j]),
                                                      haplotype="".join(fit["haplotypes"][j].get(c, "-") for c in cols))
                                                 for j in range(len(fit["haplotypes"]))]),
              open(f"{out_dir}/{family}.clusters.json", "w"), indent=1)
    n = collections.Counter("placed" if r.split("\t")[1] != "-" else "abstain" for r in open(f"{out_dir}/{family}.calls.tsv").readlines()[1:])
    print(f"[samplecopy] {family}: {len(names)} reads, {len(cols)} informative columns, {len(fit['haplotypes'])} clusters "
          f"({sum(1 for v in labels.values() if v)} labelled), error {fit['error']:.4f}; tied reads {dict(n)}", flush=True)
    return 0


def cmd_score(truth_p, out_dir, families):
    """Per family: the tied reads with a donor copy (unique / unique_elsewhere) placed right / wrong / abstained, and the
    donor-tied (own_tie) ones placed."""
    truth = {r["read"]: r for r in csv.DictReader(open(truth_p), delimiter="\t")}
    tot = collections.Counter()
    for f in families:
        n = collections.Counter()
        for r in csv.DictReader(open(f"{out_dir}/{f}.calls.tsv"), delimiter="\t"):
            t = truth.get(r["read"])
            if t is None:
                continue
            placed = r["call"] != "-"
            if t["class"] in ("unique", "unique_elsewhere"):
                n["right" if placed and r["call"] == t["truth"] else "wrong" if placed else "abstain"] += 1
            elif t["class"] == "own_tie":
                n["own_tie placed" if placed else "own_tie abstain"] += 1
            else:
                n[f"{t['class']} {'placed' if placed else 'abstain'}"] += 1
        tot.update(n)
        print(f"{f}: {dict(sorted(n.items()))}", flush=True)
    print(f"TOTAL: {dict(sorted(tot.items()))}", flush=True)
    return 0


def main(argv):
    if len(argv) == 9 and argv[0] == "run":
        return cmd_run(*argv[1:])
    if len(argv) >= 4 and argv[0] == "score":
        return cmd_score(argv[1], argv[2], argv[3:])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
