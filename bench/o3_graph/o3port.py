#!/usr/bin/env python3
"""O3 port (2026-10-23): Strainy's phasing with IsoCon's significance test, and an allele-slot certificate.
Design: docs/superpowers/specs/2026-10-23-o3-strainy-isocon-design.md; decision rules: docs/PREREG_o3_port_2026-10-23.md.

  o3port.py run --bam B --copies COPIES.tsv --genome G.fa --net net_reads.tsv --out PREFIX [--families F1,F2]

Per family (the reads of its net with used = 1, as `candidate_copies` built it):
  1. piles: a read joins the pile of every copy holding one of its alignments within AS_TIE of its best (O2's tie rule);
  2. informative columns (Strainy's T-informative SNPs, IsoCon's test): >= 2 alleles each carried by >= MIN_READS reads at
     base quality >= MIN_Q, a minor allele improbable as sequencing error (Poisson tail of the carriers' error probabilities,
     Bonferroni over the tested columns, ALPHA); a base that is no supported allele counts as missing;
  3. read graph (clusteringDivergence 0): two reads are linked when they share >= 1 informative column and disagree at none;
  4. clusters: deterministic label propagation (each distinct allele pattern votes with its read count);
  5. refinement: each cluster is re-clustered on the columns informative inside it, until nothing splits;
  6. certificate: per pile, the largest set of pairwise-conflicting clusters (distinct haplotypes present together); over the
     family, a maximum b-matching of those haplotypes into allele slots (2 per copy) of the copies whose pile holds >= half of
     their reads; an unmatched haplotype is a certified extra copy;
  7. its contig: the consensus of its reads (majority base, deletion or insertion at every position >= MIN_READS of them cover)
     along their alignments on the copy whose pile holds most of them.
Products (the `candidate_copies` contract the acceptance harness reads): PREFIX.candidates.tsv (one row per certified haplotype,
flagged = 1), PREFIX.contigs.fa, PREFIX.clusters.tsv, PREFIX.reads.tsv, PREFIX.families.tsv.
"""
import argparse
import collections
import csv
import math
import sys

MIN_READS = 2      # the project's ">= 2 independent observations"
MIN_Q = 10         # base quality below which a base counts as missing
AS_TIE = 10        # O2's tie rule (Eichler): alignments within 10 AS of the read's best
ALPHA = 0.01       # IsoCon-style test, Bonferroni over the tested columns
MAX_ROUNDS = 100   # label propagation stops earlier when no label changes
MIN_SUPPORT = 6    # reads a haplotype needs to take a slot or be certified: the shipped stage's own flag floor
                   # (candidate_copies --min-support 6, ruling R1), so the port is judged on the stage's terms

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ---------------------------------------------------------------- statistics and columns (pure; tested)

def poisson_tail(k, lam):
    """P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    term = math.exp(-lam)
    cdf = term
    for j in range(1, k):
        term *= lam / j
        cdf += term
    return max(0.0, 1.0 - cdf) if cdf < 1.0 - 1e-12 else _small_tail(k, lam)


def _small_tail(k, lam):
    """The tail by its leading terms when 1 - cdf underflows (lam small); log space, so a large k cannot overflow."""
    lt = -lam + k * math.log(lam) - math.lgamma(k + 1)
    s, j = 0.0, k
    while lt > -700 and j < k + 50:
        s += math.exp(lt)
        j += 1
        lt += math.log(lam) - math.log(j)
    return s


def err(x):
    """A read's per-base error probability at a column: its empirical rate on the copy when known (IsoCon estimates error
    rates from the reads' own alignments; PacBio base qualities are capped at Q93 and overconfident), else its base quality."""
    return x[3] if len(x) > 3 else 10 ** (-x[2] / 10)


def informative(obs, n_tests=None):
    """obs: {column: [(read, base, qual[, error rate]), ...]}. Returns {column: set of supported alleles} for the columns where
    >= 2 alleles have >= MIN_READS carriers at quality >= MIN_Q and a minor allele passes the error test (Bonferroni over the
    columns tested here, i.e. those with >= 2 such alleles, unless n_tests is given)."""
    cand = {}
    for col, xs in obs.items():
        good = [x for x in xs if x[2] >= MIN_Q and x[1] in "ACGT"]
        cnt = collections.Counter(x[1] for x in good)
        alleles = [b for b, n in cnt.items() if n >= MIN_READS]
        if len(alleles) >= 2:
            cand[col] = (good, cnt, alleles)
    m = n_tests if n_tests is not None else max(1, len(cand))
    out = {}
    for col, (good, cnt, alleles) in cand.items():
        lam = sum(err(x) / 3 for x in good)
        major = max(alleles, key=lambda b: (cnt[b], b))
        real = {major} | {b for b in alleles if b != major and poisson_tail(cnt[b], lam) < ALPHA / m}
        if len(real) >= 2:
            out[col] = real
    return out


def patterns(reads, obs, cols):
    """{read: {col: allele}} restricted to informative columns, a base that is no supported allele counted as missing."""
    pat = {r: {} for r in reads}
    for col, real in cols.items():
        for x in obs.get(col, ()):
            r, b, q = x[0], x[1], x[2]
            if r in pat and q >= MIN_Q and b in real:
                pat[r][col] = b
    return pat


def compatible(a, b):
    shared = a.keys() & b.keys()
    return bool(shared) and all(a[c] == b[c] for c in shared)


def label_propagation(pat):
    """Deterministic label propagation on the zero-disagreement read graph. Reads with the same allele pattern have the same
    neighbourhood, so the graph is built on the distinct patterns, each voting with its read count (a read also hears the other
    reads of its own pattern). Returns a list of read sets (reads with no informative column are left out)."""
    groups = collections.defaultdict(list)
    for r, p in pat.items():
        if p:
            groups[tuple(sorted(p.items()))].append(r)
    keys = sorted(groups, key=lambda k: (-len(groups[k]), k))
    idx = {k: i for i, k in enumerate(keys)}
    nbr = {i: [] for i in range(len(keys))}
    dicts = [dict(k) for k in keys]
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if compatible(dicts[i], dicts[j]):
                nbr[i].append(j)
                nbr[j].append(i)
    label = list(range(len(keys)))
    w = [len(groups[k]) for k in keys]
    for _ in range(MAX_ROUNDS):
        changed = False
        for i in range(len(keys)):
            votes = collections.Counter({label[i]: w[i] - 1})
            for j in nbr[i]:
                votes[label[j]] += w[j]
            top = max(votes.values())
            if votes[label[i]] == top:
                continue
            new = min(l for l, v in votes.items() if v == top)
            label[i] = new
            changed = True
        if not changed:
            break
    clusters = collections.defaultdict(set)
    for k, i in idx.items():
        clusters[label[i]].update(groups[k])
    return sorted(clusters.values(), key=lambda s: (-len(s), sorted(s)))


def cluster_pile(reads, obs):
    """Strainy's clustering with refinement: label propagation on the pile's informative columns, then each cluster again on
    the columns informative inside it, until nothing splits. Returns (clusters of >= MIN_READS reads, every column used, the
    pile-level columns). Only the pile-level columns may prove two haplotypes distinct (the certificate): a column that passes
    the test only inside one cluster (smaller Bonferroni family, lower error mass) can come from a recurrent error or chimera."""
    cols = informative({c: [x for x in xs if x[0] in reads] for c, xs in obs.items()})
    used = dict(cols)
    todo = [c for c in label_propagation(patterns(reads, obs, cols)) if len(c) >= MIN_READS]
    done = []
    while todo:
        c = todo.pop()
        inner = informative({col: [x for x in xs if x[0] in c] for col, xs in obs.items()})
        if not inner:
            done.append(c)
            continue
        parts = [p for p in label_propagation(patterns(c, obs, inner)) if len(p) >= MIN_READS]
        if len(parts) >= 2:
            for col, real in inner.items():
                used.setdefault(col, set()).update(real)
            todo.extend(parts)
        else:
            done.append(c)
    return sorted(done, key=lambda s: (-len(s), sorted(s))), used, cols


def haplotype(cluster, obs, cols):
    """A cluster's allele at each column: the strict majority of its members' supported alleles, carried by >= MIN_READS."""
    h = {}
    for col, real in cols.items():
        cnt = collections.Counter(x[1] for x in obs.get(col, ()) if x[0] in cluster and x[2] >= MIN_Q and x[1] in real)
        if cnt:
            b, n = cnt.most_common(1)[0]
            if n >= MIN_READS and 2 * n > sum(cnt.values()):
                h[col] = b
    return h


def conflicting_set(haps):
    """The largest set of pairwise-conflicting haplotypes: distinct haplotypes present together. Two haplotypes conflict when
    they differ at >= MIN_READS shared columns (two independent observations of a difference, the same >= 2 rule as for reads;
    a single differing column is one observation, e.g. a recurrent error splitting off a few reads). Exact by branch and
    bound; piles hold few clusters."""
    n = len(haps)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            shared = haps[i].keys() & haps[j].keys()
            adj[i][j] = adj[j][i] = sum(haps[i][c] != haps[j][c] for c in shared) >= MIN_READS
    best = []

    def grow(cur, rest):
        nonlocal best
        if len(cur) > len(best):
            best = list(cur)
        for k, v in enumerate(rest):
            if len(cur) + len(rest) - k <= len(best):
                return
            grow(cur + [v], [u for u in rest[k + 1:] if adj[v][u]])

    grow([], list(range(n)))
    return sorted(best)


def b_matching(options, capacity, order=None):
    """options: {haplotype: [copies]}; capacity: {copy: slots}. Maximum b-matching by augmenting paths, haplotypes taken in
    `order` (default: fewest options first). An augmenting path never unmatches a matched haplotype, so the ones left
    unmatched are the last in the order. Returns {haplotype: copy or None}."""
    load = collections.defaultdict(list)
    match = {}

    def augment(h, seen):
        for c in options[h]:
            if c in seen:
                continue
            seen.add(c)
            if len(load[c]) < capacity[c]:
                load[c].append(h)
                match[h] = c
                return True
            for other in list(load[c]):
                if augment_from(other, c, seen):
                    load[c].remove(other)
                    load[c].append(h)
                    match[h] = c
                    return True
        return False

    def augment_from(h, avoid, seen):
        for c in options[h]:
            if c == avoid or c in seen:
                continue
            seen.add(c)
            if len(load[c]) < capacity[c]:
                load[c].append(h)
                match[h] = c
                return True
            for other in list(load[c]):
                if augment_from(other, c, seen):
                    load[c].remove(other)
                    load[c].append(h)
                    match[h] = c
                    return True
        return False

    for h in (order if order is not None else sorted(options, key=lambda h: (len(options[h]), h))):
        if not augment(h, set()):
            match[h] = None
    return match


def consensus(alignments):
    """alignments: {read: [(ref_pos, base or '-', inserted string)]} on one copy. The majority at every position >= MIN_READS
    of them cover (a base, or '-' for a deletion), followed by the majority insertion when most of the covering reads carry
    one. Positions in genome order."""
    at = collections.defaultdict(list)
    for read, steps in alignments.items():
        for pos, base, ins in steps:
            at[pos].append((base, ins))
    out = []
    for pos in sorted(at):
        xs = at[pos]
        if len(xs) < MIN_READS:
            continue
        b, _ = collections.Counter(b for b, _ in xs).most_common(1)[0]
        if b != "-":
            out.append(b)
        ins = collections.Counter(i for _, i in xs if i)
        if ins:
            s, n = ins.most_common(1)[0]
            if 2 * sum(ins.values()) > len(xs) and n >= MIN_READS:
                out.append(s)
    return "".join(out)


# ---------------------------------------------------------------- BAM side

def load_copies(path):
    by_fam = collections.defaultdict(list)
    for r in csv.DictReader(open(path), delimiter="\t"):
        blocks = sorted(tuple(map(int, b.split("-"))) for b in r["exons"].split(",") if b)
        by_fam[r["family_id"]].append((r["tid"], r["chrom"], blocks[0][0], blocks[-1][1]))
    return by_fam


def load_net(path, families):
    net = collections.defaultdict(set)
    n_net = collections.Counter()
    for r in csv.DictReader(open(path), delimiter="\t"):
        if families and r["family"] not in families:
            continue
        n_net[r["family"]] += 1
        if r["used"] == "1":
            net[r["family"]].add(r["read"])
    return net, n_net


def read_primaries(bam, wanted):
    """{read: (seq, quals, is_reverse, best AS over every record of the read)} in one pass over the BAM."""
    prim, best = {}, collections.defaultdict(lambda: -10 ** 9)
    for rec in bam.fetch(until_eof=True):
        n = rec.query_name
        if n not in wanted or rec.is_unmapped:
            continue
        if rec.has_tag("AS"):
            best[n] = max(best[n], rec.get_tag("AS"))
        if not rec.is_secondary and not rec.is_supplementary and rec.query_sequence:
            q = rec.query_qualities
            prim[n] = (rec.query_sequence, list(q) if q is not None else [40] * len(rec.query_sequence), rec.is_reverse)
    return {n: v + (best[n],) for n, v in prim.items()}


def oriented(prim, rec):
    """The read's bases and qualities in this record's orientation, and the offset of its first CIGAR query base."""
    seq, qual, rev = prim[:3]
    if rec.is_reverse != rev:
        seq, qual = seq.translate(COMP)[::-1], qual[::-1]
    lead = rec.cigartuples[0][1] if rec.cigartuples and rec.cigartuples[0][0] == 5 else 0
    return seq, qual, lead


def pile_records(bam, copies, reads, prim):
    """{copy: {read: record}} for the records within AS_TIE of the read's best that fall on the copy's span."""
    piles = {}
    for cid, chrom, s, e in copies:
        p = {}
        for rec in bam.fetch(chrom, s, e):
            n = rec.query_name
            if n not in reads or n not in prim or rec.is_unmapped or not rec.has_tag("AS"):
                continue
            if rec.get_tag("AS") < prim[n][3] - AS_TIE:
                continue
            if n not in p or rec.get_tag("AS") > p[n].get_tag("AS"):
                p[n] = rec
        piles[cid] = p
    return piles


def walk(rec):
    """(kind, ref position, record query index) along the CIGAR: M/=/X bases, D positions, I runs (attached to the previous
    reference position, with the run's length); N (introns), S and H advance without emitting."""
    qp, rp = 0, rec.reference_start
    for op, n in rec.cigartuples:
        if op in (0, 7, 8):
            for k in range(n):
                yield "M", rp + k, qp + k
            qp += n
            rp += n
        elif op == 1:
            yield "I", rp - 1, (qp, n)
            qp += n
        elif op == 2:
            for k in range(n):
                yield "D", rp + k, None
            rp += n
        elif op == 3:
            rp += n
        elif op == 4:
            qp += n


def observations(pile, prim, span):
    """{column: [(read, base, qual)]} and {read: [(pos, base or '-', insertion)]} for a pile, restricted to the copy's span."""
    obs = collections.defaultdict(list)
    steps = {}
    s0, e0 = span
    for n, rec in pile.items():
        seq, qual, lead = oriented(prim[n], rec)
        aligned = sum(k for op, k in rec.cigartuples if op in (0, 1, 7, 8))
        eps = rec.get_tag("NM") / aligned if rec.has_tag("NM") and aligned else None
        st = []
        for kind, rp, q in walk(rec):
            if not (s0 <= rp < e0):
                continue
            if kind == "M":
                b, bq = seq[q + lead], qual[q + lead]
                st.append((rp, b, ""))
                obs[rp].append((n, b, bq) if eps is None else (n, b, bq, eps))
            elif kind == "D":
                st.append((rp, "-", ""))
            elif st and st[-1][0] == rp:
                qs, ln = q
                pos, b, ins = st[-1]
                st[-1] = (pos, b, ins + seq[qs + lead:qs + lead + ln])
        steps[n] = st
    return obs, steps


def explanation(hp, parents):
    """DADA2's bimera rule on a haplotype's columns. parents: [(name, {col: base})], each more abundant than the haplotype or a
    reference copy. Returns ("single", name) when one parent carries every allele, ("join", (left, right)) when a left part
    matches one parent and the rest another (a chimeric molecule), else None. A column a parent does not cover counts as a
    mismatch for that parent."""
    cols = sorted(hp)
    if not cols:
        return None
    pre, suf = [], []
    for name, pb in parents:
        m = [pb.get(c) == hp[c] for c in cols]
        if all(m):
            return "single", name
        lp = next((k for k, x in enumerate(m) if not x), len(m))
        ls = next((k for k, x in enumerate(reversed(m)) if not x), len(m))
        pre.append((lp, name))
        suf.append((ls, name))
    pre.sort(reverse=True)
    suf.sort(reverse=True)
    for lp, a in pre[:2]:
        for ls, b in suf[:2]:
            if a != b and lp + ls >= len(cols):
                return "join", (a, b)
    return None


def homology(home, others, genome, tmpdir):
    """{other copy: {home position: (other position, strand)}} from minimap2 -c --eqx -x asm20 of each other copy's locus
    against the home copy's locus (best-matching alignment first; a position keeps its first mapping)."""
    import os
    import re
    import subprocess
    hid, hch, hs, he = home
    t = os.path.join(tmpdir, "t.fa")
    q = os.path.join(tmpdir, "q.fa")
    with open(t, "w") as fh:
        fh.write(f">{hid}\n{genome.fetch(hch, hs, he)}\n")
    with open(q, "w") as fh:
        for cid, ch, s, e in others:
            fh.write(f">{cid}\n{genome.fetch(ch, s, e)}\n")
    paf = subprocess.run(["minimap2", "-c", "--eqx", "-x", "asm20", "-t", "1", t, q], capture_output=True, text=True).stdout
    start = {cid: s for cid, _, s, _ in others}
    maps = {}
    lines = [ln.split("\t") for ln in paf.splitlines()]
    for p in sorted(lines, key=lambda p: -int(p[9])):
        cid, strand, qs, qe, ts = p[0], p[4], int(p[2]), int(p[3]), int(p[7])
        cg = next(x[5:] for x in p[12:] if x.startswith("cg:Z:"))
        mp = maps.setdefault(cid, {})
        tpos, qo = ts, 0
        for n, op in re.findall(r"(\d+)([=XIDM])", cg):
            n = int(n)
            if op in "=XM":
                for k in range(n):
                    qq = qs + qo + k if strand == "+" else qe - 1 - (qo + k)
                    mp.setdefault(hs + tpos + k, (start[cid] + qq, strand))
                tpos += n
                qo += n
            elif op == "I":
                qo += n
            elif op == "D":
                tpos += n
    return maps


def pile_consensus(obs):
    return {c: collections.Counter(x[1] for x in xs).most_common(1)[0][0] for c, xs in obs.items() if xs}


CMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def private_alleles(i, merged, match, options):
    """Alleles of unmatched haplotype i that no slot-holding haplotype of its option copies carries, counted per pile where
    i has a pattern (the largest count). A haplotype whose every allele some slot holder carries is a mosaic of them (a
    chimeric molecule), not sequence the reference copies cannot hold: an extra copy needs >= MIN_READS private alleles."""
    best = 0
    for cid, hp in merged[i]["hap"].items():
        holders = [merged[j]["hap"][cid] for j in match if match[j] is not None and j != i
                   and cid in merged[j]["hap"] and cid in options[i]]
        best = max(best, sum(1 for col, b in hp.items() if all(h.get(col) != b for h in holders)))
    return best


def divergence(h, ref):
    """Fraction of a haplotype's columns where its allele is not the reference base."""
    cols = [c for c in h if ref.get(c, "N") != "N"]
    return sum(h[c] != ref[c] for c in cols) / len(cols) if cols else 0.0


def run_family(fam, copies, reads, bam, prim, genome, slot_order="divergence", tmpdir="."):
    piles = pile_records(bam, copies, reads, prim)
    span = {cid: (s, e) for cid, _, s, e in copies}
    chrom = {cid: ch for cid, ch, _, _ in copies}
    haps, per_copy, allhaps, cons, refs = [], {}, [], {}, {}
    for cid, pile in piles.items():
        if len(pile) < 2 * MIN_READS:
            continue
        obs, steps = observations(pile, prim, span[cid])
        clusters, cols, top = cluster_pile(set(pile), obs)
        hs = [haplotype(c, obs, top) for c in clusters]
        keep = conflicting_set(hs) if hs else []
        per_copy[cid] = (len(pile), len(clusters), len(keep))
        ref = {col: genome.fetch(chrom[cid], col, col + 1).upper() for col in top}
        cons[cid] = pile_consensus(obs)
        for i, c in enumerate(clusters):
            allhaps.append((cid, len(c), hs[i]))
        for i in keep:
            haps.append((cid, clusters[i], steps, divergence(hs[i], ref), hs[i]))
    # one haplotype found in several piles (tied reads) is one haplotype: merge by read overlap (>= half of the smaller)
    merged = []
    for cid, reads_h, steps, dv, hp in sorted(haps, key=lambda x: (-len(x[1]), x[0])):
        for m in merged:
            if len(m["reads"] & reads_h) * 2 >= min(len(m["reads"]), len(reads_h)):
                m["piles"].add(cid)
                m["div"] = min(m["div"], dv)
                m["hap"].setdefault(cid, hp)
                break
        else:
            merged.append({"reads": set(reads_h), "piles": {cid}, "div": dv, "hap": {cid: hp}})
    options = {}
    for i, m in enumerate(merged):
        options[i] = sorted(c for c, p in piles.items() if len(m["reads"] & p.keys()) * 2 >= len(m["reads"]))
    # DADA2's rule before the slots: a duplicate of a more abundant haplotype shares its slot, a chimera (left part one parent,
    # right part another) needs none; a haplotype equal to copy Z's reference sequence is Z's allele (Z becomes an option)
    maps_of, why = {}, collections.Counter()
    demand = []
    for i, m in enumerate(merged):
        if len(m["reads"]) < MIN_SUPPORT:
            why["below the floor"] += 1
            continue
        home = max(m["hap"], key=lambda c: (len(m["reads"] & piles[c].keys()), c))
        hp = m["hap"][home]
        if home not in maps_of:
            maps_of[home] = homology((home, chrom[home]) + span[home], [c for c in copies if c[0] != home], genome, tmpdir)
        mp = maps_of[home]
        parents = [(f"ref:{home}", {c: genome.fetch(chrom[home], c, c + 1).upper() for c in hp})]
        for cid, ch, _, _ in copies:
            if cid == home or cid not in mp:
                continue
            pb = {}
            for c in hp:
                if c in mp[cid]:
                    z, st = mp[cid][c]
                    b = genome.fetch(ch, z, z + 1).upper()
                    pb[c] = b if st == "+" else CMP.get(b, "N")
            parents.append((f"ref:{cid}", pb))
        for k, (pid, n, g) in enumerate(allhaps):
            if n <= len(m["reads"]):
                continue
            pb = {}
            for c in hp:
                if pid == home:
                    pb[c] = g.get(c, cons[pid].get(c))
                elif pid in mp and c in mp[pid]:
                    z, st = mp[pid][c]
                    b = g.get(z, cons[pid].get(z, genome.fetch(chrom[pid], z, z + 1).upper()))
                    pb[c] = b if st == "+" else CMP.get(b, "N")
            parents.append((f"hap{k}", pb))
        ex = explanation(hp, parents)
        if ex is None:
            demand.append(i)
            why["unexplained"] += 1
        elif ex[0] == "single" and ex[1].startswith("ref:"):
            z = ex[1][4:]
            if z not in options[i] and z in piles:
                options[i] = sorted(options[i] + [z])
            demand.append(i)
            why["equals a reference copy"] += 1
        else:
            why["duplicate" if ex[0] == "single" else "chimera"] += 1
    capacity = {c: 2 for c in piles}
    if slot_order == "support":
        order = sorted(demand, key=lambda i: (-len(merged[i]["reads"]), merged[i]["div"], i))
    else:
        order = sorted(demand, key=lambda i: (merged[i]["div"], -len(merged[i]["reads"]), i))
    match = b_matching({i: options[i] for i in demand}, capacity, order)
    certified = []
    for i in demand:
        m = merged[i]
        if match[i] is None and private_alleles(i, merged, match, options) >= MIN_READS:
            home = max(m["piles"], key=lambda c: (len(m["reads"] & piles[c].keys()), c))
            _, steps = observations({r: piles[home][r] for r in m["reads"] if r in piles[home]}, prim, span[home])
            certified.append((home, m["reads"], consensus(steps)))
    return certified, per_copy, len(merged), why


def isocon_side(prefix, net_path):
    """The IsoCon side (2026-10-23, design step 'reads that align to no copy'): the reference-free candidates of the shipped stage
    (`candidate_copies`, in-house IsoCon: clusters of the whole net -> consensus -> significance merge -> flag), kept when at
    least half of their member reads joined the net by attribution (net_reads.tsv via = B: unmapped or poorly placed reads
    with no record on a copy), the reads the Strainy side never sees. Returns {family: [(stage candidate, members, contig)]}."""
    via = {(r["family"], r["read"]): r["via"] for r in csv.DictReader(open(net_path), delimiter="\t")}
    cand = {r["candidate"]: r for r in csv.DictReader(open(f"{prefix}.candidates.tsv"), delimiter="\t")}
    owner = {r["cluster"]: r["candidate"] for r in csv.DictReader(open(f"{prefix}.clusters.tsv"), delimiter="\t")
             if r["candidate"] in cand}
    members = collections.defaultdict(set)
    for r in csv.DictReader(open(f"{prefix}.reads.tsv"), delimiter="\t"):
        if r["cluster"] in owner:
            members[owner[r["cluster"]]].add((r["family"], r["read"]))
    seqs, name = {}, None
    for ln in open(f"{prefix}.contigs.fa"):
        ln = ln.strip()
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(ln)
    out = collections.defaultdict(list)
    for c, r in sorted(cand.items()):
        m = members[c]
        if r["flagged"] == "1" and m and c in seqs and 2 * sum(via.get(x) == "B" for x in m) >= len(m):
            out[r["family"]].append((c, {x[1] for x in m}, "".join(seqs[c])))
    return out


def cmd_run(a):
    import pysam
    fams = set(a.families.split(",")) if a.families else None
    by_fam = load_copies(a.copies)
    net, n_net = load_net(a.net, fams)
    bam = pysam.AlignmentFile(a.bam)
    wanted = set().union(*net.values()) if net else set()
    prim = read_primaries(bam, wanted)
    genome = pysam.FastaFile(a.genome)
    iso = isocon_side(a.isocon, a.net) if a.isocon else {}
    import os
    import tempfile
    tmpdir = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(a.out)))
    with open(f"{a.out}.candidates.tsv", "w") as ct, open(f"{a.out}.contigs.fa", "w") as cf, \
            open(f"{a.out}.clusters.tsv", "w") as cl, open(f"{a.out}.reads.tsv", "w") as rd, \
            open(f"{a.out}.families.tsv", "w") as ft:
        ct.write("family\tcandidate\tn_clusters\tn_reads\tflagged\tunion_len\tnearest_locus\td\tn_net\tn_used\n")
        cl.write("family\tcluster\tcandidate\tn_reads\tconsensus_len\tfate\tlinked_to\td\n")
        rd.write("read\tfamily\tcluster\n")
        ft.write("family\tn_net\tn_used\tpiles\tpile_reads\tclusters\tconflicting\thaplotypes\tcertified\n")
        for fam in sorted(net):
            reads = net[fam]
            certified, per_copy, n_h, why = run_family(fam, by_fam[fam], reads, bam, prim, genome, a.slot_order, tmpdir)
            certified = certified + [(f"isocon:{c}", rs, seq) for c, rs, seq in iso.get(fam, [])]
            ft.write(f"{fam}\t{n_net[fam]}\t{len(reads)}\t{len(per_copy)}\t{sum(v[0] for v in per_copy.values())}\t"
                     f"{sum(v[1] for v in per_copy.values())}\t{sum(v[2] for v in per_copy.values())}\t{n_h}\t{len(certified)}\n")
            for k, (home, rs, seq) in enumerate(certified):
                cand = f"cand_{fam}_{k}"
                ct.write(f"{fam}\t{cand}\t1\t{len(rs)}\t1\t{len(seq)}\t{home}\tNA\t{n_net[fam]}\t{len(reads)}\n")
                cf.write(f">{cand}\n{seq}\n")
                cl.write(f"{fam}\t{fam}:h{k}\t{cand}\t{len(rs)}\t{len(seq)}\tnew_copy\t-\tNA\n")
                for r in sorted(rs):
                    rd.write(f"{r}\t{fam}\t{fam}:h{k}\n")
            print(f"[o3port] {fam}: reads {len(reads)}, piles {len(per_copy)}, haplotypes {n_h} {dict(why)}, "
                  f"certified {len(certified)}", file=sys.stderr, flush=True)


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--bam", required=True)
    r.add_argument("--copies", required=True)
    r.add_argument("--genome", required=True)
    r.add_argument("--net", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--families", default="")
    r.add_argument("--slot-order", choices=["divergence", "support"], default="divergence")
    r.add_argument("--isocon", default="", help="the stage's product prefix (candidate_copies) for the IsoCon side")
    a = ap.parse_args(argv)
    cmd_run(a)


if __name__ == "__main__":
    main(sys.argv[1:])
