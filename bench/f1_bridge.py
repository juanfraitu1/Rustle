#!/usr/bin/env python3
"""f1_bridge.py -- F1, bridge-aware regrouping of an assembled GTF (locus_fix_design; DESIGN PHASE, dev evidence only).

A post-processor of the shipped assembler's final GTF (after every polish). It never adds, drops or edits a transcript's
exons: it only rewrites `gene_id` (and appends two attributes to bridge transcripts' `transcript` lines), so intron
chains, `=` counts and chain precision are unchanged by construction.

RULE (binding for this design; every constant is inherited, source in brackets)
Transcripts: `transcript` lines of IN.gtf with their `exon` lines (1-based closed), strand, `reads`, the input gene_id.
Per input gene_id g and strand, for every intron J = [s, e] (1-based closed) used by >= 1 transcript of g:
  T_J = g's transcripts that use J; R_J = g's other transcripts on J's strand.
  Components of R_J under same-strand exon overlap (>= 1 shared base) [RG3, r1127].
  A component is UP when it starts before J's donor and ends before J's acceptor, DOWN when it starts at or after the
  donor side and ends beyond the acceptor, STRADDLE when it starts before the donor and ends beyond the acceptor
  (transcript orientation; '-' mirrored); a component wholly inside the intron is neither.
  STRUCTURAL(J) :<=> no STRADDLE component, >= 1 UP and >= 1 DOWN component (T_J is then the only link).
  UP-proof(J)   :<=> the U population of J [readthrough filter R, r1117: spliced primaries on J's strand whose 5' end
                    is upstream of J's donor and whose 3' end lies inside J's intron], deduplicated on (strand, start,
                    end, intron chain) [tss_measure / --polish-tes evidence], holds >= 1 PAS-PROVEN 3' cluster
                    [--polish-tes, r1132/r1133: single linkage of oriented 3' ends, gap <= 21 bp; >= 2 reads; mode =
                    most-ended position, ties 3'-most; proven = AATAAA or ATTAAA wholly inside oriented
                    [mode-35, mode-10] AND the mode not internally primed (>= 60% A in the 20 bp downstream, or A6)].
                    Reading of "at or before the bridge junction's intron": the upstream gene terminates INSIDE the
                    intron the bridge skips (the U population); exonic 3' piles are excluded because gorilla
                    internal-exon piles are PAS-proven at .125 (r1132).
  DOWN-proof(J) :<=> V1(J) >= 1 [readthrough filter's V1 = Q1's count, readthrough_rules.py verbatim: spliced
                    primaries, NOT deduplicated, whose 5' end is inside J's intron and in a real start cluster (5' ends
                    per strand, a gap > 100 bp opens a new cluster, >= 3 reads), whose 3' end is beyond J's acceptor,
                    and whose first donor (any read with the same strand and both ends) lies inside J's intron].
  BRIDGE(J)     :<=> STRUCTURAL(J) AND UP-proof(J) AND DOWN-proof(J).  (--mode struct|up|down: ablations that drop
                    the missing clauses; --mode nopas: the UP clause without the PAS (an unprimed cluster suffices);
                    --mode full is F1. The junctions table marks clusters P = proven, u = unprimed without PAS.)
Bridge transcripts B_g = union of T_J over g's bridge junctions (every J judged on the ORIGINAL locus; order-free).
New gene_ids: the non-bridge transcripts of g are split into same-strand exon-overlap components and named exactly as
RG3 names pieces (keeper = the piece whose representative max(reads, span, -index) is largest keeps g; the others
"<g>.rg<k>", k = 2.. in representative-index order). So without any bridge F1 = RG3 (rg3.py ec17e540).
Bridge transcripts: grouped by exon overlap among themselves; each group gets "<g>.fus<k>" (k = 1.. in index order);
their `transcript` lines get `fusion_of "<piece>,<piece>,..."` (the non-bridge pieces of g whose exons they overlap,
transcript order) and `fusion_junction "<chrom>:<s>-<e>:<strand>[,...]"`.
The bridges are RELATIONS, not loci: OUT.families.gtf = OUT.gtf minus every bridge transcript is what the families
stage reads, so a bridge is never a locus representative and never a family node.

Outputs: OUT.gtf, OUT.families.gtf, OUT.junctions.tsv (every STRUCTURAL junction with its evidence), OUT.stats.json.
usage: f1_bridge.py IN.gtf BAM FASTA OUT_PREFIX [--mode full|struct|up|down|none|nopas] [--contigs c1,c2]
(--mode none: no bridge at all, which must reproduce rg3.py byte for byte.)
"""
from __future__ import annotations

import bisect
import collections
import json
import sys

import pysam

W = 21                  # --polish-tes single-linkage gap (TSS_W = 2*TOL+1, TOL = 10)
PAS = ('AATAAA', 'ATTAAA')
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def attr(s: str, key: str):
    pat = f'{key} "'
    i = s.find(pat)
    if i < 0:
        return None
    i += len(pat)
    j = s.find('"', i)
    return s[i:j] if j >= 0 else None


class Tx:
    __slots__ = ('tid', 'gene', 'chrom', 'strand', 'reads', 'span', 'idx', 'exons')

    def introns(self):
        return [(a[1] + 1, b[0] - 1) for a, b in zip(self.exons, self.exons[1:])]


def parse(path: str) -> list:
    txs, by_id = [], {}
    for line in open(path):
        if not line or line[0] == '#':
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        tid = attr(f[8], 'transcript_id')
        if tid is None:
            continue
        if f[2] == 'transcript':
            t = Tx()
            t.tid, t.gene, t.chrom, t.strand = tid, attr(f[8], 'gene_id'), f[0], f[6]
            rv = attr(f[8], 'reads')
            t.reads = int(rv) if rv is not None and rv.lstrip('-').isdigit() else 0
            t.span = int(f[4]) - int(f[3]) + 1
            t.idx, t.exons = len(txs), []
            assert tid not in by_id, tid
            by_id[tid] = t
            txs.append(t)
        elif f[2] == 'exon':
            by_id[tid].exons.append((int(f[3]), int(f[4])))
    for t in txs:
        t.exons.sort()
    return txs


class UF:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.p[a] = b


def components(ts: list) -> list:
    """same-strand exon-overlap components (lists of Tx) of ts (one chrom); rg3._join_exon_overlap sweep."""
    uf = UF(len(ts))
    by = collections.defaultdict(list)
    for i, t in enumerate(ts):
        for a, b in t.exons:
            by[t.strand].append((a, b, i))
    for ex in by.values():
        ex.sort()
        mx, own = None, None
        for a, b, i in ex:
            if mx is not None and a <= mx:
                uf.union(i, own)
            if mx is None or b > mx:
                mx, own = b, i
    comp = collections.defaultdict(list)
    for i in range(len(ts)):
        comp[uf.find(i)].append(ts[i])
    return list(comp.values())


# ------------------------------------------------------------------------------------------------ read evidence
class Evidence:
    """one contig's spliced primaries: the U population (dedup) indexed by 3' end, and V1 per readthrough_rules.py."""

    def __init__(self, bam: pysam.AlignmentFile, fa: pysam.FastaFile, contig: str):
        self.seq = fa.fetch(contig).upper()
        seen = set()
        ends = collections.defaultdict(list)          # strand -> (e3, e5)  (all spliced primaries, V1)
        first_donor = collections.defaultdict(list)  # (strand, e5, e3) -> first donors
        u3 = collections.defaultdict(list)            # strand -> (e3, e5) deduplicated (U population)
        n = 0
        for rd in bam.fetch(contig):
            if rd.flag & 2308:
                continue
            pos, introns = rd.reference_start, []
            for op, ln in rd.cigartuples:
                if op == 3:
                    introns.append((pos + 1, pos + ln))
                    pos += ln
                elif op in (0, 2, 7, 8):
                    pos += ln
            if not introns:
                continue
            ts = rd.get_tag('ts') if rd.has_tag('ts') else '+'
            st = '+' if (ts == '+') != rd.is_reverse else '-'
            e3, e5 = (rd.reference_end, rd.reference_start + 1) if st == '+' else (rd.reference_start + 1, rd.reference_end)
            ends[st].append((e3, e5))
            first_donor[(st, e5, e3)].append(introns[0][0] - 1 if st == '+' else introns[-1][1] + 1)
            key = (st, rd.reference_start, rd.reference_end, tuple(introns))
            if key not in seen:
                seen.add(key)
                u3[st].append((e3, e5))
            n += 1
        self.n_spliced = n
        self.u3 = {st: sorted(v) for st, v in u3.items()}
        # V1's start clusters (readthrough_rules.py): 5' ends per strand, a gap > 100 opens a cluster, >= 3 = real
        self.v1 = {}
        for st, E in ends.items():
            S5 = sorted((e5, e3) for e3, e5 in E)
            lab, grp = [], []
            for i, (e5, _) in enumerate(S5):
                if grp and e5 - S5[grp[-1]][0] > 100:
                    lab += [len(grp) >= 3] * len(grp)
                    grp = []
                grp.append(i)
            lab += [len(grp) >= 3] * len(grp)
            self.v1[st] = (S5, lab)
        self.first_donor = first_donor

    # --- --polish-tes predicates (oriented coordinates: o = g on '+', o = -g on '-')
    def oseq(self, minus: bool, olo: int, ohi: int) -> str:
        a, b = (-ohi, -olo) if minus else (olo, ohi)
        i, j = max(0, a - 1), min(len(self.seq), b)
        if i >= j:
            return ''
        s = self.seq[i:j]
        return s.translate(COMP)[::-1] if minus else s

    def pas(self, minus: bool, c: int) -> bool:
        w = self.oseq(minus, c - 35, c - 10)
        return any(p in w for p in PAS)

    def primed(self, minus: bool, c: int) -> bool:
        d = self.oseq(minus, c + 1, c + 20)
        return bool(d) and (d.count('A') >= 0.6 * len(d) or 'AAAAAA' in d)

    def clusters(self, ends_o: list, minus: bool) -> list:
        ends_o = sorted(ends_o)
        out, i = [], 0
        while i < len(ends_o):
            j = i + 1
            while j < len(ends_o) and ends_o[j] - ends_o[j - 1] <= W:
                j += 1
            if j - i >= 2:
                mode, best, k = ends_o[i], 0, i
                while k < j:
                    e = k + 1
                    while e < j and ends_o[e] == ends_o[k]:
                        e += 1
                    if e - k >= best:
                        best, mode = e - k, ends_o[k]
                    k = e
                pr = self.primed(minus, mode)
                out.append((mode, j - i, self.pas(minus, mode) and not pr, not pr))
            i = j
        return out

    def up(self, s: int, e: int, st: str):
        """U population of intron [s, e] and its 3' clusters -> (U, clusters [(genomic mode, n, proven)])."""
        E = self.u3.get(st, [])
        lo, hi = bisect.bisect_left(E, (s, -1)), bisect.bisect_right(E, (e, 1 << 62))
        up = [e3 for e3, e5 in E[lo:hi] if (e5 < s if st == '+' else e5 > e)]
        minus = st == '-'
        cl = self.clusters([-x if minus else x for x in up], minus)
        return len(up), [(-m if minus else m, n, p, u) for m, n, p, u in cl]

    def v1_count(self, s: int, e: int, st: str) -> int:
        S5, lab = self.v1.get(st, ([], []))
        a, b = bisect.bisect_left(S5, (s, -1)), bisect.bisect_right(S5, (e, 1 << 62))
        n = 0
        for i in range(a, b):
            e5, e3 = S5[i]
            if not lab[i]:
                continue
            if (st == '+' and e3 > e) or (st == '-' and e3 < s):
                if any(s <= d <= e for d in self.first_donor[(st, e5, e3)]):
                    n += 1
        return n


# ------------------------------------------------------------------------------------------------ the rule
def side(comp: list, s: int, e: int, st: str) -> str:
    lo = min(t.exons[0][0] for t in comp)
    hi = max(t.exons[-1][1] for t in comp)
    before, after = lo < s, hi > e             # genomic: has bases left of the intron / right of it
    if before and after:
        return 'STRADDLE'
    if not before and not after:
        return 'INSIDE'
    if st == '+':
        return 'UP' if before else 'DOWN'
    return 'DOWN' if before else 'UP'


def decide(txs: list, ev_of, mode: str = 'full'):
    """-> (bridge tids, junction rows)."""
    by_g = collections.defaultdict(list)
    for t in txs:
        by_g[(t.gene, t.chrom)].append(t)
    bridges, rows = set(), []
    for (g, chrom), ts in by_g.items():
        if len(ts) < 3:           # a bridge needs >= 1 UP, >= 1 DOWN and >= 1 linking transcript
            continue
        J = collections.defaultdict(list)
        for t in ts:
            for iv in t.introns():
                J[(iv, t.strand)].append(t)
        for ((s, e), st), TJ in sorted(J.items()):
            tj = {t.tid for t in TJ}
            R = [t for t in ts if t.strand == st and t.tid not in tj]
            if len(R) < 2:
                continue
            sides = collections.Counter(side(c, s, e, st) for c in components(R))
            if sides['STRADDLE'] or not sides['UP'] or not sides['DOWN']:
                continue
            ev = ev_of(chrom)
            U, cl = ev.up(s, e, st)
            up_ok = any(p for _, _, p, _ in cl)
            up_np = any(u for _, _, _, u in cl)       # ablation noPAS: an unprimed cluster, PAS not required
            v1 = ev.v1_count(s, e, st)
            down_ok = v1 >= 1
            is_b = {'full': up_ok and down_ok, 'struct': True, 'up': up_ok, 'down': down_ok, 'none': False,
                    'nopas': up_np and down_ok}[mode]
            if is_b:
                bridges |= tj
            rows.append(dict(gene=g, chrom=chrom, s=s, e=e, strand=st, n_TJ=len(TJ), reads_TJ=sum(t.reads for t in TJ),
                             up=sides['UP'], down=sides['DOWN'], inside=sides['INSIDE'], U=U,
                             clusters=';'.join(f'{m}:{n}:{"P" if p else ("u" if u else "-")}' for m, n, p, u in cl) or '.',
                             up_proof=up_ok, V1=v1, down_proof=down_ok, bridge=is_b,
                             TJ=','.join(sorted(tj))))
    return bridges, rows


def regroup(txs: list, bridges: set):
    """-> ({tid: new gene_id}, {tid: (fusion_of, fusion_junction)}), RG3 naming on the non-bridge transcripts."""
    by_g = collections.defaultdict(list)
    for t in txs:
        by_g[t.gene].append(t)
    new, rel = {}, {}
    for g, ts in by_g.items():
        nb = [t for t in ts if t.tid not in bridges]
        comps = []
        by_c = collections.defaultdict(list)
        for t in nb:
            by_c[t.chrom].append(t)
        for c, cts in by_c.items():
            comps += components(cts)
        reps = [max(c, key=lambda t: (t.reads, t.span, -t.idx)) for c in comps]
        if comps:
            keep = max(range(len(comps)), key=lambda i: (reps[i].reads, reps[i].span, -reps[i].idx))
            order = sorted(range(len(comps)), key=lambda i: reps[i].idx)
            k = 2
            name = {}
            for i in order:
                if i == keep:
                    name[i] = g
                else:
                    name[i] = f'{g}.rg{k}'
                    k += 1
            for i, c in enumerate(comps):
                for t in c:
                    new[t.tid] = name[i]
        bt = [t for t in ts if t.tid in bridges]
        if bt:
            bcomps = sorted(components(bt), key=lambda c: min(t.idx for t in c))
            for k, c in enumerate(bcomps, 1):
                bname = f'{g}.fus{k}'
                for t in c:
                    new[t.tid] = bname
                    pieces = []
                    for i, pc in enumerate(comps):
                        if any(x.strand == t.strand and any(a <= d and c0 <= b for a, b in x.exons for c0, d in t.exons)
                               for x in pc):
                            pieces.append((min(x.exons[0][0] for x in pc), name[i]))
                    pieces.sort(reverse=(t.strand == '-'))
                    jn = ','.join(f'{t.chrom}:{s}-{e}:{t.strand}' for s, e in t.introns() if (s, e, t.strand) in BJ)
                    rel[t.tid] = (','.join(n for _, n in pieces), jn)
    return new, rel


BJ = set()


def rewrite(src: str, dst: str, dst_fam: str, new: dict, rel: dict) -> dict:
    n_changed = n_fam_dropped = 0
    with open(src) as fi, open(dst, 'w') as fo, open(dst_fam, 'w') as ff:
        for line in fi:
            out = line
            tid = None
            if line and line[0] != '#':
                f = line.rstrip('\n').split('\t')
                if len(f) >= 9:
                    tid = attr(f[8], 'transcript_id')
                    old = attr(f[8], 'gene_id')
                    nv = new.get(tid) if tid else None
                    if nv is not None and old is not None and nv != old:
                        f[8] = f[8].replace(f'gene_id "{old}"', f'gene_id "{nv}"', 1)
                        n_changed += 1
                    if tid in rel and f[2] == 'transcript':
                        fo_, fj = rel[tid]
                        f[8] = f[8].rstrip() + f' fusion_of "{fo_}"; fusion_junction "{fj}";'
                    out = '\t'.join(f) + '\n'
            fo.write(out)
            if tid is not None and tid in rel:
                n_fam_dropped += 1
                continue
            ff.write(out)
    return dict(lines_changed=n_changed, family_lines_dropped=n_fam_dropped)


def main(argv):
    if len(argv) < 4:
        print(__doc__)
        return 2
    src, bam_p, fa_p, out = argv[:4]
    mode = argv[argv.index('--mode') + 1] if '--mode' in argv else 'full'
    contigs = set(argv[argv.index('--contigs') + 1].split(',')) if '--contigs' in argv else None
    txs = parse(src)
    if contigs is not None:
        assert all(t.chrom in contigs for t in txs), 'the GTF holds a contig outside --contigs'
    bam, fa = pysam.AlignmentFile(bam_p), pysam.FastaFile(fa_p)
    cache = {}

    def ev_of(c):
        if c not in cache:
            cache.clear()
            cache[c] = Evidence(bam, fa, c)
        return cache[c]
    txs_sorted = sorted(txs, key=lambda t: t.chrom)          # one contig's evidence at a time
    bridges, rows = decide(txs_sorted, ev_of, mode)
    for r in rows:
        if r['bridge']:
            BJ.add((r['s'], r['e'], r['strand']))
    new, rel = regroup(txs, bridges)
    st = rewrite(src, out + '.gtf', out + '.families.gtf', new, rel)
    cols = ['gene', 'chrom', 's', 'e', 'strand', 'n_TJ', 'reads_TJ', 'up', 'down', 'inside', 'U', 'clusters', 'up_proof',
            'V1', 'down_proof', 'bridge', 'TJ']
    with open(out + '.junctions.tsv', 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')
    genes_in = {t.gene for t in txs}
    split = {t.gene for t in txs if new[t.tid] != t.gene and t.tid not in bridges}
    bgenes = {t.gene for t in txs if t.tid in bridges}
    stats = dict(rule=f'F1 mode={mode}', transcripts=len(txs), gene_ids=len(genes_in),
                 structural_junctions=len(rows), up_proof=sum(r['up_proof'] for r in rows),
                 down_proof=sum(r['down_proof'] for r in rows), bridge_junctions=sum(r['bridge'] for r in rows),
                 bridge_transcripts=len(bridges), gene_ids_with_bridge=len(bgenes),
                 gene_ids_split=len(split | bgenes), gene_ids_after=len(set(new.values())),
                 families_gene_ids=len({new[t.tid] for t in txs if t.tid not in bridges}), **st)
    json.dump(stats, open(out + '.stats.json', 'w'), indent=1)
    print(json.dumps(stats))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
