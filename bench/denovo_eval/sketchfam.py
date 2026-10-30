#!/usr/bin/env python3
"""Exon-only FracMinHash sketch families (PENDING 1c, 2026-10-29, user "lets try the sketch"). For every de novo locus of a
`<p>.fam.loci.gff3` take the EXON sequences only, sketch their canonical k-mers (FracMinHash: keep hash < 2^64 / scale), discard every
hash seen in more than --max-loci loci (the k-mer form of the dispersion rule: TE k-mers sit in hundreds of loci, true families in a few
dozen), join two loci when they share >= --min-shared retained hashes AND containment (shared / smaller sketch) >= --min-cont, families =
connected components of that graph (>= 2 loci). No minimap2, no gate. Writes OUT.fam.copies.tsv (famtruth/pairjudge format),
OUT.edges.tsv (a b shared containment masked_share_of_shared) and OUT.sketch.log.

  sketchfam.py LOCI_GFF3 GENOME_FA OUT [--k 21] [--scale 20] [--max-loci 0] [--max-chrom 0 --fai GENOME_FAI] [--min-shared 5] [--min-cont 0.2] [--unmasked-only] [--allow-overlap] [--cluster components|greedy|tree --family-cont 0.3 --root-cont 0.1 --cut fixed|relative|gap --rel-frac 0.5 --gap-min 0.15] [--weak-cont 0.1 --intron-cont 0.2 --intron-min-hashes 10] [--merge-overlap] [--metric containment|jaccard]
Needs numpy and samtools (per-chromosome fetch).
"""
import sys, re, collections, subprocess, time
import numpy as np

a = sys.argv[1:]
def opt(name, default, cast=float):
    return cast(a[a.index(name) + 1]) if name in a else default
GFF, FA, OUT = a[0], a[1], a[2]
K = opt('--k', 21, int); SCALE = opt('--scale', 20, int); MAXL = opt('--max-loci', 0, int)
MAXC = opt('--max-chrom', 0, int); FAI = a[a.index('--fai') + 1] if '--fai' in a else None; MINBP = opt('--min-chrom-bp', 10_000_000, int)
MINS = opt('--min-shared', 5, int); MINC = opt('--min-cont', 0.2); CLUSTER = a[a.index('--cluster') + 1] if '--cluster' in a else 'components'; FAMC = opt('--family-cont', 0.3); ROOTC = opt('--root-cont', 0.1); CUT = a[a.index('--cut') + 1] if '--cut' in a else 'fixed'; RELF = opt('--rel-frac', 0.5); GAPMIN = opt('--gap-min', 0.15); WEAKC = opt('--weak-cont', 0.0); INTRC = opt('--intron-cont', 0.2); INTRMIN = opt('--intron-min-hashes', 10, int); MERGE = '--merge-overlap' in a; METRIC = a[a.index('--metric') + 1] if '--metric' in a else 'containment'; UNMASKED = '--unmasked-only' in a; OVERLAP_OK = '--allow-overlap' in a
t0 = time.time()
log = open(OUT + '.sketch.log', 'w')
def say(*x):
    print(f"[sketchfam {time.time()-t0:6.0f}s]", *x, file=sys.stderr); print(*x, file=log); log.flush()

# loci: tid -> (chrom, strand, [(s0, e)])
loci = {}; exons = collections.defaultdict(list)
for line in open(GFF):
    if line[0] == '#': continue
    f = line.rstrip('\n').split('\t')
    if f[2] == 'gene':
        loci[re.search(r'ID=gene-([^;]+)', f[8]).group(1)] = (f[0], f[6] if f[6] in '+-' else '+')
    elif f[2] == 'exon':
        exons[re.search(r'gene=([^;]+)', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
if MERGE:
    # LOCUS MERGE (2026-10-29, user "can we ensure that all mainstays are correctly and specifically recalled?"): on A119b 36% of the
    # de novo loci overlap an earlier same-strand locus (TBC1D3: a 14-exon locus and a 2-exon intron-retaining locus at every copy,
    # each set forming its own family; NPIPB4: six overlapping 2-exon loci). Same-strand loci whose spans overlap are one locus:
    # union of their exons (merged intervals), name of the first. Sub-gene families and duplicate families disappear.
    by_cs = collections.defaultdict(list)
    for t, (c, s_) in loci.items():
        if exons[t]: by_cs[(c, s_)].append((min(a for a, _ in exons[t]), max(b for _, b in exons[t]), t))
    n_before = len(loci); loci2 = {}; exons2 = {}
    for (c, s_), v in by_cs.items():
        v.sort(); cur = None
        for a, b, t in v:
            if cur and a < cur[1]:
                cur[1] = max(cur[1], b); cur[2] += exons[t]
            else:
                if cur: loci2[cur[3]] = (c, s_); exons2[cur[3]] = cur[2]
                cur = [a, b, list(exons[t]), t]
        if cur: loci2[cur[3]] = (c, s_); exons2[cur[3]] = cur[2]
    for t, ex in exons2.items():
        m = []
        for a, b in sorted(ex):
            if m and a <= m[-1][1]: m[-1][1] = max(m[-1][1], b)
            else: m.append([a, b])
        exons2[t] = [tuple(x) for x in m]
    loci = loci2; exons = collections.defaultdict(list, exons2)
    say(f"merge-overlap: {n_before} loci -> {len(loci)}")
say(f"loci {len(loci)}, exons {sum(len(v) for v in exons.values())}, k {K} scale {SCALE} max-loci {MAXL} min-shared {MINS} min-cont {MINC}")

LUT = np.full(256, 4, dtype=np.uint8)
for c, v in zip(b'ACGTacgt', [0, 1, 2, 3, 0, 1, 2, 3]): LUT[c] = v
LOW = np.zeros(256, dtype=bool)
for c in b'acgtn': LOW[c] = True
MASK = np.uint64((1 << (2 * K)) - 1); THRESH = np.uint64((1 << 64) // SCALE)
def mix(h):  # splitmix64 finaliser
    h = (h ^ (h >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    h = (h ^ (h >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    return h ^ (h >> np.uint64(31))
def kmer_hashes(seq_u8, low):
    """Per k-mer START position of one 2-bit array: splitmix hash of the canonical k-mer (0 where the k-mer has an N, or is
    soft-masked under --unmasked-only), and the masked flag (k-mer touches a soft-masked base). Chromosome-wide, once."""
    n = len(seq_u8) - K + 1
    if n <= 0: return np.empty(0, dtype=np.uint64), np.empty(0, dtype=bool)
    codes = seq_u8.astype(np.uint64)
    fwd = np.zeros(n, dtype=np.uint64); rev = np.zeros(n, dtype=np.uint64); bad = np.zeros(n, dtype=bool); msk = np.zeros(n, dtype=bool)
    for j in range(K):
        w = codes[j:j + n]
        bad |= w > 3; msk |= low[j:j + n]
        fwd = (fwd << np.uint64(2)) | (w & np.uint64(3))
        rev |= ((np.uint64(3) - (w & np.uint64(3))) << np.uint64(2 * j))
    if UNMASKED: bad |= msk
    h = mix(np.minimum(fwd, rev)); h[bad] = np.uint64(0); h[h >= THRESH] = np.uint64(0)  # 0 = not in the sketch
    return h, msk
SEP = np.full(K - 1, 4, dtype=np.uint8); SEPL = np.zeros(K - 1, dtype=bool)
def sketch_pieces(arr, lowarr, ranges):
    """FracMinHash of the k-mers inside each [s, e) piece, ONE numpy pass per locus: the pieces are concatenated with K-1 N's
    between them, so every k-mer crossing a piece boundary contains an N and is dropped. Only the pieces' bp are hashed."""
    ranges = [(s0, e) for s0, e in ranges if e - s0 >= K]
    if not ranges: return np.empty(0, dtype=np.uint64), np.empty(0, dtype=bool)
    parts = []; lows = []
    for s0, e in ranges: parts += [arr[s0:e], SEP]; lows += [lowarr[s0:e], SEPL]
    h, m = kmer_hashes(np.concatenate(parts[:-1]), np.concatenate(lows[:-1]))
    nz = h != 0; h, m = h[nz], m[nz]
    u, first = np.unique(h, return_index=True)
    return u, m[first]
def sketch_ranges(hpos, mpos, ranges):
    """FracMinHash of the k-mers fully inside the given [s, e) ranges (k-mers never cross a range boundary): gather the
    chromosome-wide hashes at positions s .. e-K, drop zeros, unique."""
    idx = np.concatenate([np.arange(s0, e - K + 1) for s0, e in ranges if e - s0 >= K]) if ranges else np.empty(0, dtype=np.int64)
    if len(idx) == 0: return np.empty(0, dtype=np.uint64), np.empty(0, dtype=bool)
    h = hpos[idx]; m = mpos[idx]; nz = h != 0; h, m = h[nz], m[nz]
    u, first = np.unique(h, return_index=True)
    return u, m[first]

with np.errstate(over='ignore'):
    by_chrom = collections.defaultdict(list)
    for tid, (c, s) in loci.items(): by_chrom[c].append(tid)
    sk = {}; masked = {}; isk = {}
    for chrom, tids in by_chrom.items():
        seq = subprocess.run(['samtools', 'faidx', FA, chrom], capture_output=True, text=True).stdout
        seq = ''.join(seq.split('\n')[1:]).encode()
        raw = np.frombuffer(seq, dtype=np.uint8); arr = LUT[raw]; lowarr = LOW[raw]
        for tid in tids:
            # hash the LOCUS SPAN once (chromosome-wide hashing is too slow and memory-heavy on 3 Gb), then gather each exon's
            # range: k-mers never cross an exon boundary
            ex = sorted(exons[tid])
            if not ex: sk[tid] = np.empty(0, dtype=np.uint64); masked[tid] = np.empty(0, dtype=bool); continue
            sk[tid], masked[tid] = sketch_pieces(arr, lowarr, ex)
            if WEAKC:
                # UNSPLICED signal (2026-10-29 user "could there be a method for unspliced and one for spliced?"): sketch of the
                # unmasked intron k-mers (gaps between consecutive exons); a WEAK exon edge (containment in [weak_cont, min_cont))
                # is accepted only when the two loci's intron sketches share >= intron_cont (young segmental duplicates copy their
                # introns; old paralogs and unrelated adjacent genes do not).
                h, m = sketch_pieces(arr, lowarr, [(e0, s1) for (s0, e0), (s1, e1) in zip(ex, ex[1:]) if s1 - e0 >= K])
                isk[tid] = h[~m]
    say(f"sketched {len(sk)} loci, mean hashes/locus {np.mean([len(v) for v in sk.values()]):.1f}, empty {sum(1 for v in sk.values() if len(v)==0)}")
    # high-copy removal: by locus count (--max-loci) and/or by CHROMOSOME DISPERSION (--max-chrom, r1471 lesson: a 46-copy array on
    # one chromosome is a family; a hash on 8 chromosomes is a TE). Sequences < --min-chrom-bp count as one 'unplaced' chromosome.
    allh = np.concatenate([v for v in sk.values() if len(v)])
    uh, cnt = np.unique(allh, return_counts=True)
    hot = uh[cnt > MAXL] if MAXL else np.empty(0, dtype=np.uint64)
    say(f"distinct hashes {len(uh)}, in > {MAXL} loci: {len(hot)} (max copies {cnt.max()})")
    if MAXC:
        big = {l.split('\t')[0] for l in open(FAI) if int(l.split('\t')[1]) >= MINBP} if FAI else set(by_chrom)
        cid = {c: (i if c in big else -1) for i, c in enumerate(by_chrom)}
        hs = []; cs = []
        for tid, v in sk.items():
            if len(v): hs.append(v); cs.append(np.full(len(v), cid[loci[tid][0]], dtype=np.int64))
        hs = np.concatenate(hs); cs = np.concatenate(cs)
        pairs_hc = np.unique(np.stack([hs, cs.astype(np.uint64)], axis=1), axis=0)
        uh2, nchrom = np.unique(pairs_hc[:, 0], return_counts=True)
        hot2 = uh2[nchrom > MAXC]
        say(f"hashes on > {MAXC} chromosomes: {len(hot2)} (max {nchrom.max()})")
        hot = np.union1d(hot, hot2)
    if len(hot):
        order = [t for t in sk if len(sk[t])]; allv = np.concatenate([sk[t] for t in order]); keep = ~np.isin(allv, hot)
        off = 0
        for t in order:
            n = len(sk[t]); k = keep[off:off + n]; sk[t] = sk[t][k]; masked[t] = masked[t][k]; off += n
    say(f"masked share of retained hashes {np.mean(np.concatenate([masked[t] for t in sk if len(masked[t])])):.2f}")
    # pair counts, VECTORISED (2026-10-30, user "lets do it" on speed: the Python inverted-index loop took > 1 h on 76k human loci at
    # scale 10). Exact, same pairs and counts as the loop: sort (hash, locus) rows, emit every within-hash pair by group size.
    tids = sorted(sk); tix = {t: i for i, t in enumerate(tids)}
    H = np.concatenate([sk[t] for t in tids]); L = np.concatenate([np.full(len(sk[t]), tix[t], dtype=np.int64) for t in tids])
    M = np.concatenate([masked[t] for t in tids])
    o = np.lexsort((L, H)); H, L, M = H[o], L[o], M[o]
    starts = np.flatnonzero(np.r_[True, H[1:] != H[:-1]]); sizes = np.diff(np.r_[starts, len(H)])
    keys = []; mflag = []
    for sz in np.unique(sizes[sizes >= 2]):
        st = starts[sizes == sz]
        ii, jj = np.triu_indices(sz, 1)
        a = (st[:, None] + ii[None, :]).ravel(); b = (st[:, None] + jj[None, :]).ravel()
        keys.append(L[a] * len(tids) + L[b]); mflag.append(M[a] | M[b])
    pair = collections.Counter(); pairm = collections.Counter()
    if keys:
        keys = np.concatenate(keys); mflag = np.concatenate(mflag)
        uk, inv, cnt = np.unique(keys, return_inverse=True, return_counts=True)
        mcnt = np.bincount(inv, weights=mflag).astype(np.int64)
        for k, n, mc in zip(uk.tolist(), cnt.tolist(), mcnt.tolist()):
            p = (tids[k // len(tids)], tids[k % len(tids)]); pair[p] = n
            if mc: pairm[p] = mc
    span = {t: (loci[t][0], min(s0 for s0, _ in exons[t]), max(e for _, e in exons[t])) for t in sk}
    def overlapping(x, y):
        a, b = span[x], span[y]
        return a[0] == b[0] and a[1] < b[2] and b[1] < a[2]
    edges = []; n_overlap = 0; n_weak = 0
    for (x, y), n in pair.items():
        if n < MINS: continue
        if not OVERLAP_OK and overlapping(x, y):  # the base's reject_overlapping: antisense / nested loci share k-mers, not ancestry
            n_overlap += 1; continue
        c = n / max(1, min(len(sk[x]), len(sk[y]))) if METRIC == 'containment' else n / max(1, len(sk[x]) + len(sk[y]) - n)  # jaccard (2026-10-30, user: guided mode, full-length templates)
        if c >= MINC: edges.append((x, y, n, c, pairm[(x, y)] / n))
        elif WEAKC and c >= WEAKC:
            ix, iy = isk.get(x), isk.get(y)
            if ix is not None and iy is not None and min(len(ix), len(iy)) >= INTRMIN:
                ic = len(np.intersect1d(ix, iy, assume_unique=True)) / min(len(ix), len(iy))
                if ic >= INTRC: edges.append((x, y, n, c, pairm[(x, y)] / n)); n_weak += 1
    say(f"pairs sharing >= 1 hash {len(pair)}, overlapping pairs rejected {n_overlap}, edges {len(edges)} (weak exon edges rescued by introns {n_weak})")
    if CLUSTER == 'greedy':
        # Family = UNION SKETCH (linclust / CD-HIT style, 2026-10-29 user: "faithfully encapsulates gene families, ok to lose the most
        # divergent"): loci in decreasing sketch size; a locus joins the first family whose union sketch contains >= --family-cont of
        # its hashes (ties: largest containment), among families it has an edge into; else it founds one. No transitive chaining.
        nbr = collections.defaultdict(set)
        for x, y, _, _, _ in edges: nbr[x].add(y); nbr[y].add(x)
        fam_of = {}; union = {}; members = {}
        for t in sorted(sk, key=lambda t: -len(sk[t])):
            if len(sk[t]) == 0: continue
            best = None; bc = 0.0; mine = set(sk[t].tolist())
            for f in {fam_of[n] for n in nbr[t] if n in fam_of}:
                c = len(mine & union[f]) / len(mine)
                if c >= FAMC and c > bc: best, bc = f, c
            if best is None:
                best = t; union[best] = set(); members[best] = []
            fam_of[t] = best; union[best] |= mine; members[best].append(t)
        fams = sorted([v for v in members.values() if len(v) >= 2], key=lambda v: -len(v))
    elif CLUSTER == 'tree':
        # Family = a TREE (2026-10-29 user: "one that respects clades"): per connected component, average-linkage (UPGMA) on
        # distance 1 - containment (unlinked pairs = 1), families = cut at 1 - --root-cont, clades = the subtrees; Newick per family
        # in OUT.trees.nwk so NPIPA / NPIPB are one family with two clades instead of two families.
        from scipy.cluster.hierarchy import linkage, fcluster, to_tree
        from scipy.spatial.distance import squareform
        parent = {t: t for t in sk}
        def find(x):
            while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for x, y, _, _, _ in edges: parent[find(x)] = find(y)
        comp = collections.defaultdict(list)
        for t in sk: comp[find(t)].append(t)
        cedges = collections.defaultdict(list)  # per-component edge list (iterating all edges per component was O(C x E))
        for x, y, _, c, _ in edges: cedges[find(x)].append((x, y, c))
        fams = []; nwk = open(OUT + '.trees.nwk', 'w')
        def newick(root, names):  # iterative (a 10k-leaf component overflowed the recursive version, 2026-10-30)
            out = []; stack = [(root, 0)]
            while stack:
                node, state = stack.pop()
                if node.is_leaf(): out.append(names[node.id]); continue
                if state == 0: out.append('('); stack.append((node, 1)); stack.append((node.get_left(), 0))
                elif state == 1: out.append(','); stack.append((node, 2)); stack.append((node.get_right(), 0))
                else: out.append(f"):{node.dist:.3f}")
            return ''.join(out)
        for root, members_ in comp.items():
            if len(members_) < 2: continue
            m = sorted(members_); n = len(m); D = np.ones((n, n))
            np.fill_diagonal(D, 0.0)
            ix = {t: i for i, t in enumerate(m)}
            for x, y, c in cedges[root]:
                D[ix[x], ix[y]] = D[ix[y], ix[x]] = 1.0 - c
            Z = linkage(squareform(D, checks=False), method='average')
            # Per-family relative cut (2026-10-29 goal "replace the global containment cut with a per-family relative cut"): nested
            # families (NPIP, GOLGA8, KRAB-ZNF) span containment 0.2-1.0 inside one component while flat arrays sit at 0.9-1.0, so a
            # fixed root containment splits the former and is irrelevant to the latter.
            #   fixed:    cut at 1 - root_cont (the old behaviour)
            #   relative: cut at containment rel_frac x the component's median edge containment, never below root_cont
            #   gap:      cut at the largest gap between consecutive merge heights if it exceeds gap_min, else at 1 - root_cont
            if CUT == 'relative':
                med = float(np.median([c for _, _, c in cedges[root]])) if cedges[root] else ROOTC
                t_cut = 1.0 - max(ROOTC, RELF * med)
            elif CUT == 'gap':
                h = np.sort(Z[:, 2]); h = h[h >= 0.2]; t_cut = 1.0 - ROOTC  # gaps among the deeper merges only (identical copies at 0 are not a gap)
                if len(h) >= 2:
                    gaps = np.diff(h); i = int(np.argmax(gaps))
                    if gaps[i] >= GAPMIN: t_cut = min(t_cut, float((h[i] + h[i + 1]) / 2))
            else:
                t_cut = 1.0 - ROOTC
            lab = fcluster(Z, t=t_cut, criterion='distance')
            groups = collections.defaultdict(list)
            for t, l in zip(m, lab): groups[l].append(t)
            for g in groups.values():
                if len(g) >= 2: fams.append(g)
            if n >= 3: nwk.write(newick(to_tree(Z), m) + ';\n')
        nwk.close()
        fams = sorted(fams, key=lambda v: -len(v))
    else:
        parent = {t: t for t in sk}
        def find(x):
            while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for x, y, _, _, _ in edges: parent[find(x)] = find(y)
        comp = collections.defaultdict(list)
        for t in sk: comp[find(t)].append(t)
        fams = sorted([v for v in comp.values() if len(v) >= 2], key=lambda v: -len(v))
    say(f"families {len(fams)}, members {sum(len(v) for v in fams)}, largest {[len(v) for v in fams[:5]]}")
    with open(OUT + '.edges.tsv', 'w') as w:
        for x, y, n, c, m in edges: w.write(f"{x}\t{y}\t{n}\t{c:.3f}\t{m:.2f}\n")
    with open(OUT + '.fam.copies.tsv', 'w') as w:
        w.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
        for i, v in enumerate(fams):
            for j, t in enumerate(sorted(v)):
                c, s = loci[t]; ex = sorted(exons[t])
                w.write(f"SK{i}\t{j}\t{t}\t{c}\t{ex[0][0]}\t{ex[-1][1]}\t{len(ex)}\t{s}\t0\t{','.join(f'{a}-{b}' for a, b in ex)}\n")
say("done")
