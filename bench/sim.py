#!/usr/bin/env python3
"""Read simulators with truth in the read names (wave 7, 2026-09-24; old scripts at git tag `notebook-2026-09-24`).

Old -> new (every old command keeps its arguments):
  ideal_chromosome_sim.py --gff G --genome FA --chrom C --out OUT.fa --arm ideal|trunc|rt [...]
                                         -> sim.py chromosome (same flags)
  missing_copy_sim.py MODE REF.gtf REF.fa CHROM K DIV OUT SEED
                                         -> sim.py missing-copy MODE REF.gtf REF.fa CHROM K DIV OUT SEED [--threads 4]
  /mnt/linuxdisk/tmp/gw22/o3/simB.py REF.gtf FA K DIV OUT SEED   (out of repo; chr20 hard-coded)
                                         -> sim.py missing-copy genomic REF.gtf FA chr20 K DIV OUT SEED
  tandem_copy_sim.py --fasta ... [flags] -> sim.py tandem (same flags)
  copy_assign_read_truth.py sim COPIES.tsv COPIES.fa INDEX.mmi OUT SEED
                                         -> sim.py copies COPIES.tsv COPIES.fa INDEX.mmi OUT SEED [--threads 4]
  copy_assign_excision.py FAM_DIR X OUT  (usage line said o2_excision.py; BAM/FASTA/BIN hard-coded)
                                         -> sim.py excise FAM_DIR X OUT [--bam ... --fasta ... --bin ... --threads 4]
  from sim_reads import simulate_reads   -> from sim import simulate_reads   (sim_reads.write_fastq -> sim.write_fastq)
  copy_assign_read_truth.py score ...    -> score.py reads (not here)

⚠ B2 (fixed here): `missing-copy` and `copies` used to seed per-gene / per-copy RNGs with Python's `hash(str)`, which is
salted per process (PYTHONHASHSEED unset), so no two runs produced the same reads and no earlier run can be
reproduced. They now use `stable_seed()` (zlib.crc32 of the parts joined by an unambiguous separator; the old `fam + ci`
concatenation also collided: "FAM1"+"12" == "FAM11"+"2"). **Their reads therefore differ from every earlier run**;
numbers measured on those runs (REPRODUCE.md O2 read-truth "Expected", the O3 40/40 positive control) must be
re-measured, not compared read for read. `chromosome` and `tandem` never used hash() and are unchanged.

Only the standard library is imported at module top; pysam is imported inside the subcommands.
"""
import argparse
import collections
import csv
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402

BASES = "ACGT"
# the shipped read mapping (every simulator's reads are mapped exactly as the pipeline's real reads are)
MM2 = "minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes"


# ================================================================ read model (was sim_reads.py)
def simulate_reads(seq, n, err=0.003, indel=0.001, seed=0, trunc_frac=0.0):
    """Shared full-length HiFi transcript-read simulator (deterministic). IsoSeq reads are full-length,
    high-accuracy transcripts, so each read = one transcript sequence + HiFi errors (substitutions +
    rare short indels), with optional 5'/3' degradation truncation. A KNOWN error rate is what the
    identifiability theorem needs (PSVs must clear the error floor). No wall-clock RNG.

    n reads from `seq` with per-base substitution rate `err` and indel rate `indel` each.
    trunc_frac: max fraction trimmed from a random end (models IsoSeq 5'/3' degradation).
    Returns list of (read_seq, qual_str). Deterministic given (seq, n, seed)."""
    out = []
    for i in range(n):
        rng = random.Random((seed * 1_000_003) ^ (i * 2_654_435_761) ^ len(seq))
        s = seq
        if trunc_frac > 0:
            t = int(rng.random() * trunc_frac * len(s))
            if t:
                if rng.random() < 0.5:
                    s = s[t:]
                else:
                    s = s[:len(s) - t]
        buf = []
        for ch in s:
            r = rng.random()
            if r < err:
                buf.append(rng.choice([b for b in BASES if b != ch]))   # substitution
            elif r < err + indel:
                continue                                                # deletion
            elif r < err + 2 * indel:
                buf.append(ch); buf.append(rng.choice(BASES))           # insertion
            else:
                buf.append(ch)
        rd = "".join(buf)
        out.append((rd, "~" * len(rd)))   # '~' = Q93 placeholder (HiFi-grade)
    return out


def write_fastq(fh, name, read_qual):
    rd, q = read_qual
    fh.write(f"@{name}\n{rd}\n+\n{q}\n")


def stable_seed(*parts):
    """A per-process-STABLE integer from string parts (B2 fix; replaces hash(), which Python salts per process).
    zlib.crc32 of the parts joined by U+001F, so ('FAM1', '12') and ('FAM11', '2') differ."""
    return zlib.crc32('\x1f'.join(parts).encode())


def jitter(rd, q, rng, maxj=30):
    """Trim 0..maxj bases from each end (MANDATORY: identical reads collapse under dedup, §6n0)."""
    j5, j3 = rng.randint(0, maxj), rng.randint(0, maxj)
    return rd[j5:len(rd) - j3], q[j5:len(q) - j3]


def mutate(seq, rate, rng, protect=()):
    """round(rate * len) substitutions at distinct positions; positions in `protect` (donor/acceptor dinucleotides)
    are never touched. With protect=() it draws exactly what missing_copy_sim.mutate drew (random.sample over
    range(L) and over list(range(L)) pick the same positions)."""
    s = list(seq); cand = [i for i in range(len(s)) if i not in protect]; n = int(round(rate * len(s)))
    for i in rng.sample(cand, min(n, len(cand))):
        s[i] = rng.choice([b for b in 'ACGT' if b != s[i]])
    return ''.join(s)


# ================================================================ chromosome (was ideal_chromosome_sim.py)
# ⚠RefSeq spells it `lnc_RNA`, not `lncRNA` (919 records on chr16 alone), and 534 chr16 pseudogene
# exons carry `Parent=gene-...` with NO transcript record at all -- both are handled below, or 59% of
# chr16's genes silently produce no reads.
TX = {'mRNA', 'transcript', 'ncRNA', 'lnc_RNA', 'lncRNA', 'pseudogenic_transcript',
      'primary_transcript', 'tRNA', 'rRNA', 'snRNA', 'snoRNA', 'miRNA', 'misc_RNA',
      'V_gene_segment', 'C_gene_segment', 'J_gene_segment', 'ncRNA_gene'}


def load_transcripts(gff, chrom):
    """transcript_id -> (gene_name, strand, [exons]) for one chromosome."""
    gene_of, t2g, ex, strand = {}, {}, collections.defaultdict(list), {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom:
            continue
        if f[2] in ('gene', 'pseudogene', 'ncRNA_gene'):
            n = re.search(r'Name=([^;]+)', f[8]); i = re.search(r'ID=([^;]+)', f[8])
            if n and i:
                gene_of[i.group(1)] = n.group(1)
        elif f[2] in TX:
            i = re.search(r'ID=([^;]+)', f[8]); p = re.search(r'Parent=([^;,]+)', f[8])
            if i and p and p.group(1) in gene_of:
                t2g[i.group(1)] = gene_of[p.group(1)]; strand[i.group(1)] = f[6]
        elif f[2] == 'exon':
            p = re.search(r'Parent=([^;,]+)', f[8])
            if not p:
                continue
            par = p.group(1)
            if par in t2g:
                ex[par].append((int(f[3]), int(f[4])))
            elif par in gene_of:
                # exon parented straight to the gene (pseudogenes): the gene IS the transcript
                t2g.setdefault(par, gene_of[par]); strand.setdefault(par, f[6])
                ex[par].append((int(f[3]), int(f[4])))
    out = {}
    for t, e in ex.items():
        e.sort()
        if e:
            out[t] = (t2g[t], strand[t], e)
    return out


def cmd_chromosome(a):
    """Ideal-scenario chromosome simulation, per `docs/PREREG_ideal_chromosome_sim_2026-09-21.md` (md5 `ff226f41`).

    Establishes the CEILING of node construction under ideal input, with single-copy genes as the negative control.
    Arms:
      ideal  full-length reads, ends jittered, no readthrough      -> the ceiling
      trunc  + 5' degradation                                      -> cost of truncation
      rt     + readthrough molecules at the MEASURED 7.51%         -> cost of readthrough
    ⚠Jitter is mandatory: identical (chrom,pos,CIGAR) reads collapse under dedup (§6n0).
    ⚠The readthrough rate is measured, not invented: 33,058 of 439,985 real primary MAPQ-60 chr16 reads (7.51%) hit
     the exons of >= 2 distinct annotated genes at >= 25 bp. Writes FASTA reads `>tx|gene|fl or rt|k` to --out."""
    import pysam
    fa = pysam.FastaFile(a.genome)
    tx = load_transcripts(a.gff, a.chrom)
    print(f'{a.chrom}: {len(tx)} transcripts over '
          f'{len({v[0] for v in tx.values()})} genes', file=sys.stderr)

    # neighbouring same-strand transcript, for readthrough molecules
    order = sorted(tx.items(), key=lambda kv: kv[1][2][0][0])
    nxt = {}
    for i, (t, (g, s, e)) in enumerate(order):
        for j in range(i + 1, min(i + 25, len(order))):
            t2, (g2, s2, e2) = order[j]
            if g2 != g and s2 == s and e2[0][0] > e[-1][1]:
                nxt[t] = t2; break

    n_rt = 0
    seqs = {}
    with open(a.out, 'w') as fh:
        for idx, (t, (g, s, e)) in enumerate(order):
            if t not in seqs:
                seqs[t] = lib.spliced1(fa, a.chrom, e, s)
            base = seqs[t]
            if len(base) < 120:
                continue
            for k in range(a.reads_per_tx):
                rng = random.Random((a.seed * 7919) ^ (idx * 104729) ^ (k * 1299709))
                body = base
                tag = 'fl'
                if a.arm == 'rt' and nxt.get(t) and rng.random() < a.rt_rate:
                    t2 = nxt[t]
                    if t2 not in seqs:
                        seqs[t2] = lib.spliced1(fa, a.chrom, tx[t2][2], tx[t2][1])
                    if len(seqs[t2]) >= 120:
                        body = base + seqs[t2]; tag = 'rt'; n_rt += 1
                # MANDATORY end jitter -- identical reads collapse under dedup (§6n0)
                lo = rng.randint(0, a.jitter); hi = rng.randint(0, a.jitter)
                body = body[lo:len(body) - hi] if len(body) - hi > lo + 100 else body
                tf = 0.30 if a.arm == 'trunc' else 0.0
                for ri, (rd, _q) in enumerate(simulate_reads(body, 1, err=a.err, indel=a.err / 3,
                                                             seed=(idx * 131 + k), trunc_frac=tf)):
                    fh.write(f'>{t}|{g}|{tag}|{k}\n{rd}\n')
    print(f'arm={a.arm}  wrote {a.out}  (readthrough molecules: {n_rt})', file=sys.stderr)


# ================================================================ missing-copy (was missing_copy_sim.py)
def cmd_missing_copy(a):
    """Missing-copy simulations with truth (thesis objective O3; docs/PREREG_o3_reference_bias_2026-09-23.md arm A,
    docs/PREREG_o3_rna_only_2026-09-23.md control + addendum 2): reads from K annotated genes of CHROM plus an EXTRA
    COPY that is absent from the reference, mapped to the unmodified reference with the shipped minimap2 settings.

    modes
      transcript  each transcript of the gene mutated independently at divergence d (arm A: where do the copy's reads go?)
      genomic     the gene SPAN mutated once, every transcript read off it (one real genomic copy: the missing_copy_flag
                  positive control); also writes OUT.copies.fa/.mmi, the mutated spans, usable as a --confirm genome
      shuffled    the extra copy's transcript has exons 2 and 3 swapped (+ divergence d): how minimap2 represents an
                  exon-order rearrangement (insertion vs clip vs supplementary), the structural detector's control

    Read names carry the truth: `template|gene|tx|i` / `extra|gene|tx|i` (`shuffled|gene|i` in shuffled mode).
    ⚠ Seeds are stable_seed() since wave 7 (B2): reads differ from every run made with the old hash() seeds."""
    import pysam
    MODE, gtf, fa_p, CHROM, K, div, out, seed = a.mode, a.gtf, a.fasta, a.chrom, a.k, a.divergence, a.out, a.seed
    fa = pysam.FastaFile(fa_p); rng = random.Random(seed)

    ex = collections.defaultdict(list); gene = {}; st = {}
    for ln in open(gtf):
        if ln[0] == '#':
            continue
        r = ln.rstrip('\n').split('\t')
        if len(r) < 9 or r[0] != CHROM or r[2] != 'exon':
            continue
        tm = re.search(r'transcript_id "([^"]+)"', r[8]); gm = re.search(r'gene_id "([^"]+)"', r[8]) or re.search(r'gene_name "([^"]+)"', r[8])
        if not tm or not gm:
            continue
        ex[tm.group(1)].append((int(r[3]), int(r[4]))); gene[tm.group(1)] = gm.group(1); st[tm.group(1)] = r[6]
    min_exons = 5 if MODE == 'shuffled' else 3
    by_gene = collections.defaultdict(list)
    for t, e in ex.items():
        if len(e) >= min_exons:
            by_gene[gene[t]].append(t)
    genes = sorted(g for g, ts in by_gene.items() if (len(ts) == 1 if MODE == 'shuffled' else 1 <= len(ts) <= 6))
    chosen = rng.sample(genes, K)

    def tx_exons(t):
        e = sorted(ex[t]); exs = [fa.fetch(CHROM, x - 1, y).upper() for x, y in e]
        if st[t] == '-':
            exs = [lib.rc(x) for x in exs][::-1]
        return e, exs

    truth = {}; n = 0
    copies_fa = open(out + '.copies.fa', 'w') if MODE == 'genomic' else None
    with open(out + '.fq', 'w') as fq:
        for g in chosen:
            grng = random.Random(seed * 7919 + stable_seed(g) % 100000); ts = by_gene[g]
            span = None
            if MODE == 'genomic':
                gs = min(x for t in ts for x, y in ex[t]); ge = max(y for t in ts for x, y in ex[t])
                span = mutate(fa.fetch(CHROM, gs - 1, ge).upper(), div, grng)
                copies_fa.write(f'>{g}_copy {CHROM}:{gs}-{ge} d={div}\n{span}\n')
            for t in (ts[:1] if MODE == 'shuffled' else ts):
                e, exs = tx_exons(t); body = ''.join(exs)
                if MODE == 'transcript':
                    extra = mutate(body, div, grng)
                elif MODE == 'genomic':
                    cp = [span[x - gs:y - gs + 1] for x, y in e]
                    if st[t] == '-':
                        cp = [lib.rc(x) for x in cp][::-1]
                    extra = ''.join(cp)
                else:
                    sh = exs[:]; sh[1], sh[2] = sh[2], sh[1]; extra = mutate(''.join(sh), div, random.Random(seed + stable_seed(g) % 100000))
                kinds = (('template', body), ('shuffled' if MODE == 'shuffled' else 'extra', extra))
                for kind, seq in kinds:
                    for i, (rd, q) in enumerate(simulate_reads(seq, 10, err=0.001, indel=0.0003, seed=seed * 31 + (1 if kind != 'template' else 0) + stable_seed(t) % 100000)):
                        rd, q = jitter(rd, q, rng)
                        rid = f'{kind}|{g}|{t}|{i}' if MODE != 'shuffled' else f'{kind}|{g}|{i}'
                        fq.write(f'@{rid}\n{rd}\n+\n{q}\n'); truth[rid] = (kind, g, t, e[0][0], e[-1][1], len(e)); n += 1
    if copies_fa:
        copies_fa.close()
    with open(out + '.genes.txt', 'w') as fh:
        for g in chosen:
            fh.write(g + '\n')
    subprocess.run(f"{MM2} -t {a.threads} {fa_p} {out}.fq 2>/dev/null | samtools sort -@2 -o {out}.bam - && samtools index {out}.bam", shell=True, check=True)
    if MODE == 'genomic':
        subprocess.run(f"minimap2 -x splice:hq -t {a.threads} -d {out}.copies.mmi {out}.copies.fa 2>/dev/null", shell=True, check=True)
    print(f'[missing_copy_sim {MODE} d={div}] {len(chosen)} genes, {n} reads -> {out}.bam', flush=True)

    # --- classify the reads' placements
    recs = collections.defaultdict(list)
    for r in pysam.AlignmentFile(out + '.bam'):
        if r.is_supplementary and MODE != 'shuffled':
            continue
        recs[r.query_name].append(r)
    S = collections.defaultdict(collections.Counter); de = collections.defaultdict(list); nint = collections.defaultdict(list); ins = collections.defaultdict(list)
    for rid, (kind, g, t, gs, ge, ne) in truth.items():
        rs = recs.get(rid, []); prim = [r for r in rs if not r.is_secondary and not r.is_supplementary and not r.is_unmapped]
        if not prim:
            S[kind]['unmapped'] += 1; continue
        p = prim[0]; at = p.reference_start + 1 <= ge and p.reference_end >= gs
        S[kind]['absorbed at template locus' if at else 'primary elsewhere'] += 1
        secs = [r for r in rs if r.is_secondary]; pas = p.get_tag('AS') if p.has_tag('AS') else 0
        if any((r.get_tag('AS') if r.has_tag('AS') else 0) >= 0.98 * pas for r in secs) or p.mapping_quality == 0:
            S[kind]['AS-tied / MAPQ 0'] += 1
        if p.has_tag('de'):
            de[kind].append(p.get_tag('de'))
        if MODE == 'shuffled':
            S[kind]['with supplementary'] += any(r.is_supplementary for r in rs)
            nint[kind].append(sum(1 for op, l in p.cigartuples if op == 3)); ins[kind].append(max([l for op, l in p.cigartuples if op == 1] or [0]))
            S[kind]['insertion >= 50 bp'] += ins[kind][-1] >= 50
    for kind in sorted(S):
        T = sum(1 for v in truth.values() if v[0] == kind); c = S[kind]
        line = f"[{MODE} d={div:.3f}] {kind:9s} reads {T:,}: " + ', '.join(f"{k} {v} ({100*v/T:.1f}%)" for k, v in c.items())
        if de[kind]:
            line += f" | median de {statistics.median(de[kind]):.4f}"
        if MODE == 'shuffled' and nint[kind]:
            line += f" | introns in primary median {statistics.median(nint[kind])} | largest insertion median {statistics.median(ins[kind])}"
        print(line)


# ================================================================ tandem (was tandem_copy_sim.py)
POLISH = "--assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 --polish-retained-ratio 10"


def which_copy(copies, s, e):
    """copy whose span covers >= 50% of [s,e) (used for catalog copies and assembled transcripts)"""
    best = None
    for i, c in enumerate(copies):
        o = min(e, c['end']) - max(s, c['start'])
        if o >= 0.5 * (e - s) and (best is None or o > best[1]):
            best = (i, o)
    return best[0] if best else None


def which_exon(copies, s, e):
    """(copy, exon) of the planted exon that overlaps [s,e) most (>= 50% of the shorter), else None"""
    best = None
    for i, c in enumerate(copies):
        for x, (a_, b_) in ((1, c['e1']), (2, c['e2'])):
            o = min(e, b_) - max(s, a_)
            if o >= 0.5 * min(e - s, b_ - a_) and (best is None or o > best[2]):
                best = (i, x, o)
    return best[:2] if best else None


def tandem_classify(bam, copies):
    import pysam
    recs = collections.defaultdict(list)
    for r in pysam.AlignmentFile(bam):
        recs[r.query_name].append(r)
    rows = []
    for name, rs in recs.items():
        src = int(name.split('|')[0][4:])
        prim = [r for r in rs if not r.is_secondary and not r.is_supplementary and not r.is_unmapped]
        sup = [r for r in rs if r.is_supplementary]; sec = [r for r in rs if r.is_secondary]
        if not prim:
            rows.append((name, src, 'unmapped', -1, 0, len(sup) > 0)); continue
        p = prim[0]; groups = [[]]; rp = p.reference_start
        for op, l in p.cigartuples:
            if op in (0, 7, 8, 2):
                groups[-1].append((rp, rp + l)); rp += l
            elif op == 3:
                rp += l; groups.append([])
        groups = [g for g in groups if g]
        clip = sum(l for op, l in p.cigartuples if op == 4)
        tied = p.mapping_quality == 0 or any(r.get_tag('AS') >= 0.98 * p.get_tag('AS') for r in sec if r.has_tag('AS'))
        if len(groups) == 1:
            cls = 'unspliced'
        else:
            g1, g2 = groups[0], groups[-1]
            h1 = which_exon(copies, g1[0][0], g1[-1][1]); h2 = which_exon(copies, g2[0][0], g2[-1][1])
            if h1 is None or h2 is None:
                cls = 'outside'
            elif h1[0] == h2[0] == src:
                cls = 'same_copy'
            elif h1[0] == h2[0]:
                cls = 'other_copy'
            else:
                cls = f'cross:E1@copy{h1[0]},E2@copy{h2[0]}'
        if clip >= 50 and cls in ('same_copy', 'other_copy', 'unspliced'):
            cls = 'partial'
        # a supplementary whose position is upstream of the primary while the read continues = backward join
        for s_ in sup:
            if s_.reference_start < p.reference_start and cls != 'cross_forward':
                cls = 'cross_backward'
        rows.append((name, src, cls, p.mapping_quality, int(tied), len(sup) > 0))
    return rows


def _run(cmd, log):
    with open(log, 'w') as lg:
        return subprocess.run(cmd, shell=True, stdout=lg, stderr=subprocess.STDOUT).returncode


def tandem_pipeline(a, prefix, bam, fasta, contig_len, copies, truth_src, placed):
    out = {}
    B = a.bin
    rc_ = _run(f"{B}/copy_assign --assemble-only --assembly-junctions strict {POLISH} --bam {bam} --fasta {fasta} --region sim:1-{contig_len} --out {prefix}.asm", f"{prefix}.asm.log")
    chim = ntx = 0; per_copy = collections.Counter()
    if rc_ == 0 and os.path.exists(f"{prefix}.asm.gtf"):
        tx = collections.defaultdict(list)
        for ln in open(f"{prefix}.asm.gtf"):
            f = ln.split('\t')
            if len(f) > 8 and f[2] == 'exon':
                tx[re.search(r'transcript_id "([^"]+)"', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
        for t, ex in tx.items():
            ntx += 1; hits = [which_exon(copies, s_, e_) for s_, e_ in ex]; cs = {h[0] for h in hits if h}
            if len(cs) >= 2:
                chim += 1
            elif len(cs) == 1:
                per_copy[next(iter(cs))] += 1
    out.update(asm_rc=rc_, asm_tx=ntx, asm_chimeric=chim, asm_per_copy=','.join(f'{i}:{per_copy[i]}' for i in range(len(copies))))
    rc_ = _run(f"{B}/gw_family_catalog --bam {bam} --fasta {fasta} --threads {a.threads} --out {prefix}.cat", f"{prefix}.cat.log")
    ncop = nfam = 0
    if rc_ == 0 and os.path.exists(f"{prefix}.cat.copies.tsv"):
        ncop = sum(1 for _ in open(f"{prefix}.cat.copies.tsv")) - 1
        nfam = sum(1 for l in open(f"{prefix}.cat.families.tsv") if not l.startswith('family_id') and int(l.split('\t')[1]) >= 2)
    out.update(cat_rc=rc_, cat_copies=ncop, cat_multicopy_families=nfam)
    acc = {'correct': 0, 'wrong': 0, 'abstain': 0}
    if ncop >= 2:
        with open(f"{prefix}.regions.txt", 'w') as f:
            f.write(f"sim:1-{contig_len}\n")
        rc_ = _run(f"{B}/copy_assign --bam {bam} --fasta {fasta} --regions {prefix}.regions.txt --families {prefix}.cat.copies.tsv --copies-fa {prefix}.cat.copies.fa --out {prefix}.o2", f"{prefix}.o2.log")
        if rc_ == 0 and os.path.exists(f"{prefix}.o2.assignments.tsv"):
            cat = {}
            for l in open(f"{prefix}.cat.copies.tsv"):
                f = l.rstrip('\n').split('\t')
                if f[0] != 'family_id':
                    cat[(f[0], f[1])] = (int(f[4]), int(f[5]))
            by = collections.defaultdict(list)
            for l in open(f"{prefix}.o2.assignments.tsv"):
                f = l.rstrip('\n').split('\t')
                if f[0] == 'read_name':
                    hdr = f; continue
                by[f[0]].append(dict(zip(hdr, f)))
            for name, src in truth_src.items():
                if name not in by:  # not a tied read: the aligner placed it; score that placement
                    acc['unique_' + ('correct' if placed.get(name) == src else 'wrong')] = acc.get('unique_' + ('correct' if placed.get(name) == src else 'wrong'), 0) + 1; continue
                rows = [r for r in by.get(name, []) if r['status'] == 'assigned' and r['origin_rejected'] == '0']
                if not rows:
                    acc['abstain'] += 1; continue
                ok = any(which_copy(copies, *cat.get((r['family_id'], r['catalog_copy_idx']), (0, 0))) == src for r in rows)
                acc['correct' if ok else 'wrong'] += 1
        out['o2_rc'] = rc_
    out.update(o2_tied_correct=acc['correct'], o2_tied_wrong=acc['wrong'], o2_tied_abstain=acc['abstain'], o2_unique_correct=acc.get('unique_correct', 0), o2_unique_wrong=acc.get('unique_wrong', 0))
    return out


def cmd_tandem(a):
    """Tandem-copy simulator (docs/PREREG_tandem_copy_sim_2026-09-24.md): plant k copies of a real two-exon gene into a
    chr20 background at identity p and spacing D, simulate reads from every copy, map with the shipped minimap2
    settings, classify each read's exon placements, and optionally run the pipeline (assembler -> catalog ->
    assignment) on it.

    Per-read classes (primary alignment):
      same_copy       both exons inside the read's source copy
      other_copy      both exons inside one other copy
      cross:E1@copyI,E2@copyJ  exon 1 in one copy, exon 2 in another (one colinear chain across copies)
      cross_backward  exon 2 upstream of exon 1 — impossible as one alignment; counted from supplementary records
      unspliced       no intron in the primary alignment
      partial         >= 50 bp soft-clipped
      unmapped
    plus MAPQ, AS-tied (a secondary within 0.98 of the primary AS) and a supplementary flag.
    Outputs: PREFIX.<layout>.k<k>.p<p>.fa/.fq/.bam/.reads.tsv per condition and PREFIX.<layout>.k<k>.summary.tsv."""
    import pysam
    # ---- the gene: exon 1, an 800 bp intron with its own donor/acceptor ends, exon 2
    fa = pysam.FastaFile(a.fasta)
    exons = []
    for ln in open(a.gtf):
        f = ln.split('\t')
        if len(f) > 8 and f[0] == a.chrom and f[2] == 'exon' and f'transcript_id "{a.transcript}"' in f[8]:
            exons.append((int(f[3]), int(f[4]))); strand = f[6]
    exons = sorted(set(exons)); assert len(exons) == 2, f'{a.transcript}: need exactly 2 exons, got {exons}'
    (e1s, e1e), (e2s, e2e) = exons
    rc = lib.rc
    ex1 = fa.fetch(a.chrom, e1s - 1, e1e).upper(); ex2 = fa.fetch(a.chrom, e2s - 1, e2e).upper(); intr = fa.fetch(a.chrom, e1e, e2s - 1).upper()
    half = a.intron // 2; intr = intr[:half] + intr[-half:]
    if strand == '-':  # keep everything in the transcript's own orientation and plant it on the + strand of the contig
        ex1, ex2, intr = rc(ex2), rc(ex1), rc(intr)
    unit = ex1 + intr + ex2
    bchrom, brange = a.background.split(':'); bs, be = [int(x) for x in brange.split('-')]
    background = fa.fetch(bchrom, bs - 1, be).upper()

    L1, LI, L2 = len(ex1), len(intr), len(ex2)
    SPLICE = {L1, L1 + 1, L1 + LI - 2, L1 + LI - 1}  # GT..AG of the planted intron

    def build(p, seed):
        """contig with k copies; copies[i] = dict(e1, e2 (0-based half-open exon intervals), tx (spliced transcript))."""
        rng = random.Random(seed); offset = 50000; copies = []; seq = background
        units = [unit] + [mutate(unit, 1 - p, rng, SPLICE) for _ in range(a.copies - 1)]
        if a.layout == 'tandem':
            pos = offset
            for i, u in enumerate(units):
                seq = seq[:pos] + u + seq[pos + len(u):]
                copies.append({'e1': (pos, pos + L1), 'e2': (pos + L1 + LI, pos + len(u)), 'tx': u[:L1] + u[L1 + LI:]})
                pos += len(u) + a.distance
        else:
            # E1_0 +donor [d] E1_1 +donor [d] ... | intron body | [d] acceptor+ E2_0 [d] acceptor+ E2_1 ...
            # every planted exon keeps 50 bp of its own intron flank (donor after E1, acceptor before E2), so ANY chain
            # E1_i -> E2_j is canonical; spacers are background; the intron body is copy 0's
            d = a.distance; F = 50
            head, body_mid, tail = intr[:F], intr[F:-F], intr[-F:]
            e1_block = ''; e1_pos = []
            for i, u in enumerate(units):
                e1_pos.append(offset + len(e1_block)); e1_block += u[:L1] + head + background[offset + 1000 * (i + 1): offset + 1000 * (i + 1) + d]
            e2_block = ''; e2_pos = []
            for i, u in enumerate(units):
                e2_block += background[offset + 30000 + 1000 * (i + 1): offset + 30000 + 1000 * (i + 1) + d] + tail
                e2_pos.append(offset + len(e1_block) + len(body_mid) + len(e2_block)); e2_block += u[L1 + LI:]
            body = e1_block + body_mid + e2_block
            seq = seq[:offset] + body + seq[offset + len(body):]
            for i, u in enumerate(units):
                copies.append({'e1': (e1_pos[i], e1_pos[i] + L1), 'e2': (e2_pos[i], e2_pos[i] + L2), 'tx': u[:L1] + u[L1 + LI:]})
        for c in copies:
            c['start'], c['end'] = c['e1'][0], c['e2'][1]
        return seq, copies

    sweep = [float(x) for x in a.sweep.split(',')] if a.sweep else [a.identity]
    summary = []
    for p in sweep:
        prefix = f"{a.out}.{a.layout}.k{a.copies}.p{p}"
        seq, copies = build(p, a.seed)
        with open(prefix + '.fa', 'w') as f:
            f.write(f">sim\n{seq}\n")
        subprocess.run(f"samtools faidx {prefix}.fa", shell=True, check=True)
        truth_src = {}
        with open(prefix + '.fq', 'w') as fq:
            for i, c in enumerate(copies):
                rng = random.Random(a.seed * 101 + i)
                for n, (rd, q) in enumerate(simulate_reads(c['tx'], a.reads, err=0.001, indel=0.0003, seed=a.seed * 31 + i, trunc_frac=0.10)):
                    rd, q = jitter(rd, q, rng)
                    name = f"copy{i}|{n}"; fq.write(f"@{name}\n{rd}\n+\n{q}\n"); truth_src[name] = i
        subprocess.run(f"{MM2} -t {a.threads} {prefix}.fa {prefix}.fq 2>/dev/null | samtools sort -o {prefix}.bam - 2>/dev/null && samtools index {prefix}.bam", shell=True, check=True)
        rows = tandem_classify(prefix + '.bam', copies)
        with open(prefix + '.reads.tsv', 'w') as f:
            f.write('read\tsource_copy\tclass\tmapq\ttied\tsupplementary\n')
            for r in rows:
                f.write('\t'.join(str(x) for x in r) + '\n')
        n = len(rows); cls = collections.Counter(r[2] for r in rows)
        by_src = {i: collections.Counter(r[2] for r in rows if r[1] == i) for i in range(len(copies))}
        row = {'copies': a.copies, 'identity': p, 'distance': a.distance, 'reads': n, 'mapq0': sum(1 for r in rows if r[3] == 0), 'tied': sum(r[4] for r in rows), 'supplementary': sum(1 for r in rows if r[5])}
        for k in ('same_copy', 'other_copy', 'cross_backward', 'unspliced', 'partial', 'outside', 'unmapped'):
            row[k] = cls[k]
        crosses = {k: v for k, v in cls.items() if k.startswith('cross:')}
        row['cross_copy'] = sum(crosses.values()); row['cross_kinds'] = ';'.join(f'{k[6:]}={v}' for k, v in sorted(crosses.items()))
        row['cross_from_copy0'] = sum(v for k, v in by_src[0].items() if k.startswith('cross:')); row['cross_from_others'] = sum(v for i in range(1, len(copies)) for k, v in by_src[i].items() if k.startswith('cross:'))
        if a.pipeline:
            # aligner placement per read: the copy holding the primary's first exon (None when cross/unspliced/unmapped)
            recs = {}
            for r_ in pysam.AlignmentFile(prefix + '.bam'):
                if not (r_.is_secondary or r_.is_supplementary or r_.is_unmapped):
                    recs[r_.query_name] = which_copy(copies, r_.reference_start, r_.reference_end)
            row.update(tandem_pipeline(a, prefix, prefix + '.bam', prefix + '.fa', len(seq), copies, truth_src, recs))
        summary.append(row)
        print(f"[{a.layout} k={a.copies} p={p}] reads {n}: " + ', '.join(f"{k} {v}" for k, v in row.items() if k not in ('copies', 'identity', 'distance', 'reads') and v), flush=True)
    keys = list(summary[0].keys())
    with open(f"{a.out}.{a.layout}.k{a.copies}.summary.tsv", 'w') as f:
        f.write('\t'.join(keys) + '\n')
        for r in summary:
            f.write('\t'.join(str(r.get(k, '')) for k in keys) + '\n')
    print(f"wrote {a.out}.{a.layout}.k{a.copies}.summary.tsv")


# ================================================================ copies (was copy_assign_read_truth.py sim)
def cmd_copies(a):
    """O2 read-level truth simulation (docs/PREREG_o2_read_truth_2026-09-23.md): every copy of every multi-copy family
    (copies.tsv; >= 2 copies, spliced sequence >= 300 bp) gets min(100, max(10, n_reads)) HiFi-model reads (jittered,
    <= 10% trimmed) named `family|copy|i`, mapped with the shipped minimap2 settings against INDEX (genome-wide for the
    real experiment); also writes OUT.copies_used.tsv and per-copy closest-sibling identity (OUT.sibling.tsv).
    ⚠ copies.tsv is read POSITIONALLY (family_id, copy_idx, tid, chrom, start, end, n_exon, strand, n_reads).
    ⚠ Per-copy seeds are stable_seed(family, copy) since wave 7 (B2): reads differ from every earlier run."""
    tsv, fa, idx, out, seed = a.copies_tsv, a.copies_fa, a.index, a.out, a.seed
    rng = random.Random(seed)
    # copies.tsv: family_id copy_idx tid chrom start end n_exon strand n_reads
    rows = [l.rstrip('\n').split('\t') for l in open(tsv)][1:]
    nfam = collections.Counter(r[0] for r in rows)
    seqs = {}
    name = None
    for l in open(fa):
        if l.startswith('>'):
            h = l[1:].split('|'); name = (h[0], h[1]); seqs[name] = []
        else:
            seqs[name].append(l.strip())
    seqs = {k: ''.join(v).upper() for k, v in seqs.items()}
    n = 0; ncopies = 0
    with open(out + '.fq', 'w') as fq, open(out + '.copies_used.tsv', 'w') as cu:
        cu.write('family_id\tcopy_idx\tchrom\tstart\tend\tn_reads_real\tn_sim\tlen\n')
        for r in rows:
            fam, ci = r[0], r[1]
            if nfam[fam] < 2:
                continue
            s = seqs.get((fam, ci))
            if not s or len(s) < 300:
                continue
            k = min(100, max(10, int(r[8])))
            ncopies += 1
            cu.write(f'{fam}\t{ci}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[8]}\t{k}\t{len(s)}\n')
            for i, (rd, q) in enumerate(simulate_reads(s, k, err=0.001, indel=0.0003, seed=seed * 131 + stable_seed(fam, ci) % 1000003, trunc_frac=0.10)):
                rd, q = jitter(rd, q, rng)
                fq.write(f'@{fam}|{ci}|{i}\n{rd}\n+\n{q}\n'); n += 1
    print(f'[simO2] {ncopies} copies, {n} reads', flush=True)
    subprocess.run(f"{MM2} -t {a.threads} {idx} {out}.fq 2> {out}.mm2.log | samtools sort -@2 -o {out}.bam - && samtools index {out}.bam", shell=True, check=True)
    # closest sibling identity per copy: all-vs-all of the multi-copy families' sequences
    with open(out + '.copies.fa', 'w') as f:
        for (fam, ci), s in seqs.items():
            if nfam[fam] >= 2:
                f.write(f'>{fam}|{ci}\n{s}\n')
    subprocess.run(f"minimap2 -x asm20 -c -X -N 200 -p 0 -t {a.threads} {out}.copies.fa {out}.copies.fa > {out}.ava.paf 2>/dev/null", shell=True, check=True)
    best = {}
    for l in open(out + '.ava.paf'):
        f = l.split('\t'); x, y = f[0], f[5]
        if x == y or x.split('|')[0] != y.split('|')[0]:
            continue
        idn = int(f[9]) / int(f[10]); cov = (int(f[3]) - int(f[2])) / int(f[1])
        if cov < 0.5:
            continue
        for z in (x, y):
            if idn > best.get(z, (0, ''))[0]:
                best[z] = (idn, y if z == x else x)
    with open(out + '.sibling.tsv', 'w') as f:
        f.write('family_id\tcopy_idx\tclosest_identity\tclosest_copy\n')
        for (fam, ci) in seqs:
            if nfam[fam] >= 2:
                idn, nb = best.get(f'{fam}|{ci}', (float('nan'), 'NA'))
                f.write(f'{fam}\t{ci}\t{idn:.5f}\t{nb}\n')
    print(f'[simO2] wrote {out}.bam, {out}.sibling.tsv ({len(best)} copies with a sibling hit)')


# ================================================================ excise (was copy_assign_excision.py)
def cmd_excise(a):
    """PREREG adj/excise: remove copy X from a family, rerun copy_assign (genomic read-star default; skipped when
    OUT/A.assignments.tsv exists), follow X's MAPQ-60 reads, and look for the missing-copy signature: consistent
    mismatch sites among the origin-rejected reads that share a best candidate; reconstruct the copy from them.
    FAM_DIR holds copies.tsv, copies.fa, regions, forecast.tsv (and optionally A.assignments.tsv of the un-excised run)."""
    import pysam
    fam_dir, X, out = a.fam_dir, a.x, a.out
    BAM_P, FA_P, BIN = a.bam, a.fasta, a.bin
    B = pysam.AlignmentFile(BAM_P); FA = pysam.FastaFile(FA_P)
    os.makedirs(out, exist_ok=True)
    cp = list(csv.DictReader(open(f'{fam_dir}/copies.tsv'), delimiter='\t')); fid = cp[0]['family_id']
    keep = [r for r in cp if r['copy_idx'] != X]; remap = {r['copy_idx']: str(i) for i, r in enumerate(keep)}
    fa = {}; name = None
    for l in open(f'{fam_dir}/copies.fa'):
        if l.startswith('>'):
            name = l[1:].split('|')[1]; fa[name] = [l]
        else:
            fa[name].append(l)
    with open(f'{out}/copies.tsv', 'w') as o:
        w = csv.DictWriter(o, fieldnames=list(cp[0].keys()), delimiter='\t'); w.writeheader()
        for r in keep:
            r2 = dict(r); r2['copy_idx'] = remap[r['copy_idx']]; w.writerow(r2)
    with open(f'{out}/copies.fa', 'w') as o:
        for r in keep:
            h = fa[r['copy_idx']][0].split('|'); h[1] = remap[r['copy_idx']]; o.write('|'.join(h)); o.write(''.join(fa[r['copy_idx']][1:]))
    shutil.copy(f'{fam_dir}/regions', out)
    with open(f'{out}/forecast.tsv', 'w') as o:
        for l in open(f'{fam_dir}/forecast.tsv'):
            f = l.split('\t')
            if f[0] == 'copy_idx':
                o.write(l)
            elif f[0] in remap:
                o.write('\t'.join([remap[f[0]]] + f[1:]))
    if not os.path.exists(f'{out}/A.assignments.tsv'):
        subprocess.run(f'ulimit -v 10000000; {BIN} --bam {BAM_P} --fasta {FA_P} --families {out}/copies.tsv --copies-fa {out}/copies.fa --regions {out}/regions --out {out}/A > {out}/A.log 2> {out}/A.err', shell=True)
    A = {r['read_name']: r for r in csv.DictReader(open(f'{out}/A.assignments.tsv'), delimiter='\t')}
    A0 = {r['read_name']: r for r in csv.DictReader(open(f'{fam_dir}/A.assignments.tsv'), delimiter='\t')} if os.path.exists(f'{fam_dir}/A.assignments.tsv') else {}
    xr = [r for r in cp if r['copy_idx'] == X][0]; xc, xs, xe = xr['chrom'], int(xr['start']), int(xr['end'])
    xreads = {}
    for al in B.fetch(xc, xs, xe):
        if al.flag & 2308 or al.mapping_quality < 60:
            continue
        if any(b1 > xs and b0 < xe for b0, b1 in al.get_blocks()):
            xreads[al.query_name] = al
    fate = collections.Counter(); best_of_rej = collections.Counter()
    for n in xreads:
        r = A.get(n)
        if r is None:
            fate['O3 orphan (touches no remaining candidate; not in O2\'s output — the admission prototype\'s class)'] += 1; continue
        if r.get('n_candidates') == '0':
            fate['O2 orphan row (overlaps a copy, no locus aligns it)'] += 1
        elif r['status'] == 'assigned':
            fate['SILENT misassignment'] += 1
        elif r.get('origin_rejected') == '1':
            fate['origin-rejected'] += 1; best_of_rej[r['catalog_copy_idx']] += 1
        else:
            fate[r['status']] += 1
    n = len(xreads)
    print(f'== {fid} excise copy {X} ({xc}:{xs}-{xe}, nearest_ident {[l.split()[1] for l in open(fam_dir+"/forecast.tsv") if l.split()[0]==X][0]}): {n} MAPQ-60 reads')
    print('   fate:', {k: f'{v} ({100*v/n:.0f}%)' for k, v in fate.items()}, '| rejected reads\' best candidate (excised index):', dict(best_of_rej.most_common(3)))

    # ---- detector: consistent mismatch sites among origin-rejected reads sharing a best candidate (excised run) vs controls
    def consistent_sites(reads, y, tag):
        """align reads to Y's padded locus; return (n_reads, covered_kb, consistent_sites, top allele per site)"""
        yc, ys, ye = y['chrom'], int(y['start']), int(y['end']); pad = max((len(al.query_sequence) for al in reads), default=0)
        ls, le = max(0, ys - pad), ye + pad
        open(f'{out}/{tag}_Y.fa', 'w').write(f'>Y\n{FA.fetch(yc, ls, le)}\n'); open(f'{out}/{tag}_reads.fa', 'w').write(''.join(f'>{al.query_name}\n{al.query_sequence}\n' for al in reads))
        paf = subprocess.run(['minimap2', '-x', 'splice', '-c', '--eqx', '-N', '1', '-t', str(a.threads), f'{out}/{tag}_Y.fa', f'{out}/{tag}_reads.fa'], capture_output=True, text=True).stdout
        mism = collections.defaultdict(collections.Counter); cov = collections.Counter(); qseq = {al.query_name: al.query_sequence for al in reads}
        for l in paf.splitlines():
            f = l.split('\t'); cg = [t for t in f[12:] if t.startswith('cg:Z:')][0][5:]; t = int(f[7]); q = int(f[2]) if f[4] == '+' else int(f[1]) - int(f[3]); s = qseq[f[0]]
            if f[4] == '-':
                s = s.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]
            for num, op in re.findall(r'(\d+)([=XIDNS])', cg):
                num = int(num)
                if op == '=':
                    for k in range(num):
                        cov[t + k] += 1
                    t += num; q += num
                elif op == 'X':
                    for k in range(num):
                        cov[t + k] += 1; mism[t + k][s[q + k].upper()] += 1
                    t += num; q += num
                elif op in 'DN':
                    t += num
                elif op in 'IS':
                    q += num
        cons = {p: c.most_common(1)[0] for p, c in mism.items() if sum(c.values()) >= 3 and sum(c.values()) >= 0.5 * cov[p]}
        covered = sorted(p for p, c in cov.items() if c >= 3); covered_kb = len(covered) / 1000
        return len(reads), covered_kb, cons, (ls, covered)
    for yidx, cnt in best_of_rej.most_common(1):
        y = keep[int(yidx)]
        rej = [xreads[n] for n in xreads if A.get(n, {}).get('origin_rejected') == '1' and A[n]['catalog_copy_idx'] == yidx]
        nr, kb, cons, (ls, covered) = consistent_sites(rej, y, 'exc')
        dens = len(cons) / kb if kb else 0
        # control 1: Y's own MAPQ-60 reads vs Y ; control 2: rejected reads of the UN-excised run with best Y (if any)
        yc, ys, ye = y['chrom'], int(y['start']), int(y['end']); yreads = [al for al in B.fetch(yc, ys, ye) if not al.flag & 2308 and al.mapping_quality >= 60 and any(b1 > ys and b0 < ye for b0, b1 in al.get_blocks())][:200]
        nr1, kb1, cons1, _ = consistent_sites(yreads, y, 'ctl1') if yreads else (0, 0, {}, 0)
        orig_idx = [k for k, v in remap.items() if v == yidx][0]
        rej0 = [al for al in B.fetch(yc, max(0, ys - 200000), ye + 200000) if not al.flag & 2308 and A0.get(al.query_name, {}).get('origin_rejected') == '1' and A0[al.query_name]['catalog_copy_idx'] == orig_idx][:200]
        nr2, kb2, cons2, _ = consistent_sites(rej0, y, 'ctl2') if rej0 else (0, 0, {}, 0)
        print(f'   detector at Y=excised idx {yidx} (orig {orig_idx}): rejected reads {nr}, covered {kb:.1f} kb, consistent sites {len(cons)} = {dens:.1f}/kb | control Y-own reads: {nr1} reads, {len(cons1)/kb1 if kb1 else 0:.1f}/kb | control un-excised rejected@Y: {nr2} reads, {len(cons2)/kb2 if kb2 else 0:.1f}/kb')
        # recovery: patch Y with the consistent alleles, compare to X and to Y
        if cons and covered:
            # reconstruct ONLY the covered stretch (the largest run of read-covered positions), Y's bases with the
            # consistent sites patched by the majority read allele; score it against X's locus and Y's locus
            runs = []
            for p in covered:
                if runs and p - runs[-1][1] <= 50:
                    runs[-1][1] = p
                else:
                    runs.append([p, p])
            lo, hi = max(runs, key=lambda r: r[1] - r[0]); lo += ls; hi += ls + 1
            seq = list(FA.fetch(yc, lo, hi).upper())
            for p, (allele, _) in cons.items():
                q = p + ls - lo
                if 0 <= q < len(seq):
                    seq[q] = allele
            yorig = FA.fetch(yc, lo, hi).upper()
            open(f'{out}/reconstructed.fa', 'w').write('>reconstructed\n' + ''.join(seq) + '\n>Y_unpatched\n' + yorig + '\n')
            pad = max((len(al.query_sequence) for al in rej), default=0)
            open(f'{out}/targets.fa', 'w').write(f'>X_true\n{FA.fetch(xc, max(0, xs - pad), xe + pad)}\n>Y_locus\n{FA.fetch(yc, max(0, ys - pad), ye + pad)}\n')
            paf = subprocess.run(['minimap2', '-x', 'asm20', '-c', '-N', '5', '-p', '0', '-t', str(a.threads), f'{out}/targets.fa', f'{out}/reconstructed.fa'], capture_output=True, text=True).stdout
            best = {}
            for l in paf.splitlines():
                f = l.split('\t'); ident = int(f[9]) / max(1, int(f[10])); k = (f[0], f[5])
                if k not in best or int(f[10]) > best[k][1]:
                    best[k] = (ident, int(f[10]))
            g = lambda q, t: best.get((q, t), (0, 0))
            print(f'   recovery over the covered {hi-lo} bp: reconstructed -> X_true {g("reconstructed","X_true")[0]:.4f} | -> Y {g("reconstructed","Y_locus")[0]:.4f} ; unpatched Y -> X_true {g("Y_unpatched","X_true")[0]:.4f} | -> Y {g("Y_unpatched","Y_locus")[0]:.4f}')


# ================================================================ CLI
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    def add(name, func, help_):
        p = sub.add_parser(name, help=help_, description=func.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
        p.set_defaults(func=func)
        return p

    p = add('chromosome', cmd_chromosome, 'ideal-chromosome read simulation (was ideal_chromosome_sim.py)')
    for x in ('--gff', '--genome', '--chrom', '--out'):
        p.add_argument(x, required=True)
    p.add_argument('--arm', required=True, choices=['ideal', 'trunc', 'rt'])
    p.add_argument('--reads-per-tx', type=int, default=10)
    p.add_argument('--err', type=float, default=0.001)
    p.add_argument('--jitter', type=int, default=30, help='max bp trimmed from EACH end (mandatory)')
    p.add_argument('--rt-rate', type=float, default=0.0751, help='MEASURED on real chr16 reads')
    p.add_argument('--seed', type=int, default=17)

    p = add('missing-copy', cmd_missing_copy, 'O3 missing-copy simulations (was missing_copy_sim.py)')
    p.add_argument('mode', choices=('transcript', 'genomic', 'shuffled'))
    p.add_argument('gtf'); p.add_argument('fasta'); p.add_argument('chrom')
    p.add_argument('k', type=int); p.add_argument('divergence', type=float); p.add_argument('out'); p.add_argument('seed', type=int)
    p.add_argument('--threads', type=int, default=4, help='minimap2 threads (the old script hard-coded 4)')

    p = add('tandem', cmd_tandem, 'tandem / interleaved copy simulation (was tandem_copy_sim.py)')
    p.add_argument('--fasta', required=True); p.add_argument('--gtf', required=True); p.add_argument('--out', required=True)
    p.add_argument('--transcript', default='rna-NR_161305.1'); p.add_argument('--chrom', default='chr20')
    p.add_argument('--background', default='chr20:20000000-20200000')
    p.add_argument('--copies', type=int, default=2); p.add_argument('--identity', type=float, default=0.99); p.add_argument('--sweep', default='')
    p.add_argument('--layout', choices=['tandem', 'interleaved'], default='tandem', help='tandem: whole-gene copies spaced --distance apart; interleaved: exon-level tandem duplication (E1 E1p ... E2 E2p, spacer --distance), where the cross-copy chain E1p->E2 has the SHORTER intron')
    p.add_argument('--distance', type=int, default=8000); p.add_argument('--reads', type=int, default=50); p.add_argument('--intron', type=int, default=800)
    p.add_argument('--seed', type=int, default=1); p.add_argument('--pipeline', action='store_true'); p.add_argument('--bin', default=os.environ.get('RUSTLE_BIN', 'target/release'))
    p.add_argument('--threads', type=int, default=3)

    p = add('copies', cmd_copies, 'O2 read-truth simulation from a catalog (was copy_assign_read_truth.py sim)')
    p.add_argument('copies_tsv'); p.add_argument('copies_fa'); p.add_argument('index', help='minimap2 index (.mmi) or FASTA to map to')
    p.add_argument('out'); p.add_argument('seed', type=int)
    p.add_argument('--threads', type=int, default=4, help='minimap2 threads (the old script hard-coded 4)')

    p = add('excise', cmd_excise, 'excise one copy and look for the missing-copy signature (was copy_assign_excision.py)')
    p.add_argument('fam_dir'); p.add_argument('x', help='copy_idx to excise'); p.add_argument('out')
    p.add_argument('--bam', default='/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam')
    p.add_argument('--fasta', default='/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3_contigs.fa')
    p.add_argument('--bin', default='/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign')
    p.add_argument('--threads', type=int, default=4, help='minimap2 threads (the old script hard-coded 4)')

    a = ap.parse_args(argv)
    a.func(a)


if __name__ == '__main__':
    main()
