#!/usr/bin/env python3
"""PREREG ISOCON-ONLY TRANSCRIPTS (register 1515, docs/PREREG_isocon_only_2026-11-03.md).
  isocon_only.py set     SQ_DIR REF.gtf GFFC_DIR OUT.tsv        -> the Only set: IsoCon-FSM reference transcripts with no ours FSM
  isocon_only.py copies  TARGETS.tsv REF.gtf GENOME.fa OUT      -> OUT.copies.tsv + OUT.copies.fa: one single-copy family per gene
  isocon_only.py gtf     ALL_CLUSTERS.fa GENOME.mmi OUT.gtf      -> primary splice:hq placements as a GTF (the r1482 IsoCon mapping)
  isocon_only.py recover ONLY.tsv ARM_classification.txt OUT.tsv -> per Only transcript: recovered by an arm FSM or not
SQ_DIR holds sq3/{ours,isocon}/{arm}_classification.txt (unfiltered); GFFC_DIR the gffcompare *.tracking of the end-polish variants."""
import sys, os, re, csv, collections, subprocess
cmd = sys.argv[1]
FSM = 'full-splice_match'
def sh(c): return subprocess.run(c, shell=True, capture_output=True, text=True).stdout
def classification(path):
    return list(csv.DictReader(open(path), delimiter='\t'))
def ref_tx(gtf):
    """transcript id -> (gene_id, gene_name, n_exon); gene_id -> sorted exon list"""
    tx, gx = {}, collections.defaultdict(set)
    for line in open(gtf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        a = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
        t = tx.setdefault(a['transcript_id'], [a.get('gene_id', ''), a.get('gene_name', ''), 0])
        t[2] += 1
        gx[a.get('gene_id', '')].add((int(f[3]), int(f[4])))
    return {k: tuple(v) for k, v in tx.items()}, {g: sorted(e) for g, e in gx.items()}
if cmd == 'set':
    sq, ref, gffc, out = sys.argv[2:6]
    tx, _ = ref_tx(ref)
    cls = {a: classification(f"{sq}/sq3/{a}/{a}_classification.txt") for a in ('ours', 'isocon')}
    fsm = {a: {r['associated_transcript'] for r in cls[a] if r['structural_category'] == FSM} for a in cls}
    only = sorted(fsm['isocon'] - fsm['ours'])
    fl = collections.defaultdict(int)
    for r in cls['isocon']:
        if r['structural_category'] == FSM:
            try: fl[r['associated_transcript']] = max(fl[r['associated_transcript']], int(float(r['FL'])))
            except ValueError: pass
    # ours' best category on the Only transcript's gene
    rank = {'ISM_this': 0, 'ISM_other': 1, 'NIC': 2, 'NNC': 3, 'FSM_other': 4}
    by_gene = collections.defaultdict(list)
    for r in cls['ours']: by_gene[r['associated_gene']].append(r)
    def best_ours(t, g):
        cats = set()
        for r in by_gene.get(g, []):
            c = r['structural_category']
            if c == 'incomplete-splice_match': cats.add('ISM_this' if r['associated_transcript'] == t else 'ISM_other')
            elif c == 'novel_in_catalog': cats.add('NIC')
            elif c == 'novel_not_in_catalog': cats.add('NNC')
            elif c == FSM: cats.add('FSM_other')
        return min(cats, key=rank.get) if cats else 'none'
    variants = ['ours_tssr', 'ours_med', 'ours_ext', 'ours_extmed', 'ours_all']
    eq = {v: set() for v in variants}
    for v in variants:
        p = f"{gffc}/{v}.tracking"
        if not os.path.exists(p): continue
        for line in open(p):
            f = line.split('\t')
            if len(f) > 3 and f[3] == '=': eq[v].add(f[2].split('|')[-1])
    with open(out, 'w') as w:
        w.write("transcript\tgene\tname\tn_exon\tstratum\tisocon_fl\tours_best\t" + "\t".join(variants) + "\n")
        for t in only:
            g, name, n = tx.get(t, ('?', '?', 0))
            w.write(f"{t}\t{g}\t{name}\t{n}\t{'multi' if n > 1 else 'mono'}\t{fl[t]}\t{best_ours(t, g)}\t"
                    + "\t".join('1' if t in eq[v] else '0' for v in variants) + "\n")
    rows = [l.split('\t') for l in open(out).read().splitlines()[1:]]
    multi = [r for r in rows if r[4] == 'multi']
    print(f"S_ours {len(fsm['ours'])} S_isocon {len(fsm['isocon'])} Only {len(only)} (multi-exon {len(multi)}, mono {len(rows) - len(multi)}); "
          f"ours_only {len(fsm['ours'] - fsm['isocon'])}")
    print("ours_best on Only (multi):", dict(collections.Counter(r[6] for r in multi)))
    for i, v in enumerate(variants):
        if eq[v]: print(f"drift {v}: matches {sum(int(r[7 + i]) for r in multi)} of {len(multi)} multi-exon Only")
elif cmd == 'copies':
    targets, ref, genome, out = sys.argv[2:6]
    import pysam
    fa = pysam.FastaFile(genome)
    _, gx = ref_tx(ref)
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    n = 0
    with open(f"{out}.copies.tsv", 'w') as wt, open(f"{out}.copies.fa", 'w') as wf:
        wt.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tmax_family_identity\tlocus_start\tlocus_end\n")
        for r in csv.DictReader(open(targets), delimiter='\t'):
            g, chrom, strand = r['gene'], r['chrom'], r['strand']
            ex = []  # union of the gene's CAT exons (1-based inclusive)
            for a, b in gx.get(g, []):
                if ex and a <= ex[-1][1] + 1: ex[-1][1] = max(ex[-1][1], b)
                else: ex.append([a, b])
            if not ex: print(f"{g}\tno exons in {ref}", file=sys.stderr); continue
            seq = ''.join(fa.fetch(chrom, a - 1, b).upper() for a, b in ex)
            if strand == '-': seq = seq.translate(comp)[::-1]
            s, e = ex[0][0], ex[-1][1]
            wt.write(f"{g}\t0\t{g}\t{chrom}\t{s}\t{e}\t{len(ex)}\t{strand}\t{r['records']}\t{','.join(f'{a}-{b}' for a, b in ex)}\t0\t{s}\t{e}\n")
            wf.write(f">{g}|0|{chrom}:{s}-{e}|{strand}|nexon={len(ex)}\n{seq}\n"); n += 1
    print(f"{n} single-copy families -> {out}.copies.{{tsv,fa}}")
elif cmd == 'gtf':
    fa, mmi, out = sys.argv[2:5]
    sam = sh(f"minimap2 -ax splice:hq -uf --eqx -N 50 -p 0.1 --secondary=yes -t 16 {mmi} {fa} 2>/dev/null")
    n = 0
    with open(out, 'w') as w:
        for line in sam.splitlines():
            if line.startswith('@'): continue
            f = line.split('\t', 11)
            if int(f[1]) & 2308 or f[2] == '*': continue
            pos = int(f[3]); exons = []; cur = pos
            for ln, op in re.findall(r'(\d+)([MIDNSHP=X])', f[5]):
                ln = int(ln)
                if op in 'M=XD':
                    if not exons or exons[-1][1] != cur: exons.append([cur, cur + ln])
                    else: exons[-1][1] = cur + ln
                    cur += ln
                elif op == 'N': cur += ln
            exons = [(a, b - 1) for a, b in exons]
            strand = '-' if int(f[1]) & 16 else '+'
            fam = f[0].split(':')[0]; tid = re.sub(r'[;,()]', '_', f[0]); n += 1
            w.write(f"{f[2]}\tcand\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\tgene_id \"{fam}\"; transcript_id \"{tid}\";\n")
            for i, (a, b) in enumerate(exons): w.write(f"{f[2]}\tcand\texon\t{a}\t{b}\t.\t{strand}\t.\tgene_id \"{fam}\"; transcript_id \"{tid}\"; exon_number \"{i + 1}\";\n")
    print(f"{n} placed consensuses -> {out}")
elif cmd == 'recover':
    only, cls, out = sys.argv[2:5]
    fsm = {r['associated_transcript'] for r in classification(cls) if r['structural_category'] == FSM}
    rows = list(csv.DictReader(open(only), delimiter='\t'))
    with open(out, 'w') as w:
        w.write("transcript\tgene\tname\tstratum\tours_best\trecovered\n")
        for r in rows: w.write(f"{r['transcript']}\t{r['gene']}\t{r['name']}\t{r['stratum']}\t{r['ours_best']}\t{int(r['transcript'] in fsm)}\n")
    for st in ('multi', 'mono'):
        rs = [r for r in rows if r['stratum'] == st]
        k = sum(r['transcript'] in fsm for r in rs)
        print(f"{st}: recovered {k} of {len(rs)} ({100 * k / max(1, len(rs)):.1f}%)")
else:
    sys.exit(__doc__)
