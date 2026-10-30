#!/usr/bin/env python3
"""SEMI-GUIDED FAMILIES (2026-10-30, user: "Lets focus on implementing the semi-guided mode, notice that in the artifact for CDY2B for
example, if we used the secondaries we would be able to find at least part of the locus"). De novo families (minimap2 all-vs-all + MCL
over assembled loci) plus the SILENT COPIES that the reads' secondary alignments reveal (`semisec`: a candidate = merged secondaries
of >= 2 seed reads with AS >= 0.8 x the read's primary AS, silent = fewer than 2 primaries over it). CDY2B on A119b: 24 primary reads,
none spliced, so no assembled locus; 709 secondaries from CDY1 / CDY1B reads, 166 of 171 traced at AS >= 0.8 x primary, 52 spliced.

Rule (SEMIFAM-2/-3 lessons, LEDGER_O1 r1461-r1465: the seed read's primary names the right family only when the copy is CONFIRMED by
aligning the seed locus itself; unconfirmed candidates are repeat hubs and land in the wrong family 2 times in 3):
  1. seed locus  = the de novo locus (<base>.loci.gff3) overlapping the candidate's best read's seed pile (`best_seed`), most overlap;
  2. template    = the seed locus's spliced exon sequence (its copies.tsv row when it is in a family, else its loci.gff3 exons);
  3. confirm     = minimap2 -x splice of the template onto the candidate span (+ --pad each side): template coverage >= --cov (0.5)
                   and identity >= --min-id (0.7), the SEMIFAM-2 lift rule; the copy's exons = the aligned target blocks;
  4. attach      = the confirmed copy joins the seed locus's family (source `semisec`, member_status `silent`, n_reads = candidate reads);
                   a seed locus in no family forms a new family `SEMI<n>` with its confirmed copies (--no-new-families keeps base only).
A candidate overlapping an existing copy of the run (same strand, checked on the aligned strand too) is skipped. The seed pile resolves to a
copy already in a family when one overlaps it (loci.gff3 holds overlapping records of one locus), else to the loci.gff3 record. Families are not re-clustered.
  semifam.py BASE.fam.copies.tsv BASE.fam.loci.gff3 SEMISEC_DIR GENOME_FA OUT [--cov 0.5] [--min-id 0.7] [--pad 2000] [--threads 8]
             [--no-new-families] [--all-candidates]
--all-candidates also attaches candidates that are NOT silent (>= 2 primaries: expressed copies the assembler dropped).
Writes OUT.fam.copies.tsv (base rows + attached copies), OUT.semi.tsv (one row per candidate: seed locus, family, gate numbers, verdict)."""
import sys, os, csv, collections, subprocess, tempfile, bisect
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
a = sys.argv[1:]
opt = lambda k, d, f=float: f(a[a.index(k) + 1]) if k in a else d
COV = opt('--cov', 0.5); MINID = opt('--min-id', 0.7); PAD = opt('--pad', 2000, int); T = opt('--threads', 8, int)
NEWFAM = '--no-new-families' not in a; ALL = '--all-candidates' in a
pos = [x for i, x in enumerate(a) if not x.startswith('--') and (i == 0 or a[i - 1] not in ('--cov', '--min-id', '--pad', '--threads'))]
COPIES, LOCI, SS, GENOME, OUT = pos[:5]
rows = list(csv.DictReader(open(COPIES), delimiter='\t'))
hdr = list(rows[0].keys())
# de novo loci (all, families or not): exons per locus from the GFF3
loc_ex = collections.defaultdict(list); loc_meta = {}
for line in open(LOCI):
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon': continue
    g = f[8].split('gene=')[1].split(';')[0]
    loc_ex[g].append((int(f[3]) - 1, int(f[4]))); loc_meta[g] = (f[0], f[6] if f[6] in '+-' else '+')
by_chrom = collections.defaultdict(list)
for g, ex in loc_ex.items():
    ex.sort(); by_chrom[loc_meta[g][0]].append((ex[0][0], ex[-1][1], g))
for c in by_chrom: by_chrom[c].sort()
starts = {c: [x[0] for x in v] for c, v in by_chrom.items()}
maxlen = {c: max(x[1] - x[0] for x in v) for c, v in by_chrom.items()}
def loci_over(chrom, s, e):
    v = by_chrom.get(chrom, []); i = bisect.bisect_left(starts.get(chrom, []), s - maxlen.get(chrom, 0)); out = []
    while i < len(v) and v[i][0] < e:
        if v[i][1] > s: out.append((min(e, v[i][1]) - max(s, v[i][0]), v[i][2]))
        i += 1
    return out
fam_of = {r['tid']: r['family_id'] for r in rows}
row_of = {r['tid']: r for r in rows}
copy_iv = collections.defaultdict(list)  # (chrom, strand) -> [(start, end)]
for r in rows: copy_iv[(r['chrom'], r['strand'])].append((int(r['start']), int(r['end'])))
for k in copy_iv: copy_iv[k].sort()
def overlaps_copy(chrom, strand, s, e):
    v = copy_iv.get((chrom, strand), []); i = bisect.bisect_left(v, (s - 2_000_000, 0))
    while i < len(v) and v[i][0] < e:
        if v[i][1] > s: return True
        i += 1
    return False
cands = []
for r in csv.DictReader(open(f'{SS}/cands.tsv'), delimiter='\t'):
    if r['silent'] != 'true' and not ALL: continue
    cands.append(r)
print(f"{len(rows)} base copies in {len(set(fam_of.values()))} families; {len(loc_ex)} de novo loci; {len(cands)} semisec candidates", file=sys.stderr)
# seed locus per candidate
plan = []  # (cand, seed_locus, template exons (chrom, strand, exons), reason)
for i, c in enumerate(cands):
    sc, iv = c['best_seed'].rsplit(':', 1); ss, se = map(int, iv.split('-'))
    hits = loci_over(sc, ss, se)
    s, e = int(c['start']), int(c['end'])
    if overlaps_copy(c['chrom'], c['strand'], s, e): plan.append((c, None, None, 'overlaps_copy')); continue
    if not hits: plan.append((c, None, None, 'seed_no_locus')); continue
    # the seed pile's locus: a copy already in a family first (overlapping loci records are one locus: STRAINY-A audit r1463),
    # else the loci.gff3 record with most overlap
    inc = [(ov, g) for ov, g in hits if g in row_of]
    if not inc:
        best_l = max(hits)[1]; ch, st = loc_meta[best_l]
        inc = [(min(se, ce) - max(ss, cs), t) for t, r0 in row_of.items() if r0['chrom'] == ch and r0['strand'] == st
               for cs, ce in [(int(r0['start']), int(r0['end']))] if ce > ss and cs < se]
    seed = max(inc)[1] if inc else max(hits)[1]
    if seed in row_of: tex = (row_of[seed]['chrom'], row_of[seed]['strand'], F.blocks(row_of[seed]['exons']))
    else: tex = (loc_meta[seed][0], loc_meta[seed][1], loc_ex[seed])
    plan.append((c, seed, tex, 'align'))
# sequences: templates (spliced) and candidate regions (padded genomic span)
rc = str.maketrans('ACGTacgt', 'TGCAtgca')
chrom_seq = {}
def seq_of(chrom):
    if chrom not in chrom_seq:
        chrom_seq[chrom] = ''.join(subprocess.run(['samtools', 'faidx', GENOME, chrom], capture_output=True, text=True).stdout.split('\n')[1:])
    return chrom_seq[chrom]
need = sorted({p[0]['chrom'] for p in plan if p[3] == 'align'} | {p[2][0] for p in plan if p[3] == 'align'})
for ch in need: seq_of(ch)
tmp = tempfile.mkdtemp(prefix='semifam_')
region = {}
with open(f'{tmp}/tpl.fa', 'w') as wt, open(f'{tmp}/reg.fa', 'w') as wr:
    for i, (c, seed, tex, why) in enumerate(plan):
        if why != 'align': continue
        ch, st, ex = tex; sq = seq_of(ch); t = ''.join(sq[s0:e0] for s0, e0 in ex).upper()
        if st == '-': t = t.translate(rc)[::-1]
        rs, re_ = max(0, int(c['start']) - PAD), min(len(seq_of(c['chrom'])), int(c['end']) + PAD)
        region[i] = (rs, re_)
        wt.write(f">c{i}\n{t}\n"); wr.write(f">c{i}\n{seq_of(c['chrom'])[rs:re_].upper()}\n")
paf = subprocess.run(['minimap2', '-c', '-x', 'splice', '--secondary=yes', '-N', '200', '-p', '0.05', '-t', str(T), f'{tmp}/reg.fa', f'{tmp}/tpl.fa'],
                     capture_output=True, text=True).stdout
best = {}  # cand index -> (cov, identity, strand, target blocks)
for line in paf.split('\n'):
    f = line.split('\t')
    if len(f) < 12 or f[0] != f[5]: continue
    i = int(f[0][1:]); ql = int(f[1]); cov = (int(f[3]) - int(f[2])) / ql; idn = int(f[9]) / int(f[10]) if int(f[10]) else 0
    cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), '')
    blocks = []; pos_ = int(f[7]); bs = pos_; n = ''
    for ch in cg:
        if ch.isdigit(): n += ch; continue
        k = int(n); n = ''
        if ch in 'MDX=': pos_ += k
        elif ch == 'N': blocks.append((bs, pos_)); pos_ += k; bs = pos_
    blocks.append((bs, pos_))
    if i not in best or cov * idn > best[i][0] * best[i][1]: best[i] = (cov, idn, f[4], blocks)
new_rows = []; newfam = {}; n_new = 0; verdicts = collections.Counter()
with open(OUT + '.semi.tsv', 'w') as w:
    w.write('chrom\tstrand\tstart\tend\treads\tsilent\tclass\tbest_read\tseed_locus\tfamily\ttpl_cov\tidentity\tn_exon\tverdict\n')
    for i, (c, seed, tex, why) in enumerate(plan):
        cov = idn = 0.0; fam = ''; nex = 0; v = why
        if why == 'align':
            if i not in best: v = 'no_alignment'
            else:
                cov, idn, pst, bl = best[i]
                if cov < COV or idn < MINID: v = 'below_gate'
                else:
                    v = 'attached'
                    rs = region[i][0]; ex = [(rs + s0, rs + e0) for s0, e0 in bl]
                    strand = pst  # the template is in transcript sense, the region is + genomic: the PAF strand IS the copy's strand
                    if overlaps_copy(c['chrom'], strand, ex[0][0], ex[-1][1]): v = 'overlaps_copy'
                    fam = fam_of.get(seed)
                    if v != 'attached': pass
                    elif fam is None:
                        if not NEWFAM: v = 'seed_singleton'
                        else:
                            if seed not in newfam:
                                n_new += 1; newfam[seed] = f"SEMI{n_new}"
                                r0 = {k: 'NA' for k in hdr}; ch, st, sex = tex
                                r0.update(family_id=newfam[seed], copy_idx='0', tid=seed, chrom=ch, start=str(sex[0][0]), end=str(sex[-1][1]), n_exon=str(len(sex)),
                                          strand=st, n_reads='0', exons=','.join(f"{s0}-{e0}" for s0, e0 in sex), source='locus_rep', gene_id=seed,
                                          member_status='ungated', locus_start=str(sex[0][0]), locus_end=str(sex[-1][1]))
                                new_rows.append(r0)
                            fam = newfam[seed]
                    if v == 'attached':
                        nex = len(ex); tid = f"SS_{c['chrom']}_{ex[0][0]}_{strand}"
                        r1 = {k: 'NA' for k in hdr}
                        r1.update(family_id=fam, copy_idx='0', tid=tid, chrom=c['chrom'], start=str(ex[0][0]), end=str(ex[-1][1]), n_exon=str(nex), strand=strand,
                                  n_reads=c['reads'], exons=','.join(f"{s0}-{e0}" for s0, e0 in ex), max_family_identity=f"{idn:.4f}", source='semisec',
                                  gene_id=seed, member_status='silent', locus_start=str(ex[0][0]), locus_end=str(ex[-1][1]))
                        new_rows.append(r1)
        verdicts[v] += 1
        w.write(f"{c['chrom']}\t{c['strand']}\t{c['start']}\t{c['end']}\t{c['reads']}\t{c['silent']}\t{c.get('class', '')}\t{c['best']}\t{seed or ''}\t{fam or ''}\t{cov:.2f}\t{idn:.3f}\t{nex}\t{v}\n")
# copy_idx within family
cnt = collections.Counter(r['family_id'] for r in rows)
for r in new_rows: r['copy_idx'] = str(cnt[r['family_id']]); cnt[r['family_id']] += 1
with open(OUT + '.fam.copies.tsv', 'w') as w:
    w.write('\t'.join(hdr) + '\n')
    for r in rows + new_rows: w.write('\t'.join(r[k] for k in hdr) + '\n')
print(f"verdicts: {dict(verdicts)}; new families {n_new}; copies {len(rows)} -> {len(rows) + len(new_rows)}", file=sys.stderr)
