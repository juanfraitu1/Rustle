#!/usr/bin/env python3
"""O2 read-level truth (docs/PREREG_o2_read_truth_2026-09-23.md), one script, two modes.

  sim    O2 read-level truth simulation (PREREG_o2_read_truth_2026-09-23): reads from every copy of every multi-copy
         usage: copy_assign_read_truth.py sim copies.tsv copies.fa INDEX.mmi OUT_PREFIX SEED
  score  Per-read scoring of copy_assign --families output (one row per read x family). Three readings of the table:
         usage: CATALOG_TSV=cat.copies.tsv copy_assign_read_truth.py score PREFIX O2PREFIX

sim: every copy of every multi-copy family gets HiFi-model reads (jittered, trimmed) named `family|copy|i`, mapped
genome-wide with the shipped minimap2 settings; also writes per-copy closest-sibling identity.
score: per read, correct / wrong / conflict / abstain / lost under three readings of copy_assign's per-family table
(OWN = the read's true family's row, PRIMARY = rows with primary_local=1, ANY = any assigned row), by MAPQ-0 stratum
and by the source copy's divergence bin.
"""
import sys, csv, collections, subprocess, os, math
import sys, random, subprocess, collections, os


def sim(argv):
    sys.path.insert(0, '/mnt/c/Users/jfris/Desktop/Rustle/bench'); from sim_reads import simulate_reads
    tsv, fa, idx, out, seed = argv[0], argv[1], argv[2], argv[3], int(argv[4])
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
            if nfam[fam] < 2: continue
            s = seqs.get((fam, ci))
            if not s or len(s) < 300: continue
            k = min(100, max(10, int(r[8])))
            ncopies += 1
            cu.write(f'{fam}\t{ci}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[8]}\t{k}\t{len(s)}\n')
            for i, (rd, q) in enumerate(simulate_reads(s, k, err=0.001, indel=0.0003, seed=seed * 131 + hash(fam + ci) % 1000003, trunc_frac=0.10)):
                j5, j3 = rng.randint(0, 30), rng.randint(0, 30)
                rd = rd[j5:len(rd) - j3]; q = q[j5:len(q) - j3]
                fq.write(f'@{fam}|{ci}|{i}\n{rd}\n+\n{q}\n'); n += 1
    print(f'[simO2] {ncopies} copies, {n} reads', flush=True)
    subprocess.run(f"minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 4 {idx} {out}.fq 2> {out}.mm2.log | samtools sort -@2 -o {out}.bam - && samtools index {out}.bam", shell=True, check=True)
    # closest sibling identity per copy: all-vs-all of the multi-copy families' sequences
    with open(out + '.copies.fa', 'w') as f:
        for (fam, ci), s in seqs.items():
            if nfam[fam] >= 2: f.write(f'>{fam}|{ci}\n{s}\n')
    subprocess.run(f"minimap2 -x asm20 -c -X -N 200 -p 0 -t 4 {out}.copies.fa {out}.copies.fa > {out}.ava.paf 2>/dev/null", shell=True, check=True)
    best = {}
    for l in open(out + '.ava.paf'):
        f = l.split('\t'); a, b = f[0], f[5]
        if a == b or a.split('|')[0] != b.split('|')[0]: continue
        idn = int(f[9]) / int(f[10]); cov = (int(f[3]) - int(f[2])) / int(f[1])
        if cov < 0.5: continue
        for x in (a, b):
            if idn > best.get(x, (0, ''))[0]: best[x] = (idn, b if x == a else a)
    with open(out + '.sibling.tsv', 'w') as f:
        f.write('family_id\tcopy_idx\tclosest_identity\tclosest_copy\n')
        for (fam, ci) in seqs:
            if nfam[fam] >= 2:
                idn, nb = best.get(f'{fam}|{ci}', (float('nan'), 'NA'))
                f.write(f'{fam}\t{ci}\t{idn:.5f}\t{nb}\n')
    print(f'[simO2] wrote {out}.bam, {out}.sibling.tsv ({len(best)} copies with a sibling hit)')



def score(argv):
    P, O = argv[0], argv[1]; CAT = os.environ['CATALOG_TSV']
    cat = {}; div = {}
    for r in csv.DictReader(open(CAT), delimiter='\t'):
        cat[(r['family_id'], r['copy_idx'])] = (r['chrom'], int(r['start']), int(r['end']))
        if r.get('max_family_identity') not in (None, '', 'NA'): div[(r['family_id'], r['copy_idx'])] = 1 - float(r['max_family_identity'])
    def dbin(d):
        if d is None: return 'NA'
        return '<0.5%' if d < 0.005 else '0.5-1%' if d < 0.01 else '1-2%' if d < 0.02 else '2-5%' if d < 0.05 else '>=5%'
    def same_locus(a, b):
        if a is None or b is None or a[0] != b[0]: return False
        o = min(a[2], b[2]) - max(a[1], b[1]); return o >= 0.5 * min(a[2] - a[1], b[2] - b[1])
    prim = {}
    for ln in subprocess.run(['samtools', 'view', '-F', '2308', P + '.bam'], capture_output=True, text=True).stdout.splitlines():
        f = ln.split('\t'); prim[f[0]] = int(f[4])
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(O + '.assignments.tsv'), delimiter='\t'): by[r['read_name']].append(r)
    def tr(n): return tuple(n.split('|')[:2])
    def judge(rows_assigned, t):
        loci = {(r['family_id'], r['catalog_copy_idx']) for r in rows_assigned}
        if not loci: return 'abstain'
        ok = [k == t or same_locus(cat.get(k), cat.get(t)) for k in loci]
        if len(loci) > 1 and not all(ok): return 'conflict' if any(ok) else 'wrong'
        return 'correct' if all(ok) else 'wrong'
    S = {v: collections.defaultdict(collections.Counter) for v in ('OWN', 'PRIMARY', 'ANY')}
    n_mapq0 = 0; own_status = collections.Counter(); nprim = collections.Counter()
    for name, mq in prim.items():
        if mq != 0: continue
        n_mapq0 += 1
        t = tr(name); b = dbin(div.get(t)); rows = by.get(name, [])
        asg = lambda rs: [r for r in rs if r['status'] == 'assigned' and r['origin_rejected'] == '0']
        own = [r for r in rows if r['family_id'] == t[0]]
        own_status[tuple(sorted(r['status'] for r in own)) or ('no_row',)] += 1
        o_own = 'lost' if not rows else ('no_own_row' if not own else judge(asg(own), t))
        pr = [r for r in rows if r['primary_local'] == '1']; nprim[len(pr)] += 1
        o_pr = 'lost' if not rows else ('no_primary_row' if not pr else judge(asg(pr), t))
        o_any = 'lost' if not rows else judge(asg(rows), t)
        for view, o in (('OWN', o_own), ('PRIMARY', o_pr), ('ANY', o_any)):
            for key in ('ALL', b): S[view][key][o] += 1
    print(f'MAPQ-0 reads {n_mapq0}; own-family row status combos: {own_status.most_common(6)}; primary_local rows per read: {dict(nprim)}')
    order = ['ALL', '<0.5%', '0.5-1%', '1-2%', '2-5%', '>=5%', 'NA']
    for view in ('OWN', 'PRIMARY', 'ANY'):
        print(f'== {view}')
        for k in order:
            c = S[view].get(k)
            if not c: continue
            n = sum(c.values()); a = c['correct'] + c['wrong'] + c['conflict']
            acc = c['correct'] / a if a else float('nan')
            print(f"  {k:8s} n={n:5d} correct {c['correct']:4d} wrong {c['wrong']:4d} conflict {c['conflict']:4d} abstain {c['abstain']:4d} other {n - a - c['abstain']:4d} | acc_assigned {acc:.3f} coverage {a / n:.3f}")



if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in ('sim', 'score'):
        sys.exit(__doc__)
    {'sim': sim, 'score': score}[sys.argv[1]](sys.argv[2:])
