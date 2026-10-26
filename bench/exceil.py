"""EXCEIL (LEDGER_O3 prereg b7894f54): the most RNA could see of exon gains. Run with miniforge python (edlib).
usage: exceil.py EXHAP_DIR SDCN_DIR"""
import collections, random, sys
import edlib, pysam
E, S = sys.argv[1], sys.argv[2]
BAM = '/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam'
FA = {'ref': '/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta', 'mat': '/mnt/linuxdisk/tmp/p18/mat.fa',
      'pat': '/mnt/linuxdisk/tmp/p18/pat.fa'}
SEXREF = {'NC_073247.2', 'NC_073248.2', 'NC_011120.1'}
COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]
def sex(g, t): return t in SEXREF if g == 'ref' else t.split('_')[0] in ('chrX', 'chrY', 'chrM')
H = {g: collections.defaultdict(list) for g in FA}
for g in FA:
    for l in open(f'{S}/ex_sr_{g}.paf'):
        f = l.split('\t')
        if (int(f[3]) - int(f[2])) / int(f[1]) >= 0.8 and int(f[9]) / int(f[10]) >= 0.90 and not sex(g, f[5]):
            H[g][f[0]].append((f[5], int(f[7]), int(f[8]), f[4]))
def merged(xs):
    out = []
    for t, s, e, st in sorted(xs):
        if out and out[-1][0] == t and s < out[-1][2]: out[-1][2] = max(out[-1][2], e)
        else: out.append([t, s, e, st])
    return out
# the same unit keys as bench/exhap_units.py; keep one exon per unit
uid = {}
units = collections.defaultdict(dict)
for u in H['ref']:
    fam = u.split('|')[1]; loci = merged(H['ref'][u])
    key = tuple((t, s // 50, e // 50) for t, s, e, _ in loci)
    units[fam].setdefault(key, u)
for fam in units:
    for k, key in enumerate(sorted(units[fam])):
        uid[f'{fam}x{k}'] = units[fam][key]
gain = {}; split = {}
for l in open(f'{E}/units.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0] != 'unit': gain[f[0]] = f[6] == '1'; split[f[0]] = f[2]
nets = collections.defaultdict(list)
for p in ('net_dev.tsv', 'net_heldout_0.tsv', 'net_heldout_1.tsv', 'net_heldout_2.tsv'):
    for l in open(f'{E}/{p}'):
        f = l.split('\t')
        if f[0] != 'family': nets[f[0]].append(f[1])
fa = {g: pysam.FastaFile(p) for g, p in FA.items()}
def seqs(g, ex):
    return [(rc if st == '-' else str)(fa[g].fetch(t, s, e).upper()) for t, s, e, st in merged(H[g][ex])]
bam = pysam.AlignmentFile(BAM)
loci_of = collections.defaultdict(list)
for l in open(f'{E}/copies.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0] != 'family_id': loci_of[f[0]].append((f[3], int(f[4]), int(f[5]), f[7]))
def level(u):
    ex = uid[u]; R = seqs('ref', ex); rs = set(R)
    hap_only = set()
    for g in ('mat', 'pat'):
        X = seqs(g, ex)
        if len(X) > len(R) or not gain[u]:
            hap_only |= {x for x in X if x not in rs}
    if not hap_only: return 'none', 0
    want = set(nets[u]); hits = 0
    for c, s, e, st in loci_of[u]:
        for r in bam.fetch(c, s, e):
            if r.query_name not in want or r.is_secondary or r.is_supplementary: continue
            want.discard(r.query_name)
            qp = [q for q, p in r.get_aligned_pairs(matches_only=True) if s - 20 <= p < e + 20]
            if not qp: continue
            seg = r.query_sequence[min(qp):max(qp) + 1]
            seg = rc(seg) if st == '-' else seg
            d = lambda x: edlib.align(x, seg, mode='HW')['editDistance']
            if min(d(x) for x in hap_only) < min(d(x) for x in rs): hits += 1
    return ('B' if hits >= 2 else 'A'), hits
G = sorted(u for u in uid if gain.get(u))
rng = random.Random(20261025); C = sorted(rng.sample(sorted(u for u in uid if u in gain and not gain[u]), 300))
for name, us in (('GAIN units', G), ('control (non-gaining)', C)):
    c = collections.Counter(); bys = collections.defaultdict(collections.Counter)
    for u in us:
        lv, _ = level(u); c[lv] += 1; bys[split[u]][lv] += 1
    n = len(us)
    print(f"{name} n={n}: level A (sequence) {c['A'] + c['B']} ({100 * (c['A'] + c['B']) / n:.1f}%), level B (reads) {c['B']} "
          f"({100 * c['B'] / n:.1f}%), none {c['none']}; by split {dict((k, dict(v)) for k, v in bys.items())}")
