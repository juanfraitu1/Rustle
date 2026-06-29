#!/usr/bin/env bash
# Adversarial read-backing verification of a candidate VG win (the ZFY-real vs
# DAZ3-phantom discriminant). Reuses the locus.bam already extracted by
# winloci_eval.sh (G's window + sibling window), so NO extra genome pass.
#
# For the reads whose alignment AT GENE G supports G's annotated intron chain,
# compare each read's edit distance (NM) at G vs its best NM at any OTHER copy
# present in the bam. A copy is read-backed ("real") only if >= MIN_ANCHOR reads
# STRICTLY prefer it (NM_G < NM_other - T) or align only here; it is a phantom
# ("pure_spillover", the DAZ3 pattern) if no read prefers it and >=1 read's
# sequence favors a sibling.
#
# Usage: winloci_verify.sh GENE CHROM START END STRAND
# Env:   GI=<gene_introns.tsv>  T=2 (NM margin)  MIN_ANCHOR=2  K=0 (0=adaptive)
# Emits: one JSON object on stdout.
set -uo pipefail
GENE="${1:?gene}"; CHROM="${2:?chrom}"; START="${3:?start}"; END="${4:?end}"; STRAND="${5:?strand}"
ROOT=/mnt/c/Users/jfris/Desktop/Rustle
GI="${GI:-$ROOT/bench/paralog_secondary_scan/scan_out/gene_introns.tsv}"
BAM="/tmp/winloci/$GENE/locus.bam"
T="${T:-2}"; MIN_ANCHOR="${MIN_ANCHOR:-2}"; K="${K:-0}"
if [ ! -s "$BAM" ]; then echo "{\"gene\":\"$GENE\",\"verdict\":\"error\",\"reason\":\"no locus.bam (run winloci_eval first)\"}"; exit 0; fi

samtools view "$BAM" 2>/dev/null > "/tmp/winloci/$GENE/aln.sam" || true

python3 - "$GENE" "$CHROM" "$START" "$END" "$STRAND" "$GI" "/tmp/winloci/$GENE/aln.sam" "$T" "$MIN_ANCHOR" "$K" <<'PY'
import sys, re, collections, json
gene, chrom, start, end, strand, gi, sam, T, MIN_ANCHOR, K = sys.argv[1:11]
start, end, T, MIN_ANCHOR, K = int(start), int(end), float(T), int(MIN_ANCHOR), int(K)

# gene's annotated intron chain
introns = set()
for line in open(gi):
    if line.startswith('gene_id'): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 7 or f[0] != gene: continue
    if f[6]:
        for pr in f[6].split(','):
            d, a = pr.split(':'); introns.add((int(d), int(a)))
    break
n_int = len(introns)
Kg = K if K > 0 else max(2, min(8, (n_int + 2)//3))

def cig_parse(c): return re.findall(r'(\d+)([MIDNSH=X])', c)
def junctions(pos, c):
    cur = pos; js = []
    for ln, op in cig_parse(c):
        ln = int(ln)
        if op == 'N': js.append((cur-1, cur+ln)); cur += ln
        elif op in 'MD=X': cur += ln
    return js
def alen(c): return sum(int(l) for l, o in cig_parse(c) if o in 'M=XD') or 1

# group placements per read
reads = collections.defaultdict(list)   # qname -> [(chrom,pos,nm,alen,nmatch_G)]
for line in open(sam):
    f = line.rstrip('\n').split('\t')
    if len(f) < 11: continue
    flag = int(f[1])
    if flag & 0x800: continue   # supplementary
    rstrand = '-' if flag & 16 else '+'
    c, pos, cig = f[2], int(f[3]), f[5]
    if cig == '*': continue
    nm = next((int(t[5:]) for t in f[11:] if t.startswith('NM:i:')), None)
    if nm is None: continue
    al = alen(cig)
    at_G = (c == chrom and start-10 <= pos <= end and rstrand == strand)
    nmatch = sum(1 for j in junctions(pos, cig) if j in introns) if at_G else 0
    reads[f[0]].append((c, pos, nm, al, nmatch, at_G))

ga = sa = amb = pure = 0
sib = collections.Counter()
for qn, placs in reads.items():
    # best chain-supporting alignment AT G
    gbest = None
    for (c, pos, nm, al, nmatch, at_G) in placs:
        if at_G and nmatch >= Kg:
            if gbest is None or nm < gbest[0]: gbest = (nm, al)
    if gbest is None: continue
    nm_G, al_G = gbest
    # best alignment elsewhere (other locus), comparable extent
    other = None
    for (c, pos, nm, al, nmatch, at_G) in placs:
        if at_G: continue
        if al < 0.7 * al_G: continue
        if other is None or nm < other[0]: other = (nm, (c, pos))
    if other is None:
        pure += 1; ga += 1; continue
    dM = other[0] - nm_G       # positive => G fits better
    if dM >= T: ga += 1
    elif dM <= -T: sa += 1; sib[other[1][0]] += 1
    else: amb += 1

n_chain = ga + sa + amb
if ga >= MIN_ANCHOR: verdict = 'real'
elif ga == 0 and sa > 0: verdict = 'phantom'
elif amb >= max(sa, ga) and amb >= 3: verdict = 'ambiguous'
else: verdict = 'weak'
own_frac = round(ga / n_chain, 3) if n_chain else 0.0
print(json.dumps(dict(
    gene=gene, n_introns=n_int, K=Kg, n_chain_reads=n_chain,
    n_strict_owner=ga, n_pure_owner=pure, n_tied=amb, n_sibling_favored=sa,
    owner_frac=own_frac, dom_sibling_chrom=(sib.most_common(1)[0][0] if sib else None),
    verdict=verdict,
    reason=f"ga={ga}(>= {MIN_ANCHOR}?) sa={sa} amb={amb} pure={pure}")))
PY
