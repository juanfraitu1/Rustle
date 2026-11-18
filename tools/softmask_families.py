#!/usr/bin/env python3
"""Families whose copies are MOSTLY soft-masked sequence: one family id per line (register 1690, opt-in).
usage: softmask_families.py COPIES_TSV GENOME_FASTA FRAC > softmask.txt
A copy counts when at least FRAC of its exon bases are lowercase in the FASTA (needs FASTA.fai); a family is listed when
strictly more than half of its copies count (the majority rule of annotated_families.py --majority-in).
Satellite and transposon-class families are lowercase in any RepeatMasker / WindowMasker genome; rRNA often is not."""
import collections, csv, sys

def die(msg):
    print(f"[softmask_families] {msg}", file=sys.stderr); sys.exit(2)

def listed(copies, fasta, frac):
    fai = {}
    try:
        for line in open(fasta + '.fai'):
            f = line.split('\t'); fai[f[0]] = tuple(map(int, f[1:5]))
    except OSError as e:
        die(f"{e} (run samtools faidx)")
    n = collections.Counter(); hit = collections.Counter()
    with open(fasta, 'rb') as fa:
        def fetch(chrom, a, b):  # 0-based half-open
            length, off, lb, lw = fai[chrom]
            if not 0 <= a < b <= length: die(f"{chrom}:{a}-{b} is outside the contig ({length} bp)")
            fa.seek(off + a // lb * lw + a % lb)
            return fa.read((b - a) + (b - a) // lb + 2).replace(b'\n', b'').replace(b'\r', b'')[:b - a]
        for r in csv.DictReader(open(copies), delimiter='\t'):
            if r['chrom'] not in fai: die(f"contig {r['chrom']} of {copies} is not in {fasta}.fai")
            low = tot = 0
            for e in r['exons'].split(','):
                a, b = e.split('-'); s = fetch(r['chrom'], int(a), int(b))  # the exons column is 0-based half-open (mcl_families unit_seq)
                tot += len(s); low += sum(1 for c in s if 97 <= c <= 122)
            n[r['family_id']] += 1; hit[r['family_id']] += tot > 0 and low >= frac * tot
    return sorted(f for f in n if hit[f] * 2 > n[f]), len(n)

if __name__ == '__main__':
    if len(sys.argv) != 4: die(__doc__)
    try:
        frac = float(sys.argv[3])
    except ValueError:
        frac = -1
    if not 0 < frac <= 1: die(f"FRAC must be a number in (0, 1], got '{sys.argv[3]}' (0.9 is the tested value)")
    out, total = listed(sys.argv[1], sys.argv[2], frac)
    print('\n'.join(out))
    print(f"[softmask_families] {total} families, {len(out)} listed (most copies >= {frac:g} lowercase exon bases)", file=sys.stderr)
