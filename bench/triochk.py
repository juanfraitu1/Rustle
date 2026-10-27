"""TRIOCHK (docs/LEDGER_O3.md): trio WGS dosage of haplotype-specific 21-mers of two maternal-only copies.

  triochk.py query OUT_DIR MAT_FA PAT_FA   -> OUT_DIR/qc.u64 (sorted canonical 21-mers, ks's encoding) + OUT_DIR/labels.json
  triochk.py score OUT_DIR                 -> reads OUT_DIR/{mat,pat,Jim,Trib,Dolly}.i32 (ks outputs), prints the table
"""
import json
import subprocess
import sys

import numpy as np

K = 21
CODE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
REGIONS = {  # name: (haplotype, contig, start, end), 1-based inclusive as samtools faidx
    'N': ('mat', 'CM054600.2', 6_720_000, 6_820_000),
    'S': ('mat', 'CM054594.2', 95_909_997, 96_031_621),
    'P': ('pat', 'CM054570.2', 91_200_203, 91_204_230),
    'C': ('mat', 'CM054594.2', 95_700_000, 95_900_000),
}
ANIMALS = ['Jim', 'Trib', 'Dolly']
TRANSMIT = {'N': 'Dolly', 'S': 'Dolly', 'P': 'Trib'}


def kmers(seq):
    """Canonical 21-mers as ks computes them (A0 C1 G2 T3, min of forward and reverse complement); N resets."""
    mask = (1 << (2 * K)) - 1
    out, fw, rv, n = [], 0, 0, 0
    for ch in seq.upper():
        c = CODE.get(ch)
        if c is None:
            fw = rv = n = 0
            continue
        fw = ((fw << 2) | c) & mask
        rv = (rv >> 2) | ((3 - c) << (2 * (K - 1)))
        n += 1
        if n >= K:
            out.append(min(fw, rv))
    return out


def fetch(fa, contig, s, e):
    r = subprocess.run(['samtools', 'faidx', fa, f'{contig}:{s}-{e}'], capture_output=True, text=True, check=True)
    return ''.join(r.stdout.split('\n')[1:])


def query(out, mat, pat):
    fa = {'mat': mat, 'pat': pat}
    sets = {name: sorted(set(kmers(fetch(fa[h], c, s, e)))) for name, (h, c, s, e) in REGIONS.items()}
    allk = np.array(sorted(set().union(*map(set, sets.values()))), dtype=np.uint64)
    allk.tofile(f'{out}/qc.u64')
    json.dump({n: [int(x) for x in v] for n, v in sets.items()}, open(f'{out}/labels.json', 'w'))
    print({n: len(v) for n, v in sets.items()}, 'total', len(allk))


def score(out):
    qc = np.fromfile(f'{out}/qc.u64', dtype=np.uint64)
    cnt = {s: np.fromfile(f'{out}/{s}.i32', dtype=np.int32) for s in ['mat', 'pat'] + ANIMALS}
    sets = json.load(open(f'{out}/labels.json'))
    idx = {n: np.searchsorted(qc, np.array(v, dtype=np.uint64)) for n, v in sets.items()}
    own = {n: REGIONS[n][0] for n in sets}
    keep = {}
    for n, i in idx.items():
        if n == 'C':
            keep[n] = i[(cnt['mat'][i] == 1) & (cnt['pat'][i] == 1)]
        else:
            other = 'pat' if own[n] == 'mat' else 'mat'
            keep[n] = i[(cnt[own[n]][i] == 1) & (cnt[other][i] == 0)]
    lam = {a: float(np.median(cnt[a][keep['C']])) for a in ANIMALS}
    print('k-mers kept:', {n: len(v) for n, v in keep.items()}, 'of', {n: len(v) for n, v in idx.items()})
    print('lambda (2 copies):', lam)
    res = {'lambda': lam, 'kept': {n: len(v) for n, v in keep.items()}, 'targets': {}}
    for n in ['N', 'S', 'P']:
        row = {}
        for a in ANIMALS:
            c = cnt[a][keep[n]]
            row[a] = {'dosage': float(np.median(c)) / (lam[a] / 2), 'present': float((c >= 3).mean())}
        j, t = row['Jim']['dosage'], row[TRANSMIT[n]]['dosage']
        verdict = 'CONFIRMED' if 0.5 <= j <= 1.5 and t >= 0.5 else 'ASSEMBLY ERROR' if j < 0.25 else 'UNRESOLVED'
        other = [a for a in ('Trib', 'Dolly') if a != TRANSMIT[n]][0]
        od = row[other]['dosage']
        parents = 'parents differ' if od < 0.25 else 'both parents carry it' if od >= 0.5 else 'unclear'
        res['targets'][n] = {'animals': row, 'verdict': verdict, 'non_transmitting': parents}
        print(n, verdict, '|', parents, '|', ' '.join(f"{a} {row[a]['dosage']:.2f} ({row[a]['present']:.2f})" for a in ANIMALS))
    json.dump(res, open(f'{out}/result.json', 'w'), indent=1)


if __name__ == '__main__':
    if sys.argv[1] == 'query':
        query(*sys.argv[2:5])
    else:
        score(sys.argv[2])
