#!/usr/bin/env python3
"""PREREG r1829 (docs/PREREG_near_tie_census_2026-11-19.md), layer 2: the cascade of what could separate a near-tied read.

  near_tie_classes.py READS.tsv SEL.bam TIE_SPLIT.tsv OWN_CALLS.tsv ANNOT.gff3.gz REF.fa FAMILY_NAMES.fa OUT_PREFIX

READS.tsv = neutral_truth.py score (NT_REPAIR=1); TIE_SPLIT.tsv = tie_split.py (TS_STRATA=S1,S2); OWN_CALLS = own_psv.py's
OP_CALLS (read, own call, loose call, any certified site); ANNOT = CAT/Liftoff slim gff3 (gene, transcript, exon).
Each read in the first class that applies:
  1 truth limit (outside / no near-perfect hit / competitor unlocatable / tie split ONE PLACE or UNSPLIT or not split)
  2 undecidable (tie, SEVERAL COPIES), split by span: the CHM13 genes at the two best-ms tied places have longest transcripts
    whose spliced sequences differ / are identical / a place has no gene (or both places are one gene)
  3 reference lead (decidable; the best-ms place leads by >= 1 and is the origin)
  4 own variant (decidable; own_psv's strict call is the origin)
  5 structure (decidable; origin and best competitor differ in intron count and the one with more introns is the origin)
  6 decidable, no signal
Beside: reference lead pointing to the wrong place; structure difference that would pick the wrong place.
Places are the truth's loci (reads.tsv, AS order); a place's ms / introns = its best-ms record (records overlapping it, same
contig and strand)."""
import bisect, collections, gzip, sys
import pysam

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def truth(path):
    T = {}
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if f[0] == 'read':
            continue
        loci = []
        for x in f[6].split(';'):
            i, c, se, st, a = x.rsplit(':', 4); s, e = se.split('-'); loci.append((c, int(s), int(e), st))
        T[f[0]] = dict(st=f[1], status=f[3], origin={int(i) for i in f[5].split(',') if i}, loci=loci)
    return T

def annot(path):
    tx = {}; ex = collections.defaultdict(list)
    for l in gzip.open(path, 'rt'):
        if l[0] == '#':
            continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        at = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if f[2] == 'transcript':
            tx[at['ID']] = (f[0], f[6], at.get('Parent', at.get('gene')))
        elif f[2] == 'exon':
            ex[at['Parent']].append((int(f[3]) - 1, int(f[4])))
    spans = collections.defaultdict(list); longest = {}
    for t, (c, st, g) in tx.items():
        e = sorted(ex.get(t, []))
        if not e:
            continue
        spans[c].append((e[0][0], e[-1][1], t))
        L = sum(b - a for a, b in e)
        if g not in longest or L > longest[g][0]:
            longest[g] = (L, t)
    idx = {}
    for c, v in spans.items():
        v.sort(); idx[c] = ([x[0] for x in v], v, max(x[1] - x[0] for x in v))
    return tx, ex, idx, longest

def gene_at(A, c, s, e, strand):
    tx, ex, idx, _ = A
    if c not in idx:
        return None
    starts, v, mx = idx[c]
    best = {}
    for k in range(bisect.bisect_left(starts, s - mx), bisect.bisect_left(starts, e)):
        ts, te, t = v[k]
        if te <= s:
            continue
        ov = sum(max(0, min(b, e) - max(a, s)) for a, b in ex[t])
        if ov <= 0:
            continue
        g = tx[t][2]; same = tx[t][1] == strand
        key = (same, ov)
        if g not in best or key > best[g]:
            best[g] = key
    if not best:
        return None
    return max(best, key=lambda g: best[g])

def tx_seq(A, fa, g):
    tx, ex, _, longest = A
    t = longest[g][1]; c, st, _ = tx[t]
    s = ''.join(fa.fetch(c, a, b) for a, b in sorted(ex[t])).upper()
    return s.translate(COMP)[::-1] if st == '-' else s

def main(reads_tsv, sel, split_tsv, own_tsv, gff, ref, famfa, out):
    T = truth(reads_tsv)
    fam = {l[1:].split()[0] for l in open(famfa) if l[0] == '>'}
    split = {}
    for l in open(split_tsv):
        f = l.rstrip('\n').split('\t')
        if f[0] != 'read':
            split[f[0]] = (f[2], [int(x) for x in f[7].split(',') if x])
    own = {}
    for l in open(own_tsv):
        f = l.rstrip('\n').split('\t')
        own[f[0]] = None if f[1] == 'NA' else int(f[1])
    R = collections.defaultdict(list)
    for r in pysam.AlignmentFile(sel, check_sq=False):
        if r.is_unmapped or r.is_supplementary or r.query_name not in T:
            continue
        nI = sum(1 for op, n in r.cigartuples if op == 3)
        R[r.query_name].append((r.reference_name, '-' if r.is_reverse else '+', r.reference_start, r.reference_end, r.get_tag('ms'), nI))
    A = annot(gff); fa = pysam.FastaFile(ref)
    C = collections.Counter(); rows = []
    for q, d in T.items():
        P = []
        for c, s, e, st in d['loci']:
            rr = [x for x in R[q] if x[0] == c and x[1] == st and x[2] < e and s < x[3]]
            b = max(rr, key=lambda x: x[4]) if rr else None
            P.append((b[4], b[5]) if b else (-1, -1))
        status = d['status']; sub = ''
        if status in ('outside', 'no near-perfect hit', 'competitor unlocatable'):
            cls = '1 truth limit'; sub = status
        elif status == 'tie':
            sp = split.get(q)
            if sp is None or sp[0] != 'SEVERAL COPIES':
                cls = '1 truth limit'; sub = 'tie ' + (sp[0] if sp else 'not split')
            else:
                cls = '2 undecidable'
                # a tie place outside the truth's listed loci (its extra places, negative ids) has no record here: left out
                tied = sorted((i for i in sp[1] if 0 <= i < len(P)), key=lambda i: -P[i][0])[:2]
                if len(tied) < 2:
                    sub = 'no gene'
                else:
                    gs = [gene_at(A, *d['loci'][i]) for i in tied]
                    if None in gs or gs[0] == gs[1]:
                        sub = 'no gene' if None in gs else 'one gene'
                    else:
                        sub = 'copies differ elsewhere' if tx_seq(A, fa, gs[0]) != tx_seq(A, fa, gs[1]) else 'identical transcripts'
        elif status == 'decidable' and next(iter(d['origin'])) >= len(P):
            # origin is a place the truth table does not list (a pick beyond its first loci): no score to compare
            o = next(iter(d['origin']))
            if own.get(q) is not None and own[q] == o:
                cls = '4 own variant'
            else:
                cls = '6 decidable, no signal'; sub = 'origin beyond the listed places'
        elif status == 'decidable':
            o = next(iter(d['origin']))
            order = sorted(range(len(P)), key=lambda i: -P[i][0])
            top = order[0]; lead = P[order[0]][0] - P[order[1]][0] if len(order) > 1 else 10 ** 6
            comp = next((i for i in order if i != o), None)
            if lead >= 1 and top != o:
                C['beside: reference lead points to the wrong place'] += 1
            sdiff = comp is not None and P[o][1] != P[comp][1]
            if sdiff and P[o][1] < P[comp][1]:
                C['beside: structure would pick the wrong place'] += 1
            if lead >= 1 and top == o:
                cls = '3 reference lead'
            elif own.get(q) is not None and own[q] == o:
                cls = '4 own variant'
            elif sdiff and P[o][1] > P[comp][1]:
                cls = '5 structure'
            else:
                cls = '6 decidable, no signal'
        else:
            cls = '1 truth limit'; sub = status
        C[cls] += 1; C[(cls, sub)] += 1; C[('stratum', d['st'])] += 1
        if q in fam:
            C[('family', cls)] += 1
            if status == 'tie' and q in split:
                C[('family split', split[q][0])] += 1
        rows.append((q, d['st'], status, cls, sub, int(q in fam)))
    with open(out + '.classes.tsv', 'w') as w:
        w.write('read\tstratum\tstatus\tclass\tsub\tfamily_read\n')
        for r in rows:
            w.write('\t'.join(map(str, r)) + '\n')
    n = len(rows); L = [f"reads {n}; strata " + ', '.join(f"{k[1]} {v}" for k, v in sorted(C.items(), key=lambda kv: str(kv[0])) if isinstance(k, tuple) and k[0] == 'stratum')]
    for cls in ('1 truth limit', '2 undecidable', '3 reference lead', '4 own variant', '5 structure', '6 decidable, no signal'):
        L.append(f"{cls}\t{C[cls]}\t{C[cls] / max(1, n):.3f}")
        for k, v in sorted(C.items(), key=lambda kv: str(kv[0])):
            if isinstance(k, tuple) and k[0] == cls and k[1]:
                L.append(f"    {k[1]}\t{v}\t{v / max(1, C[cls]):.3f}")
    for k in ('beside: reference lead points to the wrong place', 'beside: structure would pick the wrong place'):
        L.append(f"{k}\t{C[k]}")
    sc, op = C[('family split', 'SEVERAL COPIES')], C[('family split', 'ONE PLACE')]
    L.append(f"K2 family reads: SEVERAL COPIES {sc} / (SEVERAL + ONE PLACE {op}) = {sc / max(1, sc + op):.3f} (row 1791: 0.694; bar 0.60-0.80)")
    L.append('family reads by class: ' + ', '.join(f"{k[1]} {v}" for k, v in sorted(C.items(), key=lambda kv: str(kv[0])) if isinstance(k, tuple) and k[0] == 'family'))
    open(out + '.summary.txt', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))

if __name__ == '__main__':
    main(*sys.argv[1:9])
