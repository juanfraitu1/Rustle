#!/usr/bin/env python3
"""Row 1573: attach locus pieces to the UNCHANGED families (docs/PREREG_locus_pieces_attach_2026-11-07.md).

usage: piece_attach.py A_COPIES.tsv C_PREFIX CM_COPIES.tsv OUT_COPIES.tsv [--keep-all]
  A_COPIES  the families run on the uncut input (its families are kept as they are)
  C_PREFIX  the run on the input with pieces (its .fam.loci.gated.paf links pieces)
  CM_COPIES C's copy table after piece_merge.py
       piece_attach.py --selftest
"""
import collections, re, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import piece_merge as PM

orig = lambda g: re.sub(r'\.p\d+$', '', g)


def read(path):
    with open(path) as fh:
        hdr = fh.readline()
        cols = hdr.rstrip('\n').split('\t')
        return hdr, [(l, dict(zip(cols, l.rstrip('\n').split('\t')))) for l in fh]


def attach(a_rows, c_rows, linked):
    """a_rows/c_rows: [(line, dict)]. Returns the output lines (A's families, pieces attached)."""
    a_members = collections.defaultdict(set)
    for _, r in a_rows:
        a_members[r['family_id']].add(orig(r['gene_id']))
    c_members = collections.defaultdict(set)
    for _, r in c_rows:
        c_members[r['family_id']].add(orig(r['gene_id']))
    placed = collections.defaultdict(list)          # (A family, parent locus) -> [(line, dict)] of pieces
    for line, r in c_rows:
        if not re.search(r'\.p\d+$', r['gene_id']):
            continue
        parent = orig(r['gene_id'])
        mates = c_members[r['family_id']] - {parent}
        best = max(sorted(a_members), key=lambda f: len(a_members[f] & mates), default=None)
        if best is None or not a_members[best] & mates:
            continue
        placed[(best, parent)].append((line, r))
    replaced = set(placed)
    out = [line for line, r in a_rows if (r['family_id'], orig(r['gene_id'])) not in replaced]
    for (fam, parent), ps in placed.items():
        rows = [dict(i=i, fam=fam, orig=parent, key=f"{r['chrom']}:{int(r['locus_start']) + 1}-{r['locus_end']}",
                     reads=float(r['n_reads'])) for i, (_, r) in enumerate(ps)]
        for i in sorted(PM.merge(rows, linked)):
            line, r = ps[i]
            out.append(line.replace(r['family_id'], fam, 1) if line.startswith(r['family_id'] + '\t') else line)
    return out


def main(a_copies, c_prefix, cm_copies, out_path, keep_all=False):
    hdr, a_rows = read(a_copies)
    _, c_rows = read(cm_copies)
    keys = {f"{r['chrom']}:{int(r['locus_start']) + 1}-{r['locus_end']}" for _, r in c_rows if re.search(r'\.p\d+$', r['gene_id'])}
    pairs = set()
    for line in open(c_prefix + '.fam.loci.gated.paf'):
        f = line.split('\t', 12)
        if f[0] in keys and f[5] in keys and f[0] != f[5] and int(f[10]) >= PM.MIN_ALEN and int(f[9]) / int(f[10]) >= PM.MIN_ID:
            pairs.add(frozenset((f[0], f[5])))
    # --keep-all (row 1577): exon-disjoint groups are separate transcription units; keep every attached group
    out = attach(a_rows, c_rows, (lambda a, b: True) if keep_all else (lambda a, b: frozenset((a, b)) in pairs))
    with open(out_path, 'w') as oh:
        oh.write(hdr)
        oh.writelines(out)
    print(f"[piece_attach] A copies {len(a_rows)} -> {len(out)}", file=sys.stderr)


def selftest():
    mk = lambda fam, gid, s, reads=5: (f"{fam}\t{gid}\tchr1\t{s}\n",
                                       dict(family_id=fam, gene_id=gid, chrom='chr1', locus_start=str(s), locus_end=str(s + 10), n_reads=str(reads)))
    a = [mk('A1', 'L', 0), mk('A1', 'M', 100), mk('A1', 'N', 200), mk('A2', 'Q', 300)]
    # L's piece p2 sits with M, N in C (-> A1, replacing L's row); p1 with Q (-> A2); a piece with no A mate is dropped
    c = [mk('C1', 'L.p2', 20), mk('C1', 'M', 100), mk('C1', 'N', 200), mk('C2', 'L.p1', 0), mk('C2', 'Q', 300),
         mk('C3', 'Z.p1', 900), mk('C3', 'Y', 950)]
    out = attach(a, c, lambda x, y: False)
    fams = sorted(l.split('\t')[:2] for l in out)
    assert fams == [['A1', 'L.p2'], ['A1', 'M'], ['A1', 'N'], ['A2', 'L.p1'], ['A2', 'Q']], fams
    # two unlinked pieces of L in one target family: one kept by default, both with keep_all (row 1577)
    c2 = [mk('C1', 'L.p1', 0, 9), mk('C1', 'L.p2', 20, 3), mk('C1', 'M', 100)]
    assert sorted(l.split('\t')[1] for l in attach(a, c2, lambda x, y: False)) == ['L.p1', 'M', 'N', 'Q']
    assert sorted(l.split('\t')[1] for l in attach(a, c2, lambda x, y: True)) == ['L.p1', 'L.p2', 'M', 'N', 'Q']
    print('selftest ok')


if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        selftest()
    else:
        main(*[x for x in sys.argv[1:] if x != '--keep-all'][:4], keep_all='--keep-all' in sys.argv)
