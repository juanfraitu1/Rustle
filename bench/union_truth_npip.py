#!/usr/bin/env python3
"""Build a UNION truth for the chr16 NPIP family, per `docs/PREREG_union_truth_npip_2026-09-22.md` (§6x1).

`bench/mode_family_score.py` intersects every prediction with the truth universe, so a predicted member
carrying no Soto label is deleted from BOTH numerator and denominator -- it is invisible, not a false
positive. That INFLATES precision (register 770). This builds two larger truths so the arms can be
re-scored, and emits them in Soto's own two-column format so the scorer takes them as `--soto` unchanged.

  U0  Soto NPIP families on the chromosome (the shipped comparator)
  U1  U0 + genes whose RefSeq `Name` contains the family token
      ⚠ register 902-adjacent: symbol ROOTS are void as a genome-wide truth GENERATOR. Here the family is
        already defined by Soto and RefSeq only adds curated members to it. Kept separable on purpose.
  U2  U1 + genes where >= --min-frac of the gene's OWN length aligns to a U0 member, with a LENGTH FLOOR.
      Prediction-independent: only U0 member sequence and the annotation are consulted, never an arm's
      clusters -- scoring an arm against a truth derived from that arm is circular.

⚠⚠ TWO traps this script exists to avoid, both of which bit during development:

1. GENE NAMES ARE NOT UNIQUE (2 on chr16: CLN3, LOC102724181). Keying an alignment to an annotation by
   Name silently merges records: CLN3's 253 bp record aligns 100% while the real 16,713 bp CLN3 aligns 0%,
   and the merge gave the short element the long gene's length, carrying it past the floor. Every record is
   therefore addressed as `Name|chrom:start-end` end to end.
2. WITHOUT A LENGTH FLOOR the truth fills with short elements sitting INSIDE a duplicated block -- CLN3
   (253 bp), UBL5P4 (236 bp), PAWRP2 (634 bp) all reach 100% containment. That is register 913's "short
   genes become hubs" arriving through the truth instead of through the metric. The floor is the SHORTEST
   U0 member ("a family member is at least as long as the shortest known member"), so it is fixed by the
   truth before any score is seen, not chosen from the admitted list.

Added genes join the family of the U0 member they align to best (by merged aligned bp) -- the same
alignment, so no second criterion is introduced.

Usage:
  union_truth_npip.py --gff chr16.genes.gff --soto soto_famCN_S1C.tsv --fasta chr16.fa \
      --chrom chr16 --token NPIP --out-dir DIR
"""
import argparse, collections, csv, os, re, subprocess, sys


def merged_bp(spans):
    spans.sort(); out = []
    for s, e in spans:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return sum(e - s for s, e in out)


def gene_records(gff, chrom):
    """[(uid, name, start, end)] -- uid is unique, name is NOT (see trap 1)."""
    out = []
    for ln in open(gff):
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene', 'ncRNA_gene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            s, e = int(f[3]) - 1, int(f[4])
            out.append((f"{m.group(1)}|{chrom}:{s}-{e}", m.group(1), s, e))
    return out


def faidx(fasta, bed, fa):
    subprocess.run(['bedtools', 'getfasta', '-fi', fasta, '-bed', bed, '-nameOnly', '-fo', fa], check=True)


def main():
    ap = argparse.ArgumentParser()
    for a in ('--gff', '--soto', '--fasta', '--out-dir'):
        ap.add_argument(a, required=True)
    ap.add_argument('--chrom', default='chr16')
    ap.add_argument('--token', default='NPIP', help='family token matched against RefSeq Name')
    ap.add_argument('--min-frac', type=float, default=0.95)
    ap.add_argument('--threads', type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    D = lambda n: os.path.join(a.out_dir, n)

    fam = {}
    for r in csv.DictReader(open(a.soto), delimiter='\t'):
        f, g = (r.get('Family ID') or '').strip(), (r.get('Gene Name') or '').strip()
        if f and f != 'N/A' and g:
            fam.setdefault(g, f)

    recs = gene_records(a.gff, a.chrom)
    length = {uid: e - s for uid, _, s, e in recs}
    on_chrom = {n for _, n, _, _ in recs}
    fams = {f for g, f in fam.items() if g in on_chrom and a.token in g}
    U0 = {g for g, f in fam.items() if f in fams and g in on_chrom}
    named = {n for n in on_chrom if a.token in n}
    U1 = U0 | named
    floor = min(length[uid] for uid, n, _, _ in recs if n in U0)
    print(f"U0 {len(U0)} in {len(fams)} families | +{len(U1 - U0)} named | length floor {floor} bp "
          f"(shortest U0 member)", file=sys.stderr)

    with open(D('q.bed'), 'w') as q, open(D('r.bed'), 'w') as r:
        for uid, n, s, e in recs:
            q.write(f"{a.chrom}\t{s}\t{e}\t{uid}\n")
            if n in U0:
                r.write(f"{a.chrom}\t{s}\t{e}\t{uid}\n")
    faidx(a.fasta, D('q.bed'), D('q.fa'))
    faidx(a.fasta, D('r.bed'), D('r.fa'))
    with open(D('qr.paf'), 'w') as o:
        subprocess.run(['minimap2', '-c', '-N', '50', '-p', '0.1', '--secondary=yes', '-x', 'asm20',
                        '-t', str(a.threads), D('r.fa'), D('q.fa')], stdout=o, check=True)

    spans = collections.defaultdict(list)
    per_target = collections.defaultdict(lambda: collections.defaultdict(list))
    qlen = {}
    for ln in open(D('qr.paf')):
        f = ln.split('\t')
        qn, tn = f[0].split('|')[0], f[5].split('|')[0]
        if qn == tn:                                   # self, by NAME -- a gene never seeds itself
            continue
        qlen[f[0]] = int(f[1])
        spans[f[0]].append((int(f[2]), int(f[3])))
        per_target[qn][tn].append((int(f[2]), int(f[3])))
    admitted = {uid.split('|')[0] for uid, sp in spans.items()
                if length.get(uid, 0) >= floor and merged_bp(sp) / qlen[uid] >= a.min_frac}
    U2 = U1 | admitted
    print(f"U2 {len(U2)} (+{len(U2 - U1)} aligned at >= {a.min_frac:.0%} of own length)", file=sys.stderr)

    for tag, U in (('U1', U1), ('U2', U2)):
        with open(D(f'{tag}_truth.tsv'), 'w') as o:
            o.write('Gene Name\tFamily ID\n')
            for g in sorted(U):
                f = fam[g] if g in U0 else fam[max(per_target[g].items(),
                                                   key=lambda kv: merged_bp(kv[1]))[0]]
                o.write(f'{g}\t{f}\n')
        print(f"wrote {D(f'{tag}_truth.tsv')}", file=sys.stderr)


if __name__ == '__main__':
    main()
