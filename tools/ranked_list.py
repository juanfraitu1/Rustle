#!/usr/bin/env python3
"""ranked_list.py FINAL.gtf ONEREAD_RANKED.gtf OUT.gtf [--sensitive-share 0.1]

One ranked list of an assembly (register rows 1880 / 1884; user decision 2026-11-22):
  tier 1  every transcript of FINAL.gtf, most `reads` first (ties keep the file's order);
  tier 2  every one-read chain of ONEREAD_RANKED.gtf whose intron chain FINAL.gtf lacks, best `oneread_score` first.
Each transcript line gets  list_tier "1|2"; list_rank "N";  and  list_sensitive "1";  when it belongs to the named cut
"sensitive" = tier 1 plus the best 10% of tier 2 (--sensitive-share). No transcript is changed, added or removed; the list is
a reading order. Cut at another assembler's output size it held at least that assembler's number of annotated intron chains
on 65 of 67 comparisons over twenty datasets (row 1880). Tier 2 is a recall tier: one read each.

Prints one summary line. Exit 2 on unusable input."""
import re
import sys

TID = re.compile(r'transcript_id "([^"]*)"')


def blocks(path):
    """Transcripts of a GTF in file order: [transcript line or None, [other lines]], header comments apart."""
    order, by_id, head = [], {}, []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                head.append(line)
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            m = TID.search(cols[8])
            if not m:
                continue
            key = (cols[0], cols[6], m.group(1))
            if key not in by_id:
                by_id[key] = [None, []]
                order.append(key)
            if cols[2] == 'transcript':
                by_id[key][0] = line
            else:
                by_id[key][1].append(line)
    return head, [(k, by_id[k]) for k in order]


def chain(key, block):
    """(contig, strand, introns) of a transcript; introns = () for a single exon."""
    ex = sorted((int(c[3]), int(c[4])) for c in (l.split('\t') for l in block[1]) if c[2] == 'exon')
    return key[0], key[1], tuple((ex[i][1], ex[i + 1][0] - 1) for i in range(len(ex) - 1))


def number(line, name):
    m = re.search(r'\b' + name + r' "([^"]*)"', line or '')
    try:
        return float(m.group(1)) if m else 0.0
    except ValueError:
        return 0.0


def tag(line, tier, rank, sensitive):
    extra = f' list_tier "{tier}"; list_rank "{rank}";' + (' list_sensitive "1";' if sensitive else '')
    return line.rstrip('\n').rstrip() + extra + '\n'


def main(argv):
    args = [a for a in argv if not a.startswith('--')]
    share = 0.1
    if '--sensitive-share' in argv:
        try:
            share = float(argv[argv.index('--sensitive-share') + 1])
        except (IndexError, ValueError):
            share = -1.0
    if len(args) != 3 or not 0.0 <= share <= 1.0:
        sys.stderr.write(__doc__)
        return 2
    final, ranked, out = args
    _, t1 = blocks(final)
    _, t2 = blocks(ranked)
    t1 = [(k, b) for k, b in t1 if b[0] is not None]
    t2 = [(k, b) for k, b in t2 if b[0] is not None]
    if not t1:
        sys.stderr.write(f'ranked_list: no transcript line in {final}\n')
        return 2
    have = {chain(k, b) for k, b in t1}
    first = sorted(range(len(t1)), key=lambda i: (-number(t1[i][1][0], 'reads'), i))
    extra = [(k, b) for k, b in t2 if chain(k, b) not in have]
    second = sorted(range(len(extra)), key=lambda i: (-number(extra[i][1][0], 'oneread_score'), i))
    n_sens = int(share * len(second))
    with open(out, 'w') as f:
        f.write(f'# ranked list (tools/ranked_list.py): tier 1 = {len(first)} transcripts of the default output by reads; '
                f'tier 2 = {len(second)} one-read chains by oneread_score; list_sensitive = tier 1 + the best {n_sens} of tier 2\n')
        rank = 0
        for i in first:
            rank += 1
            f.write(tag(t1[i][1][0], 1, rank, True))
            f.writelines(t1[i][1][1])
        for j, i in enumerate(second):
            rank += 1
            f.write(tag(extra[i][1][0], 2, rank, j < n_sens))
            f.writelines(extra[i][1][1])
    print(f'{len(first)} transcripts in tier 1, {len(second)} one-read chains in tier 2 '
          f'({len(t2) - len(extra)} already in tier 1), sensitive cut = tier 1 + {n_sens}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
