#!/usr/bin/env python3
"""loop_home.py -- the closed loop's read-home table and its development scoring.

docs/archive/2026-09/PREREG_tied_read_loop_2026-09-25.md: families -> copy assignment of tied reads (union certificate) -> each
assigned read given to its copy -> re-assembly with RUSTLE_READ_HOME_TABLE (src/rustle/vg_family/denovo_assemble.rs,
`ReadHomeTable`). This script builds the table the assembler reads and the counts the pre-registration names.

    loop_home.py union   --union U.union_certificate.tsv --families UNITS.tsv --bam B --out HOME.tsv [--summary S.tsv]
        one row per molecule the union certificate ASSIGNED (verdict assigned_family: winner `family_id:copy_idx` ->
        that unit's span; assigned_outside: winner `outside:chrom:start1-end`). tied / ambiguous / no_result rows are
        never moved (prereg §1.2). A molecule with no mapped non-supplementary record whose aligned block overlaps its
        home is left OUT of the table and counted as H (its placement would need a lifted record; prereg §1.3).
    loop_home.py oracle  --sim-dir DIR --out HOME.tsv [--summary S.tsv]
        the ORACLE arm (simulation only): every simulated molecule `family|copy_idx|i` of DIR/sim.bam gets its true
        source copy (DIR/cat.copies.tsv) as home; same presence check and exclusion as `union`.
    loop_home.py g0      --log PREFIX.assemble.log [--home HOME.tsv]
        the G0 counts (prereg §2) from the filter's own log: kept at home, R_add, R_drop; bounded iff R_drop + R_add < 20.
    loop_home.py sim-truth-gtf --sim-dir DIR --out TRUTH.gtf
        the simulated source copies (>= 1 simulated read) as a GTF, for M5.
    loop_home.py score-sim --sim-dir DIR --arm NAME=GTF [...] --out M5.tsv [--gffcompare BIN] [--work DIR]
        M5 (prereg §3): per arm, multi-exon copies whose exact intron chain an arm transcript has (gffcompare `=`),
        single-exon copies covered >= 0.5 of their exon bases by an arm locus; precision = arm transcripts matching a
        copy (`=`) / arm transcripts overlapping a simulated copy's span.

Coordinates: copies.tsv / units.tsv `start`, `end`, `exons` are 0-based half-open (the Rust contract); the table is too.
"""
from __future__ import annotations

import argparse
import collections
import csv
import os
import re
import shutil
import subprocess
import sys

G0_FLOOR = 20  # prereg §2: fewer than 20 records moved = the sample is bounded


def read_units(path):
    """copies.tsv / units.tsv (the `--families` contract), by header name: {(family_id, copy_idx): row}."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            r['start'], r['end'] = int(r['start']), int(r['end'])
            r['exon_list'] = [tuple(map(int, b.split('-'))) for b in r['exons'].split(',') if b.strip()]
            out[(r['family_id'], str(int(r['copy_idx'])))] = r
    return out


def parse_winner(winner, units):
    """(chrom, start, end, kind) of a union winner label, or None when it cannot be resolved."""
    if winner.startswith('outside:'):
        chrom, _, span = winner[len('outside:'):].rpartition(':')
        m = re.fullmatch(r'(\d+)-(\d+)', span)
        if not chrom or not m:
            return None
        return chrom, int(m.group(1)) - 1, int(m.group(2)), 'outside'   # label start is 1-based
    fam, _, idx = winner.rpartition(':')
    if not fam or not idx.isdigit():          # the `fid:#ci` fallback label names no catalog row
        return None
    u = units.get((fam, idx))
    return (u['chrom'], u['start'], u['end'], 'family') if u else None


def present_at_home(bam_path, homes):
    """{read_name} of the molecules with >= 1 mapped non-supplementary record whose aligned block (M/=/X) overlaps
    their home span: the rule the assembler applies (`aligned_blocks_overlap`). homes: {read: [(chrom, s, e)]}."""
    import pysam
    by_home = collections.defaultdict(set)
    for name, hs in homes.items():
        for h in hs:
            by_home[h].add(name)
    found = set()
    with pysam.AlignmentFile(bam_path) as bam:
        contigs = set(bam.references)
        for (chrom, s, e), names in sorted(by_home.items()):
            if chrom not in contigs:
                continue
            for a in bam.fetch(chrom, s, e):
                if a.is_unmapped or a.is_supplementary or a.query_name not in names or a.query_name in found:
                    continue
                if any(b0 < e and b1 > s for b0, b1 in a.get_blocks()):
                    found.add(a.query_name)
    return found


def write_home(out, rows):
    with open(out, 'w') as fh:
        fh.write('read_name\tchrom\tstart\tend\tsource\twinner\n')
        for r in sorted(rows):
            fh.write('\t'.join(map(str, r)) + '\n')


def write_summary(path, counts):
    lines = ''.join(f'{k}\t{v}\n' for k, v in counts.items())
    sys.stderr.write(lines)
    if path:
        with open(path, 'w') as fh:
            fh.write('item\tn\n' + lines)


def cmd_union(a):
    units = read_units(a.families)
    verdicts = collections.Counter()
    homes, meta, unresolved = {}, {}, 0
    with open(a.union) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            verdicts[r['verdict']] += 1
            if r['verdict'] not in ('assigned_family', 'assigned_outside'):
                continue
            w = parse_winner(r['winner'], units)
            if w is None:
                unresolved += 1
                continue
            homes[r['read_name']] = [w[:3]]
            meta[r['read_name']] = (w[3], r['winner'])
    ok = set(homes) if a.no_bam_check else present_at_home(a.bam, homes)
    rows = [(n, *homes[n][0], 'union_' + meta[n][0], meta[n][1]) for n in homes if n in ok]
    write_home(a.out, rows)
    counts = collections.OrderedDict()
    for v in ('assigned_family', 'assigned_outside', 'tied', 'ambiguous', 'no_result'):
        counts[f'verdict_{v}'] = verdicts.get(v, 0)
    counts['A_assigned'] = verdicts.get('assigned_family', 0) + verdicts.get('assigned_outside', 0)
    counts['winner_unresolved'] = unresolved
    counts['H_no_record_at_home'] = len(homes) - len(ok & set(homes))
    counts['home_table_molecules'] = len(rows)
    counts['bam_checked'] = 'no' if a.no_bam_check else 'yes'
    write_summary(a.summary, counts)


def cmd_oracle(a):
    import pysam
    d = a.sim_dir
    units = read_units(a.families or os.path.join(d, 'cat.copies.tsv'))
    bam = a.bam or os.path.join(d, 'sim.bam')
    homes, bad = {}, 0
    with pysam.AlignmentFile(bam) as fh:
        for rec in fh.fetch(until_eof=True):
            name = rec.query_name
            if name in homes:
                continue
            parts = name.split('|')
            u = units.get((parts[0], parts[1])) if len(parts) == 3 else None
            if u is None:
                bad += 1
                continue
            homes[name] = [(u['chrom'], u['start'], u['end'])]
    ok = set(homes) if a.no_bam_check else present_at_home(bam, homes)
    rows = [(n, *homes[n][0], 'oracle', n.rsplit('|', 1)[0].replace('|', ':')) for n in homes if n in ok]
    write_home(a.out, rows)
    write_summary(a.summary, collections.OrderedDict(
        [('simulated_molecules', len(homes)), ('name_not_a_catalog_copy', bad),
         ('H_no_record_at_home', len(homes) - len(ok)), ('home_table_molecules', len(rows))]))


G0_RE = re.compile(r'totals kept_at_home (\d+) R_add (\d+) R_drop (\d+) dropped_away_outside_pool (\d+)')


def g0_counts(log):
    last = None
    with open(log) as fh:
        for line in fh:
            if line.startswith('[read-home]'):
                m = G0_RE.search(line)
                if m:
                    last = tuple(map(int, m.groups()))
    return last or (0, 0, 0, 0)


def cmd_g0(a):
    kept, add, drop, inert = g0_counts(a.log)
    n_home = sum(1 for _ in open(a.home)) - 1 if a.home else None
    moved = add + drop
    print(f'home_table_molecules\t{n_home if n_home is not None else "NA"}\nkept_at_home\t{kept}\nR_add\t{add}\n'
          f'R_drop\t{drop}\ndropped_away_outside_pool\t{inert}\nR_drop_plus_R_add\t{moved}\n'
          f'G0\t{"bounded (< %d records moved)" % G0_FLOOR if moved < G0_FLOOR else "open"}')


# ---------------------------------------------------------------- M5: scoring against the simulated copies
def simulated_copies(d):
    units = read_units(os.path.join(d, 'cat.copies.tsv'))
    out = []
    with open(os.path.join(d, 'sim.copies_used.tsv')) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if int(r['n_sim']) >= 1:
                u = units.get((r['family_id'], str(int(r['copy_idx']))))
                if u:
                    out.append(u)
    return out


def cmd_sim_truth_gtf(a):
    with open(a.out, 'w') as fh:
        for u in simulated_copies(a.sim_dir):
            tid = f"{u['family_id']}:{int(u['copy_idx'])}"
            attr = f'gene_id "{tid}"; transcript_id "{tid}";'
            ex = u['exon_list']
            fh.write(f"{u['chrom']}\tsim\ttranscript\t{ex[0][0] + 1}\t{ex[-1][1]}\t.\t{u['strand']}\t.\t{attr}\n")
            for s, e in ex:
                fh.write(f"{u['chrom']}\tsim\texon\t{s + 1}\t{e}\t.\t{u['strand']}\t.\t{attr}\n")


def gtf_transcripts(path):
    """{tid: (chrom, [(s0, e)] sorted)} from a GTF's exon lines (0-based half-open)."""
    tx = collections.defaultdict(list)
    chrom = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'exon':
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if m:
                tx[m.group(1)].append((int(f[3]) - 1, int(f[4])))
                chrom[m.group(1)] = f[0]
    return {t: (chrom[t], sorted(v)) for t, v in tx.items()}


def _cov(blocks, of):
    """bases of `of` covered by the union of `blocks`"""
    tot = 0
    for s, e in of:
        pts = sorted((max(s, b0), min(e, b1)) for b0, b1 in blocks if b0 < e and b1 > s)
        cur = s
        for x0, x1 in pts:
            x0 = max(x0, cur)
            if x1 > x0:
                tot += x1 - x0
                cur = x1
    return tot


def cmd_score_sim(a):
    d = a.sim_dir
    work = a.work or os.path.join(d, 'loop_score')
    os.makedirs(work, exist_ok=True)
    truth = os.path.join(work, 'truth.gtf')
    cmd_sim_truth_gtf(argparse.Namespace(sim_dir=d, out=truth))
    copies = simulated_copies(d)
    multi = {f"{u['family_id']}:{int(u['copy_idx'])}" for u in copies if len(u['exon_list']) >= 2}
    single = [u for u in copies if len(u['exon_list']) == 1]
    spans = collections.defaultdict(list)
    for u in copies:
        spans[u['chrom']].append((u['start'], u['end']))
    gffc = a.gffcompare or shutil.which('gffcompare') or '/home/juanfra/miniforge3/bin/gffcompare'
    rows = []
    for spec in a.arm:
        name, _, gtf = spec.partition('=')
        q = os.path.join(work, f'{name}.gtf')
        shutil.copyfile(gtf, q)
        subprocess.run([gffc, '-r', truth, '-o', os.path.join(work, f'cmp_{name}'), q], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        tmap = os.path.join(work, f'cmp_{name}.{name}.gtf.tmap')
        matched_ref, matched_q = set(), set()
        with open(tmap) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if r['class_code'] == '=':
                    matched_ref.add(r['ref_id'])
                    matched_q.add(r['qry_id'])
        tx = gtf_transcripts(q)
        n_multi_hit = len(multi & matched_ref)
        by_chrom = collections.defaultdict(list)
        for c, ex in tx.values():
            by_chrom[c].extend(ex)
        n_single_hit = sum(1 for u in single
                           if _cov(by_chrom.get(u['chrom'], []), u['exon_list']) >= 0.5 * sum(e - s for s, e in u['exon_list']))
        overl = [t for t, (c, ex) in tx.items()
                 if any(ex[0][0] < e and ex[-1][1] > s for s, e in spans.get(c, []))]
        n_prec = sum(1 for t in overl if t in matched_q)
        rows.append((name, len(tx), len(multi), n_multi_hit, len(single), n_single_hit, len(overl), n_prec))
    with open(a.out, 'w') as fh:
        fh.write('arm\tn_transcripts\tn_multi_exon_copies\tmulti_exon_chain_found\tn_single_exon_copies\t'
                 'single_exon_covered\tn_transcripts_at_copies\tn_transcripts_matching_a_copy\n')
        for r in rows:
            fh.write('\t'.join(map(str, r)) + '\n')
    sys.stderr.write(open(a.out).read())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    u = sub.add_parser('union')
    u.add_argument('--union', required=True)
    u.add_argument('--families', required=True, help='the units/copies TSV the assignment ran on')
    u.add_argument('--bam', help='the BAM the assignment ran on (presence-at-home check)')
    u.add_argument('--no-bam-check', action='store_true', help='skip the presence check (tests only)')
    u.add_argument('--out', required=True)
    u.add_argument('--summary')
    o = sub.add_parser('oracle')
    o.add_argument('--sim-dir', required=True)
    o.add_argument('--families')
    o.add_argument('--bam')
    o.add_argument('--no-bam-check', action='store_true')
    o.add_argument('--out', required=True)
    o.add_argument('--summary')
    g = sub.add_parser('g0')
    g.add_argument('--log', required=True)
    g.add_argument('--home')
    t = sub.add_parser('sim-truth-gtf')
    t.add_argument('--sim-dir', required=True)
    t.add_argument('--out', required=True)
    s = sub.add_parser('score-sim')
    s.add_argument('--sim-dir', required=True)
    s.add_argument('--arm', action='append', required=True, metavar='NAME=GTF')
    s.add_argument('--out', required=True)
    s.add_argument('--gffcompare')
    s.add_argument('--work')
    a = ap.parse_args(argv)
    if a.cmd == 'union' and not a.no_bam_check and not a.bam:
        ap.error('union needs --bam (or --no-bam-check)')
    {'union': cmd_union, 'oracle': cmd_oracle, 'g0': cmd_g0, 'sim-truth-gtf': cmd_sim_truth_gtf,
     'score-sim': cmd_score_sim}[a.cmd](a)


if __name__ == '__main__':
    main()
