#!/usr/bin/env python3
"""family_container.py -- the CONTAINER of a family member's extra pieces (accessory exon blocks) and their relations.

docs/archive/2026-09/PREREG_fusion_container_sim_2026-09-28.md §1 (the binding definition), implemented as a post-processor of the
driver's `families` stage (`tools/rustle_pipeline.sh families` = `mcl_families --from-gtf ... --out PREFIX.fam`).
It reads the families products and never changes a family.

    family_container.py --gtf PREFIX.gtf --fam PREFIX.fam --out OUT
        inputs  PREFIX.fam.clusters.tsv, PREFIX.fam.loci.gff3, PREFIX.fam.loci.tsv, PREFIX.fam.loci.paf
                (each overridable: --clusters / --loci-gff3 / --loci-tsv / --paf; a PAF may be .gz)
        outputs OUT.blocks.tsv     one row per (clustered locus, exon block)
                OUT.relations.tsv  one row per directed family relation F -> F'
                OUT.summary.tsv    counts (also printed on stderr)

Definition (prereg §1):
  * LOCUS m = a member row of clusters.tsv (a graph node key CONTIG:START-END) plus every annotation record that
    loci.tsv folds into it (`--fold-within-clusters`, the shipped default); FAMILY F = its cluster_id.
    A record is one assembled gene_id: loci.gff3's `gene` line gives its key and `Name=` its gene_id (two gene_ids
    with the same span share one key and are both taken).
  * EXON BLOCKS of m = the union of the exons of ALL transcripts of ALL gene_ids of m's records (not only the
    representative that loci.gff3 lists), merged where they OVERLAP (share >= 1 base; abutting exons stay separate).
  * Block b of m is CORE iff some PAF record between a record of m and a record of another member m' of F has an
    aligned column (CIGAR M/=/X) whose m-side base lies in b and whose m'-side base is an exon base of m'. The exon
    bases of m' are m''s own exon blocks (the same all-transcript union). No other constant: every PAF record counts,
    whatever its identity, length or primary/secondary flag. ACCESSORY = not core.
  * RELATION: for each accessory block, the same test against the members of every OTHER family F'; each hit
    records (m, b, F'). The family relation F -> F' exists iff some member of F carries such a block; `reciprocal`
    says whether F' -> F exists too (the undirected F~F' is the union of both directions). Core blocks are not
    tested for relations (their relation columns are '.').
  * Unclustered loci get no rows and are never partners. Records between two records of the same locus are skipped.

⚠ The shipped edge rule (`--exonic-both-sides`, `annotation_families::graph_from_paf_loci`) tests a record's whole
query and target INTERVALS against exons; the prereg asks for the test "projected through the CIGAR", so this
container tests COLUMNS: a record whose query interval touches an exon of m and whose target interval touches an
exon of m' is still not evidence unless one aligned column joins the two exon bases.

Coordinate contract (checked against src/bin/mcl_families.rs `loci_from_gtf`):
  * loci.fa / PAF sequence names are record keys `CONTIG:START-END`, GFF 1-based closed; the sequence is the genome's
    FORWARD strand from START to END whatever the locus strand (`fetch_sequence(c, START-1, END)`), so the 0-based PAF
    offset o is the 1-based genome base START + o (0-based START - 1 + o).
  * PAF `+`: the CIGAR walks target [ts,te) and query [qs,qe) both ascending. PAF `-`: the CIGAR walks target [ts,te)
    ascending against the REVERSE COMPLEMENT of query [qs,qe), i.e. forward query offsets from qe-1 downwards.
  * M/=/X consume both sides (aligned columns); I consumes the query only; D and N consume the target only. Any other
    op, a CIGAR whose consumed lengths disagree with the PAF columns, or a projected record without cg:Z is an error.
  * Output start/end are GFF 1-based closed (like clusters.tsv and the GTF); block numbers are 1-based in GENOMIC
    order (not transcription order). Internally everything is 0-based half-open genome coordinates.

Tests: bench/test_family_container.py (hand-made GTF + PAF + clusters fixtures). Only the standard library is used.
"""
import argparse
import bisect
import collections
import gzip
import os
import re
import sys

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

BLOCK_HEADER = ['family_id', 'locus', 'chrom', 'strand', 'gene_ids', 'n_records', 'block', 'n_blocks', 'start', 'end',
                'bp', 'class', 'core_bp', 'core_partners', 'rel_families', 'rel_bp', 'rel_partners']
REL_HEADER = ['family_id', 'related_family', 'n_members', 'n_blocks', 'bp', 'reciprocal', 'members']


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def gtf_attr(s, key):
    """The value of `key "..."` in a GTF attribute column (the same substring rule as mcl_families `gtf_loci`)."""
    pat = key + ' "'
    i = s.find(pat)
    if i < 0:
        return None
    i += len(pat)
    j = s.find('"', i)
    return s[i:j] if j >= 0 else None


def parse_key(name):
    """`CONTIG:START-END` -> (contig, start, end) (annotation_families::parse_gene_key), None if malformed."""
    c, sep, r = name.rpartition(':')
    if not sep:
        return None
    a, sep, b = r.partition('-')
    if not sep:
        return None
    try:
        return (c, int(a), int(b))
    except ValueError:
        return None


def key_str(k):
    return f'{k[0]}:{k[1]}-{k[2]}'


def merge_blocks(intervals):
    """Merge 0-based half-open intervals that OVERLAP (share >= 1 base); abutting intervals stay separate."""
    out = []
    for s, e in sorted(intervals):
        if out and s < out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def union_len(intervals):
    """Bases covered by a list of 0-based half-open intervals (overlapping or abutting)."""
    n, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                n += cur_e - cur_s
            cur_s, cur_e = s, e
        elif e > cur_e:
            cur_e = e
    if cur_e is not None:
        n += cur_e - cur_s
    return n


def aligned_runs(cigar, strand, qs, qe, ts, te):
    """The aligned runs of one PAF record, in the record's own offset frame.

    Yields (t0, q0, n): n aligned columns; column i (0 <= i < n) joins target offset t0 + i with query offset q0 + i on
    `+`, and with query offset q0 + n - 1 - i on `-` (so [q0, q0 + n) is the run's query interval on both strands).
    """
    ops = CIGAR_RE.findall(cigar)
    if ''.join(n + op for n, op in ops) != cigar:
        raise ValueError(f'unparseable CIGAR {cigar[:60]!r}')
    if strand not in '+-' or len(strand) != 1:
        raise ValueError(f'PAF strand {strand!r}')
    t, qc = ts, 0  # target position; query bases consumed so far (from qs on '+', from qe downwards on '-')
    for n, op in ops:
        n = int(n)
        if op in 'M=X':
            yield (t, qs + qc if strand == '+' else qe - qc - n, n)
            t += n
            qc += n
        elif op == 'I':
            qc += n
        elif op in 'DN':
            t += n
        else:
            raise ValueError(f'CIGAR op {op} not allowed in a PAF record')
    if t != te or qc != qe - qs:
        raise ValueError(f'CIGAR consumes target {t - ts} / query {qc} but the record spans {te - ts} / {qe - qs}')


def mask_ranges(starts, ends, a, b):
    """(lo, hi, block_index) for every block of a locus (sorted, disjoint) intersecting [a, b)."""
    out = []
    i = bisect.bisect_right(ends, a)
    while i < len(starts) and starts[i] < b:
        lo, hi = max(starts[i], a), min(ends[i], b)
        if hi > lo:
            out.append((lo, hi, i))
        i += 1
    return out


def overlaps_any(starts, ends, a, b):
    i = bisect.bisect_right(ends, a)
    return i < len(starts) and starts[i] < b


def exon_columns(cigar, strand, qs, qe, ts, te, q_off, t_off, q_blocks, t_blocks):
    """Aligned columns of one record joining an exon base of the query locus to an exon base of the target locus.

    q_off / t_off: the genome 0-based coordinate of offset 0 of the query / target record. q_blocks / t_blocks:
    (starts, ends) of the query / target LOCUS blocks, genome 0-based half-open. Yields
    (q_lo, q_hi, q_block, t_lo, t_hi, t_block) in genome coordinates; the two intervals have the same length and are
    joined column by column (reversed on `-`).
    """
    qst, qen = q_blocks
    tst, ten = t_blocks
    for t0, q0, n in aligned_runs(cigar, strand, qs, qe, ts, te):
        tg, qg = t_off + t0, q_off + q0
        tr = [(lo - tg, hi - tg, bi) for lo, hi, bi in mask_ranges(tst, ten, tg, tg + n)]
        if not tr:
            continue
        qm = mask_ranges(qst, qen, qg, qg + n)
        if not qm:
            continue
        if strand == '+':
            qr = [(lo - qg, hi - qg, bi) for lo, hi, bi in qm]
        else:  # column i <-> query qg + n - 1 - i, so query [lo, hi) <-> i in [qg + n - hi, qg + n - lo)
            qr = [(qg + n - hi, qg + n - lo, bi) for lo, hi, bi in reversed(qm)]
        x = y = 0
        while x < len(qr) and y < len(tr):
            lo, hi = max(qr[x][0], tr[y][0]), min(qr[x][1], tr[y][1])
            if hi > lo:
                if strand == '+':
                    ql, qh = qg + lo, qg + hi
                else:
                    ql, qh = qg + n - hi, qg + n - lo
                yield (ql, qh, qr[x][2], tg + lo, tg + hi, tr[y][2])
            if qr[x][1] <= tr[y][1]:
                x += 1
            else:
                y += 1


# ---------------------------------------------------------------------------------------------------------- inputs

def read_clusters(path):
    """clusters.tsv -> (family of each member key, members per family in file order, family order)."""
    fam_of, members, order = {}, collections.OrderedDict(), {}
    with open_text(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        col = {h: i for i, h in enumerate(header)}
        for need in ('cluster_id', 'chrom', 'start', 'end'):
            if need not in col:
                raise ValueError(f'{path}: no {need} column (header {header})')
        for line in fh:
            r = line.rstrip('\n').split('\t')
            if len(r) < len(header):
                continue
            fid = r[col['cluster_id']]
            key = (r[col['chrom']], int(r[col['start']]), int(r[col['end']]))
            if key in fam_of:
                if fam_of[key] != fid:
                    raise ValueError(f'{path}: {key_str(key)} is in {fam_of[key]} and {fid} (not a strict partition)')
                continue
            fam_of[key] = fid
            order.setdefault(fid, len(order))
            members.setdefault(fid, []).append(key)
    return fam_of, members, order


def read_folds(path):
    """loci.tsv (`annotation representative`, keys CONTIG:START-END) -> {annotation key: representative key}."""
    folds = {}
    with open_text(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        if header[:2] != ['annotation', 'representative']:
            raise ValueError(f'{path}: header {header} is not `annotation representative`')
        for line in fh:
            r = line.rstrip('\n').split('\t')
            if len(r) < 2:
                continue
            a, b = parse_key(r[0]), parse_key(r[1])
            if a is None or b is None:
                raise ValueError(f'{path}: malformed keys {r[:2]}')
            if a != b:
                folds[a] = b
    return folds


def read_loci_gff3(path):
    """loci.gff3 `gene` lines -> ({key: [gene_id, ...]}, {key: strand of its first gene line})."""
    genes, strand = collections.defaultdict(list), {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            r = line.rstrip('\n').split('\t')
            if len(r) < 9 or r[2] != 'gene':
                continue
            key = (r[0], int(r[3]), int(r[4]))
            name = next((kv[5:] for kv in r[8].split(';') if kv.startswith('Name=')), None)
            if name is None:
                raise ValueError(f'{path}: gene line without Name=: {line.strip()}')
            genes[key].append(name)
            strand.setdefault(key, r[6])
    return genes, strand


def read_gtf(path, wanted_genes):
    """Assembled GTF -> {gene_id: [(chrom, start1, end1), ...]} over all its transcripts, for the wanted gene_ids.

    Transcript -> gene comes from `transcript` lines and exons join by transcript_id, as mcl_families `gtf_loci` does.
    """
    gene_of, exons = {}, collections.defaultdict(list)
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            r = line.rstrip('\n').split('\t')
            if len(r) < 9:
                continue
            t = gtf_attr(r[8], 'transcript_id')
            if t is None:
                continue
            if r[2] == 'transcript':
                g = gtf_attr(r[8], 'gene_id') or t
                if g in wanted_genes:
                    gene_of[t] = g
            elif r[2] == 'exon':
                exons[t].append((r[0], int(r[3]), int(r[4])))
    out, n_tx = collections.defaultdict(list), collections.Counter()
    for t, g in gene_of.items():
        out[g].extend(exons.get(t, []))
        n_tx[g] += 1
    return out, n_tx


# ------------------------------------------------------------------------------------------------------------ core

def build_loci(fam_of, folds, loci_genes, gene_exons):
    """The clustered loci: records, gene_ids, blocks. Returns (loci dict, locus_of_record, counters)."""
    cnt = collections.Counter()
    records = {m: [m] for m in fam_of}
    for ann, rep in folds.items():
        if rep in fam_of:
            if ann in fam_of:
                raise ValueError(f'loci.tsv folds {key_str(ann)} into {key_str(rep)} but it is itself a cluster member')
            records[rep].append(ann)
            cnt['folded_records'] += 1
    locus_of_record = {}
    loci = {}
    for m, recs in records.items():
        gids, exs = [], []
        for rec in recs:
            if rec in locus_of_record:
                raise ValueError(f'record {key_str(rec)} belongs to two loci')
            locus_of_record[rec] = m
            if rec not in loci_genes:
                raise ValueError(f'{key_str(rec)} has no gene line in loci.gff3 (was the families stage run --from-gtf?)')
            rec_ex = []
            for g in loci_genes[rec]:
                if g not in gene_exons or not gene_exons[g]:
                    raise ValueError(f'gene_id {g} of {key_str(rec)} has no exons in the GTF')
                rec_ex.extend(gene_exons[g])
                gids.append(g)
            if len(loci_genes[rec]) > 1:
                cnt['records_with_key_collision'] += 1
            if any(c != rec[0] for c, _, _ in rec_ex):
                raise ValueError(f'{key_str(rec)} has exons on another contig')
            if (min(s for _, s, _ in rec_ex), max(e for _, _, e in rec_ex)) != (rec[1], rec[2]):
                cnt['records_span_mismatch'] += 1  # the key is not the gene's all-transcript span
            exs.extend(rec_ex)
        blocks = merge_blocks([(s - 1, e) for _, s, e in exs])
        loci[m] = {'records': recs, 'gene_ids': gids, 'starts': [s for s, _ in blocks], 'ends': [e for _, e in blocks]}
    return loci, locus_of_record, cnt


def project_paf(path, loci, locus_of_record, cnt):
    """Stream the PAF; hits[m][block][partner] = genome intervals of m's block joined to an exon base of partner."""
    hits = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))

    def add(lst, lo, hi):
        if lst and lo <= lst[-1][1] and hi >= lst[-1][0]:  # coalesce with the last interval when they touch
            lst[-1] = (min(lo, lst[-1][0]), max(hi, lst[-1][1]))
        else:
            lst.append((lo, hi))

    with open_text(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 12:
                cnt['paf_malformed'] += 1
                continue
            cnt['paf_records'] += 1
            if f[0] == f[5]:
                cnt['paf_self'] += 1
                continue
            qk, tk = parse_key(f[0]), parse_key(f[5])
            lq, lt = locus_of_record.get(qk), locus_of_record.get(tk)
            if lq is None or lt is None:
                cnt['paf_unclustered'] += 1
                continue
            if lq == lt:
                cnt['paf_same_locus'] += 1
                continue
            qs, qe, strand, ts, te = int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8])
            q_off, t_off = qk[1] - 1, tk[1] - 1
            Lq, Lt = loci[lq], loci[lt]
            if not (overlaps_any(Lq['starts'], Lq['ends'], q_off + qs, q_off + qe)
                    and overlaps_any(Lt['starts'], Lt['ends'], t_off + ts, t_off + te)):
                cnt['paf_no_exon_interval'] += 1
                continue
            cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
            if cg is None:
                raise ValueError(f'{path}: record {f[0]} -> {f[5]} has no cg:Z CIGAR (the families PAF is run with -c)')
            cnt['paf_projected'] += 1
            any_col = False
            for ql, qh, qb, tl, th, tb in exon_columns(cg, strand, qs, qe, ts, te, q_off, t_off,
                                                       (Lq['starts'], Lq['ends']), (Lt['starts'], Lt['ends'])):
                add(hits[lq][qb][lt], ql, qh)
                add(hits[lt][tb][lq], tl, th)
                any_col = True
            if any_col:
                cnt['paf_exon_exon'] += 1
    return hits


def classify(fam_of, members, order, loci, hits, strand_of):
    """Rows of blocks.tsv and relations.tsv, plus counters."""
    cnt = collections.Counter()
    rows = []
    rel = collections.defaultdict(lambda: {'members': [], 'n_blocks': 0, 'bp': 0})
    fam_key = lambda fid: order[fid]
    for fid in members:
        for m in members[fid]:
            L = loci[m]
            nb = len(L['starts'])
            has_acc = False
            for bi in range(nb):
                s, e = L['starts'][bi], L['ends'][bi]
                partners = hits.get(m, {}).get(bi, {})
                same = sorted(p for p in partners if fam_of[p] == fid)
                row = {'family_id': fid, 'locus': key_str(m), 'chrom': m[0], 'strand': strand_of.get(m, '.'),
                       'gene_ids': ','.join(L['gene_ids']), 'n_records': len(L['records']), 'block': bi + 1,
                       'n_blocks': nb, 'start': s + 1, 'end': e, 'bp': e - s}
                cnt['blocks'] += 1
                if same:
                    row['class'] = 'core'
                    row['core_bp'] = union_len([iv for p in same for iv in partners[p]])
                    row['core_partners'] = ','.join(f'{key_str(p)}={union_len(partners[p])}' for p in same)
                    row['rel_families'] = row['rel_bp'] = row['rel_partners'] = '.'
                    cnt['core_blocks'] += 1
                    cnt['core_block_bp'] += e - s
                else:
                    has_acc = True
                    other = sorted((p for p in partners if fam_of[p] != fid), key=lambda p: (fam_key(fam_of[p]), p))
                    row['class'] = 'accessory'
                    row['core_bp'] = 0
                    row['core_partners'] = '.'
                    cnt['accessory_blocks'] += 1
                    cnt['accessory_bp'] += e - s
                    if other:
                        fams = sorted({fam_of[p] for p in other}, key=fam_key)
                        row['rel_families'] = ','.join(fams)
                        row['rel_bp'] = union_len([iv for p in other for iv in partners[p]])
                        row['rel_partners'] = ','.join(f'{fam_of[p]}|{key_str(p)}={union_len(partners[p])}'
                                                       for p in other)
                        cnt['accessory_blocks_related'] += 1
                        for f2 in fams:
                            r = rel[(fid, f2)]
                            if not r['members'] or r['members'][-1] != m:
                                r['members'].append(m)
                            r['n_blocks'] += 1
                            r['bp'] += union_len([iv for p in other if fam_of[p] == f2 for iv in partners[p]])
                    else:
                        row['rel_families'] = row['rel_partners'] = '.'
                        row['rel_bp'] = 0
                rows.append(row)
            cnt['loci'] += 1
            cnt['loci_with_accessory' if has_acc else 'loci_all_core'] += 1
    rel_rows = []
    for (f1, f2) in sorted(rel, key=lambda k: (fam_key(k[0]), fam_key(k[1]))):
        r = rel[(f1, f2)]
        rel_rows.append({'family_id': f1, 'related_family': f2, 'n_members': len(r['members']),
                         'n_blocks': r['n_blocks'], 'bp': r['bp'], 'reciprocal': 'yes' if (f2, f1) in rel else 'no',
                         'members': ','.join(key_str(m) for m in r['members'])})
    cnt['family_relations_directed'] = len(rel_rows)
    cnt['family_relations_reciprocal'] = sum(1 for r in rel_rows if r['reciprocal'] == 'yes')
    cnt['families'] = len(members)
    cnt['families_with_relation'] = len({f1 for f1, _ in rel})
    return rows, rel_rows, cnt


def run(gtf, clusters, loci_gff3, loci_tsv, paf):
    """The whole container on in-memory results: (block rows, relation rows, summary Counter)."""
    fam_of, members, order = read_clusters(clusters)
    folds = read_folds(loci_tsv) if loci_tsv else {}
    loci_genes, strand_of = read_loci_gff3(loci_gff3)
    wanted = set()
    for m in fam_of:
        wanted.update(loci_genes.get(m, []))
    for ann, rep in folds.items():
        if rep in fam_of:
            wanted.update(loci_genes.get(ann, []))
    gene_exons, _ = read_gtf(gtf, wanted)
    loci, locus_of_record, cnt = build_loci(fam_of, folds, loci_genes, gene_exons)
    cnt['loci_gff3_keys'] = len(loci_genes)
    cnt['gene_key_collisions'] = sum(1 for v in loci_genes.values() if len(v) > 1)
    hits = project_paf(paf, loci, locus_of_record, cnt)
    rows, rel_rows, c2 = classify(fam_of, members, order, loci, hits, strand_of)
    cnt.update(c2)
    return rows, rel_rows, cnt


SUMMARY_KEYS = ['families', 'loci', 'folded_records', 'loci_gff3_keys', 'gene_key_collisions',
                'records_with_key_collision', 'records_span_mismatch', 'blocks', 'core_blocks', 'accessory_blocks',
                'core_block_bp', 'accessory_bp', 'accessory_blocks_related', 'loci_all_core', 'loci_with_accessory',
                'families_with_relation', 'family_relations_directed', 'family_relations_reciprocal', 'paf_records',
                'paf_malformed', 'paf_self', 'paf_unclustered', 'paf_same_locus', 'paf_no_exon_interval',
                'paf_projected', 'paf_exon_exon']


def write_tsv(path, header, rows):
    with open(path, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[h]) for h in header) + '\n')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument('--gtf', required=True, help='the assembled GTF the families stage ran on (PREFIX.gtf)')
    ap.add_argument('--fam', help='the families prefix (PREFIX.fam); fills the four paths below')
    ap.add_argument('--clusters')
    ap.add_argument('--loci-gff3')
    ap.add_argument('--loci-tsv', help='the fold table; a missing file under --fam means no folded records')
    ap.add_argument('--paf')
    ap.add_argument('--out', required=True, help='output prefix')
    a = ap.parse_args(argv)
    fam = a.fam
    clusters = a.clusters or (fam and fam + '.clusters.tsv')
    loci_gff3 = a.loci_gff3 or (fam and fam + '.loci.gff3')
    paf = a.paf or (fam and fam + '.loci.paf')
    loci_tsv = a.loci_tsv
    if loci_tsv is None and fam:
        loci_tsv = fam + '.loci.tsv'
        if not os.path.exists(loci_tsv):
            print(f'[family_container] no {loci_tsv}: no folded records', file=sys.stderr)
            loci_tsv = None
    if not (clusters and loci_gff3 and paf):
        ap.error('give --fam, or all of --clusters --loci-gff3 --paf')
    rows, rel_rows, cnt = run(a.gtf, clusters, loci_gff3, loci_tsv, paf)
    write_tsv(a.out + '.blocks.tsv', BLOCK_HEADER, rows)
    write_tsv(a.out + '.relations.tsv', REL_HEADER, rel_rows)
    with open(a.out + '.summary.tsv', 'w') as fh:
        fh.write('key\tvalue\n')
        for k in SUMMARY_KEYS:
            fh.write(f'{k}\t{cnt.get(k, 0)}\n')
    print('[family_container] ' + ' '.join(f'{k}={cnt.get(k, 0)}' for k in SUMMARY_KEYS), file=sys.stderr)
    print(f'[family_container] wrote {a.out}.blocks.tsv {a.out}.relations.tsv {a.out}.summary.tsv', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
