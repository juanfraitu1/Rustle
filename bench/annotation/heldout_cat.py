#!/usr/bin/env python3
"""Step 3 of docs/CAT_RERUN_PROTOCOL_2026-10-01.md: the held-out family test (docs/PREREG_heldout_families_2026-09-20.md,
Test 2) scored with the Soto S1C join by `Gene ID` (rule R5), so it can be run on CAT/Liftoff v2.0.

`score.py heldout --soto` (the registered scorer) keys genes on the GFF `Name=` and joins Soto by `Gene Name`. Under the
CAT slim GFF `Name=` is the CAT gene id, which is the S1C `Gene ID`, so the only change needed is the join key. This
script reuses `score.heldout_load_genes`, `score.heldout_predicted_clusters`, `lib.families_on` and
`lib.bipartite_families` unchanged; score.py and lib.py are not edited.

    score     one scorer, two joins:
                --join name : the registered path (genes by `Name=`, Soto by `Gene Name`, lib.soto_gene_family);
                              on the original RefSeq outputs it must reproduce the registered numbers.
                --join id   : Soto by `Gene ID` (soto_gene_family_by_id: same exclusion rule, a Gene ID carrying more
                              than one distinct Family ID is excluded; `N/A` / empty dropped).
              Per chromosome: score.py heldout's summary line and JSON (same keys). Per arm: the pooled numbers of
              docs/HELDOUT_FAMILIES_RESULT_2026-09-20.md (mean over all the arm's truth families of F / sens / prec;
              an unmatched family scores 0 and is kept). A .gz GFF is read with the same logic as
              score.heldout_load_genes.
    universe  which Soto genes / families each join puts in each chromosome's truth, and why they differ
              (RefSeq lacks the symbol on that chromosome; a Gene Name shared by several Gene IDs; a name carrying
              more than one Family ID; a RefSeq symbol on this chromosome whose Soto gene lies elsewhere).

    python3 bench/annotation/heldout_cat.py score --gff G --soto bench/soto/soto_famCN_S1C.tsv --join id \\
        --clusters chr2=D/chr2_fam.clusters.tsv ... --arm heldout=chr2,chr8,chr10 --arm development=chr5,chr7,chr21 \\
        --json-dir D --families-tsv D/families.tsv
"""
import argparse
import collections
import csv
import gzip
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib  # noqa: E402
import score  # noqa: E402


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def load_genes_all(gff):
    """{contig: {(start1, end): Name}} for gene/pseudogene records. A plain file goes through
    score.heldout_load_genes(gff, 'ALL') unchanged; a .gz file is read with the same rule (last record wins on a
    repeated (contig, start, end), as in score.py)."""
    if gff.endswith('.gz'):
        flat = {}
        name_re = re.compile(r'Name=([^;]+)')
        with _open(gff) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 9 or f[2] not in ('gene', 'pseudogene'):
                    continue
                m = name_re.search(f[8])
                if m:
                    flat[(f[0], int(f[3]), int(f[4]))] = m.group(1)
    else:
        flat = score.heldout_load_genes(gff, 'ALL')
    per = collections.defaultdict(dict)
    for (c, s, e), nm in flat.items():
        per[c][(s, e)] = nm
    return per


def soto_gene_family_by_id(s1c):
    """S1C -> {Gene ID: Family ID}; lib.soto_gene_family's rule with `Gene ID` as the key: a Gene ID carrying more than
    one distinct Family ID is EXCLUDED; `N/A` and empty IDs are dropped."""
    ids = collections.defaultdict(set)
    for r in csv.DictReader(open(s1c), delimiter='\t'):
        fid = (r.get('Family ID') or '').strip(); gid = (r.get('Gene ID') or '').strip()
        if fid and fid != 'N/A' and gid:
            ids[gid].add(fid)
    return {g: next(iter(f)) for g, f in ids.items() if len(f) == 1}


def label_map(s1c, join):
    return lib.soto_gene_family(s1c) if join == 'name' else soto_gene_family_by_id(s1c)


def translated_clusters(clusters_tsv, genes, tmap, exact_only=False):
    """score.heldout_predicted_clusters with each member re-keyed to its R1 CAT gene: the member is resolved to a
    record of `genes` by the same probe ((start+1, end) then (start, end); (start, end) only with exact_only), and
    that record's (start, end) is looked up in `tmap` (the refseq_map rows of the chromosome, quality != none). A
    record with no CAT gene stays 'refseq:SYMBOL' (it can only cost precision); an unresolved row stays chrom:s-e."""
    out = collections.defaultdict(list)
    with open(clusters_tsv) as fh:
        for line in fh:
            if line.startswith('cluster_id'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 8:
                continue
            cid, s, e = p[0], int(p[6]), int(p[7])
            key = (s, e) if exact_only else ((s + 1, e) if (s + 1, e) in genes else (s, e))
            if key not in genes:
                out[cid].append(f'{p[5]}:{s}-{e}')
            else:
                out[cid].append(tmap.get(key) or f'refseq:{genes[key]}')
    return dict(out)


def score_chrom(chrom, genes, label_of, clusters, exact_only, only=None, truth_genes=None, tmap=None):
    """score.cmd_heldout's --soto body for one chromosome. Returns (summary, per_family, truth, pred) or None.
    only: a set of family ids; the truth is restricted to them BEFORE the bipartite matching (comparison on a fixed
    family set, e.g. the families the registered RefSeq run scored). truth_genes: the genes whose names define the
    truth universe (default `genes`). tmap: re-key members to CAT genes (translated_clusters)."""
    import numpy as np
    truth = lib.families_on(label_of, set((truth_genes if truth_genes is not None else genes).values()), 3)
    if only is not None:
        truth = {f: m for f, m in truth.items() if f in only}
    if tmap is None:
        pred = score.heldout_predicted_clusters(clusters, genes, exact_only, False)
    else:
        pred = translated_clusters(clusters, genes, tmap, exact_only)
    per = lib.bipartite_families(truth, pred)
    if per is None:
        return None
    fs = [v['f'] for v in per.values()]
    summary = dict(chrom=chrom, truth_families=len(truth), pred_clusters=len(pred),
                   truth_families_touched=sum(1 for v in per.values() if v['cluster']),
                   mean_F=round(float(np.mean(fs)), 4),
                   mean_sens=round(float(np.mean([v['sens'] for v in per.values()])), 4),
                   mean_prec=round(float(np.mean([v['prec'] for v in per.values()])), 4),
                   exact_recoveries=sum(1 for v in per.values() if v['f'] == 1.0))
    return summary, per, truth, pred


def cmd_score(a):
    import numpy as np
    clusters = dict(x.split('=', 1) for x in a.clusters)
    arms = [(x.split('=', 1)[0], x.split('=', 1)[1].split(',')) for x in a.arm]
    genes_all = load_genes_all(a.gff)
    label_of = label_map(a.soto, a.join)
    names = {}
    if a.names:
        for r in csv.DictReader(open(a.names), delimiter='\t'):
            names[r['gene_id']] = r['gene_name']
    only = None
    if a.only_families:
        only = collections.defaultdict(set)
        for ln in open(a.only_families):
            f = ln.split()
            if len(f) >= 2 and not ln.startswith('#'):
                only[f[0]].add(f[1])
    truth_all = load_genes_all(a.truth_gff) if a.truth_gff else None
    tmaps = None
    if a.translate:
        tmaps = collections.defaultdict(dict)
        for r in csv.DictReader(open(a.translate), delimiter='\t'):
            if r['quality'] != 'none' and r['cat_gene']:
                tmaps[r['chrom']][(int(r['start0']) + 1, int(r['end']))] = r['cat_gene']
    rows, by_chrom = [], {}
    for arm, chroms in arms:
        for c in chroms:
            res = score_chrom(c, genes_all.get(c, {}), label_of, clusters[c], a.exact_only,
                              None if only is None else only.get(c, set()),
                              None if truth_all is None else truth_all.get(c, {}),
                              None if tmaps is None else tmaps.get(c, {}))
            if res is None:
                print(f'{c}: NO TRUTH FAMILIES (>=3 members) — chromosome not scoreable')
                by_chrom[c] = None
                continue
            s, per, truth, pred = res
            by_chrom[c] = (s, per)
            print(f"{c}: truth families {s['truth_families']} | predicted clusters {s['pred_clusters']} | "
                  f"touched {s['truth_families_touched']} | mean F {s['mean_F']} "
                  f"(sens {s['mean_sens']} / prec {s['mean_prec']}) | exact {s['exact_recoveries']}")
            if a.json_dir:
                with open(os.path.join(a.json_dir, f'{c}_soto{a.json_suffix}.json'), 'w') as fh:
                    json.dump(dict(summary=s, per_family=per), fh, indent=1, sort_keys=True)
            for fam, v in sorted(per.items()):
                tm = sorted(truth[fam]); cm = sorted(pred[v['cluster']]) if v['cluster'] else []
                nm = lambda xs: ','.join(f'{x}({names[x]})' if x in names else x for x in xs)
                rows.append([arm, c, fam, v['n_truth'], v['n_pred'], v['hit'], v['sens'], v['prec'], v['f'],
                             v['cluster'] or '', nm(tm), nm(cm)])
    print()
    for arm, chroms in arms:
        per_all = [v for c in chroms if by_chrom.get(c) for v in by_chrom[c][1].values()]
        if not per_all:
            print(f'{arm}: no scoreable family'); continue
        fs = [v['f'] for v in per_all]
        print(f"{arm} ({','.join(chroms)}): families {len(per_all)} | touched {sum(1 for v in per_all if v['cluster'])} | "
              f"pooled F {np.mean(fs):.4f} (sens {np.mean([v['sens'] for v in per_all]):.4f} / "
              f"prec {np.mean([v['prec'] for v in per_all]):.4f}) | exact {sum(1 for f in fs if f == 1.0)}")
    if a.families_tsv:
        with open(a.families_tsv, 'w') as fh:
            fh.write('arm\tchrom\tfamily\tn_truth\tn_pred\thit\tsens\tprec\tf\tcluster\ttruth_members\tcluster_members\n')
            for r in rows:
                fh.write('\t'.join(map(str, r)) + '\n')


def cmd_universe(a):
    """Per chromosome and Soto family: the name-join truth (RefSeq GFF names) beside the id-join truth (CAT GFF ids),
    and one row per Soto gene saying where it lands under each join and why."""
    chroms = a.chroms.split(',')
    ref = load_genes_all(a.refseq_gff)
    cat = load_genes_all(a.cat_gff)
    by_name, by_id = lib.soto_gene_family(a.soto), soto_gene_family_by_id(a.soto)
    s1c = list(csv.DictReader(open(a.soto), delimiter='\t'))
    id2name = {r['Gene ID'].strip(): r['Gene Name'].strip() for r in s1c}
    name2ids = collections.defaultdict(set)
    for g, n in id2name.items():
        name2ids[n].add(g)
    cat_chrom = {nm: c for c, d in cat.items() for nm in d.values()}       # CAT id -> contig
    out_g = open(a.out + '.genes.tsv', 'w')
    out_g.write('chrom\tgene_id\tgene_name\tid_family\tname_family\tin_id_truth_pool\tin_name_truth_pool\tnote\n')
    out_f = open(a.out + '.families.tsv', 'w')
    out_f.write('chrom\tfamily\tn_id_join\tn_name_join\tscored_cat\tscored_refseq\n')
    summary = []
    for c in chroms:
        ids_here = set(cat.get(c, {}).values()); names_here = set(ref.get(c, {}).values())
        t_id = lib.families_on(by_id, ids_here, 1); t_nm = lib.families_on(by_name, names_here, 1)
        fams = sorted(set(t_id) | set(t_nm), key=lambda f: (int(f.split('_')[1]) if f.split('_')[1].isdigit() else 0, f))
        for f in fams:
            ni, nn = len(t_id.get(f, [])), len(t_nm.get(f, []))
            if ni >= 3 or nn >= 3:
                out_f.write(f'{c}\t{f}\t{ni}\t{nn}\t{"yes" if ni >= 3 else "no"}\t{"yes" if nn >= 3 else "no"}\n')
        # gene rows: every Soto Gene ID whose CAT gene is on c, plus every Soto name present on c in RefSeq
        seen = set()
        for g in sorted(id2name):
            nm = id2name[g]
            on_cat = cat_chrom.get(g) == c
            on_ref = nm in names_here
            if not (on_cat or on_ref):
                continue
            seen.add(g)
            fi, fn = by_id.get(g, ''), by_name.get(nm, '')
            notes = []
            if g not in cat_chrom:
                notes.append('gene_id_absent_from_CAT')
            elif not on_cat:
                notes.append(f'CAT_gene_on_{cat_chrom[g]}_but_RefSeq_has_symbol_on_{c}')
            if on_cat and not on_ref:
                notes.append('symbol_absent_from_RefSeq_on_this_chrom')
            if not fi:
                notes.append('gene_id_multi_family_excluded')
            if fi and not fn:
                notes.append('gene_name_multi_family_excluded')
            if len(name2ids[nm]) > 1:
                notes.append(f'gene_name_shared_by_{len(name2ids[nm])}_gene_ids')
            out_g.write(f'{c}\t{g}\t{nm}\t{fi}\t{fn}\t{"yes" if on_cat and fi else "no"}\t'
                        f'{"yes" if on_ref and fn else "no"}\t{";".join(notes) or "-"}\n')
        # summary over the SCORED truth (families with >= 3 members under that join)
        sc_id = {f: m for f, m in t_id.items() if len(m) >= 3}
        sc_nm = {f: m for f, m in t_nm.items() if len(m) >= 3}
        why = collections.Counter()
        named = set()
        for f, members in sc_id.items():
            for g in members:
                n = id2name[g]
                if n not in names_here:
                    why['symbol_absent_from_RefSeq_here'] += 1
                elif n not in by_name:
                    why['name_multi_family_excluded'] += 1
                elif (f, n) in named:
                    why['name_shared_with_another_id_collapsed'] += 1
                elif f not in sc_nm:
                    named.add((f, n)); why['family_below_3_by_name'] += 1
                else:
                    named.add((f, n)); why['shared'] += 1
        ref_only = 0
        for f, members in sc_nm.items():
            for n in members:
                if not any(cat_chrom.get(g) == c for g in name2ids[n]):
                    ref_only += 1
        summary.append((c, len(sc_nm), len(sc_id), len(set(sc_nm) & set(sc_id)),
                        sum(map(len, sc_nm.values())), sum(map(len, sc_id.values())), why, ref_only))
    out_g.close(); out_f.close()
    keys = ('shared', 'symbol_absent_from_RefSeq_here', 'name_shared_with_another_id_collapsed',
            'name_multi_family_excluded', 'family_below_3_by_name')
    with open(a.out + '.summary.tsv', 'w') as fh:
        fh.write('chrom\tfamilies_name_join\tfamilies_id_join\tfamilies_both\tmembers_name_join\tmembers_id_join\t'
                 + '\t'.join('id_members_' + k for k in keys) + '\tname_members_whose_soto_gene_is_elsewhere\n')
        for c, fn, fi, fb, mn, mi, why, ro in summary:
            fh.write(f'{c}\t{fn}\t{fi}\t{fb}\t{mn}\t{mi}\t' + '\t'.join(str(why[k]) for k in keys) + f'\t{ro}\n')
    print(open(a.out + '.summary.tsv').read())
    print(f'wrote {a.out}.genes.tsv + {a.out}.families.tsv + {a.out}.summary.tsv')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('score')
    p.add_argument('--gff', required=True, help='GFF (plain or .gz) whose gene/pseudogene Name= keys the members')
    p.add_argument('--soto', required=True)
    p.add_argument('--join', choices=('name', 'id'), required=True)
    p.add_argument('--clusters', action='append', required=True, metavar='CHROM=CLUSTERS_TSV')
    p.add_argument('--arm', action='append', required=True, metavar='ARM=chrA,chrB')
    p.add_argument('--exact-only', action='store_true', help="score.py heldout's --exact-only (B4)")
    p.add_argument('--json-dir'); p.add_argument('--json-suffix', default='')
    p.add_argument('--families-tsv')
    p.add_argument('--names', help='CAT genes.tsv: print CAT ids as id(gene_name) in --families-tsv')
    p.add_argument('--only-families', help='file of "chrom family" lines: restrict each chromosome\'s truth to these '
                   'families before matching (a chromosome with none listed is not scoreable)')
    p.add_argument('--truth-gff', help='GFF whose gene Name= set defines the truth universe (default --gff); with '
                   '--translate: the CAT GFF, while --gff is the RefSeq GFF the clusters were built on')
    p.add_argument('--translate', help='chm13v2.0_CAT_Liftoff.refseq_map.tsv: re-key RefSeq cluster members to their '
                   'R1 CAT gene (decomposition: RefSeq clusters scored against the CAT truth)')
    p.set_defaults(func=cmd_score)
    p = sub.add_parser('universe')
    for x in ('--refseq-gff', '--cat-gff', '--soto', '--chroms', '--out'):
        p.add_argument(x, required=True)
    p.set_defaults(func=cmd_universe)
    a = ap.parse_args(argv)
    a.func(a)


if __name__ == '__main__':
    main()
