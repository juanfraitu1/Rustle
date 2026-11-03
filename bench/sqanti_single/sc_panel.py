#!/usr/bin/env python3
"""PREREG SQANTI SINGLE-COPY (register 1482, docs/PREREG_sqanti_single_copy_2026-11-01.md).
  sc_panel.py select OUR.gtf COPIES.tsv BAM CHROM OUT_DIR   -> OUT_DIR/targets.tsv (isocon_targets.py's columns, class single)
  sc_panel.py cut    OUT_DIR/targets.tsv IN.gtf[.gz] OUT.gtf  -> transcripts overlapping a target span (same strand or '.')
  sc_panel.py table  OUT_DIR/targets.tsv [--expressed REF.gtf:BAM[:K]] ARM=classification.txt[:filter.txt] ...  -> the SQANTI summary table
      --expressed adds expr_recall: FSM-matched reference isoforms / reference multi-exon isoforms whose exact intron chain is
      carried by >= K (default 2) primary reads on the target spans (cached in OUT_DIR/expressed.kK.txt), and a footer with the
      expressed isoforms no arm matches.
Universe: CAT protein-coding, Ensembl id, no Compara paralog at any level, expressed in OUR.gtf, overlapping no copy of COPIES.tsv,
50-5,000 BAM records on the span."""
import sys, os, gzip, csv, collections, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import famtruth as F
cmd = sys.argv[1]
def sh(c): return subprocess.run(c, shell=True, capture_output=True, text=True).stdout
def targets(path):
    return [r for r in csv.DictReader(open(path), delimiter='\t')]
if cmd == 'select':
    gtf, copies, bam, chrom, out = sys.argv[2:7]; os.makedirs(out, exist_ok=True)
    G = F.load_cat(); idx = F.exon_index(G)
    levels = set()
    for v in F.LEVELS.values(): levels |= set(v)
    pairs, universe = F.truth_pairs(G, F.PARALOGS, levels)
    para = {g for p in pairs for g in p}
    expressed = F.expressed_from_gtf(gtf, idx)
    cp = [(int(r['locus_start']), int(r['locus_end'])) for r in csv.DictReader(open(copies), delimiter='\t') if r['chrom'] == chrom]
    n_all = n_pc = n_univ = n_para = n_expr = n_copy = n_depth = 0; rows = []
    for g, d in sorted(G.items()):
        if d['chrom'] != chrom or not d['exons']: continue
        n_all += 1
        if d['biotype'] != 'protein_coding': continue
        n_pc += 1
        if g not in universe or not d['ens']: continue
        n_univ += 1
        if g in para: n_para += 1; continue
        if g not in expressed: continue
        n_expr += 1
        s, e = d['exons'][0][0], d['exons'][-1][1]
        if any(a < e and s < b for a, b in cp): n_copy += 1; continue
        n = int(sh(f"samtools view -c {bam} {chrom}:{s + 1}-{e}").strip() or 0)
        if not 50 <= n <= 5000: n_depth += 1; continue
        rows.append((d['name'], g, s, e, d['strand'], n))
    with open(f"{out}/targets.tsv", 'w') as w:
        w.write("target\tclass\tgene\tname\tchrom\tstart0\tend\tstrand\trecords\n")
        for name, g, s, e, st, n in rows: w.write(f"{name}\tsingle\t{g}\t{name}\t{chrom}\t{s}\t{e}\t{st}\t{n}\n")
    print(f"{chrom}: genes {n_all}, protein_coding {n_pc}, with Ensembl id in the Compara universe {n_univ}, dropped: paralog {n_para}, "
          f"not expressed {n_univ - n_para - n_expr}, overlaps a family copy {n_copy}, outside 50-5000 records {n_depth}; universe {len(rows)}")
elif cmd == 'cut':
    tpath, gin, gout = sys.argv[2:5]
    T = targets(tpath); spans = collections.defaultdict(list)
    for r in T: spans[r['chrom']].append((int(r['start0']), int(r['end']), r['strand']))
    op = gzip.open if gin.endswith('.gz') else open
    tx = collections.defaultdict(list)
    for line in op(gin, 'rt'):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] not in spans: continue
        i = f[8].find('transcript_id "'); tid = f[8][i + 15:f[8].find('"', i + 15)]
        tx[tid].append(f)
    kept = 0
    with open(gout, 'w') as w:
        for tid, fs in tx.items():
            ex = [f for f in fs if f[2] == 'exon'] or fs
            s, e, st, c = min(int(f[3]) for f in ex) - 1, max(int(f[4]) for f in ex), ex[0][6], ex[0][0]
            if any(a < e and s < b and (st == '.' or st == t) for a, b, t in spans[c]):
                kept += 1
                for f in fs: w.write('\t'.join(f) + '\n')
    print(f"{gin}: {kept} of {len(tx)} transcripts on the universe spans -> {gout}")
elif cmd == 'table':
    import re
    T = targets(sys.argv[2]); spans = [(r['chrom'], int(r['start0']), int(r['end']), r['strand'], r['gene'], r['name']) for r in T]
    ART = {'genic', 'antisense', 'intergenic', 'genic_intron', 'fusion'}
    args = sys.argv[3:]; expressed = None
    if args and args[0] == '--expressed':
        ref, bam, *k = args[1].split(':'); k = int(k[0]) if k else 2; args = args[2:]
        cache = os.path.join(os.path.dirname(os.path.abspath(sys.argv[2])), f'expressed.k{k}.txt')
        if not os.path.exists(cache):
            import pysam
            ex = collections.defaultdict(list)
            for l in open(ref):
                f = l.rstrip('\n').split('\t')
                if len(f) > 8 and f[2] == 'exon':
                    ex[re.search(r'transcript_id "([^"]+)"', f[8]).group(1)].append((f[0], int(f[3]) - 1, int(f[4])))
            chain = {}
            for t, e in ex.items():
                e.sort()
                if len(e) > 1: chain[t] = (e[0][0], tuple((a[2], b[1]) for a, b in zip(e, e[1:])))
            cnt = collections.Counter(); seen = set(); B = pysam.AlignmentFile(bam)
            for c, a, b, st, g, name in spans:
                for r in B.fetch(c, a, b):
                    if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_name in seen: continue
                    seen.add(r.query_name); pos = r.reference_start; ints = []
                    for op, l in r.cigartuples:
                        if op == 3: ints.append((pos, pos + l))
                        if op in (0, 2, 3, 7, 8): pos += l
                    if ints: cnt[(c, tuple(ints))] += 1
            with open(cache, 'w') as w:
                for t, ch in chain.items():
                    if cnt[ch] >= k: w.write(t + '\n')
        expressed = set(l.strip() for l in open(cache)); matched_any = set()
    print("arm\ttranscripts\tgenes_modeled/" + str(len(T)) + "\tFSM\tISM\tNIC\tNNC\tartifact\tFSM+ISM%\tartifact%\trules_PASS%\ttx_per_gene"
          + (f"\texpr_recall/{len(expressed)}" if expressed is not None else ""))
    for spec in args:
        arm, files = spec.split('='); cls = files.split(':')[0]; flt = files.split(':')[1] if ':' in files else None
        C = list(csv.DictReader(open(cls), delimiter='\t'))
        cat = collections.Counter(); hit = set()
        for r in C:
            cat['artifact' if r['structural_category'] in ART else r['structural_category']] += 1
            s, e = int(r['start' if 'start' in r else 'tts']), int(r['end' if 'end' in r else 'tss'])
            lo, hi = min(s, e), max(s, e)
            for c, a, b, st, g, name in spans:
                if r['chrom'] == c and a < hi and lo < b: hit.add(g)
        n = len(C); fi = cat['full-splice_match'] + cat['incomplete-splice_match']
        ps = 'NA'
        if flt:
            Fr = list(csv.DictReader(open(flt), delimiter='\t'))
            ps = f"{100 * sum(1 for r in Fr if r.get('filter_result') == 'Isoform') / max(len(Fr), 1):.1f}"
        print(f"{arm}\t{n}\t{len(hit)}\t{cat['full-splice_match']}\t{cat['incomplete-splice_match']}\t{cat['novel_in_catalog']}\t"
              f"{cat['novel_not_in_catalog']}\t{cat['artifact']}\t{100 * fi / max(n, 1):.1f}\t{100 * cat['artifact'] / max(n, 1):.1f}\t{ps}\t{n / max(len(hit), 1):.2f}"
              + (f"\t{len(m)} ({100 * len(m) / max(len(expressed), 1):.1f}%)" if expressed is not None and (m := {r['associated_transcript'] for r in C if r['structural_category'] == 'full-splice_match'} & expressed) is not None else ""))
        if expressed is not None: matched_any |= m
    if expressed is not None:
        print(f"# expressed reference isoforms matched by no arm: {len(expressed - matched_any)} of {len(expressed)}")
