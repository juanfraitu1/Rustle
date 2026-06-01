#!/usr/bin/env python3
"""Enrich the secondary_dependent hits from run_scan with the signals needed
for adversarial verification:
  - n_unique   : chain-supporting reads whose PRIMARY is at the gene AND that
                 have NO secondary alignment anywhere (unambiguously this
                 gene's reads -> strongest 'genuinely expressed' anchor)
  - gene/sibling biotype, name, description (from gene_meta.tsv)
  - same_family: gene Name and sibling Name share a family stem
  - locus_rel  : trans (diff contig) | overlapping | tandem (<1Mb) | distal

Usage: enrich_hits.py SCAN.tsv GENE_INTRONS.tsv GENE_META.tsv MULTINAMES.txt BAM OUT.json
"""
import sys, re, subprocess, json, collections

def junctions(pos, cig):
    cur = pos; js = []
    for ln, op in re.findall(r'(\d+)([MIDNSH=X])', cig):
        ln = int(ln)
        if op == 'N':
            js.append((cur - 1, cur + ln)); cur += ln
        elif op in 'MD=X':
            cur += ln
    return js

def stem(name):
    if not name or name == 'None':
        return ''
    s = re.sub(r'[-_]?\d.*$', '', name)   # strip trailing digits/suffix
    return s

def main():
    scan, gi, meta, mm, bam, out = sys.argv[1:7]
    genes = {}
    for line in open(gi):
        if line.startswith('gene_id'):
            continue
        gid, chrom, s, e, strand, ntx, chain = line.rstrip('\n').split('\t')
        ints = set()
        if chain:
            for pr in chain.split(','):
                d, a = pr.split(':'); ints.add((int(d), int(a)))
        genes[gid] = dict(chrom=chrom, start=int(s), end=int(e),
                          strand=strand, introns=ints)
    M = {}
    for line in open(meta):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 4:
            M[f[0]] = dict(biotype=f[1], name=f[2], desc=f[3])
        elif len(f) >= 1:
            M[f[0]] = dict(biotype=f[1] if len(f) > 1 else '',
                           name=f[2] if len(f) > 2 else '', desc='')
    mmset = set(l.rstrip('\n') for l in open(mm))
    print(f"loaded {len(mmset)} multimapper names", file=sys.stderr)
    rows = []
    header = None
    for line in open(scan):
        f = line.rstrip('\n').split('\t')
        if header is None:
            header = f; continue
        r = dict(zip(header, f))
        if r['klass'] not in ('secondary_dependent_ambiguous',
                              'secondary_dependent_identifiable'):
            continue
        rows.append(r)
    print(f"enriching {len(rows)} secondary_dependent hits", file=sys.stderr)
    out_rows = []
    for r in rows:
        gid = r['gene']; g = genes[gid]; K = int(r['K'])
        region = f"{g['chrom']}:{g['start']}-{g['end']}"
        try:
            txt = subprocess.run(['samtools', 'view', bam, region],
                                 capture_output=True, text=True,
                                 timeout=120).stdout
        except Exception:
            txt = ''
        n_unique = 0
        seen = set()
        for line in txt.splitlines():
            ff = line.split('\t')
            if len(ff) < 11:
                continue
            flag = int(ff[1])
            if flag & 0x900:                 # secondary or supplementary
                continue
            if ('-' if flag & 16 else '+') != g['strand']:
                continue
            js = junctions(int(ff[3]), ff[5])
            if sum(1 for j in js if j in g['introns']) < K:
                continue
            qn = ff[0]
            if qn in seen:
                continue
            seen.add(qn)
            if qn not in mmset:              # primary here AND no secondary anywhere
                n_unique += 1
        sib = r['dom_sibling']
        gm = M.get(gid, {}); sm = M.get(sib, {})
        gname = gm.get('name', ''); sname = sm.get('name', '')
        same_family = bool(stem(gname) and stem(gname) == stem(sname))
        # locus relationship
        locus_rel = 'no_sibling'
        if sib in genes:
            sg = genes[sib]
            if sg['chrom'] != g['chrom']:
                locus_rel = 'trans'
            elif not (g['end'] < sg['start'] or sg['end'] < g['start']):
                locus_rel = 'overlapping'
            else:
                d = min(abs(g['start'] - sg['end']), abs(sg['start'] - g['end']))
                locus_rel = 'tandem' if d < 1_000_000 else 'distal'
        o = dict(r)
        o.update(n_unique=n_unique, gene_name=gname, gene_desc=gm.get('desc', ''),
                 gene_biotype=gm.get('biotype', ''),
                 sibling_name=sname, sibling_biotype=sm.get('biotype', ''),
                 sibling_desc=sm.get('desc', ''),
                 same_family=same_family, locus_rel=locus_rel)
        out_rows.append(o)
    json.dump(out_rows, open(out, 'w'), indent=1)
    # quick deterministic provisional taxonomy
    tax = collections.Counter()
    for o in out_rows:
        nu = o['n_unique']; pg = int(o['primary_at_gene'])
        if nu >= 3 or pg >= 5:
            t = 'expressed_copy'           # genuine independent expression
        elif nu == 0 and pg == 0:
            t = 'pure_spillover'           # fabrication risk
        else:
            t = 'weak_expression'
        o['provisional_tax'] = t
        tax[(o['klass'], t)] += 1
    json.dump(out_rows, open(out, 'w'), indent=1)
    print("provisional taxonomy (klass, tax):", file=sys.stderr)
    for k, v in sorted(tax.items()):
        print(f"   {k}: {v}", file=sys.stderr)

if __name__ == '__main__':
    main()
