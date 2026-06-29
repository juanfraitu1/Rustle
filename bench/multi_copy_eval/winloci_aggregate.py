#!/usr/bin/env python3
"""Aggregate the authoritative sweep into the final portfolio.
A GENUINE VG win transcript must be ALL of:
  (1) in VG output (vg_eq) and NOT in rustle primary-only baseline (vg_gain_vs_base)  -> VG/multimapping-specific
  (2) annotated at the CANDIDATE (secondary-dependent) locus                          -> the copy ST structurally misses
  (3) NOT matched by whole-genome StringTie                                            -> not a region artifact
  (4) the candidate copy is read-backed (verify verdict == 'real')                     -> not a phantom
Usage: winloci_aggregate.py SWEEPDIR candidates.json genomeST.tmap GGO_genomic.gff
"""
import json, re, glob, sys, os
sweep, cand_f, gst_tmap, gff = sys.argv[1:5]
cands = {c['gene_id']: c for c in json.load(open(cand_f))}
stg = set()
for l in open(gst_tmap):
    f = l.rstrip('\n').split('\t')
    if len(f) > 2 and f[2] == '=': stg.add(f[1])
def st_has(t):
    b = t.replace('rna-', '').replace('gene-', '')
    return ('rna-'+b in stg) or (b in stg) or (t in stg)

# load evals
ev = {}
for p in glob.glob(f'{sweep}/*.eval.json'):
    g = os.path.basename(p)[:-len('.eval.json')]
    try: ev[g] = json.load(open(p))
    except Exception: pass
ver = {}
for p in glob.glob(f'{sweep}/*.verify.json'):
    g = os.path.basename(p)[:-len('.verify.json')]
    try: ver[g] = json.load(open(p))
    except Exception: pass

# tx coords for at-candidate test
want = set()
for g, e in ev.items():
    for t in e.get('vg_gain_vs_st', []): want.add(t.replace('rna-', '').replace('gene-', ''))
txloc = {}
for l in open(gff):
    if l[0] == '#': continue
    f = l.split('\t')
    if len(f) < 9: continue
    m = re.search(r'ID=(?:rna-|gene-)?([^;]+)', f[8])
    if m and m.group(1) in want and f[2] in ('mRNA', 'transcript', 'lnc_RNA', 'ncRNA'):
        txloc[m.group(1)] = (f[0], int(f[3]), int(f[4]))
def at_cand(g, t):
    c = cands.get(g); loc = txloc.get(t.replace('rna-', '').replace('gene-', ''))
    return bool(c and loc and loc[0] == c['chrom'] and loc[1] <= c['end'] and loc[2] >= c['start'])

from collections import Counter
cls = Counter(e.get('classification', 'NA') for e in ev.values())
print(f"=== SWEEP: {len(ev)} evals ===")
print("class distribution:", dict(cls))

portfolio = []
for g, e in ev.items():
    if e.get('classification') != 'win_vs_st': continue
    gainbase = set(e.get('vg_gain_vs_base', []))   # VG-specific over rustle baseline
    genuine = [t for t in e.get('vg_gain_vs_st', [])
               if t in gainbase and at_cand(g, t) and not st_has(t)]
    v = ver.get(g, {})
    readbacked = v.get('verdict') == 'real'
    if genuine and readbacked:
        portfolio.append(dict(gene=g, region=e.get('region'), genuine=genuine,
                              owner_frac=v.get('owner_frac', 0), pure=v.get('n_pure_owner', 0),
                              strict=v.get('n_strict_owner', 0), reads=v.get('n_chain_reads', 0),
                              gene_name=cands.get(g, {}).get('gene_name', g)))
portfolio.sort(key=lambda r: (-len(r['genuine']), -r['owner_frac']))
ntx = len(set(t for r in portfolio for t in r['genuine']))
print(f"\n=== FINAL PORTFOLIO: {len(portfolio)} loci, {ntx} unique VG-specific transcripts "
      f"(VG-only ∧ at-candidate ∧ genome-ST-missed ∧ read-backed) ===")
print(f"{'gene':22}{'own_frac':>9}{'pure':>5}{'strict':>7}{'reads':>6}{'ntx':>4}  transcripts")
for r in portfolio:
    print(f"{r['gene'].replace('gene-',''):22}{r['owner_frac']:>9.3f}{r['pure']:>5}{r['strict']:>7}{r['reads']:>6}"
          f"{len(r['genuine']):>4}  {','.join(t.replace('rna-','') for t in r['genuine'])}")
json.dump(portfolio, open(f'{sweep}/PORTFOLIO.json', 'w'), indent=1)
print(f"\nwrote {sweep}/PORTFOLIO.json")
# also report regressions (where VG HURTS)
reg = [(g, e.get('vg_loss_vs_st', [])) for g, e in ev.items() if e.get('classification') == 'regress']
print(f"\n=== VG REGRESSIONS (loses real RefSeq tx ST keeps): {len(reg)} loci ===")
for g, lost in sorted(reg, key=lambda x: -len(x[1])):
    print(f"  {g.replace('gene-',''):22} lost {len(lost)}: {','.join(t.replace('rna-','') for t in lost[:5])}")
