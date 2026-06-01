#!/usr/bin/env python3
"""Genome-wide scan for the DAZ pattern: annotated multi-exon genes whose
expected splice chain is covered by reads that are PREDOMINANTLY SECONDARY
alignments (their primary alignment was claimed by a sibling paralog), so a
primary-only assembler (StringTie) misses the gene.

Pipeline:
  1. Candidate generation: multi-exon genes with >= MIN_SEC secondary
     alignments starting inside the gene locus (multi-copy / paralog hotspots).
  2. Pass A (per candidate, samtools region query): chain-supporting reads
     (alignment on the gene strand sharing >= K of the gene's annotated
     junctions); record NM-at-gene and whether THIS alignment is primary.
  3. Global primary lookup: one BAM pass over primaries, filtered to the read
     names whose chain-supporting alignment at a candidate is SECONDARY.
  4. Pass B (classify): for each candidate, primary_at_gene fraction; for
     primary-elsewhere reads, resolve the sibling gene at the primary locus,
     test whether the read ALSO matches the sibling's chain, and compare NM.

Classes:
  primary_sufficient            -- >25% of chain evidence is primary-at-gene
  secondary_dependent_identif.  -- mostly secondary, but sequence favors THIS
                                   gene (NM_gene <= NM_primary) or no sibling
                                   chain -> clean recovery via secondaries
  secondary_dependent_ambiguous -- DAZ-like: mostly secondary, reads ALSO match
                                   a sibling's chain AND sequence favors the
                                   sibling -> genuinely non-identifiable

Usage: run_scan.py GENE_INTRONS.tsv SECONDARY_POS.tsv BAM OUT_PREFIX
Env: SCAN_MIN_SEC(10) SCAN_MINEXPR(5) SCAN_K(0) SCAN_SECFRAC(0.75) SCAN_NM_MARGIN(0.0)
"""
import sys, os, re, subprocess, json, bisect, collections, statistics as st

MIN_SEC   = int(os.environ.get('SCAN_MIN_SEC', '10'))
MIN_EXPR  = int(os.environ.get('SCAN_MINEXPR', '5'))
K_FIXED   = int(os.environ.get('SCAN_K', '0'))     # 0 => adaptive
SECFRAC   = float(os.environ.get('SCAN_SECFRAC', '0.75'))
NM_MARGIN = float(os.environ.get('SCAN_NM_MARGIN', '0.0'))

def junctions(pos, cig):
    cur = pos; js = []
    for ln, op in re.findall(r'(\d+)([MIDNSH=X])', cig):
        ln = int(ln)
        if op == 'N':
            js.append((cur - 1, cur + ln)); cur += ln
        elif op in 'MD=X':
            cur += ln
    return js

def aligned_len(cig):
    return sum(int(l) for l, o in re.findall(r'(\d+)([MIDNSH=X])', cig)
               if o in 'M=XD')

def load_genes(path):
    genes = {}
    by_chrom = collections.defaultdict(list)   # chrom -> [(start,end,gid)]
    for line in open(path):
        if line.startswith('gene_id'):
            continue
        gid, chrom, start, end, strand, ntx, chain = line.rstrip('\n').split('\t')
        ints = set()
        if chain:
            for pr in chain.split(','):
                d, a = pr.split(':'); ints.add((int(d), int(a)))
        genes[gid] = dict(chrom=chrom, start=int(start), end=int(end),
                          strand=strand, introns=ints)
        by_chrom[chrom].append((int(start), int(end), gid))
    for c in by_chrom:
        by_chrom[c].sort()
    return genes, by_chrom

def gene_at(by_chrom, chrom, pos):
    """return gid whose locus contains pos (first hit), or None."""
    arr = by_chrom.get(chrom)
    if not arr:
        return None
    i = bisect.bisect_right(arr, (pos, float('inf'), '')) - 1
    # scan back a few in case of nested/overlapping genes
    for j in range(i, max(-1, i - 6), -1):
        if j < 0:
            break
        s, e, gid = arr[j]
        if s <= pos <= e:
            return gid
    return None

def main():
    gi, secpos, bam, outpre = sys.argv[1:5]
    genes, by_chrom = load_genes(gi)
    # ---- candidate generation ----
    sec_by_chrom = collections.defaultdict(list)
    for line in open(secpos):
        c, p = line.rstrip('\n').split('\t')
        sec_by_chrom[c].append(int(p))
    for c in sec_by_chrom:
        sec_by_chrom[c].sort()
    cands = []
    for gid, g in genes.items():
        if not g['introns']:
            continue
        arr = sec_by_chrom.get(g['chrom'])
        if not arr:
            continue
        lo = bisect.bisect_left(arr, g['start'])
        hi = bisect.bisect_right(arr, g['end'])
        nsec = hi - lo
        if nsec >= MIN_SEC:
            cands.append((nsec, gid))
    cands.sort(reverse=True)
    print(f"[scan] {len(cands)} candidate multi-exon genes "
          f"(>= {MIN_SEC} secondaries in locus)", file=sys.stderr)

    # ---- Pass A: per-candidate region query ----
    def kmin(g):
        n = len(g['introns'])
        return K_FIXED if K_FIXED > 0 else max(2, min(8, (n + 2) // 3))
    candinfo = {}                 # gid -> dict
    need_primary = set()          # read names whose G-alignment is secondary
    for nsec, gid in cands:
        g = genes[gid]
        region = f"{g['chrom']}:{g['start']}-{g['end']}"
        try:
            out = subprocess.run(['samtools', 'view', bam, region],
                                 capture_output=True, text=True,
                                 timeout=120).stdout
        except Exception:
            continue
        K = kmin(g)
        reads = {}                # qname -> dict(nm, is_primary_here, nmatch)
        for line in out.splitlines():
            f = line.split('\t')
            if len(f) < 11:
                continue
            flag = int(f[1])
            if flag & 0x800:      # supplementary
                continue
            strand = '-' if flag & 16 else '+'
            if strand != g['strand']:
                continue
            js = junctions(int(f[3]), f[5])
            nmatch = sum(1 for j in js if j in g['introns'])
            if nmatch < K:
                continue
            nm = next((int(t[5:]) for t in f[11:] if t.startswith('NM:i:')), None)
            al = aligned_len(f[5]) or 1
            is_prim = not (flag & 0x100)
            r = reads.get(f[0])
            cand = dict(nm=nm, alen=al, is_primary_here=is_prim, nmatch=nmatch)
            # keep the best (most junctions) alignment of this read at G
            if r is None or nmatch > r['nmatch']:
                reads[f[0]] = cand
        if len(reads) < MIN_EXPR:
            continue
        candinfo[gid] = dict(nsec=nsec, K=K, reads=reads)
        for qn, r in reads.items():
            if not r['is_primary_here']:
                need_primary.add(qn)

    print(f"[scan] {len(candinfo)} expressed candidates; "
          f"{len(need_primary)} read names need primary lookup", file=sys.stderr)

    # ---- global primary lookup (one pass over primaries) ----
    primary = {}                  # qname -> (chrom,pos,strand,cigar,nm,alen)
    proc = subprocess.Popen(['samtools', 'view', '-F', '0x900', bam],
                            stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        i = line.find('\t')
        qn = line[:i]
        if qn not in need_primary:
            continue
        f = line.split('\t')
        flag = int(f[1]); strand = '-' if flag & 16 else '+'
        nm = next((int(t[5:]) for t in f[11:] if t.startswith('NM:i:')), None)
        primary[qn] = (f[2], int(f[3]), strand, f[5], nm, aligned_len(f[5]) or 1)
    proc.wait()
    print(f"[scan] resolved primaries for {len(primary)} reads", file=sys.stderr)

    # ---- Pass B: classify ----
    rows = []
    summary = collections.Counter()
    for gid, ci in candinfo.items():
        g = genes[gid]; reads = ci['reads']
        n_chain = len(reads)
        prim_here = sum(1 for r in reads.values() if r['is_primary_here'])
        prim_frac = prim_here / n_chain
        # analyze primary-elsewhere reads
        sib_counts = collections.Counter()
        n_sib_chain = 0; n_seq_favors_sib = 0; n_elsewhere = 0
        nm_g = []; nm_p = []
        for qn, r in reads.items():
            if r['is_primary_here']:
                continue
            n_elsewhere += 1
            p = primary.get(qn)
            if p is None:
                continue
            pc, pp, pstr, pcig, pnm, palen = p
            sib = gene_at(by_chrom, pc, pp)
            sib_counts[sib if sib else f"{pc}:intergenic"] += 1
            # does the read's PRIMARY match the sibling gene's chain?
            if sib and genes[sib]['introns']:
                pjs = junctions(pp, pcig)
                pm = sum(1 for j in pjs if j in genes[sib]['introns'])
                if pm >= max(2, ci['K']):
                    n_sib_chain += 1
            # sequence preference: NM-rate at primary vs at G
            if r['nm'] is not None and pnm is not None:
                rg = r['nm'] / r['alen']; rp = pnm / palen
                nm_g.append(rg); nm_p.append(rp)
                if rp <= rg + NM_MARGIN:
                    n_seq_favors_sib += 1
        dom_sib, dom_n = (sib_counts.most_common(1)[0]
                          if sib_counts else (None, 0))
        # classification
        if prim_frac > (1 - SECFRAC):
            cls = 'primary_sufficient'
        else:
            denom = max(1, n_elsewhere)
            sib_chain_frac = n_sib_chain / denom
            seq_sib_frac = (n_seq_favors_sib / len([x for x in nm_g]) ) if nm_g else 0.0
            if sib_chain_frac >= 0.5 and seq_sib_frac >= 0.5:
                cls = 'secondary_dependent_ambiguous'
            else:
                cls = 'secondary_dependent_identifiable'
        summary[cls] += 1
        rows.append(dict(
            gene=gid, chrom=g['chrom'], start=g['start'], end=g['end'],
            strand=g['strand'], n_introns=len(g['introns']), K=ci['K'],
            n_chain_reads=n_chain, primary_at_gene=prim_here,
            primary_frac=round(prim_frac, 3), n_secondary=n_elsewhere,
            n_sib_chain_match=n_sib_chain,
            n_seq_favors_sibling=n_seq_favors_sib,
            nm_rate_gene=round(st.median(nm_g), 4) if nm_g else None,
            nm_rate_primary=round(st.median(nm_p), 4) if nm_p else None,
            dom_sibling=dom_sib, dom_sibling_n=dom_n, klass=cls))
    rows.sort(key=lambda r: (r['klass'] != 'secondary_dependent_ambiguous',
                             -r['n_chain_reads']))
    with open(outpre + '.tsv', 'w') as o:
        cols = ['gene', 'chrom', 'start', 'end', 'strand', 'n_introns', 'K',
                'n_chain_reads', 'primary_at_gene', 'primary_frac',
                'n_secondary', 'n_sib_chain_match', 'n_seq_favors_sibling',
                'nm_rate_gene', 'nm_rate_primary', 'dom_sibling',
                'dom_sibling_n', 'klass']
        o.write('\t'.join(cols) + '\n')
        for r in rows:
            o.write('\t'.join(str(r[c]) for c in cols) + '\n')
    with open(outpre + '.summary.json', 'w') as o:
        json.dump(dict(n_candidates=len(cands), n_expressed=len(candinfo),
                       classes=dict(summary),
                       params=dict(MIN_SEC=MIN_SEC, MIN_EXPR=MIN_EXPR,
                                   SECFRAC=SECFRAC, NM_MARGIN=NM_MARGIN)),
                  o, indent=2)
    print(f"[scan] classes: {dict(summary)}", file=sys.stderr)
    print(f"[scan] wrote {outpre}.tsv + {outpre}.summary.json", file=sys.stderr)

if __name__ == '__main__':
    main()
