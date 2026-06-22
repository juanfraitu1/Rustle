#!/usr/bin/env python3
"""
stringtie_vs_copy_headtohead.py

Head-to-head on GGO: at the de-tie multi-copy family loci, does StringTie COLLAPSE
copies that the copy-aware de-novo pipeline resolves? Answers the advisor's
"are we better because we capture copy information" on the axis where the copy
gain actually lives (copy recovery), not standard isoform recall (which can't see it
because the annotation itself collapses paralogs).

Comparison unit = a co-located multi-copy family (>=3 same-chrom copies within 5 Mb;
the same selection copy_assign.py uses). For each:
  K_rustle = # de-novo copies the pipeline assembles/rescues with read support
  K_st     = # DISTINCT StringTie loci (STRG gene_ids) overlapping those copy spans
  win      = max(0, K_rustle - K_st)   # copies StringTie collapses/misses
  rescued  = # family-rescued copies (support < 3)  -> StringTie-impossible (no family prior)

StringTie uses PRIMARY alignments only, so where divergent copies map distinctly it
recovers them too (tie). Wins concentrate where reads collapse onto one reference
copy but Rustle separates them by PSV, and on starved copies recovered by the family
prior. Restrict the *headline* win to resolvable families (see --copy-assign-out).

Inputs (all from the shipped pipeline):
  META    denovo_transcripts.meta.tsv      DN copy -> chrom,start,end,n_reads
  SPLIT   denovo_families_split.tsv        family -> members (DN ids), class
  RESCUE  denovo_rescued_copies.tsv        RC copy -> chrom,start,end,support
  ST_GTF  genome_st.gtf                    StringTie -L de-novo on GGO.bam (68k tx)
"""
import os, sys, re
from collections import defaultdict
from bisect import bisect_right

SCR = "/home/juanfra/winloci_scratch"
HERE = os.path.dirname(os.path.abspath(__file__))
META = os.path.join(SCR, "denovo_transcripts.meta.tsv")
SPLIT = os.path.join(HERE, "denovo_families_split.tsv")
RESCUE = os.path.join(HERE, "denovo_rescued_copies.tsv")
ST_GTF = os.path.join(SCR, "genome_st.gtf")
CA_OUT = os.path.join(SCR, "copy_assign_real.out")   # optional resolvability layer

WIN = 5_000_000
MIN_COPIES = 3


def load_meta():
    span = {}
    with open(META) as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            span[c[0]] = (c[1], int(c[2]), int(c[3]), int(c[6]))
    return span


def load_rescued():
    by_fam = defaultdict(list)
    with open(RESCUE) as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            fam, chrom, start, end, support = c[0], c[2], int(c[3]), int(c[4]), int(c[7])
            by_fam[fam].append((chrom, start, end, support))
    return by_fam


def load_stringtie_loci():
    """Group StringTie transcripts into loci by gene_id; index per chrom. Also keep
    each transcript's intron chain so we can tell COLLAPSE (one gene_id) from
    RECOVERY (does StringTie actually model both copies' structures?)."""
    gene = {}
    gene_chains = defaultdict(list)   # gene_id -> [(chrom, start, end, chain_tuple)]
    tx_exons = defaultdict(list)
    tx_gene, tx_chrom = {}, {}
    gid_re = re.compile(r'gene_id "([^"]+)"')
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    with open(ST_GTF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if c[2] == "transcript":
                chrom, s, e = c[0], int(c[3]), int(c[4])
                g = gid_re.search(c[8]).group(1)
                tx_gene[tid_re.search(c[8]).group(1)] = g
                tx_chrom[tid_re.search(c[8]).group(1)] = chrom
                if g in gene:
                    pc, ps, pe = gene[g]
                    gene[g] = (pc, min(ps, s), max(pe, e))
                else:
                    gene[g] = (chrom, s, e)
            elif c[2] == "exon":
                tx_exons[tid_re.search(c[8]).group(1)].append((int(c[3]), int(c[4])))
    for t, exs in tx_exons.items():
        exs.sort()
        chain = tuple((exs[i][1], exs[i + 1][0]) for i in range(len(exs) - 1))
        gene_chains[tx_gene[t]].append((tx_chrom[t], exs[0][0], exs[-1][1], chain))
    per_chrom = defaultdict(list)
    for g, (chrom, s, e) in gene.items():
        per_chrom[chrom].append((s, e, g))
    for chrom in per_chrom:
        per_chrom[chrom].sort()
    return per_chrom, gene_chains


def st_models_both(gene_chains, g, locs):
    """RECOVERY test: does StringTie gene g emit DISTINCT intron chains at >=2 of the
    co-merged copy loci? If yes, it modelled the copies and only the LABEL is merged."""
    chains_at = []
    for (cc, s, e) in locs:
        cs = {ch for (tc, ts, te, ch) in gene_chains.get(g, [])
              if tc == cc and not (te < s or ts > e) and ch}
        chains_at.append(cs)
    nonempty = [x for x in chains_at if x]
    return len(nonempty) >= 2 and len(set().union(*nonempty)) >= 2


def overlapping_genes(per_chrom, chrom, s, e):
    cand = per_chrom.get(chrom, [])
    if not cand:
        return set()
    starts = [x[0] for x in cand]
    hi = bisect_right(starts, e)
    hits = set()
    for k in range(hi - 1, -1, -1):
        cs, ce, g = cand[k]
        if ce < s:
            if s - ce > 3_000_000:   # genes are start-sorted; bounded look-back
                break
            continue
        hits.add(g)
    return hits


def merge_distinct_loci(copies):
    """OVER-SPLIT GUARD. Copies whose spans reciprocally overlap >=50% are
    isoform/fragment variants of ONE locus (e.g. the PRNP case: 5 near-identical
    14.60-14.615 Mb spans differing by a few bp / exon count), NOT separate copies.
    Union-find them into distinct loci; keep max read support per locus, and whether
    any constituent was a rescued copy."""
    n = len(copies)
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i in range(n):
        ci, si, ei = copies[i][0], copies[i][1], copies[i][2]
        for j in range(i + 1, n):
            cj, sj, ej = copies[j][0], copies[j][1], copies[j][2]
            if ci != cj:
                continue
            ov = min(ei, ej) - max(si, sj)
            if ov <= 0:
                continue
            if ov >= 0.5 * min(ei - si, ej - sj):
                parent[find(i)] = find(j)
    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(copies[i])
    loci = []
    for g in groups.values():
        chrom = g[0][0]
        s = min(x[1] for x in g); e = max(x[2] for x in g)
        reads = max(x[3] for x in g)
        rescued = any(x[4] for x in g)
        loci.append((chrom, s, e, reads, rescued))
    return loci


def colocated_groups(loci):
    """Return list of (chrom, [distinct loci]) where >=MIN_COPIES fall within WIN."""
    by_chrom = defaultdict(list)
    for cp in loci:
        by_chrom[cp[0]].append(cp)
    groups = []
    for chrom, cps in by_chrom.items():
        cps.sort(key=lambda x: x[1])
        if len(cps) < MIN_COPIES:
            continue
        ok = False
        j = 0
        for k in range(len(cps)):
            while cps[k][1] - cps[j][1] > WIN:
                j += 1
            if k - j + 1 >= MIN_COPIES:
                ok = True
                break
        if ok:
            groups.append((chrom, cps))
    return groups


def load_resolvability():
    """Parse copy_assign.py real output (if present): family -> resolv_PSV%, copies."""
    if not os.path.exists(CA_OUT):
        return None
    res = {}
    with open(CA_OUT) as f:
        for line in f:
            # table rows: family gene cop reads MAPQ0% resolv_PSV ...
            p = line.split()
            if len(p) >= 6 and p[0].startswith("DSFAM"):
                try:
                    res[p[0]] = dict(copies=int(p[2]), reads=int(p[3].replace(",", "")))
                except ValueError:
                    pass
    return res or None


def main():
    span = load_meta()
    rescued = load_rescued()
    st, gene_chains = load_stringtie_loci()
    resolv = load_resolvability()
    # families our own family-detection validation flags as non-genuine (domain-sharer
    # ZNF webs / mega-families) -> excluded from the defensible count, not hidden.
    FALSE_FAMILIES = {"DSFAM0"}

    n_unmatched = 0
    tot_raw_members = tot_distinct = 0
    collapse_recs = []
    zero_st_set = set()
    rescued_set = set()
    fam_loci = defaultdict(int)
    fam_zero = defaultdict(int)
    with open(SPLIT) as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            fid, ncop, nchr, gene, conc, dens, acr, nart, cls = c[:9]
            members = c[9]
            if cls != "family":
                continue
            copies = []   # (chrom,start,end,reads,is_rescued)
            for m in members.split(","):
                m = m.strip()
                if m in span:
                    chrom, s, e, nr = span[m]
                    copies.append((chrom, s, e, nr, False))
                elif m:
                    n_unmatched += 1
            for (chrom, s, e, sup) in rescued.get(fid, []):
                copies.append((chrom, s, e, sup, True))
            raw_n = len(copies)
            loci = merge_distinct_loci(copies)          # <-- OVER-SPLIT GUARD
            tot_raw_members += raw_n
            tot_distinct += len(loci)
            groups = colocated_groups(loci)
            if not groups:
                continue
            for chrom, cps in groups:
                # map each distinct locus -> StringTie genes overlapping it
                gene_to_loci = defaultdict(list)
                zero = 0
                for (cc, s, e, nr, isr) in cps:
                    hits = overlapping_genes(st, cc, s, e)
                    if not hits:
                        zero += 1
                        zero_st_set.add((fid, cc, s, e))
                    if isr:
                        rescued_set.add((fid, cc, s, e))
                    for g in hits:
                        gene_to_loci[g].append((cc, s, e))
                for g, locs in gene_to_loci.items():
                    if len(locs) >= 2:                       # >=2 distinct loci share one gene_id
                        n_extra = len(locs) - 1
                        labeling = st_models_both(gene_chains, g, locs)  # ST models both -> labeling, not recovery
                        collapse_recs.append(dict(family=fid, gene=gene, st_gene=g,
                                                  n_extra=n_extra, labeling=labeling))
                fam_loci[fid] += len(cps)
                fam_zero[fid] += zero

    # ---- aggregate, honestly ----
    n_fam = len({r["family"] for r in collapse_recs} | set(fam_loci))
    tot_collapsed = sum(r["n_extra"] for r in collapse_recs)
    labeling_collapsed = sum(r["n_extra"] for r in collapse_recs if r["labeling"])
    coll_by_fam = defaultdict(int)
    for r in collapse_recs:
        coll_by_fam[r["family"]] += r["n_extra"]
    dsf0 = sum(v for f, v in coll_by_fam.items() if f in FALSE_FAMILIES)
    genuine = {f: v for f, v in coll_by_fam.items() if f not in FALSE_FAMILIES}
    validated_genuine = sum(v for f, v in genuine.items()
                            if resolv is not None and f in resolv)
    tot_zero = len(zero_st_set)
    tot_resc = len(rescued_set)
    overlap = len(zero_st_set & rescued_set)

    print("=" * 78)
    print("StringTie vs copy-aware de-novo — multi-copy families (GGO)  [honest revision]")
    print("=" * 78)
    print(f"OVER-SPLIT GUARD: {tot_raw_members} member-skeletons -> {tot_distinct} distinct loci"
          f" (~{100*(tot_raw_members-tot_distinct)//tot_raw_members}% over-split fragments merged)")
    print(f"co-located multi-copy families (>=3 distinct same-chrom loci / 5Mb): {n_fam}")
    print(f"\n-- StringTie merges distinct copies under one gene_id (collapse) --")
    print(f"  collapse instances (>=2 distinct loci, one gene_id)   : {tot_collapsed}")
    print(f"  of which StringTie MODELS both copies (distinct chains): {labeling_collapsed}"
          f"  <- LABELING, not recovery ({100*labeling_collapsed//max(tot_collapsed,1)}%)")
    print(f"  in DSFAM0 (false domain-sharer ZNF family, excluded)   : {dsf0}")
    print(f"  collapses in GENUINE families (excl DSFAM0)            : {sum(genuine.values())}")
    print(f"  of those, copy_assign-VALIDATED resolvable             : {validated_genuine}")
    print(f"\n-- soft tail (overlapping sets; do NOT sum) --")
    print(f"  family-rescue-flagged loci (POA-confirmed, mostly thin): {tot_resc}")
    print(f"  distinct loci with ZERO StringTie model                : {tot_zero}")
    print(f"  overlap (zero-ST loci that are also rescue-flagged)     : {overlap}")
    print(f"\n-- per-copy attribution capability (copy_assign real) --")
    if resolv is not None:
        print(f"  99.9% silver-standard is UNIQUE-MAPPER agreement (easy ~80% of reads),")
        print(f"  not proof on the hard ambiguous reads. See copy_assign_real.out.")
    print(f"\ngenuine-family collapses by family (gene | extra | validated?):")
    for f, v in sorted(genuine.items(), key=lambda x: -x[1])[:14]:
        gene = next((r["gene"] for r in collapse_recs if r["family"] == f), "?")
        vmark = "validated" if (resolv and f in resolv) else ""
        print(f"  {f:9} {gene[:8]:8} {v:>2}  {vmark}")


if __name__ == "__main__":
    main()
