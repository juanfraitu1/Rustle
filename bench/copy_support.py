#!/usr/bin/env python3
"""Spliced support per annotated copy (docs/PREREG_spliced_copy_support_2026-10-04.md): is a copy FOUND because reads are spliced
transcripts of it, or only because something overlaps it?

    copy_support.py --copies copies.tsv --truth truth.gtf --bam reads.bam --family NPIP --out PREFIX
                    [--loci NAME=loci.gff3 ...] [--nodes pagedata.json]

Per copy (rows of --copies with the given family): the primary reads (-F 2308) overlapping its annotated exon union on its strand;
junctions = N ops >= 50 bp (exact donor/acceptor); supported junction = carried by >= 3 reads at the copy; k = min(2, annotated introns
>= 50 bp); structural-support read = >= k supported junctions (k = 0: aligned blocks cover >= 50% of the exon union); spliced-expressed
= >= 2 such reads. For each --loci set: OLD = a same-strand locus whose rep exons overlap the exon union; STRICT = spliced-expressed and a
same-strand locus whose rep junctions include >= k supported junctions (k = 0: rep exons cover >= 50% of the union). --nodes restricts
the loci of each arm to the page's NPIP-cluster nodes (pagedata.json rows: cid -> node[arm]) and reports found-within-NPIP-clusters too.
Amendment A (2026-10-04): the FOUND verdict is annotation-anchored — `ann_support_reads` = reads with >= k of the copy's ANNOTATED introns,
`ann_expressed` = >= 2 of them, `<arm>_ann_found` = ann_expressed and the representative carries >= k annotated introns (`<arm>_locus_ann_found`
at the locus level); the read-defined columns (`support_reads`, `<arm>_strict_found`, ...) stay as the annotation-free reading, reported beside.
Writes PREFIX.copies.tsv and PREFIX.json.
"""
import argparse
import collections
import statistics
import csv
import json
import re

import pysam

MIN_INTRON = 50
MIN_JUNCTION_READS = 3
FLOOR = 2
COVER = 0.5


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def ilen(iv):
    return sum(e - s for s, e in iv)


def inter(a, b):
    i = j = 0
    tot = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            tot += e - s
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def gtf_transcripts(path, gene_ids):
    """gene_id -> {transcript_id: sorted exon list [[s0, e]]} for the given gene ids (GTF, 1-based closed)."""
    tx = collections.defaultdict(lambda: collections.defaultdict(list))
    for ln in open(path):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        g = re.search(r'gene_id "([^"]+)"', f[8])
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        if not g or g.group(1) not in gene_ids:
            continue
        tx[g.group(1)][t.group(1) if t else g.group(1)].append([int(f[3]) - 1, int(f[4])])
    return {g: {t: sorted(ex) for t, ex in d.items()} for g, d in tx.items()}


def introns_of(exons):
    return [(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1) if exons[i + 1][0] - exons[i][1] >= MIN_INTRON]


def read_blocks_junctions(rd):
    """aligned reference blocks [[s0, e]] and junctions [(donor_end, acceptor_start)] from the CIGAR (N >= MIN_INTRON)."""
    blocks, juncs = [], []
    pos = rd.reference_start
    cur_s = pos
    for op, ln in rd.cigartuples:
        if op in (0, 7, 8, 2):          # M, =, X, D consume the reference
            pos += ln
        elif op == 3:                   # N
            if ln >= MIN_INTRON:
                blocks.append([cur_s, pos])
                juncs.append((pos, pos + ln))
                cur_s = pos + ln
            pos += ln
    blocks.append([cur_s, pos])
    return merge([b for b in blocks if b[1] > b[0]]), juncs


def chain_class(chain, models):
    """FSM / ISM / NIC / NNC of a junction chain against model intron chains (Amendment C's class column; `noann` without models)."""
    if not models:
        return "noann"
    chain = tuple(chain)
    if any(chain == tuple(m) for m in models):
        return "FSM"
    for m in models:
        n, L = len(m), len(chain)
        if L < n and any(tuple(m[i:i + L]) == chain for i in range(n - L + 1)):
            return "ISM"
    sites = {s for m in models for j in m for s in j}
    return "NIC" if all(j[0] in sites and j[1] in sites for j in chain) else "NNC"


def chain_match(js, chains, k):
    """Amendment B: `js` (a read's or representative's junctions inside the copy span, in order) is a contiguous sub-chain of >= k junctions of
    one of `chains` (each a transcript's intron chain, in order). k = 0 never matches here (the caller uses the coverage rule)."""
    if k == 0 or len(js) < k:
        return False
    for ch in chains:
        n = len(ch)
        if len(js) > n:
            continue
        for i in range(n - len(js) + 1):
            if ch[i:i + len(js)] == js:
                return True
    return False


def load_loci(path):
    """locus name -> dict(chrom, strand, exons (merged rep exons), juncs set) from a loci GFF3 (gene + exon rows)."""
    loci = {}
    for ln in open(path):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "gene":
            loci[at["Name"]] = dict(chrom=f[0], strand=f[6], exons=[])
        elif f[2] == "exon":
            loci[at["gene"]]["exons"].append([int(f[3]) - 1, int(f[4])])
    for v in loci.values():
        ex = sorted(v["exons"])
        v["exons"] = merge(ex)
        v["juncs"] = set(introns_of(ex))
    return loci


def locus_transcript_exons(gtf):
    """locus (gene_id) -> list of its transcripts' exon lists."""
    ex = collections.defaultdict(lambda: collections.defaultdict(list))
    for ln in open(gtf):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        ex[g][t.group(1) if t else g].append([int(f[3]) - 1, int(f[4])])
    return {g: [sorted(e) for e in d.values()] for g, d in ex.items()}


def locus_transcript_chains(gtf):
    """locus (gene_id) -> list of its transcripts' junction chains (each in order)."""
    ex = collections.defaultdict(lambda: collections.defaultdict(list))
    for ln in open(gtf):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        ex[g][t.group(1) if t else g].append([int(f[3]) - 1, int(f[4])])
    return {g: [introns_of(sorted(e)) for e in d.values()] for g, d in ex.items()}


def locus_transcript_junctions(gtf):
    """locus (gene_id) -> set of junctions over ALL its transcripts (the locus-level reading, reported beside the representative's)."""
    ex = collections.defaultdict(lambda: collections.defaultdict(list))
    for ln in open(gtf):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        ex[g][t.group(1) if t else g].append([int(f[3]) - 1, int(f[4])])
    return {g: set().union(*(set(introns_of(sorted(e))) for e in d.values())) for g, d in ex.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--loci", action="append", default=[], help="NAME=loci.gff3[,transcripts.gtf] (repeatable); with the GTF the locus-level reading is reported beside")
    ap.add_argument("--nodes", help="pagedata.json of the read-pool page (its node[arm] per cid)")
    ap.add_argument("--chain-floor", type=int, default=3, help="Amendment C: unique reads an identical chain needs to be an expressed chain")
    ap.add_argument("--tss-tol", type=int, default=150, help="Amendment D: |read 5' end - transcript TSS| tolerance (bp); 50 and 300 reported beside")
    ap.add_argument("--truth2", help="Amendment B: a second annotation's truth GTF (e.g. RefSeq) whose models of the same gene also count")
    ap.add_argument("--copies2", help="the copies table of --truth2 (its cid keys the GTF; genes matched to --copies by name / refseq_name)")
    a = ap.parse_args(argv)

    copies = [r for r in csv.DictReader(open(a.copies), delimiter="\t") if r["family"] == a.family]
    # the truth GTF keys its genes by the copy id (cid) in the recovery benchmarks' truth files, by the annotation id elsewhere
    genes = {r["isoform_gene"] for r in copies} | {r["cid"] for r in copies}
    tx = gtf_transcripts(a.truth, genes)
    tx2 = {}
    if a.truth2 and a.copies2:
        c2 = list(csv.DictReader(open(a.copies2), delimiter="\t"))
        by_name = {r["name"]: r["cid"] for r in c2}
        t2 = gtf_transcripts(a.truth2, {r["cid"] for r in c2})
        for c in copies:
            nm = c.get("refseq_name") or c["name"]
            if nm in by_name and by_name[nm] in t2:
                tx2[c["cid"]] = t2[by_name[nm]]
    arms = {}
    arm_tx, arm_tx_chains, arm_tx_exons = {}, {}, {}
    for spec in a.loci:
        name, path = spec.split("=", 1)
        path, _, gtf = path.partition(",")
        arms[name] = load_loci(path)
        if gtf:
            arm_tx[name] = locus_transcript_junctions(gtf)
            arm_tx_chains[name] = locus_transcript_chains(gtf)
            arm_tx_exons[name] = locus_transcript_exons(gtf)
    nodes = {}
    if a.nodes:
        pd = json.load(open(a.nodes))
        for row in pd.get("rows", pd if isinstance(pd, list) else []):
            nodes[row["cid"]] = row.get("node", {})

    bam = pysam.AlignmentFile(a.bam)
    rows, summary = [], {"copies": len(copies), "spliced_expressed": 0, "ann_expressed": 0, "arms": {}}
    for arm in arms:
        summary["arms"][arm] = {"old_overlap": 0, "strict_found": 0, "old_overlap_in_npip_nodes": 0, "strict_found_in_npip_nodes": 0,
                                "locus_level_found": 0, "locus_level_found_in_npip_nodes": 0}
    for c in copies:
        chrom, strand = c["chrom"], c["strand"]
        texons = tx.get(c["cid"]) or tx.get(c["isoform_gene"], {})
        texons2 = tx2.get(c["cid"], {})
        label = c.get("refseq_name") or c.get("cat_name") or c["name"]
        # Amendment B: intron chains of every model (annotation 1, annotation 2, and both)
        chains1 = [introns_of(ex) for ex in texons.values() if len(introns_of(ex)) >= 1]
        chains2 = [introns_of(ex) for ex in texons2.values() if len(introns_of(ex)) >= 1]
        chains = chains1 + chains2
        n_intron_b = max((len(ch) for ch in chains), default=0)
        kb = min(2, n_intron_b)
        # Amendment D: models as (TSS, introns oriented 5'->3', exon union) per annotation
        def oriented(exs):
            ex = sorted(exs)
            intr = introns_of(ex)
            if strand == "-":
                return ex[-1][1], list(reversed(intr)), ex
            return ex[0][0], intr, ex
        models_d1 = [oriented(ex) for ex in texons.values() if ex]
        models_d2 = [oriented(ex) for ex in texons2.values() if ex]
        models_d = models_d1 + models_d2

        def tss_support(p5, js, mdl, tol):
            """read 5' end p5 and junctions js (sorted by coordinate) vs one model: within tol of the TSS and first m junctions = first m introns."""
            tss, intr, ex = mdl
            if abs(p5 - tss) > tol:
                return False
            m = min(3, len(intr))
            if m == 0:
                return True  # intronless model: the caller checks coverage
            rj = list(reversed(js)) if strand == "-" else list(js)
            return len(rj) >= m and rj[:m] == intr[:m]

        def td_count(models, tol, uniq=False):
            n = 0
            for (bl, js), q, p5 in zip(reads, read_mapq, read_5p):
                if uniq and q == 0:
                    continue
                for mdl in models:
                    if tss_support(p5, js, mdl, tol) and (len(mdl[1]) > 0 or inter(bl, mdl[2]) >= COVER * ilen(mdl[2])):
                        n += 1
                        break
            return n

        union = merge([e for ex in texons.values() for e in ex]) or [[int(c["terr_lo0"]), int(c["terr_hi"])]]
        ann_introns = {t: set(introns_of(ex)) for t, ex in texons.items()}
        n_intron = max((len(v) for v in ann_introns.values()), default=0)
        k = min(2, n_intron)
        ann_all = set().union(*ann_introns.values()) if ann_introns else set()
        lo, hi = union[0][0], union[-1][1]
        reads, read_mapq, read_5p = [], [], []
        for rd in bam.fetch(chrom, lo, hi):
            if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                continue
            if ("-" if rd.is_reverse else "+") != strand:
                continue
            blocks, juncs = read_blocks_junctions(rd)
            if inter(blocks, union) == 0:
                continue
            juncs = [j for j in juncs if j[0] >= lo and j[1] <= hi]
            reads.append((blocks, juncs))
            read_mapq.append(rd.mapping_quality)
            read_5p.append(rd.reference_end if strand == "-" else rd.reference_start)
        jcount = collections.Counter(j for _, js in reads for j in js)
        supported = {j for j, n in jcount.items() if n >= MIN_JUNCTION_READS}
        n_unspliced = sum(1 for _, js in reads if not js)
        n_one = sum(1 for _, js in reads if len(js) == 1)
        if k == 0:
            n_support = sum(1 for bl, _ in reads if inter(bl, union) >= COVER * ilen(union))
        else:
            n_support = sum(1 for _, js in reads if sum(1 for j in js if j in supported) >= k)
        n_ann2 = sum(1 for _, js in reads if sum(1 for j in js if j in ann_all) >= min(2, n_intron)) if n_intron else 0
        # Amendment A: support = the copy's OWN annotated introns (k of them; k = 0: coverage of the exon union)
        n_ann_support = n_ann2 if k else n_support
        ann_expressed = n_ann_support >= FLOOR
        # Amendment B: chain support (contiguous sub-chain of a model; both annotations, and each alone)
        if kb == 0:
            n_chain = n_chain1 = n_chain2 = n_support
        else:
            n_chain = sum(1 for _, js in reads if chain_match(js, chains, kb))
            n_chain1 = sum(1 for _, js in reads if chain_match(js, chains1, kb))
            n_chain2 = sum(1 for _, js in reads if chain_match(js, chains2, kb))
        chain_expressed = n_chain >= FLOOR
        # Amendment C: expressed chains = identical >= 2-junction chains of >= chain_floor UNIQUE reads; support = equal or contiguous sub-chain
        cc_u, cc_all = collections.Counter(), collections.Counter()
        for (_, js), q in zip(reads, read_mapq):
            if len(js) >= 2:
                cc_all[tuple(js)] += 1
                if q > 0:
                    cc_u[tuple(js)] += 1
        expressed_chains = [list(ch) for ch, n in cc_u.items() if n >= a.chain_floor]
        expressed_all = [list(ch) for ch, n in cc_all.items() if n >= a.chain_floor]      # tied reads included (beside)
        xc_support = sum(1 for (_, js), q in zip(reads, read_mapq) if q > 0 and len(js) >= 2 and chain_match(js, expressed_chains, 2))
        xc_support_all = sum(1 for (_, js) in reads if len(js) >= 2 and chain_match(js, expressed_all, 2))
        xc_expressed = len(expressed_chains) >= 1 if kb else n_support >= FLOOR
        n_exp2 = sum(1 for ch, n in cc_u.items() if n >= 2)
        n_exp5 = sum(1 for ch, n in cc_u.items() if n >= 5)
        dom = max(cc_u.items(), key=lambda kv: (kv[1], len(kv[0])), default=None)
        dom_cls = chain_class(dom[0], chains) if dom else "-"
        dom_reads = dom[1] if dom else 0
        exp_classes = collections.Counter(chain_class(ch, chains) for ch in expressed_chains)
        td_150 = td_count(models_d, a.tss_tol)
        td_50, td_300 = td_count(models_d, 50), td_count(models_d, 300)
        td_cat, td_rs = td_count(models_d1, a.tss_tol), td_count(models_d2, a.tss_tol)
        td_uniq = td_count(models_d, a.tss_tol, uniq=True)
        td_expressed = td_150 >= FLOOR
        # D' (reported beside): the same test against the reads' own expressed chains; TSS_e = modal 5' end (20-bp bins) of the chain's unique carriers
        te_models = []
        for ch in expressed_chains:
            p5s = [p5 for (_, js), q, p5 in zip(reads, read_mapq, read_5p) if q > 0 and tuple(js) == tuple(ch)]
            if not p5s:
                continue
            binned = collections.Counter(p // 20 for p in p5s)
            b = max(binned.items(), key=lambda kv: (kv[1], -abs(kv[0])))[0]
            tss_e = int(statistics.median([p for p in p5s if p // 20 == b]))
            intr = list(reversed(ch)) if strand == "-" else list(ch)
            te_models.append((tss_e, intr, union))
        te_150 = td_count(te_models, a.tss_tol) if te_models else 0
        te_uniq = td_count(te_models, a.tss_tol, uniq=True) if te_models else 0
        te_expressed = te_150 >= FLOOR
        n_exact = sum(1 for _, js in reads if js and any(set(js) == s for s in ann_introns.values() if s))
        expressed = n_support >= FLOOR
        summary["spliced_expressed"] += expressed
        summary["ann_expressed"] += ann_expressed
        summary["chain_expressed"] = summary.get("chain_expressed", 0) + chain_expressed
        summary["xc_expressed"] = summary.get("xc_expressed", 0) + xc_expressed
        summary["td_expressed"] = summary.get("td_expressed", 0) + td_expressed
        summary["te_expressed"] = summary.get("te_expressed", 0) + te_expressed
        row = dict(cid=c["cid"], name=label, chrom=chrom, strand=strand, span=f"{lo}-{hi}", n_tx=len(texons), ann_introns=n_intron, k=k,
                   reads=len(reads), unspliced=n_unspliced, one_junction=n_one, supported_junctions=len(supported), support_reads=n_support,
                   ann2_reads=n_ann2, exact_chain_reads=n_exact, spliced_expressed=int(expressed),
                   ann_support_reads=n_ann_support, ann_expressed=int(ann_expressed),
                   chain_support_reads=n_chain, chain_support_ann1=n_chain1, chain_support_ann2=n_chain2, chain_expressed=int(chain_expressed),
                   models_ann1=len(chains1), models_ann2=len(chains2),
                   xc_expressed_chains=len(expressed_chains), xc_expressed_chains_ge2=n_exp2, xc_expressed_chains_ge5=n_exp5,
                   xc_support_reads=xc_support, xc_support_reads_incl_tied=xc_support_all, xc_expressed=int(xc_expressed),
                   xc_dominant_reads=dom_reads, xc_dominant_junctions=len(dom[0]) if dom else 0, xc_dominant_class=dom_cls,
                   xc_expressed_FSM_ISM=exp_classes["FSM"] + exp_classes["ISM"], xc_expressed_NIC=exp_classes["NIC"], xc_expressed_NNC=exp_classes["NNC"],
                   td_support_reads=td_150, td_support_50=td_50, td_support_300=td_300, td_support_cat=td_cat, td_support_refseq=td_rs,
                   td_support_unique=td_uniq, td_expressed=int(td_expressed),
                   te_support_reads=te_150, te_support_unique=te_uniq, te_expressed=int(te_expressed), te_chains=len(te_models),
                   te_tss=";".join(str(m[0]) for m in te_models[:3]))
        for arm, loci in arms.items():
            same = [L for L in loci.values() if L["chrom"] == chrom and L["strand"] == strand and inter(L["exons"], union) > 0]
            if k == 0:
                strict_loci = [L for L in same if inter(L["exons"], union) >= COVER * ilen(union)]
            else:
                strict_loci = [L for L in same if len(L["juncs"] & supported) >= k]
            old = len(same) > 0
            strict = expressed and len(strict_loci) > 0
            row[f"{arm}_old_overlap_loci"] = len(same)
            row[f"{arm}_strict_loci"] = len(strict_loci)
            row[f"{arm}_strict_found"] = int(strict)
            row[f"{arm}_rep_supported_junctions_max"] = max((len(L["juncs"] & supported) for L in same), default=0)
            # Amendment A: FOUND = annotation-anchored — the representative carries >= k of the copy's ANNOTATED introns
            rep_ann = max((len(L["juncs"] & ann_all) for L in same), default=0)
            ann_found = ann_expressed and ((len(strict_loci) > 0) if k == 0 else rep_ann >= k)
            row[f"{arm}_rep_ann_junctions_max"] = rep_ann
            row[f"{arm}_ann_found"] = int(ann_found)
            summary["arms"][arm]["ann_found"] = summary["arms"][arm].get("ann_found", 0) + ann_found
            # Amendment B: the representative's junction chain inside the span must be a contiguous sub-chain of a model
            def rep_chain_ok(L):
                js = sorted(j for j in L["juncs"] if j[0] >= lo and j[1] <= hi)
                return chain_match(js, chains, kb)
            chain_loci = [L for L in same if (inter(L["exons"], union) >= COVER * ilen(union)) if kb == 0] if kb == 0 else [L for L in same if rep_chain_ok(L)]
            chain_found = chain_expressed and len(chain_loci) > 0
            row[f"{arm}_chain_found"] = int(chain_found)
            summary["arms"][arm]["chain_found"] = summary["arms"][arm].get("chain_found", 0) + chain_found
            # Amendment C: the representative's in-span chain equals / is a sub-chain of an EXPRESSED chain of the copy
            def rep_xc_ok(L):
                js = sorted(j for j in L["juncs"] if j[0] >= lo and j[1] <= hi)
                return chain_match(js, expressed_chains, 2)
            xc_loci = [L for L in same if inter(L["exons"], union) >= COVER * ilen(union)] if not kb else [L for L in same if rep_xc_ok(L)]
            xc_found = xc_expressed and len(xc_loci) > 0
            row[f"{arm}_xc_found"] = int(xc_found)
            summary["arms"][arm]["xc_found"] = summary["arms"][arm].get("xc_found", 0) + xc_found
            # Amendment D: the representative starts at a model's TSS (± tol) and carries its first m introns
            def rep_td_ok(L):
                ex = L["exons"]
                p5 = ex[-1][1] if strand == "-" else ex[0][0]
                js = sorted(L["juncs"])
                return any(tss_support(p5, js, mdl, a.tss_tol) and (len(mdl[1]) > 0 or inter(ex, mdl[2]) >= COVER * ilen(mdl[2])) for mdl in models_d)
            td_loci = [L for L in same if rep_td_ok(L)]
            td_found = td_expressed and len(td_loci) > 0
            row[f"{arm}_td_found"] = int(td_found)
            summary["arms"][arm]["td_found"] = summary["arms"][arm].get("td_found", 0) + td_found
            def rep_te_ok(L):
                ex = L["exons"]
                p5 = ex[-1][1] if strand == "-" else ex[0][0]
                return any(tss_support(p5, sorted(L["juncs"]), mdl, a.tss_tol) for mdl in te_models)
            te_found = te_expressed and any(rep_te_ok(L) for L in same)
            row[f"{arm}_te_found"] = int(te_found)
            summary["arms"][arm]["te_found"] = summary["arms"][arm].get("te_found", 0) + te_found
            summary["arms"][arm]["old_overlap"] += old
            summary["arms"][arm]["strict_found"] += strict
            # locus level (reported beside): any transcript of a same-strand overlapping locus carries >= k supported junctions
            ltx = arm_tx.get(arm)
            locus_found = None
            if ltx is not None:
                names_same = [n for n, L in loci.items() if L in same]
                best = max((len(ltx.get(n, set()) & supported) for n in names_same), default=0)
                row[f"{arm}_locus_supported_junctions_max"] = best
                locus_found = expressed and (best >= k if k else len(strict_loci) > 0)
                row[f"{arm}_locus_level_found"] = int(locus_found)
                summary["arms"][arm]["locus_level_found"] += locus_found
                best_ann = max((len(ltx.get(n, set()) & ann_all) for n in names_same), default=0)
                locus_ann_found = ann_expressed and ((len(strict_loci) > 0) if k == 0 else best_ann >= k)
                row[f"{arm}_locus_ann_junctions_max"] = best_ann
                row[f"{arm}_locus_ann_found"] = int(locus_ann_found)
                summary["arms"][arm]["locus_ann_found"] = summary["arms"][arm].get("locus_ann_found", 0) + locus_ann_found
                # Amendment B at the locus level: any transcript of a same-strand overlapping locus chain-matches a model
                if kb == 0:
                    locus_chain_found = chain_expressed and len(chain_loci) > 0
                else:
                    locus_chain_found = chain_expressed and any(
                        chain_match(sorted(j for j in ltx_t if j[0] >= lo and j[1] <= hi), chains, kb)
                        for n in names_same for ltx_t in [arm_tx_chains.get(arm, {}).get(n, [])] for ltx_t in ltx_t)
                row[f"{arm}_locus_chain_found"] = int(locus_chain_found)
                summary["arms"][arm]["locus_chain_found"] = summary["arms"][arm].get("locus_chain_found", 0) + locus_chain_found
                if not kb:
                    locus_xc_found = xc_found
                else:
                    locus_xc_found = xc_expressed and any(
                        chain_match(sorted(j for j in t if j[0] >= lo and j[1] <= hi), expressed_chains, 2)
                        for n in names_same for t in arm_tx_chains.get(arm, {}).get(n, []))
                row[f"{arm}_locus_xc_found"] = int(locus_xc_found)
                summary["arms"][arm]["locus_xc_found"] = summary["arms"][arm].get("locus_xc_found", 0) + locus_xc_found
                # Amendment D at the locus level: any transcript of the locus starts at a model's TSS and carries its first m introns
                def tx_td_ok(exs):
                    ex = sorted(exs)
                    p5 = ex[-1][1] if strand == "-" else ex[0][0]
                    js = introns_of(ex)
                    return any(tss_support(p5, js, mdl, a.tss_tol) and (len(mdl[1]) > 0 or inter(ex, mdl[2]) >= COVER * ilen(mdl[2])) for mdl in models_d)
                locus_td_found = td_expressed and any(tx_td_ok(exs) for n in names_same for exs in arm_tx_exons.get(arm, {}).get(n, []))
                row[f"{arm}_locus_td_found"] = int(locus_td_found)
                summary["arms"][arm]["locus_td_found"] = summary["arms"][arm].get("locus_td_found", 0) + locus_td_found
                def tx_te_ok(exs):
                    ex = sorted(exs)
                    p5 = ex[-1][1] if strand == "-" else ex[0][0]
                    return any(tss_support(p5, introns_of(ex), mdl, a.tss_tol) for mdl in te_models)
                locus_te_found = te_expressed and any(tx_te_ok(exs) for n in names_same for exs in arm_tx_exons.get(arm, {}).get(n, []))
                row[f"{arm}_locus_te_found"] = int(locus_te_found)
                summary["arms"][arm]["locus_te_found"] = summary["arms"][arm].get("locus_te_found", 0) + locus_te_found
            if nodes:
                own = bool(nodes.get(c["cid"], {}).get(arm, False))
                row[f"{arm}_page_own_node"] = int(own)
                summary["arms"][arm]["old_overlap_in_npip_nodes"] += own
                summary["arms"][arm]["strict_found_in_npip_nodes"] += own and strict
                summary["arms"][arm]["ann_found_in_npip_nodes"] = summary["arms"][arm].get("ann_found_in_npip_nodes", 0) + (own and ann_found)
                summary["arms"][arm]["chain_found_in_npip_nodes"] = summary["arms"][arm].get("chain_found_in_npip_nodes", 0) + (own and chain_found)
                summary["arms"][arm]["xc_found_in_npip_nodes"] = summary["arms"][arm].get("xc_found_in_npip_nodes", 0) + (own and xc_found)
                summary["arms"][arm]["td_found_in_npip_nodes"] = summary["arms"][arm].get("td_found_in_npip_nodes", 0) + (own and td_found)
                summary["arms"][arm]["te_found_in_npip_nodes"] = summary["arms"][arm].get("te_found_in_npip_nodes", 0) + (own and te_found)
                if locus_found is not None:
                    summary["arms"][arm]["locus_level_found_in_npip_nodes"] += own and locus_found
        rows.append(row)
    with open(a.out + ".copies.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    json.dump(summary, open(a.out + ".json", "w"), indent=1)
    print(json.dumps(summary))
    hdr = ["name", "reads", "td_support_reads", "td_support_50", "td_support_300", "td_support_cat", "td_support_refseq", "td_support_unique", "td_expressed", "xc_support_reads"] + \
          [f"{arm}_{x}" for arm in arms for x in ("td_found", "xc_found") + (("locus_td_found",) if arm in arm_tx else ()) + (("page_own_node",) if nodes else ())]
    print("\t".join(hdr))
    for r in rows:
        print("\t".join(str(r.get(h, "")) for h in hdr))


if __name__ == "__main__":
    main()
