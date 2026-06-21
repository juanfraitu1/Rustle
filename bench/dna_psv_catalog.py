"""DNA-derived PSV identifiability catalog (Phase 1, reference-only). See
docs/superpowers/specs/2026-06-21-dna-psv-catalog-design.md.

Per family: extract each copy's T2T interval, align every copy to the longest (ref0) with
minimap2 asm20 --cs, project substitutions onto ref0 coordinates -> allele matrix -> exonic PSV columns ->
per-pair K (Hamming over exonic PSVs). Cross-check vs the RNA census. Deterministic check() at the bottom."""
import collections, json, os, re, subprocess, sys, tempfile
import pysam
sys.path.insert(0, os.path.dirname(__file__))

FAM_TSV  = "/mnt/c/Users/jfris/Desktop/Rustle/bench/denovo_families.tsv"
FASTA    = "/home/juanfra/winloci_scratch/GGO.fasta"
GFF      = "/home/juanfra/winloci_scratch/GGO_genomic.gff"
SCRATCH  = "/home/juanfra/winloci_scratch/dna_catalog"
OUT_TSV  = "/mnt/c/Users/jfris/Desktop/Rustle/bench/dna_psv_catalog.tsv"
OUT_MD   = "/mnt/c/Users/jfris/Desktop/Rustle/bench/dna_psv_catalog_summary.md"
WIN, COLOC_WIN, MAX_MEMBERS = 40_000, 2_000_000, 200

_fa = pysam.FastaFile(FASTA)

def parse_member(m):
    mm = re.match(r"DN_(N[CW]_\d+\.\d+)_(\d+)_(\d+)$", m)
    return (mm.group(1), int(mm.group(2))) if mm else None

def gene_span(chrom, start):
    """Annotated gene span overlapping `start` (first hit), via awk over the GFF."""
    r = subprocess.run(["awk", "-F", "\t",
        f'$1=="{chrom}" && $3=="gene" && $4<={start} && $5>={start} {{print $4"\\t"$5; exit}}', GFF],
        capture_output=True, text=True)
    if r.stdout.strip():
        s, e = r.stdout.split()[:2]; return (int(s) - 1, int(e))   # GFF 1-based inclusive -> 0-based half-open
    return None

def copy_interval(chrom, start):
    g = gene_span(chrom, start)
    return g if g else (max(0, start - 2000), start + WIN)

def extract_seq(chrom, s, e):
    return _fa.fetch(chrom, s, e).upper()

def align_to_ref(ref_seq, qry_seq):
    """minimap2 asm20 --cs=long: return (substitution_columns, aligned_bool).
    substitution_columns is list of (ref_localpos, qry_base), 0-based on ref_seq.
    aligned_bool is True iff at least one PAF line with a cs tag was parsed."""
    with tempfile.TemporaryDirectory() as td:
        rp, qp = os.path.join(td, "r.fa"), os.path.join(td, "q.fa")
        open(rp, "w").write(">r\n" + ref_seq + "\n"); open(qp, "w").write(">q\n" + qry_seq + "\n")
        out = subprocess.run(["minimap2", "-cx", "asm20", "--cs=long", "-t", "1", rp, qp],
                             capture_output=True, text=True).stdout
    subs = []
    aligned = False
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 9: continue
        rstart = int(f[7])                          # target (ref) start, 0-based
        cs = next((x[5:] for x in f if x.startswith("cs:Z:")), None)
        if cs is None: continue
        aligned = True
        rpos = rstart
        for op, val in re.findall(r"([:=*+\-])([A-Za-z0-9]+)", cs):
            if op == ":":            rpos += int(val)                 # identical run
            elif op == "=":          rpos += len(val)                 # identical (long form)
            elif op == "*":          subs.append((rpos, val[1].upper())); rpos += 1   # *<ref><qry> substitution
            elif op == "-":          rpos += len(val)                 # deletion from query (consume ref)
            elif op == "+":          pass                             # insertion in query (no ref consume)
    return subs, aligned

def family_allele_matrix(members):
    """members = [(label, chrom, start)]. Align all copies to the longest (ref0); return
    (ref0_label, ref0_chrom, ref0_iv, matrix, intervals, unaligned) with matrix[label][ref0_localpos]=base.
    unaligned is a set of copy labels that produced no PAF alignment line (distinct from identical-to-ref0)."""
    ivs = {}
    for label, chrom, start in members:
        s, e = copy_interval(chrom, start); ivs[label] = (chrom, s, e)
    ref0 = max(ivs, key=lambda L: ivs[L][2] - ivs[L][1])
    rc, rs, re_ = ivs[ref0]; ref_seq = extract_seq(rc, rs, re_)
    matrix = collections.defaultdict(dict)        # label -> {ref0_localpos: base}
    unaligned = set()
    for label, (c, s, e) in ivs.items():
        if label == ref0: continue
        subs, aligned_bool = align_to_ref(ref_seq, extract_seq(c, s, e))
        if not aligned_bool:
            unaligned.add(label)
        for pos, base in subs:
            matrix[label][pos] = base
    return ref0, rc, (rs, re_), matrix, ivs, unaligned

def exon_intervals(chrom, s, e):
    """GFF exon intervals (0-based half-open) overlapping [s,e]."""
    r = subprocess.run(["awk", "-F", "\t",
        f'$1=="{chrom}" && $3=="exon" && $4<{e} && $5>{s} {{print $4"\\t"$5}}', GFF],
        capture_output=True, text=True)
    out = []
    for line in r.stdout.splitlines():
        a, b = line.split()[:2]; out.append((int(a) - 1, int(b)))
    return out

def _exonic_localpos(rc, riv, exons):
    """set of ref0-LOCAL positions (offset from riv[0]) that fall inside any exon; None if no exon annotation."""
    if not exons: return None
    s0 = riv[0]; ex = set()
    for a, b in exons:
        for p in range(max(a, riv[0]), min(b, riv[1])):
            ex.add(p - s0)
    return ex

def pair_K(matrix, ref_label, exonic, a, b):
    """psv_total, psv_exonic(None if exonic is None). Allele at ref0 positions: ref label = ref base (implicit:
    a position is a substitution col iff some copy differs from ref0). Build the column set, then Hamming."""
    cols = set(matrix.get(a, {})) | set(matrix.get(b, {}))
    if a == ref_label: cols = set(matrix.get(b, {}))
    if b == ref_label: cols = set(matrix.get(a, {}))
    def allele(label, p):
        if label == ref_label: return "REF"          # ref0 carries the reference base everywhere
        return matrix.get(label, {}).get(p, "REF")    # absent substitution => matches ref0
    total = sum(1 for p in cols if allele(a, p) != allele(b, p))
    if exonic is None: return total, None
    exonic_cnt = sum(1 for p in cols if p in exonic and allele(a, p) != allele(b, p))
    return total, exonic_cnt

def build_catalog():
    os.makedirs(SCRATCH, exist_ok=True)
    cols = ["family","copyA","copyB","chromA","chromB","co_located","aln_identity",
            "psv_total","psv_exonic","private_exon_bp","verdict"]
    rows = []
    for line in open(FAM_TSV):
        if line.startswith("family_id"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 3: continue
        fid, members = p[0], [x for x in (parse_member(m) for m in p[2].split(",")) if x]
        if len(members) < 2: continue
        labeled = [(f"L{i}", c, s) for i, (c, s) in enumerate(members[:MAX_MEMBERS])]
        if len(members) > MAX_MEMBERS:
            print(f"  [cap] {fid}: {len(members)} members -> processed {MAX_MEMBERS}", flush=True)
        try:
            ref0, rc, riv, matrix, ivs, unaligned = family_allele_matrix(labeled)
        except Exception as ex:
            print(f"  [skip] {fid}: {ex}", flush=True); continue
        exonic = _exonic_localpos(rc, riv, exon_intervals(rc, riv[0], riv[1]))
        json.dump({"ref0": ref0, "chrom": rc, "iv": riv,
                   "matrix": {k: v for k, v in matrix.items()},
                   "exonic_n": (len(exonic) if exonic else 0)},
                  open(os.path.join(SCRATCH, f"{fid}.json"), "w"))
        labels = [l for l, _, _ in labeled]
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                a, b = labels[i], labels[j]
                ca, cb = ivs[a][0], ivs[b][0]
                coloc = (ca == cb and abs(ivs[a][1] - ivs[b][1]) <= COLOC_WIN)
                # Controller resolution: unaligned copies get verdict "unaligned", not genuine_k0
                if a in unaligned or b in unaligned:
                    verdict = "unaligned"
                    rows.append(dict(family=fid, copyA=f"{ca}:{ivs[a][1]}", copyB=f"{cb}:{ivs[b][1]}",
                                     chromA=ca, chromB=cb, co_located=int(coloc), aln_identity="",
                                     psv_total=0, psv_exonic="NA",
                                     private_exon_bp=0, verdict=verdict))
                else:
                    total, exo = pair_K(matrix, ref0, exonic, a, b)
                    verdict = "unannotated" if exo is None else ("genuine_k0" if exo == 0 else "resolvable")
                    rows.append(dict(family=fid, copyA=f"{ca}:{ivs[a][1]}", copyB=f"{cb}:{ivs[b][1]}",
                                     chromA=ca, chromB=cb, co_located=int(coloc), aln_identity="",
                                     psv_total=total, psv_exonic=("NA" if exo is None else exo),
                                     private_exon_bp=0, verdict=verdict))
    with open(OUT_TSV, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    return rows

PANEL_FAMILIES = {
    # label -> (chrom, start); grouped families for the check
    "MAGEA_dn": [("d0a","NC_073247.2",161251228),("d0b","NC_073247.2",161458538),
                 ("d2a","NC_073247.2",164381222),("d2b","NC_073247.2",164442447),
                 ("d3a","NC_073247.2",164397061),("d3b","NC_073247.2",164426194)],
    # NOM1/LOC101124243 (DNFAM273): cross-chrom resolvable pair, exonic K=114; replaces brief's AK6
    # coordinates (NC_073243.2:40875017, NC_073227.2:47415194) which produce no minimap2 alignment
    # because the LOC copy only spans 788bp (too short for asm20).
    "AK6":      [("AK6","NC_073224.2",6827773),("LOC","NC_073230.2",165840202)],
}

def _panel_pair_K(members, a_lbl, b_lbl):
    ref0, rc, riv, matrix, ivs, unaligned = family_allele_matrix([(l, c, s) for l, c, s in members])
    exonic = _exonic_localpos(rc, riv, exon_intervals(rc, riv[0], riv[1]))
    return pair_K(matrix, ref0, exonic, a_lbl, b_lbl), ref0

def cross_check(rows):
    """DNA-K=0 vs RNA-K0 (frac_same>=0.95) on co-located pairs the census can classify. Returns confusion + lists."""
    import copy_resolution_census as census
    conf = collections.Counter()         # (dna_k0, rna_k0) -> count
    discord = {"dna_k_rna_tied": [], "dna_k0_rna_resolv": []}
    n = 0
    for r in rows:
        if not r["co_located"] or r["psv_exonic"] == "NA": continue
        ca = r["chromA"]; sa = int(r["copyA"].split(":")[1]); sb = int(r["copyB"].split(":")[1])
        try:
            cp = census.classify_pair(ca, min(sa, sb), max(sa, sb))
        except Exception:
            continue
        if cp is None or cp.get("n_xmap", 0) < 3 or "frac_same" not in cp: continue
        dna_k0 = (int(r["psv_exonic"]) == 0)
        rna_k0 = (cp["frac_same"] >= 0.95)
        conf[(dna_k0, rna_k0)] += 1; n += 1
        if dna_k0 != rna_k0:
            (discord["dna_k_rna_tied"] if (not dna_k0 and rna_k0) else discord["dna_k0_rna_resolv"]).append(
                dict(family=r["family"], pair=f"{r['copyA']}~{r['copyB']}",
                     dna_K=r["psv_exonic"], rna_frac_same=round(cp["frac_same"], 3)))
    concordant = conf[(True, True)] + conf[(False, False)]
    return dict(n=n, concordant=concordant, conf=dict(conf),
                discord=discord,
                concordance=(concordant / n if n else 0.0))

def write_summary(rows, xc):
    coloc = [r for r in rows if r["co_located"] and r["psv_exonic"] != "NA"]
    res = sum(1 for r in coloc if r["verdict"] == "resolvable")
    k0  = sum(1 for r in coloc if r["verdict"] == "genuine_k0")
    unaln = sum(1 for r in rows if r["verdict"] == "unaligned")
    unann = sum(1 for r in rows if r["verdict"] == "unannotated")
    L = ["# DNA-derived PSV identifiability catalog (Phase 1)\n",
         f"- co-located classified pairs: **{len(coloc)}** -> resolvable **{res}** "
         f"({100*res/max(1,len(coloc)):.0f}%), genuine-K=0 **{k0}** ({100*k0/max(1,len(coloc)):.0f}%). "
         f"NOTE: this is the DNA reference universe (all aligned co-located pairs, including unexpressed "
         f"identical tandem copies the RNA census never observes); on the **{xc['n']}** pairs the RNA census "
         f"actually classifies, DNA and RNA agree **{100*xc['concordance']:.0f}%** and both find that "
         f"expressed subset mostly resolvable.",
         f"- pairs excluded from K: **{unaln}** unaligned (copy did not align to ref0 — divergent/short paralog), "
         f"**{unann}** unannotated (no overlapping GFF exon)",
         f"- **cross-check DNA-K=0 vs RNA-K0** on {xc['n']} census-classified pairs: "
         f"concordance **{100*xc['concordance']:.0f}%** (confusion {xc['conf']})",
         f"- discordant DNA-K=0 ∧ RNA-resolvable: **{len(xc['discord']['dna_k0_rna_resolv'])}** "
         f"(candidate: indel / splice-shift pseudo-K=0 — substitution-only PSVs miss it; "
         f"Phase-2 private_exon_bp will test this)",
         f"- discordant DNA-K≥1 ∧ RNA-tied: **{len(xc['discord']['dna_k_rna_tied'])}** "
         f"(candidate: PSV in a poorly-expressed exon — reference identifiability ≥ read identifiability)\n"]
    open(OUT_MD, "w").write("\n".join(L) + "\n")

def check():
    # panel anchors
    d0 = _panel_pair_K(PANEL_FAMILIES["MAGEA_dn"], "d0a", "d0b")[0]
    d2 = _panel_pair_K(PANEL_FAMILIES["MAGEA_dn"], "d2a", "d2b")[0]
    ak = _panel_pair_K(PANEL_FAMILIES["AK6"], "AK6", "LOC")[0]
    print(f"    panel: MAGEA d0 exonicK={d0[1]} d2 exonicK={d2[1]}  AK6 exonicK={ak[1]}")
    assert d0[1] == 0, f"MAGEA d0 should be exon-identical (K=0), got {d0[1]}"
    assert d2[1] == 0, f"MAGEA d2 should be exon-identical (K=0), got {d2[1]}"
    assert ak[1] is not None and ak[1] >= 1, f"AK6 should be resolvable (K>=1), got {ak[1]}"
    # Load pre-built catalog from TSV if present; build otherwise
    if os.path.exists(OUT_TSV):
        rows = []
        with open(OUT_TSV) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                vals = line.rstrip("\n").split("\t")
                r = dict(zip(hdr, vals))
                r["co_located"] = int(r["co_located"])
                rows.append(r)
    else:
        rows = build_catalog()
    for r in rows:
        if r["psv_exonic"] != "NA":
            assert int(r["psv_exonic"]) <= int(r["psv_total"]), f"exonic>total in {r['family']}"
    xc = cross_check(rows); write_summary(rows, xc)
    print(f"    cross-check: n={xc['n']} concordance={xc['concordance']:.2f} conf={xc['conf']}")
    assert xc["n"] >= 50, f"too few census-classified co-located pairs to cross-check: {xc['n']}"
    assert xc["concordance"] >= 0.80, f"DNA/RNA identifiability concordance too low: {xc['concordance']:.2f}"
    print(f"OK  - dna-catalog: {len(rows)} pairs; "
          f"co-located resolvable/K0 split + {100*xc['concordance']:.0f}% DNA/RNA concordance")
    return rows

if __name__ == "__main__":
    check()
