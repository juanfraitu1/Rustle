"""DNA-derived PSV identifiability catalog (Phase 1, reference-only). See
docs/superpowers/specs/2026-06-21-dna-psv-catalog-design.md.

Per family: extract each copy's T2T interval, align every copy to the longest (ref0) with
minimap2 asm20 --cs, project substitutions onto ref0 coordinates -> allele matrix -> exonic PSV columns ->
per-pair K (Hamming over exonic PSVs). Cross-check vs the RNA census. Deterministic check() at the bottom."""
import collections, os, re, subprocess, tempfile
import pysam

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
    """minimap2 asm20 --cs=long: return substitution columns as (ref_localpos, qry_base), 0-based on ref_seq."""
    with tempfile.TemporaryDirectory() as td:
        rp, qp = os.path.join(td, "r.fa"), os.path.join(td, "q.fa")
        open(rp, "w").write(">r\n" + ref_seq + "\n"); open(qp, "w").write(">q\n" + qry_seq + "\n")
        out = subprocess.run(["minimap2", "-cx", "asm20", "--cs=long", "-t", "1", rp, qp],
                             capture_output=True, text=True).stdout
    subs = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 9: continue
        rstart = int(f[7])                          # target (ref) start, 0-based
        cs = next((x[5:] for x in f if x.startswith("cs:Z:")), None)
        if cs is None: continue
        rpos = rstart
        for op, val in re.findall(r"([:=*+\-])([A-Za-z0-9]+)", cs):
            if op == ":":            rpos += int(val)                 # identical run
            elif op == "=":          rpos += len(val)                 # identical (long form)
            elif op == "*":          subs.append((rpos, val[1].upper())); rpos += 1   # *<ref><qry> substitution
            elif op == "-":          rpos += len(val)                 # deletion from query (consume ref)
            elif op == "+":          pass                             # insertion in query (no ref consume)
    return subs

def family_allele_matrix(members):
    """members = [(label, chrom, start)]. Align all copies to the longest (ref0); return
    (ref0_label, ref0_chrom, ref0_iv, matrix, intervals) with matrix[label][ref0_localpos]=base."""
    ivs = {}
    for label, chrom, start in members:
        s, e = copy_interval(chrom, start); ivs[label] = (chrom, s, e)
    ref0 = max(ivs, key=lambda L: ivs[L][2] - ivs[L][1])
    rc, rs, re_ = ivs[ref0]; ref_seq = extract_seq(rc, rs, re_)
    matrix = collections.defaultdict(dict)        # label -> {ref0_localpos: base}
    for label, (c, s, e) in ivs.items():
        if label == ref0: continue
        subs = align_to_ref(ref_seq, extract_seq(c, s, e))
        for pos, base in subs:
            matrix[label][pos] = base
    return ref0, rc, (rs, re_), matrix, ivs

if __name__ == "__main__":
    # smoke: two MAGEA copies expected to differ (MAGd1 pair)
    members = [("A", "NC_073247.2", 161266095), ("B", "NC_073247.2", 161413981)]
    ref0, rc, riv, matrix, ivs = family_allele_matrix(members)
    nonref = [L for L in ivs if L != ref0][0]
    print(f"ref0={ref0} {rc}:{riv[0]}-{riv[1]} ({riv[1]-riv[0]}bp); "
          f"substitution columns {nonref} vs ref0: {len(matrix[nonref])}")
