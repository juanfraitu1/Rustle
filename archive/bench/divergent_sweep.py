#!/usr/bin/env python3
"""divergent_sweep.py — re-refine the RAW conflict catalog with a configurable homology detector, to ask
"can a looser/different aligner include MORE divergent paralog families, and at what precision cost?"

Refinement logic is identical to the shipped one (homology connected-component AND >=2 spatially-DISTINCT
loci); only the EDGE detector changes:
  --method asm20            minimap2 -cx asm20            (baseline, identity gate --id)
  --method ont              minimap2 -cx map-ont
  --method sensitive        minimap2 -c -k11 -w5 -g5000 -A1 -B2 -O2,32 -E1,0 (smaller seeds = more divergent)
  --method protein          longest-ORF 6-frame translate exon-sum -> mmseqs all-vs-all (fident gate --id)

Per refined family it reports the levers that separate a real divergent paralog from a repeat/domain merge:
  min_identity   (lowest passing pair identity = how divergent the family is)
  mean_aln_frac  (aligned length / shorter length; LOW = the copies share only a fragment = repeat/domain)
Writes <out>.tsv (one row per refined family). Run on the RAW catalog (gw_xchrom_raw).
"""
import sys, csv, subprocess, collections, os, tempfile, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--raw", default="/home/juanfra/winloci_scratch/gw_xchrom_raw")
ap.add_argument("--method", default="asm20", choices=["asm20", "ont", "sensitive", "protein"])
ap.add_argument("--id", type=float, default=0.80)
ap.add_argument("--cov", type=float, default=0.50)
ap.add_argument("--mm2", default="/home/juanfra/miniforge3/bin/minimap2")
ap.add_argument("--mmseqs", default="mmseqs")
ap.add_argument("--out", default="/home/juanfra/winloci_scratch/divsweep")
A = ap.parse_args()

# raw families: exon-sum seqs grouped by family, plus loci
fam_seq = collections.defaultdict(dict)   # fid -> {ci: seq}
fam_loc = collections.defaultdict(dict)   # fid -> {ci: (chrom,start,end,strand)}
fid = ci = None; buf = []
def flush():
    if fid is not None and buf: fam_seq[fid][ci] = "".join(buf)
for line in open(f"{A.raw}.copies.fa"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        flush()
        h = line[1:].split("|"); fid, ci = h[0], h[1]
        chrom, span = h[2].rsplit(":", 1); s, e = span.split("-")
        fam_loc[fid][ci] = (chrom, int(s), int(e), h[3]); buf = []
    else:
        buf.append(line)
flush()
fams = {f["family_id"]: f for f in csv.DictReader(open(f"{A.raw}.families.tsv"), delimiter="\t")}

CODON = {}
_bases = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, c1 in enumerate(_bases):
    for j, c2 in enumerate(_bases):
        for k, c3 in enumerate(_bases):
            CODON[c1 + c2 + c3] = _aa[i * 16 + j * 4 + k]
def revcomp(s):
    return s.translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))[::-1]
def longest_orf_aa(seq):
    seq = seq.upper()
    best = ""
    for st in (seq, revcomp(seq)):
        for fr in range(3):
            aa = "".join(CODON.get(st[i:i+3], "X") for i in range(fr, len(st) - 2, 3))
            for part in aa.split("*"):
                if len(part) > len(best):
                    best = part
    return best

def edges_nucleotide(seqs, mm_args):
    cis = list(range(len(seqs)))
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        p = fh.name
        for i in cis: fh.write(f">{i}\n{seqs[i]}\n")
    try:
        out = subprocess.run([A.mm2, *mm_args, "-X", "--no-long-join", p, p],
                             capture_output=True, text=True).stdout
    finally:
        os.unlink(p)
    lens = [len(s) for s in seqs]
    edges = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 12: continue
        try: q, t = int(f[0]), int(f[5])
        except ValueError: continue
        if q == t: continue
        qs, qe = int(f[2]), int(f[3])
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        ident = (1 - de) if de is not None else int(f[9]) / max(1, int(f[10]))
        frac = (qe - qs) / max(1, min(lens[q], lens[t]))
        if ident >= A.id and frac >= A.cov:
            k = (min(q, t), max(q, t))
            if k not in edges or ident > edges[k][0]: edges[k] = (ident, frac)
    return edges

def edges_protein(seqs):
    prots = [longest_orf_aa(s) for s in seqs]
    cis = [i for i, p in enumerate(prots) if len(p) >= 40]
    if len(cis) < 2: return {}
    d = tempfile.mkdtemp()
    qf = os.path.join(d, "q.fa")
    with open(qf, "w") as fh:
        for i in cis: fh.write(f">{i}\n{prots[i]}\n")
    res = os.path.join(d, "res.m8")
    subprocess.run([A.mmseqs, "easy-search", qf, qf, res, os.path.join(d, "tmp"),
                    "-s", "7.5", "--format-output", "query,target,fident,qcov,tcov,evalue"],
                   capture_output=True, text=True)
    edges = {}
    if os.path.exists(res):
        for ln in open(res):
            f = ln.split("\t")
            if len(f) < 6: continue
            q, t, fident, qcov, tcov = int(f[0]), int(f[1]), float(f[2]), float(f[3]), float(f[4])
            if q == t: continue
            if fident >= A.id and min(qcov, tcov) >= A.cov:
                k = (min(q, t), max(q, t))
                if k not in edges or fident > edges[k][0]: edges[k] = (fident, min(qcov, tcov))
    subprocess.run(["rm", "-rf", d])
    return edges

MM = {"asm20": ["-cx", "asm20"], "ont": ["-cx", "map-ont"],
      "sensitive": ["-c", "-k", "11", "-w", "5", "-g", "5000", "-A", "1", "-B", "2", "-O", "2,32", "-E", "1,0"]}

def refine(fid):
    cis = sorted(fam_seq[fid], key=int)
    seqs = [fam_seq[fid][c] for c in cis]
    if len(seqs) < 2: return []
    edges = edges_protein(seqs) if A.method == "protein" else edges_nucleotide(seqs, MM[A.method])
    # components
    par = list(range(len(cis)))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for (i, j) in edges: par[find(i)] = find(j)
    comps = collections.defaultdict(list)
    for i in range(len(cis)): comps[find(i)].append(i)
    out = []
    for comp in comps.values():
        if len(comp) < 2: continue
        # distinct loci (same chrom+strand overlap OR minority-antisense collapse omitted; coordinate-disjoint)
        locs = [fam_loc[fid][cis[i]] for i in comp]
        keep = []
        for li in range(len(locs)):
            dup = any(locs[li][0] == locs[lj][0] and min(locs[li][2], locs[lj][2]) > max(locs[li][1], locs[lj][1])
                      for lj in keep)
            if not dup: keep.append(li)
        if len(keep) < 2: continue
        idents = [edges[(min(comp[a], comp[b]), max(comp[a], comp[b]))][0]
                  for a in range(len(comp)) for b in range(a + 1, len(comp))
                  if (min(comp[a], comp[b]), max(comp[a], comp[b])) in edges]
        fracs = [edges[(min(comp[a], comp[b]), max(comp[a], comp[b]))][1]
                 for a in range(len(comp)) for b in range(a + 1, len(comp))
                 if (min(comp[a], comp[b]), max(comp[a], comp[b])) in edges]
        loci_str = ";".join(f"{fam_loc[fid][cis[i]][0]}:{fam_loc[fid][cis[i]][1]}-{fam_loc[fid][cis[i]][2]}" for i in comp)
        chroms = set(fam_loc[fid][cis[i]][0] for i in comp)
        out.append((len([k for k in keep]), len(chroms) > 1, min(idents) if idents else 0,
                    sum(fracs)/len(fracs) if fracs else 0, loci_str))
    return out

refined = []
for fid in fam_seq:
    refined += [(fid,) + r for r in refine(fid)]
xc = sum(1 for r in refined if r[2])
with open(f"{A.out}.tsv", "w") as fh:
    fh.write("orig_family\tn_loci\tcross_chrom\tmin_identity\tmean_aln_frac\tloci\n")
    for r in sorted(refined, key=lambda x: x[3]):
        fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.3f}\t{r[4]:.3f}\t{r[5]}\n")
print(f"method={A.method} id>={A.id} cov>={A.cov}: {len(refined)} refined families ({xc} cross-chrom) -> {A.out}.tsv")
# divergence histogram of the families this method admits
import statistics
ids = [r[3] for r in refined]
if ids:
    print(f"  family min-identity: median={statistics.median(ids):.3f} min={min(ids):.3f} "
          f"(<0.85: {sum(1 for x in ids if x<0.85)}, <0.75: {sum(1 for x in ids if x<0.75)})")
    lowfrac = sum(1 for r in refined if r[3] < 0.30)
    print(f"  families with mean_aln_frac<0.30 (repeat/domain-share risk): {lowfrac}")
