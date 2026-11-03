"""r1513 (docs/PREREG_secondary_pool_2026-11-03.md): flagged candidate contigs of a candidate_copies arm on Jim's haplotypes.
Usage: pool_score.py ARM_PREFIX HITS_TSV MAT.mmi PAT.mmi OUT_PREFIX   -> OUT_PREFIX.contigs.tsv + summary on stdout.
Class per flagged contig (best haplotype hit by matches): absent_from_ref (identity >= 0.998, coverage >= 0.9, outside every
catalog copy's haplotype locus from HITS_TSV `C|family|copy` rows at g = 0), allele_of_known_copy (>= 0.998 on a copy locus),
not_in_genome (< 0.995), else `between` (0.995-0.998: reported, no class)."""
import sys, os, subprocess, collections, tempfile
arm, hits_p, mat, pat, out = sys.argv[1:6]
flagged = set()
with open(arm + ".candidates.tsv") as fh:
    h = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(h, line.rstrip("\n").split("\t")))
        if r["flagged"] == "1": flagged.add((r["family"], r["candidate"]))
seqs, name = {}, None
for line in open(arm + ".contigs.fa"):
    if line.startswith(">"): name = line[1:].split()[0]
    else: seqs[name] = seqs.get(name, "") + line.strip()
fam_of = {c: f for f, c in flagged}
loci = collections.defaultdict(set)   # haplotype -> {(chrom, s, e, family)}
best = {}
rows = collections.defaultdict(list)
with open(hits_p) as fh:
    next(fh)
    for line in fh:
        q, hap, ch, s, e, nm, al = line.rstrip("\n").split("\t")
        if q.startswith("C|"): rows[q].append((hap, ch, int(s), int(e), int(nm)))
for q, hs in rows.items():
    _, f, _ = q.split("|")
    bst = {}
    for h in hs: bst[h[0]] = min(bst.get(h[0], 10 ** 9), h[4])
    for hap, ch, s, e, nm in hs:
        if nm <= bst[hap]: loci[hap].add((ch, s, e, f))
with tempfile.TemporaryDirectory() as d:
    fa = os.path.join(d, "c.fa")
    with open(fa, "w") as w:
        for c in sorted(fam_of): w.write(f">{c}\n{seqs[c]}\n")
    for hap, mmi in (("mat", mat), ("pat", pat)):
        if not fam_of: break
        paf = subprocess.run(["minimap2", "-c", "-x", "splice:hq", "-uf", "--secondary=no", "-t", "8", mmi, fa], capture_output=True, text=True).stdout
        for line in paf.splitlines():
            f = line.split("\t")
            q, ql, qs, qe, ch, ts, te, m, bl = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
            idn, cov = m / bl if bl else 0.0, (qe - qs) / ql if ql else 0.0
            if q not in best or m > best[q][0]: best[q] = (m, hap, ch, ts, te, idn, cov)
def on_copy(hap, ch, s, e, fam):
    return any(c == ch and s < le and ls < e for (c, ls, le, f) in loci[hap])
n = collections.Counter()
with open(out + ".contigs.tsv", "w") as w:
    w.write("family\tcandidate\tlen\thap\tchrom\tstart\tend\tidentity\tcoverage\ton_copy_locus\tclass\n")
    for c in sorted(fam_of):
        fam = fam_of[c]
        if c not in best:
            cls, rec = "not_in_genome", ("NA", "NA", 0, 0, 0.0, 0.0); oc = False
        else:
            _, hap, ch, ts, te, idn, cov = best[c]; rec = (hap, ch, ts, te, idn, cov); oc = on_copy(hap, ch, ts, te, fam)
            cls = "not_in_genome" if idn < 0.995 else "between" if idn < 0.998 or cov < 0.9 else "allele_of_known_copy" if oc else "absent_from_ref"
        n[cls] += 1
        w.write(f"{fam}\t{c}\t{len(seqs[c])}\t{rec[0]}\t{rec[1]}\t{rec[2]}\t{rec[3]}\t{rec[4]:.4f}\t{rec[5]:.3f}\t{int(oc)}\t{cls}\n")
print(f"{os.path.basename(arm)}\tflagged {len(fam_of)}\t" + "\t".join(f"{k} {n[k]}" for k in ("absent_from_ref", "allele_of_known_copy", "not_in_genome", "between")))
