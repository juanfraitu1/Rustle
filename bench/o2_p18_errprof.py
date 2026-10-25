# P18 scorer parameters: per-event error rates of real fibroblast FLNC reads (MAPQ 60 primaries; real allele differences included)
import pysam, collections, json, sys
BAM = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
FA = pysam.FastaFile("/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta")
CHROMS = ["NC_073224.2", "NC_073227.2", "NC_073229.2", "NC_073233.2"]
N_PER = 6000
def hp(ref, p, b):
    "ref run length of base b touching position p (indel between p-1 and p)"
    n = 0; i = p
    while i < len(ref) and ref[i] == b: n += 1; i += 1
    i = p - 1
    while i >= 0 and ref[i] == b: n += 1; i -= 1
    return n
cnt = collections.Counter()
bam = pysam.AlignmentFile(BAM)
for c in CHROMS:
    k = 0
    for r in bam.fetch(c):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or r.mapping_quality < 60: continue
        ref = FA.fetch(c, r.reference_start, r.reference_end).upper(); q = r.query_sequence.upper()
        rp, qp = 0, 0
        for op, L in r.cigartuples:
            if op == 7: cnt["M"] += L; rp += L; qp += L
            elif op == 8: cnt["X"] += L; cnt["M"] += 0; rp += L; qp += L
            elif op == 1:
                s = q[qp:qp + L]; ctx = "hp" if len(set(s)) == 1 and hp(ref, rp, s[0]) >= 1 else "nhp"
                cnt[f"I_{ctx}_{min(L, 3)}"] += 1; qp += L
            elif op == 2:
                s = ref[rp:rp + L]; ctx = "hp" if len(set(s)) == 1 and hp(ref, rp, s[0]) >= L + 1 else "nhp"
                cnt[f"D_{ctx}_{min(L, 3)}"] += 1; rp += L
            elif op == 3: rp += L
            elif op == 4: qp += L
        k += 1
        if k >= N_PER: break
cols = cnt["M"] + cnt["X"]
rates = {k: v / cols for k, v in cnt.items() if k != "M"}
rates["n_reads"] = len(CHROMS) * N_PER; rates["aligned_columns"] = cols
json.dump(rates, open("errprof.json", "w"), indent=1)
for k in sorted(rates): print(k, rates[k])
