#!/usr/bin/env python3
"""Non-exon-signal rescue POC (per family). For one Soto family, align its reads to the
FULL-GENOMIC sequences of ALL its copies (UTR/intron/flank included) and measure how many reads
become uniquely assignable (MAPQ>0) — i.e. how many exon-homogenized (K=0) members are rescuable
by non-exon signal the exon-space PSV search never sees.

Usage: nonexon_rescue.py <FAMILY_ID>   e.g. nonexon_rescue.py ID_431
Emits one JSON object on stdout. Deterministic (fixed FLANK, fixed minimap2 flags)."""
import sys, json, subprocess, tempfile, os

BED   = "/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/80_fams.chr.bed"
FLOOR = "/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/soto_floor_decomposition.tsv"
D     = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
FASTA = f"{D}/chm13v2.0.fa"
BAM   = f"{D}/A119b.t2t.bam"
MM    = "/home/juanfra/miniforge3/bin/minimap2"
ST    = "samtools"
FLANK = 1500
MIN_READS = 3          # existing seed floor: a member is "rescued" iff >=3 distinguishable reads

def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=900)

def family_copies(fid):
    """[(gene, chrom, start, end)] for every member of this family."""
    out = []
    for ln in open(BED):
        c = ln.rstrip("\n").split("\t")
        tag = c[3]  # GENE|ID_k
        if tag.endswith("|" + fid) or ("|" + fid) in tag:
            gene = tag.split("|")[0]
            out.append((gene, c[0], int(c[1]), int(c[2])))
    return out

def k0_members(fid):
    """(chrom,start,end,gene) of the K=0 (exon-homogenized/young-identical) members of this family."""
    out = []
    for ln in open(FLOOR):
        c = ln.rstrip("\n").split("\t")
        if c[0] == fid and ("exon-homogenized" in c[-1] or "young identical" in c[-1]):
            out.append((c[2], int(c[3]), int(c[4]), c[1]))
    return out

def run(fid):
    copies = family_copies(fid)
    k0 = k0_members(fid)
    res = {"family_id": fid, "n_copies": len(copies), "n_k0_members": len(k0),
           "genes": sorted({g for (g, _, _, _) in copies})}
    if len(copies) < 2:
        res.update(status="n/a", reason="single-copy family (no multi-copy reference)")
        return res

    with tempfile.TemporaryDirectory(dir="/home/juanfra/winloci_scratch") as td:
        ref = os.path.join(td, "ref.fa")
        # multi-copy full-genomic reference (each copy + flank), named copy0..copyN with its span
        copy_spans = []
        with open(ref, "w") as fo:
            for i, (g, c, s, e) in enumerate(copies):
                lo = max(0, s - FLANK); hi = e + FLANK
                fa = sh(f'{ST} faidx {FASTA} {c}:{lo}-{hi}').stdout
                if not fa or ">" not in fa:
                    continue
                seq = "".join(fa.splitlines()[1:])
                fo.write(f">copy{i}\n{seq}\n")
                copy_spans.append((i, g, c, s, e, lo, hi))
        if len(copy_spans) < 2:
            res.update(status="n/a", reason="could not extract >=2 copy sequences")
            return res

        # family reads: all copy loci, deduped by read name (primary preferred), to fastq
        fq = os.path.join(td, "reads.fq")
        regions = " ".join(f"{c}:{s}-{e}" for (_, _, c, s, e, _, _) in copy_spans)
        sh(f'''{ST} view {BAM} {regions} 2>/dev/null | awk '!seen[$1]++{{print "@"$1"\\n"$10"\\n+\\n"$11}}' > {fq}''')
        n_reads = sum(1 for _ in open(fq)) // 4 if os.path.exists(fq) else 0
        res["n_reads"] = n_reads
        if n_reads < MIN_READS:
            res.update(status="n/a", reason=f"too few family reads ({n_reads})")
            return res

        # align reads to the multi-copy full-genomic reference
        sam = sh(f'{MM} -ax map-hifi --eqx {ref} {fq} 2>/dev/null | {ST} view 2>/dev/null').stdout
        # per copy: count distinguishable primaries (MAPQ>0), tie (MAPQ0), soft-clip flank
        per_copy = {}       # copy_ref -> distinguishable read count
        distinguishable = tie = softclip = 0
        for row in sam.splitlines():
            f = row.split("\t")
            if len(f) < 6:
                continue
            flag = int(f[1])
            if flag & 0x100 or flag & 0x800:   # skip secondary/supplementary
                continue
            if int(f[3]) == 0 and f[2] == "*":  # unmapped
                continue
            mapq = int(f[4]); cig = f[5]
            if mapq > 0:
                distinguishable += 1
                per_copy[f[2]] = per_copy.get(f[2], 0) + 1
                # soft-clip >=20bp?
                import re
                if any(int(m) >= 20 for m in re.findall(r'(\d+)S', cig)):
                    softclip += 1
            else:
                tie += 1
        res["n_distinguishable"] = distinguishable
        res["n_tie"] = tie
        aligned = distinguishable + tie
        res["pct_distinguishable"] = round(100 * distinguishable / aligned, 1) if aligned else 0.0
        res["n_softclip_flank_of_distinguishable"] = softclip

        # which K=0 members are rescued: their copy has >=MIN_READS distinguishable reads.
        # map each K=0 member to its copy index by span overlap.
        rescued = []
        for (kc, ks, ke, kg) in k0:
            for (i, g, c, s, e, lo, hi) in copy_spans:
                if c == kc and min(ke, e) > max(ks, s):   # overlaps this copy
                    if per_copy.get(f"copy{i}", 0) >= MIN_READS:
                        rescued.append({"gene": kg, "locus": f"{kc}:{ks}-{ke}",
                                        "distinguishable_reads": per_copy.get(f"copy{i}", 0)})
                    break
        res["k0_members_rescued"] = rescued
        res["n_k0_rescued"] = len(rescued)
        res["status"] = "ok"
        return res

if __name__ == "__main__":
    print(json.dumps(run(sys.argv[1])))
