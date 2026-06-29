#!/usr/bin/env python
"""F1 — O4 positive control on REAL data: mask known copies out of the reference, show O4 re-admits them.

The bounded O4 run on the T2T-complete gorilla admitted 0 reference-absent copies. That is ambiguous: either
the detector works and the reference is simply complete, OR the detector is inert. This resolves it: we DELETE
(N-mask) several copies of a real co-located multi-copy family from the reference, REALIGN that locus's reads
to the masked reference (the masked copies' reads now collapse onto the surviving paralogs), and run
`--absent-copies`. If O4 RE-ADMITS the masked copies as `AC_` copies, the detector demonstrably fires when a
copy is truly absent — so "0 admitted on real GGO" means the reference is complete, not that O4 is inert.

Family GWFAM71 = a 9-copy tandem on NC_073247.2 (~59.70-59.86 Mb). We mask 3 well-supported copies (idx 2,3,6)
and keep the rest as the reference host. O4's gate-1 needs >=3 haplotypes at the locus (host + >=2 divergent
copies), which masking 3 copies provides.

CONTROL: the identical pipeline against the UNMASKED contig must admit 0 AC_ (the copies are present -> remap
>=98% -> DnaNeeds), isolating "absence" as the cause of admission.

Run: /home/juanfra/miniforge3/bin/python bench/o4_mask_readmit.py
"""
import os
import subprocess
import sys

import pysam

SCRATCH = "/home/juanfra/winloci_scratch"
GENOME = f"{SCRATCH}/GGO.fasta"
BAM = f"{SCRATCH}/GGO_mm.bam"
BIN = "/mnt/c/Users/jfris/Desktop/Rustle/target/release/copy_assign"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
SAM = "/home/juanfra/miniforge3/bin/samtools"
WORK = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/03f0156d-63b6-4d0f-bf09-862bddee491c/scratchpad/o4_mask"

# GWFAM47 = a 6-copy co-located family on NC_073239.2 whose copies are ~95% identical (min 0.948, max 0.959) —
# the ADMISSION regime (divergent enough for remap<98%, similar enough to collapse, >=3-way). We mask 3
# well-supported copies (idx 0,1,5 = 47,63,61 reads) and keep copy2 (220 reads) as the reference host; the
# masked copies' reads collapse onto copy2 as ~95% divergent, mutually-distinct clusters -> >=3 haplotypes.
CONTIG = "NC_073239.2"
REGION = f"{CONTIG}:122615409-123378628"
MASK = [(122615409, 122617135), (122634080, 122635806), (123333518, 123378628)]


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, text=True, capture_output=True, **kw)


def build_refs():
    # clear any stale artifacts from a previous family's run (esp. contig.fa.fai pointing at the wrong contig)
    if os.path.isdir(WORK):
        import shutil
        shutil.rmtree(WORK)
    os.makedirs(WORK, exist_ok=True)
    # extract the whole contig
    ctg = f"{WORK}/contig.fa"
    r = sh(f"{SAM} faidx {GENOME} {CONTIG} > {ctg} && {SAM} faidx {ctg}")
    if r.returncode:
        sys.exit(f"faidx failed: {r.stderr}")
    fa = pysam.FastaFile(ctg)
    seq = list(fa.fetch(CONTIG))
    # masked copy: N out the masked intervals
    masked = seq[:]
    nmask = 0
    for (s, e) in MASK:
        for i in range(s, e):
            masked[i] = "N"
            nmask += 1
    unmasked_fa = f"{WORK}/unmasked.fa"
    masked_fa = f"{WORK}/masked.fa"
    with open(unmasked_fa, "w") as fh:
        fh.write(f">{CONTIG}\n{''.join(seq)}\n")
    with open(masked_fa, "w") as fh:
        fh.write(f">{CONTIG}\n{''.join(masked)}\n")
    sh(f"{SAM} faidx {unmasked_fa}")
    sh(f"{SAM} faidx {masked_fa}")
    print(f"masked {nmask:,} bp across {len(MASK)} copies on {CONTIG}")
    return unmasked_fa, masked_fa


def extract_reads():
    reads_fa = f"{WORK}/reads.fa"
    # primary, non-supplementary region reads -> fasta (their original sequences)
    r = sh(f"{SAM} view -F 0x900 {BAM} {REGION} | "
           f"awk '{{print \">\"$1\"\\n\"$10}}' > {reads_fa}")
    n = sum(1 for _ in open(reads_fa)) // 2
    print(f"extracted {n} region reads")
    return reads_fa


def realign(ref_fa, reads_fa, tag):
    bam = f"{WORK}/{tag}.bam"
    cmd = (f"{MM2} -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 4 {ref_fa} {reads_fa} 2>/dev/null "
           f"| {SAM} sort -o {bam} - && {SAM} index {bam}")
    r = sh(cmd)
    if r.returncode:
        sys.exit(f"realign {tag} failed: {r.stderr}")
    return bam


def run_o4(ref_fa, bam, tag):
    out = f"{WORK}/{tag}_o4"
    env = dict(os.environ)
    env.pop("RUSTLE_ABSENT_MMI", None)          # remap target = ref_fa (the MASKED/unmasked contig)
    env["RUSTLE_MINIMAP2"] = MM2
    r = subprocess.run(
        [BIN, "--bam", bam, "--fasta", ref_fa, "--region", REGION,
         "--min-copies", "2", "--absent-copies", "--skip-poa-diagnostic", "--out", out],
        text=True, capture_output=True, env=env,
    )
    if r.returncode:
        sys.exit(f"copy_assign {tag} failed:\n{r.stderr[-1500:]}")
    assigns = f"{out}.assignments.tsv"
    dna = f"{out}.dna_needs.tsv"
    ac = 0
    if os.path.exists(assigns):
        ac = sum(1 for line in open(assigns) if "AC_" in line)
    dn = (sum(1 for _ in open(dna)) - 1) if os.path.exists(dna) else 0
    return ac, dn, out


def main():
    unmasked_fa, masked_fa = build_refs()
    reads_fa = extract_reads()
    print("\n--- realigning + running O4 on the MASKED reference (copies deleted) ---")
    mbam = realign(masked_fa, reads_fa, "masked")
    m_ac, m_dn, m_out = run_o4(masked_fa, mbam, "masked")
    print("--- CONTROL: identical pipeline on the UNMASKED reference (copies present) ---")
    cbam = realign(unmasked_fa, reads_fa, "unmasked")
    c_ac, c_dn, c_out = run_o4(unmasked_fa, cbam, "unmasked")

    print("\n================= RESULT =================")
    print(f"MASKED  reference: AC_ copies admitted = {m_ac:>3} | DNA-needs = {m_dn}")
    print(f"CONTROL (unmasked): AC_ copies admitted = {c_ac:>3} | DNA-needs = {c_dn}")
    verdict = "PASS — O4 re-admits the absent copy(ies); inert on the complete control" \
        if (m_ac > 0 and c_ac == 0) else \
        ("PARTIAL — masked admits but control also admits (specificity issue)" if m_ac > 0 else
         "NO ADMISSION — masked copies routed to DNA-needs; see reasons below (gate is conservative on real data)")
    print(f"VERDICT: {verdict}")
    dnf = f"{m_out}.dna_needs.tsv"
    if os.path.exists(dnf):
        print(f"\nMASKED-run DNA-needs reasons (top):")
        import collections
        reasons = collections.Counter(line.rstrip("\n").split("\t")[-1] for line in open(dnf).readlines()[1:])
        for reason, k in reasons.most_common():
            print(f"  {k:>3}  {reason}")
    print(f"\noutputs under {WORK}/")


if __name__ == "__main__":
    main()
