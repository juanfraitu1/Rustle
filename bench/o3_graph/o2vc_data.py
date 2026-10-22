#!/usr/bin/env python3
"""The caller side of O2 from the sample's own variant calls (design 2026-10-22): run longcallR on the teacher reads of a copy
and read its VCF into the Calls that o2vc_rule.genotype() turns into a copy's alleles at a position. Nothing here decides a
placement (o2vc_rule.py does) and nothing reads truth.

What longcallR 2.0.0 does, observed on SMN2 (HG002, chr5:70808816-70837687, primary records, `-p hifi-masseq --min-mapq 0
--no-bam-output`), and what the code below rests on:
- `-o PREFIX` writes PREFIX.vcf and nothing else: VCF 4.3, plain text, one sample ("Sample"). It does not make the directory of
  PREFIX (exit 101). A region without reads still exits 0 and leaves a header-only VCF, so a missing file is an error and an
  empty one is not.
- FILTER is one of PASS, LowQual, HomRef, RnaEdit, Multiallelic, Dense (the header's six; SMN2 shows the first four). FORMAT is
  GT:GQ:PS:DP:AF:PQ on a PASS or LowQual record (PS "." when unphased) and GT:GQ:DP:AF on the others, so FORMAT keys are read by
  name, never by position. REF is upper-cased though the CHM13 FASTA is soft-masked (46 of 47 SMN2 records sit on a lowercase
  base): compare it with the genome base upper-cased.
- `-r CHROM:START-END`: START is 1-based and must be >= 1 (0 panics: exit 101, no VCF), END is exclusive, so a record is in the
  region when START <= POS < END.
- It needs the BAM sorted and indexed and a `.fai` beside the FASTA (it does not make one: exit 101), and it reads the whole FASTA:
  16 s and 3.4 GB RSS against the 3.1 GB CHM13 file, 0.7 s and 0.5 GB against a one-contig FASTA of chr5, with the same calls.
- It is not deterministic: on identical inputs the phase labels of phased records flip (0|1 / 1|0) and in about 1 run in 10 two weak
  T>C records move from HomRef / LowQual to RnaEdit (SMN2, 41 runs at -t 1 and -t 2). The rule reads a GT as a set of alleles, so
  the labels never matter, but a re-call can differ at a weak RnaEdit-like site: keep the VCF a run used, do not re-call to
  reproduce it.
"""
import gzip
import importlib.util
import os
import re
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
LONGCALLR = "/home/juanfra/miniforge3/envs/longcallr/bin/longcallR"
PRESET = "hifi-masseq"                                          # transcript-strand HiFi reads: no strand-bias filter
MIN_MAPQ = 0                                                    # the teachers are chosen by AS margin, not by MAPQ


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_vc", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


Call = _sibling("o2vc_rule").Call


def longcallr_version(binary=LONGCALLR):
    """The caller's own version line ("longcallR 2.0.0"), kept in run.json so a result names the caller that made it."""
    return subprocess.run([binary, "--version"], capture_output=True, text=True, check=True).stdout.strip()


def run_longcallr(teacher_bam, genome, out_prefix, binary=LONGCALLR, threads=2, region=None):
    """Run longcallR (PRESET, MIN_MAPQ, no BAM output, every other option its default) on teacher_bam against genome and return
    the VCF it wrote, out_prefix + ".vcf". teacher_bam must be sorted and indexed, genome a FASTA with its .fai. region is
    "CHROM:START-END" with START 1-based and >= 1 and END exclusive (the records with START <= POS < END), or None for the whole
    BAM. The directory of out_prefix is made first (longcallR does not). A non-zero exit raises RuntimeError with the last 20
    lines of stderr, as does an exit 0 that wrote no VCF."""
    cmd = [binary, "-b", teacher_bam, "-f", genome, "-p", PRESET, "--min-mapq", str(MIN_MAPQ), "--no-bam-output",
           "-t", str(threads), "-o", out_prefix] + (["-r", region] if region else [])
    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    tail = "\n".join(r.stderr.splitlines()[-20:])
    if r.returncode:
        raise RuntimeError(f"longcallR exited {r.returncode} ({' '.join(cmd)}):\n{tail}")
    vcf = out_prefix + ".vcf"
    if not os.path.exists(vcf):
        raise RuntimeError(f"longcallR wrote no VCF at {vcf} ({' '.join(cmd)}):\n{tail}")
    return vcf


def parse_vcf(path):
    """{(chrom, POS - 1): Call} for every record of a longcallR VCF, .vcf or .vcf.gz (POS is 1-based in the file, 0-based in the
    key). FORMAT keys are read by name. gt is the tuple of the GT's allele indices (0 = REF), () when it calls none ("./.",
    ".|.", "."): a GT missing only one allele counts as calling none too, for half a genotype would make a copy look homozygous.
    The record is kept either way, and o2vc_rule.genotype() leaves such a copy unknown. dp is the record's own DP (0 when absent):
    the caller's count, not the teacher depth the rule's MIN_DEPTH uses."""
    out = {}
    with (gzip.open if path.endswith(".gz") else open)(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\r\n").split("\t")
            f = dict(zip(c[8].split(":"), c[9].split(":")))
            idx = re.split(r"[/|]", f["GT"])
            gt = tuple(int(i) for i in idx) if all(i.isdigit() for i in idx) else ()
            alts = () if c[4] == "." else tuple(c[4].split(","))
            out[c[0], int(c[1]) - 1] = Call(c[3], alts, gt, c[6], int(f["DP"]) if f.get("DP", ".").isdigit() else 0)
    return out
