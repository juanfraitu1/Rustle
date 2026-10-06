"""align + run: the shipped minimap2 command and the Rustle stages on a scenario directory.

Stages (each optional, each a documented recipe, logs in DIR/logs/):
  assemble  tools/rustle_pipeline.sh assemble   (as_table + copy_assign --assemble-only --genome-wide, shipped polish, f1v2)
  denovo    tools/rustle_pipeline.sh families   -> DIR/run.fam.clusters.tsv + run.fam.copies.tsv/.fa
  guided    truth.gff3 gene regions -> samtools faidx -> minimap2 -x asm20 -c --eqx -P -> mcl_families --paf --gff
            (figures/_o1_recovery.py::guided_families, the documented guided recipe) -> DIR/guided.clusters.tsv
  assign    copy_assign --families run.fam.copies.tsv (or the legacy catalog with --copy-table catalog) -> run.assign.*
  catalog   gw_family_catalog (legacy copy table) -> run.cat.*
  flag      missing_copy_flag --scan-only + --from-scan on a splice .mmi of genome.ref.fa -> run.flag.missing_copy.tsv
Binaries from --bin / RUSTLE_BIN. Nothing is rebuilt here.
"""
import os
import shutil
import subprocess

from sim import MM2  # noqa: E402  the shipped read mapping

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DRIVER = os.path.join(REPO, "tools", "rustle_pipeline.sh")
DEFAULT_BIN = os.environ.get("RUSTLE_BIN", "/mnt/linuxdisk/home/juanfraitu/rustle_target/release")
STAGES = ("assemble", "denovo", "guided", "assign", "flag")


def sh(cmd, log_path, env=None):
    """Run a shell command, stdout+stderr to log_path; raise with the log's tail on failure."""
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "w") as lg:
        rc = subprocess.run(cmd, shell=True, stdout=lg, stderr=subprocess.STDOUT, env=env).returncode
    if rc != 0:
        tail = "".join(open(log_path).readlines()[-15:])
        raise RuntimeError(f"command failed ({rc}): {cmd}\n--- {log_path} tail ---\n{tail}")
    return rc


def align(out_dir, threads=2, log=print):
    ref, fq, bam = (os.path.join(out_dir, x) for x in ("genome.ref.fa", "reads.fq", "reads.bam"))
    sh(f"{MM2} -t {threads} {ref} {fq} | samtools sort -@1 -o {bam} - && samtools index {bam}", os.path.join(out_dir, "logs", "align.log"))
    n = subprocess.run(f"samtools view -c -F 2308 {bam}", shell=True, stdout=subprocess.PIPE, text=True).stdout.strip()
    log(f"align: {n} primary alignments -> reads.bam  ({MM2})")
    return bam


def run(out_dir, stages=STAGES, bin_dir=DEFAULT_BIN, threads=2, copy_table="families", log=print, env_extra=None):
    ref, bam = os.path.join(out_dir, "genome.ref.fa"), os.path.join(out_dir, "reads.bam")
    prefix = os.path.join(out_dir, "run")
    logs = os.path.join(out_dir, "logs")
    env = dict(os.environ)
    env.setdefault("TMPDIR", os.environ.get("TMPDIR", "/mnt/linuxdisk/tmp") if os.path.isdir("/mnt/linuxdisk/tmp") else "/tmp")
    env.update(env_extra or {})
    for b in ("copy_assign", "mcl_families", "as_table"):
        if not os.path.exists(os.path.join(bin_dir, b)):
            raise FileNotFoundError(f"{b} not found in {bin_dir} (--bin / RUSTLE_BIN)")
    drv = f"bash {DRIVER} {{stage}} --bam {bam} --fasta {ref} --out {prefix} --bin {bin_dir} --threads {threads} --no-cache"
    done = {}
    if "assemble" in stages:
        sh(drv.format(stage="assemble"), os.path.join(logs, "assemble.log"), env)
        n = sum(1 for l in open(prefix + ".gtf") if "\ttranscript\t" in l)
        log(f"assemble: {n} transcripts -> run.gtf"); done["assemble"] = prefix + ".gtf"
    if "denovo" in stages:
        sh(drv.format(stage="families"), os.path.join(logs, "families.log"), env)
        n = len({l.split('\t')[0] for l in open(prefix + ".fam.clusters.tsv") if not l.startswith("cluster_id")})
        log(f"denovo: {n} clusters -> run.fam.clusters.tsv"); done["denovo"] = prefix + ".fam.clusters.tsv"
    if "guided" in stages:
        done["guided"] = guided(out_dir, bin_dir, threads, log)
    if "catalog" in stages or (copy_table == "catalog" and "assign" in stages):
        sh(drv.format(stage="catalog"), os.path.join(logs, "catalog.log"), env)
        log("catalog: run.cat.copies.tsv"); done["catalog"] = prefix + ".cat.copies.tsv"
    if "assign" in stages:
        tab = prefix + (".cat" if copy_table == "catalog" else ".fam")
        if not os.path.exists(tab + ".copies.tsv"):
            raise FileNotFoundError(f"{tab}.copies.tsv missing: run the {'catalog' if copy_table == 'catalog' else 'denovo'} stage first")
        n_cop = sum(1 for l in open(tab + ".copies.tsv")) - 1
        if n_cop == 0:
            log("assign: no copies in the copy table (no multi-copy family found) -> nothing to assign"); done["assign"] = None
        else:
            regions = prefix + ".regions.txt"
            sh(f"samtools view -H {bam} | awk '$1==\"@SQ\"{{sub(\"SN:\",\"\",$2); sub(\"LN:\",\"\",$3); print $2\":1-\"$3}}' > {regions}",
               os.path.join(logs, "regions.log"))
            sh(f"{bin_dir}/copy_assign --bam {bam} --fasta {ref} --regions {regions} --families {tab}.copies.tsv --copies-fa {tab}.copies.fa "
               f"--threads {threads} --out {prefix}.assign", os.path.join(logs, "assign.log"), env)
            n = sum(1 for l in open(prefix + ".assign.assignments.tsv")) - 1
            log(f"assign: {n} read x family rows -> run.assign.assignments.tsv"); done["assign"] = prefix + ".assign.assignments.tsv"
    if "flag" in stages:
        done["flag"] = flag(out_dir, bin_dir, threads, log, env)
    return done


def guided(out_dir, bin_dir, threads, log=print):
    """The guided O1 recipe on the truth annotation's gene bodies."""
    ref, gff = os.path.join(out_dir, "genome.ref.fa"), os.path.join(out_dir, "truth.gff3")
    prefix = os.path.join(out_dir, "guided")
    logs = os.path.join(out_dir, "logs")
    ref_contigs = {l.split("\t")[0] for l in open(ref + ".fai")}
    regs = sorted({f"{f[0]}:{f[3]}-{f[4]}" for f in (l.rstrip("\n").split("\t") for l in open(gff) if not l.startswith("#"))
                   if len(f) > 4 and f[2] in ("gene", "pseudogene") and f[0] in ref_contigs}, key=lambda s: s.encode())
    with open(prefix + ".regions", "w") as fh:
        fh.write("\n".join(regs) + "\n")
    sh(f"samtools faidx {ref} -r {prefix}.regions > {prefix}.bodies.fa", os.path.join(logs, "guided.faidx.log"))
    sh(f"minimap2 -x asm20 -c --eqx -P -t {threads} {prefix}.bodies.fa {prefix}.bodies.fa > {prefix}.paf", os.path.join(logs, "guided.mm2.log"))
    # mcl_families --gff reads gene/pseudogene records by Name= and exons by gene=; absent-contig genes are not nodes
    sh(f"{bin_dir}/mcl_families --paf {prefix}.paf --gff {gff} --min-exonic-bp 1 --min-shared-exon-frac 0.60 --out {prefix}",
       os.path.join(logs, "guided.mcl.log"))
    n = len({l.split('\t')[0] for l in open(prefix + ".clusters.tsv") if not l.startswith("cluster_id")})
    log(f"guided: {len(regs)} nodes, {n} clusters -> guided.clusters.tsv")
    return prefix + ".clusters.tsv"


def flag(out_dir, bin_dir, threads, log=print, env=None):
    """missing_copy_flag on the de novo loci (run.families.gtf, else run.gtf) with a splice index of the reference."""
    ref, bam = os.path.join(out_dir, "genome.ref.fa"), os.path.join(out_dir, "reads.bam")
    prefix = os.path.join(out_dir, "run")
    logs = os.path.join(out_dir, "logs")
    mmi = os.path.join(out_dir, "genome.ref.splice.mmi")
    if not os.path.exists(mmi):
        sh(f"minimap2 -x splice:hq -d {mmi} {ref}", os.path.join(logs, "index.log"))
    loci = prefix + ".families.gtf" if os.path.exists(prefix + ".families.gtf") else prefix + ".gtf"
    if not os.path.exists(loci):
        raise FileNotFoundError("flag needs the assembled loci: run the assemble stage first")
    confirm = ""
    tmmi = os.path.join(out_dir, "genome.truth.splice.mmi")
    if not os.path.exists(tmmi):
        sh(f"minimap2 -x splice:hq -d {tmmi} {os.path.join(out_dir, 'genome.truth.fa')}", os.path.join(logs, "index_truth.log"))
    confirm = f"--confirm truth={tmmi}"
    sh(f"{bin_dir}/missing_copy_flag --bam {bam} --fasta {ref} --loci {loci} --index x --threads {threads} --out {prefix}.flag_scan --scan-only",
       os.path.join(logs, "flag_scan.log"), env)
    sh(f"{bin_dir}/missing_copy_flag --bam {bam} --fasta {ref} --loci {loci} --index {mmi} --threads {threads} {confirm} "
       f"--out {prefix}.flag --from-scan {prefix}.flag_scan", os.path.join(logs, "flag.log"), env)
    out = prefix + ".flag.missing_copy.tsv"
    n = sum(1 for l in open(out)) - 1 if os.path.exists(out) else 0
    log(f"flag: {n} loci -> run.flag.missing_copy.tsv")
    return out
