"""Shared helpers for the Phase 0 flow-recall diagnostic."""
import json, os, re, subprocess  # subprocess used by ensure_locus_logs (Task 2)

# Repo root derived from this file's location (bench/flow_recall_phase0/lib.py).
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RG = f"{ROOT}/bench/copy_recovery_eval/results_genomewide"
PERCHROM = f"{RG}/perchrom"
CACHE = f"{ROOT}/bench/flow_recall_phase0/cache"
RUSTLE = f"{ROOT}/target/release/rustle"
ST = f"{ROOT}/tools/stringtie/stringtie"

def parse_gtf(path):
    """transcript_id -> {strand, chrom, exons, introns(tuple), span, splice_sites(set)}"""
    tx = {}
    if not os.path.exists(path):
        return tx
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m:
                continue
            d = tx.setdefault(m.group(1), {"strand": f[6], "chrom": f[0], "exons": []})
            d["exons"].append((int(f[3]), int(f[4])))
    for d in tx.values():
        d["exons"].sort()
        ex = d["exons"]
        d["introns"] = tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
        d["span"] = (ex[0][0], ex[-1][1]) if ex else (0, 0)
        d["splice_sites"] = {c for iv in d["introns"] for c in iv}
    return tx

def chain_str(introns):
    return ",".join(f"{a}-{b}" for a, b in introns)

def parse_chain_str(s):
    return tuple(tuple(int(x) for x in p.split("-")) for p in s.split(",") if "-" in p)

def read_parity(path, step):
    out = []
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            try:
                o = json.loads(line)
            except ValueError:
                continue
            if o.get("step") == step:
                out.append(o)
    return out

def ref_chain(chrom, tid):
    """Return (span_lo, span_hi, introns_tuple) for an annotated transcript, or None."""
    ex = []
    path = f"{PERCHROM}/{chrom}/ref.gtf"
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[2] != "exon" or tid not in f[8]:
                continue
            ex.append((int(f[3]), int(f[4])))
    if not ex:
        return None
    ex.sort()
    introns = tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
    return ex[0][0], ex[-1][1], introns

def _run_tool(binpath, log_var, slice_bam, chrom, out_gtf, log_file):
    env = dict(os.environ)
    env[log_var] = log_file
    env[log_var.replace("LOG", "FILTER_CHROM")] = chrom
    env[log_var.replace("LOG", "FILTER_STEPS")] = "path_extracted,junction_accept"
    subprocess.run([binpath, "-L", slice_bam, "-o", out_gtf], env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def ensure_locus_logs(chrom, tid, pad=8000):
    """Idempotently slice the locus and run BOTH tools with parity logging. Cached."""
    wd = f"{CACHE}/{chrom}/{tid}"
    if os.path.exists(f"{wd}/.done"):
        return wd
    ch = ref_chain(chrom, tid)
    if ch is None:
        return None
    lo, hi, _ = ch
    os.makedirs(wd, exist_ok=True)
    sb = f"{wd}/slice.bam"
    reg = f"{chrom}:{max(1, lo - pad)}-{hi + pad}"
    with open(sb, "wb") as fh:
        subprocess.run(["samtools", "view", "-b", f"{PERCHROM}/{chrom}/c.bam", reg],
                       stdout=fh, stderr=subprocess.DEVNULL)
    subprocess.run(["samtools", "index", sb], stderr=subprocess.DEVNULL)
    _run_tool(RUSTLE, "RUSTLE_PARITY_LOG", sb, chrom, f"{wd}/r.gtf", f"{wd}/r.jsonl")
    _run_tool(ST, "STRINGTIE_PARITY_LOG", sb, chrom, f"{wd}/s.gtf", f"{wd}/s.jsonl")
    # No check=True: a tool failure should not abort the 668-locus resumable sweep.
    # Analyses tolerate empty logs; a missing-.done locus is simply re-run next pass.
    open(f"{wd}/.done", "w").close()
    return wd
