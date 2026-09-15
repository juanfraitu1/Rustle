#!/usr/bin/env python3
"""Evidence-agnostic families (spec 2026-09-14 §1a): duplication evidence D1-D4 in one internal format.

Pair table (TSV, header = PAIR_COLS): 0-based half-open intervals of side A and side B, strands, identity, CIGAR (M/I/D;
M both sides, D side A, I side B; side B reverse-complemented when strand_b is '-'), source (D1:sedef, D2:biser, D3:selfaln).
"""
import argparse
import csv
import glob
import os
import re
import shutil
import signal
import subprocess
import sys
import time

PAIR_COLS = ("chrom_a", "start_a", "end_a", "chrom_b", "start_b", "end_b", "strand_a", "strand_b", "identity", "cigar", "source")
MIN_SD_BP, MIN_SD_ID = 1000, 0.90


def normalize_cigar(cg):
    ops = [(int(n), "M" if o in "=X" else o) for n, o in re.findall(r"(\d+)([MIDNSHP=X])", cg)]
    for n, o in ops:
        if o not in ("M", "I", "D"):
            raise ValueError(f"Invalid CIGAR operation: {o}")
    out = []
    for n, o in ops:
        if out and out[-1][1] == o:
            out[-1] = (out[-1][0] + n, o)
        else:
            out.append((n, o))
    return "".join(f"{n}{o}" for n, o in out)


def from_sedef(line, fmt):
    f = line.rstrip("\n").split("\t")
    if len(f) < 23 or not f[1].isdigit():
        return None
    if fmt == "gorilla":
        ident, cig = float(f[20]), (f[32] if len(f) > 32 else "")
    elif fmt == "human":
        ident, cig = float(f[22]), ""
    else:
        raise ValueError(fmt)
    return (f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]), f[8], f[9], ident, cig, "D1:sedef")


def from_biser(line):
    f = line.rstrip("\n").split("\t")
    if len(f) < 13 or not f[1].isdigit():
        return None
    m = re.search(r"X=([\d.]+)", f[13]) if len(f) > 13 else None
    err = float(m.group(1)) if m else float(f[7])
    return (f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]), f[8], f[9], round(1 - err / 100, 9),
            normalize_cigar(f[12]), "D2:biser")


def from_selfpaf(line):
    """minimap2 self-alignment record of a chunk (query `chrom@offset`) against the genome. Side A = target (forward),
    side B = query with the record's strand, so minimap2's CIGAR (D = target, I = query) already has SEDEF semantics.
    Keeps canonical non-self records >= MIN_SD_BP at identity >= MIN_SD_ID."""
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        return None
    if "@" not in f[0]:
        return None
    qc, off = f[0].rsplit("@", 1)
    qs, qe = int(off) + int(f[2]), int(off) + int(f[3])
    tc, ts, te = f[5], int(f[7]), int(f[8])
    if qc == tc and qs < te and ts < qe:
        return None
    nm, bl = int(f[9]), int(f[10])
    if bl < MIN_SD_BP or nm / bl < MIN_SD_ID:
        return None
    if (tc, ts) > (qc, qs):
        return None
    cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
    return (tc, ts, te, qc, qs, qe, "+", f[4], round(nm / bl, 9), normalize_cigar(cg), "D3:selfaln")


def to_sedef_gorilla(p):
    ca, a1, a2, cb, b1, b2, sa, sb, ident, cig, _src = p
    L = sum(int(n) for n, o in re.findall(r"(\d+)([MID])", cig) if o == "M") or max(a2 - a1, b2 - b1)
    m = round(ident * L)
    mm = L - m
    frac = m / (m + mm) if L else 0.0
    row = [""] * 34
    row[0:6] = [ca, str(a1), str(a2), cb, str(b1), str(b2)]
    row[6], row[7], row[8], row[9] = "S", f"{(1 - ident) * 100:.1f}", sa, sb
    row[16], row[17], row[20], row[32], row[33] = str(m), str(mm), f"{frac:.6f}", cig, f"{frac:.6f}"
    return "\t".join(row)


def write_pairs(pairs, path):
    with open(path, "w") as fh:
        fh.write("\t".join(PAIR_COLS) + "\n")
        for p in pairs:
            fh.write("\t".join(map(str, p)) + "\n")


def read_pairs(path):
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        out.append((r["chrom_a"], int(r["start_a"]), int(r["end_a"]), r["chrom_b"], int(r["start_b"]), int(r["end_b"]),
                    r["strand_a"], r["strand_b"], float(r["identity"]), r["cigar"], r["source"]))
    return out


MERYL = "/home/juanfra/miniforge3/envs/phasing_eval/bin/meryl"
MERYL_LOOKUP = "/home/juanfra/miniforge3/envs/phasing_eval/bin/meryl-lookup"


def valley(hist, max_count=100000):
    cs = sorted(c for c in hist if 2 <= c <= max_count)
    for i in range(1, len(cs) - 1):
        c = cs[i]
        if hist[c] < hist[cs[i - 1]] and hist[c] < hist[cs[i + 1]] and any(hist[d] > hist[c] for d in cs[i + 1:]):
            return c
    return None


def _run(cmd, out=None, log=None):
    """Run a command, capturing stdout and stderr atomically.

    If out is given, writes stdout to out (using a .tmp file and os.replace on success).
    If log is given, appends command line and stderr to log.
    Raises RuntimeError on non-zero exit.
    """
    out_tmp = f"{out}.tmp" if out else None

    # Remove stale .tmp files before starting
    if out_tmp and os.path.exists(out_tmp):
        if os.path.isdir(out_tmp):
            shutil.rmtree(out_tmp)
        else:
            os.remove(out_tmp)

    # Write command header to log
    if log:
        with open(log, "a") as fh:
            fh.write(f"$ {' '.join(cmd)}\n")

    # Run command, capturing stderr
    with (open(out_tmp, "w") if out_tmp else open(os.devnull, "w")) as out_fh:
        result = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True)

    # Append stderr to log
    if log and result.stderr:
        with open(log, "a") as fh:
            fh.write(result.stderr)

    # Check for errors
    if result.returncode != 0:
        # Clean up tmp file on failure
        if out_tmp and os.path.exists(out_tmp):
            if os.path.isdir(out_tmp):
                shutil.rmtree(out_tmp)
            else:
                os.remove(out_tmp)
        if log:
            raise RuntimeError(f"{cmd[0]} failed (exit {result.returncode}); see {log}")
        else:
            raise RuntimeError(f"{cmd[0]} failed (exit {result.returncode})")

    # Atomically rename temp file to final location on success
    if out_tmp and os.path.exists(out_tmp):
        os.replace(out_tmp, out)


def _run_to(cmd_template, final, log=None):
    """Run a command with a tool-named output, atomically producing the final output.

    cmd_template: list with "TMPPATH" placeholder replaced by final.tmp
    final: final output path (directory or file)
    log: log file path
    """
    final_tmp = f"{final}.tmp"

    # Remove stale .tmp files before starting
    if os.path.exists(final_tmp):
        if os.path.isdir(final_tmp):
            shutil.rmtree(final_tmp)
        else:
            os.remove(final_tmp)

    # Replace placeholder with tmp path
    cmd = [final_tmp if x == "TMPPATH" else x for x in cmd_template]

    # Write command header to log
    if log:
        with open(log, "a") as fh:
            fh.write(f"$ {' '.join(cmd)}\n")

    # Run command
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Append stderr to log
    if log and result.stderr:
        with open(log, "a") as fh:
            fh.write(result.stderr)

    # Check for errors
    if result.returncode != 0:
        # Clean up tmp file on failure
        if os.path.exists(final_tmp):
            if os.path.isdir(final_tmp):
                shutil.rmtree(final_tmp)
            else:
                os.remove(final_tmp)
        if log:
            raise RuntimeError(f"{cmd[0]} failed (exit {result.returncode}); see {log}")
        else:
            raise RuntimeError(f"{cmd[0]} failed (exit {result.returncode})")

    # Atomically rename temp file/dir to final location on success
    if os.path.exists(final_tmp):
        os.replace(final_tmp, final)


def cmd_meryl(genome, outdir, threads=4):
    import repeat_evidence as rep
    os.makedirs(outdir, exist_ok=True)
    log = f"{outdir}/meryl.log"

    db = f"{outdir}/kmers.meryl"
    if not os.path.exists(db):
        _run_to([MERYL, "count", "k=31", f"threads={threads}", "memory=12", str(genome), "output", "TMPPATH"], db, log=log)

    hist_path = f"{outdir}/hist.tsv"
    if not os.path.exists(hist_path):
        _run([MERYL, "histogram", db], out=hist_path, log=log)

    hist = {int(a): int(b) for a, b in (l.split()[:2] for l in open(hist_path) if l.strip() and l.split()[0].isdigit())}
    c = valley(hist)
    open(f"{outdir}/cmax.txt", "w").write(f"{c if c is not None else 'NA'}\n")
    if c is None:
        print(f"[meryl] no histogram valley: D4 and R4 not available for {genome}")
        return None

    for name, ops in (("low", ["less-than", str(c + 1), "[", "greater-than", "1", db, "]"]), ("high", ["greater-than", str(c), db])):
        sub = f"{outdir}/{name}.meryl"
        if not os.path.exists(sub):
            _run_to([MERYL] + ops + ["output", "TMPPATH"], sub, log=log)
        bed = f"{outdir}/{name}_copy.runs.bed"
        if not os.path.exists(bed):
            _run_to([MERYL_LOOKUP, "-bed-runs", "-sequence", str(genome), "-mers", sub, "-output", "TMPPATH"], bed, log=log)

    for name, src in (("low", "D4:meryl"), ("high", "R4:meryl")):
        ivs = [(f[0], int(f[1]), int(f[2]), ".") for f in (l.split("\t") for l in open(f"{outdir}/{name}_copy.runs.bed")) if len(f) >= 3]
        rep.write_bed(ivs, f"{outdir}/{name}_copy.bed", src)

    print(f"[meryl] C_max = {c}")
    return c


BISER = "/home/juanfra/miniforge3/envs/biser/bin/biser"


def run_budget(cmd, budget, stdout=None):
    """Run cmd in its own process group; on budget expiry kill the WHOLE group (BISER workers, Liftoff's minimap2) so no
    orphan survives (WSL crash rule). Returns True if it finished."""
    proc = subprocess.Popen(cmd, stdout=stdout or subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
    try:
        rc = proc.wait(timeout=budget)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()
        return False
    if rc != 0:
        raise subprocess.CalledProcessError(rc, cmd)
    return True


def nonrepeat_aligned(p, merged):
    import repeat_evidence as rep
    ca, a1, a2, cb, b1, b2, sa, sb, ident, cig, _ = p
    oa, ob, tot = 0, 0, 0
    for n, o in re.findall(r"(\d+)([MID])", cig):
        n = int(n)
        if o == "M":
            ga = (a1 + oa, a1 + oa + n)
            gb = (b2 - ob - n, b2 - ob) if sb == "-" else (b1 + ob, b1 + ob + n)
            tot += min(n - rep.masked_bases(merged, ca, *ga), n - rep.masked_bases(merged, cb, *gb))
            oa += n
            ob += n
        elif o == "D":
            oa += n
        else:
            ob += n
    return tot


def cmd_biser(a):
    import pysam
    os.makedirs(a.outdir, exist_ok=True)
    bed, tmp = f"{a.outdir}/biser.bed", f"{a.outdir}/biser_tmp"
    if not os.path.exists(bed):
        cmd = [BISER, "-t", str(a.threads), "-o", bed, "--keep-temp", "-T", tmp]
        if os.path.isdir(tmp):
            cmd += ["--resume", tmp]
        if not os.path.exists(str(a.genome) + ".fai"):
            pysam.faidx(str(a.genome))
        if not run_budget(cmd + [str(a.genome)], a.budget):
            print("[biser] budget spent; rerun the same command to resume")
            return
    pairs = [p for p in (from_biser(l) for l in open(bed)) if p]
    write_pairs(pairs, f"{a.outdir}/pairs.D2.tsv")
    print(f"[biser] {len(pairs)} pairs -> {a.outdir}/pairs.D2.tsv")


def cmd_selfaln(a):
    import pysam
    import repeat_evidence as rep
    d = f"{a.outdir}/selfaln"
    os.makedirs(d, exist_ok=True)
    g = pysam.FastaFile(str(a.genome))
    mmi = f"{d}/genome.mmi"
    if not os.path.exists(mmi):
        subprocess.run(["minimap2", "-x", "asm20", "-t", str(a.threads), "-d", mmi, str(a.genome)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    chunks = []
    for c, L in zip(g.references, g.lengths):
        for off in range(0, L, a.chunk_bp):
            chunks.append((c, off, min(L, off + a.chunk_bp)))
    t0 = time.time()
    done = 0
    for i, (c, s, e) in enumerate(chunks):
        paf = f"{d}/c{i:05d}.paf"
        if os.path.exists(paf):
            done += 1
            continue
        if time.time() - t0 > a.budget:
            break
        q = f"{d}/c{i:05d}.fa"
        open(q, "w").write(f">{c}@{s}\n{g.fetch(c, s, e).upper()}\n")
        with open(paf + ".tmp", "w") as fh:
            ok = run_budget(["minimap2", "-x", "asm20", "-c", "-N", "50", "-p", "0.1", "-t", str(a.threads), mmi, q],
                            max(30, a.budget - (time.time() - t0)), stdout=fh)
        if not ok:
            os.remove(paf + ".tmp")
            print(f"[selfaln] chunk {i} did not finish inside the budget; rerun (a chunk that never fits needs -f 0.001)")
            break
        os.replace(paf + ".tmp", paf)
        os.remove(q)
        done += 1
    print(f"[selfaln] {done}/{len(chunks)} chunks aligned")
    if done < len(chunks):
        return
    merged = rep.merge(rep.read_bed(a.repeats))
    pairs = []
    for paf in sorted(glob.glob(f"{d}/c*.paf")):
        for line in open(paf):
            p = from_selfpaf(line)
            if p and nonrepeat_aligned(p, merged) >= MIN_SD_BP:
                pairs.append(p)
    write_pairs(pairs, f"{a.outdir}/pairs.D3.tsv")
    print(f"[selfaln] {len(pairs)} pairs -> {a.outdir}/pairs.D3.tsv")


def cmd_d1(a):
    contigs = set(a.contigs.split(","))
    pairs = [p for p in (from_sedef(l, a.fmt) for l in open(a.sedef) if not l.startswith("#")) if p and p[0] in contigs and p[3] in contigs]
    write_pairs(pairs, a.out)
    print(f"[d1] {len(pairs)} pairs -> {a.out}")


def cmd_to_sedef(a):
    with open(a.out, "w") as fh:
        for p in read_pairs(a.pairs):
            if p[9]:
                fh.write(to_sedef_gorilla(p) + "\n")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("biser")
    p.add_argument("--genome", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--budget", type=int, default=560)
    p = sub.add_parser("selfaln")
    p.add_argument("--genome", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--repeats", required=True)
    p.add_argument("--chunk-bp", type=int, default=2_000_000)
    p.add_argument("--budget", type=int, default=540)
    p.add_argument("--threads", type=int, default=4)
    p = sub.add_parser("d1")
    p.add_argument("--sedef", required=True)
    p.add_argument("--fmt", required=True)
    p.add_argument("--contigs", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("meryl")
    p.add_argument("--genome", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--threads", type=int, default=4)
    p = sub.add_parser("to-sedef")
    p.add_argument("--pairs", required=True)
    p.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "meryl":
        cmd_meryl(a.genome, a.outdir, a.threads)
    else:
        {"biser": cmd_biser, "selfaln": cmd_selfaln, "d1": cmd_d1, "to-sedef": cmd_to_sedef}[a.cmd](a)


if __name__ == "__main__":
    main()
