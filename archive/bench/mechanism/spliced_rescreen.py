#!/usr/bin/env python3
"""L1 spliced identity re-screen.

The V4b candidate shortlist was selected on GENOMIC identity (whole copy span, introns
included).  The pipeline's admission gate 5 rejects a reconstructed copy whose best remap
identity is >= 0.98, and it operates on SPLICED (transcript) alignments.  Those are
different axes.  This script measures each family on the gate's axis.

L1(c_i)  -- for copy c_i of family F:
    reads  = primary (-F 2308) alignments in GGO_mm.bam whose reference interval is fully
             CONTAINED in [start, end) of c_i   (containment, NOT overlap)
    target = the GENOMIC sequence of the OTHER copies c_j, j != i   (family-restricted)
    map    = minimap2 -c --eqx -x splice -N 50 -p 0.1 targets_i.fa reads_i.fa
    per read: best hit by number of matching bases (PAF col 10); identity = 1 - de:f:
    L1(c_i)= median identity over that copy's reads THAT HAVE A HIT
L1(F)    = median over copies of L1(c_i)

Identity is ALWAYS 1 - de (minimap2's gap-compressed per-base divergence), never
nmatch/blocklen.  Reads with no hit are counted separately and are NOT folded in as 0.

Interpretation: L1 is a CEILING on what any reconstruction of that copy could score against
gate 5, because it uses the real reads rather than a reconstructed consensus.  L1 >= 0.98
=> the family can never be recovered by this method.  Well below 0.98 => gate-5 headroom.
"""
import argparse
import collections
import json
import os
import pathlib
import re
import statistics
import subprocess
import sys

MM2 = "/home/juanfra/miniforge3/bin/minimap2"
SAMTOOLS = "/home/juanfra/miniforge3/bin/samtools"
GGO_FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
GGO_BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
COPIES = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_gwcat.copies.tsv"
WORK = pathlib.Path("/home/juanfra/winloci_scratch/spliced_rescreen")
GATE = 0.98

_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def ref_len_of_cigar(cig):
    """Reference bases consumed.  N = intron spliced OUT -- it DOES consume reference."""
    return sum(int(n) for n, o in _CIG.findall(cig) if o in "MDN=X")


def run(cmd, **kw):
    p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)
    if p.returncode != 0:
        sys.stderr.write(" ".join(str(c) for c in cmd) + "\n" + p.stderr[-4000:] + "\n")
        p.check_returncode()
    return p


# ----------------------------------------------------------------------------- catalog
def load_copies(path, families):
    want = set(families)
    fams = collections.defaultdict(list)
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if want and r["family_id"] not in want:
                continue
            fams[r["family_id"]].append({
                "copy_idx": int(r["copy_idx"]), "tid": r["tid"], "chrom": r["chrom"],
                "start": int(r["start"]), "end": int(r["end"]),      # 0-based half-open
                "n_exon": int(r["n_exon"]), "strand": r["strand"],
                "n_reads_catalog": int(r["n_reads"]),
            })
    for f in fams:
        fams[f].sort(key=lambda c: c["copy_idx"])
    return dict(fams)


# ------------------------------------------------------------------------------- reads
def contained_reads(chrom, start, end, bam=GGO_BAM):
    """Primary alignments (-F 2308) fully contained in [start, end).  0-based half-open in,
    1-based inclusive region out.  Returns {read_name: SEQ} (deterministic by name)."""
    p = run([SAMTOOLS, "view", "-F", "2308", bam, f"{chrom}:{start + 1}-{end}"])
    out = {}
    for ln in p.stdout.splitlines():
        f = ln.split("\t")
        if f[5] == "*" or f[9] == "*":
            continue
        st = int(f[3]) - 1                      # SAM POS is 1-based -> 0-based
        en = st + ref_len_of_cigar(f[5])        # 0-based half-open end
        if st >= start and en <= end:
            out[f[0]] = f[9]
    return out


def write_fa(seqs, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        for n in sorted(seqs):
            fh.write(f">{n}\n{seqs[n]}\n")
    os.replace(tmp, path)
    return len(seqs)


def read_fa(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = ln[1:].strip().split()[0], []
            else:
                buf.append(ln.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def cached_reads(fam, copy, cache_dir):
    fa = cache_dir / fam / f"reads_{copy['copy_idx']}.fa"
    meta = fa.with_suffix(".meta.json")
    key = {"chrom": copy["chrom"], "start": copy["start"], "end": copy["end"],
           "bam": GGO_BAM, "rule": "contained_primary_F2308_v1"}
    if fa.exists() and meta.exists():
        try:
            if json.loads(meta.read_text()) == key:
                return read_fa(fa), fa
        except Exception:
            pass
    seqs = contained_reads(copy["chrom"], copy["start"], copy["end"])
    write_fa(seqs, fa)
    meta.write_text(json.dumps(key, sort_keys=True))
    return seqs, fa


def sibling_targets(fam, copies, i, cache_dir, fasta=GGO_FASTA):
    """Genomic sequence of every copy except index i, one multi-fasta."""
    fa = cache_dir / fam / f"targets_{copies[i]['copy_idx']}.fa"
    fa.parent.mkdir(parents=True, exist_ok=True)
    regions = []
    for j, c in enumerate(copies):
        if j == i:
            continue
        regions.append((c["copy_idx"], f"{c['chrom']}:{c['start'] + 1}-{c['end']}"))
    if not regions:
        return None, []
    tmp = fa.with_suffix(".fa.tmp")
    with open(tmp, "w") as fh:
        for cidx, reg in regions:
            r = run([SAMTOOLS, "faidx", fasta, reg])
            seq = "".join(r.stdout.splitlines()[1:])
            fh.write(f">copy{cidx}|{reg}\n")
            for k in range(0, len(seq), 80):
                fh.write(seq[k:k + 80] + "\n")
    os.replace(tmp, fa)
    return fa, [c for c, _ in regions]


# --------------------------------------------------------------------------- alignment
def best_hits(reads_fa, target_fa, preset="splice", extra=("-N", "50", "-p", "0.1"),
              threads=3):
    """Returns {read: (nmatch, identity, qcov, target)}.  identity = 1 - de:f:.
    Best hit chosen by PAF col 10 = number of matching bases."""
    cmd = [MM2, "-c", "--eqx", "-x", preset, *extra, "-t", str(threads),
           str(target_fa), str(reads_fa)]
    p = run(cmd)
    best = {}
    for ln in p.stdout.splitlines():
        f = ln.split("\t")
        if len(f) < 13:
            continue
        q, qlen, qs, qe, tgt = f[0], int(f[1]), int(f[2]), int(f[3]), f[5]
        nmatch = int(f[9])
        de = None
        for t in f[12:]:
            if t.startswith("de:f:"):
                de = float(t[5:])
                break
        if de is None:                 # requires -c; without it there is no de
            continue
        rec = (nmatch, 1.0 - de, (qe - qs) / qlen if qlen else 0.0, tgt)
        prev = best.get(q)
        # deterministic tie-break: higher nmatch, then higher identity, then target name
        if prev is None or (nmatch, rec[1], tgt) > (prev[0], prev[1], prev[3]):
            best[q] = rec
    return best


def stats(ids, qcovs):
    if not ids:
        return {"median_identity": None, "min_identity": None, "max_identity": None,
                "frac_ge_098": None, "n_ge_098": 0, "n_lt_098": 0, "median_qcov": None}
    return {"median_identity": round(statistics.median(ids), 6),
            "min_identity": round(min(ids), 6), "max_identity": round(max(ids), 6),
            "frac_ge_098": round(sum(i >= GATE for i in ids) / len(ids), 6),
            "n_ge_098": sum(i >= GATE for i in ids),
            "n_lt_098": sum(i < GATE for i in ids),
            "median_qcov": round(statistics.median(qcovs), 6)}


# ------------------------------------------------------------------------------ driver
def measure_family(fam, copies, cache_dir, threads=3, preset="splice"):
    out = {"family_id": fam, "n_copies": len(copies), "copies": [], "notes": []}
    if len(copies) < 2:
        out["notes"].append("fewer than 2 copies -- family-restricted L1 undefined")
        out.update({"L1_family": None, "L1_min": None, "L1_max": None})
        return out
    pooled_ids, pooled_qc = [], []
    for i, c in enumerate(copies):
        seqs, reads_fa = cached_reads(fam, c, cache_dir)
        tgt_fa, tgt_idx = sibling_targets(fam, copies, i, cache_dir)
        rec = {"copy_idx": c["copy_idx"], "tid": c["tid"], "chrom": c["chrom"],
               "start": c["start"], "end": c["end"], "span_bp": c["end"] - c["start"],
               "strand": c["strand"], "n_exon": c["n_exon"],
               "n_reads_catalog": c["n_reads_catalog"],
               "n_reads": len(seqs), "target_copies": tgt_idx,
               "reads_fa": str(reads_fa), "targets_fa": str(tgt_fa)}
        if not seqs:
            rec.update({"n_reads_with_hit": 0, "n_reads_no_hit": 0,
                        "target_preference": {}, **stats([], [])})
            rec["L1"] = None
            rec["note"] = "no reads fully contained in this copy's span"
            out["copies"].append(rec)
            continue
        best = best_hits(reads_fa, tgt_fa, preset=preset, threads=threads)
        ids = [v[1] for v in best.values()]
        qcs = [v[2] for v in best.values()]
        pref = collections.Counter(v[3] for v in best.values())
        rec.update({"n_reads_with_hit": len(best),
                    "n_reads_no_hit": len(seqs) - len(best),
                    "target_preference": dict(sorted(pref.items())), **stats(ids, qcs)})
        rec["L1"] = rec["median_identity"]
        pooled_ids += ids
        pooled_qc += qcs
        out["copies"].append(rec)
    per_copy = [c["L1"] for c in out["copies"] if c["L1"] is not None]
    out["n_copies_measured"] = len(per_copy)
    out["L1_family"] = round(statistics.median(per_copy), 6) if per_copy else None
    out["L1_min"] = round(min(per_copy), 6) if per_copy else None
    out["L1_max"] = round(max(per_copy), 6) if per_copy else None
    out["n_reads_total"] = sum(c["n_reads"] for c in out["copies"])
    out["n_reads_with_hit_total"] = sum(c.get("n_reads_with_hit", 0) for c in out["copies"])
    out["n_reads_no_hit_total"] = sum(c.get("n_reads_no_hit", 0) for c in out["copies"])
    out["pooled"] = stats(pooled_ids, pooled_qc)
    if per_copy:
        out["verdict"] = ("NO gate-5 headroom (L1 >= 0.98)" if out["L1_family"] >= GATE
                          else "gate-5 headroom (L1 < 0.98)")
    else:
        out["verdict"] = "UNMEASURABLE"
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", required=True, help="comma-separated family ids")
    ap.add_argument("--out", required=True, help="output json path")
    ap.add_argument("--copies", default=COPIES)
    ap.add_argument("--cache", default=str(WORK / "cache"))
    ap.add_argument("--preset", default="splice",
                    help="minimap2 -x preset (default splice; :hq gives DIFFERENT numbers)")
    ap.add_argument("--threads", type=int, default=3)
    a = ap.parse_args()

    fams = [f for f in a.families.split(",") if f]
    cache_dir = pathlib.Path(a.cache)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cat = load_copies(a.copies, fams)

    results = {"_meta": {"copies_tsv": a.copies, "bam": GGO_BAM, "fasta": GGO_FASTA,
                         "preset": a.preset, "gate": GATE,
                         "identity": "1 - de:f: (gap-compressed), best hit by PAF col 10",
                         "read_rule": "primary (-F 2308), fully CONTAINED in copy span",
                         "geometry": "family-restricted: reads vs sibling copy genomic seqs"},
               "families": {}}
    for f in fams:
        if f not in cat:
            results["families"][f] = {"family_id": f, "error": "not in copies.tsv"}
            print(f"{f}\tERROR not in copies.tsv")
            continue
        r = measure_family(f, cat[f], cache_dir, threads=a.threads, preset=a.preset)
        results["families"][f] = r
        print(f"{f}\tn_copies={r['n_copies']}\tL1={r['L1_family']}\t"
              f"[{r['L1_min']},{r['L1_max']}]\treads={r.get('n_reads_total')}\t"
              f"hit={r.get('n_reads_with_hit_total')}\tnohit={r.get('n_reads_no_hit_total')}\t"
              f"pooled_frac_ge_098={r['pooled']['frac_ge_098'] if 'pooled' in r else None}\t"
              f"{r.get('verdict')}", flush=True)

    outp = pathlib.Path(a.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(results, indent=2, sort_keys=True))
    print(f"# wrote {outp}")


if __name__ == "__main__":
    main()
