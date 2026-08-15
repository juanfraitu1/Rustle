"""V4b removal-recovery for 3-COPY families, at a lowered gate-1 cluster floor.

Why this driver exists (and why it is NOT real_tandem.py)
--------------------------------------------------------
`real_tandem.py` requires >=4 supported copies because it deletes TWO copies so that
host+2 = 3 clusters can reach the SHIPPING gate-1 floor (min_clusters=3). The two loci that
stopped architecturally in the prior V4b run (GSTM, GWFAM227) are bare 3-copy families:
deleting one copy leaves at most host + deleted = 2 clusters, so the shipping gate is
MATHEMATICALLY UNREACHABLE there, for a reason that has nothing to do with divergence.

This driver deletes ONE copy of a 3-copy family and runs the admission gate at
RUSTLE_ABSENT_MIN_CLUSTERS=2 (`--absent-min-clusters 2`), which is the value that makes the
gate reachable for this arithmetic. THE DEFAULT REMAINS 3 -- see absent_copy.rs.

THE CONTROL IS NOT OPTIONAL
---------------------------
Lowering a gate that guards copy-vs-allele is exactly how a false positive is manufactured.
So every ablation here is PAIRED with an identical run against the INTACT (unmasked)
region-local assembly, same reads, same min_clusters. Both are always reported. If the
no-deletion control also admits an absent copy, the "recovery" is an artefact of the lowered
gate and the result is void.

Everything except the mask is held identical between the two arms: same extracted read set
(one fastq, aligned twice), same padded window, same contig length, same copy_assign command.
"""
import argparse
import json
import os
import pathlib
import shutil
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from sim_tandem import set_match  # noqa: E402

GGO_FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
GGO_BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
MM2 = shutil.which("minimap2") or "/home/juanfra/miniforge3/bin/minimap2"
SAMTOOLS = shutil.which("samtools") or "/home/juanfra/miniforge3/bin/samtools"
COPY_ASSIGN = str(pathlib.Path(__file__).resolve().parents[2] / "target/release/copy_assign")


def _run(cmd, env=None, log=None):
    print("+", " ".join(str(c) for c in cmd), file=sys.stderr, flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if log is not None:
        pathlib.Path(log).write_text((r.stdout or "") + "\n===STDERR===\n" + (r.stderr or ""))
    if r.stderr:
        print(r.stderr[-3000:], file=sys.stderr, flush=True)
    r.check_returncode()
    return r


def _read_tsv(path):
    p = pathlib.Path(path)
    if not p.exists():
        return []
    with open(p) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(dict(zip(header, line.split("\t"))))
        return rows


def run_copy_assign(bam, fasta, region, out_prefix, extra_args=(), min_clusters=None, log=None):
    cmd = [
        COPY_ASSIGN, "--bam", str(bam), "--fasta", str(fasta), "--region", region,
        "--out", str(out_prefix), "--homology-primary", "--min-copies", "2",
        *extra_args,
    ]
    env = dict(os.environ)
    env.pop("RUSTLE_ABSENT_MIN_CLUSTERS", None)  # never inherit; set explicitly or not at all
    if min_clusters is not None:
        cmd += ["--absent-min-clusters", str(min_clusters)]
    r = _run(cmd, env=env, log=log)
    return r.stdout + r.stderr


def load_psv(prefix):
    cols = {int(r["col_index"]): int(r["genome_pos"]) for r in _read_tsv(f"{prefix}.psv_cols.tsv")}
    rows = _read_tsv(f"{prefix}.psv_copies.tsv")
    copies = {int(r["copy_index"]): r["alleles"] for r in rows}
    tids = {int(r["copy_index"]): r["copy_tid"] for r in rows}
    return cols, copies, tids


def identity(a, b):
    n = same = 0
    for x, y in zip(a, b):
        if x == "." or y == ".":
            continue
        n += 1
        same += (x == y)
    return (same / n if n else 0.0), n


def pick_single_deletion(quant_rows, copies):
    """3-copy design: choose (keep, host, hidden) with ONE hidden copy that is strictly closer
    (PSV-allele identity) to `host` than to the other copy that stays in the reference, so its
    reads pool onto host -> host + hidden = 2 clusters at the surviving locus."""
    idxs = sorted(int(r["copy_index"]) for r in quant_rows)
    best = None
    for host in idxs:
        for hidden in idxs:
            if hidden == host:
                continue
            keep = [i for i in idxs if i != hidden]
            id_host, n_host = identity(copies[host], copies[hidden])
            ok = n_host > 0
            for k in keep:
                if k == host:
                    continue
                id_k, _ = identity(copies[k], copies[hidden])
                if id_k >= id_host:
                    ok = False
            if ok and (best is None or id_host > best[0]):
                best = (id_host, host, hidden, keep, n_host)
    if best is None:
        return None
    return {"host": best[1], "hidden": [best[2]], "keep": best[3],
            "psv_identity_hidden_to_host": best[0], "n_psv_cols_compared": best[4]}


def build_local_ref(fasta, contig, pad_start, pad_end, mask_spans, out_fa, new_contig):
    r = _run([SAMTOOLS, "faidx", fasta, f"{contig}:{pad_start}-{pad_end}"])
    seq = list("".join(r.stdout.splitlines()[1:]))
    n_masked = 0
    for (s, e) in mask_spans:
        rel_s, rel_e = max(0, s - pad_start), min(len(seq), e - pad_start + 1)
        for i in range(rel_s, rel_e):
            if seq[i] != "N":
                n_masked += 1
            seq[i] = "N"
    with open(out_fa, "w") as f:
        f.write(f">{new_contig}\n")
        s = "".join(seq)
        for i in range(0, len(s), 80):
            f.write(s[i:i + 80] + "\n")
    _run([SAMTOOLS, "faidx", str(out_fa)])
    return len(seq), n_masked


def extract_reads(bam, contig, pad_start, pad_end, work, prefix):
    raw = work / f"{prefix}.extract.bam"
    _run([SAMTOOLS, "view", "-b", "-F", "2308", bam, f"{contig}:{pad_start}-{pad_end}", "-o", str(raw)])
    fq = work / f"{prefix}.fastq"
    r = _run([SAMTOOLS, "fastq", str(raw)])
    fq.write_text(r.stdout)
    return fq, sum(1 for _ in open(fq)) // 4


def align(fq, ref_fa, work, prefix):
    p = subprocess.run([MM2, "-ax", "splice:hq", "-N", "50", str(ref_fa), str(fq)],
                       capture_output=True, text=True)
    print(p.stderr[-1500:], file=sys.stderr, flush=True)
    p.check_returncode()
    sam = work / f"{prefix}.sam"
    sam.write_text(p.stdout)
    srt = work / f"{prefix}.sorted.bam"
    _run([SAMTOOLS, "sort", "-o", str(srt), str(sam)])
    _run([SAMTOOLS, "index", str(srt)])
    n_mapped = int(subprocess.run([SAMTOOLS, "view", "-c", "-F", "2308", str(srt)],
                                  capture_output=True, text=True).stdout.strip() or 0)
    return srt, n_mapped


def arm(label, ref_fa, new_contig, ref_len, fq, work, min_clusters):
    """Run one arm (ablation or control) and summarise its absent-copy outcome."""
    srt, n_mapped = align(fq, ref_fa, work, label)
    prefix = work / f"{label}.ca"
    log = work / f"{label}.ca.log"
    run_copy_assign(srt, ref_fa, f"{new_contig}:1-{ref_len}", prefix,
                    extra_args=["--absent-copies", "--dump-psv"],
                    min_clusters=min_clusters, log=str(log))
    fams = _read_tsv(f"{prefix}.families.tsv")
    dna = _read_tsv(f"{prefix}.dna_needs.tsv")
    psv_path = pathlib.Path(f"{prefix}.psv_copies.tsv")
    admitted = []
    cols = copies = tids = None
    if psv_path.exists():
        cols, copies, tids = load_psv(str(prefix))
        admitted = sorted(i for i, t in tids.items() if t.startswith("AC_"))
    notice = "[absent_copy] " in pathlib.Path(log).read_text()
    return {
        "label": label,
        "ref": str(ref_fa),
        "n_reads_realigned_primary": n_mapped,
        "n_families": len(fams),
        "families_n_copies": [int(f["n_copies"]) for f in fams],
        "n_dna_needs": len(dna),
        "dna_needs": [{"n_clusters": int(r["n_clusters"]), "reason": r["reason"],
                       "read_count": int(r["read_count"]),
                       "span": f'{r["chrom"]}:{r["start"]}-{r["end"]}'} for r in dna],
        "n_clusters_seen": sorted({int(r["n_clusters"]) for r in dna}),
        "max_n_clusters_seen": max([int(r["n_clusters"]) for r in dna], default=None),
        "n_admitted_absent_copies": len(admitted),
        "admitted_indices": admitted,
        "admitted_tids": [tids[i] for i in admitted] if admitted else [],
        "gate1_rejections": [r["reason"] for r in dna if "clusters (copy-vs-allele" in r["reason"]],
        "min_clusters_notice_printed": notice,
        "_psv": (cols, copies, tids) if psv_path.exists() else None,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--locus", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--pad", type=int, default=3000)
    ap.add_argument("--min-clusters", type=int, default=2)
    ap.add_argument("--scratch", default="/home/juanfra/winloci_scratch/v4b_mincl2")
    a = ap.parse_args()

    work = pathlib.Path(a.scratch) / a.label
    work.mkdir(parents=True, exist_ok=True)
    contig, span = a.locus.split(":")
    start, end = (int(x) for x in span.split("-"))
    pad_start, pad_end = start - a.pad, end + a.pad

    res = {"label": a.label, "locus": a.locus, "min_clusters": a.min_clusters,
           "pad": a.pad, "scratch": str(work)}

    # ---- Step 1: same-read-set truth against the INTACT full assembly -------------------
    intact = work / "intact"
    run_copy_assign(GGO_BAM, GGO_FASTA, a.locus, intact, extra_args=["--dump-psv"],
                    log=str(work / "intact.log"))
    fams = _read_tsv(f"{intact}.families.tsv")
    quant = _read_tsv(f"{intact}.quant.tsv")
    res["intact_n_families"] = len(fams)
    res["intact_n_copies_supported"] = int(fams[0]["n_copies"]) if fams else 0
    res["intact_copies"] = [{"copy_index": int(r["copy_index"]), "copy_tid": r.get("copy_tid"),
                             "start": int(r["copy_start"]), "end": int(r["copy_end"]),
                             "n_reads_hard": r.get("n_reads_hard"),
                             "anchored_reads": r.get("anchored_reads")} for r in quant]
    if len(quant) < 3:
        res["stop"] = f"intact run supports {len(quant)} copies (<3); no 1-of-3 deletion design"
        print(json.dumps(res, indent=2, default=str))
        return

    cols, copies, tids = load_psv(str(intact))
    design = pick_single_deletion(quant, copies)
    if design is None:
        res["stop"] = ("no valid 1-hidden design: no copy is strictly closer (PSV-allele identity) "
                       "to one host than to the other surviving copy")
        print(json.dumps(res, indent=2, default=str))
        return
    res["design"] = design
    span_by_idx = {int(r["copy_index"]): (int(r["copy_start"]), int(r["copy_end"])) for r in quant}
    mask_spans = [span_by_idx[h] for h in design["hidden"]]
    res["mask_spans"] = mask_spans

    # ---- Step 2/3: two region-local references, identical but for the mask -------------
    new_contig = f"{contig}_local"
    abl_fa = work / "ref_ablation.fa"
    ctl_fa = work / "ref_control.fa"
    ref_len, n_masked = build_local_ref(GGO_FASTA, contig, pad_start, pad_end, mask_spans,
                                        abl_fa, new_contig)
    ref_len_c, n_masked_c = build_local_ref(GGO_FASTA, contig, pad_start, pad_end, [],
                                            ctl_fa, new_contig)
    assert ref_len == ref_len_c and n_masked_c == 0
    res["local_ref_len"] = ref_len
    res["bases_masked"] = n_masked

    # ---- Step 4: ONE read extraction, aligned to both refs ------------------------------
    fq, n_reads = extract_reads(GGO_BAM, contig, pad_start, pad_end, work, "reads")
    res["n_reads_extracted"] = n_reads
    if n_reads == 0:
        res["stop"] = "zero reads extracted -- BUG, not a finding"
        print(json.dumps(res, indent=2, default=str))
        return

    # ---- Step 5: ablation, then paired control, at the SAME min_clusters ----------------
    abl = arm("ablation", abl_fa, new_contig, ref_len, fq, work, a.min_clusters)
    ctl = arm("control", ctl_fa, new_contig, ref_len, fq, work, a.min_clusters)

    # ---- Step 6: identity_to_intact of any admitted copy vs the DELETED copy ------------
    def score(armres):
        if not armres["_psv"] or not armres["admitted_indices"]:
            return None, []
        d_cols, d_copies, _ = armres["_psv"]
        offset = pad_start - 1
        d_abs = {ci: pos + offset for ci, pos in d_cols.items()}
        intact_pos_to_col = {pos: ci for ci, pos in cols.items()}
        shared = sorted({p for p in d_abs.values() if p in intact_pos_to_col})
        if not shared:
            return None, []
        d_pos_to_col = {v: k for k, v in d_abs.items()}
        truth = ["".join(copies[h][intact_pos_to_col[p]] for p in shared) for h in design["hidden"]]
        rec = ["".join(d_copies[i][d_pos_to_col[p]] for p in shared) for i in armres["admitted_indices"]]
        m = set_match(rec, truth) if rec and truth else []
        return (max((x[2] for x in m), default=None), m)

    for armres in (abl, ctl):
        idv, m = score(armres)
        armres["identity_to_intact_deleted_copy"] = idv
        armres["set_match"] = m
        armres["n_shared_psv_positions"] = len(
            {p + pad_start - 1 for p in armres["_psv"][0].values()} & set(cols.values())
        ) if armres["_psv"] else 0
        armres.pop("_psv", None)

    res["ablation"] = abl
    res["control"] = ctl
    res["control_also_admitted"] = ctl["n_admitted_absent_copies"] > 0
    res["recovered"] = bool(
        abl["n_admitted_absent_copies"] > 0
        and not res["control_also_admitted"]
        and (abl["identity_to_intact_deleted_copy"] or 0) > 0.90
    )
    (work / "result.json").write_text(json.dumps(res, indent=2, default=str))
    print(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    main()
