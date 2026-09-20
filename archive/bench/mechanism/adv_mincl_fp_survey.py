"""ADVERSARIAL: what does min_clusters=2 COST on loci where NO copy was deleted?

The V4b ablation lowered gate 1 (`n_clusters >= min_clusters`) from 3 to 2 on the argument that
in a removal ablation copy status is established BY CONSTRUCTION. That argument does not license
the knob anywhere else. This driver measures the knob's false-positive cost directly:

    for many INTACT loci (no deletion, nothing masked), run the absent-copy admission gate at
    min_clusters=3 (the shipping default) and at min_clusters=2, on the SAME alignment, and count
    admissions that appear only at 2.

Everything except the flag is held identical: one region-local reference (UNMASKED), one read
extraction, ONE alignment, two copy_assign invocations over the same sorted BAM. Any difference
in output is therefore attributable to the flag alone.

An admission on an intact assembly is not automatically a false positive (a genuinely collapsed
copy may exist there), so the reported quantity is the DELTA -- admissions that gate 1 was
holding back at the shipping default -- plus, for each such admission, how host-like it is
(PSV-allele identity to its own host copy). Gate 1 is named for copy-vs-allele confusion; an
extra admission whose alleles are near-identical to its host is exactly the failure mode.

Serial, foreground, resumable (`--limit N` per invocation), outputs under winloci_scratch.
"""
import argparse
import csv
import filecmp
import json
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import v4b_3copy  # noqa: E402
from v4b_3copy import (  # noqa: E402
    GGO_BAM, GGO_FASTA, _read_tsv, align, build_local_ref, extract_reads,
    identity, load_psv, run_copy_assign,
)

# PIN the binary. `target/release/copy_assign` is a SHARED path that another process can (and did,
# mid-run on 2026-08-08 11:22) replace with a build from different sources -- a survey whose two
# arms ran on two different binaries measures nothing. The pinned copy is verified to (a) accept
# --absent-min-clusters and (b) come from a source tree whose md5s equal the current ones.
PINNED = "/home/juanfra/winloci_scratch/adv_mincl_fp/copy_assign.pinned"
if pathlib.Path(PINNED).exists():
    v4b_3copy.COPY_ASSIGN = PINNED

SUFFIXES = ["assignments", "quant", "families", "psv_copies", "psv_cols", "dna_needs",
            "famcn_readonly", "psv_reads"]


def arm(label, ref_fa, new_contig, ref_len, srt, work, min_clusters):
    prefix = work / f"{label}.ca"
    log = work / f"{label}.ca.log"
    run_copy_assign(srt, ref_fa, f"{new_contig}:1-{ref_len}", prefix,
                    extra_args=["--absent-copies", "--dump-psv"],
                    min_clusters=min_clusters, log=str(log))
    fams = _read_tsv(f"{prefix}.families.tsv")
    dna = _read_tsv(f"{prefix}.dna_needs.tsv")
    admitted, host_ident = [], []
    psv_path = pathlib.Path(f"{prefix}.psv_copies.tsv")
    if psv_path.exists():
        _cols, copies, tids = load_psv(str(prefix))
        admitted = sorted(i for i, t in tids.items() if t.startswith("AC_"))
        # An AC_ tid is AC_<host_tid_locus>; its host is the DN_ copy at the same locus string.
        for i in admitted:
            stem = tids[i][len("AC_"):]
            hosts = [j for j, t in tids.items()
                     if t.startswith("DN_") and t[len("DN_"):].rsplit("_", 1)[0] == stem]
            best = None
            for j in hosts:
                idv, n = identity(copies[i], copies[j])
                if n and (best is None or idv > best[0]):
                    best = (idv, n, tids[j])
            host_ident.append({"ac_tid": tids[i], "copy_index": i,
                               "identity_to_host": (best[0] if best else None),
                               "n_psv_cols": (best[1] if best else 0),
                               "host_tid": (best[2] if best else None)})
    return {
        "label": label,
        "min_clusters": min_clusters,
        "n_families": len(fams),
        "families_n_copies": [int(f["n_copies"]) for f in fams],
        "collapsed_copies": [int(f["collapsed_copies"]) for f in fams],
        "n_dna_needs": len(dna),
        "dna_needs": [{"n_clusters": int(r["n_clusters"]), "reason": r["reason"],
                       "read_count": int(r["read_count"]),
                       "span": f'{r["chrom"]}:{r["start"]}-{r["end"]}'} for r in dna],
        "n_admitted_absent_copies": len(admitted),
        "admitted_tids": [h["ac_tid"] for h in host_ident],
        "admitted_host_identity": host_ident,
        "notice_printed": "[absent_copy] " in pathlib.Path(log).read_text(),
    }


def one_locus(fid, chrom, start, end, pad, scratch):
    work = pathlib.Path(scratch) / fid
    work.mkdir(parents=True, exist_ok=True)
    out = work / "result.json"
    if out.exists():
        return json.loads(out.read_text()), True
    t0 = time.time()
    pad_start, pad_end = max(1, start - pad), end + pad
    new_contig = f"{chrom}_local"
    ref_fa = work / "ref_intact.fa"
    ref_len, n_masked = build_local_ref(GGO_FASTA, chrom, pad_start, pad_end, [], ref_fa,
                                        new_contig)
    assert n_masked == 0, "control reference must mask nothing"
    fq, n_reads = extract_reads(GGO_BAM, chrom, pad_start, pad_end, work, "reads")
    res = {"family_id": fid, "locus": f"{chrom}:{start}-{end}", "pad": pad,
           "local_ref_len": ref_len, "n_reads_extracted": n_reads, "deletion": "NONE (intact)"}
    if n_reads == 0:
        res["stop"] = "zero reads extracted"
        out.write_text(json.dumps(res, indent=2))
        return res, False
    srt, n_mapped = align(fq, ref_fa, work, "reads")
    res["n_reads_realigned_primary"] = n_mapped
    res["binary"] = v4b_3copy.COPY_ASSIGN
    # SAME alignment, two flag values. default arm first (knob unset -> no env var at all).
    try:
        res["default3"] = arm("default3", ref_fa, new_contig, ref_len, srt, work, None)
        res["mincl2"] = arm("mincl2", ref_fa, new_contig, ref_len, srt, work, 2)
    except Exception as exc:  # record and move on; a crashed locus is not a silent zero
        res["stop"] = f"copy_assign failed: {type(exc).__name__}: {str(exc)[:300]}"
        out.write_text(json.dumps(res, indent=2))
        return res, False
    res["byte_identical"] = {
        s: filecmp.cmp(str(work / f"default3.ca.{s}.tsv"), str(work / f"mincl2.ca.{s}.tsv"),
                       shallow=False)
        if (work / f"default3.ca.{s}.tsv").exists() and (work / f"mincl2.ca.{s}.tsv").exists()
        else None
        for s in SUFFIXES
    }
    res["extra_admissions_at_2"] = sorted(
        set(res["mincl2"]["admitted_tids"]) - set(res["default3"]["admitted_tids"]))
    res["lost_admissions_at_2"] = sorted(
        set(res["default3"]["admitted_tids"]) - set(res["mincl2"]["admitted_tids"]))
    res["seconds"] = round(time.time() - t0, 1)
    out.write_text(json.dumps(res, indent=2))
    return res, False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci-tsv", required=True)
    ap.add_argument("--pad", type=int, default=3000)
    ap.add_argument("--limit", type=int, default=5)
    ap.add_argument("--scratch", default="/home/juanfra/winloci_scratch/adv_mincl_fp")
    a = ap.parse_args()

    todo = list(csv.DictReader(open(a.loci_tsv), delimiter="\t"))
    done_n = 0
    for row in todo:
        if done_n >= a.limit:
            break
        fid = row["family_id"]
        res, cached = one_locus(fid, row["chrom"], int(row["start"]), int(row["end"]),
                                a.pad, a.scratch)
        if cached:
            continue
        done_n += 1
        d3 = res.get("default3", {}).get("n_admitted_absent_copies")
        m2 = res.get("mincl2", {}).get("n_admitted_absent_copies")
        print(f"{fid}\treads={res.get('n_reads_extracted')}\tadmit3={d3}\tadmit2={m2}"
              f"\textra={res.get('extra_admissions_at_2')}\t{res.get('stop','')}"
              f"\t{res.get('seconds','')}s", flush=True)


if __name__ == "__main__":
    main()
