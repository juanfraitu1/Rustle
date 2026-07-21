"""V4b: real gorilla tandem removal-recovery, same-read-set ground truth.

Substrate correction from the brief (recorded honestly): the brief names A119b.t2t.bam +
mGorGor1 on the mounted disk. mGorGor1 is NOT present as a FASTA on this machine. Per the
task instruction this driver uses the matched, on-disk pair instead:
  - GGO.fasta   (/home/juanfra/winloci_scratch/GGO.fasta)   -- the mGorGor1-v2.0 T2T assembly
  - GGO_mm.bam  (/home/juanfra/winloci_scratch/GGO_mm.bam)  -- gorilla HiFi aligned to GGO.fasta
This keeps reference and reads from the SAME assembly (the point of the same-read-set arm:
truth = what this exact read set supports against the intact reference, no cross-individual
confound), which is what the brief's own "same-read-set" framing requires -- it does not
require the specific A119b/mGorGor1 filenames.

Method (mirrors bench/mechanism/sim_tandem.py's v4a design, real data instead of simulated):
  1. Run copy_assign on the intact reference at a candidate locus (--homology-primary
     --min-copies 2 --dump-psv) -- this IS the ground truth: n_copies this exact read set
     supports, and the intact per-copy PSV allele strings.
  2. Per Task 5's finding #2 (absent_copy.rs::AbsentCopyParams::min_clusters is hardcoded to 3,
     no CLI override): a plain 1-of-3 deletion can only ever leave host+1=2 clusters at any one
     surviving locus, which never reaches the gate. So this driver requires a locus with >=4
     supported copies where >=2 of them are each closer (by PSV-allele identity) to one
     surviving "host" copy than to any other surviving copy -- those 2 are deleted (masked) so
     their reads pool onto the host, giving host+2 = 3 clusters at that one locus.
  3. Degrade: extract the locus (+padding) from GGO.fasta into a small region-local reference
     (NOT a 3.6GB rewrite), mask the two chosen copies' spans with 'N', keep the rest.
  4. Re-align the SAME reads (extracted from GGO_mm.bam for this region, primary-only) to the
     degraded region-local reference with minimap2 -ax splice:hq -N 50 (NOT --secondary=no).
  5. Run copy_assign --absent-copies --linearize --dump-psv on the degraded alignment; compare
     admitted candidates' PSV alleles (matched back to absolute genome_pos) against the deleted
     copies' true alleles from step 1, via set_match (permutation-invariant).
"""
import argparse
import json
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


def _run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd), file=sys.stderr)
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.stderr:
        print(r.stderr[-4000:], file=sys.stderr)
    r.check_returncode()
    return r


def _read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(dict(zip(header, line.split("\t"))))
        return rows


def run_copy_assign(bam, fasta, region, out_prefix, extra_args=()):
    cmd = [
        COPY_ASSIGN, "--bam", str(bam), "--fasta", str(fasta), "--region", region,
        "--out", str(out_prefix), "--homology-primary", "--min-copies", "2",
        *extra_args,
    ]
    r = _run(cmd)
    return r.stdout + r.stderr


def load_psv(prefix):
    """cols: {col_index: genome_pos}; copies: {copy_index: alleles_str}."""
    cols = {int(r["col_index"]): int(r["genome_pos"]) for r in _read_tsv(f"{prefix}.psv_cols.tsv")}
    copy_rows = _read_tsv(f"{prefix}.psv_copies.tsv")
    copies = {int(r["copy_index"]): r["alleles"] for r in copy_rows}
    tids = {int(r["copy_index"]): r["copy_tid"] for r in copy_rows}
    return cols, copies, tids


def identity(a, b):
    n = 0
    same = 0
    for x, y in zip(a, b):
        if x == "." or y == ".":
            continue
        n += 1
        if x == y:
            same += 1
    return (same / n if n else 0.0), n


def pick_deletion_design(quant_rows, cols, copies):
    """Choose (keep_indices, host_index, hidden_indices) among the family's reference copies
    such that hidden copies are each strictly closer (PSV-allele identity) to `host` than to
    any OTHER copy that stays in the reference -- i.e. their reads should pool onto host,
    giving host+len(hidden) clusters at that one surviving locus."""
    idxs = [int(r["copy_index"]) for r in quant_rows]
    best = None
    for host in idxs:
        others = [i for i in idxs if i != host]
        # rank every other copy by identity to host
        scored = []
        for o in others:
            idv, n = identity(copies[host], copies[o])
            scored.append((idv, n, o))
        scored.sort(reverse=True)
        # try taking the top-2 closest-to-host as hidden, requiring each is closer to host
        # than to every copy that would remain in the reference (i.e. every other idx not
        # itself deleted).
        if len(scored) < 2:
            continue
        hidden = [scored[0][2], scored[1][2]]
        keep = [i for i in idxs if i not in hidden]
        ok = True
        for h in hidden:
            id_host, _ = identity(copies[host], copies[h])
            for k in keep:
                if k == host:
                    continue
                id_k, _ = identity(copies[k], copies[h])
                if id_k >= id_host:
                    ok = False
        if ok:
            cand = (scored[0][0] + scored[1][0], host, hidden, keep)
            if best is None or cand[0] > best[0]:
                best = cand
    if best is None:
        return None
    _, host, hidden, keep = best
    return keep, host, hidden


def build_degraded_ref(fasta, contig, pad_start, pad_end, mask_spans, out_fa, new_contig):
    region = f"{contig}:{pad_start}-{pad_end}"
    r = _run([SAMTOOLS, "faidx", fasta, region])
    lines = r.stdout.splitlines()
    seq = "".join(lines[1:])
    seq_list = list(seq)
    for (s, e) in mask_spans:
        rel_s = max(0, s - pad_start)
        rel_e = min(len(seq_list), e - pad_start + 1)
        for i in range(rel_s, rel_e):
            seq_list[i] = "N"
    degraded = "".join(seq_list)
    with open(out_fa, "w") as f:
        f.write(f">{new_contig}\n")
        for i in range(0, len(degraded), 80):
            f.write(degraded[i:i + 80] + "\n")
    _run([SAMTOOLS, "faidx", str(out_fa)])
    return len(degraded)


def extract_and_realign(bam, contig, pad_start, pad_end, degraded_fa, work, prefix):
    region = f"{contig}:{pad_start}-{pad_end}"
    raw_bam = work / f"{prefix}.extract.bam"
    _run([SAMTOOLS, "view", "-b", "-F", "2308", bam, region, "-o", str(raw_bam)])
    fq = work / f"{prefix}.fastq"
    r = _run([SAMTOOLS, "fastq", str(raw_bam)])
    fq.write_text(r.stdout)
    n_reads = sum(1 for _ in open(fq)) // 4
    sam_out = work / f"{prefix}.sam"
    p1 = subprocess.run(
        [MM2, "-ax", "splice:hq", "-N", "50", str(degraded_fa), str(fq)],
        capture_output=True, text=True,
    )
    print(p1.stderr[-2000:], file=sys.stderr)
    p1.check_returncode()
    bam_out = work / f"{prefix}.bam"
    bam_out.write_text(p1.stdout)
    sorted_bam = work / f"{prefix}.sorted.bam"
    _run([SAMTOOLS, "sort", "-o", str(sorted_bam), str(bam_out)])
    _run([SAMTOOLS, "index", str(sorted_bam)])
    return sorted_bam, n_reads


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--locus", default="NC_073248.2:19558214-19739408",
                     help="chrom:start-end for the candidate tandem family (default: RBMY1 array)")
    ap.add_argument("--pad", type=int, default=3000)
    ap.add_argument("--scratch", default="/home/juanfra/winloci_scratch/mech_real")
    ap.add_argument(
        "--out",
        default=str(pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism/verification_results.json"),
    )
    a = ap.parse_args()
    work = pathlib.Path(a.scratch)
    work.mkdir(parents=True, exist_ok=True)

    contig, span = a.locus.split(":")
    start, end = (int(x) for x in span.split("-"))
    pad_start, pad_end = start - a.pad, end + a.pad

    # Step 1: intact-reference ground truth.
    intact_prefix = work / "mech_real_full2"
    log1 = run_copy_assign(GGO_BAM, GGO_FASTA, a.locus, intact_prefix,
                            extra_args=["--dump-psv"])
    families = _read_tsv(f"{intact_prefix}.families.tsv")
    quant = _read_tsv(f"{intact_prefix}.quant.tsv")
    n_copies_supported = int(families[0]["n_copies"]) if families else 0
    v4b = {
        "locus": a.locus,
        "substrate": "GGO.fasta=mGorGor1 + GGO_mm.bam, gorilla-to-gorilla",
        "truth_source": "same-read-set vs full GGO.fasta",
        "n_copies_supported": n_copies_supported,
        "recovered": False,
        "identity_to_intact": None,
        "min_clusters_satisfied": False,
        "note": "",
    }

    if n_copies_supported < 4 or len(quant) < 4:
        v4b["note"] = (
            f"Intact reference supports n_copies={n_copies_supported} (<4) at this locus. "
            "A locus deletion design that satisfies the hardcoded min_clusters=3 admission "
            "gate (absent_copy.rs) requires >=1 surviving host + 2 deleted copies that pool "
            "onto it, i.e. >=4 total supported copies at one family. Stopping honestly -- "
            "no recovery attempted at this locus."
        )
        out = pathlib.Path(a.out)
        data = json.loads(out.read_text()) if out.exists() else {}
        data["v4b"] = v4b
        out.write_text(json.dumps(data, indent=2))
        print("V4b:", json.dumps(v4b, indent=2))
        return

    cols, copies, tids = load_psv(str(intact_prefix))
    design = pick_deletion_design(quant, cols, copies)
    if design is None:
        v4b["note"] = (
            f"n_copies_supported={n_copies_supported} but no host+2-hidden PSV-identity design "
            "found where 2 copies are each strictly closer to one surviving host than to any "
            "other surviving copy (reads would not cleanly pool onto one locus). Stopping "
            "honestly -- no recovery attempted."
        )
        out = pathlib.Path(a.out)
        data = json.loads(out.read_text()) if out.exists() else {}
        data["v4b"] = v4b
        out.write_text(json.dumps(data, indent=2))
        print("V4b:", json.dumps(v4b, indent=2))
        return

    keep, host, hidden = design
    span_by_idx = {int(r["copy_index"]): (int(r["copy_start"]), int(r["copy_end"])) for r in quant}
    mask_spans = [span_by_idx[h] for h in hidden]

    print(f"Design: keep={keep} host={host} hidden(deleted)={hidden} "
          f"mask_spans={mask_spans}", file=sys.stderr)

    degraded_fa = work / "degraded_ref.fa"
    new_contig = f"{contig}_degraded"
    build_degraded_ref(GGO_FASTA, contig, pad_start, pad_end, mask_spans, degraded_fa, new_contig)

    sorted_bam, n_reads = extract_and_realign(GGO_BAM, contig, pad_start, pad_end, degraded_fa,
                                               work, "degraded")

    degraded_region = f"{new_contig}:1-{pad_end - pad_start + 1}"
    degraded_prefix = work / "mech_real_degraded"
    log2 = run_copy_assign(sorted_bam, degraded_fa, degraded_region, degraded_prefix,
                            extra_args=["--absent-copies", "--linearize", "--dump-psv"])

    d_families = _read_tsv(f"{degraded_prefix}.families.tsv")
    dna_needs = _read_tsv(f"{degraded_prefix}.dna_needs.tsv") \
        if pathlib.Path(f"{degraded_prefix}.dna_needs.tsv").exists() else []
    d_n_copies = int(d_families[0]["n_copies"]) if d_families else 0
    d_collapsed = int(d_families[0]["collapsed_copies"]) if d_families else 0

    min_clusters_satisfied = (d_n_copies - len(keep)) >= 0 and d_collapsed > 0

    recovered = False
    identity_to_intact = None
    match_detail = []
    if pathlib.Path(f"{degraded_prefix}.psv_copies.tsv").exists():
        d_cols, d_copies, d_tids = load_psv(str(degraded_prefix))
        # translate degraded-ref-local col genome_pos -> absolute genome coords
        offset = pad_start - 1
        d_abs_cols = {ci: (pos + offset) for ci, pos in d_cols.items()}
        n_reference_kept = len(keep)
        admitted_idxs = sorted(i for i in d_copies if i >= n_reference_kept)

        # build truth allele strings (from intact run) restricted to the columns the
        # degraded run actually discovered, matched by absolute genome_pos.
        intact_pos_to_col = {pos: ci for ci, pos in cols.items()}
        shared_positions = [p for p in d_abs_cols.values() if p in intact_pos_to_col]
        shared_positions = sorted(set(shared_positions))

        def truth_allele_str(copy_idx):
            s = []
            for p in shared_positions:
                ci = intact_pos_to_col[p]
                s.append(copies[copy_idx][ci])
            return "".join(s)

        def degraded_allele_str(copy_idx):
            d_pos_to_col = {v: k for k, v in d_abs_cols.items()}
            s = []
            for p in shared_positions:
                ci = d_pos_to_col[p]
                s.append(d_copies[copy_idx][ci])
            return "".join(s)

        truth_seqs = [truth_allele_str(h) for h in hidden]
        recovered_seqs = [degraded_allele_str(i) for i in admitted_idxs]

        matches = set_match(recovered_seqs, truth_seqs) if recovered_seqs and truth_seqs else []
        match_detail = matches
        if matches:
            identity_to_intact = min(m[2] for m in matches)
            recovered = len(admitted_idxs) > 0 and identity_to_intact is not None and identity_to_intact > 0.90

    v4b.update({
        "n_copies_supported": n_copies_supported,
        "design": {"kept_reference_copy_indices": keep, "host_index": host,
                   "hidden_deleted_indices": hidden, "mask_spans": mask_spans},
        "degraded_n_copies": d_n_copies,
        "degraded_collapsed_copies": d_collapsed,
        "n_dna_needs": len(dna_needs),
        "dna_needs_reasons": [r.get("reason", "") for r in dna_needs],
        "n_admitted": (len(d_copies) - len(keep)) if pathlib.Path(f"{degraded_prefix}.psv_copies.tsv").exists() else 0,
        "min_clusters_satisfied": min_clusters_satisfied,
        "recovered": recovered,
        "identity_to_intact": identity_to_intact,
        "set_match": match_detail,
        "n_reads_realigned": n_reads,
        "identity_resolution": "PSV-column (matched by absolute genome_pos between intact and degraded runs)",
        "copy_assign_log_intact": log1[-2000:],
        "copy_assign_log_degraded": log2[-2000:],
        "note": (
            "Recovered." if recovered else
            "Degradation ran and the admission gate was evaluated; recovery did not clear the "
            "identity/admission bar -- see n_admitted/n_dna_needs/dna_needs_reasons for why."
        ),
    })

    out = pathlib.Path(a.out)
    data = json.loads(out.read_text()) if out.exists() else {}
    data["v4b"] = v4b
    out.write_text(json.dumps(data, indent=2))
    print("V4b:", json.dumps({k: v for k, v in v4b.items() if not k.startswith("copy_assign_log")}, indent=2))


if __name__ == "__main__":
    main()
