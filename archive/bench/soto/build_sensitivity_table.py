#!/usr/bin/env python3
"""Augment the per-family Soto sensitivity/precision table with this session's findings:
for every missed member, WHY it's missed (from the floor decomposition), which are recoverable
by the unspliced-seeding fix vs the confirmed K=0 identifiability floor, and an honest
post-seeding-fix sensitivity PROJECTION. Emits soto_sensitivity_precision.{tsv,md}."""
import csv, os
B = "/mnt/c/Users/jfris/Desktop/Rustle/bench/soto"

# per-family current sensitivity/precision (real, from the eval)
fam = {}
with open(f"{B}/soto_family_detection.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        fam[r["family_id"]] = r

# per missed member: cause from the floor decomposition
# classes: seeding-fixable (MERGED / unseeded)  vs  K=0 floor (exon-homogenized / young identical / silent)
cause_of = {}   # (family_id, chrom, start, end) -> cause bucket
raw_cause = {}
def bucket(cause):
    if "MERGED" in cause or "unseeded" in cause:
        return "seeding-fixable"
    if "exon-homogenized" in cause or "young identical" in cause:
        return "K0-floor"
    if "silent" in cause:
        return "silent(artifact-corrected)"
    return "other"
with open(f"{B}/soto_floor_decomposition.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        key = (r["family_id"], r["chrom"], r["start"], r["end"])
        cause_of[key] = bucket(r["cause"]); raw_cause[key] = r["cause"]

# per family: tally missed members by bucket
missed = {}  # fid -> {bucket: count, genes:{bucket:[gene]}}
with open(f"{B}/soto_member_detection.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["detected"].strip().upper() == "N":
            fid = r["family_id"]
            key = (fid, r["chrom"], r["start"], r["end"])
            b = cause_of.get(key, "other")
            m = missed.setdefault(fid, {"seeding-fixable":0, "K0-floor":0, "silent(artifact-corrected)":0, "other":0, "genes":{}})
            m[b] += 1
            m["genes"].setdefault(b, []).append(r["gene"])

rows = []
for fid, fr in fam.items():
    n_mem = int(fr["n_members"]); n_det = int(fr["n_detected"])
    mm = missed.get(fid, {})
    sfix = mm.get("seeding-fixable", 0)
    k0 = mm.get("K0-floor", 0)
    sil = mm.get("silent(artifact-corrected)", 0)
    oth = mm.get("other", 0)
    # projection: seeding fix recovers the seeding-fixable + the (corrected) silent seeding case
    proj_det = n_det + sfix + sil
    rows.append({
        "family_id": fid,
        "genes": fr.get("members","").split(";")[0].split(":")[0] if fr.get("members") else "",
        "n_members": n_mem,
        "sens_now": fr["sens_pct"],
        "precision": fr["precision"],
        "missed": n_mem - n_det,
        "seeding_fixable": sfix,
        "K0_floor": k0,
        "silent_corrected": sil,
        "other": oth,
        "proj_sens_postfix": f"{proj_det}/{n_mem} = {100*proj_det/n_mem:.0f}%",
    })

rows.sort(key=lambda r: (r["missed"] == 0, -(r["missed"])))  # incomplete first, most-missed first

# TSV
cols = ["family_id","genes","n_members","sens_now","precision","missed","seeding_fixable","K0_floor","silent_corrected","other","proj_sens_postfix"]
with open(f"{B}/soto_sensitivity_precision.tsv","w") as f:
    f.write("\t".join(cols)+"\n")
    for r in rows:
        f.write("\t".join(str(r[c]) for c in cols)+"\n")

# aggregate
n_fam = len(rows)
full_now = sum(1 for r in rows if r["missed"] == 0)
tot_mem = sum(r["n_members"] for r in rows)
tot_missed = sum(r["missed"] for r in rows)
tot_sfix = sum(r["seeding_fixable"] for r in rows)
tot_k0 = sum(r["K0_floor"] for r in rows)
tot_sil = sum(r["silent_corrected"] for r in rows)
tot_oth = sum(r["other"] for r in rows)
det_now = tot_mem - tot_missed
det_proj = det_now + tot_sfix + tot_sil

md = [f"""# Soto per-family sensitivity / precision — updated {os.popen('date +%F').read().strip()}

## Aggregate

| | members | %  |
|---|---|---|
| **detected now** (current catalog) | {det_now}/{tot_mem} | {100*det_now/tot_mem:.1f}% |
| **projected post-seeding-fix** | {det_proj}/{tot_mem} | {100*det_proj/tot_mem:.1f}% |
| families fully recovered (now) | {full_now}/{n_fam} | {100*full_now/n_fam:.0f}% |

Missed members ({tot_missed}) by cause:
- **seeding-fixable** (unspliced pooled+deleted; the seeding fix recovers): **{tot_sfix}**
- **K=0 identifiability floor** (exon-homogenized / young-identical; non-exon rescue verified = 0 → needs DNA): **{tot_k0}**
- silent (classification artifact; ≥1 is a seeding case): {tot_sil}
- other: {tot_oth}

**Honest reading:** current sensitivity {100*det_now/tot_mem:.1f}%; the seeding fix projects to
~{100*det_proj/tot_mem:.1f}% (recovering the {tot_sfix} seeding-fixable members); the residual
{tot_k0} are the genuine K=0 floor (rigorously confirmed irreducible from RNA — the aligner already
uses UTR/intron/flank and they still tie). The projection is a PROJECTION until the genome-wide
catalog is rebuilt with the seeding-fix binary (`bench/soto/seeding_validation.txt` proves the
per-locus mechanism: OLD 0 reps → NEW 7 reps).

## Per-family table (incomplete families first)

| family | gene | members | sens now | prec | missed | seeding-fix | K=0 | proj post-fix |
|---|---|---|---|---|---|---|---|---|"""]
for r in rows:
    if r["missed"] == 0:
        continue
    md.append(f"| {r['family_id']} | {r['genes']} | {r['n_members']} | {r['sens_now']} | {r['precision']} | {r['missed']} | {r['seeding_fixable']} | {r['K0_floor']} | {r['proj_sens_postfix']} |")
md.append(f"\n*({full_now} families already at 100% sensitivity / 100% precision are omitted from the table above; full data in `soto_sensitivity_precision.tsv`.)*")
open(f"{B}/soto_sensitivity_precision.md","w").write("\n".join(md)+"\n")
print(f"wrote soto_sensitivity_precision.tsv + .md  ({n_fam} families; now {100*det_now/tot_mem:.1f}% -> proj {100*det_proj/tot_mem:.1f}%; K=0 floor {tot_k0})")
