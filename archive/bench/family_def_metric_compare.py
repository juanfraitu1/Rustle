#!/usr/bin/env python3
"""family_def_metric_compare.py — settle the ~B-metric question for the family definition.

The NOTE defines ~B as identity(>=iota=0.80) AND reciprocal-coverage(>=tau=0.30) — TWO
parameters. The production Rust gate (family_graph.rs::contiguous_core_coverage) uses ONE:
the longest contiguous aligned block / shorter length (T=0.13, Compara-validated 12/12).

For every within-family copy pair (validated families, cached exon-union models — no BAM),
compute all three from a single minimap2 asm20 --eqx alignment:
  identity   = matched / (matched + mismatched)
  recip_cov  = min(aligned_q/len_q, aligned_t/len_t)
  core_cov   = longest run of consecutive aligned cols (=/X, no indel) / min(len)
and answer: (a) is iota inert given tau (do all tau-passing pairs already have id>=0.80)?
(b) do tau (>=0.30) and core (>=0.13) keep the SAME pairs? (c) where does core_cov land for
real copies vs the bridge range (<=0.055) from core_gate_atscale.md?
Run: /home/juanfra/miniforge3/bin/python bench/family_def_metric_compare.py
"""
import collections
import json
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
SEQS = "/home/juanfra/winloci_scratch/vg_reinforce_copies.fa"
FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"
TMP = "/home/juanfra/winloci_scratch/metriccmp"
OVERLAP_MERGE = 0.5
IOTA, TAU, TCORE = 0.80, 0.30, 0.13
MAX_COPIES = 8   # skip isoform-noise mega-families for the pairwise check


def load_models():
    models, name, buf = {}, None, []
    for line in open(SEQS):
        if line.startswith(">"):
            if name is not None:
                models[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name is not None:
        models[name] = "".join(buf)
    return models


def dedup(loci):
    by = collections.defaultdict(list)
    for lid, c, s, e in loci:
        by[c].append((s, e, lid))
    regs = []
    for c, it in by.items():
        it.sort()
        cur, ids = None, []
        for s, e, lid in it:
            if cur is not None and s < cur[1] and (min(cur[1], e) - s) > OVERLAP_MERGE * min(cur[1] - cur[0], e - s):
                cur = (cur[0], max(cur[1], e)); ids.append(lid); continue
            if cur is not None:
                regs.append(ids)
            cur, ids = (s, e), [lid]
        if cur is not None:
            regs.append(ids)
    return regs


def align(a, b):
    """minimap2 asm20 --eqx a vs b -> (identity, recip_cov, core_cov) or None."""
    open(f"{TMP}/a.fa", "w").write(f">a\n{a}\n")
    open(f"{TMP}/b.fa", "w").write(f">b\n{b}\n")
    out = subprocess.run(f"minimap2 -c --eqx -x asm20 {TMP}/a.fa {TMP}/b.fa 2>/dev/null",
                         shell=True, capture_output=True, text=True).stdout
    best = None
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qlen, qs, qe, tlen, ts, te = int(f[1]), int(f[2]), int(f[3]), int(f[6]), int(f[7]), int(f[8])
        cg = next((t[5:] for t in f[12:] if t.startswith("cg:Z:")), None)
        if cg is None:
            continue
        m = x = 0
        run = longest = 0
        import re
        for n, op in re.findall(r"(\d+)([=XID])", cg):
            n = int(n)
            if op == "=":
                m += n; run += n
            elif op == "X":
                x += n; run += n
            else:  # I or D -> break the contiguous aligned block
                longest = max(longest, run); run = 0
        longest = max(longest, run)
        ident = m / (m + x) if (m + x) else 0.0
        recip = min((qe - qs) / max(qlen, 1), (te - ts) / max(tlen, 1))
        core = longest / max(min(qlen, tlen), 1)
        score = m
        if best is None or score > best[3]:
            best = (round(ident, 4), round(recip, 4), round(core, 4), score)
    return best[:3] if best else None


def main():
    os.makedirs(TMP, exist_ok=True)
    models = load_models()
    fams = collections.defaultdict(list)
    for line in open(FAM_TSV).read().splitlines()[1:]:
        fi, lid, c, s, e, nr = line.split("\t")
        fams[int(fi)].append((lid, c, int(s), int(e)))

    pairs = []
    for fi in sorted(fams):
        regs = dedup(fams[fi])
        if not (2 <= len(regs) <= MAX_COPIES):
            continue
        # one representative model per copy-region = its longest exon-union
        reps = []
        for ids in regs:
            cand = [(len(models[l]), models[l]) for l in ids if models.get(l)]
            if cand:
                reps.append(max(cand)[1])
        for i in range(len(reps)):
            for j in range(i + 1, len(reps)):
                r = align(reps[i], reps[j])
                if r:
                    pairs.append(dict(family=fi, id=r[0], recip=r[1], core=r[2]))

    # analysis
    tau_pass = [p for p in pairs if p["recip"] >= TAU]
    core_pass = [p for p in pairs if p["core"] >= TCORE]
    both = [p for p in pairs if p["recip"] >= TAU and p["core"] >= TCORE]
    tau_only = [p for p in pairs if p["recip"] >= TAU and p["core"] < TCORE]
    core_only = [p for p in pairs if p["core"] >= TCORE and p["recip"] < TAU]
    iota_binds = [p for p in tau_pass if p["id"] < IOTA]   # tau passes but iota fails -> iota not inert
    cores = sorted(p["core"] for p in pairs)
    recips = sorted(p["recip"] for p in pairs)

    out = dict(
        n_pairs=len(pairs),
        tau_pass=len(tau_pass), core_pass=len(core_pass),
        agree_both=len(both), tau_only_not_core=len(tau_only), core_only_not_tau=len(core_only),
        iota_binds_on_tau_pass=len(iota_binds),
        min_id_among_tau_pass=round(min((p["id"] for p in tau_pass), default=0), 4),
        core_min=round(cores[0], 4) if cores else None,
        core_p10=round(cores[len(cores)//10], 4) if cores else None,
        core_median=round(cores[len(cores)//2], 4) if cores else None,
        recip_min=round(recips[0], 4) if recips else None,
        bridge_ref="core_gate_atscale.md: bridges core<=0.055, copies core 0.174-0.608",
    )
    json.dump(dict(summary=out, pairs=pairs),
              open(os.path.join(HERE, "family_def_metric_compare.json"), "w"), indent=2)
    print("=== ~B metric comparison over validated-family copy pairs (no BAM) ===")
    for k, v in out.items():
        print(f"  {k:28}: {v}")
    if tau_only:
        print("\n  tau-passes-but-core-fails (would change between metrics):")
        for p in tau_only[:8]:
            print(f"    fam{p['family']} id={p['id']} recip={p['recip']} core={p['core']}")
    print("\n[+] wrote family_def_metric_compare.json")


if __name__ == "__main__":
    main()
