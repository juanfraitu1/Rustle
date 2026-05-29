#!/usr/bin/env python3
"""Transfrag-level parity diff: compare read-to-transfrag abundance assignment Rustle vs ST.

Foundation infrastructure for the read-to-isoform/transfrag-assignment fix. Joins each tool's
transfrag events by (strand, intron-chain) and reports abundance divergences — pinpoints where
Rustle attributes a different read mass to a chain than StringTie (e.g. chain 16958940-17026851:
Rustle abund 1211 vs ST 1), which is the root of the longcov/cov/over-enumeration gaps.

Usage: python3 bench/transfrag_parity_diff.py [RU_STEP] [ST_STEP]
  RU_STEP default transfrag_define (constructed, pre-flow). ST_STEP default transfrag_pre_depl.
Inputs: /tmp/ru_tf.jsonl /tmp/st_tf.jsonl (capture with PARITY_FILTER_STEPS=transfrag_define,transfrag_pre_depl).
"""
import json, sys, collections, statistics as st

# Default to transfrag_pre_depl on BOTH sides — apples-to-apples (same stage, same convention).
# transfrag_define (Rustle graph-build) inflates the Rustle-only count vs ST's pre_depl stage.
RU_STEP = sys.argv[1] if len(sys.argv) > 1 else "transfrag_pre_depl"
ST_STEP = sys.argv[2] if len(sys.argv) > 2 else "transfrag_pre_depl"

def load(path, step):
    """(strand, intron_chain) -> max abundance (collapse dup transfrags of same chain)."""
    d = collections.defaultdict(float)
    rc = collections.defaultdict(float)
    for line in open(path):
        try: e = json.loads(line)
        except: continue
        if e.get("step") != step: continue
        p = e.get("payload", e)
        ich = p.get("introns", "")
        if not ich: continue
        k = (e.get("strand"), ich.rstrip(","))
        ab = float(p.get("abund", 0))
        d[k] = max(d[k], ab)
        if "read_count" in p: rc[k] = max(rc[k], float(p["read_count"]))
    return d, rc

ru, ru_rc = load("/tmp/ru_tf.jsonl", RU_STEP)
stt, _ = load("/tmp/st_tf.jsonl", ST_STEP)
print(f"Rustle {RU_STEP}: {len(ru)} multi-intron chains | ST {ST_STEP}: {len(stt)} chains")

common = set(ru) & set(stt)
ru_only = set(ru) - set(stt)
st_only = set(stt) - set(ru)
print(f"  in both: {len(common)}  Rustle-only: {len(ru_only)}  ST-only: {len(st_only)}")

ratios = sorted(ru[k] / stt[k] for k in common if stt[k] > 0 and ru[k] > 0)
if ratios:
    print(f"\n=== abundance ru/st on common chains ===")
    print(f"  median={st.median(ratios):.3f} p10={ratios[len(ratios)//10]:.3f} "
          f"p90={ratios[9*len(ratios)//10]:.3f} within2x={100*sum(1 for x in ratios if .5<=x<=2)/len(ratios):.0f}% "
          f"exact={100*sum(1 for x in ratios if abs(x-1)<.01)/len(ratios):.0f}%")
    div = [(ru[k]/stt[k], ru[k], stt[k], ru_rc.get(k,0), k) for k in common if stt[k]>0 and ru[k]/stt[k]>5]
    div.sort(reverse=True)
    print(f"\n=== worst over-attributions (Rustle abund >> ST), {len(div)} chains >5x ===")
    for r, ra, sa, rrc, k in div[:12]:
        print(f"  ru_abund={ra:.1f} st_abund={sa:.1f} ratio={r:.0f} ru_read_count={rrc:.0f} strand={k[0]} chain={k[1][:50]}")
