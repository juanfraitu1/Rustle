#!/usr/bin/env bash
# Task 7 validation: the rigorous beat-StringTie cut for the GATED read-coherence layer.
# For each chrom, gffcompare {st, rcg(gated), rc(ungated)} vs RefSeq (ref_$C.gff3) and collect
# the FSM(=) ref transcript ids from each .refmap. Then genome-wide:
#   - rcg_FSM \ st_FSM   = exact annotated tx the GATED layer finds that StringTie misses (target ~580)
#   - rc_FSM  \ st_FSM   = same for the UNGATED layer (the original 580 baseline)
#   - retention          = how much of the ungated win the gate keeps
#   - total tx counts rcg vs rc  = noise dropped (raw ~10,144 -> ~4-5k)
set -uo pipefail
OUT=/tmp/gw
CMP=/tmp/gw/rcgcmp
GFFCOMPARE=/home/juanfra/miniforge3/bin/gffcompare
mkdir -p "$CMP"
CONTIGS=$(ls "$OUT"/rcg_*.gtf 2>/dev/null | sed -E 's#.*/rcg_(.*)\.gtf#\1#')
[[ -n "$CONTIGS" ]] || { echo "no rcg_*.gtf found (run rcg_gen.sh first)"; exit 1; }
for C in $CONTIGS; do
  REF="$OUT/ref_$C.gff3"
  [[ -s "$REF" ]] || { echo "[$C] no ref, skip"; continue; }
  for ARM in st rcg rc gd; do
    Q="$OUT/${ARM}_$C.gtf"
    [[ -s "$Q" ]] || continue
    [[ -s "$CMP/${ARM}_$C".*.refmap ]] 2>/dev/null && continue
    "$GFFCOMPARE" -r "$REF" -o "$CMP/${ARM}_$C" "$Q" >/dev/null 2>&1
  done
done
echo "===== GATED READ-COHERENCE — RIGOROUS BEAT-STRINGTIE CUT ====="
python3 - "$OUT" "$CMP" <<'PY'
import sys,glob,os,re
OUT,CMP=sys.argv[1],sys.argv[2]
def fsm_set(prefix):
    """Union over chroms of FSM(=) ref ids (keyed chrom:ref_id) from <prefix>_<C>.*.refmap."""
    s=set()
    for f in glob.glob(f"{CMP}/{prefix}_*.refmap"):
        b=os.path.basename(f)
        m=re.match(rf"{prefix}_(NC_\d+\.\d+)\.",b)
        C=m.group(1) if m else b
        with open(f) as fh:
            fh.readline()
            for line in fh:
                p=line.rstrip("\n").split("\t")
                if len(p)<3: continue
                if p[2]=="=": s.add(f"{C}:{p[1]}")
    return s
def tx_count(prefix):
    n=0
    for f in glob.glob(f"{OUT}/{prefix}_*.gtf"):
        with open(f) as fh:
            for l in fh:
                if "\ttranscript\t" in l: n+=1
    return n
def tagged_count():
    n=0
    for f in glob.glob(f"{OUT}/rcg_*.gtf"):
        with open(f) as fh:
            for l in fh:
                if 'source "read_coherence"' in l: n+=1
    return n
st,rcg,rc,gd=fsm_set("st"),fsm_set("rcg"),fsm_set("rc"),fsm_set("gd")
rcg_only=rcg-st; rc_only=rc-st
kept=rcg_only & rc_only
ntx_rcg, ntx_rc = tx_count("rcg"), tx_count("rc")
print(f"chroms (rcg): {len(glob.glob(f'{OUT}/rcg_*.gtf'))}")
print(f"FSM(=) ref tx matched:  StringTie={len(st)}  flow(gd)={len(gd)}  gated-rc={len(rcg)}  ungated-rc={len(rc)}")
print()
print("--- BEAT-STRINGTIE (the headline) ---")
print(f"GATED   cut (rcg_FSM \\ st_FSM): {len(rcg_only)}   <-- exact annotated tx StringTie misses")
print(f"UNGATED cut (rc_FSM  \\ st_FSM): {len(rc_only)}   (the ~580 baseline)")
if rc_only:
    print(f"retention of ungated win by the gate: {len(kept)}/{len(rc_only)} = {100*len(kept)//len(rc_only)}%")
    lost=rc_only-rcg_only
    print(f"  (gate dropped {len(lost)} of the ungated FSM wins; gained {len(rcg_only-rc_only)} new)")
print()
print("--- CRITICAL-1 COST: flow replacement (rcg drops de-novo flow in favor of read-chain) ---")
gd_only_st = gd - st                      # what the pure flow baseline beats ST by
flow_lost  = gd - rcg                     # flow FSM the replacement gives up
flow_lost_vs_st = (gd - st) - (rcg - st)  # of the flow's beat-ST wins, how many rcg loses
print(f"flow(gd) beat-ST cut (gd_FSM \\ st_FSM): {len(gd_only_st)}")
print(f"flow FSM LOST by replacement (gd_FSM \\ rcg_FSM): {len(flow_lost)}")
print(f"  of which were flow's beat-ST wins now lost: {len(flow_lost_vs_st)}")
union_only = (rcg|gd) - st
print(f"TRULY-ADDITIVE ceiling (rcg ∪ flow) \\ st_FSM: {len(union_only)}   (what option B would deliver)")
print()
print(f"total tx:  gated-rc={ntx_rcg}   ungated-rc={ntx_rc}   (noise dropped by gate: {ntx_rc-ntx_rcg})")
print(f"read_coherence-tagged tx in gated output: {tagged_count()}")
PY
echo "RCG_VALIDATE_DONE"
