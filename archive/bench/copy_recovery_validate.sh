#!/usr/bin/env bash
# Validate genome-wide copy recovery against the 27 strand-aware genuine targets.
# Per chrom: gffcompare rec_$C (recovery ON) and gd_$C (baseline) vs RefSeq ref_$C.
# Then compare each of the 27 genuine copies' class baseline-vs-recovery, + authenticity.
# Run AFTER bench/copy_recovery_gw.sh and after /tmp/cre_guided/genuine_strand.tsv exists.
set -uo pipefail
OUT=/tmp/gw
GFFCMP=/home/juanfra/miniforge3/bin/gffcompare
GENUINE=/tmp/cre_guided/genuine_strand.tsv
[[ -s "$GENUINE" ]] || { echo "missing $GENUINE (run save_genuine.py)"; exit 1; }

# Per-chrom gffcompare rec + gd vs RefSeq (only chroms that have rec output).
for f in "$OUT"/rec_NC_*.gtf; do
  [[ -s "$f" ]] || continue
  C=$(basename "$f" | sed 's/^rec_//; s/\.gtf$//')
  [[ -s "$OUT/ref_$C.gff3" ]] || continue
  [[ -f "$OUT/gc_rec_$C".*.tmap ]] 2>/dev/null || "$GFFCMP" -r "$OUT/ref_$C.gff3" -o "$OUT/gc_rec_$C" "$f" 2>/dev/null
  [[ -s "$OUT/gd_$C.gtf" ]] && "$GFFCMP" -r "$OUT/ref_$C.gff3" -o "$OUT/gc_gd_v_$C" "$OUT/gd_$C.gtf" 2>/dev/null || true
done

python3 - <<'PY'
import glob, csv, collections, os
OUT="/tmp/gw"
def cls_map(prefix, chrom):
    m={}
    for p in glob.glob(f"{OUT}/{prefix}_{chrom}.*.tmap"):
        with open(p) as fh:
            fh.readline()
            for line in fh:
                c=line.rstrip("\n").split("\t")
                if len(c)>=3 and c[1] not in ("-",""):
                    rid = c[1][4:] if c[1].startswith("rna-") else c[1]   # gffcompare stores rna-XM_...
                    # keep best class per ref (=, c, j priority)
                    pr={"=":0,"c":1,"j":2,"k":3,"m":4,"n":5,"o":6,"e":7,"i":8,"u":9}.get(c[2],10)
                    if rid not in m or pr < m[rid][1]: m[rid]=(c[2],pr)
    return {k:v[0] for k,v in m.items()}
rows=list(csv.DictReader(open("/tmp/cre_guided/genuine_strand.tsv"),delimiter="\t"))
by_chrom=collections.defaultdict(list)
for r in rows: by_chrom[r["chrom"]].append(r)
rec_cache={}; gd_cache={}
def get(prefix,chrom,cache):
    if chrom not in cache: cache[chrom]=cls_map(prefix,chrom)
    return cache[chrom]
hdr=f"{'family':<16}{'tx':<18}{'strand':<7}{'own':>5}  base->rec"
print("="*72); print("27 GENUINE COPIES: baseline class -> recovery class"); print("="*72)
print(hdr); print("-"*72)
imp=0; rec_hit=collections.Counter(); base_hit=collections.Counter()
for r in rows:
    ch=r["chrom"]; tx=r["tx"]
    b=get("gc_gd_v",ch,gd_cache).get(tx,"MISS")
    rc=get("gc_rec",ch,rec_cache).get(tx,"MISS")
    base_hit[b]+=1; rec_hit[rc]+=1
    arrow = "  <-- IMPROVED" if (rc in ("=","c","j") and b in ("MISS","i","u","o","e") and rc!=b) else ("" if rc==b else "  (changed)")
    if "IMPROVED" in arrow: imp+=1
    print(f"{r['family']:<16}{tx:<18}{r['strand']:<7}{int(r['own']):>5}  {b:>4} -> {rc:<4}{arrow}")
print("-"*72)
print(f"baseline class dist: {dict(base_hit)}")
print(f"recovery class dist: {dict(rec_hit)}")
print(f"copies improved by recovery (MISS/u/i->=/c/j): {imp} / {len(rows)}")
print(f"  recovery exact (=): {rec_hit.get('=',0)}   contained (c): {rec_hit.get('c',0)}   near (j): {rec_hit.get('j',0)}")
PY
echo "COPY_RECOVERY_VALIDATE_DONE"
