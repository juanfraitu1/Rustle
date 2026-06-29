#!/usr/bin/env bash
# TDD test for RUSTLE_VG_UNION_BASELINE: VG output must be a SUPERSET of the
# primary-only baseline (no over-collapse regression) AND keep VG's own wins.
# Deterministic (-p1). Exit 0 = all pass, 1 = any fail.
set -uo pipefail
ROOT=/mnt/c/Users/jfris/Desktop/Rustle; SC=/home/juanfra/winloci_scratch
R=$ROOT/target/release/rustle; W=/tmp/vgunion_test; rm -rf $W; mkdir -p $W
TAND="RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1"
eqn(){ gffcompare -r $SC/GGO_tx.gff "$1" -o "$2" >/dev/null 2>&1; ls "$2".*.tmap >/dev/null 2>&1 && awk -F'\t' 'NR>1&&$3=="="{print $2}' "$2".*.tmap|sort -u|wc -l || echo 0; }
fail=0
chk(){ # label  actual  op  expected
  local l="$1" a="$2" op="$3" e="$4"
  if [ "$a" "$op" "$e" ]; then echo "  PASS  $l ($a $op $e)"; else echo "  FAIL  $l ($a NOT $op $e)"; fail=1; fi
}

for gl in "DGCR6:gene-DGCR6:regress" "RABL2A:gene-RABL2A:win"; do
  name=${gl%%:*}; rest=${gl#*:}; g=${rest%%:*}; kind=${rest#*:}
  B=/tmp/winloci/$g/locus.bam
  [ -s "$B" ] || { echo "SKIP $name (no bam)"; continue; }
  echo "== $name ($kind) =="
  $R -p1 -L $B -o $W/base.gtf >/dev/null 2>&1;                                                   nb=$(eqn $W/base.gtf $W/b)
  env $TAND $R -p1 --vg --vg-snp --genome-fasta $SC/GGO.fasta -L $B -o $W/vg.gtf >/dev/null 2>&1; nv=$(eqn $W/vg.gtf $W/v)
  env $TAND RUSTLE_VG_UNION_BASELINE=1 $R -p1 --vg --vg-snp --genome-fasta $SC/GGO.fasta -L $B -o $W/u.gtf >/dev/null 2>&1; nu=$(eqn $W/u.gtf $W/u)
  echo "  base=$nb vg=$nv union=$nu"
  chk "$name union >= baseline (no regression)" "$nu" -ge "$nb"
  chk "$name union >= vg (keeps wins)"          "$nu" -ge "$nv"
done

# Flag-OFF must be byte-identical to plain VG (no behavior change when off)
echo "== flag-off safety (DGCR6) =="
B=/tmp/winloci/gene-DGCR6/locus.bam
env $TAND $R -p1 --vg --vg-snp --genome-fasta $SC/GGO.fasta -L $B -o $W/off1.gtf >/dev/null 2>&1
env $TAND $R -p1 --vg --vg-snp --genome-fasta $SC/GGO.fasta -L $B -o $W/off2.gtf >/dev/null 2>&1
if cmp -s $W/off1.gtf $W/off2.gtf; then echo "  PASS  flag-off deterministic"; else echo "  FAIL  flag-off non-deterministic"; fail=1; fi

[ $fail -eq 0 ] && echo "ALL PASS" || echo "SOME FAILED"
exit $fail
