#!/bin/bash
# run_e2e.sh — end-to-end check of the pipeline's candidates stage on the candidate_copies fixture (spec 2026-10-02 §7,
# ruling R13: O2's scope is AS-tied molecules, so a read that realigns uniquely to a candidate is PLACED there by the
# aligner and never enters the certificate; the candidate family is assigned as a 2-copy family with the gate unchanged).
#
# usage: bash tests/fixtures/candidate_copies/driver/run_e2e.sh --bin DIR --out SCRATCH_DIR [--threads N]
#   --bin  the cargo release directory (candidate_copies, copy_assign), passed to tools/rustle_pipeline.sh --bin
#   --out  a scratch directory (created; earlier products in it are overwritten)
# It runs minimap2 and the driver and takes no lock itself: run it under the machine's heavy lock,
#   bash tools/rlock.sh heavy bash tests/fixtures/candidate_copies/driver/run_e2e.sh --bin <target>/release --out <scratch>
# Needs minimap2, samtools and python3 on PATH (as the driver does).
#
# INPUTS. The fixture has no O1 products (its one locus gives the families stage no family), so the families stage's
# three products are built here, at OUT/fx and identically at OUT/ctl/fx:
#   fx.fam.copies.tsv / .fa   the fixture's copy table: family MCL0 = copy A (DN_chrT_10000_A); copy B, 3% from A, is the
#                             copy deleted from genome.fa. Plus a CONTROL family MCL1 without a candidate: one copy at
#                             chrT:30000-30500 with n_reads 0 (the fixture has reads only at copy A, so no control family
#                             WITH reads exists).
#   fx.fam.copies.regions     mcl_families' layout, the hull +- 5 kb: MCL0 chrT:5000-18100, MCL1 chrT:25000-35500.
# Which reads are copy A's and copy B's: by the generator's names (make_fixture.py: A_00..A_59 and B_00..B_59, "the
# names carry the truth"). The fixture BAM cannot tell them apart by position: all 120 primaries lie on copy A's span.
#
# RUNS. The driver's `candidates` then `assign --candidates` on OUT/fx (the stage and the use of its products are opt-in,
# ruling R14), and a plain `assign` on OUT/ctl/fx, which holds no candidates products (--no-cache, so every run computes).
# Last, a plain `assign` on OUT/fx itself, after (a)-(c) (it replaces fx.assign.*).
#
# ASSERTIONS (exit 1 at the first that fails, naming it; exit 0 when all hold):
#   (pre) the candidates stage flags exactly one candidate, cand_MCL0_0 of MCL0
#   (a)   OUT/fx.aug.bam: the primaries of the 60 copy-B reads lie on cand_MCL0_0 with MAPQ >= 60, and the primaries of
#         the 60 copy-A reads are unchanged from the fixture BAM (flag, contig, position, MAPQ, CIGAR)
#   (b)   MCL0 is assigned as a 2-copy family that includes the candidate (fx.assign.families.tsv n_copies 2;
#         fx.assign.quant.tsv lists cand_MCL0_0 and DN_chrT_10000_A), with the AS-tied gate at its default (the candidate
#         run's log line `AS-TIED GATE (ratio 1.00)`); no assigned row puts a copy-B read on the real copy or a copy-A
#         read on the candidate (vacuous when MCL0 has no row, as here: every molecule is a clear-best mapper)
#   (c)   the split ran (candidate families on the augmented inputs, MCL1 on the originals), each concatenated table has
#         one header, and MCL1's rows are the same in every per-family table of the split and of the plain run on
#         OUT/ctl/fx. MCL1 has no reads: (c) checks the split's selection, regions and concatenation, not read assignment.
#   (d)   opt-in (ruling R14): a plain `assign` on OUT/fx, whose candidates products are present, does not use them — the
#         driver says so in one line, makes no split (no fx.assign_cand.*), and MCL0 is assigned as its 1-copy self
set -euo pipefail
usage() { echo "usage: bash $0 --bin DIR --out SCRATCH_DIR [--threads N]" >&2; exit 2; }
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
fx=$repo/tests/fixtures/candidate_copies
BIN=""; OUT=""; THREADS=2
while [ $# -gt 0 ]; do
  case "$1" in
    --bin) BIN=${2:-}; shift 2 || usage;; --out) OUT=${2:-}; shift 2 || usage;; --threads) THREADS=${2:-}; shift 2 || usage;;
    *) usage;;
  esac
done
[ -n "$BIN" ] && [ -n "$OUT" ] || usage
for b in candidate_copies copy_assign; do
  [ -x "$BIN/$b" ] || { echo "run_e2e: $BIN/$b is missing (cargo build --release)" >&2; exit 2; }
done
for t in minimap2 samtools python3; do
  command -v "$t" > /dev/null || { echo "run_e2e: $t is not on PATH" >&2; exit 2; }
done
mkdir -p "$OUT/ctl"
OUT=$(cd "$OUT" && pwd)
P=$OUT/fx
Q=$OUT/ctl/fx
pass() { echo "PASS $*"; }
fail() { echo "FAIL $*" >&2; exit 1; }
# rows of table $2 whose column named $1 equals $3 (columns located by header name, as the binaries write them)
rows_of() { awk -F'\t' -v col="$1" -v want="$3" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } $c[col] == want' "$2"; }
# column $2 (by name) of the rows of table $1 whose family_id is $3
col_of() { awk -F'\t' -v col="$2" -v fam="$3" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } $c["family_id"] == fam { print $c[col] }' "$1"; }

stage_inputs() {
  local p=$1
  { cat "$fx/copies.tsv"
    printf 'MCL1\t0\tDN_chrT_30000_ctl\tchrT\t30000\t30500\t1\t+\t0\t30000-30500\tNA\tcontrol\tgeneCtl\tNA\t1\t500\tNA\tkept\t30000\t30500\n'
  } > "$p.fam.copies.tsv"
  { cat "$fx/copies.fa"
    printf '>MCL1|0|chrT:30000-30500|+|nexon=1\n'
    samtools faidx "$fx/genome.fa" chrT:30001-30500 | grep -v '^>' | tr -d '\n'
    echo
  } > "$p.fam.copies.fa"
  printf 'MCL0\tchrT:5000-18100\nMCL1\tchrT:25000-35500\n' > "$p.fam.copies.regions"
}
stage_inputs "$P"
stage_inputs "$Q"
minimap2 -x splice -d "$OUT/genome.splice.mmi" "$fx/genome.fa" 2> "$OUT/genome.splice.mmi.log"
drv() {
  bash "$repo/tools/rustle_pipeline.sh" "$@" --bam "$fx/reads.bam" --fasta "$fx/genome.fa" --index "$OUT/genome.splice.mmi" \
    --bin "$BIN" --threads "$THREADS" --no-cache
}
echo "== driver: candidates, assign --candidates (OUT/fx); assign (OUT/ctl/fx)"
drv candidates --out "$P"
drv assign --candidates --out "$P"
drv assign --out "$Q"
echo "== assertions"

# (pre) exactly one flagged candidate
flagged=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } $c["flagged"] == 1 { print $c["family"] ":" $c["candidate"] }' "$P.cand.candidates.tsv" | paste -sd,)
[ "$flagged" = "MCL0:cand_MCL0_0" ] || fail "(pre) flagged candidates are '$flagged', expected MCL0:cand_MCL0_0"
pass "(pre) one flagged candidate: cand_MCL0_0 of MCL0"

# (a) placement by the patch realignment
[ -s "$P.aug.bam" ] || fail "(a) $P.aug.bam was not written"
b_all=$(samtools view -F 2308 "$P.aug.bam" | awk -F'\t' '$1 ~ /^B_/' | wc -l)
b_cand=$(samtools view -F 2308 "$P.aug.bam" | awk -F'\t' '$1 ~ /^B_/ && $3 == "cand_MCL0_0" && $5 >= 60' | wc -l)
[ "$b_all" = 60 ] && [ "$b_cand" = 60 ] || fail "(a) copy-B primaries on cand_MCL0_0 with MAPQ >= 60: $b_cand of $b_all (expected 60 of 60)"
pass "(a) the 60 copy-B reads' primaries lie on cand_MCL0_0 with MAPQ >= 60"
a_primaries() { samtools view -F 2308 "$1" | awk -F'\t' -v OFS='\t' '$1 ~ /^A_/ { print $1, $2, $3, $4, $5, $6 }' | LC_ALL=C sort; }
a_fx=$(a_primaries "$fx/reads.bam" | wc -l)
[ "$a_fx" = 60 ] || fail "(a) the fixture BAM holds $a_fx copy-A primaries, expected 60"
if ! diff <(a_primaries "$fx/reads.bam") <(a_primaries "$P.aug.bam") > "$OUT/a_primaries.diff"; then
  fail "(a) the copy-A primaries changed in the realignment (see $OUT/a_primaries.diff)"
fi
pass "(a) the 60 copy-A reads' primaries are unchanged from the fixture BAM (flag, contig, position, MAPQ, CIGAR; on chrT)"

# (b) the candidate family is a 2-copy family, gate unchanged
n_copies=$(col_of "$P.assign.families.tsv" n_copies MCL0)
[ "$n_copies" = 2 ] || fail "(b) MCL0's n_copies in fx.assign.families.tsv is '$n_copies', expected 2"
tids=$(col_of "$P.assign.quant.tsv" copy_tid MCL0 | LC_ALL=C sort | paste -sd,)
[ "$tids" = "DN_chrT_10000_A,cand_MCL0_0" ] || fail "(b) MCL0's copies in fx.assign.quant.tsv are '$tids', expected DN_chrT_10000_A and cand_MCL0_0"
pass "(b) MCL0 is assigned as a 2-copy family with its candidate (n_copies 2; quant copies $tids)"
grep -q 'AS-TIED GATE (ratio 1.00)' "$P.assign_cand.log" || fail "(b) no 'AS-TIED GATE (ratio 1.00)' line in fx.assign_cand.log: the gate is not at its default"
pass "(b) the AS-tied gate is at its default in the candidate run (AS-TIED GATE (ratio 1.00))"
n_rows=$(rows_of family_id "$P.assign.assignments.tsv" MCL0 | wc -l)
false_moves=$(awk -F'\t' 'FNR == 1 { for (k in c) delete c[k]; for (i = 1; i <= NF; i++) c[$i] = i; next }
  FILENAME == ARGV[1] { if ($c["family_id"] == "MCL0") tid[$c["copy_index"]] = $c["copy_tid"]; next }
  $c["family_id"] == "MCL0" && $c["status"] == "assigned" {
    t = tid[$c["assigned_copy"]]
    if (($c["read_name"] ~ /^B_/ && t != "cand_MCL0_0") || ($c["read_name"] ~ /^A_/ && t == "cand_MCL0_0")) n++
  }
  END { print n + 0 }' "$P.assign.family_join.tsv" "$P.assign.assignments.tsv")
[ "$false_moves" = 0 ] || fail "(b) $false_moves assigned rows put a copy-B read on the real copy or a copy-A read on the candidate"
if [ "$n_rows" = 0 ]; then
  pass "(b) no false move among MCL0's assigned rows (vacuous: MCL0 has 0 rows — every molecule is a clear-best mapper, outside the AS-tied gate)"
else
  pass "(b) no false move among MCL0's $n_rows rows"
fi

# (c) the split: the control family's tables are those of a run without candidates
[ -s "$P.assign_cand.families.tsv" ] && [ -s "$P.assign_rest.families.tsv" ] || fail "(c) the split did not make both runs (fx.assign_cand.* and fx.assign_rest.*)"
cand_fams=$(awk -F'\t' 'NR > 1 { print $1 }' "$P.assign_cand.families.tsv" | paste -sd,)
rest_fams=$(awk -F'\t' 'NR > 1 { print $1 }' "$P.assign_rest.families.tsv" | paste -sd,)
[ "$cand_fams" = MCL0 ] && [ "$rest_fams" = MCL1 ] || fail "(c) the runs hold '$cand_fams' (augmented) and '$rest_fams' (originals), expected MCL0 and MCL1"
pass "(c) the split ran: MCL0 on the augmented inputs, MCL1 on the originals"
for t in assignments families quant family_join famcn_readonly; do
  [ -e "$P.assign.$t.tsv" ] && [ -e "$Q.assign.$t.tsv" ] || fail "(c) fx.assign.$t.tsv is missing from the split or the plain run"
  head_line=$(head -1 "$P.assign.$t.tsv")
  [ "$head_line" = "$(head -1 "$Q.assign.$t.tsv")" ] || fail "(c) fx.assign.$t.tsv: the header differs from the plain run's"
  n_head=$(awk -v h="$head_line" '$0 == h' "$P.assign.$t.tsv" | wc -l)
  [ "$n_head" = 1 ] || fail "(c) fx.assign.$t.tsv holds its header $n_head times"
  if ! diff <(rows_of family_id "$P.assign.$t.tsv" MCL1) <(rows_of family_id "$Q.assign.$t.tsv" MCL1) > "$OUT/c_$t.diff"; then
    fail "(c) fx.assign.$t.tsv: MCL1's rows differ from the plain run's (see $OUT/c_$t.diff)"
  fi
done
fams_split=$(awk -F'\t' 'NR > 1 { print $1 }' "$P.assign.families.tsv" | LC_ALL=C sort | paste -sd,)
[ "$fams_split" = "MCL0,MCL1" ] || fail "(c) fx.assign.families.tsv lists '$fams_split', expected MCL0 and MCL1"
pass "(c) one header per concatenated table; MCL1's rows identical to the plain run in assignments, families, quant, family_join, famcn_readonly (MCL1 has no reads: this checks the split's plumbing, not read assignment)"

# (d) opt-in: without --candidates, assign leaves the candidates products alone (and says so)
drv assign --out "$P" 2> "$OUT/d_assign.stderr"
grep -q "present but unused" "$OUT/d_assign.stderr" || fail "(d) a plain assign on OUT/fx did not say that the candidates products are unused (see $OUT/d_assign.stderr)"
! ls "$P".assign_cand.* > /dev/null 2>&1 || fail "(d) a plain assign on OUT/fx made the candidate run (fx.assign_cand.* exist)"
n_copies=$(col_of "$P.assign.families.tsv" n_copies MCL0)
[ "$n_copies" = 1 ] || fail "(d) a plain assign gives MCL0 n_copies '$n_copies', expected 1 (the candidate is not used)"
pass "(d) a plain assign on OUT/fx does not use the candidates products: one line says so, no split, MCL0 has 1 copy"
echo "run_e2e: all assertions hold"
