#!/bin/bash
# bench/seed_pool/test_run_sh.sh — tests of run.sh's resume, refusal and budget-stop logic with STUB binaries (no assembly, no aligner, no BAM read).
# Written after the independent verification of 2026-10-07 found that `arm` with a mistyped name assembled the default pool and exited 0, that a stale shard
# log made an unrelated failure look like "run again", and that a non-empty file was taken for a finished step.
#   bash bench/seed_pool/test_run_sh.sh        (about 20 seconds; one rlock heavy slot per driver call)
set -u
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
T=$(mktemp -d /mnt/linuxdisk/tmp/test_run_sh.XXXXXX); trap 'rm -rf "$T"' EXIT
BIN=$T/bin; W=$T/work; mkdir -p "$BIN" "$W/mol"
BAM=/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam
printf '#as_table\tbam=%s\trecords=1\tmolecules=1\n' "$(readlink -f "$BAM")" > "$W/mol/mol.tsv"      # a table the driver accepts without building one
cat > "$BIN/copy_assign" <<'EOS'
#!/bin/bash
out=""; while [ $# -gt 0 ]; do case "$1" in --out) out=$2; shift 2;; *) shift;; esac; done
echo x >> "$STUB_DIR/assemble_calls"
printf 'chr16\tstub\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' > "$out.gtf"
sleep 1
printf 'chr16\tstub\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' > "$out.families.gtf"
EOS
cat > "$BIN/mcl_families" <<'EOS'
#!/bin/bash
if [ "${1:-}" = "--help" ]; then echo "writes <out>.copies.tsv"; exit 0; fi
out=""; gtf=""; while [ $# -gt 0 ]; do case "$1" in --out) out=$2; shift 2;; --from-gtf) gtf=$2; shift 2;; *) shift;; esac; done
echo x >> "$STUB_DIR/families_calls"
loci="${gtf%.families.gtf}.fam.loci.fa"; : > "$loci"
if [ -n "${STUB_WRAPPER:-}" ]; then      # act like the shard wrapper: leave a log that names THIS loci FASTA, then fail
  d="$MM2_SHARD_DIR/key/aaa/qkey/t4.B1"; mkdir -p "$d"
  { echo "2026-10-07 00:00:00 start mode=wrap target=$loci (key) query=$loci flags='x' dir=$d deadline_in=500s"
    [ "$STUB_WRAPPER" = progress ] && echo "2026-10-07 00:00:01 shard 0/2 done: 5s, 10 records"
    echo "2026-10-07 00:00:02 budget: shard 1/2 not started (remaining 3s, estimate 9s); 1/2 shards done"; } >> "$d/wrapper.log"
  exit 1
fi
[ -z "${STUB_FAIL_FAM:-}" ] || exit 1
printf 'cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\nMCL0\t1\t1\t0\tNA\tchr16\t1\t100\n' > "$out.clusters.tsv"
printf '##gff-version 3\nchr16\t.\tgene\t1\t100\t.\t+\t.\tID=gene-L;Name=L\n' > "$out.loci.gff3"
printf 'family\tcopy\n' > "$out.copies.tsv"
EOS
chmod +x "$BIN"/*
export STUB_DIR=$T
fail=0
run() { RS_BIN=$BIN RS_WORK=$W RS_CALL_S=100 bash "$REPO/bench/seed_pool/run.sh" "$@" > "$T/out.log" 2>&1; echo $?; }
calls() { wc -l < "$T/$1_calls" 2>/dev/null || echo 0; }
check() { if [ "$2" = "$3" ]; then echo "PASS  $1"; else echo "FAIL  $1: got '$2', wanted '$3'"; fail=1; fi; }
P=$W/human_chr16/P/run

check "a pool arm runs"                      "$(run arm human_chr16 P)" 0
check "  assemble and families ran once"     "$(calls assemble) $(calls families)" "1 1"
check "  the markers exist"                  "$(ls $P.assemble.done $P.families.done 2>/dev/null | wc -l)" 2
check "a finished arm is skipped"            "$(run arm human_chr16 P)" 0
check "  no binary was called again"         "$(calls assemble) $(calls families)" "1 1"

: > $P.families.gtf
check "a truncated product is redone"        "$(run arm human_chr16 P)" 0
check "  the assembly ran again"             "$(calls assemble)" 2

echo "arm=G98 flags=--seed-pool good --seed-as-ratio 0.98 substrate=human_chr16 contig=chr16" > $P.assemble.done
check "a marker of another arm is not trusted" "$(run arm human_chr16 P)" 0
check "  the assembly ran again"             "$(calls assemble)" 3

check "an unknown pool arm is refused"       "$(run arm human_chr16 G97)" 2
check "  nothing ran"                        "$(calls assemble)" 3
check "a Rule-1 name is refused by arm"      "$(run arm human_chr16 G98_R1)" 2
check "a Rule-1 name is refused by r1"       "$(run r1 human_chr16 P_R1)" 2
check "  nothing ran"                        "$(calls assemble)" 3

rm -f $P.assemble.done $P.families.done
check "adopt writes the markers"             "$(run adopt human_chr16)" 0
check "  an adopted arm is skipped"          "$(run arm human_chr16 P)" 0
check "  no binary was called"               "$(calls assemble)" 3

# budget stops: only a call that made progress is 'run again' (75); a stale log of another arm is ignored
mkdir -p $W/mm2_shard/stale/aaa/qstale/t4.B1
echo "2026-10-07 00:00:00 start mode=wrap target=/some/other/arm/run.fam.loci.fa (x)
2026-10-07 00:00:02 budget: shard 1/2 not started (remaining 3s, estimate 9s); 1/2 shards done" > $W/mm2_shard/stale/aaa/qstale/t4.B1/wrapper.log
rm -rf $W/human_chr16/A
check "a failure is not mistaken for a budget stop by a stale log" "$(STUB_FAIL_FAM=1 run arm human_chr16 A)" 1
rm -rf $W/mm2_shard $W/human_chr16/A
check "a budget stop with progress is 'run again'"                 "$(STUB_WRAPPER=progress run arm human_chr16 A)" 75
rm -rf $W/mm2_shard $W/human_chr16/A
check "a budget stop without progress is a failure"                "$(STUB_WRAPPER=none run arm human_chr16 A)" 1
exit $fail
