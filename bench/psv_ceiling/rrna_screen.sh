#!/bin/bash
# Annotation-free rRNA screen of a copy catalog: megablast of copies.fa against the MATURE rRNA sequences (human 18S/5.8S/28S cut
# from U13369.1 at 3657-5527 / 6623-6779 / 7935-12969, and 5S NR_023363.1). Mature rRNAs are a universal, species-independent
# sequence class (99 % human-gorilla), so this uses no annotation of the target genome. NEVER blast against the whole 45S unit: its
# intergenic spacer carries Alu-like repeats that hit >100 unrelated families at ~83 % identity.
#   bench/psv_ceiling/rrna_screen.sh copies.fa > hits.tsv      (qseqid sseqid pident length qlen evalue; E <= 1e-20)
# On the gorilla KB3781 catalog (667 families) exactly three families hit: SM5 (49/54 copies, 18S/5.8S/28S), SM7 (51/51, 5S),
# SM577 (2/2, 5.8S) — the same three the RefSeq rRNA biotype flags — and no other.
B=/home/juanfra/miniforge3/envs/blast/bin; D=$(dirname "$0"); T=$(mktemp -d)
$B/makeblastdb -in $D/rrna_mature.fa -dbtype nucl -out $T/db >/dev/null
$B/blastn -task megablast -query "$1" -db $T/db -outfmt "6 qseqid sseqid pident length qlen evalue" -evalue 1e-20 -num_threads 2 -max_target_seqs 5 2>/dev/null
rm -rf $T
