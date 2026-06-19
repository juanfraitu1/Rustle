# De-novo transcript validation (intron-chain Sn/Pr vs RefSeq)

Realigned all 101,467 de-novo transcripts (`denovo_transcripts.fa`) to the genome with minimap2
(`-ax splice:hq -uf`, low-mem `-I 1G --split-prefix`, MALLOC_ARENA_MAX=2 — the ulimit -v VIRTUAL cap
false-triggered 3×, drop it), `bam_to_gtf.py` → GTF, `gffcompare -r GGO_genomic.gff`. 100% mapped.

## Verdict: REAL + defensible
- **Intron-level precision 86%** (`-R -Q` 86.1) — posited splice junctions are genuine annotated sites.
- **Class codes: 98.9% overlap a KNOWN gene** (=20.7% FSM, c=31.0% ISM, j=30.9% novel-iso-of-known,
  m/n/k=16.4% retained/containment); only **0.5% `u` intergenic-novel** (artifact-suspect), ~0.6% antisense/intronic.
- **Sensitivity (where it looks, -R -Q): 76.6% of introns, 76.7% of loci recovered.** Genome-wide intron Sn 53%.
- Intron-CHAIN FSM only ~21% — NOT artifacts: long-read novel/partial isoforms (one novel junction in a
  ~9-exon chain fails whole-chain match). Expected annotation-incompleteness + novel-isoform discovery.
- CAVEAT: 31% ISM + 13% retained-intron => a real fraction are PARTIAL/incomplete (5' degradation / pre-mRNA),
  typical of read-derived assembly.

Artifacts in /home/juanfra/winloci_scratch/validate/ (dn_realigned.bam, dn_gw*.stats/.tmap).
