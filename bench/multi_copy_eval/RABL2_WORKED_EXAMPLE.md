# Worked example: graph-level multimapper rescue at the RABL2 paralog family

**One locus that ties the whole arc together:** RABL2 was the flagship VG win-locus, is the
single largest reservoir in the genome-wide multimapper "treasure" (~53 isoforms among the 424
recoverable-only-from-multimappers transcripts), and demonstrates the exact capability — *use long
multi-mapping reads at the graph level to rescue copy-specific transcripts a primary-only
assembler discards.*

## The setup

RABL2A (NC_073235.2) and its paralog RABL2B (NC_086018.1, + further copies on NC_073227.2 /
NC_073231.2) are near-identical recent duplicates. HiFi reads multi-map across the copies:

| copy | primary reads | secondary reads |
|---|---:|---:|
| **RABL2A** (NC_073235.2) | 32 | **183** |
| RABL2B sibling (NC_086018.1) | 54 | 237 |

minimap2 must pick one alignment as *primary*; for ~85% of RABL2A's reads it picks the **sibling**,
leaving RABL2A with only 32 primaries. A primary-only assembler therefore sees almost no
copy-specific evidence at RABL2A.

## The evidence hiding in the discarded pile

Of the **216 reads touching RABL2A**, edit distance shows where they truly originate:

| | reads |
|---|---:|
| **strictly own RABL2A** (NM here < NM at the sibling) | **164** |
| tied (NM equal — non-identifiable) | 0 |
| prefer the sibling | 52 |

**164 reads decisively originate at RABL2A** — they fit RABL2A *better* than any other copy — yet
most are filed as secondary alignments. That is the "treasure": real, copy-specific evidence sitting
in the reads a primary-only tool throws away. (Note the clean separation — 0 ties — is why this copy
is *decidable*; near-identical copies with many ties are the identifiability floor and are not
claimed.)

## Result — three tools on the same RABL2 reads (vs RefSeq)

| tool | RABL2 transcripts recovered (`=`) | predictions |
|---|---:|---:|
| **StringTie** `-L` (discards multimappers) | **4** | 15 |
| rustle baseline (primary-only) | 5 | 16 |
| **rustle VG** (multimapper graph) | **9** | 20 |

The **6 VG-unique recoveries are all the RABL2A copy's isoforms**
(`XM_055379899/901/906/919/920/922`) — precisely the copy StringTie loses by discarding the
secondaries.

## The mechanism (graph level)

1. **Family discovery** links RABL2A and its siblings into one group via shared multi-mapping reads.
2. **Family / variation graph**: homologous exons across copies become shared `ExonClass` nodes; the
   copies' distinct splice structures are retained as per-copy spans.
3. **Read→copy assignment** by edit distance (+ optional `--vg-snp` PSV fingerprints): the 164
   strictly-owning reads are credited to RABL2A rather than the sibling minimap2 favored.
4. **Per-copy assembly** then reconstructs RABL2A's isoforms from its own (mostly secondary) reads.

## Why it's honest, and where it stops

- **Verified real, not fabricated:** the recovery rests on 164 reads with *strictly* lower edit
  distance at RABL2A — the same strict-NM test that, genome-wide, classifies ~half of apparent
  multimapper recoveries as borrowed-secondary phantoms and *rejects* them. RABL2A passes decisively.
- **Bounded by identifiability:** had the copies been so similar that reads tied at both, no test
  (and no assembler) could assign them — that is the measured floor, not a tooling gap.

## The claim this supports

> Long multi-mapping reads carry copy-specific transcript evidence that primary-only assemblers
> (StringTie et al.) discard. Rustle retains them and resolves them **at the graph level** —
> family graph + edit-distance / SNP read-to-copy assignment — to rescue paralog-copy transcripts,
> here recovering **9 vs StringTie's 4** at the RABL2 family, with the recovery verified by strict
> edit-distance ownership (164/216 reads). Genome-wide this rescues ~424 transcripts in
> recently-duplicated families (RABL2 the largest), each backed by decisive own-evidence.
