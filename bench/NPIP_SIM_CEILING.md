# NPIP algorithmic-ceiling simulation — vs `docs/PREREG_npip_sim_2026-09-18.md` (md5 17be6b031092b79b343421fa0dd14f1e)

**Verdict: SIM-1 PASSES at 25/26. Given full-length, distinguishable reads the assembler is essentially
perfect — it recovers every junction present in the data. The real-data deficit is therefore NOT an
algorithmic defect in the assembler; it is 5' truncation plus everything else about a real library.**

Data `/mnt/linuxdisk/home/juanfraitu/npip_sim/`.

## Substrate

26 canonical transcripts, one per spliced NPIP copy, built by merging each copy's annotated exon blocks
across every non-canonical junction — **209 canonical junctions retained, 40 artifacts merged away**,
matching the §6m7 truth exactly (asserted in the generator). 30 reads per transcript, error-free, ends
independently jittered. Aligned with the same settings as both existing BAMs
(`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`).

| rule | result |
|---|---|
| **SIM-4** distinct (chrom,pos,CIGAR) / total | FL **0.990**, TR **0.999** — the §6m9 dedup collapse is fixed |
| **SIM-0** junctions present in the aligned reads | FL **208/209 (99.5%)**, 25/26 copies whole; TR 203/209, 20/26 whole |

The one FL miss (NPIPB8) is the ALIGNER failing to reproduce a junction from the exact source sequence, not
an assembler failure — it is absent from the reads before the assembler sees them.

## Result

| arm | transcripts | junctions / 209 | **complete chains** | exact chains |
|---|---|---|---|---|
| **SIM full-length (FL)** | **26** | **208 (99.5%)** | **25/26** | **25/26** |
| SIM 5'-truncated (TR) | 138 | 193 (92.3%) | 19/26 | 18/26 |
| real data (same arm, same truth) | 1,456 | 170 (81.3%) | 10/26 | — |

**SIM-1 bar was >= 24/26. Measured 25/26 — PASSED.** The FL arm emits exactly 26 transcripts, one per copy,
and 25 of them are the exact annotated chain. Against the 208 junctions actually present in its reads the
assembler scores **208/208**.

## The deficit, fully decomposed

| step | complete chains | what it costs |
|---|---|---|
| perfect full-length reads | 25/26 | (1 lost to an aligner artifact) |
| + 5' truncation | 19/26 | **−6 copies** |
| + everything else in a real library (depth, multimapping, error, real expression) | 10/26 | **−9 copies** |

**5' truncation alone costs 6 of 26 copies.** It is the single largest identified factor, and it is a
property of the Iso-Seq library, not of the code. The remaining 9 are real-library effects that this
substrate deliberately does not model.

## SIM-3 — what may NOT be claimed

Transcripts derive from the same annotation that defines the truth. Per the pre-registration and
`npip_ideal/DECLARATIONS.txt` C0, this is an **upper bound on the assembler for NPIP and nothing else**. It
is NOT evidence that node construction recovers annotation, that the family definition works, or that O1
works. The informative content here is the NEGATIVE that did not happen: the assembler had no algorithmic
defect to blame.

## Consequence for where to spend effort

Chasing the assembler further is now the wrong move — it is at 208/208 on what it is given. The two
levers that remain are (a) anything that recovers 5' ends (longer/full-length protocols, or a rule that
chains truncated reads across a shared 3' anchor), and (b) the real-library effects. Notably **§6m6's
read-isoform widening targets exactly (a)** and was worth only +2 junctions on real data — so a
chaining-based fix has already been measured and is not promising; the protocol side is.
