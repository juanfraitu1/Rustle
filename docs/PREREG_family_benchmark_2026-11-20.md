# PREREG: one family definition on every assembler's output (row 1855; filed as 1847 until that number turned out taken by `r1847_sirv_rules`, 2026-11-20)

User (2026-11-20): show the advisor that Rustle assembles multi-copy gene families, with an easy post hoc family definition
applied to each tool's assembly. Bar proposed by me and accepted by the user ("ok lets do that"): Rustle within 0.02 of the
best tool on both pair precision and pair recall, on A119b and OR6737. Written before any family run on a tool's GTF.

## What is compared
The shipped `families` stage (driver `tools/rustle_pipeline.sh families`, binaries of commit 3a1005d1): loci of a GTF
(`gene_id` groups), all-vs-all alignment of the locus sequences, exon-share edge rule, SD damping with the genome's SEDEF BED,
MCL. Same command for every arm; `RUSTLE_FAMREPS=0` (the extra representatives are a Rustle-side addition).

| Arm | Input GTF |
|---|---|
| `rustle` | `r1817_transfer_rerun/{a1,or}/ours/ours.gtf` with its `ours.families.gtf` (the pipeline's default family input, defaults of 2026-11-19) |
| `rustle_final` | `ours.final.gtf` of the same runs |
| `rustle_default` | as `rustle` with the extra representatives on (the pipeline as shipped); reported, not judged |
| `stringtie`, `flair` | `benchmark_collapse/{A119b,GGO}.{stringtie.gtf.gz,flair.isoforms.gtf.gz}` (de novo) |
| `isoquant` | `r1646_isoquant/{a,g}.dn_pt/x/x.transcript_models.gtf` (de novo, `--polya_trimmed all`) |
| `bambu` | `r1641_bambu/{a,g}.dn.gtf` (de novo, NDR 1) |

Tool arms run with `RUSTLE_BRIDGE_REGROUP=off` (their GTF is the family input as it is). Known asymmetry: the locus
representative is the transcript with most `reads`; other tools' GTFs carry no `reads` attribute, so theirs is the longest.

## Score
`bench/denovo_eval/pairscore.py`, default level, same-symbol pairs set aside (the standing scorer). A119b: human truth.
OR6737: protein-hit truth (`DENOVO_GGO_ENS_MAP=r1611/ens_map.tsv`), symbol truth beside it.
- Precision: as the scorer reports it.
- Recall: on a FIXED truth set per sample = truth pairs among genes with a transcript in ANY of the five assemblies
  (`--gtf` = the five input GTFs concatenated). Recall on each arm's own expressed genes is reported beside it.
- Beside, no bar: predicted pairs, hit pairs, genes in at least one hit pair (members), families and copies per arm.

## Rule
Judged arm = `rustle`. Best tool per axis among `stringtie`, `flair`, `isoquant`, `bambu`; on the precision axis only tools
whose fixed-set recall is at least half of Rustle's count.
- PASS: on both samples, precision >= best tool - 0.02 AND fixed-set recall >= best tool - 0.02.
- FAIL: otherwise, said per sample and axis. No change of the definition or its parameters after seeing the scores.

## Predictions (mine)
Parity with StringTie on both axes; FLAIR ahead on recall on A119b (it emits several times more transcripts), so the recall
bar may fail there. Rows 1670 (clean members at parity) and 1613 are the basis.

## Limits
Both samples were development samples for the family definition (edge rule, damping), tuned on Rustle's loci: the
definition may fit Rustle's locus shapes better than a tool's. Tool GTFs are from October / November and were made with each
tool's defaults; Rustle's are at the defaults of 2026-11-19.

## Result (2026-11-21): FAIL, on OR6737 recall alone
Register row 1855 (workspace `r1847_family_benchmark`; the number 1847 was also taken by `r1847_sirv_rules`, so the row is
1855 and `r1855_family_benchmark/README` points here). Jobs 4491151 (FLAIR A119b) and the earlier ones of `run.sbatch`.

Fixed truth set, precision / recall (true pairs found):

| arm | A119b (1,936 pairs) | OR6737 protein-hit truth (487 pairs) | OR6737 symbol truth (102 pairs) |
|---|---|---|---|
| `rustle` (judged) | 0.846 / 0.524 (1,015) | 0.260 / 0.320 (156) | 0.529 / 0.176 |
| `rustle_default` | 0.822 / 0.545 (1,056) | 0.263 / 0.335 (163) | |
| `rustle_final` | 0.842 / 0.500 (968) | 0.251 / 0.310 (151) | |
| `stringtie` | 0.846 / 0.427 (827) | 0.227 / 0.353 (172) | 0.586 / 0.167 |
| `flair` | 0.188 / 0.472 (914) | 0.220 / 0.279 (136) | 0.484 / 0.147 |
| `isoquant` | 0.836 / 0.333 (644) | 0.229 / 0.234 (114) | 0.556 / 0.147 |
| `bambu` | 0.886 / 0.245 (475) | 0.217 / 0.232 (113) | 0.517 / 0.147 |

- A119b: PASS on both axes. Precision equals the best eligible tool (StringTie 0.846; bambu's 475 hits are below half of
  Rustle's 1,015, so it is not eligible); recall is the highest, 0.052 above FLAIR and 0.097 above StringTie.
- OR6737: precision PASS (highest of the five); recall FAIL, 0.320 against StringTie 0.353 (-0.033, bar -0.02).
- Recall on each arm's own expressed genes, A119b / OR6737: rustle 0.704 / 0.430, stringtie 0.477 / 0.383, flair 0.573 /
  0.386, isoquant 0.495 / 0.322, bambu 0.765 / 0.395.
- Loci, A119b / OR6737: rustle 64,805 / 19,079, stringtie 75,993 / 19,424, flair 514,452 / 47,069, isoquant 108,676 /
  25,733, bambu 27,572 / 16,422. FLAIR on A119b predicts 4,860 pairs, 3,946 false.
- Predictions: parity with StringTie held on precision; on recall Rustle is ahead on A119b and behind on OR6737. FLAIR
  does not lead A119b recall (wrong).

Deviations: FLAIR transcripts with a zero-length exon were dropped so the stage could read the GTF (14 on A119b, 3 on
OR6737); A119b FLAIR needed the 1 TB node (out of memory at 250 GB); the member count (genes in a hit pair) was not
computed, the scorer's dump files overwrote each other.
