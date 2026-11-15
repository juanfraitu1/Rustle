# PREREG row 1666: which reads judge a single-exon model as internally primed (2026-11-13)

Assembly-layer row (base output, not the read-assignment layer). Step 0 only so far; no code change.

## Observation
A119b default `PREFIX.gtf` has no model at single-exon genes with hundreds of clean unspliced reads (PPP1R2B 1,107 reads,
1 primed; FAM47A 396, 13 primed). With `RUSTLE_POLISH_MONO_PRIMED=0` (whole-contig runs, workspace `r1666_mono_root`)
both return, as 5.7 kb / 8.4 kb models on 2.2 / 2.6 kb genes. `mono_primed_count` counts a read only when it overlaps
>= half of the MODEL (the rule registered in `bench/stgaps/PREREG_mono_primed_2026-10-26.md`). Hypothesis: a few long
primed reads stretch the model and are then the only reads counted. Not yet shown directly.
TTC30B / CTAGE15 are a different cause (45 kb unspliced blob removed by the shadow rule): not this row.

## Step 0 (offline, frozen before running; `bench/stgaps/mono_primed_overlap.py`)
Contigs chr2, chr5, chr7, chrX. Dropped models = single-exon transcripts of the filter-off run absent from the default
run. Per model, (counted, primed) over unspliced primary records (`-F 2820`, no N) under:
- `model` (current): overlap >= half the model;
- `shorter`: overlap >= half the shorter of read and model;
- `any`: any overlap.
Validity gate: `model` must reproduce >= 0.98 of the drop decisions (primed/counted >= 0.8, counted > 0).
Report: drops decided by <= 3 reads; flips to kept per candidate; flips by class: one of the 291 no-model members /
another annotated gene (model covers >= half the gene's exon bases, CAT) / nothing annotated.
Stop rule: go on to a switch-guarded Rust change only if, for one candidate, flips regain >= 10 annotated genes on the
four contigs AND the annotated share among flips is >= 2x the annotated share among models that stay dropped.
Otherwise the row closes as a diagnosis. Gorilla KB3781 drops only 118 single-exon models in total, so a transfer test
there is near-vacuous; HG002 would be the second sample.

## Step 0 result: GO (job 4468807)
3,449 dropped models on the four contigs; the recount reproduces all 3,449 decisions (gate passed). 1,502 drops were
decided by <= 3 reads. `shorter` and `any` flip the same 304 models to kept: 44 annotated (2 of them no-model family
members), 260 nothing annotated. Annotated share 0.145 among flips vs 0.042 (133 / 3,145) among models that stay
dropped = 3.4x. Both conditions met.

## Step 1 (frozen before the build)
Switch `RUSTLE_MONO_PRIMED_OVERLAP=shorter` (unset = today's rule, byte-identical): a read counts when it overlaps
>= half of the shorter of read and model. Own target directory. Checks, A119b genome-wide `assemble` + `families`,
against a same-binary switch-off run:
1. switch off reproduces the default GTF byte for byte (validity);
2. fixed truth set = the baseline run's funnel: no-model members that gain a model >= 10 (ceiling unknown; 2 on four
   contigs), reported separately from members placed (the restored models are stretched, so placed may lag);
3. families hits (true pairs, pairscore hit dump) not below baseline, false pairs + <= 10;
4. cost: restored single-exon models genome-wide and the share that hit nothing annotated, reported;
   multi-exon transcripts identical.
PASS (bars 1-3) = proposed to the user as an opt-in or default; the choice is the user's. HG002 is the second sample,
run only after the user sees A119b.

## Step 1 result: FAIL as registered (bar 2), job 4468867
Bar 1 passed (switch off = the default GTF byte for byte); multi-exon transcripts identical; families unchanged
(1,768 hits both, 192 false both, 1 true pair gained and 1 lost; placed 482 -> 484).
Bar 2 failed: 2 of the 291 no-model members gain a model (need >= 10).
The rule also removes far more than it restores: single-exon models 1,266 restored, 7,682 newly dropped (filter drops
13,884 -> 20,300). Step 0 missed that direction: under today's rule 10,332 single-exon models have NO counted read and
are therefore never judged; `shorter` gives 9,367 of them reads, and most of those are primed. So the half-of-model
rule mostly acts as an exemption for stretched models, not as a trap for good ones. The mechanism at PPP1R2B / FAM47A
is real but rare. Switch removed from the source; `bench/stgaps/mono_primed_overlap.py` kept.
