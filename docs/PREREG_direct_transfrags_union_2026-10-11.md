# Direct transfrags + exact chains (union), 2026-10-11

**Pre-registered before the code exists and before any arm is scored.** Follows `docs/PREREG_direct_transfrags_2026-10-11.md`
(register 1265: intron sensitivity +3.5 held-out, chain precision -1.2 to -1.7). User (2026-10-11): are the remaining StringTie-only
transcripts the flow-only ones, and can the rest be included without flow?

**Dev diagnosis (gorilla OR6737 NC_073244.2 only).** Of 195 multi-exon reference transcripts StringTie matches exactly and direct
transfrags does not: 9 (5%) are carried whole by no read (flow-only); **66 are in the DEFAULT output** — >= 2 reads carry exactly that
chain, but direct transfrags absorbs it into a longer transfrag that contains it (an alternative start or end that looks like a
truncation); 71 are carried whole by one read; 29 by >= 2 reads of which one exact; 16 exact >= 2 but dropped by gate or polish in both
modes; 4 only inside longer reads. The 66 are a cost of the absorption rule, not of the missing flow.

## The rule (`copy_assign --direct-transfrags-union`, default off = byte-identical; exclusive with the other two)

Admit every group the default admits (`>= 2` identical reads) AND every maximal transfrag whose direct-transfer abundance reaches the
floor (`--direct-transfrags`). A group in both carries the transferred abundance; a contained group admitted on its own reads keeps them
(its reads still count toward its containers). Everything else unchanged. No new parameter.

## Test (fresh held-out: none of these samples was used by any assembler test)

Arms OFF (the cached default runs) / ON, pools `primary` and `good`. **Held-out (decides):** gorilla KB3781 fibroblast (gorilla reference
of figure 1), chimpanzee PTR and orangutan PPY (their annotation GTFs, `assembly/<sample>/annotation.gff_to_gtf.gtf`). Dev: OR6737
NC_073244.2 (sanity only, no amendment unless written before any held-out number). Per sample and pool: **B1** intron sensitivity +1.0
point or more; **B2** intron-chain precision down by <= 0.5 point; **B3** intron-chain sensitivity not lower. **Verdict:** B1-B3 in all
three samples for both pools => ADOPT as the default; primary pool only => default for the benchmark arm; else opt-in with the numbers.
Reported: transcripts, transcript and locus levels; on OR6737 and A119b (already seen) the comparison to StringTie / FLAIR is descriptive.
