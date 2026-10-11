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

## Result (2026-10-11): NOT ADOPTED as the default (B2 fails on every held-out arm); kept opt-in

**Dev (NC_073244.2, primary):** 4,169 transcripts; intron 58.8 / 81.2, chain 29.0 / 38.7, transcript 27.2 / 38.7, locus 41.1 / 78.1
(default 54.8 / 83.3, 28.3 / 40.8, 26.5 / 40.7, 39.9 / 77.8; StringTie 58.9 / 83.9, 24.7 / 37.0, 23.2 / 36.8, 40.0 / 71.4).
**Fresh held-out** (sensitivity / precision; OFF = the cached default runs):

| sample, pool | transcripts OFF -> ON | intron | intron chain | transcript | locus |
|---|---|---|---|---|---|
| KB3781 primary | 71,855 -> 74,958 | 53.9 / 81.7 -> 55.9 / 80.2 | 27.0 / 36.1 -> 27.6 / 35.3 | 24.3 / 36.0 -> 24.8 / 35.2 | 28.3 / 76.8 -> 28.7 / 77.9 |
| KB3781 good | 72,622 -> 75,773 | 54.1 / 81.5 -> 56.0 / 80.0 | 27.1 / 35.9 -> 27.7 / 35.1 | 24.4 / 35.8 -> 24.9 / 35.0 | 28.5 / 75.4 -> 28.9 / 76.5 |
| chimp primary | 53,909 -> 59,370 | 47.1 / 85.8 -> 52.2 / 83.2 | 16.0 / 39.4 -> 17.1 / 38.1 | 14.7 / 38.8 -> 15.7 / 37.7 | 26.4 / 63.3 -> 27.6 / 65.2 |
| chimp good | 54,821 -> 60,320 | 47.4 / 85.6 -> 52.4 / 83.0 | 16.1 / 39.1 -> 17.2 / 37.8 | 14.8 / 38.5 -> 15.8 / 37.3 | 26.7 / 62.4 -> 27.8 / 64.3 |
| orangutan primary | 106,467 -> 121,865 | 56.3 / 70.9 -> 63.3 / 65.4 | 22.1 / 19.9 -> 23.8 / 18.4 | 20.1 / 19.0 -> 21.6 / 17.8 | 30.5 / 39.9 -> 32.0 / 43.1 |
| orangutan good | 111,629 -> 127,292 | 56.7 / 69.9 -> 63.6 / 64.4 | 22.3 / 19.3 -> 24.0 / 17.9 | 20.3 / 18.2 -> 21.8 / 17.2 | 30.8 / 36.8 -> 32.3 / 40.0 |

**B1 holds** (intron sensitivity +1.9 to +7.0); **B3 holds** (chain sensitivity +0.6 to +1.7); **B2 fails** (chain precision -0.8 to
-1.5, bar 0.5). Locus sensitivity and precision both rise in all six. Opt-in: `--direct-transfrags-union`, `RUSTLE_DIRECT_TRANSFRAGS=union`.
**Against the tools (benchmark samples, genome-wide, primary pool, descriptive — both samples were seen before):**

| | base | intron | intron chain | transcript | locus |
|---|---|---|---|---|---|
| gorilla default | 45.1 / 68.1 | 59.3 / 82.8 | 26.9 / 35.6 | 24.2 / 35.4 | 30.7 / 71.1 |
| gorilla union | 46.9 / 67.5 | 63.1 / 80.2 | 27.9 / 34.1 | 25.1 / 34.0 | 31.5 / 71.9 |
| gorilla StringTie | 50.3 / 63.1 | 63.9 / 83.1 | 24.4 / 34.6 | 22.0 / 34.3 | 31.2 / 64.5 |
| gorilla FLAIR | 48.5 / 46.0 | 62.1 / 61.5 | 25.5 / 17.6 | 23.1 / 16.0 | 28.8 / 49.2 |
| human default | 48.2 / 17.5 | 60.8 / 64.4 | 21.9 / 18.6 | 20.4 / 16.2 | 34.6 / 22.9 |
| human union | 49.2 / 19.2 | 64.4 / 59.2 | 22.9 / 17.4 | 21.3 / 15.7 | 35.1 / 25.3 |
| human StringTie | 59.6 / 14.5 | 63.5 / 57.8 | 18.1 / 15.0 | 16.9 / 12.8 | 33.8 / 20.3 |
| human FLAIR | 66.4 / 8.8 | 67.7 / 30.2 | 22.4 / 6.3 | 21.1 / 3.6 | 31.9 / 6.3 |

Human: the union is above StringTie on intron, chain, transcript and locus (both measures); StringTie keeps base sensitivity. Gorilla: the
union is above StringTie on chain and transcript sensitivity and locus, about level on intron sensitivity (-0.8), below on intron
precision (-2.9) and chain precision (-0.5), and StringTie keeps base sensitivity. Flow-only transcripts are a small part of the
remaining gap (dev: 9 of 195).
