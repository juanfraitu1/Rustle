# Pre-registration — should guided mode have SPLICED and UNSPLICED variants?

**Written 2026-09-22, §6y8, before either node set is built.** User: *"should we have a spliced and
unspliced guided mode?"*

## What the register already says, and why it does not settle it

Guided mode today = annotated gene **BODIES** (genomic span, introns included), `minimap2 -x asm20`.
A spliced variant = the same genes as **concatenated exons**.

- **r272 argues NO**: *"the extra 91 edges come from jointness (RNA supplying what DNA cannot see)"* was
  REFUTED — a DNA-only window cut to the RNA node's length gives **symmetric difference 0** and recovers
  **91/91**. *"The gain is CONTIGUITY, and DNA has it free."*
- **r279 argues YES**: a true DNA edge's covered bases are only **13.2% exonic**, and **134/440 true edges
  are <10% exonic** — 13.2% is the genome-wide exon fraction, so DNA edges are *"BLIND to transcription"*.
- **r277**: the DNA-vs-RNA gap is **fragmentation + substrate choice (one arbitrary isoform)**, not splicing.

⚠ These are about the older E_r/DNA comparison, not about today's `mcl_families` guided mode, and r272's
length control is exactly what a spliced node set does NOT have. So the question is open on this layer.

## The arms

Same chromosome, same gene set, same downstream (`mcl_families --min-exonic-bp 1
--min-shared-exon-frac 0.60`), same `minimap2 -x asm20 -c --eqx -P -t 4`. **Only the node sequence changes.**

| arm | node sequence |
|---|---|
| **U** unspliced (shipped guided) | gene body, genomic span |
| **S** spliced | the gene's exons concatenated (longest transcript) |

## What decides it

This is not "which scores higher" — a mode is worth having if it is **not redundant**. Reported:
- **edge-set overlap**: \|U ∩ S\|, U-only, S-only, and the **Jaccard**;
- family-level pairwise precision / recall / F against the protein referee and Soto, per arm;
- **exonic fraction of the aligned bases** in each arm's edges — r279's 13.2% is the null; a spliced arm
  is 100% by construction, so the informative number is U's, re-measured on this layer.

Development **chr16**; held out **chr2 / chr8 / chr10**.

| outcome | verdict |
|---|---|
| S-only edges are **≥20%** of S's edges AND S's held-out F is within 0.05 of U's | ⭐ **TWO MODES** — they see different things and both work |
| S-only < 20% but S's F ≥ U's | ⚠ **REPLACE, don't add** — one mode, spliced |
| S-only < 20% **and** S's F below U's | ⛔ **ONE MODE** — spliced is a strict subset that also scores worse |

**Predicted, before looking — ⚠/⛔.** r272 says the substrates are not complementary once length is
controlled, and §6u1/u2 adds that **half of family pairs have a genuinely intronless member**, whose
spliced and unspliced sequences are nearly identical — so S should be largely a subset of U. I expect
S-only well under 20% and S's F at or below U's, because concatenating exons destroys the intronic
homology that r279 shows carries 87% of a true DNA edge's aligned bases.

I will not change the arms, the metrics or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⛔ **ONE MODE. Spliced is a strict SUBSET of unspliced and scores far worse.**

Same 930 genes, same headers, same `minimap2 -x asm20 -c --eqx -P` and `mcl_families` config. Only the
node sequence differs (spliced 3.6 MB vs unspliced 41.9 MB = **8.6%**).

| | unspliced (shipped) | spliced |
|---|---|---|
| PAF records | 32,654 | 2,186 |
| graph nodes / edges / clusters | 143 / **293** / 29 | 20 / **25** / 5 |
| referee: sens · prec · **F** | 0.170 · 0.981 · **0.290** | 0.042 · 1.000 · **0.082** |
| Soto: sens · prec · **F** | 0.521 · 0.925 · **0.667** | 0.085 · 1.000 · **0.156** |

⛔**S-only edges = 0 (0.0% of S).** All 25 spliced edges are already among the unspliced 293; U-only is
**268 (91.5% of U)**; Jaccard **0.085**. Both ⛔ conditions fire: S-only far below 20%, and S's F is well
below U's on both truths.

⭐**This reproduces r272 on a new layer and strengthens it.** r272 found symmetric difference **0** once
node length was controlled and concluded *"the gain is CONTIGUITY, and DNA has it free."* Here, with no
length control at all, spliced is not merely non-complementary — it is a **strict subset**.

⭐**And the spliced edge yield tracks the sequence fraction almost exactly**: spliced sequence is 8.6% of
unspliced, and spliced recovers **25/293 = 8.5%** of the edges. The substrate is not finding different
homology; it is finding proportionally less of the same homology.

## r279 re-measured on this layer, and it is the reason

**23.8%** of an unspliced edge's aligned bases are exonic (r279 measured 13.2% on the older layer, the
genome-wide exon fraction). So **76.2% of what a true edge aligns on is intronic**, and concatenating exons
throws that away — which is exactly the 91.5% of edges that vanish. ⚠This also sharpens r279's warning:
the unspliced arm IS partly "blind to transcription" in the sense that most of its evidence is non-exonic,
but the shipped `--min-exonic-bp 1 --min-shared-exon-frac 0.60` conjunct is what supplies the transcription
requirement — the coverage comes from introns, the *licence* comes from the exon conjunct.

⭐**What spliced is good for: precision 1.000 on both truths.** It is a clean certificate over a tiny
population (sensitivity 0.042/0.085), the same shape as every other high-precision/low-recall signal found
today. It is not a mode.

## Prediction scorecard

Predicted ⚠/⛔ with *"S should be largely a subset of U ... because concatenating exons destroys the
intronic homology that r279 shows carries 87% of a true DNA edge's aligned bases."* **Correct, including
the mechanism** — the measured non-exonic share on this layer is 76.2% rather than 87%, and S turned out to
be not "largely" but *entirely* a subset.

> **Scorer (2026-09-22 port):** every `bench/mode_family_score.py` number above is reproduced byte-for-byte by `target/release/family_score` (same flags; 732/732 parity, scipy's assignment tie-breaking included — r1045/r1046). The Python was retired in §6z2.
