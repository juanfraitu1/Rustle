# Gene-conversion EMISSION — from a buried attribute to a first-class isoform

Following NOVEL_COMBINATION_PROOF.md, the one genuinely graph-unique cross-copy signal rustle has
is the **gene-conversion / cross-copy sequence recombinant** (a read whose 5′ exons carry copy A's
SNPs and 3′ exons carry copy B's). Previously this was only *detected* and written as a GTF
*attribute* on the native transcripts (`gtf.rs`). This change *emits* it as its own transcript.

## What it does

Opt-in `RUSTLE_VG_MOSAIC_EMIT=1` (requires the detector, `RUSTLE_VG_MOSAIC_ON=1`). For each
**family-confirmed** conversion event, emit one recombinant transcript:

- **structure** = the host copy's intron chain (a gene conversion is a *sequence* mosaic — no
  novel junction, per the proof), so it shares the host's exons;
- **distinct species** — `source = "gene_conversion"`, coverage = the supporting recombinant
  molecule count, GTF attribute `conversion_isoform "recombinant"` plus the inherited
  `gene_conversion` / `conversion_copies` / `conversion_breakpoint` / `conversion_reads`.

Default-off → byte-identical when the flag is unset.

## Implementation (TDD)

- `vg.rs::build_conversion_transcript(host, event) -> Option<Transcript>` — pure, 6 unit tests
  (`conversion_emit_*`). Builds the recombinant from a CONFIRMED event; returns `None` for
  unconfirmed events (don't-fabricate), a host that is **not a participant copy** (`copy_a`/
  `copy_b` — else it would clone the wrong chain), or a breakpoint **bracket that overlaps no
  exon**. The bracket endpoints are exonic SNP sites but legitimately span the intron between
  two exons, so the check is bracket↔exon **overlap**, not midpoint-in-exon (a regression test
  pins the bracket-spans-intron case).
- `vg.rs::emit_family_recombinants(transcripts, events, fam_id)` — pure, 3 unit tests. Emits
  **one recombinant per distinct confirmed event** in a family (not just the best one), each
  carrying its OWN event in `family_verdict.conversion` so its GTF breakpoint/copies/reads are
  its own. **Deterministic** host selection: prefers `copy_a`, then lowest `(start, end)` —
  independent of rayon cross-bundle ordering.
- `pipeline.rs` — a function-level `family_conversions` map carries ALL confirmed events per
  family from detection to the `RUSTLE_VG_MOSAIC_EMIT` block (the verdict only holds the best
  one). The block iterates families in sorted order, calling `emit_family_recombinants`.
- `gtf.rs` — emits `conversion_isoform "recombinant"` for the tagged transcript.
- `transcript_filter.rs::is_rescue_protected` — `source=="gene_conversion"` is protected from
  the opt-in subset dedup (the recombinant shares the host chain, so it would otherwise be
  dropped as a same-chain subset).

## Adversarial review (2-agent workflow)

Both lenses returned *ship-with-minor-fixes*; all real findings fixed: **(major)** host was
selected by span-containment alone → could clone a non-participant copy's chain (fabrication) and
rode on nondeterministic cross-bundle order → now copy-constrained + deterministic; **(minor)**
intronic-midpoint → bracket-overlap; **(minor)** subset-dedup → protected; **(nit)** dead
`tpm/fpkm` stores → removed. The review's *one-recombinant-per-family* bound was subsequently
**lifted** (see "Multiple events" below). `RUSTLE_VG_MOSAIC_EMIT` without `RUSTLE_VG_MOSAIC_ON`
is a silent no-op (documented).

## Multiple events per family

Originally only the family's single "best" event reached the verdict, so a family with several
real conversions emitted one recombinant. Now `family_conversions` carries **every** confirmed
event to the emission, and `emit_family_recombinants` produces one recombinant per distinct
`(copy_a, copy_b, breakpoint)`, each hosted on its participant copy and carrying its own
attributes. Verified on a synthetic family with two planted breakpoints (exon-2 and exon-3
switches) → **2 recombinants**, `6643-9485`/6 reads and `204080-206643`/7 reads.

## Verification

| check | result |
|---|---|
| unit tests (`conversion_emit_*`, `emit_family_*`) | RED→GREEN, 9/9 |
| full lib suite | 262/0 |
| **two-event** family (badread HiFi, 2 planted breakpoints) | **2** recombinants, each with its own breakpoint/reads |
| single-event regression | **1** recombinant, host = copy0 (=copy_a), `cov 11`, breakpoint 9655 |
| **don't-fabricate** — unconfirmed (0.997) / pure no-recombinant | **0 / 0** emitted |
| content-determinism (4 runs) | identical sorted content + identical recombinant lines |
| default (no flags) | unaffected (4 transcripts, gated) |

*(Note: overall GTF line-order has a pre-existing rayon cross-bundle nondeterminism — base VG
shows it too, with identical sorted content. The emission's content is deterministic.)*

## Honest scope

This is the **sequence channel**, not a novel exon combination — by the proof, a gene conversion
has no novel junction, so the recombinant's *structure* equals its host copy's. The value is
making the cross-copy molecular species a pullable GTF record with full provenance (which copy
5′, which 3′, breakpoint, read support) instead of an attribute — useful for studying gene
conversion, and the one cross-copy event a linear assembler cannot represent at all. On a
standard transcript benchmark it is an *additional* (structurally-duplicate) prediction, which is
why it is opt-in and confirmed-only.
