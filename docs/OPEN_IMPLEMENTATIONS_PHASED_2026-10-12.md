# Open implementations — phased data (2026-10-12)

**Why this exists.** A collaborator suggested working with phased data. The toy
(`docs/PHASED_TOY_2026-10-12.md`, prereg c5431f41, P0-P3 PASS) shows the principle in the EASY case: on a collapsed
reference a heterozygous allele and an extra copy give identical evidence; with two haplotype assemblies a truth-free
rule separates them 100% (dispersed duplication, unique flanks, K = 2, perfect assemblies, 0.43% divergence).
Nothing below is built. **Build an item only once its trigger is met.**

## Step 0 — ask before building (blocks everything)

> "When you say phased data, which form, and for which individual: (a) haplotype-resolved assemblies (two FASTAs,
> e.g. hifiasm trio/Hi-C), (b) a phased VCF (GT `0|1`, `PS` blocks) against which reference, or (c) a BAM with
> `HP`/`PS` haplotags? Is it the same animal as the IsoSeq (KB3781) or another one?"

| form | O2 (which copy / haplotype a read is from) | O3 (copy absent from one haplotype) | item |
|---|---|---|---|
| (a) two assemblies | yes | **yes** — the only form the toy validated | I1 |
| (b) phased VCF, collapsed reference | only in well-resolved regions | **no** — in a collapsed SD the caller reports PSVs as het sites and the phaser phases them: the toy's collapsed arm | I3 → I2 |
| (c) haplotagged BAM | as (b) — tags inherit the VCF's errors in SDs | no | I2 |

If the collaborator can provide (a), push for it; it is the form that makes O3 testable.

## I1 — two-assembly mode (`--hap-ref hap1.fa,hap2.fa`) · size M · trigger: form (a)

- **What:** promote the toy rule to a tool. Per read: map to each haplotype; NM1 == NM2 → SHARED; else home =
  lower-NM haplotype; ALLELE if the home locus has a counterpart between the same flanks on the other haplotype, else
  EXTRA_COPY (O3 candidate on the home haplotype).
- **Replaces:** the ad hoc `O3_REF=mat|pat` path in `bench/o3_maternal/` (`truth.py` already computes "locus on OTHER
  with no counterpart on REF" from lifts) and `bench/rna_allele/` chrmap logic. Reuse those, do not fork them.
- **Open design points:**
  - Counterpart test: the toy used a fixed ±2 kb flank, which works only because the flank clears the 500-bp
    duplicated pad. Real SDs are tens of kb → use the lift/synteny already in `bench/o3_maternal/truth.py`, not a
    fixed flank.
  - Tandem duplications: "which copy is the extra one" is undefined; report the copy-count difference per array, not
    a per-copy label.
  - Where it lives: new `copy_assign` mode vs `missing_copy_flag` input vs a `utilities` subcommand — decide with the
    O2/O3 owners of those binaries.
- **Tests first:** the toy's S0/SA/SB as fixtures (`bench/phased_toy/phased_toy.py` already generates them, 1 s);
  then a KB3781 family with a known mat/pat difference (GWFAM175_B0, the one copy beyond allele divergence).

## I2 — read haplotype from BAM tags (`copy_assign --hap-tag HP`) · size S · trigger: form (b) or (c)

- **What:** read `HP` (and `PS`) aux tags from the input BAM; carry them as a label into the assignments output
  (read → copy, haplotype, phase set); report per-copy haplotype counts.
- **Do NOT** feed the tag into the assignment decision until validated: in SDs the tag inherits the phaser's
  collapsed-duplication errors (form (b) row above).
- **Validation before use:** on KB3781, tags derived from its own phased variants vs the O2 certificate's mat/pat
  assignment (`bench/asj_o2/`, T1 discordance .0014 held-out) — agreement inside vs outside SDs.

## I3 — phased VCF input · size XS · trigger: form (b)

- **No VCF parser.** Convert with `whatshap haplotag --reference ref.fa phased.vcf.gz reads.bam` → I2. Document the
  command (and the RNA caveat: spliced reads are fine, but allele-specific expression leaves some haplotypes unseen).

## I4 — naming collision with the existing `copy_assign --phase` · size XS · do anytime

- `--phase` WRITES `PS`/`HP` where **HP = the assigned COPY**, not a parental haplotype (help text at
  `src/copy_assign_bin.rs:1076-1087` already says so; register line 640). An input flag about parental haplotypes next
  to it will confuse readers and downstream tools that consume `HP`.
- **Options:** rename the output flag (e.g. `--emit-copy-tags`, keeping `--phase` as a deprecated alias) or name the
  new input flags without "phase". Pick one before I1/I2 land.

## I5 — toy rungs (harder cases) · size S each · trigger: before claiming anything beyond the easy case

Each a pre-registered rung of `bench/phased_toy/` (or a famsim two-haplotype mode, `bench/famsim/`):
1. tandem duplication (no unique flank between copies);
2. gene conversion between haplotypes / between copies;
3. divergence sweep down to the read error budget (PSVs/het sites < errors per read);
4. K = 3 copies — identifiability by separation + linkage fails at K >= 3 (register line 561);
5. imperfect assemblies (collapsed copy in one haplotype assembly);
6. real KB3781 mat/pat sequences of one family instead of simulated variants.

## Limits that travel with any claim

- Copy vs allele cannot be settled on RNA alone (register line 767) — phased **assemblies** are the referee, not the reads.
- Genome-wide whole-gene CN difference between KB3781's haplotypes was uninformative (control 5.4x the signal,
  register line 751) — expect few real mat/pat copy differences; GWFAM175_B0 is the known one.
- Phased data must be of the SAME individual as the RNA for the O3 test; another individual's haplotypes test
  population CN variation, not reference absence in this sample.
