# O3 as an RNA-only "look harder" screen: synthesis of the vault survey (2026-10-22)

**Status:** synthesis of four survey files (`docs/o3survey_2026-10-22/`, angles A–D) plus two quick checks run on the frozen O2-VC
outputs. It is design and exploration only. No screening score was computed on held-out truth, so HG00733 and HG002 stay usable for
evaluation.

## The framing

- **The advisor's view (Canzar).** "Cluster by PSV, BAM tags and find a certain print of a missing copy"; it "does not need to be
  definitive (that is what DNA WGS is for) but should suggest when to look harder".
- **The user's rule (2026-10-22):** RNA reads against one reference only. No DNA reads, no per-sample assembly.

So O3 is a **screening rank**, judged by AUC or by sensitivity at a fixed false-alarm rate. It is not a hard abstention gate for O2.

## Two statements that bound what is possible

**1. An extra copy can be certified from RNA reads** (survey C, statement 1; survey A, fingerprint 1).

- **The test.** Take the read patterns over a family's columns, where each pattern is seen in ≥ 2 reads and two patterns conflict at
  ≥ 2 columns. If the patterns that mutually conflict cannot all fit in the 2 allele slots per copy of the copies the reads tie to,
  the donor carries more sequence than the reference: an extra copy.
- **Properties.** The test is exact, checked by bipartite matching, and uses no thresholds. The "≥ 2 independent observations"
  rule handles noise. Comparing the haplotype count with 2 × the copy count is its weak special case.

**2. An absent copy cannot be certified from RNA alone** (survey C Theorem A; survey A's verdict; survey D's register evidence).

- **The argument.** Suppose the donor carries copy X as h1/h2 and lacks Y. Another donor carries X as h1/h1 and Y as h2/h2. Both
  express the same read sequences inside the shared block, so every RNA test behaves the same on both.
- **What follows.**
  - The PRAMEF7.1 case can only ever be a "consider" flag.
  - H5's lost right placements are forced, not a tuning failure.
- **The single exception.** A read that reaches the copy's own unique sequence (a flank, a unique exon, a readthrough) certifies
  that the copy is **present**. That is survey A's "anchored ownership".

## The screen: four tiers, threshold-free (ranked by tier, then by consequence)

| tier | fingerprint | what it means | source | status |
|---|---|---|---|---|
| 1 certified extra copy | haplotype excess: mutually conflicting patterns exceed the tied copies' allele slots | the donor has a copy CHM13 lacks | C1, A1, B1 | needs each copy's **teacher** alleles. Re-call with longcallR's haplotagged BAM (we ran `--no-bam-output`) |
| 2 consider: novel recurrent haplotype | a recurrent haplotype no copy's teachers show; generalises the current linked foreign-allele flag to whole haplotypes | extra copy, mosaic/fusion, or an untaught allele | C3, D3 | the current flag already fires on HG002's NPIPB2 (donor copy with a 57-bp insertion), SPDYE1 (4 flags; 28 real donor-insertion reads) and the CNV families NCF1, NOMO1, GOLGA8A, SMN2. It does not fire in HG00733's 46 absent-copy families (0/46) and false-alarms in 5 of 128 |
| 2b consider: unexplained recurrent indel | an indel no tied copy explains; a copy just outside the 10-AS tie has fewer edits | a structurally different donor copy (NPIPB2.2) | B2 | screen; also fires on simulated KRTAP21 |
| 3 consider absent | collapse-compatible pair (two copies' reads fit in two haplotypes) AND no read anchors either copy to its own unique sequence | a reference copy the donor may lack (PRAMEF7.1) | C4, A2, A3, H5 | necessary, never sufficient (statement 2). Its false-alarm set (present homozygous pairs) is computable on the truth side before running |

**Dead, do not re-propose** (survey D and the register):
- MEC used as a haplotype count;
- colouring the raw conflict graph;
- the four-gamete test (it never certifies when ≥ 2 copies are tied);
- depth or expression as a copy-number signal;
- soft clips;
- raw NM / `de` (HiFi indel error dominates);
- reference-PSV-only O3.

**Prior evidence for tier 2's lineage** (survey D):
- the linked-cluster to IsoCon chain recovered 425/426 in simulation with 0/152 false;
- real data were weaker (r1292–1296: not claimed; r1304: alleles falsely flagged in 4 of 8 intact families);
- IsoCon detected both of NPIP's maternal-only units with nothing outside (r1363).

## Quick checks run tonight (descriptive, post hoc)

- **The current flag as a family-level screen:**
  - HG00733: 0/46 absent-copy families flagged, 5/128 others flagged. That is expected, since a reference copy the donor lacks
    creates no foreign alleles, and HG00733's simulation contains no donor-only copies.
  - HG002: it fires in the families with real structural differences listed in the table.
- **A family-level haplotype-excess check on tied reads only** (survey A): 0 of 131 families in HG00733 and 0 of 50 in HG002, so
  specific but blind. That includes SPDYE1 and RASA4B, which do contain real donor insertions. The certificate needs the **teacher**
  reads' alleles per copy, not just the tied reads.

## What to build next (in order)

1. **Haplotagged teacher reads.** Re-call longcallR with BAM output, so each copy's teacher reads carry phase. This is light: about
   0.7 s per copy.
2. **The tier 1–3 screen as reported-only O3 outputs.**
   - Tier 1: bipartite-matching certificate.
   - Tier 2: novel recurrent haplotype.
   - Tier 2b: unexplained indel and tie widening by edit count.
   - Tier 3: collapse-compatible pair plus anchored ownership.
3. **Extend `donorsim.py` to simulate donor-only copies** (today it simulates only orthologs), together with absent copies and
   uneven allele expression.
4. **Pre-register the screen on a fresh HPRC donor.**
   - Metric: AUC, and sensitivity at the incumbent flag's own false-alarm rate (a matched operating point, no new constant).
   - Controls: EIF3C as a negative control (r1329); haploid X/Y copies in males as a free control (a het on a non-PAR Y copy is
     always a mixture).
   - HG002: reported, never claimed, because of its liftover-label artifacts.

## Related tonight

- **Naive O2 (`docs/O2_NAIVE_TIED_2026-10-22.md`):**
  - **Highest tied AS:** it is minimap2, adding nothing.
  - **Crediting every tied copy:** it recovers a few expressed copies that primaries miss, but invents more.
- **RNA-only O2 gates (`docs/O2O3_EXPLORE_2026-10-22.md`):** none brings the held-out donor under 2%. The hindsight H5 does, but at a
  forced cost (statement 2).
