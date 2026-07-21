> **⚠ SUPERSEDED (2026-07-20).** This reply discusses a profile-HMM assignment proposal that was abandoned. The current method (significance-certificate assignment, no HMM) is in `bench/rustle_mechanism.html`. Kept for history.

# Response to your email

> "I don't understand your motivation of using profile HMMs to assign a read to a paralog… profile HMMs capture the alignment of multiple sequences… so this might be useful for checking if a new copy belongs to a family, but I don't see how this can help you in assigning a read to a paralog. Second, I also don't understand your idea of building a profile HMM per copy, and how you would use it to solve the assignment problem. Explain these in isolation, without mixing them up with other concepts such as variation graphs etc."

You are right on all your points. Here is the honest mapping, point by point, with the data that backs each answer.

---

## Point 1.  *"Profile HMMs are for family membership, not copy assignment."*

**You are correct.** A classical profile HMM (HMMER-style) built from an MSA of multiple copies cannot discriminate between copies. At every SNP column the MSA has multiple alleles, so the match-state emission becomes degenerate (≈ 0.5 / 0.5 for two alleles in two copies). The model literally absorbs the variation it would need to use for discrimination.

**Direct evidence:** `figures/fig_hmm4_family_vs_copy.png`
- Left panel: family HMM (MSA of all copies). SNP columns become 0.50/0.50.
  ONE model, so ΔlogL is mathematically undefined → assignment impossible.
- Right panel: per-copy HMMs. SNP columns are pure (0.95 vs 0.017).
  TWO models, so ΔlogL is defined.

So we agree with you that the standard "family HMM" cannot do read-to-copy assignment. That was never our proposal.

---

## Point 2.  *"What do you mean by per-copy profile HMM?"*

We do not mean a classical profile HMM built from an MSA. We mean a much simpler object: **a per-position emission model derived from one copy's transcript sequence**.

For a copy `c` with transcript sequence `S_c`, position `i`:

> `P(base = S_c[i] | copy c)        = p_match  ≈ 0.95`
> `P(base ≠ S_c[i] | copy c)        = (1 − p_match) / 3  ≈ 0.017`

For a read `R` of length L aligned to copy `c`, with M matches:

> `logL(R | copy c)  =  M · log(0.95)  +  (L − M) · log(0.017)`

That is the entire "per-copy HMM." It is a per-base independent-emission model. No MSA, no Viterbi over a profile of multiple sequences.

**Direct evidence in isolation (no VG):** `figures/fig_hmm2_scoring.png`
- Toy example: two copies of length 6 with SNPs at positions 3 and 5.
- Per-copy emission heatmaps shown side by side.
- ΔlogL computed and shown as a bar chart per read.

**This is mathematically equivalent to running per-copy Smith-Waterman alignment with specific scoring constants** (match = +log(0.95), mismatch = log(0.017)). minimap2 already does this for us via its AS score — which directly proves you were right that the HMM framing adds nothing essential.

---

## Point 3.  *"How does this solve the assignment problem?"*

It does not, in isolation, for our use case. **That is the honest answer.**

Once we tested on real GOLGA8I FLNC IsoSeq reads, the picture was:

a) **For FLNC reads on coherent multi-copy families with distinct copy paths, minimap2 alone solves the assignment.**
   - `figures/fig_flnc_assignment.png` — 4 simulated FLNC reads from 4 different GOLGA8I copies, all aligned to a 7-transcript reference with `minimap2 -ax map-hifi`. **All 4 → MAPQ = 60, correct copy.** No HMM needed.

b) **Where minimap2 ties (MAPQ = 0), per-copy HMM scoring also ties.** Both methods agree at the equivalence-class boundary. The HMM does *not* rescue what minimap2 cannot resolve.
   - `figures/fig_class_boundary_resolution.png` — two cases (GAGE class_2 with 10 members, GOLGA8I truncated copy 171) where minimap2 reports MAPQ = 0 and HMM also ties. Neither method can break the within-class tie because of Theorem 1.

c) **The one thing that does close the residual gap is a transcript-length consistency term**, not an HMM:
   - For FLNC reads, the read length should match the source transcript length.
   - `figures/fig_em_bootstrap.png` — without length term: 171 collapses into 071 (sequence too similar due to paralog repeats); with length term: all 7 copies recovered correctly from 100 simulated reads.

So: **the HMM step in the original proposal is unnecessary for our actual data.** What does load-bearing work is alignment + transcript-length consistency.

---

## Point 4.  *"Explain in isolation, without mixing up with variation graphs."*

The HMM explanation in Points 1–3 above contains no VG concepts. The VG is a *separate* tool used for a *different* purpose:

- **What HMM scoring tries to answer:** "Given a read, which copy did it come from?" Per-position match likelihood. No graph involved.

- **What the VG framework answers:** "Which copies *can in principle* be distinguished, and which collapse into equivalence classes by sequence identity?" This is purely topological — it does not score reads, it defines the resolvability floor.

These are two separate questions. The VG is mathematically separable from per-copy HMM scoring.

---

## Summary of where we ended up

| Your concern | Our (now-honest) position |
|---|---|
| Profile HMMs from MSAs cannot assign reads to copies | We agree. We never meant a classical profile HMM. |
| You don't see how per-copy HMMs help | They are formally equivalent to per-copy alignment scoring. minimap2 already does this work for FLNC reads. |
| The HMM proposal seems redundant | It IS redundant for our data. The actual methodology is: VG → equivalence classes (resolvability floor) + minimap2 (alignment) + FLNC-length term (closes residual ties from paralog repeats) + EM bootstrap (distributes any remaining ambiguous reads). |
| Mixing in variation graphs confused the picture | The VG defines what's resolvable (a floor). The alignment + length term does the assignment. They are now presented as two separate components. |

## Supporting evidence (in order to bring to the meeting)

1. **`fig_hmm2_scoring.png`** — what we mean by per-copy HMM, in isolation (toy 6-base example).
2. **`fig_hmm4_family_vs_copy.png`** — confirms that family HMMs (MSA-based) cannot discriminate copies; per-copy HMMs can.
3. **`fig_flnc_assignment.png`** — for real FLNC reads on GOLGA8I, minimap2 alone gives MAPQ = 60 on every read.
4. **`fig_class_boundary_resolution.png`** — where minimap2 reports MAPQ = 0, HMM ALSO ties; Theorem 1 explains why.
5. **`fig_em_bootstrap.png`** — without FLNC-length term, source 171 is lost; with the term, full recovery.
6. **`fig_advisor_floor.png`** — VG equivalence classes (separate concept, the resolvability floor).
7. **`fig_advisor_pipeline.png`** — the actual pipeline we end up running: VG → minimap2 → FLNC-length + EM.
8. **`fig_cross_family_specificity.png`** — zero false positives in family clustering across 6 unrelated families.
