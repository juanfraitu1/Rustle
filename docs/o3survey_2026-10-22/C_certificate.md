# C. A combinatorial missing-copy certificate as an O3 screening score (angle C, 2026-10-22)

**Status.** Design only. Nothing was computed on any truth table, and nothing was run on the frozen outputs (so that HG00733 and HG002 stay
clean for the evaluation). Read: BRIEF; prereg Results; EXPLORE §0-§5; expected-failures §0, F1-F11, §3; memory `feedback_o3_rna_only`,
`project_samplecopy_2026-10-19`, `project_o3_npip_maternal_copy_2026-10-21`, `project_mec_hidden_copy_nogo` (08-30),
`project_rna_haplotype_count_parked`, `project_o3_psv_clusters_2026-10-15`; register lines 579/617/624/625/640/642 and r1175, r1194, r1196,
r1209-1212, r1292-1296, r1336-1337, r1363; vault notes and full text for Chaisson 2017, SDA, longcallR, pHapCompass, Paraphase 2023/2025,
Canzar 2016, IsoCon and devider. "(inf.)" marks my own inference. "(std.)" marks a standard theorem that is not in the vault.

## 0. Bottom line
1. **The extra-copy direction has a threshold-free certificate:** a **Hall deficiency** δ(Q) > 0.
   - Take a clique Q of pairwise-conflicting certified read patterns. If Q cannot be matched injectively into the diploid haplotype slots of
     its tie-set copies, the reads need ≥ δ(Q) haplotypes beyond what the reference copies can host.
   - It is valid under an explicit error model (E2) and an explicit no-extra-locus hypothesis (H0), and is computed by bipartite matching.
2. **The absent-copy direction has no certificate, and none can exist (Theorem A).**
   - Every RNA-only output is a function of the multiset of read sequences plus CHM13.
   - "One heterozygous copy, paralog absent" and "two homozygous copies" give the same multiset once expression is matched.
   - So on matched pairs P(flag | absent) = P(flag | present) for every test. H5's cost (66 / 108 right placements lost, EXPLORE §5.3) is
     therefore **forced**, not a tuning failure.
   - The best possible is a **necessary condition**: the "collapse-compatible" pair (a bipartite certified conflict graph). Its recall is
     provable, and its false-alarm set can be predicted.
3. **The other bounds the prompt lists are weak or dead:**
   - the prompt's own bound, ω > 2n (the clique number against twice the copy count), is valid but almost never fires;
   - four-gamete never certifies in O2's domain;
   - MEC cannot count haplotypes (Chaisson; the 08-30 no-go);
   - raw conflict-graph colouring is dead (line 642).
4. **Score:** lexicographic tiers, with no constant except "every fact is observed twice". This principle is already in O2-VC's flag.
5. **Evaluate it as an identifiability prediction** (the r1209 design).
   - From the truth side, compute in advance which families and copies the certificate *can* see; it must fire on exactly those.
   - HG00733 has no donor-only copies, so the certified tier must give 0 there.

## 1. Model and definitions (O2-VC objects; `calls.tsv`, `reads.tsv`, `columns.tsv`, `vcf/`)
- **Tie set and teachers.**
  - T(r) is the set of copies within 10 AS of read r's best alignment, in the O2-VC domain (no locus off the family).
  - A teacher has T(r) = {c} at margin ≥ 10.
  - A tie component is a connected set of copies linked by tied reads; **n** is its size.
- **Columns.** O2-VC columns (`separating_call`, `separating_ref`, `erased_ref`, `foreign`) minus:
  - A>G changes in the transcript frame;
  - `RnaEdit`, `Dense` and unknown (`.`) sites.

  Bases count at Q ≥ 10. Every column is caller- or reference-validated.
- **Patterns.**
  - A pattern P = (S, x) is an allele vector x on a column set S.
  - Its *support* is the set of reads that cover S and equal x on S.
  - P is **certified** when ≥ 2 distinct reads support it.
  - Vertices are the agreement restrictions of read pairs, deduplicated. Patterns that differ at exactly one column are merged by single
    linkage.
  - Two patterns **conflict** when they disagree at **≥ 2 columns**. G₂ is the certified conflict graph.
- **Slots and lists.**
  - Capacity κ(c) = 2, or 1 for male non-PAR X/Y. The slots of copy c are (c, 1..κ(c)).
  - A pattern's list L(P) is the union of the slots of T(r) over its supporting reads r.
  - **Leak closure:** a teacher pattern of c also lists the slots of every copy whose teacher patterns are compatible with it.
- **E2 (error model):** no read pair carries identical sequencing errors at two columns.
- **H0 (no extra locus):** each read's source haplotype is hosted at a copy in T(r), and each copy hosts ≤ κ(c) haplotypes.

## 2. Statements
| # | statement | proof sketch | fails when |
|---|---|---|---|
| S1 | **Clique bound.** If k certified patterns pairwise conflict, the reads come from ≥ k distinct haplotypes. | A 2-column disagreement between patterns that each have 2 supporting reads cannot come from one haplotype without a read pair sharing errors at 2 columns (E2), or two independent coincidences at different columns. The merge rule removes the latter: a pattern one column from another is not a separate vertex. | Recurrent artifacts that are not errors: editing outside A>G, template-switch chimeras, Kinnex concatemers (§5). |
| S2 | **Exact count and polynomial bounds.** With 1-column conflicts, min #haplotypes = χ(G) (Helly: pairwise-agreeing patterns agree column by column, so they have a common extension). With the 2-column rule, ω ≤ ⌈ϑ(Ḡ₂)⌉ ≤ χ (std., Lovász sandwich); ϑ is polynomial. | Colour classes are haplotypes. | χ is NP-hard and inapproximable (line 625). Instances are tiny after deduplication; use ω and Hall as the certificate and χ only as the report. |
| S3 | **Monotone and conservative.** Adding columns or reads never lowers ω or δ, and merging never raises them. | Supergraph and subgraph. | Never. Restricting to O2's columns loses power, never validity. |
| S4 | **Hall certificate.** For a clique Q, δ(Q) = \|Q\| − ν(Q → slots) = max over Q' ⊆ Q of (\|Q'\| − \|L(Q')\|) (std., Ore's deficiency). δ > 0 refutes H0 ∧ E2. | Under H0, distinct haplotypes need distinct slots. | Unfaithful ties: a donor SV pushes the true copy out of the 10-AS window (NPIPB2.2's 57-bp insertion), or teacher leakage escapes the closure (F1 with an untaught source, F5). |
| S5 | **Family tier.** With every list set to all 2n slots, δ = ω(G₂) − 2n: the prompt's bound. | S4 with full lists. | Powerless. It needs ≥ 2n + 1 mutually conflicting expressed haplotypes. |
| S6 | **Four-gamete.** Four gametes at columns (i, j) give k = 4 in S1. O2's domain has n ≥ 2, so 4 ≤ 2n and the test never certifies alone. Perfect-phylogeny violations count conversions, not copies. | Arithmetic. PP25: more than 75% of African PMS2/PMS2CL haplotypes are converted (expected-failures F7). | Not applicable. |
| S7 | **MEC cannot count.** Chaisson §3: "Since the MEC objective itself decreases monotonically with ploidy, it is not possible to estimate the ploidy using the MEC objective." Under missing data MEC also shrinks the scored cells (08-30: GWFAM115 k=6 reaches MEC 0 with 15.3% of cells scored). | Not needed. | Use only MEC = 0 feasibility (S2), never a penalised K. |
| **A** | **Absence is not identifiable from RNA.** D1: copy X is h1/h2 and Y is absent. D2: X is h1/h1 and Y is h2/h2, both inside the X/Y block. For any expression in D1 there is one in D2 with the same law of read sequences. Alignments, ties, teachers, longcallR calls and `columns.tsv` depend only on read sequence and CHM13. So every test φ has P_D1(φ) = P_D2(φ). | IsoCon says it "has no functionality … that would distinguish allele-specific variants from distinct gene copies". Theorem A says that, for absence, no method can. | A read reaches copy-unique sequence (a flank, unique exon or readthrough). That read is untied and certifies **presence** only. Expression is constrained: it is not (ASE, F4). |
| A1 | **Collapse-compatibility is necessary.** Let R(A,B) be the teachers of A and B plus reads with T ⊆ {A,B}. If A or B is absent and the other is diploid, G₂[R(A,B)] is bipartite (a linear-time check). | ≤ 2 source haplotypes. | Leakage from a third copy C makes the graph non-bipartite and lowers recall. Its false alarms are exactly the present pairs that express ≤ 2 distinct haplotypes: both homozygous, or one silent. |
| A2 | **What is certifiable.** χ(G₂[R(A,B)]) ≥ 3 means ≥ 2 loci feed A∪B: the pair is not a collapse. The prompt's "the pair's columns never co-occur discordantly" holds in both D1 and D2, so it is not a certificate. | — | A third-copy leak gives a true "≥ 2 loci" that is the wrong locus. It fails toward VC's status quo. |

PRAMEF7 (.0, .1) in HG00733 fits A and A1 exactly:
- two haplotypes, each copy homozygous (four PASS 1/1 records on .1, no het);
- collapse-compatible, so "consider absent";
- H5's present-pair losses sit in the same class, since H5 is A1 plus a "≥ 2 columns" constant (inf.).

## 3. The screening score (rank = tier first, then the key; no constants)
| tier | verdict | condition | within-tier key | meaning to report |
|---|---|---|---|---|
| 1 | **certified extra** | max over Q of δ₂(Q) > 0 (S4, κ-slot lists, leak-closed) | δ₂, then the clique's read support | ≥ δ₂ expressed haplotypes the tie-set copies cannot host, or a recurrent artifact (§5). DNA places it. |
| 2 | consider extra | δ_obs > 0: as S4, but a copy's capacity is the number of distinct leak-closed teacher haplotypes actually observed (from the VCF's GT/PS) | δ_obs, then support | Tied haplotypes no teacher carries: the haplotype-level foreign flag. It **includes mosaics**, novel combinations of known alleles that the per-column flag cannot see. The excess δ_obs − δ₂ is bounded by the free slots, i.e. untaught second alleles (F5). |
| 3 | consider absent | the pair (A, B) is collapse-compatible (A1), and VC would place reads on it | the number of tied reads VC places on A or B (prioritise by consequence) | Necessary for absence, never sufficient (Theorem A). The expected false-alarm set is computable on the truth side. |
| 4 | none | non-bipartite A∪B, δ_obs = 0 | — | Certified "not a collapse" (A2). |

## 4. Fingerprints, items (a)-(g)
| FP | (a) signal | (b) scenario | (c) literature | (d) confounders | (e) form | (f) compute | (g) register |
|---|---|---|---|---|---|---|---|
| C1 Hall-deficient clique | δ₂ > 0 | Donor-only copy with saturated tie copies (SPDYE1-like extra) | PP25: "an extra haplotype … indicates a copy number gain"; Chaisson CC (copy count as output) | Template switching, editing outside A>G, tie-set loss from SVs (NPIPB2), leaks from an untaught source | Combinatorial, exact (matching) | `columns.tsv` patterns + VCF GT/PS teacher haplotypes + `reads.tsv` tie sets | New. Raw-graph colouring is dead (line 642); this is certified and list-constrained. |
| C2 ω(G₂) > 2n | Family clique exceeds diploid capacity | Same, unsaturated lists | Paraphase CN = "number of unique haplotypes, adjusted by depth"; RNA has no depth adjustment, so it is a lower bound | As C1, minus the tie-set issue | Combinatorial | As C1 | New; the special case of C1 |
| C3 novel recurrent haplotype | δ_obs > 0 | Extra copy, untaught allele, or mosaic/fusion | SDA repulsion "PSVs are mutually exclusive across all of the sequence reads"; DM's discordant reads | F5 (untaught allele: present copy), other apes' or other loci's transcripts (r1295) | Combinatorial; ≥ 2 support | As C1 | Overlaps the O2-VC foreign flag (C1 FAIL after relabel) and r1292-1295 (real data NOT CLAIMED) |
| C4 collapse-compatible pair | bipartite G₂[R(A,B)] | Absent copy taught by its paralog (PRAMEF7.1) | IsoCon quote; r1194: "RNA cannot separate homozygous from hemizygous" | Homozygous expressed copies (forced, Theorem A), a silent copy | Combinatorial, linear time | VCF het/hom per copy + `columns.tsv` | H5 is its thresholded form (dev numbers); G4 swap is inert, consistent with A |
| C5 rejected forms | raw ω/χ; MEC-K; four-gamete; EM | — | Chaisson S7 | Error edges | — | — | Killed: line 642 (3× K, about 58 colours for K = 2), 08-30 MEC no-go, S6, line 579 (EM 0/3,081) |

## 5. Noise without thresholds (the principle: every fact is observed twice: a vertex by 2 reads, an edge by 2 columns)
| source | effect | handling | residual (stated, not tuned) |
|---|---|---|---|
| Sequencing error | Spurious patterns. With 2-read support alone, the expected number of coincident-error pairs is Σ_j C(d_j, 2)(e/3)², ≈ 4 in a 2,000-read × 20-column family at e = 10⁻³ (inf.). | 2-column edges plus 1-column merging push a false clique to E2 violations, Σ C(d, 2)(e/3)⁴ ≈ 10⁻⁴ (inf.). Measure e from the reads' own mismatch rate off the columns, and **report the budget per family**. | True haplotypes one site apart merge (DAZ2/3, r1336): the identifiability limit without a rate model. |
| RNA editing | Recurrent A>G "haplotypes" | Exclude A>G in the transcript frame and RnaEdit sites (already O2-VC policy) | Non-A>G editing (rare; F6) |
| Chimeras and template switching | A recurrent mosaic of two observed haplotypes | Report "single-switch-explainable" patterns separately; they never enter tier 1 | Real mosaic extra copies fall to tier 2 |
| Teacher leakage (≈ 2%, F1) | False saturation, i.e. false tier 1 | Leak closure (S4): a teacher pattern compatible with a rival's may be hosted there | Leak whose source haplotype is untaught (F5) |
| ASE / monoallelic expression | Hides haplotypes | Harmless for validity (S3); costs power | It is what makes Theorem A bite |
| Donor-vs-reference SVs | A true copy leaves T(r) | Report the "present-divergent" class separately (NPIPB2.2) | O3-relevant reference bias, not an extra copy |
| Liftover and chains | Truth labels only | Score with the controller's verified list | 29/133 `no_counterpart` reads sit in chain gaps (EXPLORE) |

## 6. Literature check (cited only from what I read)
| paper | supports | does not support |
|---|---|---|
| longcallR | "weighted minimum error correction problem, which is NP-hard … may miss the global minimum" | Diploid, K = 2 fixed, MAPQ ≥ 10: no haplotype counting, no multi-copy |
| pHapCompass | Valid-phasing pruning by genotype; uncertainty | K and genotypes given (note §Assumes): it cannot count |
| SDA | Attraction/repulsion PSV graph; copy number as the cluster count | Thresholds (log10 ratio ≥ 1.5, < 10% of mean coverage) and depth to tell PSVs from alleles; RNA has no such depth |
| Chaisson 2017 | MEC equals maximum likelihood; CC discovers the count | MEC monotone in ploidy (S7); CC needs "uniquely identifying variants" and unstated c_min/c_max/d |
| Paraphase 23/25 | Count unique haplotypes: an extra one is a CN gain | Depth adjustment and identical haplotypes ("two of the three haplotypes being identical", PP23) need DNA |
| Canzar 2016 | Conflicts "if and only if the corresponding intervals overlap"; interval graphs are perfect, so ω = χ (std.) | Facility location, no abstention; our G₂ is not interval or chordal (line 642: 23% chordal), so ω < χ is possible |
| IsoCon / devider | Alleles vs copies stated open; devider has no preset count and abstains on ties | devider arm REFUTED for O3 (reference-bound, memory 10-02/03) |

## 7. Register overlaps
| row | what it killed | consequence here |
|---|---|---|
| line 579 / 617 / 624 | EM / soft SDA: 0 of 3,081 decisions changed; the objective decomposes | The certificate is an infeasibility proof, not a likelihood: no EM |
| line 625 / 642 | χ(H) is inapproximable; the raw graph needs about 58 colours for K = 2 | Only the certified G₂ is used; ω and Hall are the certificates |
| 08-30 MEC no-go | k+1 MEC is degenerate; the clean negative drops hardest; about 0.2 expected true positives | No MEC. The power warning stands: positives are rare |
| r1175 | O2 is a step function in m; m ≤ 1 abstains | Same ceiling: a pattern must conflict with every slot-holder at ≥ 2 covered columns |
| r1209-1212 | Y chain: copies a few SNVs apart are invisible at delta | With κ = 1 on the Y there is no delta, but the 2-column rule makes 1-site copies invisible. r1336's two RBMY1 groups sharing one label is a κ = 1 Hall violation already seen; r1337 was UNDERPOWERED |
| r1196 | A missing copy's reads swamp the absorbing copy (native base 2-6%, under a 0.20 floor) | S4 has no allele-fraction floor: ≥ 2 reads suffice |
| r1292-1296 | Simulation PASS; real data NOT CLAIMED; 35/57 intact clusters are other loci | Keep the O2-VC domain (family-best); expect contaminant transcripts in tiers 1-2 (O3 = detection, user 10-09) |

## 8. Evaluation design (screening rank; pre-register before computing)
| arm | unit | positives | negatives | metric | pre-registered prediction |
|---|---|---|---|---|---|
| E-abs HG00733 | copy pair | Copies with no ortholog on either haplotype (`orthologs.tsv`; build.json: 61 without an ortholog, 8 silent) that have an O2-domain paralog | Pairs with orthologs on both haplotypes | Tier-3 recall; AUC of the tier-3 key | Recall = the leak-free fraction (≈ 1). **False-alarm set = present pairs that express ≤ 2 distinct donor haplotypes**, computed in advance from the haplotypes and the simulated reads. An equality test, not a bar. |
| E-ext HG00733 | family | none (donorsim simulates no donor-only copies) | all families | Tier-1 count | **0.** Each one is a falsifier, classified by mechanism. |
| E-ext HG002 | family / read | 5 families, 42 verified reads (SPDYE1 28, RASA4B 8, USP17L11 3, POTEM 2, BMS1P1 1) | Families whose tied reads all have a donor counterpart | AUC; sensitivity at the **incumbent foreign flag's own false-alarm rate** (a matched operating point, no new constant) | Tier 1 fires only where the tie copies are saturated (identifiable set computed from the donor first). NPIPB2.2 is its own class. HG002 is a used donor: report it, do not claim from it. |
| E-fresh | both | Never-used HPRC donor (HG01978 or HG02257, EXPLORE §5.4). Extend donorsim to simulate **donor-only copies** (reverse orthologs), with beta-binomial ASE and an editing/chimera injection arm. | As above | As above | Same equalities; the injection arm measures the §5 residuals |

**Falsifiers:**
- any tier-1 call on HG00733 or on E-fresh's no-extra families;
- a tier-3 false-alarm set unequal to its predicted set;
- tier 2 not beating the foreign flag's AUC on E-fresh.

## 9. Computing it from the frozen outputs (a light job, about a minute of Python; not run)
- **Tied-read patterns.** Read them from `columns.tsv`. A column's key is the tuple of copy positions; work in the read frame.
- **Teacher haplotypes.** Take them from `vcf/<copy>.vcf`: hom-only gives 1; a phased 0|1 within one PS gives 2 per block; an unphased
  LowQual 0/1 gives "2, phase unknown" and stays compatible. Teacher read pileups at the columns would be stronger (inf.).
- **Lists.** Build them from `reads.tsv` tie sets.
- **Per tie component:**
  1. deduplicate and merge patterns;
  2. build G₂;
  3. enumerate maximal cliques (Bron-Kerbosch on tens of vertices);
  4. compute δ₂ and δ_obs by Hopcroft-Karp;
  5. run the A1 bipartite check per tied pair;
  6. emit the tier, the key and the E2 budget.
