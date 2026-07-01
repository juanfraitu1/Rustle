# RNA family = transcript homology ($E_r$), not read-conflict ($E_c$) — results note

**2026-06-30.** Short numeric companion to `bench/family_definition_formal.md`. The RNA-level multi-copy family is reframed from the read-conflict component ($E_c$, an *ambiguity* oracle) to the transcript-homology component ($E_r$, a *homology* oracle, parallel to $E_a$/$E_b$/$E_p$). $E_c$ + the SUN 3-tier ladder are demoted to the within-family **O2** (copy-assignment) structure. Full argument, lattice, and provenance in `family_definition_formal.md`; this note is the measurement ledger.

## The flaw (one line)

$E_c$ links only loci `minimap2` **fails** to resolve (de-tie / MAPQ-0). A multi-copy family whose copies are divergent enough to each map uniquely produces **no de-tie edge**, fragments into singletons, and vanishes — the globin case (MB↔HBB share no read). So $E_c$-as-family-definition **silently drops the divergent, uniquely-mapping families** — the size-graded, globin-style class measured at **29.8% of multi-copy families** (the O1-drop; §Missed class). (It does *not* keep only the assignment-hard core — most retained families are easy but still de-tie; the 82–96% SUN/MAPQ figures are per-copy O2 resolvability, a distinct quantity, §O2.)

## Substrate

- $E_r$ catalog = `bench/denovo_families.tsv` (POA transcript-homology, core_recip $\ge0.13$): **1,130 multi-copy families / 3,636 expressed de-novo loci**. Excluding the one repeat/ZNF over-merge blob **DNFAM0 (728 loci**, dissolved by the shared $\gamma$-quasi-clique operator $R$): **1,129 families / 2,908 loci**; evaluable (fetchable BAM coords) **1,126 / 2,902**.
- Resolvability measured *inside* each $E_r$ family (the genome-wide $E_c$ catalog is not materialized): `GGO_mm.bam` MAPQ + `bench/sun_identifiability.tsv` (154 families / 412 copies) + `bench/psv_graph_genomewide.json` (154 co-located families). Scripts (run from repo root; outputs to `bench/rna_reframe_*_rows.tsv`, override with `RNA_REFRAME_OUT`): `bench/rna_reframe_measure_detie.py`, `bench/rna_reframe_measure_sigtie.py`, `bench/rna_reframe_measure_ec_er.py`, `bench/rna_reframe_validate.py`.

## Missed class: $E_r\setminus E_c$

| criterion | families dropped | copies dropped |
|---|---|---|
| de-tie (shipped default $\Delta$0.005/$\mathrm{DE_{max}}$0.05/min3) | **335 / 1,126 = 29.8%** | **741 / 2,902 = 25.5%** |
| sig-tie ($E_c^{\mathrm{sig}}$, $\varepsilon$1e-3/$\alpha$1e-3) | **370 / 1,126 = 32.9%** | **828 / 2,902 = 28.5%** |

Size-graded (the flaw's signature): size-2 **33%**, size-3 24%, size-4/5 ~15%, size $\ge7$ **0%**. Dropped class is genuine paralogy — roots ZNF/ZBTB/TRIM/KRAB-ZNF/TMEM/PNMA/SLC4A; **266/335** contain protein_coding; **48 unnamed/novel** de-novo families (invisible to annotation *and* to $E_c$).

Hard-tie (MAPQ=0) reading — **per-copy O2 resolvability, NOT a larger O1-drop.** The shipped $E_c$ edge is the de-tie predicate, not MAPQ=0, and a de-tie is far more permissive than a hard tie, so the O1-drop stays **29.8%**. For O2 context: $E_r$ loci **0.51% MAPQ=0** (99.49% uniquely resolvable), genome-wide **0.04%** (**99.96% MAPQ$>0$**); ambiguity concentrates **14×** in the co-located subset (7.4% MAPQ=0). (An earlier helper `rna_reframe_measure_ec_er.py` reported a **96.6% "EC_DROPPED"** from a $\ge5\%$-MAPQ0 per-locus proxy — that is O2 hard-ambiguity, **not** de-tie edge-formation, and overstates the O1-drop ~3×; the faithful measure is the 29.8% above.)

## $E_c\subseteq E_r$ (the one containment, asymmetric — machine-checked)

`bench/rna_reframe_validate.py` (§B–§C), 154 SUN families mapped onto $E_r$: among **132** with a multi-copy $E_r$ label, **129/132 = 97.7%** have all co-located ambiguous copies in ONE $E_r$ family under the shipped **symmetric** core_recip$\ge0.13$ operational oracle. The **3** splits are **operational shared-exon leaks, NOT EDGE_LINKED** (retracted): `denovo_family_edges.tsv` carries **0** cross-family `core_recip` edges touching these loci (families *are* its components), yet the loci genuinely de-tie (fam16 **20**, fam34 **121**, fam65 **65** reads) — their shared homology is a single conserved exon/domain **below the 0.13 reciprocal floor**, so the symmetric operational oracle misses them (~2.3%). They are **absorbed by the DEFINITIONAL permissive-local $E_r^{\mathrm{asym}}$** (each de-tie read is itself a significant local-alignment witness), so the theorem $E_c\subseteq E_r^{\mathrm{asym}}$ **stands** — only the operational 0.13 oracle leaks, and $R$ cannot fix it ($R$ only splits, never merges). Two-tier: **permissive $E_r^{\mathrm{asym}}$** contains every witnessed $E_c$ edge; the **operational 0.13** oracle contains 97.7%. Asymmetry is honest (permissive local homology: exonic at $\ge1$ endpoint, homologous at the other — the RNA mirror of $E_b$'s asymmetry; forced by `read_conflict.rs:86–106` checking only `de`/`aln_len`, no partner-exon overlap). Chain: $E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_r^{\mathrm{asym}}$.

## O2 per-copy resolvability inside these families (SUN, 154 families / 412 copies) — NOT the O1-drop

Tier-1 (per-copy resolvable) **338 = 82.0%** · Tier-2 **1** · Tier-3 (collapsed, $E_c$-visible) **73 = 17.7%**. **132/154 = 85.7% are fully-Tier-1** — i.e. **every copy is individually read-resolvable**, a *copy-assignment* (O2) fact that does **NOT** measure the O1 family-drop. "Fully-Tier-1" does **not** imply an edgeless $E_c$: a single distinguishing SUN moves $de$ by only ~$1/\text{read\_len}\le0.005$, so SUN-resolvable copies still de-tie. **Measured** (`rna_reframe_validate.py` §D): **132/132 = 100%** of these fully-Tier-1 co-located families still form a de-tie edge and **SURVIVE** $E_c$ — **0** vanish. So the earlier "fully-Tier-1 $\Rightarrow$ edgeless-$E_c$ $\Rightarrow$ vanish $\Rightarrow$ captures $\approx1/8$ / drops ~86%" is **retracted** (empirically inverted; on this co-located catalog the $E_c$ O1-drop is ~0%, well below the full-catalog 29.8%). The **sole** O1-drop headline is the de-tie **29.8%**; 82–96% is per-copy O2 resolvability. `psv_graph_genomewide.json`: **82.5% per-copy-resolvable / 17.5% frontier(K=0)** (again O2, not an O1-drop).

## Updated 4-level lattice

O1 homology oracles (pairwise incomparable): $E_a$ (genomic) — $E_b$ (exonic) — **$E_r$ (spliced transcript)** — $E_p$ (protein). $E_c$ leaves O1 and becomes an O2 substructure of $E_r$. Clean-containment chain moves down one level from $E_b$ to $E_r$:

$$E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_r^{\mathrm{asym}}\qquad(\text{sole unconditional containment; tightest, since }E_c,E_r\text{ share }V_R;\ E_c\subseteq E_r\subseteq_{\mathrm{proj}}E_b).$$

$E_r$ vs $E_b$ = **symmetric refinement, incomparable** (not equal): $E_r$ = symmetric "exonic-at-BOTH" core of the asymmetric $E_b$, **plus** a spliced-substrate superset direction (retrocopies: EEF1A1/CNN2, OCLN~SEPTIN7, BCAS4~CCDC30 — clean transcript-to-transcript but no contiguous genomic block).

## O1/O2 separation and the prior-unification re-read

O1 = homology (chance-match null, includes easy copies); O2 = resolvability ($E_c$/SUN/$\chi(H)\le\Lambda$, sequencing-error null $\varepsilon$, K-frontier abstains). This **partially undoes** the prior "one significance criterion for O1+O2" (`CONFLICT_EDGE_UNIFICATION`), which baked the flaw in by making the family edge = the gate's own significance-tie. **Kept:** $E_c^{\mathrm{sig}}\subseteq E_c$ as an *O2-internal* refinement. **Re-read:** the GGO 81$\to$71 "narrowing" = a *within-family O2 sharpening* (fewer confusable pairs), **not** an O1 family-count change; the 10 "lost" families are fully-resolvable O1 $E_r$ families.

## Implementation implication

Shipped family catalog = the homology / `denovo_families` grouping ($E_r$, $R$-refined to dissolve DNFAM0), **not** the read-conflict graph. `src/rustle/vg_family/read_conflict.rs` is unchanged — it correctly stays the O2 carrier (`MCC=χ(H)`, exact-decomposition unit), consumed *inside* each $E_r$ family.

**Verdict:** fix the definition. The current read-conflict family is O2's substructure, not O1's family; the reframe recovers the missed ~30% while keeping $E_c$ as the correct within-family `MCC=χ(H)` carrier.
