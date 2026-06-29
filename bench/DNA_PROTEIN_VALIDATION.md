# Dna Protein Validation (consolidated)

> Merged from 5 source docs (verbatim, git keeps the originals' history). Each section below was a separate `bench/*.md`.

**Contents:** [dna_psv_catalog_summary](#dna-psv-catalog-summary) · [dna_rna_overlay](#dna-rna-overlay) · [compara_fetch](#compara-fetch) · [compara_validation](#compara-validation) · [transcript_validation](#transcript-validation)


---

## dna_psv_catalog_summary

# DNA-derived PSV identifiability catalog (Phase 1)

- co-located classified pairs: **1387** -> resolvable **331** (24%), genuine-K=0 **1056** (76%). NOTE: this is the DNA reference universe (all aligned co-located pairs, including unexpressed identical tandem copies the RNA census never observes); on the **137** pairs the RNA census actually classifies, DNA and RNA agree **86%** and both find that expressed subset mostly resolvable.
- pairs excluded from K: **14262** unaligned (copy did not align to ref0 — divergent/short paralog), **14** unannotated (no overlapping GFF exon)
- **cross-check DNA-K=0 vs RNA-K0** on 137 census-classified pairs: concordance **86%** (confusion {(True, False): 14, (True, True): 8, (False, False): 110, (False, True): 5})
- discordant DNA-K=0 ∧ RNA-resolvable: **14** (candidate: indel / splice-shift pseudo-K=0 — substitution-only PSVs miss it; Phase-2 private_exon_bp will test this)
- discordant DNA-K≥1 ∧ RNA-tied: **5** (candidate: PSV in a poorly-expressed exon — reference identifiability ≥ read identifiability)


---

## dna_rna_overlay

# Two-tier overlay: RNA layer on DNA/protein-defined multi-copy families

**DNA tier** = protein clusters (mmseqs2 easy-cluster, ≥30% id / 50% cov, on translated CDS) →
formal multi-copy families, INCLUDING ancient/diverged families the RNA-similarity definition
missed. **RNA overlay** = per copy, is it transcribed (real GGO IsoSeq), and how well.

- DNA multi-copy families (protein clusters ≥2): **3,587** (14,545 gene copies)
- RNA-tier families (POA, for comparison): 1,337

## What actually transcribes (genome vs transcriptome)
- of 14,545 copies in DNA-defined families: **transcribed (≥5 reads): 10,490 (72%)** ; well-expressed (≥40): 6,824 (47%) ; **silent: 4,055 (28%)**
- i.e. the DNA tier enumerates every copy; the RNA layer shows which are live — many copies are
  present in the genome but transcriptionally silent.

## Ancient-family gain (the whole point of the DNA tier)
- curated families recovered by the DNA tier that the RNA tier MISSED: **DEFB*, SIGLEC***
| family | DNA tier | RNA tier | DNA copies | transcribed |
|---|---|---|---|---|
| APOBEC3 | YES | YES | 4 | 4 |
| CRYBG (ANCIENT) | no | no | 0 | 0 |
| DAZ | no | YES | 1 | 1 |
| DEFB* (ANCIENT) | YES | no | 4 | 0 |
| MAGEA* | YES | YES | 5 | 5 |
| PRAMEF* | YES | YES | 4 | 0 |
| RABL2 | YES | YES | 2 | 2 |
| RFPL | YES | YES | 4 | 3 |
| SIGLEC* (ANCIENT) | YES | no | 3 | 2 |
| TAS2R* (ANCIENT) | YES | YES | 18 | 7 |

## Per-family 3-number summary (sample: largest + curated)
| DNA family | copies | transcribed | well-expressed | example members |
|---|---|---|---|---|
| DFAM0 | 501 | 83 | 4 | LOC101123691, LOC101123789, LOC101123793, LOC101124039, LOC101124044… |
| DFAM1 | 229 | 219 | 180 | KRBOX5, LOC101123988, LOC101124084, LOC101124732, LOC101124778… |
| DFAM2 | 136 | 134 | 108 | LOC101126065, LOC101126415, LOC101127631, LOC101128844, LOC101129578… |
| DFAM3 | 93 | 82 | 50 | CDC42, IFT27, LOC101128843, LOC101133567, LOC101142457… |
| DFAM4 | 74 | 19 | 8 | BFSP2, DES, GFAP, GHAA, INA… |
| DFAM5 | 64 | 31 | 2 | ACR, CELA1, CFD, CTRL, F10… |
| DFAM6 | 47 | 12 | 1 | LOC101124648, LOC101136188, LOC101144932, LOC101146296, LOC101147007… |
| DFAM7 | 47 | 0 | 0 | LHB, LOC101123748, LOC101124600, LOC101154318, LOC109024055… |
| DFAM8 | 46 | 42 | 31 | CDK1, CDK10, CDK14, CDK15, CDK16… |
| DFAM9 | 45 | 26 | 13 | ACKR2, ACKR4, AGTR1, APLNR, CCR10… |
| DFAM10 | 44 | 0 | 0 | LOC101123968, LOC101128508, LOC115932853, LOC115932855, LOC129523578… |
| DFAM11 | 43 | 17 | 0 | LOC101126783, LOC101127524, LOC129524344, LOC129524346, LOC129524347… |

## Honest scope
- DNA tier = protein clustering of annotated CODING genes → catches ancient coding families +
  currently-silent coding copies. Non-coding/pseudogene + UNANNOTATED copies need a genomic
  self-alignment pass (next extension).
- 'transcribed' uses the REAL IsoSeq (not the ideal-coverage synthetic) — so silent = silent in
  this sample; ideal-coverage GGO would test detectability, not biology.
- copy-resolvability (which transcribed copies are distinguishable per-read) follows from the
  identifiability theorem: dispersed copies resolve by locus; co-located need PSVs (rare here).


---

## compara_fetch

# Ensembl Compara paralogy fetch — coverage report

Non-circular validation prep. This phase ONLY fetches/caches Ensembl Compara paralogy for the NAMED universe genes; the comparison vs our family grouping is the next phase.

- Cache: `bench/compara_cache.json` (keyed by `endpoint|symbol`, resumable; reruns fetch only missing).

- Relation summary: `bench/compara_paralog_relation.json`

- Source universe: `bench/copy_recovery_eval/results/universe.tsv`


## Universe gene inventory

| metric | count |
| --- | --- |
| total distinct gene_ids | 195 |
| named genes (not `^LOC[0-9]+`) | 41 |
| LOC-only genes | 154 |
| families | 62 |

## Mapping coverage (named genes)

| metric | count | of named |
| --- | --- | --- |
| got ENSG id | 40 | 41 |
| got paralogue data (HTTP ok) | 37 | 41 |
| non-empty paralogue set | 32 | 41 |
| unmapped (symbol not in Ensembl) | 0 | 41 |
| persistent fetch errors | 5 | 41 |

## Family-level checkability

| metric | count | of 62 families |
| --- | --- | --- |
| families with >=2 NAMED genes | 11 | 62 |
| families with >=2 MAPPED genes (checkable within-family pair) | 10 | 62 |

Within-universe symmetric paralog pairs (both genes mapped, Compara-linked): **5**


### Symbols with persistent fetch errors (rerun to retry)

`CREB1`, `GP1BB`, `NCAPH2`, `RABL2B`, `USP18`


Deterministic given the cache. Fetched 82 new API responses this run.


---

## compara_validation

# Non-circular validation: Ensembl Compara vs our paralog families
_Deterministic (no RNG). Minimizer k=15, w=10, canonical blake2b-64. Generated by `bench/compara_validation.py`. Compara release: `vertebrates`, `https://rest.ensembl.org`._
## Headline (the non-circular number)
> On the **10** universe families (12 within-family named gene pairs) that map to Ensembl, human Compara independently confirms **33%** (4/12) of within-family pairs as paralogs.
This **de-circularizes** the earlier `family_detection_validation.py`, which compared one sequence-similarity clustering (minimizer-Jaccard) against another (the minimap2-built `universe.tsv`) -- both reward shared subsequence, so that agreement was partly methodological. Here the truth set is Ensembl Compara paralogy: **gene-tree + species-tree reconciliation**, which never sees our minimap2 alignments. A Compara 'yes' is therefore an independent biological confirmation.
## KEY FRAMING: Compara is COARSER, so PRECISION is the metric
Ensembl Compara paralogy is gene-tree-based and far **coarser** than our families. Querying Compara for `RABL2A` returns the entire RAB / small-GTPase **superfamily** (68 paralogues genome-wide in this dump), whereas our family is the 2-copy tandem `RABL2A`/`RABL2B`. So **our families should be SUBSETS of Compara groups**: we split ancient superfamilies into recent copy-clusters, which is correct granularity, not error. The right question is therefore PRECISION -- *of the gene pairs we put in the same family, what fraction does Compara independently also call paralogs?* -- not recall. Compara lumping ancient paralogs we correctly separate is expected granularity and is **intentionally not scored as a miss**.
## (4) Coverage (stated up front -- small, honest sample)
- Universe: **195 genes / 62 families** (`universe.tsv`). Of these, **154 genes are LOC-only** (no human symbol) and cannot be mapped to Ensembl -- they are out of scope here.
- Named genes attempted: **41**. Got an ENSG id: **40/41** (1 lookup error: `GP1BB`). Non-empty paralogue set returned: **32** (4 paralogue-fetch errors: `CREB1`, `NCAPH2`, `RABL2B`, `USP18` returned an error and so contribute NO confirmations even if real).
- **Evaluable set: 10 universe families with >=2 mapped named genes -> 12 checkable within-family pairs.** This is the entire basis of the headline; it is a SMALL sample and is reported as exact counts, not extrapolated.
## (1) PRECISION vs Compara -- UNIVERSE families (the headline)
Of the **12** within-universe-family named pairs (both genes mapped to an ENSG), Compara independently confirms **4** as paralogs => **precision = 33%** (4/12).

| universe family | pair | Compara confirms? |
|---|---|---|
| `APOBEC3D` | APOBEC3D <-> APOBEC3F | **yes** |
| `ASDURF` | ASDURF <-> ASNSD1 | no |
| `CASP8` | CASP8 <-> FLACC1 | no |
| `CCDC188` | CCDC188 <-> ZDHHC8 | no |
| `CDPF1` | CDPF1 <-> PPARA | no |
| `CREB1` | CREB1 <-> METTL21A | no |
| `GCA` | GCA <-> KCNH7 | no |
| `GPR39` | GPR39 <-> LYPD1 | no |
| `LOC134758217` | RFPL1 <-> RFPL2 | **yes** |
| `LOC134758217` | RFPL1 <-> RFPL3 | **yes** |
| `LOC134758217` | RFPL2 <-> RFPL3 | **yes** |
| `RABL2A` | RABL2A <-> RABL2B | no |

**Reading the non-confirmations.** Most of the `no` rows are pairs where the universe's minimap2 step grouped two DIFFERENTLY-NAMED genes (e.g. `GCA<->KCNH7`, `CDPF1<->PPARA`, `CREB1<->METTL21A`) that Compara's gene tree does NOT place in the same paralog group. A non-confirmation can mean EITHER (a) a genuine universe over-grouping (sequence similarity from a shared domain / repeat, not orthologous-block paralogy), OR (b) a Compara gap / coarser-but-different granularity, OR (c) one partner's paralogue fetch errored (`CREB1`, `RABL2B` did -- so `CREB1<->METTL21A` and `RABL2A<->RABL2B` could not be confirmed even if true). The confirmed pairs are the unambiguous tandem/recent-duplicate families: `APOBEC3D<->APOBEC3F`, `RFPL1<->RFPL2`, `RFPL1<->RFPL3`, `RFPL2<->RFPL3`.
## (2) PRECISION vs Compara -- RUSTLE's minimizer-Jaccard grouping
Re-using the EXACT criterion from `family_detection_validation.py` (canonical k=15/w=10 minimizers, blake2b-64, union-find) over the same mapped named genes:

| grouping | within-cluster pairs | Compara-confirmed | precision |
|---|---|---|---|
| universe families | 12 | 4 | 33% |
| rustle Jaccard @ 0.3 | 2 | 0 | 0% |
| rustle Jaccard @ 0.06 | 13 | 5 | 38% |

rustle's non-trivial predicted clusters (mapped named genes only):
- @ 0.3: {CREB1, METTL21A}; {RABL2A, RABL2B}
- @ 0.06: {APOBEC3D, APOBEC3F}; {ASDURF, ASNSD1}; {CASP8, FLACC1}; {CCDC188, ZDHHC8}; {CDPF1, PPARA}; {CREB1, METTL21A}; {GCA, KCNH7}; {GGT1, GGTLC2}; {GPR39, LYPD1}; {RABL2A, RABL2B}; {RFPL1, RFPL2, RFPL3}

The restriction to the 40 mapped NAMED genes means rustle's minimizer-Jaccard rarely puts two of them in one cluster at the strict 0.30 bar -- the named genes that share a universe family are often different-symbol genes whose whole-gene Jaccard is below 0.30 (the same broad within-family Jaccard distribution documented in the prior report). Where rustle DOES cluster a named pair, that pair is the relevant precision test; counts are reported exactly above and are small.
## (3) GRANULARITY (observation, not error)
- Non-trivial Compara groups containing >=2 of our mapped named genes: **3**.
- Compara groups that MERGE >1 distinct universe family (Compara coarser): **1**. This is EXPECTED granularity -- Compara lumps at the superfamily level what we split into copy-clusters -- and is reported as an observation, NOT an error.

| Compara group (our genes) | distinct universe families merged |
|---|---|
| APOBEC3D, APOBEC3F | 1: `APOBEC3D` |
| GGT1, GGTLC2 | 2: `GGT1`, `GGTLC2` |
| RFPL1, RFPL2, RFPL3 | 1: `LOC134758217` |

- REVERSE check (real over-merge by us): universe families whose mapped named genes span >1 Compara component: **8**.
  - `ASDURF` (ASDURF, ASNSD1) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `CASP8` (CASP8, FLACC1) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `CCDC188` (CCDC188, ZDHHC8) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `CDPF1` (CDPF1, PPARA) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `CREB1` (CREB1, METTL21A) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `GCA` (GCA, KCNH7) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `GPR39` (GPR39, LYPD1) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `RABL2A` (RABL2A, RABL2B) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
## Cross-check: JSON's stated within-universe paralog relation
The JSON records **5** within-universe symmetric Compara paralog pairs (both genes mapped + linked): `APOBEC3D<->APOBEC3F`, `GGT1<->GGTLC2`, `RFPL1<->RFPL2`, `RFPL1<->RFPL3`, `RFPL2<->RFPL3`. Note `GGT1<->GGTLC2` is in this list but the two genes sit in SEPARATE universe families (`GGT1` and `GGTLC2`), so Compara confirms a paralogy that our minimap2 universe **split** -- a recall observation (granularity), not a precision failure, and consistent with the framing that we cluster more finely than Compara.
## Honest caveats
- **Human-as-proxy for gorilla.** The sequences are gorilla genes; Compara paralogy is queried on the HUMAN ortholog's ENSG. Paralogy is deeply conserved across great apes (these duplications predate the ~10-Mya gorilla/human split), so human paralogy is a strong proxy -- but it is a proxy, not the gorilla gene tree.
- **Compara is COARSER (superfamily-level).** A NON-confirmation can mean a genuine universe false-grouping OR a Compara gap / different granularity. Every non-confirmed pair is named above so the reader can judge; we do not silently treat them all as errors.
- **Only NAMED genes mapped.** 154 of 195 universe genes are LOC-only (no symbol -> no ENSG) and are out of scope. This validates the NAMED-gene backbone of the families, not all 62 families.
- **Symbol -> Ensembl mapping noise.** Genes were resolved by symbol via the Ensembl REST lookup; 40/41 got an ENSG, with 1 lookup error and 4 paralogue-fetch errors. An errored fetch yields an EMPTY paralogue set, which can only DEPRESS the confirmation rate (it cannot create a false confirmation) -- so the reported precision is a conservative lower bound for the affected pairs.
- **Recall vs Compara is intentionally NOT the metric.** Compara lumping ancient paralogs we correctly split is expected granularity, not error (see the framing section). We report granularity as an observation only.
- **Small N.** Only 10 families / 12 pairs are checkable. The headline is an exact count on a small sample, not a generalized rate; treat it as a directional, independent sanity check rather than a population estimate.
- **Determinism.** No RNG anywhere (the full mapped set is evaluated; no sampling). Minimizer hash is blake2b. Output is byte-stable.


---

## transcript_validation

# De-novo transcript validation (intron-chain Sn/Pr vs RefSeq)

Realigned all 101,467 de-novo transcripts (`denovo_transcripts.fa`) to the genome with minimap2
(`-ax splice:hq -uf`, low-mem `-I 1G --split-prefix`, MALLOC_ARENA_MAX=2 — the ulimit -v VIRTUAL cap
false-triggered 3×, drop it), `bam_to_gtf.py` → GTF, `gffcompare -r GGO_genomic.gff`. 100% mapped.

## Verdict: REAL + defensible
- **Intron-level precision 86%** (`-R -Q` 86.1) — posited splice junctions are genuine annotated sites.
- **Class codes: 98.9% overlap a KNOWN gene** (=20.7% FSM, c=31.0% ISM, j=30.9% novel-iso-of-known,
  m/n/k=16.4% retained/containment); only **0.5% `u` intergenic-novel** (artifact-suspect), ~0.6% antisense/intronic.
- **Sensitivity (where it looks, -R -Q): 76.6% of introns, 76.7% of loci recovered.** Genome-wide intron Sn 53%.
- Intron-CHAIN FSM only ~21% — NOT artifacts: long-read novel/partial isoforms (one novel junction in a
  ~9-exon chain fails whole-chain match). Expected annotation-incompleteness + novel-isoform discovery.
- CAVEAT: 31% ISM + 13% retained-intron => a real fraction are PARTIAL/incomplete (5' degradation / pre-mRNA),
  typical of read-derived assembly.

Artifacts in /home/juanfra/winloci_scratch/validate/ (dn_realigned.bam, dn_gw*.stats/.tmap).

