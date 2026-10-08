# Pre-registration: RG3 (regroup after polish) on the NPIP block of a held-out human library, with ape no-change controls

**Written 2026-09-28 (KEY=prereg) before any human_testis NPIP-block number exists.** This file binds once Amendment 1
records (1) the user's explicit acceptance of this file, of the substrate reuse in §11 and of the cost in §12, and
(2) this file's sha1. Until then nothing held-out runs.

**User decisions (2026-09-28, relayed by the orchestrating session; Amendment 1 records the user's own confirmation):**
1. **Ship RG3**: regroup after polish, splitting a `gene_id` only where its surviving transcripts share no same-strand
   exonic base. The rule is realised by the frozen `rg3.py` ec17e5408dedc82e0d7201e62f77bb181a45ff95 (§1). The Rust
   port must equal it byte for byte (gate G2, §10); the port never changes the rule.
2. **Test it held-out on the NPIP BLOCK, not genome-wide.** Substrate: **human_testis** (a different library from the
   development library A119b; same genome T2T-CHM13 v2.0, same annotation, same truths — §11 states what that means).
   Universe U = NPIP's loci plus every family linked to them by alignment, defined from BASE's own PAF and families,
   never from gene names (§2). BASE, RG3 and NULL_RG3 all run the families stage on that same U (§3).
3. **Ape NPIP contigs as no-change controls:** gorilla OR6737 and KB3781, chimp PTR, orangutan PPY (§5.5). Prediction:
   no change in NPIP-copy placement.
4. RG3 is opt-in. A default flip is the user's call whatever the outcome.

**What this author read.** `rg3.py`, `rg3_null.py`, `test_rg3.py`, `emu.py` (all frozen, sha1s in §10); the parked
`docs/archive/2026-09/PREREG_locus_representatives_2026-09-26.md` (untracked, STATUS PARKED); the scratchpad notes `npf_audit.md`,
`npf_variants.md`, `npf_critique.md`, `rt4_ghost.md`, `rep3_RG.md`; the memory files named in the task; the frozen
binaries' `SHA1SUMS`; the *file-name listing* of every `rt_arms/<sample>/` directory (names only, no size, no
content); the headers and formats of the **chimp_PTR** BASE families products (`loci.tsv`, `clusters.tsv`,
`copies.tsv`, `params.tsv`, PAF columns 1-12 of line 1); the truth tables' headers; the `--help` of the frozen
`mcl_families`, `family_score` and `copy_assign`; `tools/rustle_pipeline.sh stage_families`, `fam_call.sh`,
`tools/rlock.sh`; `figures/samples.tsv`. **No human_testis product was opened** (no GTF, families, PAF, cache, table
or log). **The NPIP block of human_testis has never been read by any author.**

**What was seen before this file, and where it came from** (design evidence; all of it is on development data
except the genome-wide testis rows listed in §11):
- Human A119b chr16 (dev): 26 Dishuck copies, 11 in fused loci, 22 in the NPIP family under BASE; RG3 recovers NPIPB2
  and NPIPB6 (both ghost links: the polish dropped the only bridging transcript), 22 → 24 of 26; Compara true pairs
  66 → 87 with 0 false; protein referee 100 → 125 with 0 false; Soto +8 and NPIP-U2 +12 "false" pairs, all involving
  NPIPB2/NPIPB6 and explained by a truth-filing disagreement plus pre-existing partner contamination (npf_critique
  §6); RG3 beats its matched null on 5/5 seeds (NULL 22 copies, 66/66 Compara); split correctness on dev 16 > 2,
  5 > 3, 4 > 0 (chr16, chr20, gorilla NC_073244.2). Register rows 1125-1127, 1134; `rep3_RG.md`, `npf_variants.md`.
- The truths disagree on the annotated readthrough copies NPIPA1, PDXDC2P-NPIPB14P and PKD1P6-NPIPP1: Soto and the
  union truth file them under PKD1 (Soto ID_149), Dishuck under NPIP (npf_critique §5.2). §4.3 handles this.
- Ape NPIP contigs were described in `npf_audit` (2026-09-28): across 108 holder-locus entries, no ape NPIP holder has
  an RG3 split or even a minority-bridge candidate. The ape NPIP contigs are therefore **no longer blind**; they serve
  only as no-change controls, and prediction P7 (§7) is a check of a known census, not a blind prediction.
- RG3's selection history (four adjacencies measured on dev after review2 X2; exon overlap named before measuring and
  kept; `rg.py`'s exact-junction rule fails split correctness 3/3) is in the parked prereg §2.1 and `rep3_RG` §1.

## 0. What this test asks, in one paragraph

On dev, the only information one locus representative loses at the family level is *other genes hidden in the same
node* (parked prereg §1, r1121-r1127), and on the NPIP block the hidden genes were two NPIP copies whose fused loci
had a partner's transcript as representative. RG3 makes those copies their own nodes without any constant. This test
asks whether that holds on a **different library** at the **same block**: does RG3, applied to human_testis's BASE
GTF, place at least as many Dishuck NPIP copies in the NPIP family as BASE and more than a matched random regroup,
without adding false pairs on the annotation-independent truths, without fragmenting genes, and without losing a
member — while leaving the ape NPIP blocks untouched. It is a **block-level, one-library** test; a genome-wide
verdict on RG3 (the parked design, ~18-24 h) is not attempted.

## 1. The rule RG3 (binding; realised by `rg3.py` ec17e5408dedc82e0d7201e62f77bb181a45ff95)

### 1.1 The rule, verbatim from the frozen file's docstring

```
RULE (the only rule the CLI runs):
  Input: an assembled GTF (the frozen BASE GTF of a sample). A transcript = a `transcript` line (GTF order = its index
  i); its strand = column 7 of that line; its exons = the `exon` lines with its transcript_id (1-based closed, sorted);
  its reads = the `reads` attribute (absent -> 0); its span = end - start + 1 of its transcript line.
  Pieces: the connected components of the transcripts of ONE gene_id, two transcripts adjacent iff they are on the
  same chrom and the same strand and share >= 1 exonic base (exons [a1,b1], [a2,b2] with max(a1,a2) <= min(b1,b2)).
  (A shared junction implies a shared exonic base, the donor base, so on one strand this adjacency CONTAINS rg.py's;
  a single-exon transcript joins any transcript it overlaps.)
  Piece representative: max (reads, span, -i) over the piece (earliest index on ties).
  Naming: a gene_id that is ONE piece is untouched. A gene_id split into m >= 2 pieces: the piece whose representative
  is max (reads, span, -i) keeps the gene_id; the other m - 1 pieces are renamed "<gene_id>.rg<k>", k = 2..m in order
  of their representative's index. (".rg<k>" and not rg.py's ".<k>": review2 L1, ".<k>" equals existing transcript_ids.)
  Output: the input GTF with `gene_id "<old>"` (its first occurrence) replaced by `gene_id "<new>"` on every line of a
  renamed transcript. No line is added, removed or reordered; nothing else changes.

ASSERTED on every run (exit 3 on failure): (1) the output equals the input once gene_id is removed from every line;
(2) every output gene_id holds exactly one piece; (3) no new name equals an input gene_id or an input transcript_id;
(4) every output gene_id maps to exactly one input gene_id (split-only); (5) no two output gene_ids of one input
gene_id share an exonic base on one strand (the pieces are the finest exon-disjoint split).
```

### 1.2 The three functions that are the rule (verbatim, except that the dev-comparator branches of `pieces` — `junction`, `donor`, `acceptor`, never run by the CLI — and the assert/stats bookkeeping of `regroup` are elided as `...`; the port reproduces these; `parse` and `rewrite` are stated below)

```python
def _join_exon_overlap(uf, txs):
    """same (gene_id, chrom, strand), >= 1 shared exonic base: sweep by start, join with the max-end owner so far."""
    groups = collections.defaultdict(list)
    for i, t in enumerate(txs):
        for a, b in t.exons:
            groups[(t.gene, t.chrom, t.strand)].append((a, b, i))
    for ex in groups.values():
        ex.sort()
        max_end, owner = None, None
        for a, b, i in ex:
            if max_end is not None and a <= max_end:
                uf.union(i, owner)
            if max_end is None or b > max_end:
                max_end, owner = b, i


def pieces(txs: list, adj: str = ADJ) -> list:          # ADJ = "exon"; the CLI never runs another adjacency
    """piece id (the representative's index) of every transcript, in txs order."""
    n = len(txs)
    uf = UF(n)
    ...  # adj == "exon": _join_exon_overlap(uf, txs)
    comp = collections.defaultdict(list)
    for i in range(n):
        comp[uf.find(i)].append(i)
    rep = [0] * n
    for mem in comp.values():
        r = max(mem, key=lambda i: (txs[i].reads, txs[i].span, -i))
        for m in mem:
            rep[m] = r
    return rep


def regroup(txs: list, adj: str = ADJ, suffix: str = ".rg") -> tuple[dict, dict]:
    rep = pieces(txs, adj)
    by_gene = collections.defaultdict(set)
    for i, t in enumerate(txs):
        by_gene[t.gene].add(rep[i])
    name = {}
    for g, reps in by_gene.items():
        if len(reps) == 1:
            name[(g, next(iter(reps)))] = g
            continue
        order = sorted(reps)                                   # ascending representative index
        best = max(order, key=lambda r: (txs[r].reads, txs[r].span, -r))
        name[(g, best)] = g
        k = 2
        for r in order:
            if r != best:
                name[(g, r)] = f"{g}{suffix}{k}"
                k += 1
    new = {t.tid: name[(t.gene, rep[i])] for i, t in enumerate(txs)}
    ...  # asserts 2, 3, 4, 5 (exit 3), stats
    return new, stats
```

**`parse` semantics** (what the port must match): lines starting with `#` are skipped for parsing but kept in the
output; lines with fewer than 9 tab fields are skipped; `attr(s, key)` returns the text between the FIRST `key "` and
the next `"`; a `transcript` line creates the transcript (index i = its order among `transcript` lines; `gene_id` from
its own attributes; `reads` = `int(v)` when `v.lstrip("-").isdigit()`, else 0; span = column 5 − column 4 + 1); a
duplicate `transcript_id` on a `transcript` line and an `exon` line before its `transcript` line each exit 3; exon
lists are sorted by (start, end). Any line with a `transcript_id` but a feature other than `transcript`/`exon` is
neither parsed nor rewritten unless its transcript was renamed (then its `gene_id` is rewritten like any other line
of that transcript). **`rewrite` semantics:** on every non-`#` line with ≥ 9 fields whose `transcript_id` was renamed
and whose `gene_id` differs from the new name, `f[8].replace('gene_id "<old>"', 'gene_id "<new>"', 1)`; the line's
trailing newline is preserved; assert 1 compares the two files line by line with `gene_id "…"` stripped and requires
equal line counts. **CLI:** `rg3.py IN.gtf OUT.gtf [--stats OUT.json]`, exit 0 / 2 (usage) / 3 (assert).

### 1.3 Properties (unchanged from the parked prereg §2.1)

- Threshold-free (≥ 1 shared base means "overlaps"), split-only, deterministic, per `gene_id` (so the result on a
  GTF restricted to a set of `gene_id`s equals the restriction of the result on the whole GTF, except that assert 3's
  name universe shrinks; gate G4 checks this on every substrate).
- It changes `gene_id` only: intron chains and every transcript are identical by construction.
- It is not r846's node cut (no transcript is cut), not r1013/r1017's coordinate split, and it merges nothing.

### 1.4 The Rust port ("ship") and its gate

- The port is written by the port agent from `rt4_ghost.md` "Proposed Rust change" with RG3's adjacency: an opt-in
  flag on `copy_assign --assemble-only` (proposed name `--gtf-regroup`, clap bool, default false; driver env
  `RUSTLE_GTF_REGROUP`, unset = byte-identical), applied to the polished GTF lines after the polish block and before
  `--gtf-tpm`, for any `--assembly-polish` mode, plus `params.tsv` rows only when set; and a **GTF-to-GTF path** so an
  existing GTF can be regrouped without re-assembly (exact CLI = the port agent's choice; Amendment 1 records it).
- **Gate G2 (§10).** Unset output cmp-identical to the frozen binary's; set output cmp-identical to `rg3.py`'s on the
  dev contigs, on the 19 fixtures ported as unit tests, and — after the freeze, light — on every substrate's BASE GTF
  through the GTF-to-GTF path. **The judged arm is `rg3.py`'s output.** A port mismatch is a port bug: fixed in the
  port, never in the rule, and it never blocks the verdict; the ship waits for G2.

## 2. The universe U (alignment-defined; per substrate; names never used)

Read only the STRUCTURE needed from the substrate's BASE products: `gene_id`s and exon coordinates (GTF), locus
keys and the fold map (`<s>.BASE.fam.loci.tsv`: `annotation → representative`, keys `chrom:start-end`), cluster
membership (`<s>.BASE.fam.clusters.tsv`: `cluster_id … chrom start end` of representatives), and PAF pairs
(`<s>.BASE.fam.loci.paf`, columns 1, 6, 10, 11). **Never NPIP membership counts, never a score.**

Notation. L = the BASE `gene_id`s with ≥ 1 exon line. K(g) = `chrom:min(exon start)-max(exon end)` (1-based) over
all exon lines of g, exactly `emu.read_gtf`'s locus key (= the `--from-gtf` node key; the PAF names loci by K). Two
`gene_id`s may share a key (chimp_PTR BASE has 3 such collisions); a key in U admits every `gene_id` with that key.
rep(g) = the `representative` column for K(g) in `loci.tsv` (K(g) itself if absent). Cl(g) = the `cluster_id` whose
row (chrom, start, end) equals rep(g), or ∅.

C = the substrate species' NPIP copy set (§4.1, §4.5), each copy with (chrom, start, end, strand, exons).

1. **H0** = { g ∈ L : chrom(K(g)) = chrom(c) and [start(K(g)), end(K(g))] ∩ [start(c), end(c)] ≠ ∅ for some c ∈ C }
   (span overlap, ≥ 1 base, strand-blind: inclusive on purpose; the strict same-strand exon-overlap holder of §5.1
   is a subset).
2. **F0** = { Cl(g) : g ∈ H0, Cl(g) ≠ ∅ }; **U1** = H0 ∪ { g ∈ L : Cl(g) ∈ F0 } (every locus whose representative is a
   member of a family that holds an NPIP-overlapping locus, folded loci included).
3. **N(U1)** = { g ∈ L : some PAF record has {query, target} = {K(g), K(h)} with h ∈ U1, query ≠ target, column 11
   (alignment block length) ≥ 300 and column 10 / column 11 ≥ 0.70 } — the record filters `mcl_families` applies
   before any exon clause (`emu.read_paf`: MIN_BP 300, MIN_ID 0.70), so "alignment" means a record the graph step
   would consider, admitted or not.
4. **U = U1 ∪ N(U1). Closed once:** N(U) is *not* added. Reported: |H0|, |F0|, |U1|, |N(U1) \ U1|, |U|, and U's contig
   composition.

**Orangutan.** Its genome-wide BASE families never finished (`orangutan_PPY.BASE.fam.clusters.tsv` does not exist; a
`loci.paf` is on disk, but the parked prereg's ops note records the run looping on shard 78/100, so that PAF is not
trusted as complete and is not used). Its BASE reference is the contig-restricted families run over the 4
NPIP-carrying contigs (`dev_small_frozen/npf_audit/ppy/fam/ppy.{clusters,loci}.tsv`, 5,608 loci, the same binary and
flags). Its PAF was not kept, so the run agent re-runs that same command once (≈ 6 min, heavy); the re-run must
reproduce `ppy.clusters.tsv` byte for byte (gate G8), and U_PPY is built from the re-run's PAF and clusters.

**Why this U.** It contains every locus that can exchange an edge with an NPIP-overlapping locus or with its family,
so an RG3 piece from any U locus can reach NPIP, and any locus NPIP's family could lose to a piece is present. It is
built from BASE only, so RG3's new pieces are added to it (§3) and nothing in it depends on any arm's result.

## 3. The arms on U (all three on the same U; frozen `mcl_families` 91ef2e1cdf8f74fe878912dfc132253faba553b5)

| arm | GTF_U | how |
|---|---|---|
| **BASE_U** | every BASE GTF line whose `gene_id` ∈ U, original order (`#` lines dropped) | restriction only |
| **RG3_U** | `rg3.py <s>.BASE.gtf` (the whole GTF), then the lines whose `gene_id` ∈ U or equals `<g>.rg<k>` with g ∈ U | RG3 on the whole GTF (assert 3 sees every name), restricted; **gate G4:** equals `rg3.py GTF_U(BASE)` byte for byte |
| **NULL_RG3_U(k)**, k = 1..5 | `rg3_null.py GTF_U(BASE) OUT --label <s>.U<k>` | the §9 protocol on the U-restricted GTF, so its donors are U loci and it matches RG3's U splits; 5 fixed labels ⇒ 5 fixed draws |

The families stage, identical for every arm (the driver's `stage_families` command verbatim, `--emit-units` included,
stock minimap2 2.30 on PATH, `RUSTLE_MINIMAP2` unset, `TMPDIR` under the scratch dir), under `bash tools/rlock.sh heavy`:

```
mcl_families --from-gtf <GTF_U> --fasta <species fasta> --threads 4 --min-exonic-bp 1 --min-shared-exon-frac 0.60 \
  --emit-units --out <scratch>/<s>.<arm>.fam > <scratch>/<s>.<arm>.families.log 2>&1
```

Human fasta `/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa`; gorilla `…/_from_wsl/winloci_scratch/GGO.fasta`;
chimp `…/winloci_data/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna`; orangutan
`…/winloci_data/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna` (the registry values, `figures/inputs.local.tsv`).

- **What a U run is and is not.** The all-vs-all on U is a *new* alignment (r1126: minimap2's output depends on the
  target set; `-N 50` caps secondaries per query, so a U run may keep records the genome-wide run dropped). That is
  why all three arms are run on U and compared with each other, never with the genome-wide BASE. Reported, not
  judged: the concordance of BASE_U with the genome-wide BASE restricted to U (loci whose cluster co-membership
  changes), and c1 (§5.1) under the genome-wide BASE clusters.
- **Identity gates per U run** (a failure is a bug; it stops that substrate, whose clauses become "not measured"):
  `rg3.py` asserts 1-5; `emu.py` fac9a560 as R0 reproduces the run's `clusters.tsv` byte for byte from its own PAF
  (`run_emu3.py` b6ab23d8 pattern); 10 random node relabellings (`relabel_null.py` 8c050970 pattern, seeds
  20260926+k) give BASE_U's and RG3_U's partitions (else MCL tie order could confound a ±1 result: C1 and C2 are then
  not measured on s); every BASE `gene_id` ∈ U appears in GTF_U(RG3) with the same transcript count (the keeper keeps
  its name).
- Apes: BASE_U, RG3_U and NULL_RG3_U(1) only (the null is reported there, never judged).

## 4. Truths

### 4.1 NPIP copies, human: Dishuck 2025 via `docs/lit_subclusters_npip_dishuck_check.tsv` (sha1 6b6a1027…)

The **26 chr16 rows** (NPIPB1P on chr18 is reported as a 27th row and never counted). A copy's coordinates are the
table's; its exons are `audit.human_truth()` (`audit.py` fd0305ca…, npf_audit): the copy's RefSeq record exons from
`families_gw/species/human/genes.tsv`, with two stated edits — PKD1P6-NPIPP1 is its NPIPP1 half (exons at or above
15,126,650) and NPIPB14P, which has no exons of its own, takes the PDXDC2P-NPIPB14P readthrough exons inside its
span. Dishuck's level-2 groups are Iso-Seq expression groups and are **never scored** (memory
`project_tbc1d3_subclusters_population_unit`).

### 4.2 Pair truths, human (frozen `family_score` 7723029bb3af134d…, `--pairwise --per-family`)

| truth | file (sha1 of the families table) | `family_score` call | role |
|---|---|---|---|
| **Ensembl Compara Primates paralogues** (genome-wide, `bench/truth.py compara`) | `families_gw/species/human/compara.Primates.families.tsv` (`Gene Name / Family ID / Contig`) | `--gff families_gw/species/human/genes_only.gff --chrom ALL` | **judged** |
| **Protein-homology families, chr16** (the §6ko rule; the genome-wide table was never built — queue 2 of the parked prereg never ran, and no `ph.families.tsv` exists for any species) | `fig7/current/human_chr16_ref.families.tsv` ebc7227c…, gff `fig7/current/gff/human_chr16.gff` 1bd6fafa… | `--chrom chr16` | **judged** (U loci off chr16 are unjudgeable here; their count is reported) |
| Soto 2025 families (a COVER; house rule "not independent") | `families_gw/species/human/soto.families.tsv` | `--chrom ALL` | reported, both filings (§4.3) |
| NPIP union truth U2 (register 990; shares O1's aligner) | `families_gw/species/human/npip_u2.families.tsv` = `npip_union_truth` | `--chrom chr16` | reported, both filings |

Pair definitions are `family_score`'s: a locus resolves to the annotated gene with the largest **span** overlap
(`gene_at`, strand-blind); a predicted pair is two members of one cluster whose genes are both in the truth
universe; TP iff both genes share a truth family; **false** iff both are in the universe and share none. The
universe intersection deletes pairs with an unlabelled gene (r770/r991): every arm's unjudgeable predicted pairs are
reported beside the counts. `gene_at` names 9 of 21 NPIP holders after a partner on dev (npf_audit §2): this bias is
identical across the arms, and it is why the primary clause C1 is per copy and exon-level, not pair-level.

**Pair sets.** "New false pairs" needs sets, not counts. Instrument: `rg_score.pair_sets` (`rep2_RG/rg_score.py`
fea2167b…) over the U run's clusters; **gate G6:** its (predicted, TP) counts equal `family_score`'s on BASE_U for each
judged truth. If G6 fails for a truth, the count form fp = predicted − TP is used for that truth and the report says
so (the count form allows a swap, which is weaker; disclosed).

### 4.3 The PKD1-vs-NPIP filing disagreement: report under both filings

Soto and U2 file NPIPA1, PDXDC2P-NPIPB14P and PKD1P6-NPIPP1 under PKD1 (Soto ID_149); U2's ID_149 also holds
PKD1P3-NPIPA1, PKD1P5-LOC105376752 and LOC131696449 (the label of the NPIPA6+A7 unit); Dishuck calls the copies these
records hold NPIP copies. Two filings, applied to **every** truth table at scoring time:
- **Filing A (as published):** the table unchanged.
- **Filing B (Dishuck):** for each of the six names above that occurs in the table, one row is ADDED giving it the
  truth's NPIP family id, where the truth's NPIP family = the `Family ID` holding the most of the 26 copies' RefSeq
  symbols. Existing rows are kept (a cover-consistent minimal edit; `family_score` accepts repeated genes). Where none
  of the six names occurs, B = A (expected for Compara and the protein referee, whose universes exclude
  pseudogene/readthrough records; if not, both filings are still computed).
- The **judged** truths must pass C2 under **both** filings. Soto and U2 are reported under both, with every gained
  false pair classified as (i) filing (both genes NPIP copies by Dishuck), (ii) an RG3 piece paired with a
  partner-labelled unit that was ALREADY in the NPIP family under BASE_U (pre-existing contamination multiplied), or
  (iii) other. On dev, all +8 / +12 were (i) or (ii).

### 4.4 Split-correctness truth: the species genes table (human, gorilla only)

`families_gw/species/{human,gorilla}/genes.tsv` + `genes_only.gff` (record `Name|start1`, strand from the GFF), as
`cls.py` 413390944e4a7e54213256ea649a31802283c135 reads them. No chimp or orangutan table exists (queue 2 never ran):
C3 is not computable there, and the apes need it only as a reported row.

### 4.5 NPIP copies, apes (controls)

- Gorilla OR6737 and KB3781: T_member, the 2026-09-17 proxy (liftoff AND an independent identity landing agree on
  the same human copy; 25 copies), as `audit.ape_member_truth('GGO')` builds it from
  `/mnt/linuxdisk/home/juanfraitu/ggo_npip/labels/GGO.truth.tsv`.
- Chimp PTR: T_member (19; `ape_member_truth('PTR')`) and T_native (46 records named NPIP;
  `/mnt/linuxdisk/home/juanfraitu/npip_membership/out/PTR.T_native.tsv`, `audit.ptr_native_truth`).
- Orangutan PPY: the landing truth (19) and the native truth (20) built by `ppy_truth.py` 57be8160… in npf_audit
  (`ppy/ppy_truth.{member,native}.tsv`), as `audit.ppy_truth(kind)` reads them.
- Each ape truth set gives its own U (§2, with C = that set) and its own control row; a species with two truths has two
  rows, both required for the control (§5.5).

## 5. Metrics and clauses (integer arithmetic; per substrate; nothing pooled across species)

### 5.1 Per-copy placement (the primary instrument; npf_audit's, frozen as SHA1_NPF in Amendment 1)

For an arm A on U with units = the `gene_id`s of GTF_U(A), families from the run's `clusters.tsv` through its fold map
(`loci.tsv`):
- **holder_A(c)** = the unit u maximising the exonic overlap (bases) between the exon union of u's transcripts on
  strand(c) and c's exons; ≥ 1 base required, else c is **absent**. Ties: the unit with the smaller K start, then the
  lexicographically smaller `gene_id`.
- **present(c)** iff holder_BASE(c) exists. **n_present** = |{ c : present(c) }| (fixed by BASE; a split-only rewrite
  cannot remove an overlap, so RG3's holders exist for the same copies — checked).
- **fam_A(u)** = Cl(rep(u)) in the U run (∅ = unclustered).
- **NPIP_A** = the `cluster_id` with the most present copies placed in it (fam_A(holder_A(c)) = it); ties: more members,
  then the lower MCL index.
- **placed_A(c)** ∈ { NPIP, other cluster, singleton (fam = ∅), absent }; **c1(A)** = |{ c present : placed_A(c) = NPIP }|.
- Reported per copy: every same-strand touching unit (≥ 1 exonic base) and its cluster (npf_critique §4: "copies in
  NPIP" is a max-overlap-holder statement; the touching-set table shows what it hides); the number of clusters
  holding ≥ 1 copy and the copies per cluster; NPIP unit precision (units of NPIP_A whose exons overlap a copy /
  units); partners pulled in (NPIP_A units whose exons hit a non-NPIP annotated gene, `npip_units.py` 6dafdd7d…);
  for each fused holder its RG3 class (split / one piece) and, for description only, the ghost / minority /
  dominant bridge class of npf_audit (which needs annotated sides and is therefore not a rule).
- **Gate G5:** SHA1_NPF reproduces npf_audit's dev table (`npf_audit/out/hsa16.dishuck.copies.tsv`: BASE 22 of 26 in
  NPIP; NPIPB2 singleton, NPIPB5 in the SMG1P family, NPIPB6 in the EIF3C family, A4 in a fragment cluster) and
  npf_variants' V1 = 24 and NULL_V1 = 22 rows (`npf_variants/out/hsa16.<arm>.score.json`), from the frozen dev families
  products, before any held-out product is opened.

### 5.2 The clauses (human_testis; zero tolerance; no new constant)

| clause | passes iff | notes |
|---|---|---|
| **C1** NPIP copies placed | (i) for every present c: placed_BASE(c) = NPIP ⇒ placed_RG3(c) = NPIP; (ii) c1(RG3) ≥ c1(BASE); (iii) c1(RG3) > max over k = 1..5 of c1(NULL_k) | (i) is per copy, stronger than the count; (iii) protects against "more, smaller loci" (§9). If every NULL draw is a total shortfall (no eligible donor in U), NULL_k = BASE_U by construction and (iii) reads c1(RG3) > c1(BASE). **A no-effect outcome (c1(RG3) = c1(BASE) = max NULL) fails (iii)** and is handled by the verdict table as KEEP OPT-IN, never REFUTE |
| **C2** true pairs | on EACH judged truth (Compara, protein-homology chr16) under BOTH filings: tp(RG3_U) ≥ tp(BASE_U) **and** FP(RG3_U) \ FP(BASE_U) = ∅ (G6; else fp(RG3_U) ≤ fp(BASE_U)) | Soto and U2 reported under both filings with the (i)/(ii)/(iii) breakdown, never judged; the per-family sign (b = truth families whose TP rise, c = those that fall) reported |
| **C3** split correctness | over RG3's splits whose input `gene_id` ∈ U: SEP > FRAG (strict; a tie fails), classes exactly as `cls.py` 41339094 defines them against the human genes table (A(P) = annotated genes sharing ≥ 1 same-strand exonic base with a transcript of piece P; FRAG = some gene in A(P) for ≥ 2 pieces; SEP = no gene cut and every A(P) non-empty; UNJ = no gene cut and some A(P) empty, reported) | **judged iff SEP + FRAG ≥ 1**; n is reported and is expected to be small (dev chr16 had 19 splits on the whole contig). PURE, UNJ, the RefSeq-readthrough-dropped form, and NULL_k's SEP/FRAG are reported |
| **C4** no NPIP-block locus lost | (a) every u ∈ members(NPIP_BASE_U) has its keeper or ≥ 1 of its RG3 pieces in members(NPIP_RG3_U); (b) every u ∈ U clustered in BASE_U has its keeper or ≥ 1 piece clustered in RG3_U (any family) | (b) can fail by MCL re-flow through a changed keeper model (a 2-member family losing its edge) — that is the recall-side harm the genome-wide P2 would catch, so it is judged here at tolerance 0 |

### 5.3 What is reported for human_testis and never judged

|U| and its construction counts (§2); BASE_U vs genome-wide BASE concordance; c1 under the genome-wide clusters; the
RG3 statistics (`--stats`: gene_ids, splits, pieces added, transcripts relabelled) on the whole GTF and within U; the
NULL statistics (splits matched, bin exact/nearest/unmatched, shortfalls, pool size); largest family (L) and π = RG3
pieces among its members, with the H form L(RG3) ≤ L(BASE) + π stated as reported; **revealed copies** (truth-labelled
genes that are family members only through an RG3 piece and share a truth family with another member of their
family; counted per piece, so family size cannot multiply them; dev: {NPIPB2, NPIPB6}); gained TP pairs by carrier
(piece / keeper / unchanged); unjudgeable predicted pairs per truth; the bipartite F rows (scipy tie policy, r1045)
for continuity only.

### 5.4 Can every judged clause fail? (negative controls)

| clause | the failing outcome | shown on dev, or plausible held-out |
|---|---|---|
| C1 | a copy leaves NPIP through re-flow (i); no copy gained (iii, the no-effect case); a random regroup gains as much (iii) | dev: NULL_V1 22 = BASE 22, so (iii) held only because RG3 gained 2; testis at 1.9 isoforms per locus (seen, §11) may hold no ghost at NPIP — (iii) then fails and the verdict is KEEP OPT-IN |
| C2 | a piece joins a large family wrongly (FP pairs scale with family size: one wrong member of a 30-member family adds up to 29 false pairs), or a true pair is lost through re-flow | dev control from arm A: R1 (Compara pair precision 1.000 → .753, largest 27 → 48); Soto/U2 on dev show what a filing artefact looks like (+8/+12) |
| C3 | SEP ≤ FRAG | dev: RG3 on the R3 chr20 GTF ties 1 = 1; 5 of RG3's 33 dev splits are FRAG (a dropped bridge of one long gene, or an annotated readthrough record spanning both parents) |
| C4 | a keeper's changed exon model drops the edge that held a 2-member family; a NPIP member's every piece lands outside NPIP | not seen on dev (0 lost); plausible on the SD-rich block |
| A (apes) | any copy's placement changes | not seen on dev; plausible through re-flow from a non-holder U locus that RG3 splits |

### 5.5 Ape no-change control (clause A; per ape sample, per truth set)

For every present copy c: holder_RG3(c), mapped to its input `gene_id`, equals holder_BASE(c); and the member set of
its cluster, with pieces mapped to their input `gene_id`s, is identical between BASE_U and RG3_U. **A(s) passes iff
this holds for every present copy under every truth set of s.** Reported: RG3 splits among U loci (holders and
non-holders separately), c1 per arm, the NULL(1) row, |U|. Prediction P7 (§7). The parked census (npf_variants §8: 0
RG3 splits among 108 holder entries) is already known; the blind part is re-flow from non-holder splits.

## 6. Verdict (dev never enters; species never pooled; human_testis carries the verdict, the apes carry the control)

Common rules. **Not judged:** C3 with SEP + FRAG = 0. **Not measured:** a clause disabled by an identity gate (§3) or
uncomputable within the machine rules; it caps the verdict at KEEP OPT-IN unless REFUTE already holds. **Undecided:**
a BASE prerequisite missing (human_testis BASE families exist; the orangutan re-run of G8 failing makes only the
orangutan control not measured).

Let F = the set of failed clauses among {C1, C2, C3, C4, A}, where A fails iff it fails on ≥ 1 ape sample and truth set.
**Harm** = any of: C1 (i) fails; c1(RG3) < c1(BASE); a new false pair on a judged truth under either filing (or
fp rises under the count form); C4 fails.

| verdict | condition |
|---|---|
| **EFFECTIVE (block-level, one human library)** | F = ∅, with c1(RG3) − c1(BASE) ≥ 1 |
| **KEEP OPT-IN** | \|F\| = 1 and no harm (this includes the no-effect case: C1 (iii) fails with c1(RG3) = c1(BASE), everything else passes — reported as "no ghost at NPIP in this library"); or F = ∅ but a clause is not measured or C3 is not judged |
| **REFUTE** | harm; or \|F\| ≥ 2 |

An identity-gate failure on human_testis makes every clause not measured there (KEEP OPT-IN at most). The control A
failing alone is |F| = 1 → KEEP OPT-IN with the ape change described; A failing together with any human clause → REFUTE.

## 7. Predictions (set before any held-out number; probabilities are this author's)

1. `rg3.py` asserts 1-5 pass on 6/6 BASE GTFs (0.95). RG3 splits ≥ 1 U `gene_id` on human_testis (0.75); pieces add
   0.2-1% loci genome-wide, testis the fewest (post-exposure: testis has 1.9 isoforms per locus, seen).
2. **C1:** c1(RG3) − c1(BASE) ≥ +1 (0.35); ≥ +2 (0.15); = 0 (0.55); < 0 (0.05). NPIPB2's GSPT1 fusion is a ghost in
   testis too (0.30); NPIPB6's EIF3CL fusion (0.30); the two are correlated (one polish, one library). The ghost SET
   differs from A119b's (0.7): ghosts are library-dependent because the polish drops a bridge on read counts.
3. max_k c1(NULL_k) = c1(BASE) (0.85); ≥ 1 NULL draw is a shortfall or bin-unmatched in U (0.5).
4. **C2:** 0 new false pairs on Compara and on the chr16 protein referee under both filings (0.80); tp(RG3) ≥ tp(BASE)
   on both (0.90). Soto/U2 show ≥ 1 "false" pair under filing A iff C1 gains a copy (0.8), fewer under filing B (0.7).
5. **C3:** judged (0.7); passes given judged (0.7); FRAG ≥ 1 within U (0.35).
6. **C4:** (a) 0.95, (b) 0.90.
7. **Apes:** all four samples unchanged under every truth set (0.65; per sample 0.90); ≥ 1 ape has ≥ 1 RG3 split among
   non-holder U loci (0.6); 0 splits among holders (0.95, a known census).
8. **Port:** G2 byte-identical on dev at the first build (0.6), after ≤ 2 fixes (0.95); the GTF-to-GTF path
   cmp-identical to `rg3.py` on all 5 BASE GTFs (0.95 given dev passes).
9. **Verdict:** EFFECTIVE 0.25; KEEP OPT-IN 0.55 (mostly the no-effect case); REFUTE 0.12; not measured / undecided 0.08.
   Calibration: P(C1 gain) 0.35 × P(C2 | gain) 0.8 × P(C3 pass or unjudged) ≈ 0.8 × P(C4) 0.9 × P(A) 0.65 ≈ 0.13 under
   independence; the components are positively correlated (a clean ghost piece passes C2-C4 together), so 0.25.

## 8. What would falsify the design reasoning (reported whatever the verdict)

- **F1 "RG3's NPIP gain is a ghost piece joining its true family."** Measured iff c1(RG3) > c1(BASE). Falsified if any
  gained copy's holder_RG3 is a keeper or an unchanged unit (the gain would then be MCL re-flow, not the rule). Dev: 2/2
  gained copies are pieces.
- **F2 "Ghosts are a library property, not a locus property."** Reported: the fused NPIP holders of testis BASE_U
  classified by RG3 (split / one piece) against A119b's list (GSPT1~NPIPB2, NPIPB6~EIF3CL ghosts; NPIPB5 real bridge;
  four dominant readthroughs). Identical sets on two libraries would falsify it (and would make the polish's role the
  next question, §13).
- **F3 "A random regroup of the same size does nothing at NPIP."** Falsified if any NULL draw raises c1 above BASE.
- **F4 "RG3 cannot fragment a gene on the block."** Falsified by FRAG ≥ 1 within U; reported with the gene named.

## 9. The NULL: `rg3_null.py` 6ece9ae088766f8c22ecc251c40a883d1451c426 (protocol verbatim; 5 fixed labels)

`rg3_null.py` is `rg_null.py` 08b5c317 verbatim with RG3's pieces and the `.rg<k>` names. Its docstring:

```
RG splits the gene_ids whose transcripts form >= 2 junction components (rg.py). NULL_RG leaves those untouched and
instead splits the same number of gene_ids that are ONE junction component, moving the same numbers of transcripts
into the same numbers of new gene_ids, so that a family / locus effect of RG can be attributed to its predicate
(junction disconnection) rather than to "more, smaller loci".

PROTOCOL (per contig c; seed string "20260926:<label>:<c>", Python random.Random, version-2 string seeding, so the
draw does not depend on PYTHONHASHSEED; every list sorted before a draw):
  S_c  = the gene_ids RG splits on c, sorted by (first exon start, gene_id). For g in S_c, its non-keeper pieces in
         RG's naming order, each with n_p transcripts; T_g = sum n_p; b_g = max over those pieces of
         floor(log2(max(1, reads of the piece representative))).
  Pool = gene_ids on c that RG leaves as ONE piece and that hold >= 2 transcripts, sorted by (first exon start,
         gene_id). For a donor h: rep(h) = max (reads, span, -index) (rg.py's rule); cand(h) = h's SPLICED transcripts
         other than rep(h), in GTF order.
  For each g in S_c in order:
    bin search b = b_g, b_g - 1, b_g + 1, b_g - 2, b_g + 2, ... (0 <= b <= 40; the lower bin first on a tie of distance);
    eligible(b) = donors still in Pool with |cand(h)| >= T_g and >= 1 transcript of cand(h) in bin b;
    the first non-empty eligible(b) is used; if none exists at any bin: donors with |cand(h)| >= T_g (bin unmatched,
    counted); if none: g is a SHORTFALL (counted, not replaced).
    h = rng.choice(eligible); h leaves the Pool (one split per donor).
    seed = rng.choice(cand(h) transcripts in bin b)  (any cand(h) transcript when the bin was unmatched);
    rest = rng.sample(cand(h) minus seed, T_g - 1);  moved = [seed] + rest.
    Pieces are filled in order: piece 1 gets moved[0 : n_1], piece 2 the next n_2, ...; they are named "<h>.2",
    "<h>.3", ... (RG's naming form; asserted not to collide with an input gene_id).
  Output: the BASE GTF with gene_id rewritten on the moved transcripts' lines only (rg.rewrite, same identity checks).
(read "rg.py" as "rg3.py", "junction" as "exon-overlap", ".2" as ".rg2")
```

- Applied to GTF_U(BASE) with labels `<s>.U1` … `<s>.U5` (human_testis) and `<s>.U1` (apes), so the seed strings are
  fixed by this file: `"20260926:human_testis.U1:chr16"` etc. It must be run from a directory holding the frozen `rg3.py`
  (it does `import rg3 as rg`).
- Matched: split `gene_id`s per contig within U, pieces per split, transcripts per piece, the depth bin of the deepest
  new piece's seed, spliced. Not matched, reported: the donor's own depth; exon disjointness (NULL pieces overlap their
  keeper by construction — which is what makes the null a predicate control).
- Dev (npf_variants): NULL_V1 × 5 gives 22 copies on 5/5 seeds and BASE's 66/66 Compara pairs.

## 10. Gates and frozen sha1s (all before the freeze except G2(f), G3, G4, G8, which need held-out products)

| id | what must hold | who |
|---|---|---|
| **G1** | `rg3.py` ec17e5408dedc82e0d7201e62f77bb181a45ff95; `rg3_null.py` 6ece9ae088766f8c22ecc251c40a883d1451c426; `test_rg3.py` 03e10572b144bf8f05fc26048cfe9a871aefae3f (19 fixtures pass); `rg3.py` reproduces `dev_small_frozen/rep3_RG/rg/*.RG3.gtf` byte for byte and is idempotent; `rg3_null.py` reproduces `rep3_RG/null/*.NULL3.gtf` under two `PYTHONHASHSEED`s | run agent, dev |
| **G2** port | (a) unset: the port binary regenerates the 3 dev BASE GTFs (A119b chr16 `chr16:0-96330374`, chr20 `chr20:0-66210255`, gorilla OR6737 `NC_073244.2:0-80312928`; `--region` on the full BAM, seeded, driver flags, `--assembly-junctions strict`, the driver polish string, `--gtf-tpm`; `fj_impl/run.sh` a23cd82c with `EXE` = the port) cmp-identical to the frozen `fj_bin_frozen/copy_assign` b13b6ae6…'s regeneration, every product; (b) set: the GTF cmp-identical to `rg3.py` on the unset GTF, 3/3, and non-GTF products identical except the new `params.tsv` rows; (c) the 19 fixtures as Rust unit tests; (d) the GTF-to-GTF path cmp-identical to `rg3.py` on the 3 dev GTFs; (e) `cargo test --release` passes (log to a file); (f) after the freeze: the GTF-to-GTF path on the 5 substrates' BASE GTFs cmp-identical to `rg3.py` (light) | port agent (a-e), run agent (f) |
| **G3** | emu R0 (`emu.py` fac9a5602f223fa94d889f204d309b4765498520) byte-exact on every U run; the 10-relabelling check on BASE_U and RG3_U | run agent |
| **G4** | `rg3.py`(whole GTF) restricted to U = `rg3.py`(GTF_U) byte for byte, per substrate | run agent |
| **G5** | SHA1_NPF reproduces the dev per-copy tables (§5.1) | run agent, dev |
| **G6** | `pair_sets` counts = `family_score` counts on BASE_U per judged truth (else the count form, disclosed) | run agent |
| **G7** | `cls.py` 41339094 reproduces `rep3_RG/cls/cls.json` (RG3 16/2/1, 5/3/2, 4/0/0; old rule 15/38/31, 4/22/11, 3/4/1) | run agent, dev |
| **G8** | the orangutan contig-restricted BASE re-run reproduces `npf_audit/ppy/fam/ppy.clusters.tsv` byte for byte | run agent |

Frozen binaries `/mnt/linuxdisk/tmp/rustle_figures/fj_bin_frozen/` (built from `main@e0cf3282` + `source.diff`
652394a8…; the task states they equal HEAD 67f86286, and `git diff e0cf3282 HEAD -- src` touches the same four `src/`
files as `source.diff`, which additionally carries `tools/rustle_pipeline.sh`):
`mcl_families` 91ef2e1cdf8f74fe878912dfc132253faba553b5, `family_score` 7723029bb3af134da0d8bc5606afde743d7474b0,
`copy_assign` b13b6ae627134d771ca19cad7a05dd387080ef96, `as_table` 35e2a17a6f595e5658878993db4d09927870a491.
Other frozen inputs: `bench/truth.py` 2c59a2d8…, `bench/lib.py` e7765bc0…, `figures/_o1.py` e3994c68…,
`figures/inputs.local.tsv` e26e49fc…, `docs/lit_subclusters_npip_dishuck_check.tsv` 6b6a1027…, `npf.py` fdc1a7d3…,
`audit.py` fd0305ca…, `score_arm.py` 258d1886… (its `rustle_figures_dev/npf_*` paths now resolve under
`rustle_figures/dev_small_frozen/`; SHA1_NPF is the frozen copy with those paths fixed and a `--u` restriction
argument, nothing else), `attr.py` 5ff32b43… (attribution only; its pair counts are never quoted, npf_critique §6),
`npip_units.py` 6dafdd7d…, `relabel_null.py` 8c050970…, `run_emu3.py` b6ab23d8…, `ppy_truth.py` 57be8160….
Amendment 1 records SHA1_NPF, the U-builder SHA1_U, the scorer SHA1_SCORE, the port's commit and flag, and every
command line.

**Any failure is fixed in the instrument, never in the rule, the universe, a null, a clause or a floor.** A fixed
instrument re-passes every dev gate before it runs again.

## 11. Exposure (what "held-out" means here, plainly)

- **human_testis is held out in LIBRARY only.** It shares the genome (CHM13 v2.0), the annotation, the Dishuck copy
  table and every pair truth with the development library A119b, whose chr16 is where RG3 was selected and where
  its NPIP gain was found. What differs is the reads, hence the assembled loci, the fused-locus set, the polish's
  dropped bridges (the ghosts) and the representatives. A copy lost to a ghost in A119b need not be lost in testis
  and vice versa. This is disclosed as the weakest form of held-out short of the dev contig itself; the user chose
  it because the apes have no fused NPIP members to test (npf_audit §2, npf_variants §8) and no other human library
  exists on this machine. Any claim from this test is "on a second library at the same block".
- **Reuse count.** The six samples (A119b, testis, OR6737, KB3781, PTR, PPY) have served verdict-bearing held-out
  tests at least three times: the readthrough R prereg (09-25, r1117), the readthrough v2/v3 prereg (09-26, r1119/r1120;
  v2 §5 says a rule chosen after seeing held-out results needs a new library or the user's explicit acceptance of a
  further reuse), and the TES `pas-end` rule (09-27, r1136), with descriptive held-out reads in r1137-r1139. **This is
  at least the fourth verdict use of human_testis, the first on node regrouping and the first on the NPIP block.** The
  parked RG prereg would have been a third reuse and never ran. After this test, human_testis's NPIP block is spent for
  regrouping work.
- **Seen of human_testis before this file (by earlier authors, genome-wide, never NPIP-block):** the published
  `rt_arms/tables_v3` BASE rows (Compara pairwise TP 222 of 3,053, 239 predicted; Soto 73 / 1,826, 91 predicted; 338
  families, largest 47; `a.fused`, TES/TSS, `found_annotated`); 1.9 isoforms per locus with 11% of exon bp and 22% of
  junctions outside the representative (`rep_loss.py`); BASE families wall time ≈ 21 min; the BASE GTF sha1 50f239d6
  (parked prereg §13); the r1136-r1139 genome-wide rows. Predictions 1 and 2's "few ghosts" lean rests on the 1.9
  isoforms figure and is marked post-exposure. This author saw none of those files.
- **Apes:** the NPIP contigs were read descriptively in npf_audit and npf_variants (holders, fused copies, RG3 split
  census) — they are controls, not verdict substrates, and P7 says so.
- **chr16 of A119b is dev; chr16 of testis is not** (memory 09-28): the dev-contig exclusion of the parked prereg §7.2
  does not apply to testis, which is the whole point of this test.

## 12. Order, stop rules, machine rules, cost

**Order.**
1. **Freeze.** The main session records this file's sha1 and the user's acceptance (§11 reuse, §12 cost) in Amendment 1.
   Until then: dev only; no testis product opened by anyone.
2. **Dev, in parallel (never a testis file):** the port agent builds the port and runs G2(a-e); the run agent runs G1,
   G5, G6 (on dev products), G7, and freezes SHA1_NPF / SHA1_U / SHA1_SCORE. The port agent never opens a testis product.
3. **Held-out (run agent only, after the freeze), per substrate in this order: human_testis, gorilla_OR6737,
   gorilla_KB3781, chimp_PTR, orangutan_PPY (G8 first).** For each: build U (structure only; write `U.tsv` with the
   construction counts); GTF_U(BASE); `rg3.py` whole GTF + restriction (+ G4); `rg3_null.py` × 5 (apes × 1); the
   families runs (heavy, one at a time; ≈ 7 for testis, 3 per ape, + the orangutan re-run); G3; G2(f); scoring
   (light): SHA1_NPF per arm, `family_score` × truths × filings, `pair_sets`, `cls.py`, C4; write the substrate's
   table before the next substrate starts. **Nothing about the testis result changes anything downstream.**
4. **Verdict** (§6), Outcome appended to this file, one register row per numbered claim, Amendment 2 with every sha1
   and command line; the port ships only after G2 passes in full.

**Stop rules.** No change after the freeze to the rule, the universe rule, a null, a clause, a filing, a truth file or a
tolerance; a bug fix re-runs every affected step and is recorded as an amendment. A failed gate stops its substrate
(clauses not measured), never touches a rule. A U families run that hits the 600 s call cap is re-run ONCE with
`RLOCK_TIMEOUT=1800` (one process, foreground); if it fails again, that substrate's family clauses are not measured.
No variant is ever substituted (no `rg.py`, no donor|acceptor, no V4s, no genome-wide run).

**Machine rules.** Heavy steps (`cargo build/test`, `mcl_families`, the orangutan re-run, any minimap2) via
`bash tools/rlock.sh heavy <cmd>`, one at a time, foreground; light steps (Python scoring, `family_score`, `rg3.py`,
`cls.py`, emu) via `bash tools/rlock.sh light <cmd>`. Build ONLY with
`CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target`, `--release`, cargo output to a FILE, `touch` every
edited `.rs` first. `TMPDIR` under `/mnt/linuxdisk`. Never `pkill -f`; kill by PID. Scratch:
`/mnt/linuxdisk/tmp/rustle_figures_dev/rg3_<key>/`; delete own `loci.fa` / `copies.fa` / PAF copies after G3 and the
scoring, keep `clusters.tsv`, `loci.tsv`, `params.tsv`, logs, GTF_U arms and every table (small). Do not commit or push
(the main session does).

**Cost (estimate; the user accepts it).** Port: build 10-20 min × 2, `cargo test --release` 15-25 min, dev
regenerations 3 × ≤ 1 min, comparisons — ≈ 1-1.5 h. Held-out: U construction = one Python pass over each PAF
(0.43-0.53 GB) ≈ 1-3 min × 5; families on U: |U| is expected in the hundreds to low thousands of loci (dev chr16,
2,802 loci, ran in 43-72 s), so ≈ 20 runs × ≤ 2 min ≈ 40 min, + 6 min orangutan re-run; scoring ≈ 20 min. **≈ 1.5-2.5 h
of the heavy lock, ≈ 3-5 h wall with gates and reviews.** Disk: < 3 GB transient, < 200 MB kept.

## 13. Not in this test (and why)

- A genome-wide RG3 verdict (the parked prereg: FU, S, P1, P2, H, T, L on six samples; ≈ 18-24 h, ≈ 55-60 GB). This
  test does not decide the default; it decides whether the block-level gain replicates on a second library.
- The polish's role in making ghosts (npf_critique §3: under the pre-09-23 polish the bridge may come back). Needs a
  re-assembly with `--polish-retained-ratio 0 --assembly-junctions majority`; it is the next question if F2 holds.
- Any cut through a real read bridge (V4/V4g/V4s, npf_variants): all failed split correctness on dev; V4s is not a
  secondary arm (npf_critique §4).
- The fusion relation R (npf_variants §6): output only; not scored (r845's `multi` policy lowers bipartite F).
- Dishuck's expression subgroups; Soto or U2 as judged truths (a cover and an aligner-sharing truth; both reported).
- A merge rule for BASE `gene_id`s that already share an exonic base; the old `rg.py`; donor|acceptor; the `.<k>` names.
- Reading any testis number before the freeze, by anyone.

## 14. Hostile self-review (the advisor's role), with the fixes applied

1. **"C1 (iii) makes 'no ghost in this library' a failure, so the test cannot distinguish 'no effect' from 'refuted'."**
   Fixed: the verdict table (§6) sends the no-effect case to KEEP OPT-IN and reserves REFUTE for harm or two failures;
   prediction 9 puts 0.55 on exactly that outcome.
2. **"Family size multiplies one unit into ~20 pairs (the 'selecting which component' trap); a pair clause rewards
   joining the biggest family."** Fixed: the primary clause C1 is per copy; C2's false-pair half is a set difference at
   tolerance 0 (a wrong member of a 30-family costs 29 pairs and fails it); revealed copies are counted per piece.
3. **"`gene_at` is strand-blind and names 9 of 21 NPIP holders after a partner; pair truths are partly a function of
   the prediction's geometry."** Acknowledged in §4.2; identical bias across arms on the same U; the holder rule of §5.1
   is exon-level and same-strand; the touching-unit table is reported per copy.
4. **"U is defined from BASE, so the arms are not symmetric."** U is the only construction that does not depend on an
   arm's outcome; RG3's pieces are added; every arm is aligned on the same target set; BASE_U vs genome-wide BASE
   concordance is reported so a U-restriction artefact is visible. The one-step closure is stated, not iterated, to
   keep U bounded; its size and composition are reported.
5. **"The null restricted to U may have no donors (shortfall) and collapse to BASE."** Stated in C1 (iii); shortfalls
   are reported; five draws are fixed by label.
6. **"The truths disagree on the readthrough copies, so any filing choice is a thumb on the scale."** Both filings are
   computed for every truth; the judged truths must pass under both; Soto/U2 gains are decomposed into filing /
   pre-existing contamination / other. Judging Soto/U2 was ruled out for stated reasons (cover; shared aligner) and,
   honestly, after seeing that they alone produced "false" pairs on dev — disclosed here so the choice is visible.
7. **"Same genome, same truths, same block: this is dev with different reads."** §11 says exactly that, names it the
   weakest held-out form, and the verdict label carries "block-level, one human library". The apes cannot supply a
   positive test (census); they are controls.
8. **"n is tiny: C3 may be judged on one split; C1 on ±1 copy."** Stated; C3 is judged at ≥ 1 with n reported, and a
   single FRAG with no SEP fails C3 but is not harm (KEEP OPT-IN unless another clause fails); the 10-relabelling gate
   protects a ±1 copy result against MCL tie order.
9. **"The genome-wide protein referee does not exist; the chr16 table cannot judge U loci off chr16."** Stated in §4.2;
   Compara ALL judges them; unjudgeable counts reported per truth.
10. **"The port is 'shipped' before it is tested on the block."** The port never produces a judged number; G2 requires
    byte-identity with `rg3.py` on dev and on every substrate's BASE GTF; the ship waits for G2, the verdict does not
    wait for the port.
11. **"The ape NPIP contigs were already read; P7 is not a prediction."** Marked as a known census; the blind residue
    (re-flow from non-holder splits) is what A tests.
12. **"A U families run is a different alignment from the genome-wide one (r1126), so 'BASE_U' is not the shipped
    BASE."** Correct and stated: the comparison is within U, and the concordance row shows how far BASE_U is from the
    shipped families on those loci.
13. **"Ghost detection depends on the polish constants (retained-ratio 10, strict junctions, 2-read floor)."** True;
    RG3 inherits them and adds none; the polish's role is the next question (§13), not this test's.
14. **"The prereg was written by an agent that read the dev results, including the +8/+12 Soto/U2 pairs."** Yes; every
    such influence is named in the header and in item 6. No testis number influenced anything.

## Amendments

### Amendment 1 — the freeze (2026-09-28, run agent KEY=run, written BEFORE any human_testis product was opened)

**Acceptance.** The user's decision was relayed by the orchestrating session on 2026-09-28 ("ok lets ship the
regroup RG3 and test"; decisions: ship RG3; test it held-out on the NPIP block of human_testis; apes as no-change
controls). The user's own sentence does not name §11 (reuse) or §12 (cost); their acceptance is carried by the
orchestrating session's decision list, and this is recorded as such rather than as the user's literal words.

**This file's sha1 before this amendment:** `b5e465d09382bd47271fe02a4f81adb8255f957b` (53,545 bytes).

**Frozen instruments** (`/mnt/linuxdisk/tmp/rustle_figures_dev/rg3_run/lib/`; the frozen inputs of §10 were copied
there with their sha1s verified: rg3.py ec17e540, rg3_null.py 6ece9ae0, test_rg3.py 03e10572, emu.py fac9a560,
cls.py 41339094, relabel_null.py 8c050970, run_emu3.py b6ab23d8, npf.py fdc1a7d3, audit.py fd0305ca, rg_score.py
fea2167b, score_arm.py 258d1886, npip_units.py 6dafdd7d, ppy_truth.py 57be8160):
- **SHA1_NPF** = `npf_score.py` `ed04211134a6ac16b11235c1be4a05bd47377c60` — §5.1 implemented verbatim over the frozen
  `npf.py` / `audit.py` truths (holder = max same-strand exonic overlap, ties by K start then gene_id; NPIP = most
  copies, then more members, then lower MCL index), with a `--u` restriction and npf_audit's description columns.
  *Deviation from §10's wording* ("score_arm.py with paths fixed and `--u`, nothing else"): `score_arm.py` is
  chr16/human-only and breaks the holder tie by reads rather than by §5.1's rule, so §5.1 was re-implemented; ties
  are counted and reported (0 on dev).
- **SHA1_U** = `build_u.py` `951b07d8a7a812cfbf808329cdd7869c461bc288` (§2 verbatim; PAF columns 1/6/10/11 only).
- **SHA1_SCORE** = `score_u.py` `3e5cf0e298b83149366fc97aec1b988ad3504418`, with `pairs.py`
  `e1e0340904ca91504288caacc60a53a9e475e032` (family_score 7723029b wrapper, filings A/B of §4.3, `pair_sets` of
  rg_score fea2167b generalised to many contigs with both the first-family and any-shared-family TP rule),
  `make_gtf_u.py` `9ebba32d14a1adb6f3fcf97b63e332946967e488` (GTF_U of §3), `g3_emu.py`
  `8a7fe6aebda9f532077a391a9abd9d3d209daa2f` (G3: emu R0 + 10 relabellings), `run_fam.sh`
  `026375a9b2e726539cd83b08bf18e21ce4c9c2a8` (the §3 `mcl_families` command verbatim, `/usr/bin/time -v` to a
  separate `.time.log`).

**Port (G2, from the port agent's report `rg3_port.md`).** Uncommitted working tree on `main@67f86286`
(`src/bin/copy_assign.rs` +511 lines, `tools/rustle_pipeline.sh`); flag `--gtf-regroup` (clap bool, default off,
beside `--gtf-tpm`, refuses a run without `--gtf`/`--assemble-only`); driver `RUSTLE_GTF_REGROUP=1` (`""|0` off, else
exit 2); frozen build `/mnt/linuxdisk/tmp/rustle_figures/rg3_bin_frozen/copy_assign` 9fb8b6b8, `source.diff`
5d44fdc1. G2(a) unset: 6/6 products cmp-identical to `fj_bin_frozen/copy_assign` b13b6ae6 on chr16 / chr20 /
NC_073244.2 (also under R3); G2(b) set: GTF cmp-identical to `rg3.py` 3/3 (+ R3 3/3), splits 19/10/4, other products
identical, `params.tsv` = off + 4 rows; G2(c) the 19 fixtures as 9 Rust unit tests, 9/9; G2(e) `cargo test --release`
21 suites, 954 passed / 0 failed / 13 ignored. **G2(d) NOT MEASURABLE — the port has no GTF-to-GTF path** (no such CLI
exists in `copy_assign.rs` or `source.diff`); hence **G2(f) is not measurable either.** Per §1.4 this is a port item:
the judged arm is `rg3.py`'s output and the verdict does not wait; the ship waits for the port agent to add the path.

**Dev gates (run agent).**
- **G1 PASS** on hsa16 / hsa20 / ggo44: `rg3.py` on the regenerated BASE GTFs is cmp-identical to the port agent's
  `rg3_port/ref/<s>.on.rg3py.gtf`, its `--stats` equal the frozen `rep3_RG/rg/<s>.BASE.RG3.json` (19/10/4 splits,
  171/30/12 relabelled, 1957/293/78 lines), it is idempotent, and its locus set (gene_id, span) equals the frozen
  families run's `rep3_RG/fam/<s>.BASE.RG3.loci.gff3` (2821 / 1892 / 1145 loci; ggo44's log says 1146 because two
  gene_ids share one key); `rg3_null.py` (labels `human_A119b` / `gorilla_OR6737`) is byte-identical under
  `PYTHONHASHSEED=0` and `=12345`, its stats equal `rep3_RG/null/<s>.BASE.NULL3.json` and its locus set equals
  `rep3_RG/fam/<s>.BASE.NULL3.loci.gff3`; `test_rg3.py`: 19 fixtures pass. *Deviation:* the rep3_RG RG3/NULL3 GTFs
  named in §10 were not kept (only their `--stats` JSON; the `rep4_incr/full` symlinks dangle) and the dev BASE GTFs
  `rep2_RG/dev_gtf/*` were deleted, so the BASE inputs are the port agent's frozen-binary regenerations
  `rg3_port/runs/frozen/<s>.off.unset.gtf` (the port report: same transcript counts and `rg3.py` stats as rep3_RG's
  inputs) and the byte targets are the three substitutes above.
- **G5 PASS**: SHA1_NPF on the dev BASE (`rt3` families `rep_edges/base/hsa16`, and npf_variants' `fam/hsa16.BASE`)
  gives 22 of 26 in MCL1 (27 rows) with exactly npf_audit's exceptions (NPIPB2 singleton, PKD1P6-NPIPP1 MCL25, NPIPB5
  MCL26, NPIPB6 MCL27), 0 holder ties, audit's holders identical on 26/26; on V1 (`fam/hsa16.V1`) 24 of 26 (29 rows),
  NPIPB2 and NPIPB6 gained as pieces. *Deviation:* npf_variants' NULL1_k GTFs and their `rg3_null.py` labels were not
  kept (≈ 80 candidate labels tried; none reproduces the frozen donors), so the NULL_V1 = 22 row is cited from the
  frozen `npf_variants/out/hsa16.NULL1_k.score.json` (22 on 5/5, per-copy rows identical to BASE's) instead of being
  regenerated; SHA1_NPF's null handling is covered by G1's NULL3 reproduction and by the end-to-end dev check below.
- **G6(dev) PASS**: `pair_sets` = family_score on the dev BASE for Compara (66 predicted / 66 TP) and the chr16
  protein referee (101 / 100); Soto (124 vs 103) and U2 (169 vs 145) do not reconcile, as expected for covers (never
  judged). Filing B per truth: Compara NPIP family CF153 (19 of the 26 symbols), protein PF57 (21), Soto ID_154 (14),
  U2 ID_154 (19); NPIPA1 occurs in every table (so B ≠ A for Compara and the referee too — both filings are computed
  everywhere, as §4.3 provides); Soto adds NPIPA1 / PDXDC2P-NPIPB14P / PKD1P6-NPIPP1, U2 all six.
- **G7 PASS**: `cls.py` reproduces every row of `rep3_RG/cls/cls.json` (RG3 16/2/1, 5/3/2, 4/0/0; junction rule
  15/38/31, 4/22/11, 3/4/1; NULL3 0/19/0, 0/9/1, 0/4/0).
- **End-to-end check of SHA1_SCORE on dev chr16** (pseudo-U = all 2,802 chr16 gene_ids; BASE = npf_variants
  `fam/hsa16.BASE`, RG3 = `fam/hsa16.V1`, NULL = rep3_RG `fam/hsa16.BASE.NULL3`): C1 22 → 24 (NPIPB2, NPIPB6; both
  holders are pieces; null 22); C2 Compara 66 → 87 TP with 0 new false pairs, protein 100 → 125 with 0 (fp 1 → 1),
  Soto +8 and U2 +12 false pairs all of class (i) or (ii), identical under both filings; C3 RG3 16/2/1 vs NULL 0/19/0;
  largest 27 → 29 with π = 2; concordance 0 changes (same run). **C4(b) would FAIL on dev**: the 2-member family
  MCL73 (`DN_chr16_20664599_5`, split by RG3, and the antisense `DN_chr16_20714140_9`, whose span lies inside the
  keeper's pre-split span) dissolves once the keeper's span shrinks — the failure mode §5.4 names; the clause is kept
  as registered and every lost locus is described (class, BASE cluster, overlap with a split locus's old span).
- **G8 note (deviation to record before the run):** the orangutan reference run `npf_audit/ppy/fam/ppy.clusters.tsv`
  was made with `rt3_bin_frozen/mcl_families` 0c4639e2 (its `time.log`), not with `fj_bin_frozen`'s 91ef2e1c; the
  byte-for-byte re-run of G8 therefore uses that binary and the same command on a `ppy4.gtf` rebuilt as the BASE GTF
  restricted to the 4 NPIP contigs (NC_072382.2, NC_072383.2, NC_072387.2, NC_072391.2); the U arms use
  `fj_bin_frozen/mcl_families` as §3 states.

**Command lines (held-out, per substrate s, run in this order):**
```
python3 lib/build_u.py --species SP --truth TK --gtf <s>.BASE.gtf --loci <s>.BASE.fam.loci.tsv \
    --clusters <s>.BASE.fam.clusters.tsv --paf <s>.BASE.fam.loci.paf --out U/<s>.<TK>.U
python3 lib/make_gtf_u.py <s>.BASE.gtf U/<s>.<TK>.U.tsv U/<s>.<TK>.BASE_U.gtf
python3 lib/rg3.py <s>.BASE.gtf U/<s>.RG3.whole.gtf --stats U/<s>.RG3.whole.json            # whole GTF (asserts 1-5)
python3 lib/make_gtf_u.py U/<s>.RG3.whole.gtf U/<s>.<TK>.U.tsv U/<s>.<TK>.RG3_U.gtf --pieces-of <s>.BASE.gtf
python3 lib/rg3.py U/<s>.<TK>.BASE_U.gtf U/<s>.<TK>.RG3_U.g4.gtf && cmp U/<s>.<TK>.RG3_U.gtf U/<s>.<TK>.RG3_U.g4.gtf   # G4
python3 lib/rg3_null.py U/<s>.<TK>.BASE_U.gtf U/<s>.<TK>.NULL<k>_U.gtf --label <s>.U<k> --stats ...   # k = 1..5 (apes: 1)
bash tools/rlock.sh heavy bash lib/run_fam.sh <species> U/<s>.<TK>.<arm>_U.gtf fam/<s>.<TK>.<arm>     # per arm
python3 lib/g3_emu.py U/<s>.<TK>.<arm>_U.gtf fam/<s>.<TK>.<arm>.fam.loci.paf fam/<s>.<TK>.<arm>.fam.clusters.tsv emu/<s>.<TK>.<arm> [--relabel]
python3 lib/score_u.py --sample <s> --species SP --truth TK --u U/<s>.<TK>.U.tsv --arms score/<s>.<TK>.arms.json \
    --out score/<s>.<TK> [--truth-tables] [--cls-species human|gorilla] --gw-clusters ... --gw-loci ... --base-whole-gtf ... \
    --rg3-stats U/<s>.RG3.whole.json --null-stats NULL<k>=...
```
Scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/rg3_run/`; tables `/mnt/linuxdisk/tmp/rustle_figures/rg3_npip/tables/`.
Machine note: a foreign `minimap2` job of another session (PID 3203385, 10-14 GB, 5 threads, not under the lock)
was running during the dev gates; the heavy U runs are started under `tools/rlock.sh heavy` regardless, one at a time.

**Amendment 1 (the freeze; before any held-out command) must contain:** the user's explicit acceptance of this file,
of the reuse (§11) and of the cost (§12); this file's sha1; SHA1_NPF, SHA1_U, SHA1_SCORE; the port's commit, flag name
and GTF-to-GTF CLI, with G2(a-e) results; G1, G5, G6(dev), G7 results; every command line.

**Amendment 2 (after the run) must contain:** per substrate |U| and its construction counts, the families runs'
wall times, G3/G4/G8 and G2(f) results, the clause table, the verdict, the register rows, and the paths of every kept
product.

## Amendment 3 (orangutan stopped by the user)

2026-09-28, written by the verdict agent after the run. This file's sha1 before this amendment:
`ffe5467e1ece053991eab5ec3a4ef7c46ecf0025` (63,055 bytes).

**The user's decision.** At about 16:30 the user said: "lets stop the orangutan one, I only need human and gorilla for my
meeting." The files on disk show where the run stopped. The G8 re-run had finished and passed: `clusters.tsv`, `loci.tsv`,
`copies.tsv` and `params.tsv` were byte-identical to `npf_audit/ppy/fam`, in 5:25 wall, run with `rt3_bin_frozen`
as Amendment 1's G8 note says. The driver had printed `== orangutan_PPY member: U` and had written nothing more. There is
no orangutan U, GTF_U, families run on U or clause product, and no process was left running (`ps` was checked; nothing
was killed by this agent). **The orangutan control is recorded as NOT RUN, stopped by the user.** It is neither a gate
failure nor a pass, and nothing is inferred from it.

**How this enters §6.** Control A is judged over the ape samples that were measured: gorilla OR6737, gorilla KB3781,
and chimp PTR under both truth sets. Orangutan PPY is *not measured*. Under §6 a clause that is not measured caps the verdict at
KEEP OPT-IN. The human clauses already put the verdict at KEEP OPT-IN (F = {C1}), so the cap changes nothing. The only
thing the orangutan run could still have changed is **REFUTE**: had A failed on PPY, F = {C1, A} and |F| = 2. It could
never have produced EFFECTIVE, because C1 (iii) fails on human_testis. The verdict below therefore reads "on human +
gorilla (+ chimp); orangutan control not run". The user's scope for the meeting is human + gorilla. Chimp finished
before the stop, is reported, and does not change the verdict.

**Amendment 2 was not appended.** The run agent's report (`rg3_run.md`) says Amendment 2 was appended after the run.
It was not: before this amendment the file ended at Amendment 1 (mtime 14:30). The Outcome below carries everything
Amendment 2 had to contain.

**Instrument findings of the independent recomputation.** No rule, universe, null, clause or truth was changed.
1. **Filing B is a no-op inside `family_score`.** `pairs.py` builds filing B by *appending* rows, and `family_score`
   (`soto_truth` / `truth_rows_gw`) keeps only the FIRST family per gene. An appended row for a gene that already has a
   row is therefore ignored. §4.3's premise, "`family_score` accepts repeated genes", is true only in the sense that it does
   not crash. **Consequence here: none.** In the judged truths the only added name is NPIPA1, which has no locus in
   testis U. The independent recount with cover semantics (a gene belongs to every family it is filed under, and a pair is
   true iff the two genes share any family) gives identical counts on all four truths under both filings. A future filing
   B must PREPEND its rows, or use cover semantics.
2. **The run's own Soto / U2 pair sets do not reconcile with `family_score`** (G6 "fails" for Soto and U2). The cause
   is the universe, not the cover. `pairs.pair_sets` admits genes whose truth family has a single gene, whereas
   `family_score`'s universe keeps only families with at least 2 genes. That produces the run's "fp 1 → 1" (Soto) and
   "fp 2 → 2" (U2). With `family_score`'s universe, Soto and U2 have **0 false pairs in every arm**. These are reported
   truths only, and there is no new false pair under either reading.
3. **The identity gate "same transcript count"** (§3) was checked as keeper + pieces. Read literally as the keeper's
   own count, it differs for the one split `gene_id` by construction (8 → 3 + 5). This is a wording defect, not a bug.
4. **The pre-Amendment-1 sha1 `b5e465d0…` cannot be re-derived.** No copy was kept, and deleting the Amendment 1 block
   leaves 53,531 bytes against the recorded 53,545 (a 14-byte placeholder difference). The ordering on disk is
   consistent with a freeze before the held-out work: frozen instruments mtime ≤ 14:29:55, this file 14:30, first testis
   product (`U/human_testis.dishuck.U.tsv`) 14:31.
5. **G2(d) and G2(f) are still NOT MEASURABLE.** The port (uncommitted, `src/bin/copy_assign.rs` +511) has no
   GTF-to-GTF path, so the ship still waits for it. As §1.4 provides, the verdict does not.

## Outcome

### O.1 Independent recomputation and provenance (verdict agent)

The verdict agent wrote its own code: GTF parsing, K(g), U (§2), RG3 (§1.1), holder / NPIP / c1 (§5.1), pair sets
(`family_score`'s `gene_at` and universe, §4.2), C3 classes (§5.2), C4, and control A (§5.5). The scripts are in
`/mnt/linuxdisk/tmp/rustle_figures/rg3_npip/verify/`: `common.py` b7d1a47f, `verify_human.py` 5b6f69af,
`verify_ape.py` e0505ab5, `cls_ggo.py` cacbfadf. JSON outputs sit beside them. The only code taken from the run's `lib/`
is the ape truth loaders `audit.ape_member_truth` / `audit.ptr_native_truth`, which are §4.5's truth definitions.
**Every judged number of the run was reproduced.** The only differences are the reported-only Soto / U2 false-pair
counts (Amendment 3, item 2).

- **U built as §2 defines, from BASE only.** The recount gives the same U, gene_id for gene_id, on all 5
  (sample, truth) rows. It uses PAF columns 1/6/10/11 with block ≥ 300 and 10·matches ≥ 7·block (integer form).
  Human: L 13,012, H0 18, F0 5, U1 23, N \ U1 271, **|U| 294**, and 622,353 of 628,150 q ≠ t records pass.
- **Arms.** GTF_U(BASE) is exactly the restriction of the whole BASE GTF (`#` lines dropped). The recount's RG3 names
  equal `rg3.py`'s on every transcript of the 4 whole GTFs. GTF_U(RG3) equals the restriction, and G4 holds. Human:
  NULL_1..5 GTF_U are byte-identical to GTF_U(BASE), and their PAF and `clusters.tsv` are byte-identical to BASE_U's.
  The shortfall was re-derived: the chrX pool holds 8 donors, the most spliced non-representative transcripts in any
  donor is 2, and T_g = 5, so no donor is eligible.
- **sha1s verified.** `rg3.py` ec17e540, `rg3_null.py` 6ece9ae0, `test_rg3.py` 03e10572, `emu.py` fac9a560, `cls.py`
  41339094, `relabel_null.py` 8c050970, `run_emu3.py` b6ab23d8, `npf.py` fdc1a7d3, `audit.py` fd0305ca,
  `rg_score.py` fea2167b, `score_arm.py` 258d1886, `npip_units.py` 6dafdd7d, `ppy_truth.py` 57be8160;
  SHA1_NPF `npf_score.py` ed042111, SHA1_U `build_u.py` 951b07d8, SHA1_SCORE `score_u.py` 3e5cf0e2 + `pairs.py`
  e1e03409, `make_gtf_u.py` 9ebba32d, `g3_emu.py` 8a7fe6ae, `run_fam.sh` 026375a9. Binaries: `mcl_families` 91ef2e1c
  (every U run's `time.log` names it, together with the §3 flags and the registry fasta), `family_score` 7723029b.
  Inputs: `lit_subclusters_npip_dishuck_check.tsv` 6b6a1027, testis BASE GTF 50f239d6 (= §11), protein referee ebc7227c
  with gff 1bd6fafa.
- **No dev number enters.** The testis and ape score inputs (`finish_testis.sh`, `run_ape.sh`, the `arms.json`
  files) point only at `rg3_run/{U,fam}` and `rt_arms/<sample>`. A search of the testis and gorilla score products for
  A119b / hsa16 / `dev_small` / `rep3_RG` / `npf_variants` finds nothing. The dev gates G1 / G5 / G6(dev) / G7 are used
  only as gates.

### O.2 human_testis (the verdict substrate)

n_present = **7** of 26 Dishuck copies: NPIPB4, B5, B6, B7, B9, B10P, B14P. Twelve copies have no span-overlapping
locus, and 7 more have loci but no same-strand exonic overlap. All present copies are tiny (1-2 transcripts, 2-18 reads).
Holder ties: 0 in every arm. RG3 splits 7 of 13,012 `gene_id`s on the whole GTF (0.05%, 33 transcripts relabelled),
and **exactly 1 inside U**: `DN_chrX_54382055_5`, 8 transcripts → keeper 3 + piece 5. G3 holds on all 7 runs (emu R0 byte-equal;
10/10 relabellings for BASE_U and RG3_U). The families runs on U took 3:19-3:55 wall and 3.14 GB each.

| clause | BASE_U | RG3_U | NULL_1..5 | result |
|---|---|---|---|---|
| **C1** c1 (NPIP = MCL1, 5 rows: B6, B7, B9, B10P, B14P) | 5 | 5 | 5 ×5 | (i) PASS (0 copies leave NPIP; holders identical), (ii) PASS, **(iii) FAIL, the no-effect case** (5 = 5 = max NULL; every NULL draw a total shortfall) |
| **C2** Compara (ALL), filings A and B | tp 4, fp 0 | tp **5**, fp 0 | tp 4 | PASS: 0 new false pairs; the gained pair is PAGE2-PAGE2B (CF400), carried by the RG3 piece; revealed copy PAGE2B |
| **C2** protein referee (chr16), filings A and B | tp 6, fp 0 | tp 6, fp 0 | tp 6 | PASS: 0 new false pairs |
| Soto / U2 (reported, both filings) | 5 / 4 tp, 0 fp | 5 / 4 tp, 0 fp | = BASE | 0 new false pairs (the run's own count: fp 1→1 / 2→2, Amendment 3 item 2) |
| **C3** split correctness (judged, n = 1) | – | SEP 1 / FRAG 0 / UNJ 0 | 0 splits | PASS: keeper = PAGE2, piece = PAGE2B, one annotated gene each |
| **C4** (a) NPIP members kept; (b) U loci still clustered | 5; 23 | 5/5; 23/23 | – | PASS |

Reported: 6 → 7 families, 23 → 25 loci clustered. The new family is MCL6 = {keeper, piece}. Largest family 7 → 7, with 0
pieces in it (the H form holds). c1 under the genome-wide BASE clusters is also 5. BASE_U vs genome-wide BASE: 0 of 294
loci change co-membership inside U (27 are clustered genome-wide against 23 in the U run; 4 lose partners outside U).
Unjudgeable predicted pairs: 36-39 per truth. Bipartite F (continuity only): Compara .009 → .012, protein .045 → .045.

**The A119b ghost links do not recur.** NPIPB2 has no locus in testis. NPIPB6's holder `DN_chr16_28623027_8` is an unfused
1-transcript locus already in NPIP. The only fused NPIP holder is NPIPB14P's, which contains the PDXDC2P-NPIPB14P readthrough
(1 fusion transcript, 2 reads). That readthrough dominates the locus's representative, the locus is one RG3 piece, and it is
placed in NPIP under every arm.

### O.3 Ape no-change controls (clause A; control only, never pooled with human)

| sample / truth | \|U\| | present | c1 BASE / RG3 / NULL1 | RG3 splits in U (holders) | A | whole U partition (pieces mapped) | C3 in U (reported) | families wall |
|---|---|---|---|---|---|---|---|---|
| gorilla OR6737 / T_member (25) | 436 | 11 | 10 / 10 / 10 (MCL0, 14) | 3 (0) | **PASS** | identical | SEP 1 / FRAG 2 (SEC14L1 pure; GPRASP1 cut) | 5:21-5:41 |
| gorilla KB3781 / T_member (25) | 256 | 8 | 8 / 8 / 8 (MCL0, 10) | 2 (0) | **PASS** | identical | SEP 2 / FRAG 0 | 3:20-3:36 |
| chimp PTR / T_member (19) | 349 | 16 | 16 / 16 / 16 (MCL0, 40) | 1 (0) | **PASS** | identical | no genes table | 1:42-1:54 |
| chimp PTR / T_native (46) | 356 | 37 | 35 / 35 / 35 (MCL0, 40) | 1 (0) | **PASS** | identical | no genes table | 1:40-1:42 |
| orangutan PPY / landing + native | – | – | – | – | **NOT RUN, stopped by the user** (G8 passed) | – | – | – |

G3 holds on all 12 ape runs (emu R0 byte-equal, and 10/10 relabellings for BASE and RG3). G4 holds. C4 holds on every
ape row. Whole-GTF RG3 splits: OR6737 49 / 19,597 (0.25%), KB3781 124 / 17,420 (0.71%), PTR 9 / 18,848 (0.05%).
Concordance between the U runs and genome-wide BASE: 2 / 2 / 0 / 0 loci change co-membership.

### O.4 Verdict (§6)

F = {C1}. C1 fails only in (iii), the no-effect case. There is no harm: C1 (i) holds, c1 did not fall, there is no new
false pair on either judged truth under either filing, and C4 holds. A passes on every measured ape sample and truth set.
Orangutan is not measured (Amendment 3).

**VERDICT: KEEP OPT-IN, on human + gorilla (chimp also unchanged); the orangutan control was not run because the user
stopped it.** It is reported as "no ghost at NPIP in this library". RG3's A119b NPIP gain (NPIPB2, NPIPB6) did **not**
replicate on human_testis because the precondition is absent: testis has no ghost-bridged NPIP copy. The rule did no harm
on the block. It made one correct separation (PAGE2 / PAGE2B), and that adds one true Compara pair. The only way orangutan
could still change the verdict is REFUTE (a PPY change would make |F| = 2). The default flip remains the user's call, and
the port's ship still waits for G2(d) and G2(f).

### O.5 Predictions (§7) and falsifiers (§8)

- P1: asserts pass 4/4 run (PPY not run). ≥ 1 split in testis U: yes (1). "Pieces 0.2-1% of loci genome-wide, testis
  the fewest": testis 0.054% and PTR 0.048% are below the range, OR6737 0.25% and KB3781 0.71% are inside it; **missed**
  (PTR, not testis, is the fewest).
- P2: c1 difference = 0 (0.55): **hit**. The NPIPB2 / NPIPB6 ghosts in testis (0.30 each): neither. The ghost set differs from
  A119b's (0.7): **hit** (empty).
- P3: max NULL = BASE (0.85): **hit**. ≥ 1 NULL draw a shortfall (0.5): **hit** (5/5 total shortfalls, so NULL = BASE
  and C1 (iii) reduced to c1(RG3) > c1(BASE)).
- P4: 0 new false pairs (0.80), tp non-decreasing (0.90): **hit**. Soto / U2 false pair iff a C1 gain: no gain and no false pair,
  consistent.
- P5: C3 judged (0.7), passes (0.7): **hit** (n = 1). FRAG ≥ 1 in human U (0.35): no.
- P6: C4 (a) and (b): **hit**. P7: all measured ape rows unchanged: **hit** on 3 of 4 samples (PPY not run). Non-holder splits
  in ≥ 1 ape: **hit** (3 / 2 / 1). 0 holder splits: **hit**.
- P8: G2 (a-c, e) passed per the port report (the build count was not re-checked). (d) and (f) are not measurable.
- P9: KEEP OPT-IN (0.55): **hit**.
- F1: not measurable (no gain). F2 ("ghosts are a library property"): **not falsified**. A119b's ghost set does not recur;
  testis has none at NPIP. F3: holds but is vacuous here (every NULL draw was a total shortfall). **F4 ("RG3 cannot fragment a
  gene on the block")**: holds on human U (0 FRAG) and on KB3781 U (0). **It is falsified on gorilla OR6737 U, with 2 FRAG
  among 3 splits**: `DN_NC_073228.2_14760123_17` (SEC14L1, pure fragmentation) and `DN_NC_073247.2_114501773_5`
  (GPRASP1 cut; the keeper also spans ARMCX5 / GPRASP2 / GPRASP3). Both are non-holder U neighbours, neither piece joined
  a family, and A holds. C3 is reported only on apes and never judged.

### O.6 Caveats (read before quoting)

- **n is tiny.** Only 7 NPIP copies are present in testis, all 1-2 transcripts with 2-18 reads, and 12 of 26 copies have
  no assembled locus. The block-level result is "no effect because the precondition is absent", not "no effect when a
  ghost is present".
- The test is held out **in library only**: same genome, annotation and truths as dev (§11). After this test,
  human_testis's NPIP block is spent for regrouping work.
- The NULL has no power on testis: all 5 draws were total shortfalls, so C1 (iii) compared RG3 with BASE directly.
- C3 was judged on 1 split. On gorilla OR6737 U, RG3 fragmented 2 of 3 split genes (reported, F4 above). This is the
  same failure class as dev's 5 / 33 FRAG, and it matters for any default flip. The genome-wide split-correctness question
  (the parked prereg) is still open.
- Filing B and the Soto / U2 own pair sets have the instrument defects of Amendment 3 (items 1-2). Neither affects a
  judged number.
- Orangutan: not run, stopped by the user (Amendment 3).

**Kept products.** Tables: `/mnt/linuxdisk/tmp/rustle_figures/rg3_npip/tables/` (human_testis.*,
gorilla_OR6737.member.*, gorilla_KB3781.member.*, chimp_PTR.{member,native}.*). Verification:
`/mnt/linuxdisk/tmp/rustle_figures/rg3_npip/verify/`. Scratch, with GTF_U arms, U runs, emu and scores:
`/mnt/linuxdisk/tmp/rustle_figures_dev/rg3_run/`. G8 re-run: `rg3_run/ppy/fam/`. Key sha1s: human U.tsv b2acdb73,
GTF_U(BASE) 14cf22cc, GTF_U(RG3) 2e948a69, BASE_U clusters 4f3e3cea, RG3_U clusters 6c2228b8; ape U.tsv OR6737 7712dd43,
KB3781 d3ab91e6, PTR member 3a515044, PTR native 1d595308.
