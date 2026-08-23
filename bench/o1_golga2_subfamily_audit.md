# GOLGA2 versus GOLGA6/8: family or false merge?

**Decision (2026-08-16):** the chr9 node is a real GOLGA2 transcript and a documented
ancestral homolog of the duplicated GOLGA6/8 family. It is not an unrelated false-positive
node. It should remain visible in the **broad RNA homology family**, but should be labelled
`RELATED_OUTGROUP` and excluded from the **recent-copy GOLGA6/8 subfamily**.

## 1. Node identity

Fresh Rustle emits `chr9:140463642-140483173`, minus strand, 26 exons and 39 reads. The
current T2T RefSeq GOLGA2 locus is `chr9:140462943-140483121`, minus strand, with 27 exons:
the prediction is a near-full-length GOLGA2 model. SWI5 begins at `chr9:140482470`, so the
old `GOLGA2,SWI5` label came from a small boundary overlap; `best_pc_gene` is GOLGA2.

- GOLGA2: <https://www.ncbi.nlm.nih.gov/gene/2801>
- SWI5: <https://www.ncbi.nlm.nih.gov/gene/375757>

## 2. Why the RNA graph includes it

The fresh node has five actual forward `E_r` edges into the chr15 GOLGA6/8 block:

| chr15 endpoint | RNA identity | coverage of shorter exon-sum |
|---|---:|---:|
| `20279485-20289994` | 0.7431 | 0.7718 |
| `20841196-20854572` | 0.7614 | 0.8342 |
| `20972619-20985986` | 0.7611 | 0.8342 |
| `26121347-26134747` (GOLGA8F) | 0.7987 | 0.5113 |
| `26296281-26309739` (GOLGA8G) | 0.7979 | 0.5113 |

This is redundant, forward, major-transcript homology—not one reverse-complement repeat edge
or a single articulation bridge. The core-duplicon literature independently reports that the
GOLGA6 and GOLGA8 subfamilies show their closest similarity to chromosome-9 GOLGA2, their
ancestral source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC6920530/>.

## 3. DNA/RNA certificate

On the same 19 emitted copies, `--joint-dna-rna` reports the five chr9–chr15 RNA edges above
as `RNA_ONLY`. It also reports one `DNA_ONLY` edge between GOLGA8T and GOLGA2. Inspection of
that genomic-span record gives identity 0.8420 and coverage 0.7636. Thus the default 0.60 DNA
certificate also sees shared genomic sequence, but it is below the DNA detector's ordinary
0.90 recent-segmental-duplication floor.

This supports a typed hierarchy:

1. **Broad homology family:** sensitive forward RNA homology; includes GOLGA2.
2. **Recent-copy subfamily:** high-identity DNA duplication structure (ordinary DNA floor
   0.90), with RNA expression placed on those loci; excludes GOLGA2 from the chr15 copy array.

The DNA/RNA intersection should not replace O1 membership globally. Existing joint-arm
experiments show that such a hard intersection loses genuine partial transcripts. Here it is
useful as a *subfamily/type certificate*.

## 4. Direct 0.80 RNA counterfactual

The current binary was rerun on the identical HSA regional BAM with
`--rna-forward-only --min-identity 0.80`.

- all five GOLGA2 RNA edges disappear (maximum identity 0.7987);
- the chr15 GOLGA6/8 block remains intact with 18 copies and 81 edges;
- the complete regional catalog shrinks from 165 to 151 family copies;
- independently named HSA targets fall from 52 to 50 emitted and from 49 to 46 in their
  modal family: two MAGEA members are lost and NBPF splits further.

Therefore 0.80 is a useful **GOLGA recent-subfamily view**, but a rejected global RNA purity
rule. It solves this label only by changing the family granularity and causes measured recall
damage elsewhere.

## 5. Implemented audit correction

The expanded and fresh audits now use three truth classes rather than conflating two:

- target GOLGA6/8 copy;
- `REVIEW_RELATED_OUTGROUP` (GOLGA2, purple); and
- independently named unrelated conflict (orange).

The purified recent-copy GFA excludes the purple outgroup; the audit/fresh broad-family GFA
retains it. No coordinate blacklist and no annotation-dependent Rustle membership rule was
added.

Production emission of both hierarchy levels is deliberately deferred. The proposed algorithm,
CLI/output contract, safeguards, and acceptance tests are specified in
[`docs/o1_duplication_provenance_model.md`](../docs/o1_duplication_provenance_model.md).

## Reproduction

First rebuild the fresh regional substrate:

```bash
python3 bench/o1_fresh_emission_validation.py
```

Then run the two counterfactuals on that unchanged HSA panel BAM:

```bash
RUSTLE_ER_EDGE_DUMP=/tmp/rustle_o1_fresh_emission/HSA.id80_er \
  /tmp/rustle_o1_target/release/gw_family_catalog \
  --bam /tmp/rustle_o1_fresh_emission/HSA.panel.bam \
  --fasta /mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa \
  --threads 4 --rna-forward-only --min-identity 0.80 \
  --out /tmp/rustle_o1_fresh_emission/HSA.id80

RUSTLE_ER_EDGE_DUMP=/tmp/rustle_o1_fresh_emission/HSA.joint_er \
  /tmp/rustle_o1_target/release/gw_family_catalog \
  --bam /tmp/rustle_o1_fresh_emission/HSA.panel.bam \
  --fasta /mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa \
  --threads 4 --rna-forward-only --joint-dna-rna \
  --out /tmp/rustle_o1_fresh_emission/HSA.joint
```
