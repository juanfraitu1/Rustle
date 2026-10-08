# bridge_units fixture: the units execution against its dev prototype

A synthetic PLAIN GTF (the shape `copy_assign --assemble-only --bridge-regroup off` writes: one `gene_id` per
junction-sharing component) with every case the unit split and the scoped native regroup must get right. The
`expected.*` files are the products of the dev prototype, not of this repository's Rust; the tests assert that the Rust
output equals them (after removing the attributes the prototype does not write: `fusion_locus`, `fusion_gene`,
`fusion_detector`, `fusion_evidence`).

| file | what it is |
|---|---|
| `plain.gtf` | 27 transcripts in 10 genes (see the cases below) |
| `cuts.tsv` | the list of the 5 transcripts to cut and their introns (`--bridge-units-list`; also the prototype's cut file) |
| `expected.units.gtf` | prototype `units2.py --regroup scoped --attach units --mode oracle --cuts cuts.tsv --label list:cuts.tsv` (9fea69a1) |
| `expected.units.tsv` | its units table (`expected.units.gtf.units.tsv`; the Rust table has no `oracle_label` column and ends with `evidence`) |
| `clusters.tsv`, `loci.tsv` | hand-written families products over the loci of `expected.units.gtf` (clusters MCL0..MCL5 and one fold) |
| `expected.relations.tsv`, `expected.members_by_locus.tsv` | prototype `relations.py --detector list:cuts.tsv` (d90a33da) on the above |

`s_f0.5.plain.gtf`, `s_f0.5.cuts.tsv`, `s_f0.5.expected.units.gtf`, `s_f0.5.expected.units.tsv` are the same test on REAL gorilla
loci: 34 transcripts of 6 `gene_id`s of the fusion simulation's f = .5 plain assembly (the four genes that hold an F1 bridge, three on
`-` and one on `+`, and the two smallest others), the four bridges F1 found (`R.cuts.tsv` of the units study), and the prototype's
scoped output for them (`units2.py ... --label f1`; its table's `label` column differs from ours by construction).

Cases (gene: what it tests): `g1` (+, m = 1, a third exon-disjoint component `g1.nat3`, a second isoform), `g2` (-, m = 2: units
in transcription order, genomic right to left), `g3` / `g4` (a single-exon unit attaches to the untouched gene `g4`'s piece,
which keeps its RG3 name; `g4` also holds an exon-disjoint piece `g4.rg2`), `g3` / `g5` / `g6` (a single-exon unit with equal
overlap on two transcripts goes to the lower line: `g5`), `g8` (untouched gene whose two strands share a junction: RG3 splits it,
the native rule would not), `g9` (untouched, a single-exon transcript overlapping a spliced one), `g10` (+, m = 3), `g11` (contig
`c2` with `g1`'s coordinates: the junction key holds the contig).

Regenerate: `make_fixture.py` (scratch), then the two prototype commands above. The prototype scripts live in the scratch
directory `container_units_v2/lib/` of the 2026-09-30 units study (`docs/archive/2026-09/PREREG_container_units_v2_dev_2026-09-30.md`).
