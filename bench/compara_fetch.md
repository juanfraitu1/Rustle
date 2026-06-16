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

