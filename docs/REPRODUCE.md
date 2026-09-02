# One catalog, reproducible today

> **Why this file exists.** Every catalog count the project quoted before 2026-08-30 is
> **unreproducible by the current binary** — `refine` was removed on 08-20 and the O1 node floor moved
> 3 → 2 on 08-30. The dossier's advice follows from that: *a number he cannot reproduce is a number
> he will assume was chosen.* This file pins **one** catalog end to end so that objection has a
> concrete answer. It is deliberately small: three contigs, ~40 min, one command.

## The pinned run

| | |
|---|---|
| **source** | `b16a019` on `dna-from-genome` (last commit touching `src/`; **0** src commits since) |
| **binary** | `gw_family_catalog`, md5 `14a850a9d446d0a2dfc65fca1ef6bbf9` |
| **build** | `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin gw_family_catalog` |
| **reads** | `npip3.bam` — 1,706,598,583 B, md5 `c1ad3cda53282bc3d78f4876fecf179c` — gorilla fibroblast IsoSeq on mGorGor1, 3 contigs |
| **genome** | `GGO.fasta` — 3,590,176,889 B, md5 `7a26fb029571b9e9450179769b2dc969` |
| **env** | `RUSTLE_GATE_MIN_READS=2`; every other `RUSTLE_*` **unset** (the arm script unsets them explicitly) |

```bash
export RUSTLE_GATE_MIN_READS=2
export RUSTLE_ER_EDGE_DUMP=$OUT/dump/e
gw_family_catalog --bam npip3.bam --fasta GGO.fasta --out $OUT/cat \
                  --homology-primary --threads 4
```

## What it must produce

| artifact | md5 | count |
|---|---|---|
| `cat.copies.tsv` | `9849dcb45b63e48e7b9b4d4358113a10` | **678 copies** |
| `cat.families.tsv` | `bac5d98e2b623400a8ad61758f3160ca` | **121 families** |
| `dump/e.nodes.tsv` | `0a3bb20798dc954b533371e38171216a` | **3,598 nodes** |
| `dump/e.edges.tsv` | — | **3,141 edges** |

⭐ **This was verified, not asserted.** `arm_f2` (built 08-30) and `arm_cfoff` (rebuilt 09-02 after
the coverage clause landed) are **md5-identical** on all three files above, and their edge dumps agree
**row-for-row on all 16 shared columns with 0 rows differing either way**. The only difference is
three *added* columns — `cov_longer`, `unaln_i`, `unaln_j` — which §6ba emits as **disclosure, not as
a gate**.

⚠ `dump/e.edges.tsv` has no md5 here precisely because of those three columns: a build predating
§6ba produces the same 16 columns and a different file hash. Compare the columns, not the hash.

## The certificate the run writes about itself

`dump/e.params.tsv` records the rule that produced the catalog, so the parameters do not have to be
taken on trust:

```
min_identity_asm20      0.800000     sensitive_identity   0.600000
min_coverage            0.500000     gate_min_reads       2
coverage_form   aligned span on the SHORTER sequence / len(shorter)
identity_metric 1-de (fallback nmatch/blocklen when de:f: absent)
alignment_orientation   forward-only (+)
edge_rule       ANY single record clearing both floors
mm_args_sensitive       minimap2 -c -X --no-long-join -t 4 -k 11 -w 5
```

followed by an `env.RUSTLE_*` row per behaviour flag, each reading `<unset>` for this run. That block
is the answer to *"which of your 135 knobs were on?"* — **the run says so itself, per run**, and a
flag added later without a certificate row is the M2 defect the project treats as a bug.

## What this file does NOT claim

- ⛔ **Not the genome-wide catalog.** Three contigs. No genome-wide catalog exists from the current
  binary; that run is outstanding (§6bp, §6bs).
- ⛔ **Not a validated catalog** — reproducibility is not correctness. What is reproduced is the
  *computation*, which is the precondition for arguing about the result, not the argument.
- ⚠ **The inputs are not in the repository** (5.3 GB). Their md5s are pinned above, so a
  reproduction can prove it started from the same bytes; obtaining them is out of band.
