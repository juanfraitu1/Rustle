# PREREG r1527: read-level junction snap (10 bp) for one-chain family misses

Registered 2026-11-04 after step 0 and before any arm number (user: "lets register the lever and test").

## Step 0 (measured before this file)
On the 97 one-chain members (r1526 class), junction support = same-strand spliced primaries within 50 kb of the gene (what the
pass-1 region sees): the existing rule of `copy_assign --read-junction-snap W` (a junction moves onto the junction within W bp
at donor and acceptor that strictly more reads carry, >= 2) makes >= 2 chains identical at **23 of 97** for W = 10, and 23 for
W = 20 (~/r1526/step0_snap.py). Ceiling 23. History: r1140 (2026-09-28, polish-level fuzzy merges) found near pairs are mostly
NAGNAG / real alternative sites, and support-ratio merges cost 5-7% of chains; the read-level snap was listed as an opt-in
lever (tools/rustle_pipeline.sh RUSTLE_READ_JUNCTION_SNAP) and never registered.

## Fixed before any number
- Arm S = one genome-wide `rustle_pipeline.sh assemble`, RUSTLE_ASSEMBLY_LAYERS=off RUSTLE_READ_JUNCTION_SNAP=10, strand vote
  default, binary ~/rustle_target_test/release (unchanged since r1526). Base = fusedgw/R1. Wall bar 45 min on the assemble.
- Scores as r1526: family_panel.py (T), chain_cat.py (a base chain whose junction moved counts as NOT kept: the precision
  guard r1140 asked for), score34.py --one-chain.
- Bars (all must hold): complete families >= 214 (base 209; half the ceiling); fused <= 49; ordinary hit rate >= 0.829;
  base spliced chains kept >= 0.99; new chains' exact CAT match >= 0.05; of the 97, >= 12 gain a same-strand model (half of 23).
- PASS -> RUSTLE_READ_JUNCTION_SNAP=10 becomes the assemble default; any FAIL -> stays opt-in, row closed.
