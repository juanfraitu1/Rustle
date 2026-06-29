# Graph vs linear copy-assignment — RCF_611 simulation (decisive test)

User hypothesis: a vg graph (modeling the tandem repeat as a cycle) aligns reads to copies BETTER than
linear alignment, because it places reads correctly through the repeat instead of mis-anchoring.

Test: RCF_611 = 4 co-located paralogs (NC_073247.2, ~2 Mb span) with tandem-repeat CNV (cDNAs
599–960 bp, pairwise id 0.857–0.996). 800 KNOWN-ORIGIN reads simulated (200/copy, full-length +
fragments, ~1.2% HiFi error). Assign by LINEAR pairwise (minimap2 -N50 -p0.1, primary target) vs
GRAPH (vg giraffe to a `vg msga` graph of the 4 copies; copy = which copy's UNIQUE graph nodes the
read traverses). `bench/rcf611_graph_vs_linear.py`.

## Result — hypothesis REFUTED

| method | correct | wrong | ambiguous |
|---|--:|--:|--:|
| LINEAR-primary (pipeline) | **98.0%** | 0.9% | 1.1% |
| LINEAR-ASdecisive | 93.4% | 0.0% | 6.6% |
| GRAPH (vg giraffe) | 94.5% | 0.1% | 5.4% |

- The graph **rescues 0** reads linear got wrong; it regresses 1.
- On linear's 13 MAPQ0 reads, GRAPH gets **2/13** correct vs linear's **6/13** (graph WORSE on the hard reads).
- Linear (minimap2) is already 98% accurate.

## Why
- minimap2 already handles these tandem repeats (chaining/scoring is not fooled here) — the earlier
  "linear is mangled" stat was a SECONDARY-alignment artifact (secondaries are MAPQ0 + clipped by
  definition); primary reads are mostly MAPQ>0 and correct.
- The copy-distinguishing signal is in the UNIQUE regions, which BOTH methods align fine. The shared
  tandem repeat carries NO copy information, so collapsing it into a graph cycle adds no discriminating
  power. Reads inside shared sequence stay ambiguous either way.
- Same disjoint-levers result the project keeps finding: copy resolution lives in sequence (unique
  regions); the structural/graph representation is faithful but not a resolution lever.

## Verdict
The vg graph's non-linear (repeat) modeling is real biology but is NOT a copy-assignment improvement
here. Do not graph-align for copy resolution; linear alignment + AS-decisive copy choice is sufficient.
(One family; representative co-located tandem-CNV case with a real MAPQ0 tail. A harder regime — many
near-identical copies, very long arrays — could differ, but RCF_611 shows no benefit.)
