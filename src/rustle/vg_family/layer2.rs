//! Layer-2 orchestrator (C2–C6). Reads Layer-1 splice graphs and the C1
//! side-index; never re-bundles, never re-splits, never shifts a coordinate.
//! Default-off behind `config.vg_layer2`.

use crate::genome::GenomeIndex;
use crate::graph::Graph;
use crate::path_extract::Transcript;
use crate::types::{DetHashMap, DetHashSet};
use crate::vg::FamilyGroup;
use crate::vg_family::family_graph::FamilyGraph;
use crate::vg_family::secondary_index::SecondaryIndex;
use anyhow::Result;

/// Clustering thresholds passed through to the family-graph merge (conservative
/// defaults matching the bundle-path `build_family_graph` invocation).
const FG_MIN_POS_RECIP: f64 = 0.0;
const FG_REFINE_MIN_JACCARD: f64 = 0.0;

/// Result of a Layer-2 pass for one chromosome.
#[derive(Debug, Default)]
pub struct Layer2Output {
    pub families: Vec<FamilyGroup>,
    pub family_graphs: Vec<Option<FamilyGraph>>,
    pub novel_transcripts: Vec<Transcript>,
}

/// One per-locus input to Layer 2: the Layer-1 splice graph + its chrom/strand.
pub struct Layer1Locus<'a> {
    pub chrom: String,
    pub strand: char,
    pub graph: &'a Graph,
}

/// Run Layer 2 over one chromosome's Layer-1 loci.
///
/// `primary_locus[read_name_hash] = bundle index of the read's primary`. The
/// side-index is consumed (pruned to family-candidate loci, capped per locus).
#[allow(clippy::too_many_arguments)]
pub fn run_layer2(
    loci: &[Layer1Locus<'_>],
    mut side_index: SecondaryIndex,
    primary_locus: &DetHashMap<u64, usize>,
    genome: Option<&GenomeIndex>,
    min_link: u32,
    min_similarity: f64,
    per_locus_cap: usize,
    kmer_len: usize,
    new_copies: bool,
) -> Result<Layer2Output> {
    // (C2) candidate links → compute similarity only for linked, same-chrom pairs.
    let links = crate::vg::build_multimap_index_from_secondary_index(&side_index, primary_locus);
    let mut similarity: DetHashMap<(usize, usize), f64> = DetHashMap::default();
    if let Some(g) = genome {
        for ((a, b), count) in &links {
            if *count < min_link {
                continue;
            }
            let (la, lb) = (&loci[*a], &loci[*b]);
            if la.chrom == lb.chrom {
                let sim = crate::vg::exon_kmer_similarity_between_graphs(
                    la.graph, lb.graph, &la.chrom, g, kmer_len,
                );
                similarity.insert((*a, *b), sim);
            }
        }
    }

    let families = crate::vg::discover_family_groups_layer2(
        &side_index, primary_locus, &similarity, min_link, min_similarity,
    );

    // (C1 prune) keep only family-candidate loci, cap per locus (logged drops).
    let keep: DetHashSet<usize> = families
        .iter()
        .flat_map(|f| f.bundle_indices.iter().copied())
        .collect();
    let pruned = side_index.prune_to_loci(&keep);
    let capped = side_index.cap_per_locus(per_locus_cap);
    if pruned + capped > 0 {
        eprintln!(
            "[layer2] side-index pruned {pruned} (non-candidate) + capped {capped} (per-locus>{per_locus_cap})"
        );
    }

    // (C3) merge each family's Layer-1 graphs into a variation graph.
    let mut family_graphs: Vec<Option<FamilyGraph>> = Vec::with_capacity(families.len());
    for fam in &families {
        let copies: Vec<(String, char, &Graph)> = fam
            .bundle_indices
            .iter()
            .map(|&i| (loci[i].chrom.clone(), loci[i].strand, loci[i].graph))
            .collect();
        let fg = crate::vg_family::family_graph::build_family_graph_from_layer1_graphs(
            fam.family_id,
            &copies,
            genome,
            FG_MIN_POS_RECIP,
            min_similarity,
            FG_REFINE_MIN_JACCARD,
        )
        .ok();
        family_graphs.push(fg);
    }

    let _ = new_copies; // consumed in M5
    Ok(Layer2Output {
        families,
        family_graphs,
        novel_transcripts: Vec::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_layer2_forms_family_and_merges_graph() {
        let g0 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160), (300, 360)],
        );
        let g1 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160), (5300, 5360)],
        );
        let genome = crate::vg_family::family_graph::tests_support::make_two_copy_genome();

        let loci = vec![
            Layer1Locus { chrom: "chrT".into(), strand: '+', graph: &g0 },
            Layer1Locus { chrom: "chrT".into(), strand: '+', graph: &g1 },
        ];

        let mut si = SecondaryIndex::new();
        si.push(crate::vg_family::secondary_index::SecondaryAlignment {
            read_name_hash: 7,
            chrom: "chrT".into(),
            ref_start: 100,
            ref_end: 160,
            introns: vec![],
            nm: 0,
            is_supplementary: false,
            locus: Some(1),
        });
        let mut primary_locus: DetHashMap<u64, usize> = DetHashMap::default();
        primary_locus.insert(7, 0);

        let out = run_layer2(
            &loci, si, &primary_locus, Some(&genome),
            1, 0.3, 1000, 15, false,
        )
        .unwrap();

        assert_eq!(out.families.len(), 1, "one family from C2");
        assert!(out.family_graphs[0].is_some(), "family graph merged (C3)");
        assert!(out.novel_transcripts.is_empty(), "no emission yet (M3)");
    }
}
