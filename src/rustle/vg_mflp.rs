//! MFLP (Minimum Flow Linear Program) solver for multi-mapping read assignment.
//!
//! Uses `good_lp` (already in Cargo.toml) to formulate a linear program that
//! assigns multi-mapping reads to gene family copies optimally.
//!
//! For each multi-mapping read r across copies k₁..kₙ:
//! - Variable: w[r][k] ∈ [0, 1] — weight at copy k
//! - Constraint: Σ_k w[r][k] = 1 — weights sum to 1
//! - Objective: maximize Σ_{r,k} w[r][k] * compatibility(r, k)
//!
//! Activated with `--vg --vg-solver mflp`.

use crate::types::Bundle;
use crate::vg::{junction_compatibility, EmResult, FamilyGroup};
use good_lp::{constraint, default_solver, variable, variables, Expression, Solution, SolverModel};
use std::collections::HashMap;

/// Solve the multi-mapping read assignment LP for one family group.
///
/// Returns updated weights: read_name_hash → Vec<(global_bundle_idx, weight)>.
pub fn solve_mflp_family(
    family: &FamilyGroup,
    bundles: &[Bundle],
) -> (HashMap<u64, Vec<(usize, f64)>>, EmResult) {
    let mut result = EmResult::default();
    let mut weight_map: HashMap<u64, Vec<(usize, f64)>> = HashMap::new();

    if family.multimap_reads.is_empty() {
        return (weight_map, result);
    }

    // Collect reads and their compatibility scores.
    struct ReadInfo {
        rnh: u64,
        locs: Vec<(usize, usize)>, // (global_bundle_idx, read_idx)
        scores: Vec<f64>,           // compatibility per copy
    }
    let mut read_infos: Vec<ReadInfo> = Vec::new();

    for (&rnh, locs) in &family.multimap_reads {
        let mut info = ReadInfo {
            rnh,
            locs: Vec::new(),
            scores: Vec::new(),
        };
        for &(fam_pos, ri) in locs {
            let global_bi = family.bundle_indices[fam_pos];
            if ri >= bundles[global_bi].reads.len() {
                // Supplementary-only link: give small base score.
                info.locs.push((global_bi, ri));
                info.scores.push(0.1);
                continue;
            }
            let read = &bundles[global_bi].reads[ri];
            let bundle = &bundles[global_bi];
            let compat = junction_compatibility(read, bundle);
            let context: f64 = bundle
                .junction_stats
                .iter()
                .map(|(_, st)| st.nreads_good)
                .sum::<f64>()
                .max(1.0);
            info.locs.push((global_bi, ri));
            info.scores.push((compat + 0.01) * context.ln().max(1.0));
        }
        if info.locs.len() >= 2 {
            read_infos.push(info);
        }
    }

    if read_infos.is_empty() {
        return (weight_map, result);
    }

    // Build LP.
    let mut vars = variables!();

    // Create variables: w[read_idx][copy_idx]
    let mut w: Vec<Vec<good_lp::Variable>> = Vec::new();
    for info in &read_infos {
        let mut read_vars = Vec::new();
        for _ in &info.locs {
            read_vars.push(vars.add(variable().min(0.0).max(1.0)));
        }
        w.push(read_vars);
    }

    // Objective: maximize sum of score * weight.
    let mut objective = Expression::from(0.0);
    for (ri, info) in read_infos.iter().enumerate() {
        for (ci, &score) in info.scores.iter().enumerate() {
            objective = objective + score * w[ri][ci];
        }
    }

    let mut problem = vars.maximise(objective).using(default_solver);

    // Constraints: weights sum to 1 for each read.
    for (ri, info) in read_infos.iter().enumerate() {
        let mut sum = Expression::from(0.0);
        for ci in 0..info.locs.len() {
            sum = sum + w[ri][ci];
        }
        problem = problem.with(constraint!(sum == 1.0));
    }

    // Solve.
    match problem.solve() {
        Ok(solution) => {
            let mut n_reweighted = 0usize;
            for (ri, info) in read_infos.iter().enumerate() {
                let mut weights: Vec<(usize, f64)> = Vec::new();
                for (ci, &(global_bi, _ri)) in info.locs.iter().enumerate() {
                    let new_w = solution.value(w[ri][ci]);
                    weights.push((global_bi, new_w));
                }
                // Normalize (LP solution may not exactly sum to 1 due to floating point).
                let total: f64 = weights.iter().map(|(_, w)| w).sum();
                if total > 0.0 {
                    for (_, w) in &mut weights {
                        *w /= total;
                    }
                }
                n_reweighted += weights.len();
                weight_map.insert(info.rnh, weights);
            }
            result.iterations = 1;
            result.converged = true;
            result.max_delta = 0.0;
            result.reads_reweighted = n_reweighted;
        }
        Err(e) => {
            eprintln!(
                "[VG] MFLP solver failed for family {}: {}",
                family.family_id, e
            );
        }
    }

    (weight_map, result)
}

/// Run MFLP for all families, applying weights to bundles.
pub fn run_mflp(
    families: &[FamilyGroup],
    bundles: &mut [Bundle],
) -> Vec<EmResult> {
    let mut results = Vec::with_capacity(families.len());

    for family in families {
        let (weight_map, result) = solve_mflp_family(family, bundles);

        // Apply weights.
        for (rnh, weights) in &weight_map {
            for &(global_bi, new_w) in weights {
                // Find the read in this bundle.
                if let Some(ri) = bundles[global_bi]
                    .reads
                    .iter()
                    .position(|r| r.read_name_hash == *rnh)
                {
                    bundles[global_bi].reads[ri].weight = new_w;
                }
            }
        }

        if result.reads_reweighted > 0 {
            eprintln!(
                "[VG] Family {} MFLP: reweighted {} read entries across {} copies",
                family.family_id,
                result.reads_reweighted,
                family.bundle_indices.len(),
            );
        }

        results.push(result);
    }

    let total: usize = results.iter().map(|r| r.reads_reweighted).sum();
    if total > 0 {
        eprintln!(
            "[VG] MFLP complete: {} reads adjusted across {} families",
            total,
            results.iter().filter(|r| r.reads_reweighted > 0).count(),
        );
    }

    results
}
