//! `mcl_refine` — apply the SHIPPED MCL clustering (`annotation_families::mcl`) to an arbitrary
//! weighted edge list. Measurement scaffolding for the §6iz proposal #4 controlled ablation
//! (`docs/o1_ledger.md`): mirrors `gamma_refine`'s existing shape exactly (same stdin contract, same
//! "vary the edges/partitioner independently" purpose) so a candidate graph built outside the MCL
//! pipeline (e.g. the de-novo/E_r homology graph) can be partitioned by the SAME rule the annotated
//! catalog uses, instead of a re-implementation. Not wired into any shipped pipeline.
//!
//! stdin : one edge per line, `u<TAB>v<TAB>weight` (weight in [0,1]; extra columns ignored), node ids
//!         are arbitrary strings.
//! stdout: one cluster block per line, `cluster_id<TAB>member` rows.
use std::collections::BTreeMap;
use std::io::{BufRead, Write};

use rustle::vg_family::annotation_families::{mcl, HomologyGraph};

fn main() -> anyhow::Result<()> {
    let inflation: f64 =
        std::env::var("RUSTLE_MCL_INFLATION").ok().and_then(|v| v.parse().ok()).unwrap_or(2.8);
    let prune: f64 = std::env::var("RUSTLE_MCL_PRUNE").ok().and_then(|v| v.parse().ok()).unwrap_or(1e-9);

    let mut id_of: BTreeMap<String, usize> = BTreeMap::new();
    let mut name: Vec<String> = Vec::new();
    let mut edges: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    for line in std::io::stdin().lock().lines() {
        let line = line?;
        let mut it = line.split('\t');
        let (Some(a), Some(b), Some(w)) = (it.next(), it.next(), it.next()) else { continue };
        let w: f64 = w.parse().unwrap_or(0.0);
        let mut intern = |s: &str| -> usize {
            if let Some(&i) = id_of.get(s) {
                return i;
            }
            let i = name.len();
            name.push(s.to_string());
            id_of.insert(s.to_string(), i);
            i
        };
        let (u, v) = (intern(a), intern(b));
        if u != v {
            let (u, v) = if u < v { (u, v) } else { (v, u) };
            let cur = edges.entry((u, v)).or_insert(0.0);
            if w > *cur {
                *cur = w;
            }
        }
    }
    let n = name.len();
    let g = HomologyGraph {
        genes: (0..n).map(|_| (String::new(), 0u64, 0u64)).collect(),
        edges,
        ..Default::default()
    };
    let parts = mcl(&g, inflation, 100, prune);

    let out = std::io::stdout();
    let mut w = std::io::BufWriter::new(out.lock());
    for (bi, block) in parts.iter().enumerate() {
        for &m in block {
            writeln!(w, "B{bi}\t{}", name[m])?;
        }
    }
    eprintln!(
        "[mcl_refine] inflation={inflation} prune={prune} nodes={n} edges={} -> clusters={}",
        g.n_edges(),
        parts.len()
    );
    Ok(())
}
