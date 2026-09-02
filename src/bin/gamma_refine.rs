//! `gamma_refine` — apply the SHIPPED gamma-quasi-clique refinement to an arbitrary edge list.
//!
//! Exists so a candidate graph built outside the RNA pipeline (e.g. the annotation-seeded DNA
//! projection, ledger §6be) can be split by the SAME rule the read catalog uses, rather than by a
//! re-implementation. It calls `family_definition::refine_component` directly — connected
//! components, then density gate `2|E|/(n(n-1)) >= GAMMA`, then the deterministic splitter,
//! recursively — so any comparison of the two catalogs varies the EDGES, not the partitioner.
//!
//! stdin : one edge per line, `u<TAB>v`, node ids are arbitrary strings.
//! stdout: one refined block per line, `block_id<TAB>member` rows.
use std::collections::{BTreeSet, HashMap};
use std::io::{BufRead, Write};

use rustle::vg_family::family_definition::{induced_density, refine_component, GAMMA};

fn main() -> anyhow::Result<()> {
    let gamma: f64 = std::env::var("RUSTLE_GAMMA").ok().and_then(|v| v.parse().ok()).unwrap_or(GAMMA);
    let mut id_of: HashMap<String, usize> = HashMap::new();
    let mut name: Vec<String> = Vec::new();
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for line in std::io::stdin().lock().lines() {
        let line = line?;
        let mut it = line.split('\t');
        let (Some(a), Some(b)) = (it.next(), it.next()) else { continue };
        let mut intern = |s: &str| -> usize {
            *id_of.entry(s.to_string()).or_insert_with(|| {
                name.push(s.to_string());
                name.len() - 1
            })
        };
        let (u, v) = (intern(a), intern(b));
        if u != v {
            edges.push((u, v));
        }
    }
    let n = name.len();
    let mut adj: HashMap<usize, BTreeSet<usize>> = HashMap::new();
    for &(u, v) in &edges {
        adj.entry(u).or_default().insert(v);
        adj.entry(v).or_default().insert(u);
    }
    // connected components of the whole graph = the raw families the refinement starts from
    let mut seen = vec![false; n];
    let mut comps: Vec<BTreeSet<usize>> = Vec::new();
    for s in 0..n {
        if seen[s] {
            continue;
        }
        let mut stack = vec![s];
        let mut c = BTreeSet::new();
        while let Some(x) = stack.pop() {
            if !seen[x] {
                seen[x] = true;
                c.insert(x);
                for &y in adj.get(&x).into_iter().flatten() {
                    if !seen[y] {
                        stack.push(y);
                    }
                }
            }
        }
        comps.push(c);
    }
    let out = std::io::stdout();
    let mut w = std::io::BufWriter::new(out.lock());
    let (mut blocks, mut split) = (0usize, 0usize);
    for c in &comps {
        let refined = refine_component(c, &adj, gamma);
        if refined.len() > 1 {
            split += 1;
        }
        for b in refined {
            for m in &b {
                writeln!(w, "B{}\t{}", blocks, name[*m])?;
            }
            blocks += 1;
        }
    }
    eprintln!(
        "[gamma_refine] gamma={gamma} nodes={n} edges={} components={} -> blocks={blocks} \
         (components split: {split})",
        edges.len(),
        comps.len()
    );
    let _ = induced_density; // referenced for parity with the shipped density definition
    Ok(())
}
