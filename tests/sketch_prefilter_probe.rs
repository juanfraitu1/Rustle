//! Sketch-prefilter probe (measurement harness, not a shipped test).
//!
//! Question: for the `mcl_families --from-gtf` all-vs-all (one batched `minimap2 -x asm20 -c -X` over
//! `<out>.loci.fa`, edges filtered to identity >= 0.7, cov_longer >= 0.3, >=300 bp), can the
//! candidates stage's minimizer sketch (`rustle::candidate_copies::detect::minimizer_sketch`, k=31 w=5,
//! containment = `sketch_share`) predict which locus pairs become edges — and how many pairs it would
//! let a prefilter drop before minimap2?
//!
//! Run with:
//!   SKETCH_PROBE_LOCI=/path/to/slice.fam.loci.fa SKETCH_PROBE_GRAPH=/path/to/slice.fam.graph.tsv \
//!     cargo test --profile dev-opt --test sketch_prefilter_probe -- --ignored --nocapture
//!
//! The graph file comes from `mcl_families --dump-graph` on the SAME loci (edge rows:
//! `chrom:start-end<TAB>chrom:start-end<TAB>weight`, headers matching the FASTA verbatim).

use std::collections::HashMap;

fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read_to_string(path).expect("read loci fasta");
    let mut out = Vec::new();
    let mut name: Option<String> = None;
    let mut seq = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some(n) = name.take() {
                out.push((n, std::mem::take(&mut seq)));
            }
            name = Some(h.split_whitespace().next().unwrap_or("").to_string());
        } else {
            seq.extend(line.trim().as_bytes());
        }
    }
    if let Some(n) = name.take() {
        out.push((n, seq));
    }
    out
}

#[test]
#[ignore = "measurement harness: needs SKETCH_PROBE_LOCI and SKETCH_PROBE_GRAPH"]
fn sketch_containment_vs_graph_edges() {
    let loci_path = std::env::var("SKETCH_PROBE_LOCI").expect("SKETCH_PROBE_LOCI");
    let graph_path = std::env::var("SKETCH_PROBE_GRAPH").expect("SKETCH_PROBE_GRAPH");
    let loci = read_fasta(&loci_path);
    let n = loci.len();
    eprintln!("[probe] {n} loci");

    // Sketches (k=31, w=5 — the candidates stage's constants).
    let t0 = std::time::Instant::now();
    let sketches: Vec<Vec<u64>> = loci
        .iter()
        .map(|(_, s)| rustle::candidate_copies::detect::minimizer_sketch(s, 31, 5))
        .collect();
    eprintln!(
        "[probe] sketches in {:.2?}; sketch sizes: min {} median {} max {}",
        t0.elapsed(),
        sketches.iter().map(|s| s.len()).min().unwrap_or(0),
        sketches.iter().map(|s| s.len()).nth(n / 2).unwrap_or(0),
        sketches.iter().map(|s| s.len()).max().unwrap_or(0)
    );

    // Ground-truth edges from --dump-graph.
    let idx: HashMap<&str, usize> = loci
        .iter()
        .enumerate()
        .map(|(i, (nm, _))| (nm.as_str(), i))
        .collect();
    let mut edges = std::collections::HashSet::new();
    for line in std::fs::read_to_string(&graph_path)
        .expect("read graph")
        .lines()
    {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let (Some(&a), Some(&b)) = (idx.get(f[0]), idx.get(f[1])) else {
            continue;
        };
        edges.insert(if a < b { (a, b) } else { (b, a) });
    }
    eprintln!("[probe] {} edges in graph", edges.len());

    // All-pairs containment; collect per-pair share and edge label.
    let t0 = std::time::Instant::now();
    let mut edge_shares: Vec<f64> = Vec::new();
    let mut nonedge_shares: Vec<f64> = Vec::new();
    let thresholds = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8];
    let mut pairs_passing = vec![0usize; thresholds.len()];
    let mut edges_kept = vec![0usize; thresholds.len()];
    for i in 0..n {
        for j in (i + 1)..n {
            let s = rustle::candidate_copies::detect::sketch_share(&sketches[i], &sketches[j]);
            let is_edge = edges.contains(&(i, j));
            if is_edge {
                edge_shares.push(s);
            } else {
                nonedge_shares.push(s);
            }
            for (k, &t) in thresholds.iter().enumerate() {
                if s >= t {
                    pairs_passing[k] += 1;
                    if is_edge {
                        edges_kept[k] += 1;
                    }
                }
            }
        }
    }
    let total_pairs = n * (n - 1) / 2;
    eprintln!("[probe] all-pairs containment in {:.2?}", t0.elapsed());

    let summary = |v: &mut Vec<f64>| {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if v.is_empty() {
            return "empty".to_string();
        }
        format!(
            "n={} min={:.4} p50={:.4} p90={:.4} p99={:.4} max={:.4}",
            v.len(),
            v[0],
            v[v.len() / 2],
            v[(v.len() * 9) / 10],
            v[(v.len() * 99) / 100],
            v[v.len() - 1]
        )
    };
    println!("edges:     {}", summary(&mut edge_shares));
    println!("non-edges: {}", summary(&mut nonedge_shares));
    println!();
    println!(
        "{:>8} {:>14} {:>14} {:>12} {:>12}",
        "thresh", "pairs_passing", "pair_reduction", "edges_kept", "edge_recall"
    );
    for (k, &t) in thresholds.iter().enumerate() {
        let recall = edges_kept[k] as f64 / edges.len().max(1) as f64;
        println!(
            "{:>8.3} {:>14} {:>13.1}x {:>12} {:>11.1}%",
            t,
            pairs_passing[k],
            total_pairs as f64 / pairs_passing[k].max(1) as f64,
            edges_kept[k],
            100.0 * recall
        );
    }
    // The verdict: the smallest threshold keeping 100% of edges, and its pair reduction.
    let full_recall_t = edges_kept
        .iter()
        .position(|&e| e == edges.len())
        .map(|k| thresholds[k]);
    if let Some(t) = full_recall_t {
        let k = thresholds.iter().position(|&x| x == t).unwrap();
        println!(
            "\n100% edge recall first holds at threshold {} ({} of {} pairs kept, {:.1}x fewer)",
            t,
            pairs_passing[k],
            total_pairs,
            total_pairs as f64 / pairs_passing[k].max(1) as f64
        );
    } else {
        println!("\nNO threshold keeps 100% of edges — a sketch prefilter would lose edges.");
    }

    // Single-linkage components at the full-recall threshold (and one notch up): if the
    // loci collapse into one giant repeat-linked component, a partition-then-align prefilter
    // cannot reduce the AVA; if they break into small groups, it can.
    for &t in &[0.005f64, 0.01, 0.05] {
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(p: &mut Vec<usize>, mut x: usize) -> usize {
            while p[x] != x {
                p[x] = p[p[x]];
                x = p[x];
            }
            x
        }
        let mut n_union = 0usize;
        for i in 0..n {
            for j in (i + 1)..n {
                if rustle::candidate_copies::detect::sketch_share(&sketches[i], &sketches[j]) >= t {
                    let (a, b) = (find(&mut parent, i), find(&mut parent, j));
                    if a != b {
                        parent[a.max(b)] = a.min(b);
                        n_union += 1;
                    }
                }
            }
        }
        let mut comp: HashMap<usize, usize> = HashMap::new();
        for i in 0..n {
            *comp.entry(find(&mut parent, i)).or_insert(0) += 1;
        }
        let mut sizes: Vec<usize> = comp.into_values().collect();
        sizes.sort_unstable_by(|a, b| b.cmp(a));
        let pairs_within: usize = sizes.iter().map(|&s| s * (s - 1) / 2).sum();
        println!(
            "t={:.3}: {} components, top sizes {:?}; AVA pairs within components = {} ({:.1}x fewer than {total_pairs})",
            t,
            sizes.len(),
            sizes.iter().take(6).collect::<Vec<_>>(),
            pairs_within,
            total_pairs as f64 / pairs_within.max(1) as f64,
        );
        let _ = n_union;
    }
}
