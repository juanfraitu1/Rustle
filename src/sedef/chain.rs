//! Stage 2: chain candidate window pairs into SD blocks.
//!
//! Deliberately cheap so it runs comfortably on a laptop:
//!
//! 1. Cells (window-index pairs) are grouped into bands of two diagonals
//!    (`diag = b - a`). Within a band, a duplicated segment appears as a
//!    monotonically increasing run even when the two copies drift out of
//!    grid phase by one window.
//! 2. A single sorted walk per band extends a chain while both coordinates
//!    advance by at most `gap` windows -- this bridges small indel gaps
//!    without any k-mer re-extraction.
//! 3. Chains whose A-intervals AND B-intervals overlap by at least
//!    `merge_frac` of the smaller span are unioned (union-find with an
//!    interval sweep line). This reconnects phase-shifted fragments while
//!    refusing to merge merely-nested distinct duplications.
//!
//! Everything is integer sorting and linear scans: ~O(n log n), low memory.

use std::collections::HashMap;
use std::error::Error;
use std::io::BufRead;

#[derive(Debug, Clone)]
pub struct PairRec {
    pub chr_a: String,
    pub s_a: u64,
    pub e_a: u64,
    pub chr_b: String,
    pub s_b: u64,
    pub e_b: u64,
    pub shared: u32,
}

#[derive(Debug, Clone)]
pub struct Chain {
    pub ca: String,
    pub cb: String,
    pub a0: u32,
    pub a1: u32, // inclusive window indices, side A
    pub b0: u32,
    pub b1: u32, // side B
    pub score: u32,
    pub cells: u32,
}

#[derive(Debug, Clone)]
pub struct Block {
    pub chr_a: String,
    pub a0: u64, // 0-based, half-open genomic coords
    pub a1: u64,
    pub chr_b: String,
    pub b0: u64,
    pub b1: u64,
    pub cells: u32,
    pub score: u32,
    pub est_id: f64,
}

pub struct ChainParams {
    pub window: u64,
    pub gap: u32,
    pub min_score: u32,
    pub merge_frac: f64,
    pub sketch: usize,
    pub k: usize,
    /// Cells below this estimated identity are discarded before chaining.
    pub min_cell_id: f64,
    /// Blocks below this estimated identity are dropped at the end.
    pub min_block_id: f64,
    /// Precomputed: min `shared` count corresponding to `min_cell_id`.
    pub min_shared_cell: u32,
}

pub fn load_pairs(
    path: &std::path::Path,
    min_shared: u32,
) -> Result<Vec<PairRec>, Box<dyn Error>> {
    let f = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(f);
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() < 7 {
            continue;
        }
        let shared: u32 = t[6].parse()?;
        if shared < min_shared {
            continue;
        }
        out.push(PairRec {
            chr_a: t[0].to_string(),
            s_a: t[1].parse()?,
            e_a: t[2].parse()?,
            chr_b: t[3].to_string(),
            s_b: t[4].parse()?,
            e_b: t[5].parse()?,
            shared,
        });
    }
    Ok(out)
}

/// Load a blocks TSV produced by `chain_and_merge` (via the `chain`
/// subcommand). Extra columns (cells/score/est_id) are preserved.
pub fn load_blocks(path: &std::path::Path) -> Result<Vec<Block>, Box<dyn Error>> {
    let f = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(f);
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() < 9 {
            continue;
        }
        out.push(Block {
            chr_a: t[0].to_string(),
            a0: t[1].parse()?,
            a1: t[2].parse()?,
            chr_b: t[3].to_string(),
            b0: t[4].parse()?,
            b1: t[5].parse()?,
            cells: t[6].parse()?,
            score: t[7].parse()?,
            est_id: t[8].parse()?,
        });
    }
    Ok(out)
}

/// Invert est_identity: the `shared` count that corresponds to a given
/// per-base identity for full-coverage 2 kb windows.
pub fn shared_for_identity(p: f64, sketch: usize, k: usize) -> u32 {
    let j = p.powi(k as i32) / (2.0 - p.powi(k as i32));
    let s = sketch as f64;
    let shared = (2.0 * s * j / (1.0 + j)).ceil();
    shared.max(1.0) as u32
}

/// Rough per-base identity estimate from mean shared sketch mass.
/// If J is the k-mer-set Jaccard and p the per-base identity, a k-mer
/// survives with probability p^k = J / (J + 1), so invert that.
pub fn est_identity(mean_shared: f64, sketch: usize, k: usize) -> f64 {
    let s = sketch as f64;
    if mean_shared <= 0.0 {
        return 0.0;
    }
    if mean_shared >= s {
        return 100.0;
    }
    let j = mean_shared / (2.0 * s - mean_shared);
    let p = (2.0 * j / (1.0 + j)).powf(1.0 / k as f64);
    100.0 * p
}

fn overlap_frac(x0: u32, x1: u32, y0: u32, y1: u32) -> f64 {
    let hi = x1.min(y1);
    let lo = x0.max(y0);
    if hi < lo {
        return 0.0;
    }
    let ov = (hi - lo + 1) as f64;
    ov / ((x1 - x0 + 1).min(y1 - y0 + 1)) as f64
}

/// Half-open u64 variant for genomic coordinates.
fn overlap_frac64(x0: u64, x1: u64, y0: u64, y1: u64) -> f64 {
    let hi = x1.min(y1);
    let lo = x0.max(y0);
    if hi <= lo {
        return 0.0;
    }
    let ov = (hi - lo) as f64;
    ov / ((x1 - x0).min(y1 - y0)) as f64
}

struct Dsu {
    parent: Vec<usize>,
}

impl Dsu {
    fn new(n: usize) -> Self {
        Dsu {
            parent: (0..n).collect(),
        }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        let mut c = x;
        while self.parent[c] != r {
            let n = self.parent[c];
            self.parent[c] = r;
            c = n;
        }
        r
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[ra] = rb;
        }
    }
}

pub struct ChainStats {
    pub n_pairs: usize,
    pub n_chains: usize,
    pub n_blocks: usize,
}

pub fn chain_and_merge(pairs: &[PairRec], p: &ChainParams) -> (Vec<Block>, ChainStats) {
    // ---- pass 1: band grouping + per-band diagonal walk ----------------
    // The band (d.div_euclid(2)) groups diagonals in pairs so a one-window
    // phase drift between copies stays inside one chain.
    let mut groups: HashMap<(String, String, i64), Vec<(u32, u32, u32)>> = HashMap::new();
    for r in pairs {
        let a = (r.s_a / p.window) as u32;
        let b = (r.s_b / p.window) as u32;
        let d = b as i64 - a as i64;
        groups
            .entry((r.chr_a.clone(), r.chr_b.clone(), d.div_euclid(2)))
            .or_default()
            .push((a, b, r.shared));
    }

    let mut chains: Vec<Chain> = Vec::new();
    for ((ca, cb, _band), mut cells) in groups {
        cells.sort_unstable();
        let mut open: Option<Chain> = None;
        let mut i = 0;
        while i < cells.len() {
            // The same `a` can appear on both diagonals of a band; keep the
            // strongest cell so the walk stays strictly increasing.
            let mut j = i;
            let mut best = cells[i];
            while j + 1 < cells.len() && cells[j + 1].0 == cells[i].0 {
                j += 1;
                if cells[j].2 > best.2 {
                    best = cells[j];
                }
            }
            i = j + 1;
            let (a, b, s) = best;
            let extend = matches!(&open, Some(c) if
                c.a1 < a && c.b1 < b
                && a - c.a1 <= p.gap && b - c.b1 <= p.gap);
            if extend {
                let c = open.as_mut().unwrap();
                c.a1 = a;
                c.b1 = b;
                c.score += s;
                c.cells += 1;
            } else {
                if let Some(c) = open.take() {
                    chains.push(c);
                }
                open = Some(Chain {
                    ca: ca.clone(),
                    cb: cb.clone(),
                    a0: a,
                    a1: a,
                    b0: b,
                    b1: b,
                    score: s,
                    cells: 1,
                });
            }
        }
        if let Some(c) = open.take() {
            chains.push(c);
        }
    }
    let n_pairs = pairs.len();
    chains.retain(|c| c.score >= p.min_score);
    let n_chains = chains.len();

    // ---- pass 2: merge chains overlapping on both sides ----------------
    let mut dsu = Dsu::new(chains.len());
    let mut order: Vec<usize> = (0..chains.len()).collect();
    order.sort_by_key(|&i| (chains[i].a0, chains[i].b0));
    let mut active: Vec<usize> = Vec::new();
    for &i in &order {
        active.retain(|&j| chains[j].a1 >= chains[i].a0);
        for &j in &active {
            if overlap_frac(chains[i].a0, chains[i].a1, chains[j].a0, chains[j].a1)
                >= p.merge_frac
                && overlap_frac(chains[i].b0, chains[i].b1, chains[j].b0, chains[j].b1)
                    >= p.merge_frac
            {
                dsu.union(i, j);
            }
        }
        active.push(i);
    }

    let mut agg: HashMap<usize, (String, String, u32, u32, u32, u32, u32, u32)> = HashMap::new();
    for (i, c) in chains.iter().enumerate() {
        let root = dsu.find(i);
        let e = agg
            .entry(root)
            .or_insert_with(|| (c.ca.clone(), c.cb.clone(), u32::MAX, 0, u32::MAX, 0, 0, 0));
        e.2 = e.2.min(c.a0);
        e.3 = e.3.max(c.a1);
        e.4 = e.4.min(c.b0);
        e.5 = e.5.max(c.b1);
        e.6 += c.score;
        e.7 += c.cells;
    }

    let w = p.window;
    let mut blocks: Vec<Block> = agg
        .into_iter()
        .map(|(_, (ca, cb, a0, a1, b0, b1, score, cells))| {
            let mean_shared = score as f64 / cells as f64;
            Block {
                chr_a: ca,
                a0: a0 as u64 * w,
                a1: (a1 as u64 + 1) * w,
                chr_b: cb,
                b0: b0 as u64 * w,
                b1: (b1 as u64 + 1) * w,
                cells,
                score,
                est_id: est_identity(mean_shared, p.sketch, p.k),
            }
        })
        .collect();
    blocks.sort_by(|x, y| y.score.cmp(&x.score));
    blocks.retain(|b| b.est_id >= p.min_block_id);

    // Chain-level merging can still leave near-duplicate blocks behind
    // (aggregate intervals overlap even though no individual chain pair
    // qualified). Blocks are few; merge once more at block level.
    let blocks = merge_blocks(blocks, p.merge_frac, p.sketch, p.k);
    let blocks = {
        let mut b = blocks;
        b.retain(|x| x.est_id >= p.min_block_id);
        b.sort_by(|x, y| y.score.cmp(&x.score));
        b
    };

    let n_blocks = blocks.len();
    (
        blocks,
        ChainStats {
            n_pairs,
            n_chains,
            n_blocks,
        },
    )
}

/// Union blocks that overlap on both sides by at least `merge_frac` of the
/// smaller span (second-level cleanup of grid-phase near-duplicates).
fn merge_blocks(blocks: Vec<Block>, merge_frac: f64, sketch: usize, k: usize) -> Vec<Block> {
    let n = blocks.len();
    if n <= 1 {
        return blocks;
    }
    let mut dsu = Dsu::new(n);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| (blocks[i].a0, blocks[i].b0));
    let mut active: Vec<usize> = Vec::new();
    for &i in &order {
        active.retain(|&j| blocks[j].a1 > blocks[i].a0);
        for &j in &active {
            if overlap_frac64(blocks[i].a0, blocks[i].a1, blocks[j].a0, blocks[j].a1)
                >= merge_frac
                && overlap_frac64(blocks[i].b0, blocks[i].b1, blocks[j].b0, blocks[j].b1)
                    >= merge_frac
            {
                dsu.union(i, j);
            }
        }
        active.push(i);
    }

    let mut agg: HashMap<usize, Block> = HashMap::new();
    for (i, b) in blocks.into_iter().enumerate() {
        let root = dsu.find(i);
        match agg.get_mut(&root) {
            None => {
                agg.insert(root, b);
            }
            Some(m) => {
                m.a0 = m.a0.min(b.a0);
                m.a1 = m.a1.max(b.a1);
                m.b0 = m.b0.min(b.b0);
                m.b1 = m.b1.max(b.b1);
                m.cells += b.cells;
                m.score += b.score;
            }
        }
    }
    let mut out: Vec<Block> = agg.into_values().collect();
    for b in &mut out {
        let mean_shared = b.score as f64 / b.cells as f64;
        b.est_id = est_identity(mean_shared, sketch, k);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> ChainParams {
        ChainParams {
            window: 2000,
            gap: 2,
            min_score: 10,
            merge_frac: 0.5,
            sketch: 256,
            k: 16,
            min_cell_id: 0.0,
            min_block_id: 0.0,
            min_shared_cell: 1,
        }
    }

    fn rec(s_a: u64, s_b: u64, shared: u32) -> PairRec {
        PairRec {
            chr_a: "c".into(),
            s_a,
            e_a: s_a + 2000,
            chr_b: "c".into(),
            s_b,
            e_b: s_b + 2000,
            shared,
        }
    }

    #[test]
    fn diagonal_run_forms_one_chain() {
        let pairs = vec![
            rec(0, 10 * 2000, 30),
            rec(2000, 11 * 2000, 30),
            rec(4000, 12 * 2000, 30),
        ];
        let (blocks, stats) = chain_and_merge(&pairs, &params());
        assert_eq!(stats.n_chains, 1);
        assert_eq!(blocks.len(), 1);
        assert_eq!((blocks[0].a0, blocks[0].a1), (0, 6000));
        assert_eq!((blocks[0].b0, blocks[0].b1), (20_000, 26_000));
    }

    #[test]
    fn phase_drift_stays_in_one_chain() {
        // copies drift by one window (delta diag = 1 -> same band)
        let pairs = vec![rec(0, 0, 30), rec(2000, 2000, 30), rec(4000, 6000, 30)];
        let (blocks, stats) = chain_and_merge(&pairs, &params());
        assert_eq!(stats.n_chains, 1);
        assert_eq!((blocks[0].a0, blocks[0].a1), (0, 6000));
        assert_eq!((blocks[0].b0, blocks[0].b1), (0, 8000));
    }

    #[test]
    fn distinct_nested_duplications_do_not_merge() {
        // chain 1: big duplication A[0..9k] ~ B[20..29k]
        // chain 2: small duplication A[2..6k] ~ B[80..83k] (different partner)
        let pairs = vec![
            rec(0, 20_000, 30),
            rec(2000, 22_000, 30),
            rec(4000, 24_000, 30),
            rec(6000, 26_000, 30),
            rec(8000, 28_000, 30),
            rec(2000, 80_000, 30),
            rec(4000, 82_000, 30),
        ];
        let (blocks, _) = chain_and_merge(&pairs, &params());
        assert_eq!(blocks.len(), 2);
    }

    #[test]
    fn phase_shifted_fragments_merge_into_one_block() {
        // Same SD in two fragments: first 2 cells on diag 0, next 2 cells
        // shifted by one window (diag 1 -> different band -> second chain).
        // Fragments overlap completely on both sides -> one block.
        let pairs = vec![
            rec(0, 0, 30),
            rec(2000, 2000, 30),
            rec(4000, 6000, 30),
            rec(6000, 8000, 30),
        ];
        let (blocks, _) = chain_and_merge(&pairs, &params());
        assert_eq!(blocks.len(), 1);
        assert_eq!((blocks[0].a0, blocks[0].a1), (0, 8000));
    }
}
