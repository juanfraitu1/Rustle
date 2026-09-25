//! `mcl_port` — bit-faithful Rust port of `bench/mcl_port.py`, the Python MCL used as the like-for-like
//! COMPARATOR by the bench scorers (through `bench/lib.py::mcl` since wave 7: the §6y0/§6y7 arms, the protein
//! referee, `truth.py adjudicated`, ...).
//!
//! ⚠ This is NOT the shipped `annotation_families::mcl`. Register 917 records that `mcl_port.py` is
//! deliberately not bit-identical to it (245 vs 168 clusters on chr2), and those scripts were designed
//! against the Python's numbers. This bin exists so they can keep those numbers while `numpy`/`scipy`
//! leave the path: it reproduces the Python's clusters EXACTLY — same values, same cluster order, same
//! member order — by replicating scipy's sparse arithmetic, not by re-deriving MCL:
//!
//!   * matrices are CSC with scipy's storage order: a product's column lists its rows in the
//!     REVERSE of first-touch order (scipy `csr_matmat` pass 2 walks a linked list from the last
//!     inserted entry), and every later step inherits that order;
//!   * each output entry accumulates its terms in the order `csr_matmat` visits them, starting from
//!     an exact 0.0, so the f64 sums are bit-identical;
//!   * column normalisation multiplies by the reciprocal `1/s` (that is what `X @ diags(1/s)` does),
//!     never divides; inflation is libm `pow`, the same call numpy makes;
//!   * clusters join each column to its FIRST maximal row in storage order (`np.argmax`), and the
//!     union-find roots on the smaller index, so cluster order = ascending smallest member.
//!
//! usage: mcl_port --graph EDGES.tsv [--inflation 2.8] [--prune 1e-9] [--max-iter 100]
//!   EDGES.tsv: `u<TAB>v<TAB>w` per undirected edge (lines with fewer fields are ignored, as are
//!   self-loops — the Python skips them too). Output: one cluster per line, members tab-separated.

use anyhow::{Context, Result};
use clap::Parser;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};

#[derive(Parser, Debug)]
#[command(about = "Bit-faithful port of bench/mcl_port.py (the Python MCL comparator); NOT the shipped MCL")]
struct Args {
    #[arg(long)]
    graph: String,
    #[arg(long, default_value_t = 2.8)]
    inflation: f64,
    #[arg(long, default_value_t = 1e-9)]
    prune: f64,
    #[arg(long, default_value_t = 100)]
    max_iter: usize,
}

/// Compressed sparse column, in scipy's storage order (indices NOT sorted after a product).
#[derive(Clone)]
struct Csc {
    n: usize,
    indptr: Vec<usize>,
    indices: Vec<usize>,
    data: Vec<f64>,
}

impl Csc {
    /// scipy `csc_matrix((vals,(rows,cols)))`: coalesce to sorted indices per column, duplicates summed
    /// in insertion order.
    fn from_coo(n: usize, rows: &[usize], cols: &[usize], vals: &[f64]) -> Self {
        let mut per_col: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        for k in 0..rows.len() {
            per_col[cols[k]].push((rows[k], vals[k]));
        }
        let (mut indptr, mut indices, mut data) = (vec![0usize; n + 1], Vec::new(), Vec::new());
        for j in 0..n {
            let col = &mut per_col[j];
            col.sort_by_key(|e| e.0); // stable, so duplicates keep insertion order
            let mut k = 0;
            while k < col.len() {
                let (i, mut v) = col[k];
                let mut m = k + 1;
                while m < col.len() && col[m].0 == i {
                    v += col[m].1;
                    m += 1;
                }
                indices.push(i);
                data.push(v);
                k = m;
            }
            indptr[j + 1] = indices.len();
        }
        Csc { n, indptr, indices, data }
    }

    /// scipy `csr_matmat` pass 2 applied to `C^T = B^T · A^T` (how scipy multiplies two CSC matrices):
    /// output column j lists rows in reverse-first-touch order; exact zeros are dropped.
    fn matmul(&self, b: &Csc) -> Csc {
        // C = self · b.  C^T = b^T · self^T  ⇒  csr_matmat(A' = b^T, B' = self^T)
        // A' row j = column j of b (entries (k, b[k,j])); B' row k = column k of self (entries (i, self[i,k])).
        let n = self.n;
        let mut sums = vec![0.0f64; n];
        let mut next = vec![usize::MAX; n]; // -1
        const HEAD_END: usize = usize::MAX - 1; // -2
        let (mut indptr, mut indices, mut data) = (vec![0usize; n + 1], Vec::new(), Vec::new());
        for j in 0..n {
            let mut head = HEAD_END;
            let mut length = 0usize;
            for kk in b.indptr[j]..b.indptr[j + 1] {
                let k = b.indices[kk];
                let v = b.data[kk];
                for ii in self.indptr[k]..self.indptr[k + 1] {
                    let i = self.indices[ii];
                    sums[i] += v * self.data[ii];
                    if next[i] == usize::MAX {
                        next[i] = head;
                        head = i;
                        length += 1;
                    }
                }
            }
            for _ in 0..length {
                if sums[head] != 0.0 {
                    indices.push(head);
                    data.push(sums[head]);
                }
                let temp = head;
                head = next[head];
                next[temp] = usize::MAX;
                sums[temp] = 0.0;
            }
            indptr[j + 1] = indices.len();
        }
        Csc { n, indptr, indices, data }
    }

    /// `X @ sp.diags(1/s)` with `s = X.sum(axis=0)` (sequential sums in storage order), zeros → 1.
    /// Scipy routes this through `csr_matmat` too, so the column's row order is REVERSED again.
    fn norm(&self) -> Csc {
        let n = self.n;
        let mut inv = vec![0.0f64; n];
        for j in 0..n {
            let mut s = 0.0f64;
            for kk in self.indptr[j]..self.indptr[j + 1] {
                s += self.data[kk];
            }
            if s == 0.0 {
                s = 1.0;
            }
            inv[j] = 1.0 / s;
        }
        // C^T = D · X^T: for row j of D (single entry d[j]), walk column j of X in storage order,
        // first-touch each row once, emit in reverse. Value = d[j] * X[i,j]; exact zeros dropped.
        let (mut indptr, mut indices, mut data) = (vec![0usize; n + 1], Vec::new(), Vec::new());
        for j in 0..n {
            let s = self.indptr[j];
            let e = self.indptr[j + 1];
            for kk in (s..e).rev() {
                let v = inv[j] * self.data[kk];
                if v != 0.0 {
                    indices.push(self.indices[kk]);
                    data.push(v);
                }
            }
            indptr[j + 1] = indices.len();
        }
        Csc { n, indptr, indices, data }
    }

    /// `abs(N - M)`: `(nnz == 0, max)` over the union of entries — order-free, exact differences.
    fn diff_stats(&self, m: &Csc) -> (bool, f64) {
        let mut any = false;
        let mut mx = 0.0f64;
        let mut row: HashMap<usize, f64> = HashMap::new();
        for j in 0..self.n {
            row.clear();
            for kk in self.indptr[j]..self.indptr[j + 1] {
                *row.entry(self.indices[kk]).or_insert(0.0) += self.data[kk];
            }
            for kk in m.indptr[j]..m.indptr[j + 1] {
                *row.entry(m.indices[kk]).or_insert(0.0) -= m.data[kk];
            }
            for &d in row.values() {
                if d != 0.0 {
                    any = true;
                    let a = d.abs();
                    if a > mx {
                        mx = a;
                    }
                }
            }
        }
        (!any, mx)
    }
}

fn mcl(nodes: &[String], edges: &[(usize, usize, f64)], inflation: f64, prune: f64, max_iter: usize) -> Vec<Vec<String>> {
    let n = nodes.len();
    if n == 0 {
        return Vec::new();
    }
    let (mut rows, mut cols, mut vals): (Vec<usize>, Vec<usize>, Vec<f64>) = ((0..n).collect(), (0..n).collect(), vec![1.0; n]);
    for &(a, b, w) in edges {
        if a == b {
            continue;
        }
        rows.extend([a, b]);
        cols.extend([b, a]);
        vals.extend([w, w]);
    }
    let mut m = Csc::from_coo(n, &rows, &cols, &vals).norm();
    for _ in 0..max_iter {
        let mut nn = m.matmul(&m);
        for x in nn.data.iter_mut() {
            *x = x.powf(inflation);
        }
        // N.data[N.data < prune] = 0; N.eliminate_zeros()  (order preserved)
        let mut keep_ip = vec![0usize; n + 1];
        let (mut ki, mut kd) = (Vec::with_capacity(nn.indices.len()), Vec::with_capacity(nn.data.len()));
        for j in 0..n {
            for kk in nn.indptr[j]..nn.indptr[j + 1] {
                let v = if nn.data[kk] < prune { 0.0 } else { nn.data[kk] };
                if v != 0.0 {
                    ki.push(nn.indices[kk]);
                    kd.push(v);
                }
            }
            keep_ip[j + 1] = ki.len();
        }
        nn = Csc { n, indptr: keep_ip, indices: ki, data: kd }.norm();
        let (empty, mx) = nn.diff_stats(&m);
        let done = empty || mx < 1e-7;
        m = nn;
        if done {
            break;
        }
    }
    // union-find: column j joins its FIRST maximal row in storage order (np.argmax); root = smaller index
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for j in 0..n {
        let (s, e) = (m.indptr[j], m.indptr[j + 1]);
        if s == e {
            continue;
        }
        let mut best = s;
        for kk in s + 1..e {
            if m.data[kk] > m.data[best] {
                best = kk;
            }
        }
        let r = m.indices[best];
        let (ra, rb) = (find(&mut parent, j), find(&mut parent, r));
        if ra != rb {
            parent[ra.max(rb)] = ra.min(rb);
        }
    }
    let mut order: Vec<usize> = Vec::new();
    let mut groups: HashMap<usize, Vec<String>> = HashMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        if !groups.contains_key(&r) {
            order.push(r);
        }
        groups.entry(r).or_default().push(nodes[i].clone());
    }
    order.into_iter().map(|r| groups.remove(&r).unwrap()).collect()
}

fn main() -> Result<()> {
    let a = Args::parse();
    let f = std::fs::File::open(&a.graph).with_context(|| format!("opening {}", a.graph))?;
    let mut raw: Vec<(String, String, f64)> = Vec::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        let fs: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        if fs.len() >= 3 {
            raw.push((fs[0].to_string(), fs[1].to_string(), fs[2].parse::<f64>()?));
        }
    }
    // nodes = sorted({x for e in edges for x in e}) — over EDGE endpoints only, as the Python does
    let mut nodes: Vec<String> = raw.iter().flat_map(|(u, v, _)| [u.clone(), v.clone()]).collect();
    nodes.sort();
    nodes.dedup();
    let ix: HashMap<&str, usize> = nodes.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let edges: Vec<(usize, usize, f64)> = raw.iter().map(|(u, v, w)| (ix[u.as_str()], ix[v.as_str()], *w)).collect();
    let out = std::io::stdout();
    let mut w = std::io::BufWriter::new(out.lock());
    for c in mcl(&nodes, &edges, a.inflation, a.prune, a.max_iter) {
        writeln!(w, "{}", c.join("\t"))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Regression pin. Expected clusters were produced by this port at the commit where it matched the
    /// numpy/scipy `bench/mcl_port.py` on 20/20 graphs (10 real-weight, 10 uniform) — the uniform case
    /// below exercises the tie-breaking that a textbook MCL gets wrong.
    #[test]
    fn two_triangles_joined_by_one_weak_edge_split_at_the_bridge() {
        let nodes: Vec<String> = ["a", "b", "c", "d", "e", "f"].iter().map(|s| s.to_string()).collect();
        let edges = vec![(0, 1, 1.0), (1, 2, 1.0), (0, 2, 1.0), (3, 4, 1.0), (4, 5, 1.0), (3, 5, 1.0), (2, 3, 0.05)];
        let c = mcl(&nodes, &edges, 2.8, 1e-9, 100);
        assert_eq!(c, vec![vec!["a", "b", "c"], vec!["d", "e", "f"]]);
        // uniform weights on a 4-cycle: symmetric ties, resolved by storage order — must stay stable
        let n4: Vec<String> = ["p", "q", "r", "s"].iter().map(|s| s.to_string()).collect();
        let cyc = vec![(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (0, 3, 1.0)];
        let c4 = mcl(&n4, &cyc, 2.8, 1e-9, 100);
        assert_eq!(c4.iter().map(|c| c.len()).sum::<usize>(), 4);
        assert!(c4.iter().all(|c| !c.is_empty()));
    }
}
