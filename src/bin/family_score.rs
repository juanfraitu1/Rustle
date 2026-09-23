//! Score one node-construction MODE's clusters against a family truth (Soto, the protein referee, or a
//! union truth), reporting sensitivity / precision / one-to-one bipartite F and the collapse count.
//!
//! Native Rust replacement for `bench/mode_family_score.py` (§6x0-§6z1), reproducing it byte for byte
//! on every (clusters, gff, truth, chrom, family) combination run in the 2026-09-22 session. Payoff, per
//! register 897, is the dependency removed: this drops `numpy` + `scipy.optimize.linear_sum_assignment`
//! from the family-scoring path — the scorer the standing reporting rule names (sens / prec / bipartite).
//!
//! ⚠ The bipartite matching here is EVALUATION ONLY. Families are never built with it (standing project
//! constraint); this only scores clusters that were already built without it.
//!
//! Semantics copied from the Python, including its resolvers' biases (deliberately, so arms stay
//! comparable — register 964/§6v9 showed multi-labelling is worse and §6v8 that winner-take-all
//! undercounts recall, but the bias is identical across arms):
//!   * locus -> ONE gene by max overlap, first maximum wins;
//!   * a predicted cluster is intersected with the truth universe (⚠ register 770: this DELETES unlabelled
//!     members from numerator and denominator, so precision is an upper bound — §6x1/r991);
//!   * collapse is computed gene -> best-covering locus (locus -> gene is 0 by construction, §6x0).
//!
//! usage: family_score --clusters X.clusters.tsv --gff chrN.genes.gff --soto TRUTH.tsv [--chrom chr16]
//!        [--family NPIP] [--label arm]

use anyhow::{Context, Result};
use clap::Parser;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, BufReader};

#[derive(Parser, Debug)]
#[command(about = "Score mode clusters against a family truth (sens / prec / bipartite F / collapse)")]
struct Args {
    /// `mcl_families` clusters TSV (needs columns cluster_id, chrom, start, end)
    #[arg(long)]
    clusters: String,
    /// GFF3 with gene / pseudogene / ncRNA_gene records carrying `Name=`
    #[arg(long)]
    gff: String,
    /// Truth TSV with header columns `Gene Name` and `Family ID` (Soto S1C format; a gene may repeat)
    #[arg(long)]
    soto: String,
    #[arg(long, default_value = "chr16")]
    chrom: String,
    /// Restrict truth to families containing a gene whose name contains this substring
    #[arg(long)]
    family: Option<String>,
    #[arg(long, default_value = "arm")]
    label: String,
}

/// `Gene Name` -> first `Family ID` seen (Python `setdefault`), skipping blanks and `N/A`; insertion-ordered.
fn soto_truth(path: &str) -> Result<Vec<(String, String)>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {path}"))?;
    let mut lines = BufReader::new(f).lines();
    let hdr: Vec<String> = match lines.next() {
        Some(h) => h?.split('\t').map(|s| s.to_string()).collect(),
        None => return Ok(Vec::new()),
    };
    let col = |name: &str| hdr.iter().position(|h| h == name);
    let (ci_f, ci_g) = (col("Family ID"), col("Gene Name"));
    let mut seen: HashSet<String> = HashSet::new();
    let mut out = Vec::new();
    for line in lines {
        let line = line?;
        let r: Vec<&str> = line.split('\t').collect();
        let get = |i: Option<usize>| i.and_then(|i| r.get(i)).map(|s| s.trim()).unwrap_or("");
        let (fam, g) = (get(ci_f), get(ci_g));
        if !fam.is_empty() && fam != "N/A" && !g.is_empty() && seen.insert(g.to_string()) {
            out.push((g.to_string(), fam.to_string()));
        }
    }
    Ok(out)
}

fn name_attr(attrs: &str) -> Option<&str> {
    let i = attrs.find("Name=")? + 5;
    let rest = &attrs[i..];
    Some(rest.split(';').next().unwrap_or(rest))
}

/// (start, end, name) for gene-like records on `chrom`, sorted as Python sorts the tuples.
fn gene_spans(gff: &str, chrom: &str) -> Result<Vec<(i64, i64, String)>> {
    let f = std::fs::File::open(gff).with_context(|| format!("opening {gff}"))?;
    let mut out = Vec::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fs: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        if fs.len() < 9 || fs[0] != chrom || !matches!(fs[2], "gene" | "pseudogene" | "ncRNA_gene") {
            continue;
        }
        if let Some(n) = name_attr(fs[8]) {
            out.push((fs[3].parse::<i64>()?, fs[4].parse::<i64>()?, n.to_string()));
        }
    }
    out.sort();
    Ok(out)
}

/// cluster_id -> members on `chrom`, in first-seen cluster order and file order within a cluster.
fn load_clusters(path: &str, chrom: &str) -> Result<Vec<(String, Vec<(i64, i64)>)>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    let mut rows = text.split('\n').map(|l| l.trim_end_matches('\n'));
    let hdr: Vec<&str> = rows.next().unwrap_or("").split('\t').collect();
    let ci: HashMap<&str, usize> = hdr.iter().enumerate().map(|(i, h)| (*h, i)).collect();
    for need in ["cluster_id", "chrom", "start", "end"] {
        if !ci.contains_key(need) {
            anyhow::bail!("{path}: expected columns cluster_id/chrom/start/end, got {:?}", &hdr[..hdr.len().min(8)]);
        }
    }
    let mut order: Vec<String> = Vec::new();
    let mut members: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    for line in rows {
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < hdr.len() || r[ci["chrom"]] != chrom {
            continue;
        }
        let cid = r[ci["cluster_id"]].to_string();
        let (s, e) = (r[ci["start"]].parse::<i64>()?, r[ci["end"]].parse::<i64>()?);
        if !members.contains_key(&cid) {
            order.push(cid.clone());
        }
        members.entry(cid).or_default().push((s, e));
    }
    Ok(order.into_iter().map(|c| { let m = members.remove(&c).unwrap_or_default(); (c, m) }).collect())
}

/// Best gene for a locus by max overlap, FIRST maximum wins (strict `>`), scanning the sorted spans.
fn gene_at(spans: &[(i64, i64, String)], s: i64, e: i64) -> Option<&str> {
    let mut best: Option<(i64, &str)> = None;
    for (gs, ge, g) in spans {
        if *ge < s {
            continue;
        }
        if *gs > e {
            break;
        }
        let ov = e.min(*ge) - s.max(*gs);
        if ov > 0 && best.map_or(true, |b| ov > b.0) {
            best = Some((ov, g));
        }
    }
    best.map(|b| b.1)
}

/// scipy's `linear_sum_assignment` on a rectangular cost matrix — a faithful port of
/// `scipy/optimize/rectangular_lsap.cpp` (Crouse 2016, shortest augmenting path with duals), INCLUDING its
/// tie-breaking: columns are scanned from a `remaining` list filled in REVERSE index order, an equal-cost
/// column is preferred when still unassigned, and a chosen column is swap-removed from the tail. Ties are
/// the whole reason this is a port rather than a textbook Hungarian: equally-optimal assignments differ in
/// WHICH clusters they match, and `prec` (Σ matched cluster sizes) is not tie-invariant even though `sens`
/// is. Costs are exact integers, so every comparison is bit-for-bit what scipy computes on doubles.
/// Requires `nr <= nc`; the caller transposes otherwise, as scipy does.
fn lsap_min_cost(nr: usize, nc: usize, cost: &[i64]) -> Vec<usize> {
    const INF: i64 = i64::MAX / 4;
    let (mut u, mut v) = (vec![0i64; nr], vec![0i64; nc]);
    let mut shortest = vec![INF; nc];
    let mut path = vec![usize::MAX; nc];
    let mut col4row = vec![usize::MAX; nr];
    let mut row4col = vec![usize::MAX; nc];
    let (mut sr, mut sc) = (vec![false; nr], vec![false; nc]);
    let mut remaining = vec![0usize; nc];
    for cur_row in 0..nr {
        // augmenting path from cur_row
        let mut min_val = 0i64;
        let mut num_remaining = nc;
        for it in 0..nc {
            remaining[it] = nc - it - 1;
        }
        sr.iter_mut().for_each(|x| *x = false);
        sc.iter_mut().for_each(|x| *x = false);
        shortest.iter_mut().for_each(|x| *x = INF);
        let mut sink = usize::MAX;
        let mut i = cur_row;
        while sink == usize::MAX {
            let mut index = usize::MAX;
            let mut lowest = INF;
            sr[i] = true;
            for it in 0..num_remaining {
                let j = remaining[it];
                let r = min_val + cost[i * nc + j] - u[i] - v[j];
                if r < shortest[j] {
                    path[j] = i;
                    shortest[j] = r;
                }
                if shortest[j] < lowest || (shortest[j] == lowest && row4col[j] == usize::MAX) {
                    lowest = shortest[j];
                    index = it;
                }
            }
            min_val = lowest;
            debug_assert!(min_val < INF, "infeasible assignment");
            let j = remaining[index];
            if row4col[j] == usize::MAX {
                sink = j;
            } else {
                i = row4col[j];
            }
            sc[j] = true;
            num_remaining -= 1;
            remaining[index] = remaining[num_remaining];
        }
        // update duals
        u[cur_row] += min_val;
        for k in 0..nr {
            if sr[k] && k != cur_row {
                u[k] += min_val - shortest[col4row[k]];
            }
        }
        for j in 0..nc {
            if sc[j] {
                v[j] -= min_val - shortest[j];
            }
        }
        // augment along the path
        let mut j = sink;
        loop {
            let i2 = path[j];
            row4col[j] = i2;
            std::mem::swap(&mut col4row[i2], &mut j);
            if i2 == cur_row {
                break;
            }
        }
    }
    col4row
}

/// scipy `linear_sum_assignment(-M)` on a rectangular overlap matrix: maximise total overlap. Returns
/// (row, col) pairs; every row is assigned when rows <= cols, every column otherwise.
fn max_overlap_assignment(m: &[Vec<i64>]) -> Vec<(usize, usize)> {
    let n = m.len();
    let k = if n == 0 { 0 } else { m[0].len() };
    if n == 0 || k == 0 {
        return Vec::new();
    }
    if n <= k {
        let cost: Vec<i64> = m.iter().flat_map(|r| r.iter().map(|x| -x)).collect();
        lsap_min_cost(n, k, &cost).into_iter().enumerate().collect()
    } else {
        // transpose, solve, map back — scipy does exactly this for nr > nc
        let cost: Vec<i64> = (0..k).flat_map(|j| (0..n).map(move |i| (i, j))).map(|(i, j)| -m[i][j]).collect();
        lsap_min_cost(k, n, &cost).into_iter().enumerate().map(|(j, i)| (i, j)).collect()
    }
}

fn main() -> Result<()> {
    let a = Args::parse();
    let fam = soto_truth(&a.soto)?;
    let spans = gene_spans(&a.gff, &a.chrom)?;
    let clusters = load_clusters(&a.clusters, &a.chrom)?;

    // locus -> ONE gene, by max overlap
    let mut locus_gene: HashMap<(String, i64, i64), String> = HashMap::new();
    for (cid, members) in &clusters {
        for &(s, e) in members {
            if let Some(g) = gene_at(&spans, s, e) {
                locus_gene.insert((cid.clone(), s, e), g.to_string());
            }
        }
    }

    // truth families on this chromosome, in first-seen family order
    let on_chrom: HashSet<&str> = spans.iter().map(|x| x.2.as_str()).collect();
    let mut truth_order: Vec<String> = Vec::new();
    let mut truth: HashMap<String, HashSet<String>> = HashMap::new();
    for (g, f) in &fam {
        if on_chrom.contains(g.as_str()) {
            if !truth.contains_key(f) {
                truth_order.push(f.clone());
            }
            truth.entry(f.clone()).or_default().insert(g.clone());
        }
    }
    let t_order: Vec<String> = truth_order
        .into_iter()
        .filter(|f| a.family.as_ref().map_or(true, |sub| truth[f].iter().any(|g| g.contains(sub.as_str()))))
        .filter(|f| truth[f].len() >= 2)
        .collect();
    let universe: HashSet<&str> = t_order.iter().flat_map(|f| truth[f].iter().map(|g| g.as_str())).collect();

    // predicted clusters as gene sets, intersected with the truth universe
    let mut p_order: Vec<String> = Vec::new();
    let mut pred: HashMap<String, HashSet<String>> = HashMap::new();
    for (cid, members) in &clusters {
        let gs: HashSet<String> = members
            .iter()
            .filter_map(|&(s, e)| locus_gene.get(&(cid.clone(), s, e)))
            .filter(|g| universe.contains(g.as_str()))
            .cloned()
            .collect();
        if !gs.is_empty() {
            p_order.push(cid.clone());
            pred.insert(cid.clone(), gs);
        }
    }
    if t_order.is_empty() || p_order.is_empty() {
        println!("{}: no scoreable truth/prediction overlap", a.label);
        return Ok(());
    }

    let m: Vec<Vec<i64>> = t_order
        .iter()
        .map(|tf| p_order.iter().map(|pc| truth[tf].intersection(&pred[pc]).count() as i64).collect())
        .collect();
    let assignment = max_overlap_assignment(&m);
    let matched: i64 = assignment.iter().map(|&(i, j)| m[i][j]).sum();
    let tot_truth: usize = t_order.iter().map(|f| truth[f].len()).sum();
    let tot_pred: usize = assignment.iter().filter(|&&(i, j)| m[i][j] > 0).map(|&(_, j)| pred[&p_order[j]].len()).sum();
    let sens = if tot_truth > 0 { matched as f64 / tot_truth as f64 } else { 0.0 };
    let prec = if tot_pred > 0 { matched as f64 / tot_pred as f64 } else { 0.0 };
    let f1 = if sens + prec > 0.0 { 2.0 * sens * prec / (sens + prec) } else { 0.0 };

    // collapse (register 817's failure mode): each TRUTH GENE -> the locus that best covers it; genes that
    // must share one locus. ⚠ gene -> locus, never locus -> gene (0 by construction).
    let truth_genes: HashSet<&str> = universe.clone();
    let all_loci: Vec<(&str, i64, i64)> =
        clusters.iter().flat_map(|(c, ms)| ms.iter().map(move |&(s, e)| (c.as_str(), s, e))).collect();
    let mut gene_span: HashMap<&str, (i64, i64)> = HashMap::new();
    for (gs, ge, g) in &spans {
        gene_span.insert(g.as_str(), (*gs, *ge)); // last wins on duplicate names, as the Python dict does
    }
    let mut share: BTreeMap<(&str, i64, i64), usize> = BTreeMap::new();
    let mut placed = 0usize;
    for g in &truth_genes {
        let Some(&(gs, ge)) = gene_span.get(g) else { continue };
        let mut best: Option<(i64, (&str, i64, i64))> = None;
        for &(cid, s, e) in &all_loci {
            let ov = ge.min(e) - gs.max(s);
            if ov > 0 && best.map_or(true, |b| ov > b.0) {
                best = Some((ov, (cid, s, e)));
            }
        }
        if let Some((_, l)) = best {
            *share.entry(l).or_insert(0) += 1;
            placed += 1;
        }
    }
    let collapsed: usize = share.values().filter(|&&n| n > 1).map(|&n| n - 1).sum();
    let missing = truth_genes.len() - placed;

    println!(
        "{:>14} | truth {} fams / {} genes | clusters {} | sens {:.3} prec {:.3} F {:.3} | collapsed {} | no-locus {}",
        a.label, t_order.len(), tot_truth, p_order.len(), sens, prec, f1, collapsed, missing
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hungarian_picks_the_max_overlap_assignment_on_a_rectangle() {
        // 2 truth families x 3 clusters; best is T0->C1 (3) + T1->C2 (2) = 5, not the greedy T0->C0.
        let m = vec![vec![2, 3, 0], vec![0, 1, 2]];
        let mut a = max_overlap_assignment(&m);
        a.sort();
        assert_eq!(a, vec![(0, 1), (1, 2)]);
        let total: i64 = a.iter().map(|&(i, j)| m[i][j]).sum();
        assert_eq!(total, 5);
        // more rows than columns: every column is assigned, rows may go unmatched
        let t = vec![vec![1], vec![5], vec![2]];
        assert_eq!(max_overlap_assignment(&t), vec![(1, 0)]);
        // scipy tie-break, checked against `linear_sum_assignment(-M)` on this matrix: a 2-way tie in
        // total overlap resolves to (0,0),(1,1), not (0,1),(1,0)
        let tie = vec![vec![1, 1], vec![1, 1]];
        assert_eq!(max_overlap_assignment(&tie), vec![(0, 0), (1, 1)]);
    }

    #[test]
    fn gene_at_is_first_max_overlap() {
        let spans = vec![(100, 200, "A".to_string()), (150, 400, "B".to_string()), (900, 950, "C".to_string())];
        assert_eq!(gene_at(&spans, 120, 260), Some("B")); // B overlaps 110, A overlaps 80
        assert_eq!(gene_at(&spans, 100, 200), Some("A")); // tie 100 vs 50 -> A
        assert_eq!(gene_at(&spans, 500, 800), None);
    }
}
