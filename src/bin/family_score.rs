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
//! GENOME-WIDE MODE (`--chrom ALL`, 2026-09-25). No chromosome restriction: every contig of the GFF and of the
//! clusters is read, and a gene is the pair (contig, Name), because RefSeq names repeat across contigs (CHM13: 66
//! names occur more than once, 49 of them on more than one chromosome, mostly the X/Y pseudoautosomal genes). A
//! locus is labelled with a gene on ITS OWN contig, and the collapse count looks for a gene's best-covering locus
//! on the gene's own contig. Truth: when the truth TSV has a contig column (`Chrom`, `chrom`, `Contig`, `contig`
//! or `seqid`), each row names one (contig, Name) gene, first family per (contig, Name); otherwise a row's Name
//! stands for every (contig, Name) gene the GFF has, first family per Name, as in the per-chromosome mode.
//! Families and pairs may therefore cross contigs. Restricted to one contig, this is the per-chromosome mode:
//! on inputs holding a single contig, `--chrom ALL` and `--chrom <that contig>` print the same numbers (a test
//! pins it). The per-chromosome mode itself is unchanged, byte for byte.
//!
//! OPT-IN OUTPUTS (the default output is unchanged): `--per-family OUT.tsv` writes one row per truth family (its
//! matched cluster, sizes, recall / precision / F / Jaccard, its truth and recovered pairs, its members), with
//! the same definitions as `figures/_o1_recovery.py::score_arm`; `--pairwise` prints a second line with the
//! pairwise counts over the truth universe (truth pairs, predicted pairs, true-positive pairs, and the pairwise
//! sensitivity / precision / F, which unlike the bipartite precision do not depend on how ties are broken).
//!
//! usage: family_score --clusters X.clusters.tsv --gff genes.gff --soto TRUTH.tsv [--chrom chr16 | --chrom ALL]
//!        [--family NPIP] [--label arm] [--per-family OUT.tsv] [--pairwise]

use anyhow::{Context, Result};
use clap::Parser;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};

/// `--chrom` value selecting the genome-wide mode.
const ALL: &str = "ALL";

#[derive(Parser, Debug)]
#[command(about = "Score mode clusters against a family truth (sens / prec / bipartite F / collapse)")]
struct Args {
    /// `mcl_families` clusters TSV (needs columns cluster_id, chrom, start, end)
    #[arg(long)]
    clusters: String,
    /// GFF3 with gene / pseudogene / ncRNA_gene records carrying `Name=`
    #[arg(long)]
    gff: String,
    /// Truth TSV with header columns `Gene Name` and `Family ID` (Soto S1C format; a gene may repeat). In
    /// `--chrom ALL` mode an optional contig column (`Chrom`/`chrom`/`Contig`/`contig`/`seqid`) keys genes by
    /// (contig, Name).
    #[arg(long)]
    soto: String,
    /// Chromosome to score, or `ALL` for the genome-wide mode (genes keyed by (contig, Name); see the module doc)
    #[arg(long, default_value = "chr16")]
    chrom: String,
    /// Restrict truth to families containing a gene whose name contains this substring
    #[arg(long)]
    family: Option<String>,
    #[arg(long, default_value = "arm")]
    label: String,
    /// Opt-in: write one row per truth family (matched cluster, members, recall, precision, F, pairs) to this TSV
    #[arg(long)]
    per_family: Option<String>,
    /// Opt-in: print a second line with the pairwise counts (truth / predicted / true-positive pairs over the
    /// truth universe) and the pairwise sensitivity / precision / F
    #[arg(long, default_value_t = false)]
    pairwise: bool,
}

/// A gene: (contig, Name). In the per-chromosome mode every gene has the same contig, so this is the Name.
type Gene = (String, String);

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

/// Contig column names a genome-wide truth may carry, in order of preference.
const TRUTH_CONTIG_COLUMNS: [&str; 5] = ["Chrom", "chrom", "Contig", "contig", "seqid"];

/// Genome-wide truth rows `(contig, Name, family)`. With a contig column: one (contig, Name) gene per row, first
/// family per (contig, Name). Without one: exactly [`soto_truth`], contig `None` (the Name stands for every
/// contig's gene of that Name).
fn truth_rows_gw(path: &str) -> Result<Vec<(Option<String>, String, String)>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {path}"))?;
    let mut lines = BufReader::new(f).lines();
    let hdr: Vec<String> = match lines.next() {
        Some(h) => h?.split('\t').map(|s| s.to_string()).collect(),
        None => return Ok(Vec::new()),
    };
    let col = |name: &str| hdr.iter().position(|h| h == name);
    let Some(ci_c) = TRUTH_CONTIG_COLUMNS.iter().find_map(|c| col(c)) else {
        return Ok(soto_truth(path)?.into_iter().map(|(g, f)| (None, g, f)).collect());
    };
    let (ci_f, ci_g) = (col("Family ID"), col("Gene Name"));
    let mut seen: HashSet<(String, String)> = HashSet::new();
    let mut out = Vec::new();
    for line in lines {
        let line = line?;
        let r: Vec<&str> = line.split('\t').collect();
        let get = |i: Option<usize>| i.and_then(|i| r.get(i)).map(|s| s.trim()).unwrap_or("");
        let (fam, g, c) = (get(ci_f), get(ci_g), get(Some(ci_c)));
        if !fam.is_empty() && fam != "N/A" && !g.is_empty() && !c.is_empty() && seen.insert((c.to_string(), g.to_string())) {
            out.push((Some(c.to_string()), g.to_string(), fam.to_string()));
        }
    }
    Ok(out)
}

fn name_attr(attrs: &str) -> Option<&str> {
    let i = attrs.find("Name=")? + 5;
    let rest = &attrs[i..];
    Some(rest.split(';').next().unwrap_or(rest))
}

/// contig -> (start, end, name) of its gene-like records, each list sorted as Python sorts the tuples. `chrom` =
/// `Some(c)` keeps contig `c` only (the per-chromosome mode); `None` keeps every contig.
fn gene_spans(gff: &str, chrom: Option<&str>) -> Result<BTreeMap<String, Vec<(i64, i64, String)>>> {
    let f = std::fs::File::open(gff).with_context(|| format!("opening {gff}"))?;
    let mut out: BTreeMap<String, Vec<(i64, i64, String)>> = BTreeMap::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fs: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        if fs.len() < 9 || chrom.is_some_and(|c| fs[0] != c) || !matches!(fs[2], "gene" | "pseudogene" | "ncRNA_gene") {
            continue;
        }
        if let Some(n) = name_attr(fs[8]) {
            out.entry(fs[0].to_string()).or_default().push((fs[3].parse::<i64>()?, fs[4].parse::<i64>()?, n.to_string()));
        }
    }
    for v in out.values_mut() {
        v.sort();
    }
    Ok(out)
}

/// A cluster member locus: (contig, start, end).
type Locus = (String, i64, i64);

/// cluster_id -> member loci (on `chrom` when given, else on every contig), in first-seen cluster order and file
/// order within a cluster.
fn load_clusters(path: &str, chrom: Option<&str>) -> Result<Vec<(String, Vec<Locus>)>> {
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
    let mut members: HashMap<String, Vec<Locus>> = HashMap::new();
    for line in rows {
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < hdr.len() || chrom.is_some_and(|c| r[ci["chrom"]] != c) {
            continue;
        }
        let cid = r[ci["cluster_id"]].to_string();
        let (s, e) = (r[ci["start"]].parse::<i64>()?, r[ci["end"]].parse::<i64>()?);
        if !members.contains_key(&cid) {
            order.push(cid.clone());
        }
        members.entry(cid).or_default().push((r[ci["chrom"]].to_string(), s, e));
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

/// Everything one scoring run computes; `main` prints it and the opt-in outputs read it.
struct Score {
    /// truth families in first-seen order (>= 2 genes, `--family` applied)
    t_order: Vec<String>,
    truth: HashMap<String, HashSet<Gene>>,
    /// scored clusters (>= 1 universe gene) in first-seen order
    p_order: Vec<String>,
    pred: HashMap<String, HashSet<Gene>>,
    /// (truth row, cluster column) pairs of the one-to-one assignment
    assignment: Vec<(usize, usize)>,
    /// overlap matrix, truth rows x cluster columns
    m: Vec<Vec<i64>>,
    collapsed: usize,
    missing: usize,
}

fn score(a: &Args) -> Result<Score> {
    let genome_wide = a.chrom == ALL;
    let only = if genome_wide { None } else { Some(a.chrom.as_str()) };
    // (contig, Name, family): the per-chromosome mode reads the truth exactly as it always has
    let rows: Vec<(Option<String>, String, String)> = if genome_wide {
        truth_rows_gw(&a.soto)?
    } else {
        soto_truth(&a.soto)?.into_iter().map(|(g, f)| (None, g, f)).collect()
    };
    let spans = gene_spans(&a.gff, only)?;
    let clusters = load_clusters(&a.clusters, only)?;
    let no_spans: Vec<(i64, i64, String)> = Vec::new();
    let spans_of = |c: &str| spans.get(c).unwrap_or(&no_spans);

    // locus -> ONE gene on its own contig, by max overlap
    let mut locus_gene: HashMap<(String, Locus), Gene> = HashMap::new();
    for (cid, members) in &clusters {
        for (c, s, e) in members {
            if let Some(g) = gene_at(spans_of(c), *s, *e) {
                locus_gene.insert((cid.clone(), (c.clone(), *s, *e)), (c.clone(), g.to_string()));
            }
        }
    }

    // truth families over the genes the GFF has, in first-seen family order
    let mut contigs_of: HashMap<&str, BTreeSet<&str>> = HashMap::new();
    for (c, v) in &spans {
        for x in v {
            contigs_of.entry(x.2.as_str()).or_default().insert(c.as_str());
        }
    }
    let mut truth_order: Vec<String> = Vec::new();
    let mut truth: HashMap<String, HashSet<Gene>> = HashMap::new();
    for (c, g, f) in &rows {
        let keys: Vec<Gene> = match c {
            Some(c) => contigs_of
                .get(g.as_str())
                .filter(|cs| cs.contains(c.as_str()))
                .map(|_| vec![(c.clone(), g.clone())])
                .unwrap_or_default(),
            None => contigs_of
                .get(g.as_str())
                .map(|cs| cs.iter().map(|c| (c.to_string(), g.clone())).collect())
                .unwrap_or_default(),
        };
        if keys.is_empty() {
            continue;
        }
        if !truth.contains_key(f) {
            truth_order.push(f.clone());
        }
        truth.entry(f.clone()).or_default().extend(keys);
    }
    let t_order: Vec<String> = truth_order
        .into_iter()
        .filter(|f| a.family.as_ref().map_or(true, |sub| truth[f].iter().any(|g| g.1.contains(sub.as_str()))))
        .filter(|f| truth[f].len() >= 2)
        .collect();
    let universe: HashSet<&Gene> = t_order.iter().flat_map(|f| truth[f].iter()).collect();

    // predicted clusters as gene sets, intersected with the truth universe
    let mut p_order: Vec<String> = Vec::new();
    let mut pred: HashMap<String, HashSet<Gene>> = HashMap::new();
    for (cid, members) in &clusters {
        let gs: HashSet<Gene> = members
            .iter()
            .filter_map(|l| locus_gene.get(&(cid.clone(), l.clone())))
            .filter(|g| universe.contains(g))
            .cloned()
            .collect();
        if !gs.is_empty() {
            p_order.push(cid.clone());
            pred.insert(cid.clone(), gs);
        }
    }
    if t_order.is_empty() || p_order.is_empty() {
        return Ok(Score { t_order, truth, p_order, pred, assignment: Vec::new(), m: Vec::new(), collapsed: 0, missing: 0 });
    }

    let m: Vec<Vec<i64>> = t_order
        .iter()
        .map(|tf| p_order.iter().map(|pc| truth[tf].intersection(&pred[pc]).count() as i64).collect())
        .collect();
    let assignment = max_overlap_assignment(&m);

    // collapse (register 817's failure mode): each TRUTH GENE -> the locus on its contig that best covers it;
    // genes that must share one locus. ⚠ gene -> locus, never locus -> gene (0 by construction).
    let mut gene_span: HashMap<Gene, (i64, i64)> = HashMap::new();
    for (c, v) in &spans {
        for (gs, ge, g) in v {
            gene_span.insert((c.clone(), g.clone()), (*gs, *ge)); // last wins on duplicate names, as the Python dict does
        }
    }
    let mut loci_by_contig: HashMap<&str, Vec<(&str, i64, i64)>> = HashMap::new();
    for (cid, ms) in &clusters {
        for (c, s, e) in ms {
            loci_by_contig.entry(c.as_str()).or_default().push((cid.as_str(), *s, *e));
        }
    }
    let mut share: BTreeMap<(&str, &str, i64, i64), usize> = BTreeMap::new();
    let mut placed = 0usize;
    for g in &universe {
        let Some(&(gs, ge)) = gene_span.get(*g) else { continue };
        let mut best: Option<(i64, (&str, &str, i64, i64))> = None;
        for &(cid, s, e) in loci_by_contig.get(g.0.as_str()).map(|v| v.as_slice()).unwrap_or(&[]) {
            let ov = ge.min(e) - gs.max(s);
            if ov > 0 && best.map_or(true, |b| ov > b.0) {
                best = Some((ov, (cid, g.0.as_str(), s, e)));
            }
        }
        if let Some((_, l)) = best {
            *share.entry(l).or_insert(0) += 1;
            placed += 1;
        }
    }
    let collapsed: usize = share.values().filter(|&&n| n > 1).map(|&n| n - 1).sum();
    let missing = universe.len() - placed;
    Ok(Score { t_order, truth, p_order, pred, assignment, m, collapsed, missing })
}

/// Pairwise view over the truth universe: genes interned to ids (sorted gene order), the predicted pair set (a
/// pair is predicted when some scored cluster holds both genes), and the totals. Truth families are disjoint (one
/// family per gene), so the true-positive pairs are also the sum of the per-family recovered pairs.
struct Pairs {
    ids: HashMap<Gene, u32>,
    predicted: HashSet<(u32, u32)>,
    n_truth: usize,
    tp: usize,
}

/// Unordered id pairs `(smaller, larger)` of a gene set.
fn id_pairs(genes: &HashSet<Gene>, ids: &HashMap<Gene, u32>) -> Vec<(u32, u32)> {
    let mut v: Vec<u32> = genes.iter().filter_map(|g| ids.get(g).copied()).collect();
    v.sort_unstable();
    let mut out = Vec::with_capacity(v.len() * v.len().saturating_sub(1) / 2);
    for i in 0..v.len() {
        for j in (i + 1)..v.len() {
            out.push((v[i], v[j]));
        }
    }
    out
}

fn pair_counts(sc: &Score) -> Pairs {
    let mut genes: Vec<&Gene> = sc.t_order.iter().flat_map(|f| sc.truth[f].iter()).collect();
    genes.sort();
    genes.dedup();
    let ids: HashMap<Gene, u32> = genes.iter().enumerate().map(|(i, g)| ((*g).clone(), i as u32)).collect();
    let predicted: HashSet<(u32, u32)> = sc.p_order.iter().flat_map(|c| id_pairs(&sc.pred[c], &ids)).collect();
    let (mut n_truth, mut tp) = (0usize, 0usize);
    for f in &sc.t_order {
        for p in id_pairs(&sc.truth[f], &ids) {
            n_truth += 1;
            tp += usize::from(predicted.contains(&p));
        }
    }
    Pairs { ids, predicted, n_truth, tp }
}

fn gene_label(g: &Gene, genome_wide: bool) -> String {
    if genome_wide { format!("{}:{}", g.0, g.1) } else { g.1.clone() }
}

/// `--per-family`: one row per truth family in first-seen order (the definitions of `score_arm`).
fn write_per_family(path: &str, sc: &Score, pairs: &Pairs, genome_wide: bool) -> Result<()> {
    let assigned: HashMap<usize, usize> = sc.assignment.iter().copied().collect();
    let mut fh = std::io::BufWriter::new(std::fs::File::create(path).with_context(|| format!("creating {path}"))?);
    writeln!(
        fh,
        "family_id\tn_truth\tcluster\tn_pred\thit\tsens\tprec\tf\tjaccard\ttruth_pairs\ttp_pairs\tn_contigs\tmembers\thit_members"
    )?;
    for (i, fam) in sc.t_order.iter().enumerate() {
        let tset = &sc.truth[fam];
        let nt = tset.len();
        let hit_j = assigned.get(&i).copied().filter(|&j| sc.m.get(i).is_some_and(|r| r[j] > 0));
        let (cluster, npred, hit, members_hit): (String, usize, usize, Vec<&Gene>) = match hit_j {
            Some(j) => {
                let c = &sc.p_order[j];
                let mut h: Vec<&Gene> = tset.intersection(&sc.pred[c]).collect();
                h.sort();
                (c.clone(), sc.pred[c].len(), sc.m[i][j] as usize, h)
            }
            None => ("-".into(), 0, 0, Vec::new()),
        };
        let (s, p) = if hit > 0 { (hit as f64 / nt as f64, hit as f64 / npred as f64) } else { (0.0, 0.0) };
        let f = if s + p > 0.0 { 2.0 * s * p / (s + p) } else { 0.0 };
        let jac = if hit > 0 { hit as f64 / (nt + npred - hit) as f64 } else { 0.0 };
        let fam_pairs = id_pairs(tset, &pairs.ids);
        let tp_pairs = fam_pairs.iter().filter(|p| pairs.predicted.contains(*p)).count();
        let mut members: Vec<&Gene> = tset.iter().collect();
        members.sort();
        let n_contigs = members.iter().map(|g| g.0.as_str()).collect::<BTreeSet<_>>().len();
        let join = |v: &[&Gene]| v.iter().map(|g| gene_label(g, genome_wide)).collect::<Vec<_>>().join(",");
        writeln!(
            fh,
            "{fam}\t{nt}\t{cluster}\t{npred}\t{hit}\t{s:.6}\t{p:.6}\t{f:.6}\t{jac:.6}\t{}\t{tp_pairs}\t{n_contigs}\t{}\t{}",
            fam_pairs.len(),
            join(&members),
            if members_hit.is_empty() { "-".to_string() } else { join(&members_hit) }
        )?;
    }
    fh.flush()?;
    Ok(())
}

fn main() -> Result<()> {
    let a = Args::parse();
    let genome_wide = a.chrom == ALL;
    let sc = score(&a)?;
    // pairs are counted only when an opt-in output asks for them
    let pairs = (a.pairwise || a.per_family.is_some()).then(|| pair_counts(&sc));
    if let (Some(p), Some(pr)) = (&a.per_family, &pairs) {
        write_per_family(p, &sc, pr, genome_wide)?;
    }
    if sc.t_order.is_empty() || sc.p_order.is_empty() {
        println!("{}: no scoreable truth/prediction overlap", a.label);
    } else {
        let matched: i64 = sc.assignment.iter().map(|&(i, j)| sc.m[i][j]).sum();
        let tot_truth: usize = sc.t_order.iter().map(|f| sc.truth[f].len()).sum();
        let tot_pred: usize = sc
            .assignment
            .iter()
            .filter(|&&(i, j)| sc.m[i][j] > 0)
            .map(|&(_, j)| sc.pred[&sc.p_order[j]].len())
            .sum();
        let sens = if tot_truth > 0 { matched as f64 / tot_truth as f64 } else { 0.0 };
        let prec = if tot_pred > 0 { matched as f64 / tot_pred as f64 } else { 0.0 };
        let f1 = if sens + prec > 0.0 { 2.0 * sens * prec / (sens + prec) } else { 0.0 };
        println!(
            "{:>14} | truth {} fams / {} genes | clusters {} | sens {:.3} prec {:.3} F {:.3} | collapsed {} | no-locus {}",
            a.label, sc.t_order.len(), tot_truth, sc.p_order.len(), sens, prec, f1, sc.collapsed, sc.missing
        );
    }
    if let (true, Some(pr)) = (a.pairwise, &pairs) {
        let (n_tpairs, n_ppairs, n_tp) = (pr.n_truth, pr.predicted.len(), pr.tp);
        let ps = if n_tpairs > 0 { n_tp as f64 / n_tpairs as f64 } else { 0.0 };
        let pp = (n_ppairs > 0).then(|| n_tp as f64 / n_ppairs as f64);
        let pf = pp.map(|p| if ps + p > 0.0 { 2.0 * ps * p / (ps + p) } else { 0.0 });
        let fmt = |x: Option<f64>| x.map_or("NA".to_string(), |v| format!("{v:.3}"));
        println!(
            "{:>14} | pairwise | truth pairs {n_tpairs} | predicted pairs {n_ppairs} | tp {n_tp} | sens {ps:.3} prec {} F {}",
            a.label,
            fmt(pp),
            fmt(pf)
        );
    }
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

    // ---- genome-wide mode and the opt-in outputs (2026-09-25) ----

    fn tmp(tag: &str, body: &str) -> String {
        let p = std::env::temp_dir().join(format!("family_score_test_{}_{tag}", std::process::id()));
        std::fs::write(&p, body).unwrap();
        p.display().to_string()
    }

    fn args(clusters: &str, gff: &str, truth: &str, chrom: &str) -> Args {
        Args {
            clusters: clusters.into(),
            gff: gff.into(),
            soto: truth.into(),
            chrom: chrom.into(),
            family: None,
            label: "t".into(),
            per_family: None,
            pairwise: false,
        }
    }

    fn gff_line(c: &str, s: i64, e: i64, name: &str) -> String {
        format!("{c}\tRefSeq\tgene\t{s}\t{e}\t.\t+\t.\tID=g-{name}-{c}-{s};Name={name}\n")
    }

    /// (families, truth genes, scored clusters, matched, matched-cluster members, collapsed, no-locus)
    fn summary(sc: &Score) -> (usize, usize, usize, i64, usize, usize, usize) {
        let matched: i64 = sc.assignment.iter().map(|&(i, j)| sc.m[i][j]).sum();
        let tot_pred: usize =
            sc.assignment.iter().filter(|&&(i, j)| sc.m[i][j] > 0).map(|&(_, j)| sc.pred[&sc.p_order[j]].len()).sum();
        let tot_truth: usize = sc.t_order.iter().map(|f| sc.truth[f].len()).sum();
        (sc.t_order.len(), tot_truth, sc.p_order.len(), matched, tot_pred, sc.collapsed, sc.missing)
    }

    #[test]
    fn genome_wide_on_one_contig_is_the_per_chromosome_mode() {
        let mut gff = String::from("##gff-version 3\n");
        for (s, e, n) in [(100, 900, "A"), (1_000, 1_900, "B"), (3_000, 3_900, "C"), (5_000, 5_900, "D"),
                          (7_000, 7_900, "E"), (9_000, 9_500, "A"), (11_000, 11_900, "F")] {
            gff += &gff_line("chr1", s, e, n);
        }
        let gff = tmp("one.gff", &gff);
        let cl = tmp(
            "one.clusters.tsv",
            "cluster_id\tchrom\tstart\tend\n\
             K1\tchr1\t150\t850\nK1\tchr1\t1050\t1850\nK2\tchr1\t3050\t3850\n\
             K3\tchr1\t5050\t5850\nK3\tchr1\t7050\t7850\nK3\tchr1\t10950\t11850\nK4\tchr1\t9050\t9450\n",
        );
        let truth = tmp("one.truth.tsv", "Gene Name\tFamily ID\nA\tf1\nB\tf1\nC\tf1\nD\tf2\nE\tf2\nZ\tf3\nF\tf2\n");
        let per = score(&args(&cl, &gff, &truth, "chr1")).unwrap();
        let gw = score(&args(&cl, &gff, &truth, ALL)).unwrap();
        assert_eq!(summary(&per), summary(&gw));
        assert_eq!(per.t_order, gw.t_order);
        assert_eq!(per.p_order, gw.p_order);
        assert_eq!(per.assignment, gw.assignment);
        assert_eq!(summary(&per), (2, 6, 4, 5, 5, 0, 0));
    }

    /// `tag` keeps the files of tests running in parallel apart.
    fn par_fixture(tag: &str) -> (String, String) {
        // P is a pseudoautosomal-like gene: the same Name on chrX and chrY
        let gff = tmp(
            &format!("{tag}.par.gff"),
            &(gff_line("chrX", 100, 200, "P") + &gff_line("chrY", 100, 200, "P") + &gff_line("chr1", 1_000, 2_000, "Q")),
        );
        let cl = tmp(
            &format!("{tag}.par.clusters.tsv"),
            "cluster_id\tchrom\tstart\tend\nK1\tchrX\t100\t200\nK1\tchr1\t1000\t2000\nK2\tchrY\t100\t200\n",
        );
        (gff, cl)
    }

    #[test]
    fn genome_wide_keys_a_repeated_name_by_contig_and_crosses_contigs() {
        let (gff, cl) = par_fixture("keys");
        // name-only truth: the Name stands for every contig's gene of that name
        let truth = tmp("par.truth.tsv", "Gene Name\tFamily ID\nP\tF\nQ\tF\n");
        let sc = score(&args(&cl, &gff, &truth, ALL)).unwrap();
        let mut fam: Vec<&Gene> = sc.truth["F"].iter().collect();
        fam.sort();
        assert_eq!(
            fam,
            vec![&("chr1".to_string(), "Q".to_string()), &("chrX".into(), "P".into()), &("chrY".into(), "P".into())]
        );
        assert_eq!(summary(&sc), (1, 3, 2, 2, 2, 0, 0), "F -> K1 (chrX:P + chr1:Q); chrY:P sits in K2");
        let pr = pair_counts(&sc);
        assert_eq!((pr.n_truth, pr.predicted.len(), pr.tp), (3, 1, 1));
        // a contig column names one (contig, Name) gene per row
        let truth_c = tmp("par.truth_c.tsv", "Gene Name\tFamily ID\tChrom\nP\tF\tchrX\nQ\tF\tchr1\n");
        let sc = score(&args(&cl, &gff, &truth_c, ALL)).unwrap();
        assert_eq!(summary(&sc), (1, 2, 1, 2, 2, 0, 0), "exactly K1; chrY:P is outside the universe");
        // the per-chromosome mode ignores the contig column and is unchanged: on chrX, P alone is no family
        let sc = score(&args(&cl, &gff, &truth_c, "chrX")).unwrap();
        assert!(sc.t_order.is_empty());
    }

    #[test]
    fn per_family_rows_carry_the_match_the_pairs_and_the_members() {
        let (gff, cl) = par_fixture("rows");
        let truth = tmp("par2.truth.tsv", "Gene Name\tFamily ID\nP\tF\nQ\tF\n");
        let sc = score(&args(&cl, &gff, &truth, ALL)).unwrap();
        let out = std::env::temp_dir().join(format!("family_score_test_{}_per_family.tsv", std::process::id()));
        write_per_family(out.to_str().unwrap(), &sc, &pair_counts(&sc), true).unwrap();
        let text = std::fs::read_to_string(&out).unwrap();
        let rows: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
        assert_eq!(rows[0][0], "family_id");
        assert_eq!(rows.len(), 2);
        assert_eq!(
            rows[1],
            vec!["F", "3", "K1", "2", "2", "0.666667", "1.000000", "0.800000", "0.666667", "3", "1", "3",
                 "chr1:Q,chrX:P,chrY:P", "chr1:Q,chrX:P"]
        );
    }

}
