//! Within-locus haplotype phasing for `--vg-phase`. Operates on one copy's
//! read set (a single VG-family bundle). Detects heterozygous sites from the
//! per-read mismatch pileup (B1), phases reads with exact Minimum-Error-
//! Correction (B2), and assigns hp/ps tags. Reference-base-free: alleles are
//! binary (Alt = read carries the site's dominant alt mismatch, Ref = spans
//! the site without it).

use crate::types::BundleRead;

/// A candidate heterozygous site within one copy's read set.
#[derive(Debug, Clone, PartialEq)]
pub struct HetSite {
    pub pos: u64,       // reference position (0-based; matches BundleRead.mismatches)
    pub alt_allele: u8, // dominant alternate base (ASCII A/C/G/T)
    pub n_ref: u32,
    pub n_alt: u32,
}

#[derive(Debug, Clone)]
pub struct PhasingConfig {
    pub min_maf: f64,          // minor-allele-fraction floor (default 0.25)
    pub max_maf: f64,          // ceiling — excludes fixed/ref-artifact sites (default 0.75)
    pub min_allele_reads: u32, // per-allele read floor (default 3)
    pub max_coverage: usize,   // MEC-DP active-read cap (default 15)
    pub ext_hp_frac: f64,      // external-HP-tag precedence threshold (default 0.5)
}

impl Default for PhasingConfig {
    fn default() -> Self {
        PhasingConfig {
            min_maf: 0.25,
            max_maf: 0.75,
            min_allele_reads: 3,
            max_coverage: 15,
            ext_hp_frac: 0.5,
        }
    }
}

/// True if `pos` falls inside any of the read's aligned exon intervals.
fn read_spans(read: &BundleRead, pos: u64) -> bool {
    read.exons.iter().any(|&(s, e)| pos >= s && pos < e)
}

/// Detect biallelic heterozygous sites from the mismatch pileup.
pub fn detect_het_sites(reads: &[BundleRead], cfg: &PhasingConfig) -> Vec<HetSite> {
    use std::collections::HashMap;
    // pos -> (alt_base -> count)
    let mut alt_counts: HashMap<u64, HashMap<u8, u32>> = HashMap::new();
    for r in reads {
        for &(pos, base) in &r.mismatches {
            *alt_counts.entry(pos).or_default().entry(base).or_default() += 1;
        }
    }
    let mut sites: Vec<HetSite> = Vec::new();
    for (&pos, bases) in &alt_counts {
        // Dominant alt base at this position (deterministic tie-break by base value).
        let (&alt_base, &n_alt) = match bases
            .iter()
            .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
        {
            Some(v) => v,
            None => continue,
        };
        let coverage = reads.iter().filter(|r| read_spans(r, pos)).count() as u32;
        if coverage < n_alt {
            continue; // defensive
        }
        let n_ref = coverage - n_alt;
        let denom = (n_ref + n_alt) as f64;
        if denom == 0.0 {
            continue;
        }
        let maf = n_alt as f64 / denom;
        if n_alt >= cfg.min_allele_reads
            && n_ref >= cfg.min_allele_reads
            && maf >= cfg.min_maf
            && maf <= cfg.max_maf
        {
            sites.push(HetSite { pos, alt_allele: alt_base, n_ref, n_alt });
        }
    }
    sites.sort_by(|a, b| a.pos.cmp(&b.pos).then(a.alt_allele.cmp(&b.alt_allele)));
    sites
}

/// A read's allele at each het site: Some(true)=Alt, Some(false)=Ref, None=not covered.
pub type AlleleRow = Vec<Option<bool>>;

/// Build the read × het-site allele matrix (row order matches `reads`).
pub fn allele_matrix(reads: &[BundleRead], sites: &[HetSite]) -> Vec<AlleleRow> {
    reads
        .iter()
        .map(|r| {
            sites
                .iter()
                .map(|s| {
                    if !read_spans(r, s.pos) {
                        None
                    } else {
                        Some(r.mismatches.iter().any(|&(p, b)| p == s.pos && b == s.alt_allele))
                    }
                })
                .collect()
        })
        .collect()
}

/// Cost of a read row against a haplotype allele vector (Hamming over covered sites).
fn row_cost(row: &AlleleRow, hap: &[bool]) -> u32 {
    row.iter()
        .zip(hap)
        .filter_map(|(a, &h)| a.map(|av| if av != h { 1 } else { 0 }))
        .sum()
}

/// Exact MEC over all 2^M haplotype-A allele assignments (M = #sites).
/// Each read greedily joins the cheaper of {hapA, complement(hapA)}.
/// Returns (hapA alleles, side per read [false=A,true=B], total cost).
/// Caller must ensure M is small (<= ~20).
pub fn mec_brute(matrix: &[AlleleRow], n_sites: usize) -> (Vec<bool>, Vec<bool>, u32) {
    let mut best: Option<(Vec<bool>, Vec<bool>, u32)> = None;
    for mask in 0u32..(1u32 << n_sites) {
        let hap_a: Vec<bool> = (0..n_sites).map(|j| (mask >> j) & 1 == 1).collect();
        let hap_b: Vec<bool> = hap_a.iter().map(|&x| !x).collect();
        let mut sides = Vec::with_capacity(matrix.len());
        let mut total = 0u32;
        for row in matrix {
            let ca = row_cost(row, &hap_a);
            let cb = row_cost(row, &hap_b);
            if ca <= cb {
                sides.push(false);
                total += ca;
            } else {
                sides.push(true);
                total += cb;
            }
        }
        if best.as_ref().map_or(true, |b| total < b.2) {
            best = Some((hap_a, sides, total));
        }
    }
    best.unwrap_or((vec![false; n_sites], vec![false; matrix.len()], 0))
}

use std::collections::HashMap;

/// Per-column active reads: indices into the matrix whose allele is Some at this site.
fn active_at(matrix: &[AlleleRow], site: usize, cap: usize) -> Vec<usize> {
    let mut v: Vec<usize> = (0..matrix.len())
        .filter(|&r| matrix[r][site].is_some())
        .collect();
    v.truncate(cap);
    v
}

/// Weighted column cost: each `(class, side)` contributes `weight[class]` to the
/// orientation it disagrees with. Min over the two allele orientations, exact because
/// each site's two haplotype alleles are chosen independently.
fn column_cost_classes(
    classes: &[AlleleRow],
    weight: &[u32],
    site: usize,
    sides: &[(usize, bool)],
) -> u32 {
    let mut c_a = 0u32;
    let mut c_b = 0u32;
    for &(c, side1) in sides {
        let av = classes[c][site].unwrap();
        let w = weight[c];
        let target_a = if side1 { false } else { true };
        if av != target_a {
            c_a += w;
        } else {
            c_b += w;
        }
    }
    c_a.min(c_b)
}

/// Exact MEC with read-assignment traceback, bounded by DIPLOID STRUCTURE (not a magic
/// read cap).
///
/// MEC assigns every read to one of two haplotypes; the cost is the number of
/// (read, site) cells disagreeing with the read's haplotype consensus. The intractable
/// part of an exact sweep is the state = side of every read still open across a column
/// boundary, which is `2^(boundary coverage)` — and full-length RNA reads all span the
/// gene, so that is `2^(hundreds)`.
///
/// The fix is **pattern collapse**: two reads with an IDENTICAL allele vector are WLOG on
/// the same haplotype (their per-column contribution is identical, so any split can be
/// merged to the cheaper side without raising cost). We collapse identical rows into
/// weighted classes and run the DP over CLASSES. A clean diploid locus has ~2 classes
/// (the two haplotypes) plus a few error-singletons, so the open-state is tiny and the DP
/// is exact and fast at any coverage. `cap` now bounds the number of distinct CLASSES
/// considered per column (truncating only the rare excess error-patterns) — far more
/// effective than capping raw reads, and it almost never fires on real diploid data.
///
/// Returns `(side, cost)` where `side[r]` is `Some(false)`=hapA, `Some(true)`=hapB, or
/// `None` if read `r` covers no het site. Provably equivalent to `mec_brute` (the collapse
/// is cost-preserving).
pub fn mec_dp(matrix: &[AlleleRow], n_sites: usize, cap: usize) -> (Vec<Option<bool>>, u32) {
    let n_reads = matrix.len();
    if n_sites == 0 || n_reads == 0 {
        return (vec![None; n_reads], 0);
    }

    // --- collapse identical-pattern reads into weighted classes ---
    let mut classes: Vec<AlleleRow> = Vec::new();
    let mut weight: Vec<u32> = Vec::new();
    let mut members: Vec<Vec<usize>> = Vec::new();
    let mut idx: HashMap<AlleleRow, usize> = HashMap::new();
    for (r, row) in matrix.iter().enumerate() {
        let c = *idx.entry(row.clone()).or_insert_with(|| {
            classes.push(row.clone());
            weight.push(0);
            members.push(Vec::new());
            classes.len() - 1
        });
        weight[c] += 1;
        members[c].push(r);
    }

    // last column at which each CLASS is active (honoring the per-column class cap).
    let mut last_active: Vec<isize> = vec![-1; classes.len()];
    for site in 0..n_sites {
        for &c in &active_at(&classes, site, cap) {
            last_active[c] = site as isize;
        }
    }

    use std::collections::BTreeMap;
    type Key = BTreeMap<usize, bool>; // open classes -> side
    type Val = (u32, BTreeMap<usize, bool>);
    let mut dp: HashMap<Key, Val> = HashMap::new();
    dp.insert(BTreeMap::new(), (0, BTreeMap::new()));

    for site in 0..n_sites {
        let active = active_at(&classes, site, cap);
        let mut next: HashMap<Key, Val> = HashMap::new();
        for (assign, (pcost, full)) in &dp {
            let mut continuing: Vec<(usize, bool)> = Vec::with_capacity(active.len());
            let mut entering: Vec<usize> = Vec::with_capacity(active.len());
            for &c in &active {
                match assign.get(&c) {
                    Some(&s) => continuing.push((c, s)),
                    None => entering.push(c),
                }
            }
            let n_combos = 1u32 << entering.len();
            for combo in 0..n_combos {
                let mut sides: Vec<(usize, bool)> = continuing.clone();
                for (i, &c) in entering.iter().enumerate() {
                    sides.push((c, (combo >> i) & 1 == 1));
                }
                let total = pcost + column_cost_classes(&classes, &weight, site, &sides);
                let mut new_assign = assign.clone();
                let mut new_full = full.clone();
                for (i, &c) in entering.iter().enumerate() {
                    let s = (combo >> i) & 1 == 1;
                    new_assign.insert(c, s);
                    new_full.insert(c, s);
                }
                new_assign.retain(|&c, _| last_active[c] > site as isize);
                match next.get_mut(&new_assign) {
                    Some(v) => {
                        if total < v.0 {
                            *v = (total, new_full);
                        }
                    }
                    None => {
                        next.insert(new_assign, (total, new_full));
                    }
                }
            }
        }
        dp = next;
    }

    let (_, full) = dp.into_values().min_by_key(|(c, _)| *c).unwrap_or((0, BTreeMap::new()));

    // expand class sides back to per-read sides
    let mut side: Vec<Option<bool>> = vec![None; n_reads];
    for (&c, &s) in &full {
        for &r in &members[c] {
            side[r] = Some(s);
        }
    }
    // recompute the (weighted) cost from the recovered class assignment
    let mut cost = 0u32;
    for site in 0..n_sites {
        let sides: Vec<(usize, bool)> = active_at(&classes, site, cap)
            .iter()
            .map(|&c| (c, full.get(&c).copied().expect("active class must be assigned")))
            .collect();
        cost += column_cost_classes(&classes, &weight, site, &sides);
    }

    (side, cost)
}

/// Cost-only convenience (used by the equivalence test).
pub fn mec_dp_cost(matrix: &[AlleleRow], n_sites: usize, cap: usize) -> u32 {
    mec_dp(matrix, n_sites, cap).1
}

/// One read's haplotype assignment.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadHaplotype {
    pub read_name_hash: u64,
    pub hp: u8,  // 1 or 2
    pub ps: u32, // phase-set id
}

pub struct PhasingResult {
    pub het_sites: Vec<HetSite>,
    pub assignments: Vec<ReadHaplotype>, // only reads covering >=1 het site
}

/// Phase a copy's reads. Detects hets, runs MEC, assigns every het-covering read
/// to a side by best agreement, derives phase sets from read-overlap connectivity,
/// and applies canonical HP labels (HP1 = larger side; ties -> side with the
/// smallest read_name_hash). Reads covering 0 het sites are omitted (unphased).
pub fn phase_reads(reads: &[BundleRead], cfg: &PhasingConfig) -> PhasingResult {
    let sites = detect_het_sites(reads, cfg);
    if sites.is_empty() {
        return PhasingResult { het_sites: sites, assignments: Vec::new() };
    }
    let matrix = allele_matrix(reads, &sites);
    let n_sites = sites.len();

    // false=A, true=B; None = read covers no het site (unphased).
    let side: Vec<Option<bool>> = if n_sites <= 20 {
        // Small site count: exact brute MEC, then assign each read to its cheaper
        // haplotype (the existing derivation).
        let hap_a: Vec<bool> = mec_brute(&matrix, n_sites).0;
        let hap_b: Vec<bool> = hap_a.iter().map(|&x| !x).collect();
        matrix
            .iter()
            .map(|row| {
                if row.iter().all(|a| a.is_none()) {
                    None
                } else {
                    let ca = row_cost(row, &hap_a);
                    let cb = row_cost(row, &hap_b);
                    Some(cb < ca)
                }
            })
            .collect()
    } else {
        // Many sites: exact coverage-bounded MEC-DP returns the per-read side
        // directly (the brute 2^M is intractable here).
        mec_dp(&matrix, n_sites, cfg.max_coverage).0
    };

    // Phase sets: hets linked iff co-covered by a read. Union-find over sites.
    let mut parent: Vec<usize> = (0..n_sites).collect();
    fn find(p: &mut Vec<usize>, x: usize) -> usize {
        let mut r = x;
        while p[r] != r { r = p[r]; }
        let mut c = x;
        while p[c] != r { let n = p[c]; p[c] = r; c = n; }
        r
    }
    for row in &matrix {
        let covered: Vec<usize> = row.iter().enumerate().filter_map(|(j, a)| a.map(|_| j)).collect();
        for w in covered.windows(2) {
            let (a, b) = (find(&mut parent, w[0]), find(&mut parent, w[1]));
            if a != b { parent[a] = b; }
        }
    }
    let mut comp_id: HashMap<usize, u32> = HashMap::new();
    let mut next_ps: u32 = 1;
    let mut site_ps: Vec<u32> = vec![0; n_sites];
    for j in 0..n_sites {
        let root = find(&mut parent, j);
        let id = *comp_id.entry(root).or_insert_with(|| { let v = next_ps; next_ps += 1; v });
        site_ps[j] = id;
    }

    let mut read_ps: Vec<Option<u32>> = vec![None; reads.len()];
    for (ri, row) in matrix.iter().enumerate() {
        if let Some(j) = row.iter().position(|a| a.is_some()) {
            read_ps[ri] = Some(site_ps[j]);
        }
    }

    let mut ps_sideA: HashMap<u32, Vec<u64>> = HashMap::new();
    let mut ps_sideB: HashMap<u32, Vec<u64>> = HashMap::new();
    for (ri, s) in side.iter().enumerate() {
        let (Some(sd), Some(ps)) = (s, read_ps[ri]) else { continue };
        let h = reads[ri].read_name_hash;
        if *sd { ps_sideB.entry(ps).or_default().push(h); }
        else   { ps_sideA.entry(ps).or_default().push(h); }
    }
    let mut a_is_hp1: HashMap<u32, bool> = HashMap::new();
    let empty: Vec<u64> = Vec::new();
    for ps in 1..next_ps {
        let a = ps_sideA.get(&ps).unwrap_or(&empty);
        let b = ps_sideB.get(&ps).unwrap_or(&empty);
        let a_first = a.iter().min().copied();
        let b_first = b.iter().min().copied();
        let a_hp1 = if a.len() != b.len() {
            a.len() > b.len()
        } else {
            match (a_first, b_first) {
                (Some(x), Some(y)) => x <= y,
                (Some(_), None) => true,
                (None, Some(_)) => false,
                (None, None) => true,
            }
        };
        a_is_hp1.insert(ps, a_hp1);
    }

    let mut assignments: Vec<ReadHaplotype> = Vec::new();
    for (ri, s) in side.iter().enumerate() {
        let (Some(sd), Some(ps)) = (s, read_ps[ri]) else { continue };
        let a_hp1 = *a_is_hp1.get(&ps).unwrap_or(&true);
        let is_side_a = !*sd;
        let hp = if is_side_a == a_hp1 { 1u8 } else { 2u8 };
        assignments.push(ReadHaplotype { read_name_hash: reads[ri].read_name_hash, hp, ps });
    }
    assignments.sort_by(|x, y| x.read_name_hash.cmp(&y.read_name_hash));

    PhasingResult { het_sites: sites, assignments }
}

/// Return a clone of `reads` with `hp_tag`/`ps_tag` populated, ready for
/// `split_bundle_by_phase`. If at least `ext_hp_frac` of reads already carry an
/// external `hp_tag` (e.g. a phased BAM's `HP`/`PS` tags), those are used as-is
/// (external precedence); otherwise internal MEC phasing runs. Reads that end up
/// unphased get `hp_tag = None` (they will not be split off).
pub fn assign_haplotypes(reads: &[BundleRead], cfg: &PhasingConfig) -> Vec<BundleRead> {
    use std::collections::HashMap;
    let n = reads.len();
    let mut out: Vec<BundleRead> = reads.to_vec();
    let n_ext = reads.iter().filter(|r| r.hp_tag.is_some()).count();
    if n > 0 && (n_ext as f64) / (n as f64) >= cfg.ext_hp_frac {
        return out; // external HP tags dominate -> trust the BAM's phasing
    }
    let res = phase_reads(reads, cfg);
    let by_hash: HashMap<u64, (u8, u32)> =
        res.assignments.iter().map(|a| (a.read_name_hash, (a.hp, a.ps))).collect();
    for r in out.iter_mut() {
        match by_hash.get(&r.read_name_hash) {
            Some(&(hp, ps)) => {
                r.hp_tag = Some(hp);
                r.ps_tag = Some(ps);
            }
            None => {
                r.hp_tag = None;
                r.ps_tag = None;
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::BundleRead;
    use std::sync::Arc;

    // Minimal BundleRead builder: one exon span [start,end), given mismatches.
    fn mk_read(name_hash: u64, start: u64, end: u64, mismatches: Vec<(u64, u8)>) -> BundleRead {
        BundleRead {
            read_uid: name_hash,
            read_name: Arc::from(""),
            read_name_hash: name_hash,
            ref_id: None,
            mate_ref_id: None,
            mate_start: None,
            hi: 0,
            ref_start: start,
            ref_end: end,
            exons: vec![(start, end)],
            junctions: Vec::new(),
            junction_valid: Vec::new(),
            junctions_raw: Vec::new(),
            junctions_del: Vec::new(),
            weight: 1.0,
            is_reverse: false,
            strand: '+',
            has_poly_start: false,
            has_poly_end: false,
            has_poly_start_aligned: false,
            has_poly_start_unaligned: false,
            has_poly_end_aligned: false,
            has_poly_end_unaligned: false,
            unaligned_poly_t: 0,
            unaligned_poly_a: 0,
            has_last_exon_polya: false,
            has_first_exon_polyt: false,
            query_length: None,
            clip_left: 0,
            clip_right: 0,
            nh: 1,
            nm: 0,
            de: None,
            md: None,
            insertion_sites: Vec::new(),
            unitig: false,
            unitig_cov: 0.0,
            read_count_yc: 1.0,
            countfrag_len: 0.0,
            countfrag_num: 0.0,
            junc_mismatch_weight: 0.0,
            pair_idx: Vec::new(),
            pair_count: Vec::new(),
            mapq: 60,
            mismatches,
            seq: Vec::new(),
            hp_tag: None,
            ps_tag: None,
            is_primary_alignment: true,
            em_weight_gap: -1.0,
            em_n_sites: 0,
            em_anchored: false,
            em_ev_decisive: false,
        }
    }

    #[test]
    fn detects_balanced_het() {
        let cfg = PhasingConfig::default();
        let mut reads = Vec::new();
        for i in 0..5 {
            reads.push(mk_read(i, 50, 200, vec![(100, b'A')]));
        }
        for i in 5..10 {
            reads.push(mk_read(i, 50, 200, vec![]));
        }
        let sites = detect_het_sites(&reads, &cfg);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].pos, 100);
        assert_eq!(sites[0].alt_allele, b'A');
        assert_eq!(sites[0].n_alt, 5);
        assert_eq!(sites[0].n_ref, 5);
    }

    #[test]
    fn rejects_sequencing_error_low_maf() {
        let cfg = PhasingConfig::default();
        let mut reads = Vec::new();
        reads.push(mk_read(0, 50, 200, vec![(100, b'A')]));
        for i in 1..20 {
            reads.push(mk_read(i, 50, 200, vec![]));
        }
        assert!(detect_het_sites(&reads, &cfg).is_empty());
    }

    #[test]
    fn rejects_fixed_difference_high_maf() {
        let cfg = PhasingConfig::default();
        let reads: Vec<_> = (0..10).map(|i| mk_read(i, 50, 200, vec![(100, b'A')])).collect();
        assert!(detect_het_sites(&reads, &cfg).is_empty());
    }

    #[test]
    fn mec_brute_two_clean_haplotypes() {
        let matrix = vec![
            vec![Some(true), Some(false)],
            vec![Some(true), Some(false)],
            vec![Some(false), Some(true)],
            vec![Some(false), Some(true)],
        ];
        let (_hap, sides, cost) = mec_brute(&matrix, 2);
        assert_eq!(cost, 0);
        assert_eq!(sides[0], sides[1]);
        assert_eq!(sides[2], sides[3]);
        assert_ne!(sides[0], sides[2]);
    }

    #[test]
    fn mec_brute_tolerates_one_error() {
        let matrix = vec![
            vec![Some(true), Some(true)],
            vec![Some(true), Some(false)],
            vec![Some(false), Some(true)],
            vec![Some(false), Some(true)],
        ];
        let (_hap, _sides, cost) = mec_brute(&matrix, 2);
        assert_eq!(cost, 1);
    }

    // Deterministic small pseudo-random matrices; DP must match the brute oracle.
    fn lcg(seed: &mut u64) -> u64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *seed >> 33
    }

    #[test]
    fn dp_matches_brute_on_random_small() {
        let mut seed = 0x1234_5678u64;
        for _ in 0..200 {
            let n_sites = 1 + (lcg(&mut seed) % 4) as usize; // 1..4 sites
            let n_reads = 2 + (lcg(&mut seed) % 6) as usize;  // 2..7 reads
            let mut matrix: Vec<AlleleRow> = Vec::new();
            for _ in 0..n_reads {
                let row: AlleleRow = (0..n_sites)
                    .map(|_| match lcg(&mut seed) % 3 {
                        0 => None,
                        1 => Some(false),
                        _ => Some(true),
                    })
                    .collect();
                matrix.push(row);
            }
            let (_, _, brute_cost) = mec_brute(&matrix, n_sites);
            let (dp_side, dp_cost) = mec_dp(&matrix, n_sites, 16);
            assert_eq!(dp_cost, brute_cost, "matrix={:?}", matrix);
            // The returned read assignment must itself achieve the brute optimum
            // (the DP folds the per-site orientation min, so cost equality is the
            // achievement check), and every read covering >=1 site must be placed.
            for (ri, row) in matrix.iter().enumerate() {
                let covered = row.iter().any(|a| a.is_some());
                assert_eq!(
                    dp_side[ri].is_some(),
                    covered,
                    "read {} side/coverage mismatch; matrix={:?}",
                    ri,
                    matrix
                );
            }
        }
    }

    #[test]
    fn mec_dp_scales_on_high_coverage_diploid() {
        // 200 FULL-LENGTH reads over 24 het sites from TWO haplotypes -- every read spans
        // every site, so a per-read cap cannot help and the open-set would be 2^200 -- plus
        // a handful of single-base errors at distinct sites. Pattern-collapse reduces this
        // to ~7 classes, so the EXACT DP is fast. (The previous test used random-allele
        // staggered reads, which is intractable for exact MEC -- NP-hard at 2^coverage;
        // real data has diploid structure, and that structure is the bound.)
        let n_sites = 24usize;
        let hap_a: AlleleRow = (0..n_sites).map(|j| Some(j % 2 == 0)).collect();
        let hap_b: AlleleRow = (0..n_sites).map(|j| Some(j % 2 == 1)).collect();
        let mut matrix: Vec<AlleleRow> = Vec::new();
        let mut n_err = 0u32;
        for k in 0..200usize {
            let mut row = if k % 2 == 0 { hap_a.clone() } else { hap_b.clone() };
            if k % 40 == 7 {
                let s = (k / 40) % n_sites; // distinct error site per read -> distinct patterns
                row[s] = row[s].map(|b| !b);
                n_err += 1;
            }
            matrix.push(row);
        }
        let t = std::time::Instant::now();
        let (side, cost) = mec_dp(&matrix, n_sites, 15);
        assert!(
            t.elapsed().as_secs_f64() < 0.5,
            "mec_dp slow on 200x diploid ({:?}) -- collapse is not bounding the state",
            t.elapsed()
        );
        // each single-base error costs exactly 1 cell against its haplotype consensus
        assert_eq!(cost, n_err, "MEC cost = number of error cells");
        // every read is assigned (NO read dropped, unlike the old read-index cap) ...
        assert!(side.iter().all(|s| s.is_some()), "every covered read assigned");
        // ... and the two haplotypes land on opposite sides.
        assert_ne!(side[0], side[1], "the two haplotypes separate");
    }

    #[test]
    fn assign_haplotypes_external_tags_take_precedence() {
        // 12 reads already carry external HP (6 HP1 + 6 HP2) and NO mismatches, so internal
        // phasing would find nothing — external precedence keeps both haplotypes.
        let cfg = PhasingConfig::default();
        let mut reads = Vec::new();
        for i in 0..6 {
            let mut r = mk_read(i, 50, 200, vec![]);
            r.hp_tag = Some(1);
            reads.push(r);
        }
        for i in 6..12 {
            let mut r = mk_read(i, 50, 200, vec![]);
            r.hp_tag = Some(2);
            reads.push(r);
        }
        let tagged = assign_haplotypes(&reads, &cfg);
        assert!(tagged.iter().any(|r| r.hp_tag == Some(1)));
        assert!(tagged.iter().any(|r| r.hp_tag == Some(2)));
    }

    #[test]
    fn assign_haplotypes_internal_when_no_external() {
        // No external HP -> internal MEC phasing on a clean balanced het at pos 100.
        let cfg = PhasingConfig::default();
        let mut reads = Vec::new();
        for i in 0..6 {
            reads.push(mk_read(i, 50, 200, vec![(100, b'A')]));
        }
        for i in 6..12 {
            reads.push(mk_read(i, 50, 200, vec![]));
        }
        let tagged = assign_haplotypes(&reads, &cfg);
        assert!(tagged.iter().any(|r| r.hp_tag == Some(1)), "one haplotype");
        assert!(tagged.iter().any(|r| r.hp_tag == Some(2)), "the other haplotype");
    }

    #[test]
    fn phase_reads_splits_two_haplotypes() {
        let cfg = PhasingConfig::default();
        let mut reads = Vec::new();
        for i in 0..6 {
            reads.push(mk_read(i, 50, 200, vec![(100, b'A'), (150, b'G')]));
        }
        for i in 6..12 {
            reads.push(mk_read(i, 50, 200, vec![]));
        }
        let res = phase_reads(&reads, &cfg);
        assert_eq!(res.het_sites.len(), 2);
        assert_eq!(res.assignments.len(), 12);
        let hp_alt = res.assignments.iter().find(|a| a.read_name_hash == 0).unwrap().hp;
        let hp_ref = res.assignments.iter().find(|a| a.read_name_hash == 6).unwrap().hp;
        assert_ne!(hp_alt, hp_ref);
        assert!(res.assignments.iter().all(|a| a.ps == res.assignments[0].ps));
    }

    #[test]
    fn phase_reads_unphased_when_no_hets() {
        let cfg = PhasingConfig::default();
        let reads: Vec<_> = (0..10).map(|i| mk_read(i, 50, 200, vec![])).collect();
        let res = phase_reads(&reads, &cfg);
        assert!(res.het_sites.is_empty());
        assert!(res.assignments.is_empty());
    }

    #[test]
    fn phase_reads_deterministic() {
        let cfg = PhasingConfig::default();
        let mut reads = Vec::new();
        for i in 0..6 { reads.push(mk_read(i, 50, 200, vec![(100, b'A')])); }
        for i in 6..12 { reads.push(mk_read(i, 50, 200, vec![])); }
        let a = phase_reads(&reads, &cfg).assignments;
        let b = phase_reads(&reads, &cfg).assignments;
        assert_eq!(a, b);
    }
}
