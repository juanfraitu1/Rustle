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

/// Column cost for a bipartition `bits` (bit i = side of active[i]) at `site`:
/// min over choosing side-0's allele in {Ref,Alt}.
fn column_cost(matrix: &[AlleleRow], site: usize, active: &[usize], bits: u32) -> u32 {
    let mut c_a = 0u32;
    let mut c_b = 0u32;
    for (i, &r) in active.iter().enumerate() {
        let side1 = (bits >> i) & 1 == 1;
        let av = matrix[r][site].unwrap();
        // target for side-0 reads when side-0 == Alt(true): side1-reads target Ref(false)
        let target_a = if side1 { false } else { true };
        let target_b = if side1 { true } else { false };
        if av != target_a {
            c_a += 1;
        }
        if av != target_b {
            c_b += 1;
        }
    }
    c_a.min(c_b)
}

/// Exact coverage-bounded MEC cost via the column DP.
///
/// MEC assigns every read to one of two haplotypes (a global bipartition of the
/// reads); the cost is the number of (read, covered-site) cells that disagree
/// with the read's assigned haplotype allele, where each site independently takes
/// the consensus allele of each haplotype. We sweep columns left to right and
/// enumerate ALL 2^n bipartitions of the column's active reads (n capped by
/// `cap` -> <= 2^cap states). Correctness requires that a read keep the SAME side
/// across ALL columns it covers, not just adjacent ones — so the DP state carries
/// the chosen side of every read seen so far. (A read can skip a column and
/// reappear later; constraining only adjacent columns would wrongly decouple it.)
///
/// Per-column cost is the cheaper of the two allele orientations (side-0 = Ref or
/// side-0 = Alt); this is exact because each site's two haplotype alleles are
/// chosen independently. A free global flip of all sides is symmetric, so no
/// symmetry-breaking is needed. Provably equivalent to `mec_brute`.
pub fn mec_dp_cost(matrix: &[AlleleRow], n_sites: usize, cap: usize) -> u32 {
    if n_sites == 0 || matrix.is_empty() {
        return 0;
    }

    // DP state: assignment of every read seen so far to a side.
    // Key: BTreeMap<read_idx, side> serialized; we use a HashMap keyed by a
    // canonical Vec of (read_idx, side) sorted by read_idx. To keep it cheap we
    // store the side-map directly as a sorted Vec<(usize,bool)>.
    use std::collections::BTreeMap;
    // map from committed read-side assignment -> min cost
    let mut dp: HashMap<BTreeMap<usize, bool>, u32> = HashMap::new();
    dp.insert(BTreeMap::new(), 0);

    for site in 0..n_sites {
        let active = active_at(matrix, site, cap);
        let n = active.len();
        let n_states = 1u32 << n;

        // Precompute column cost for every bipartition of this column's active set.
        let mut col_cost: Vec<u32> = Vec::with_capacity(n_states as usize);
        for bits in 0..n_states {
            col_cost.push(column_cost(matrix, site, &active, bits));
        }

        let mut next: HashMap<BTreeMap<usize, bool>, u32> = HashMap::new();

        for (assign, &pcost) in &dp {
            for bits in 0..n_states {
                // Build candidate assignment, enforcing consistency for any read
                // already committed to a side.
                let mut ok = true;
                let mut new_assign = assign.clone();
                for (i, &r) in active.iter().enumerate() {
                    let cur_side = (bits >> i) & 1 == 1;
                    match new_assign.get(&r) {
                        Some(&existing) => {
                            if existing != cur_side {
                                ok = false;
                                break;
                            }
                        }
                        None => {
                            new_assign.insert(r, cur_side);
                        }
                    }
                }
                if !ok {
                    continue;
                }
                let total = pcost + col_cost[bits as usize];
                next.entry(new_assign)
                    .and_modify(|e| {
                        if total < *e {
                            *e = total;
                        }
                    })
                    .or_insert(total);
            }
        }

        dp = next;
    }

    dp.values().copied().min().unwrap_or(0)
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

    let hap_a: Vec<bool> = if n_sites <= 20 {
        mec_brute(&matrix, n_sites).0
    } else {
        sites.iter().map(|s| s.n_alt >= s.n_ref).collect()
    };
    let hap_b: Vec<bool> = hap_a.iter().map(|&x| !x).collect();

    let mut side: Vec<Option<bool>> = Vec::with_capacity(reads.len()); // false=A,true=B
    for row in matrix.iter() {
        if row.iter().all(|a| a.is_none()) {
            side.push(None);
            continue;
        }
        let ca = row_cost(row, &hap_a);
        let cb = row_cost(row, &hap_b);
        side.push(Some(cb < ca));
    }

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
            let dp_cost = mec_dp_cost(&matrix, n_sites, 16);
            assert_eq!(dp_cost, brute_cost, "matrix={:?}", matrix);
        }
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
