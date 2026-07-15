//! Assembly-based paralog-specific copy number (parCN) orchestrator CLI. OPTIONAL assembly/DNA-side
//! supplement: projects a catalog `copies.fa` onto two phased haplotype assemblies (mat/pat) and counts
//! per-copy genomic loci, disambiguated by deterministic SUN witnesses. See
//! docs/superpowers/specs/2026-07-14-assembly-parcn-design.md and `rustle::vg_family::parcn`.
//!
//! Writes `<out>.parcn.tsv` (per-copy: SUN tier, per-haplotype locus counts, parCN, assign method) and
//! `<out>.parcn_families.tsv` (per-family: copy count, diploid famCN, unresolved-locus count).

use anyhow::{Context, Result};
use clap::Parser;
use std::collections::{BTreeMap, HashMap};
use std::io::Write;

use rustle::vg_family::genome_projection::project_with_cs;
use rustle::vg_family::parcn::{
    assign_locus, dedup_loci, format_family_row, format_parcn_row, parse_copies_fa, sun_positions,
    tabulate, Assignment, CopySun, Locus,
};

#[derive(Parser, Debug)]
#[command(
    name = "parcn",
    about = "Assembly-based paralog-specific copy number: project copies.fa onto phased mat/pat haplotypes"
)]
struct Args {
    /// Catalog `copies.fa` (`>{family}|{copy_idx}|{chrom}:{s}-{e}|{strand}|nexon={n}`), e.g. from
    /// `gw_family_catalog --enumerate-copies`.
    #[arg(long)]
    copies_fa: String,
    /// Maternal haplotype assembly FASTA.
    #[arg(long)]
    mat: String,
    /// Paternal haplotype assembly FASTA.
    #[arg(long)]
    pat: String,
    /// Output prefix; writes `<out>.parcn.tsv` and `<out>.parcn_families.tsv`.
    #[arg(long)]
    out: String,
    /// minimap2 binary (honors PATH; override for a non-standard install).
    #[arg(long, default_value = "minimap2")]
    minimap2: String,
    /// minimap2 threads.
    #[arg(long, default_value_t = 4)]
    threads: usize,
}

pub struct RunArgs {
    pub copies_fa: String,
    pub mat: String,
    pub pat: String,
    pub out: String,
    pub minimap2: String,
    pub threads: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();
    run(RunArgs {
        copies_fa: args.copies_fa,
        mat: args.mat,
        pat: args.pat,
        out: args.out,
        minimap2: args.minimap2,
        threads: args.threads,
    })
}

/// `band` for `sun_positions`: the family's copy-length spread + a fixed slack, so the banded MSA the SUN
/// scan relies on can reach every copy pair regardless of indel-driven length differences. Capped at 1024:
/// an unbounded band on a family with a wide copy-length spread (e.g. a retained-intron outlier) can blow
/// up the banded DP matrix's memory; beyond the cap `banded_msa_pair` returns `None` for that pair, which
/// `diff_offsets` already handles conservatively (empty diff set -> no spurious private positions for that
/// pair, just fewer/no confirmable SUNs for that copy).
fn band_for(copies: &[rustle::vg_family::parcn::Copy]) -> usize {
    let lens = copies.iter().map(|c| c.seq.len());
    let (lo, hi) = lens.fold((usize::MAX, 0usize), |(lo, hi), l| (lo.min(l), hi.max(l)));
    // Band = the family's copy-length spread + slack, so `banded_msa_pair` can align copies that differ in
    // length (partial/alternative-exon copies) without clipping. Capped at 8192 to bound the DP matrix on a
    // pathological huge-disparity family (~O(len*8192) ≈ a few GB): real within-family spreads are small
    // (observed max 5651 bp over the gorilla catalog's 157 families, so 8192 clips none of them), and
    // `sun_positions` runs before the projection index loads so its memory never overlaps the ~13 GB index.
    // A spread above the cap yields `None` from the banded aligner → the conservative empty-diff path.
    (if hi < lo { 64 } else { (hi - lo) + 64 }).min(8192)
}

/// Project one haplotype's `queries` onto `target`, group hits by family (`qname` = `"{family}|{copy}"`,
/// split once on the FIRST `|`), dedup overlapping loci within each family, and assign each deduped locus
/// to a copy via its SUN witness. Returns per-family `Vec<Assignment>`. Empty `project_with_cs` output
/// (missing minimap2, no hits) degrades gracefully to an empty map for every family (zero loci on this side).
fn project_and_assign(
    queries: &[(String, Vec<u8>)],
    target: &str,
    minimap2: &str,
    threads: usize,
    suns_by_fam: &BTreeMap<String, Vec<CopySun>>,
) -> Result<HashMap<String, Vec<Assignment>>> {
    let hits = project_with_cs(queries, target, 0.95, 0.90, minimap2, threads)?;
    if hits.is_empty() {
        eprintln!("[parcn] WARNING: no projection hits against haplotype target '{target}' (0 loci this side)");
        return Ok(HashMap::new());
    }
    let mut by_fam: HashMap<String, Vec<Locus>> = HashMap::new();
    for h in &hits {
        let Some((fam, copy)) = h.qname.split_once('|') else { continue };
        by_fam.entry(fam.to_string()).or_default().push(Locus {
            chrom: h.chrom.clone(),
            start: h.start,
            end: h.end,
            best_copy: copy.to_string(),
            identity: h.identity,
            runner_up_identity: 0.0,
            cs: h.cs.clone(),
            qs: h.qs,
            qe: h.qe,
            strand: h.strand,
        });
    }
    let mut out: HashMap<String, Vec<Assignment>> = HashMap::new();
    for (fam, loci) in by_fam {
        let Some(suns) = suns_by_fam.get(&fam) else { continue };
        let deduped = dedup_loci(loci);
        let mut assignments = Vec::with_capacity(deduped.len());
        for locus in &deduped {
            let Some(sun) = suns.iter().find(|s| s.copy_id == locus.best_copy) else { continue };
            assignments.push(assign_locus(locus, sun));
        }
        out.insert(fam, assignments);
    }
    Ok(out)
}

pub fn run(a: RunArgs) -> Result<()> {
    // Fail fast with a clear message: project_with_cs degrades gracefully (empty hits) on a missing
    // minimap2, but the binary should not silently emit an all-unresolved report.
    if std::process::Command::new(&a.minimap2).arg("--version").output().is_err() {
        anyhow::bail!("minimap2 ('{}') not found on PATH — parcn requires minimap2", a.minimap2);
    }

    let fams = parse_copies_fa(&a.copies_fa).with_context(|| format!("parsing {}", a.copies_fa))?;

    // Query list (once, both haplotypes reuse it), keyed by "family|copy".
    let mut queries: Vec<(String, Vec<u8>)> = Vec::new();
    for (fam_id, copies) in &fams {
        for c in copies {
            queries.push((format!("{fam_id}|{}", c.copy_id), c.seq.clone()));
        }
    }

    let suns_by_fam: BTreeMap<String, Vec<CopySun>> =
        fams.iter().map(|(f, c)| (f.clone(), sun_positions(c, band_for(c)))).collect();

    let mat_by_fam = project_and_assign(&queries, &a.mat, &a.minimap2, a.threads, &suns_by_fam)
        .with_context(|| format!("projecting onto maternal haplotype {}", a.mat))?;
    let pat_by_fam = project_and_assign(&queries, &a.pat, &a.minimap2, a.threads, &suns_by_fam)
        .with_context(|| format!("projecting onto paternal haplotype {}", a.pat))?;

    let parcn_path = format!("{}.parcn.tsv", a.out);
    let fam_path = format!("{}.parcn_families.tsv", a.out);
    let mut pw = std::fs::File::create(&parcn_path).with_context(|| format!("creating {parcn_path}"))?;
    let mut fw = std::fs::File::create(&fam_path).with_context(|| format!("creating {fam_path}"))?;
    writeln!(pw, "family_id\tcopy_id\tsun_tier\tloci_mat\tloci_pat\tparCN\tassign_method")?;
    writeln!(fw, "family_id\tn_copies\tfamCN_diploid\tn_unresolved_loci")?;

    let empty: Vec<Assignment> = Vec::new();
    for (fam_id, copies) in &fams {
        let suns = &suns_by_fam[fam_id];
        let mat_assign = mat_by_fam.get(fam_id).unwrap_or(&empty);
        let pat_assign = pat_by_fam.get(fam_id).unwrap_or(&empty);
        let (rows, n_unresolved) = tabulate(fam_id, copies, suns, mat_assign, pat_assign);
        for r in &rows {
            writeln!(pw, "{}", format_parcn_row(r))?;
        }
        writeln!(fw, "{}", format_family_row(fam_id, &rows, n_unresolved))?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parcn_end_to_end_heterozygous_mat_pat_split() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
        let dir = std::env::temp_dir();
        let tag = std::process::id();
        // Deterministic splitmix64-based non-degenerate sequence generator (SAME pattern as
        // vg_family::genome_projection's tests, see that module's `projection_enumerates_disjoint_copies`
        // note): a naive periodic index formula, or worse a homopolymer pad, aliases into a short repeated
        // motif that minimap2 (`-x splice`, no cap on the "intron" gap it will bridge) happily chains
        // across, fusing two distinct copy loci into one spurious alignment and inflating parCN.
        let splitmix = |i: u64| -> u64 {
            let mut z = i.wrapping_add(0x9E3779B97F4A7C15);
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        };
        let bases = [b'A', b'C', b'G', b'T'];
        let gen_seq = |seed: u64, len: u64| -> Vec<u8> {
            (0..len)
                .map(|i| bases[(splitmix(seed.wrapping_mul(0x2545_F491_4F6C_DD1D).wrapping_add(i)) % 4) as usize])
                .collect::<Vec<u8>>()
        };
        // Family F, two DIVERGENT copies (~12.5% substitutions, evenly spread every 8bp starting at offset
        // 4 -> 50 mismatches over 400bp) rather than the single-base difference this test used to plant.
        // `project_with_cs` (see `project_and_assign`) requires identity>=0.95 AND query coverage>=0.90 for
        // a hit to count; an EVEN 1-in-8 substitution rate keeps every contiguous >=90%-of-query window at
        // ~87.5% identity, so there is no clean sub-window real minimap2 could exploit into a spurious
        // cross-copy hit (unlike a single planted SNV, which leaves the two copies ~99.75% identical to
        // EACH OTHER and thus mutually alignable). This is what makes c1's absence from the paternal
        // haplotype (below) actually register as loci_pat=0 instead of accidentally re-discovering c1's
        // locus via a near-identical match to c0.
        // Every substituted offset is, by construction, private to whichever copy carries it (family size
        // 2, so "differs from the one sibling" IS "private" -- see `sun_positions`) -> both copies land
        // Tier-1 (SUN).
        let base = gen_seq(1, 400);
        let c0 = base.clone();
        let mut c1 = base.clone();
        let mut p = 4usize;
        while p < 400 {
            let cur_idx = bases.iter().position(|&b| b == c0[p]).unwrap();
            c1[p] = bases[(cur_idx + 1) % 4];
            p += 8;
        }
        let copies_fa = dir.join(format!("parcn_e2e_copies_{tag}.fa"));
        std::fs::write(&copies_fa, format!(">F|0|c:1-1|+|nexon=1\n{}\n>F|1|c:1-1|+|nexon=1\n{}\n", String::from_utf8_lossy(&c0), String::from_utf8_lossy(&c1))).unwrap();
        // Padding seeded distinctly (seed=99, vs. the copy bodies' seed=1) so it neither aligns to the
        // copies nor to itself in a way that would fake a third locus.
        let pad = gen_seq(99, 300);
        // HETEROZYGOUS fixture: the maternal haplotype carries BOTH copies' loci; the paternal haplotype
        // carries ONLY c0's -- a heterozygous absence of copy1 on paternal. This is exactly what the old
        // fixture never exercised: it put byte-identical content on both haplotypes, so loci_mat and
        // loci_pat -- the tool's headline per-haplotype split -- were always equal and the real orchestrator
        // path that tells them apart never ran end-to-end.
        let mut mat_seq = c0.clone(); mat_seq.extend(&pad); mat_seq.extend(&c1); mat_seq.extend(&pad);
        let mut pat_seq = c0.clone(); pat_seq.extend(&pad);
        let write_hap = |name: &str, seq: &[u8]| {
            let p = dir.join(format!("parcn_e2e_{name}_{tag}.fa"));
            std::fs::write(&p, format!(">h_{name}\n{}\n", String::from_utf8_lossy(seq))).unwrap();
            p
        };
        let mat = write_hap("mat", &mat_seq);
        let pat = write_hap("pat", &pat_seq);
        let out = dir.join(format!("parcn_e2e_out_{tag}"));
        // call the library orchestrator (factor main's body into `run(args) -> Result<()>` so it is testable)
        run(RunArgs { copies_fa: copies_fa.to_string_lossy().into(), mat: mat.to_string_lossy().into(), pat: pat.to_string_lossy().into(), out: out.to_string_lossy().into(), minimap2: "minimap2".into(), threads: 2 }).unwrap();
        let parcn = std::fs::read_to_string(format!("{}.parcn.tsv", out.to_string_lossy())).unwrap();

        // Parse the per-copy row for `cp` and return its tab-split columns:
        // family_id, copy_id, sun_tier, loci_mat, loci_pat, parCN, assign_method.
        let row = |cp: &str| -> Vec<String> {
            parcn
                .lines()
                .find(|l| { let f: Vec<&str> = l.split('\t').collect(); f.len() > 1 && f[0] == "F" && f[1] == cp })
                .unwrap_or_else(|| panic!("no row for copy {cp} in:\n{parcn}"))
                .split('\t')
                .map(|s| s.to_string())
                .collect()
        };
        // c0: present on BOTH haplotypes -> loci_mat=1, loci_pat=1, parCN=2, resolved by SUN.
        let r0 = row("0");
        assert_eq!(r0[3], "1", "c0 loci_mat: {r0:?}");
        assert_eq!(r0[4], "1", "c0 loci_pat: {r0:?}");
        assert_eq!(r0[5], "2", "c0 parCN: {r0:?}");
        assert_eq!(r0[6], "SUN", "c0 method: {r0:?}");
        // c1: maternal-only (heterozygous absence on paternal) -> loci_mat=1, loci_pat=0, parCN=1. This is
        // the headline mat/pat split the whole-branch review flagged as never exercised end-to-end -- assert
        // the individual haplotype columns explicitly, not just the parCN sum.
        let r1 = row("1");
        assert_eq!(r1[3], "1", "c1 loci_mat: {r1:?}");
        assert_eq!(r1[4], "0", "c1 loci_pat: {r1:?}");
        assert_eq!(r1[5], "1", "c1 parCN: {r1:?}");
        assert_eq!(r1[6], "SUN", "c1 method: {r1:?}");
        assert_ne!(r1[3], r1[4], "c1 must show a real mat/pat split: {r1:?}");

        for p in [copies_fa, mat, pat] { std::fs::remove_file(p).ok(); }
        std::fs::remove_file(format!("{}.parcn.tsv", out.to_string_lossy())).ok();
        std::fs::remove_file(format!("{}.parcn_families.tsv", out.to_string_lossy())).ok();
    }
}
