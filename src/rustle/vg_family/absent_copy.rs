//! Admission gate for reference-ABSENT (collapsed) copies.
//!
//! A `CollapsedCandidate` (from `copy_split`) passes through five ordered gates and either
//! becomes an admitted `Copy` (a synthetic `DenovoTranscript` carrying the copy's distinguishing
//! alleles) or a `DnaNeeds` record flagged for downstream DNA-level validation.  All five gates
//! must pass — first failure determines the reason string.
//!
//! Gate order:
//!   1. `n_clusters >= min_clusters` — enough co-varying clusters for a multi-copy call.
//!      Default 3, overridable for removal-ablation experiments ONLY — see
//!      [`AbsentCopyParams::from_env`], which documents why the default must not change.
//!   2. `min_p_distinct` — the candidate's allele vector is certifiably distinct from the host.
//!   3. `strand_symmetric_spectrum` — NOT a pure A→G editing artefact.
//!   4. `collapsed_copy_to_transcript_from_host_seq` — alleles are placeable (substitution-only).
//!   5. Genome remap identity `< remap_max_identity` — the synthetic copy does not remap
//!      back to its host locus with ≥ 98 % identity (i.e. it is genuinely divergent).
//!
//! The hermetic `admit_candidate_with_remap` takes the already-spliced host transcript and an
//! injected `remap_identity` closure so tests never shell minimap2.  The real entry point
//! `admit_candidate` rebuilds the host spliced sequence from the genome (using the copy's intron
//! chain), then delegates to the hermetic variant with a real minimap2 closure.

use std::collections::BTreeMap;

use crate::genome::GenomeIndex;
use crate::vg_family::copy_assign_pipeline::gen2off;
use crate::vg_family::copy_split::{
    collapsed_copy_to_transcript_from_host_seq, min_p_distinct, strand_symmetric_spectrum,
    CollapsedCandidate,
};
use crate::vg_family::family_detect::DenovoTranscript;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// The shipping cluster floor for gate 1. **This value must not change.** See
/// [`AbsentCopyParams::from_env`] for the (experiment-only) override and why it exists.
pub const DEFAULT_MIN_CLUSTERS: usize = 3;

/// Env var overriding [`AbsentCopyParams::min_clusters`]. Unset ⇒ [`DEFAULT_MIN_CLUSTERS`].
pub const MIN_CLUSTERS_ENV: &str = "RUSTLE_ABSENT_MIN_CLUSTERS";

/// Tunable thresholds for the admission gate.
pub struct AbsentCopyParams {
    /// Per-base sequencing error rate used in `min_p_distinct` (default 0.003 = HiFi).
    pub error_rate: f64,
    /// Bonferroni-corrected significance threshold for `min_p_distinct` (default 1e-3).
    pub alpha: f64,
    /// Minimum number of identifiable co-varying clusters at the locus
    /// (default [`DEFAULT_MIN_CLUSTERS`] = 3).
    pub min_clusters: usize,
    /// Upper bound on remap identity; above this threshold the copy maps back too cleanly
    /// to be a genuine divergent paralog (default 0.98 = 98 %).
    pub remap_max_identity: f64,
}

impl Default for AbsentCopyParams {
    fn default() -> Self {
        Self {
            error_rate: 0.003,
            alpha: 1e-3,
            min_clusters: DEFAULT_MIN_CLUSTERS,
            remap_max_identity: 0.98,
        }
    }
}

/// Parse a raw [`MIN_CLUSTERS_ENV`] value into a cluster floor.
///
/// Pure (takes the string, never reads the environment) so the policy below is testable without
/// mutating process-global state — env-mutating tests race under `cargo test`'s thread pool.
///
/// Policy is deliberately FAIL-SAFE: anything not a positive integer (absent, empty, non-numeric,
/// or `0`) yields [`DEFAULT_MIN_CLUSTERS`]. `0` is rejected rather than honoured because
/// `n_clusters < 0` is unsatisfiable, i.e. `0` would silently DISABLE gate 1 outright; a gate that
/// guards against copy-vs-allele confusion must never be switched off by a typo.
pub(crate) fn parse_min_clusters(raw: Option<&str>) -> usize {
    match raw.map(str::trim).and_then(|s| s.parse::<usize>().ok()) {
        Some(n) if n >= 1 => n,
        _ => DEFAULT_MIN_CLUSTERS,
    }
}

impl AbsentCopyParams {
    /// `Default`, but with `min_clusters` overridable via [`MIN_CLUSTERS_ENV`]
    /// (CLI: `copy_assign --absent-min-clusters N`). Every other field is the `Default` value.
    ///
    /// # Why this one threshold is configurable — and why the default MUST stay 3
    ///
    /// Gate 1 (`n_clusters >= min_clusters`) is not a tuning parameter; it is an IDENTIFIABILITY
    /// claim. With fewer than three co-varying clusters at a locus you cannot tell a second COPY
    /// from a heterozygous ALLELE of one copy using RNA alone — that separation needs DNA
    /// copy-number evidence, which is exactly what the rejection reason
    /// `"<N clusters (copy-vs-allele needs DNA)"` says. Lowering it in production would
    /// manufacture false positives by admitting heterozygous sites as reference-absent copies.
    /// **The default must not change. Do not "tune" it to raise recovery.**
    ///
    /// The override exists for ONE situation: an experiment in which copy status is established
    /// BY CONSTRUCTION rather than inferred — specifically the V4b removal-recovery ablation,
    /// where a known copy is DELETED from the assembly and the same reads are re-aligned. There
    /// the gate is answering a question the experimental design has already answered, and it is
    /// additionally UNREACHABLE by arithmetic: deleting one copy of a 3-copy family leaves at
    /// most host + deleted = 2 clusters, so `n_clusters >= 3` can never hold and the locus stops
    /// on gate 1 for a reason that has nothing to do with divergence (this is what blocked GSTM
    /// and GWFAM227).
    ///
    /// Any run using the override is only interpretable ALONGSIDE the identical run against the
    /// INTACT (undeleted) assembly at the SAME `min_clusters`. If the no-deletion control also
    /// admits a copy, the "recovery" is an artefact of the lowered gate, not a recovery.
    pub fn from_env() -> Self {
        let raw = std::env::var(MIN_CLUSTERS_ENV).ok();
        let min_clusters = parse_min_clusters(raw.as_deref());
        // Announce ONCE per process, and only when the var is set, so an unset run prints nothing
        // and stays byte-identical. A lowered gate must be visible in the run log it produced.
        if raw.is_some() {
            static NOTICE: std::sync::Once = std::sync::Once::new();
            NOTICE.call_once(|| {
                eprintln!(
                    "[absent_copy] {MIN_CLUSTERS_ENV}={:?} -> min_clusters={min_clusters} \
                     (default {DEFAULT_MIN_CLUSTERS}). Gate 1 guards copy-vs-allele; a run below \
                     the default is valid ONLY next to the intact-assembly control at the same value.",
                    raw.as_deref().unwrap_or("")
                );
            });
        }
        Self { min_clusters, ..Self::default() }
    }
}

/// A candidate that failed one or more gates and requires DNA-level validation.
#[derive(Clone, Debug)]
pub struct DnaNeedsRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_clusters: usize,
    /// Human-readable string identifying the first gate that fired.
    pub reason: String,
    pub read_count: usize,
}

/// Result of the admission gate.
pub enum Admission {
    /// All gates passed — the `DenovoTranscript` is the synthetic collapsed copy. The second
    /// field is the genome remap identity computed at gate 5 (`Some(id)`); `None` only where an
    /// admitted copy has no computed identity.
    Copy(DenovoTranscript, Option<f64>),
    /// At least one gate failed — the record is flagged for DNA validation.
    DnaNeeds(DnaNeedsRecord),
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn dna_needs(c: &CollapsedCandidate, reason: &str) -> Admission {
    Admission::DnaNeeds(DnaNeedsRecord {
        chrom: c.chrom.clone(),
        start: c.start,
        end: c.end,
        n_clusters: c.n_clusters,
        reason: reason.into(),
        read_count: c.iso.read_count,
    })
}

/// Derive the host's own bases at every PSV position from the already-spliced `host_spliced`
/// transcript (whose `gen2off` map covers exactly the exons built with the copy's intron chain).
fn host_alleles_at(host_spliced: &DenovoTranscript, psv_pos: &[u64]) -> Vec<Option<u8>> {
    let g2o: BTreeMap<u64, usize> = gen2off(host_spliced);
    psv_pos
        .iter()
        .map(|&pos| g2o.get(&pos).and_then(|&off| host_spliced.seq.get(off).copied()))
        .collect()
}

// ---------------------------------------------------------------------------
// Core gate logic (hermetic — remap is injected)
// ---------------------------------------------------------------------------

/// Hermetic admission gate.  `host_spliced` must already carry the COPY's intron chain and
/// the corresponding spliced sequence (so `gen2off` agrees with `host_spliced.seq`).
/// `remap_identity` is a closure `|seq: &[u8]| -> Option<f64>` that returns the best
/// alignment identity against the genome (or `None` if there is no alignment).  In tests,
/// pass `|_seq| Some(0.5)` to bypass minimap2 entirely.
pub(crate) fn admit_candidate_with_remap<F>(
    cand: &CollapsedCandidate,
    host_spliced: &DenovoTranscript,
    p: &AbsentCopyParams,
    remap_identity: F,
) -> Admission
where
    F: Fn(&[u8]) -> Option<f64>,
{
    // Gate 1: cluster count.
    if cand.n_clusters < p.min_clusters {
        return dna_needs(
            cand,
            &format!("<{} clusters (copy-vs-allele needs DNA)", p.min_clusters),
        );
    }

    // Gate 2: min_p distinct from host.
    let host_alleles = host_alleles_at(host_spliced, &cand.psv_pos);
    if !min_p_distinct(&cand.iso.allele_vector, &host_alleles, p.error_rate, p.alpha) {
        return dna_needs(cand, "not min_p-distinct from host");
    }

    // Gate 3: strand-symmetric spectrum (reject pure A→G editing artefact).
    if !strand_symmetric_spectrum(&host_alleles, &cand.iso.allele_vector) {
        return dna_needs(cand, "pure A->G spectrum (editing-like)");
    }

    // Gate 4: build the overlay transcript.
    let t = match collapsed_copy_to_transcript_from_host_seq(
        &cand.iso,
        &cand.psv_pos,
        host_spliced,
    ) {
        Some(t) => t,
        None => return dna_needs(cand, "consensus unplaceable (indel/intron)"),
    };

    // Gate 5: genome remap identity.
    match remap_identity(&t.seq) {
        Some(id) if id < p.remap_max_identity => Admission::Copy(t, Some(id)),
        Some(_) => dna_needs(cand, ">=98% remap identity (paralog-leak or het)"),
        None => dna_needs(cand, "no homology on remap"),
    }
}

// ---------------------------------------------------------------------------
// Real entry point
// ---------------------------------------------------------------------------

/// Admission gate with real genome-based host-seq reconstruction and minimap2 remap.
///
/// `fasta_path` is the genome FASTA file (required by minimap2; `GenomeIndex` does not expose
/// its path so it is threaded in separately).  If the host spliced sequence cannot be built
/// (e.g. non-canonical junction in the copy's intron chain), routes to `DnaNeeds`.
pub fn admit_candidate(
    cand: &CollapsedCandidate,
    host: &DenovoTranscript,
    genome: &GenomeIndex,
    fasta_path: &str,
    p: &AbsentCopyParams,
) -> Admission {
    use crate::vg_family::denovo_assemble::build_spliced_seq;

    // Rebuild host spliced seq using the COPY's intron chain (not the host's own chain), so that
    // `gen2off` inside the hermetic gate agrees with the seq bytes.
    let (seq, strand) = match build_spliced_seq(
        genome,
        &host.chrom,
        host.start,
        host.end,
        &cand.iso.intron_chain,
        // The host's strand is already KNOWN from the catalog, so pass it rather than None: with
        // RUSTLE_READ_STRAND on, an unspliced host would otherwise be rebuilt as the '+' placeholder
        // and its bytes would disagree with the reverse-complemented catalog sequence that `gen2off`
        // is computed against.
        Some(host.strand),
    ) {
        Some(v) => v,
        None => return dna_needs(cand, "host sequence unbuildable"),
    };
    let host_spliced = DenovoTranscript {
        seq,
        strand,
        introns: cand.iso.intron_chain.clone(),
        ..host.clone()
    };

    let fp = fasta_path.to_owned();
    admit_candidate_with_remap(cand, &host_spliced, p, move |seq| {
        remap_identity_minimap2(seq, &fp)
    })
}

// ---------------------------------------------------------------------------
// minimap2 remap helper
// ---------------------------------------------------------------------------

/// Shell minimap2 (`-cx splice --eqx`) to remap `query_seq` against `fasta_path`.
/// The query is a SPLICED copy consensus (exons concatenated); the target genome contains introns, so the
/// SPLICED preset (`-x splice`) is required — `asm20` (non-spliced) cannot chain across real multi-kb introns
/// and would fail to align a genuine multi-exon copy (routing real copies wrongly to DnaNeeds). With splice,
/// the `de:f:` divergence reflects the EXONIC divergence (introns are spliced out, not counted as mismatch).
/// Returns `Some(identity)` for the best alignment (by alignment block length), `None` if there
/// is no alignment.  Identity is derived from the `de:f:<x>` divergence tag when present
/// (`identity = 1 − x`), falling back to `matches / aln_block_len` (PAF cols 10/11).
///
/// PERFORMANCE (H4): a genome-wide O4 pass calls this once per candidate copy, and minimap2 re-indexes
/// the whole genome FASTA on every call (~seconds for a 3 Gb gorilla genome). Set `RUSTLE_ABSENT_MMI` to a
/// pre-built splice index (`minimap2 -x splice -d genome.splice.mmi genome.fasta`) and it is used as the
/// target instead of re-reading the FASTA — minimap2 loads the `.mmi` directly, collapsing per-call indexing
/// to a one-time cost. When the env var is unset the target is `fasta_path` exactly as before (byte-identical
/// behaviour; the `.mmi` must be built with the SAME `-x splice` preset or minimap2's k/w would mismatch).
fn remap_identity_minimap2(query_seq: &[u8], fasta_path: &str) -> Option<f64> {
    use std::io::Write;
    use std::sync::atomic::{AtomicUsize, Ordering};

    // Process-global monotonic counter so concurrent calls (even with equal-length queries) get
    // disjoint temp paths; `pid` keeps distinct processes disjoint from each other.
    static REMAP_NONCE: AtomicUsize = AtomicUsize::new(0);

    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    let nonce = REMAP_NONCE.fetch_add(1, Ordering::Relaxed);
    let qpath = dir.join(format!("rustle_absent_q_{pid}_{nonce}.fa"));

    struct Cleanup(std::path::PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }
    let _cl = Cleanup(qpath.clone());

    {
        let mut q = std::fs::File::create(&qpath).ok()?;
        q.write_all(b">q\n").ok()?;
        q.write_all(query_seq).ok()?;
        q.write_all(b"\n").ok()?;
    }

    // Use a pre-built splice `.mmi` when provided (one-time index instead of per-call genome re-read);
    // defaults to `fasta_path` so behaviour is byte-identical when the env var is unset (H4).
    let target = std::env::var("RUSTLE_ABSENT_MMI").unwrap_or_else(|_| fasta_path.to_string());
    let out = std::process::Command::new(&mm2)
        .args(["-cx", "splice", "--eqx", "--secondary=no", "-t", "1"])
        .arg(&target)
        .arg(&qpath)
        .output()
        .ok()?;

    if !out.status.success() {
        return None;
    }

    let text = String::from_utf8_lossy(&out.stdout);
    let mut best_id: Option<f64> = None;
    let mut best_span: u64 = 0;

    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        // PAF col 10 = number of matching bases; col 11 = alignment block length.
        let aln_block_len: u64 = match f[10].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if aln_block_len == 0 {
            continue;
        }
        // Prefer the de:f: divergence tag (exonic divergence under the spliced alignment).
        let id = if let Some(de_str) = f.iter().find_map(|x| x.strip_prefix("de:f:")) {
            let div: f64 = de_str.parse().unwrap_or(1.0);
            1.0 - div
        } else {
            let matches: u64 = f[9].parse().unwrap_or(0);
            matches as f64 / aln_block_len as f64
        };

        if aln_block_len > best_span {
            best_span = aln_block_len;
            best_id = Some(id);
        }
    }

    best_id
}

// ---------------------------------------------------------------------------
// Augment-and-linearize: minimap2 realign-pool closure (Task 3)
// ---------------------------------------------------------------------------

/// PAF -> per-read primary hit. `n_reads` reads were written with headers `>0..>n_reads`; targets
/// (ref contigs) with headers `>0..`. Primary = the hit carrying the `tp:A:P` tag (fall back to the
/// largest alignment-block length, PAF col 10, when no hit is tagged primary). Returns, per read,
/// `Some((target_contig_index, mapq))` or `None` if the read has no hit at all. PAF is 0-indexed
/// here: f[0]=qname, f[5]=tname, f[10]=alignment-block-len, f[11]=mapq, f[12..]=optional tags.
pub(crate) fn paf_primary_hits(paf: &str, n_reads: usize) -> Vec<Option<(usize, u32)>> {
    // per read: (is_primary, block_len, contig_idx, mapq)
    let mut best: Vec<Option<(bool, u64, usize, u32)>> = vec![None; n_reads];
    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let q: usize = match f[0].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if q >= n_reads {
            continue;
        }
        let t: usize = match f[5].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let block: u64 = f[10].parse().unwrap_or(0);
        let mapq: u32 = f[11].parse().unwrap_or(0);
        let is_p = f[12..].iter().any(|x| *x == "tp:A:P");
        let cand = (is_p, block, t, mapq);
        // Prefer primary; then larger alignment block.
        let take = match best[q] {
            None => true,
            Some((p0, b0, _, _)) => (is_p, block) > (p0, b0),
        };
        if take {
            best[q] = Some(cand);
        }
    }
    best.into_iter().map(|o| o.map(|(_, _, t, mq)| (t, mq))).collect()
}

/// Real realign-pool closure for `linearize_certificate` (Task 2): writes `ref_contigs` as a temp
/// target FASTA (`>{i}` headers) and `reads` as a temp query FASTA (`>{j}` headers), shells
/// minimap2, and reduces the resulting PAF via `paf_primary_hits` to each read's primary
/// `(target_contig_index, mapq)`. Mirrors `remap_identity_minimap2`'s temp-file/`Cleanup`/
/// `RUSTLE_MINIMAP2`/nonce pattern, but targets a POOL of contigs (many-to-many) rather than a
/// single genome FASTA, and uses an intronless preset (`map-hifi`) since both `ref_contigs` and
/// `reads` are spliced-consensus sequences here — there is no intron to chain across (unlike
/// `remap_identity_minimap2`'s genome target, which needs `-x splice`).
///
/// Any spawn/exit failure (minimap2 missing, bad args, non-zero exit) degrades gracefully to
/// `vec![None; reads.len()]`, matching this module's existing `Ok(None)`-on-failure contract.
pub fn realign_pool_minimap2(ref_contigs: &[Vec<u8>], reads: &[Vec<u8>]) -> Vec<Option<(usize, u32)>> {
    use std::io::Write;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static NONCE: AtomicUsize = AtomicUsize::new(0);

    if ref_contigs.is_empty() || reads.is_empty() {
        return vec![None; reads.len()];
    }

    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    let nonce = NONCE.fetch_add(1, Ordering::Relaxed);
    let refp = dir.join(format!("rustle_lin_ref_{pid}_{nonce}.fa"));
    let qp = dir.join(format!("rustle_lin_q_{pid}_{nonce}.fa"));

    struct Cleanup(std::path::PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }
    let _c1 = Cleanup(refp.clone());
    let _c2 = Cleanup(qp.clone());

    let write_fa = |path: &std::path::Path, seqs: &[Vec<u8>]| -> Option<()> {
        let mut fh = std::fs::File::create(path).ok()?;
        for (i, s) in seqs.iter().enumerate() {
            writeln!(fh, ">{}", i).ok()?;
            fh.write_all(s).ok()?;
            fh.write_all(b"\n").ok()?;
        }
        Some(())
    };
    if write_fa(&refp, ref_contigs).is_none() || write_fa(&qp, reads).is_none() {
        return vec![None; reads.len()];
    }

    // Intronless spliced-consensus contigs -> map-hifi (asm20-class), NOT -x splice. --secondary=no
    // so a non-unique read's only reported hit is minimap2's own primary pick; target is the ref
    // contig pool, query is the reads.
    let out = match std::process::Command::new(&mm2)
        .args(["-c", "-x", "map-hifi", "--secondary=no", "-t", "1"])
        .arg(&refp)
        .arg(&qp)
        .output()
    {
        Ok(o) if o.status.success() => o,
        _ => return vec![None; reads.len()],
    };

    paf_primary_hits(&String::from_utf8_lossy(&out.stdout), reads.len())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vg_family::copy_split::{CollapsedCandidate, CopyIsoform};
    use crate::vg_family::family_detect::DenovoTranscript;

    fn host() -> DenovoTranscript {
        DenovoTranscript {
            tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            n_reads: 9,
            strand: '+',
            introns: vec![],
            seq: b"AAAAAAAAAA".to_vec(),
         ..Default::default() }
    }

    /// The shipping cluster floor is 3 and MUST stay 3: below three co-varying clusters a second
    /// COPY is not distinguishable from a heterozygous ALLELE without DNA. This pins the default
    /// so lowering it requires deleting a test that says why it exists, not editing a literal.
    /// The override (`RUSTLE_ABSENT_MIN_CLUSTERS` / `--absent-min-clusters`) is for removal
    /// ablations, where copy status is established by construction — see `from_env`'s docs.
    #[test]
    fn default_min_clusters_is_three_and_must_stay_three() {
        assert_eq!(AbsentCopyParams::default().min_clusters, 3);
        assert_eq!(DEFAULT_MIN_CLUSTERS, 3);
    }

    /// The `MIN_CLUSTERS_ENV` parse policy, exercised WITHOUT touching the process environment
    /// (env-mutating tests race under cargo's thread pool). Fail-safe: absent/blank/non-numeric/`0`
    /// all fall back to the default. `0` in particular must NOT be honoured — `n_clusters < 0` is
    /// unsatisfiable, so it would silently switch gate 1 off entirely.
    #[test]
    fn min_clusters_env_parse_is_fail_safe() {
        assert_eq!(parse_min_clusters(None), 3, "unset -> default");
        assert_eq!(parse_min_clusters(Some("2")), 2, "ablation value honoured");
        assert_eq!(parse_min_clusters(Some(" 2 ")), 2, "whitespace trimmed");
        assert_eq!(parse_min_clusters(Some("1")), 1, "1 is a positive integer");
        assert_eq!(parse_min_clusters(Some("5")), 5, "raising the floor is allowed");
        assert_eq!(parse_min_clusters(Some("0")), 3, "0 would DISABLE gate 1 -> refuse");
        assert_eq!(parse_min_clusters(Some("")), 3, "empty -> default");
        assert_eq!(parse_min_clusters(Some("two")), 3, "non-numeric -> default");
        assert_eq!(parse_min_clusters(Some("-1")), 3, "negative -> default");
        assert_eq!(parse_min_clusters(Some("2.5")), 3, "non-integer -> default");
    }

    /// The override actually moves gate 1: the SAME 2-cluster candidate that `default()` rejects is
    /// admitted at `min_clusters = 2`. Drives the params struct directly (no env), so it documents
    /// the ablation's mechanism without process-global state.
    #[test]
    fn lowered_min_clusters_admits_the_two_cluster_candidate_default_rejects() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C'), Some(b'T'), Some(b'C')],
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            // A 3-copy family with one copy DELETED supplies at most host + deleted = 2.
            n_clusters: 2,
        };
        // Default: stops on gate 1.
        match admit_candidate_with_remap(&cand, &host(), &AbsentCopyParams::default(), |_| Some(0.5)) {
            Admission::DnaNeeds(r) => assert!(r.reason.contains("clusters"), "got: {}", r.reason),
            _ => panic!("default must reject a 2-cluster candidate"),
        }
        // Ablation floor: gate 1 passes and the remaining four gates decide.
        let p = AbsentCopyParams { min_clusters: 2, ..Default::default() };
        match admit_candidate_with_remap(&cand, &host(), &p, |_| Some(0.5)) {
            Admission::Copy(t, id) => {
                assert_eq!(t.chrom, "c1");
                assert_eq!(id, Some(0.5));
            }
            Admission::DnaNeeds(r) => panic!("expected Copy at min_clusters=2, got: {}", r.reason),
        }
    }

    /// Gate 1 fires: only 2 clusters (< min_clusters=3).  The remap closure returns 0.5
    /// (would otherwise be divergent enough to pass gate 5), so the only failure is gate 1.
    #[test]
    fn admit_rejects_two_cluster_substitution_only_as_dna_needs() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C'), Some(b'C'), Some(b'C')],
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            n_clusters: 2, // only 2 clusters
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| Some(0.5));
        match got {
            Admission::DnaNeeds(r) => assert!(
                r.reason.contains("clusters"),
                "expected 'clusters' in reason, got: {}",
                r.reason
            ),
            _ => panic!("expected DnaNeeds"),
        }
    }

    /// All five gates pass: 3 clusters, allele vector is min_p-distinct from host (3 differing
    /// columns → p = (ε/3)³ ≪ α), mixed spectrum (C/T/C not A→G), overlay is buildable
    /// (single-exon identity map), remap returns 0.5 < 0.98.
    #[test]
    fn admit_accepts_three_cluster_distinct_divergent() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C'), Some(b'T'), Some(b'C')],
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            n_clusters: 3,
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| Some(0.5));
        match got {
            Admission::Copy(t, id) => {
                assert_eq!(t.chrom, "c1");
                assert_eq!(id, Some(0.5), "gate 5 must carry through the remap identity");
            }
            _ => panic!("expected Copy"),
        }
    }

    /// Gate 2 fires: allele vector identical to host at all observed columns → not min_p-distinct.
    #[test]
    fn admit_rejects_non_distinct_from_host() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'A'), Some(b'A'), Some(b'A')], // same as host
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            n_clusters: 3,
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| Some(0.5));
        match got {
            Admission::DnaNeeds(r) => assert!(
                r.reason.contains("min_p-distinct"),
                "expected 'min_p-distinct' in reason, got: {}",
                r.reason
            ),
            _ => panic!("expected DnaNeeds"),
        }
    }

    /// Gate 3 fires: pure A→G spectrum (RNA-editing artefact).
    #[test]
    fn admit_rejects_pure_ag_spectrum() {
        // Host is all A; candidate is all G at those positions → pure A→G.
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'G'), Some(b'G'), Some(b'G')],
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            n_clusters: 3,
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| Some(0.5));
        match got {
            Admission::DnaNeeds(r) => assert!(
                r.reason.contains("A->G"),
                "expected 'A->G' in reason, got: {}",
                r.reason
            ),
            _ => panic!("expected DnaNeeds"),
        }
    }

    /// Gate 5 fires: remap identity ≥ remap_max_identity (0.99 ≥ 0.98).
    #[test]
    fn admit_rejects_high_remap_identity() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C'), Some(b'T'), Some(b'C')],
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            n_clusters: 3,
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| Some(0.99));
        match got {
            Admission::DnaNeeds(r) => assert!(
                r.reason.contains("remap identity"),
                "expected 'remap identity' in reason, got: {}",
                r.reason
            ),
            _ => panic!("expected DnaNeeds"),
        }
    }

    /// Gate 5 fires: no alignment returned (None).
    #[test]
    fn admit_rejects_no_remap_alignment() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C'), Some(b'T'), Some(b'C')],
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            n_clusters: 3,
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| None);
        match got {
            Admission::DnaNeeds(r) => assert!(
                r.reason.contains("no homology"),
                "expected 'no homology' in reason, got: {}",
                r.reason
            ),
            _ => panic!("expected DnaNeeds"),
        }
    }

    /// DnaNeedsRecord carries correct metadata from the candidate.
    #[test]
    fn dna_needs_record_fields_from_candidate() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "chr7".into(),
            start: 500,
            end: 600,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C')],
                read_count: 12,
                identifiable: true,
            },
            psv_pos: vec![502],
            n_clusters: 1, // gate 1 fires
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_| Some(0.5));
        match got {
            Admission::DnaNeeds(r) => {
                assert_eq!(r.chrom, "chr7");
                assert_eq!(r.start, 500);
                assert_eq!(r.end, 600);
                assert_eq!(r.n_clusters, 1);
                assert_eq!(r.read_count, 12);
            }
            _ => panic!("expected DnaNeeds"),
        }
    }

    /// `Admission::Copy` carries the genome remap identity as a second field (Task 1 / copy-graph
    /// v2) so downstream consumers (Task 2: `FamilyAssignment.copy_map_identity`, the `MI:f:` tag)
    /// can thread it through without recomputing minimap2.
    #[test]
    fn admission_copy_carries_identity() {
        let t = DenovoTranscript {
            tid: "t".into(),
            chrom: "c".into(),
            start: 0,
            end: 10,
            n_reads: 5,
            strand: '+',
            introns: vec![],
            seq: b"ACGTACGTAC".to_vec(),
         ..Default::default() };
        let a = Admission::Copy(t.clone(), Some(0.95));
        match a {
            Admission::Copy(_, id) => assert_eq!(id, Some(0.95)),
            _ => panic!("wrong variant"),
        }
    }

    /// Pure parse test (no minimap2 subprocess) of the PAF -> per-read-primary reducer that
    /// `realign_pool_minimap2` uses. Read 0 has two hits: the `tp:A:P` one on contig 2 (mapq 60)
    /// must win over the larger-block secondary on contig 1; read 1 has a single `tp:A:P` hit.
    #[test]
    fn primary_hit_reducer_picks_tp_p_and_reads_mapq() {
        let paf = "0\t100\t0\t100\t+\t2\t120\t0\t100\t95\t100\t60\ttp:A:P\n\
                   0\t100\t0\t100\t+\t1\t120\t0\t100\t80\t100\t0\ttp:A:S\n\
                   1\t100\t0\t100\t+\t0\t120\t0\t100\t90\t100\t0\ttp:A:P\n";
        let hits = super::paf_primary_hits(paf, 2); // 2 reads
        assert_eq!(hits[0], Some((2, 60)));
        assert_eq!(hits[1], Some((0, 0)));
    }

    /// Gate 4 fires: a PSV position (999) is NOT in the host exon frame, so
    /// `collapsed_copy_to_transcript_from_host_seq` returns None → DnaNeeds("...unplaceable...").
    /// The two in-frame columns (102, 104) carry distinct, non-A->G alleles so gates 1-3 pass;
    /// the remap closure should never be reached.
    #[test]
    fn admit_rejects_unplaceable_consensus_as_dna_needs() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                // distinct, non-A->G at the in-frame columns; an allele at the off-frame column.
                allele_vector: vec![Some(b'C'), Some(b'T'), Some(b'C')],
                read_count: 5,
                identifiable: true,
            },
            // 102 and 104 are in the host exon frame; 999 is NOT -> gate 4 None.
            psv_pos: vec![102, 104, 999],
            n_clusters: 3,
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| {
            panic!("remap must not be reached when gate 4 fails")
        });
        match got {
            Admission::DnaNeeds(r) => assert!(
                r.reason.contains("unplaceable"),
                "expected 'unplaceable' in reason, got: {}",
                r.reason
            ),
            _ => panic!("expected DnaNeeds"),
        }
    }
}
