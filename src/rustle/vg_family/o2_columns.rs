//! O2 VG-materialization step 2 (the crux): per-column paralogous-allele extraction
//! from a BAM aligned onto a single family backbone contig.
//!
//! Faithful Rust port of `bench/psv_graph_genomewide.py::column_alleles` -- the
//! pysam-pileup BAM parser that `o2_vg_visualization.materialize_family` (via the
//! genome-wide PSV-graph driver) uses to read, for every reference column of the
//! backbone `copy0`, each aligned copy's / read's base. Those columns are exactly
//! the PSV "bubbles" of `copy_assignment_theory.md`; this function is what feeds
//! [`super::o2_margin_gate`].
//!
//! ## The Python (bench/psv_graph_genomewide.py ~line 45)
//! ```python
//! def column_alleles(bampath, ref):
//!     out = collections.defaultdict(dict)        # ref_pos(int) -> {read_name -> base(char)}
//!     bam = pysam.AlignmentFile(bampath, "rb")
//!     for col in bam.pileup(ref, min_base_quality=0, truncate=False,
//!                           max_depth=300000, flag_filter=0xF04):
//!         for pr in col.pileups:
//!             if pr.query_position is not None:
//!                 out[col.reference_pos][pr.alignment.query_name] = \
//!                     pr.alignment.query_sequence[pr.query_position]
//!     bam.close(); return out
//! ```
//!
//! ## Why a plain per-read CIGAR walk is byte-identical to that pileup
//! With `min_base_quality=0` (no base-quality gate), `max_depth=300000`
//! (effectively uncapped for these single-family BAMs), `truncate=False`
//! (whole contig) and single-end long reads (no mate-overlap clipping, no BAQ
//! effect), the samtools pileup reduces to: for every kept alignment on `ref`,
//! walk its CIGAR against the stored SEQ and, at each reference position covered
//! by an aligned (M/=/X) base, record `out[ref_pos][query_name] = SEQ[query_pos]`.
//! `pr.query_position is None` for D/N (reference-only) columns, so those are
//! skipped -- exactly what advancing the reference cursor without emitting does.
//! This equivalence was verified empirically (pysam vs. this walk agree cell-for-cell)
//! on the real backbone BAMs shipped as fixtures; see `bench/gen_o2_columns_fixture.py`.
//!
//! ## flag_filter 0xF04
//! `0xF04 = 3844 = UNMAP(0x4) | SECONDARY(0x100) | QCFAIL(0x200) | DUP(0x400) |
//! SUPPLEMENTARY(0x800)`. An alignment with ANY of those bits set is EXCLUDED; the
//! kept set is primary + mapped + non-secondary + non-supplementary + non-dup +
//! non-qcfail. (Note the reverse-strand bit 0x10 is NOT in the mask, so reverse
//! alignments DO contribute -- and their SEQ is already stored forward-in-reference,
//! so no reverse-complement is applied: the base recorded is `SEQ[query_pos]` verbatim.)
//!
//! ## CIGAR op handling (BAM op codes)
//! * `M`(0) / `=`(7=SequenceMatch) / `X`(8=SequenceMismatch): aligned -- for each of
//!   the `len` covered reference positions record `SEQ[query_pos]`; advance BOTH cursors.
//! * `D`(2) / `N`(3=Skip/intron): reference-only -- advance the reference cursor, emit
//!   nothing (pysam `query_position is None` there).
//! * `I`(1) / `S`(4=SoftClip): query-only -- advance the query cursor (soft-clipped bases
//!   are present in SEQ and are stepped over; inserted bases likewise).
//! * `H`(5=HardClip) / `P`(6=Pad): advance neither (hard-clipped bases are NOT in SEQ).
//!
//! The output key is the 0-based reference position (`col.reference_pos` is 0-based;
//! `alignment_start` is 1-based so we subtract 1). Each kept alignment contributes at
//! most once per reference position; if the same `query_name` recurs it is overwritten
//! (defaultdict `dict` assignment), matching Python. `BTreeMap` gives a deterministic
//! ascending iteration order.

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::Result;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::RecordBuf;

/// pysam `flag_filter=0xF04`: exclude UNMAP|SECONDARY|QCFAIL|DUP|SUPPLEMENTARY.
pub const FLAG_FILTER: u16 = 0xF04;

/// Per-column paralogous-allele map: `ref_pos (0-based)` -> `read_name` -> `base`.
///
/// Faithful port of `bench/psv_graph_genomewide.py::column_alleles`. Reuses the
/// crate's noodles BAM reader (`crate::bam::open_bam`, the shared non-DEAD helper);
/// no `.bai` index is needed because we walk every record once and filter by the
/// backbone contig `ref_name` (the Python passes a region to `pileup`, but with
/// `truncate=False` over the whole contig the record set is identical).
///
/// `ref_name` must be a contig present in the BAM header (the shipped pipeline only
/// ever passes the backbone name `copy0`). An unknown contig yields an empty map
/// (the Python's `pileup` would raise on an unknown contig, but that never occurs;
/// an in-header contig with no alignments legitimately yields an empty map, matching
/// pysam's empty pileup iterator).
pub fn column_alleles<P: AsRef<Path>>(
    bam_path: P,
    ref_name: &str,
) -> Result<BTreeMap<i64, BTreeMap<String, char>>> {
    let mut reader = crate::bam::open_bam(bam_path, 1)?;
    let header = reader.read_header()?;

    // Resolve the backbone contig name to its numeric reference id. The BAM
    // reference id is the 0-based index into the @SQ order == IndexMap insertion
    // order, so `.position(..)` gives exactly the id `reference_sequence_id()` reports.
    let target_id: Option<usize> = header
        .reference_sequences()
        .iter()
        .position(|(name, _)| format!("{}", name) == ref_name);
    let Some(target_id) = target_id else {
        return Ok(BTreeMap::new());
    };

    let mut out: BTreeMap<i64, BTreeMap<String, char>> = BTreeMap::new();
    let mut record = RecordBuf::default();
    while reader.read_record_buf(&header, &mut record)? > 0 {
        // flag_filter 0xF04 -- ANY excluded bit set => skip (checked before touching
        // the CIGAR/SEQ, so unmapped records with no CIGAR are safely dropped here).
        if record.flags().bits() & FLAG_FILTER != 0 {
            continue;
        }
        // Only the backbone contig's alignments contribute (pysam `pileup(ref, ...)`).
        if record.reference_sequence_id() != Some(target_id) {
            continue;
        }
        let Some(start) = record.alignment_start() else {
            continue;
        };
        // pysam `reference_pos` is 0-based; `alignment_start` is 1-based.
        let ref_start = start.get() as i64 - 1;

        // SEQ (forward-in-reference; noodles decodes the 4-bit codes to ASCII bases
        // via the full 16-symbol table, matching pysam `query_sequence`). A record
        // with no SEQ (`query_sequence is None`) contributes nothing, matching the
        // Python (its `query_sequence[query_position]` would never be reached).
        let seq_bytes: &[u8] = record.sequence().as_ref();
        if seq_bytes.is_empty() {
            continue;
        }
        let name = match record.name() {
            Some(n) => n.to_string(),
            None => continue,
        };

        let mut ref_pos = ref_start;
        let mut qpos: usize = 0;
        for op in record.cigar().as_ref().iter() {
            let len = op.len();
            match op.kind() {
                // Aligned columns: emit SEQ[qpos] at each covered reference position.
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for k in 0..len {
                        let qo = qpos + k;
                        if qo >= seq_bytes.len() {
                            break; // defensive; never trips on a well-formed BAM
                        }
                        let base = seq_bytes[qo] as char;
                        out.entry(ref_pos + k as i64)
                            .or_default()
                            .insert(name.clone(), base);
                    }
                    ref_pos += len as i64;
                    qpos += len;
                }
                // Reference-only ops (query_position is None there -> nothing emitted).
                Kind::Deletion | Kind::Skip => {
                    ref_pos += len as i64;
                }
                // Query-only ops (inserted / soft-clipped bases are present in SEQ).
                Kind::Insertion | Kind::SoftClip => {
                    qpos += len;
                }
                // Neither cursor advances (hard-clipped bases are absent from SEQ).
                Kind::HardClip | Kind::Pad => {}
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use std::path::PathBuf;

    fn testdata_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/rustle/vg_family/testdata")
    }

    fn load_fixture() -> Value {
        let raw = include_str!("testdata/o2_columns_fixture.json");
        serde_json::from_str(raw).expect("parse o2_columns fixture json")
    }

    /// golden: {pos(str) -> {name -> base(str, single char)}} -> BTreeMap<i64, BTreeMap<String,char>>
    fn parse_golden(v: &Value) -> BTreeMap<i64, BTreeMap<String, char>> {
        v.as_object()
            .unwrap()
            .iter()
            .map(|(pos, inner)| {
                let p = pos.parse::<i64>().unwrap();
                let m: BTreeMap<String, char> = inner
                    .as_object()
                    .unwrap()
                    .iter()
                    .map(|(name, base)| {
                        let b = base.as_str().unwrap();
                        assert_eq!(b.chars().count(), 1, "base must be one char, got {b:?}");
                        (name.clone(), b.chars().next().unwrap())
                    })
                    .collect();
                (p, m)
            })
            .collect()
    }

    /// Compare two golden maps and, on the FIRST divergence (ascending pos, then name),
    /// return a human-readable location for the assertion message.
    fn first_mismatch(
        got: &BTreeMap<i64, BTreeMap<String, char>>,
        exp: &BTreeMap<i64, BTreeMap<String, char>>,
    ) -> Option<String> {
        let mut positions: Vec<i64> = got.keys().chain(exp.keys()).copied().collect();
        positions.sort_unstable();
        positions.dedup();
        for p in positions {
            let g = got.get(&p);
            let e = exp.get(&p);
            match (g, e) {
                (None, Some(_)) => return Some(format!("pos {p}: present in pysam golden, MISSING in rust")),
                (Some(_), None) => return Some(format!("pos {p}: present in rust, MISSING in pysam golden")),
                (Some(gm), Some(em)) => {
                    let mut names: Vec<&String> = gm.keys().chain(em.keys()).collect();
                    names.sort();
                    names.dedup();
                    for nm in names {
                        match (gm.get(nm), em.get(nm)) {
                            (Some(gb), Some(eb)) if gb == eb => {}
                            (gb, eb) => {
                                return Some(format!(
                                    "pos {p} name {nm}: rust={:?} pysam={:?}",
                                    gb, eb
                                ));
                            }
                        }
                    }
                }
                (None, None) => {}
            }
        }
        None
    }

    #[test]
    fn byte_parity_all_bams() {
        let fx = load_fixture();
        let bams = fx["bams"].as_array().expect("bams array");
        assert!(!bams.is_empty(), "no bam fixtures");

        let mut n_real = 0usize;
        let mut n_adv = 0usize;
        let mut total_cells = 0usize;

        for entry in bams {
            let file = entry["file"].as_str().unwrap();
            let ref_name = entry["ref"].as_str().unwrap();
            let kind = entry["kind"].as_str().unwrap();
            match kind {
                "real" => n_real += 1,
                _ => n_adv += 1,
            }
            let exp = parse_golden(&entry["golden"]);
            total_cells += exp.values().map(|m| m.len()).sum::<usize>();

            let path = testdata_dir().join(file);
            let got = column_alleles(&path, ref_name)
                .unwrap_or_else(|e| panic!("column_alleles({file}, {ref_name}) failed: {e}"));

            if let Some(loc) = first_mismatch(&got, &exp) {
                panic!(
                    "byte-parity MISMATCH in {file} (ref {ref_name}, kind {kind}): {loc}\n  \
                     rust cols={}, pysam cols={}",
                    got.len(),
                    exp.len()
                );
            }
            // exact whole-map equality (belt and suspenders beyond first_mismatch)
            assert_eq!(
                got, exp,
                "whole-map inequality in {file} (ref {ref_name}) not localized by first_mismatch"
            );
        }

        assert!(n_real >= 3, "expected >=3 real BAM fixtures, got {n_real}");
        assert!(n_adv >= 1, "expected >=1 adversarial BAM fixture, got {n_adv}");
        eprintln!(
            "o2_columns byte-parity: {} BAMs ({} real, {} adversarial), {} column-cells all match pysam",
            bams.len(),
            n_real,
            n_adv,
            total_cells
        );
    }

    /// Focused assertions on the adversarial BAM: the excluded-flag reads must NOT
    /// appear anywhere, reverse-strand reads contribute their stored (forward) base
    /// with no reverse-complement, and an in-header contig with no reads is empty.
    #[test]
    fn adversarial_flag_and_strand_semantics() {
        let fx = load_fixture();
        let bams = fx["bams"].as_array().unwrap();
        let adv = bams
            .iter()
            .find(|b| b["kind"].as_str().unwrap() == "adversarial" && b["ref"].as_str().unwrap() != "empty")
            .expect("adversarial (non-empty ref) entry");
        let file = adv["file"].as_str().unwrap();
        let ref_name = adv["ref"].as_str().unwrap();
        let got = column_alleles(testdata_dir().join(file), ref_name).unwrap();

        // No excluded-flag read name may appear at any column.
        let excluded = fx["excluded_read_names"]
            .as_array()
            .expect("excluded_read_names")
            .iter()
            .map(|v| v.as_str().unwrap().to_string())
            .collect::<Vec<_>>();
        assert!(!excluded.is_empty(), "fixture must list excluded read names");
        for (p, m) in &got {
            for name in m.keys() {
                assert!(
                    !excluded.contains(name),
                    "excluded-flag read {name} leaked into column {p}"
                );
            }
        }

        // The empty in-header contig yields an empty map.
        let empty = bams
            .iter()
            .find(|b| b["ref"].as_str().unwrap() == "empty")
            .expect("empty-contig entry");
        let empty_got =
            column_alleles(testdata_dir().join(empty["file"].as_str().unwrap()), "empty").unwrap();
        assert!(empty_got.is_empty(), "empty contig must yield empty map");
    }
}
