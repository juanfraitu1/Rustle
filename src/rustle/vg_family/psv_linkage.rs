//! Layer-2 "C": within-molecule PSV->junction linkage for paralog isoform recovery.
//! Task 1: PSV-column extraction (positions + per-copy alleles + per-copy genomic coords).
//!
//! A PSV (paralog-sequence variant) column is one base position, inside a
//! diagnostic exon node, where >=2 copies of the family carry distinct alleles.
//! It is the genotyping reference for later tasks: a read covering this column
//! tells us which copy it came from (its base matches that copy's allele).
//!
//! WHY this is its own primitive (pure, not yet wired): the column carries the
//! genomic coordinate IN EACH COPY's OWN FRAME (from `per_copy_spans`, never the
//! cross-copy union `ExonClass.span`), so the genotyping pass in a later task can
//! look up "the base this read has at copy X's position 502" directly. Computing
//! that frame mapping once, deterministically, keeps the downstream pass cheap and
//! the whole channel reproducible (DetHash + total-order sorts everywhere).

use crate::vg_family::family_graph::FamilyGraph;

/// One copy's allele at a PSV column: the base and its genomic position IN THAT
/// COPY's frame (the copy's exon start + the within-exon offset).
#[derive(Debug, Clone, PartialEq)]
pub struct PsvCopyAllele {
    pub copy_id: usize,
    pub genomic_pos: u64,
    pub allele: u8,
}

/// One paralog-distinguishing column: a position (within a diagnostic exon node)
/// where >= 2 copies carry distinct bases. `per_copy` is sorted by copy_id so the
/// column is canonical regardless of the order copies were pushed into the node.
#[derive(Debug, Clone, PartialEq)]
pub struct PsvColumn {
    pub node_idx: usize,
    pub per_copy: Vec<PsvCopyAllele>,
}

/// Extract the family's PSV columns. For each diagnostic node (a node whose
/// `per_copy_sequences` carry >= 2 DISTINCT sequences), for each within-exon
/// offset where the present copies' bases differ, emit a `PsvColumn`: each present
/// copy's base + the genomic coordinate `per_copy_spans[copy].0 + offset` (the
/// copy's OWN frame, never `ExonClass.span`, which is a cross-copy union and can
/// span the whole chromosome).
///
/// v1 limitation (documented): only handle nodes where ALL contributing copies'
/// per_copy_sequences are EQUAL LENGTH (no indel re-alignment) — skip a node whose
/// copies' sequences differ in length. With equal length, offset i is the same
/// genomic position relative to each copy's start, so the per-copy frame mapping is
/// exact; an indel would shift offsets between copies and require alignment, which
/// is deferred. A copy that has a sequence but no `per_copy_spans` entry is skipped
/// (we cannot place its allele in genomic coordinates without its frame).
///
/// Deterministic: the returned Vec is sorted by (node_idx, the column's first
/// genomic_pos); each column's `per_copy` is sorted by copy_id.
pub fn psv_columns_for_family(fg: &FamilyGraph) -> Vec<PsvColumn> {
    let mut columns: Vec<PsvColumn> = Vec::new();

    for node in &fg.nodes {
        // A node's contributing copies = its per_copy_sequences entries. We must be
        // able to place each in genomic coordinates, so a copy with a sequence but
        // no per_copy_spans entry is dropped (we cannot know its frame).
        let mut contributors: Vec<(usize, &[u8], u64)> = Vec::new();
        for (cid, seq) in &node.per_copy_sequences {
            if let Some((_, (start, _))) = node.per_copy_spans.iter().find(|(c, _)| c == cid) {
                contributors.push((*cid, seq.as_slice(), *start));
            }
        }

        // Diagnostic = >= 2 contributors carrying >= 2 DISTINCT sequences.
        if contributors.len() < 2 {
            continue;
        }
        let distinct = {
            let mut s: crate::types::DetHashSet<&[u8]> = crate::types::DetHashSet::default();
            for (_, seq, _) in &contributors {
                s.insert(seq);
            }
            s.len()
        };
        if distinct < 2 {
            continue;
        }

        // v1: only equal-length contributors (no indel re-alignment). If lengths
        // differ, offset i is NOT the same genomic position across copies, so we
        // skip the whole node rather than misalign columns.
        let len0 = contributors[0].1.len();
        if contributors.iter().any(|(_, seq, _)| seq.len() != len0) {
            continue;
        }

        // For each within-exon offset, if >= 2 distinct bases are present among
        // contributors, emit a column with EVERY present copy's base placed in its
        // OWN frame (start + offset). per_copy sorted by copy_id for canonical form.
        for offset in 0..len0 {
            let mut distinct_bases: crate::types::DetHashSet<u8> = crate::types::DetHashSet::default();
            for (_, seq, _) in &contributors {
                distinct_bases.insert(seq[offset]);
            }
            if distinct_bases.len() < 2 {
                continue;
            }
            let mut per_copy: Vec<PsvCopyAllele> = contributors
                .iter()
                .map(|(cid, seq, start)| PsvCopyAllele {
                    copy_id: *cid,
                    genomic_pos: start + offset as u64,
                    allele: seq[offset],
                })
                .collect();
            per_copy.sort_by_key(|p| p.copy_id);
            columns.push(PsvColumn { node_idx: node.idx.0, per_copy });
        }
    }

    // Total order: (node_idx, the column's first genomic_pos). per_copy is already
    // sorted by copy_id, so "first genomic_pos" is the lowest-copy_id allele's pos
    // — a stable, deterministic key.
    columns.sort_by(|a, b| {
        let a_first = a.per_copy.first().map(|p| p.genomic_pos).unwrap_or(0);
        let b_first = b.per_copy.first().map(|p| p.genomic_pos).unwrap_or(0);
        a.node_idx.cmp(&b.node_idx).then(a_first.cmp(&b_first))
    });

    columns
}

/// (E) family identifiability gate. A family is identifiable enough to attempt
/// PSV-linkage iff it has >= `min_psv_columns` distinguishable PSV columns on its
/// exons. The `error_rate` argument is RESERVED for the coverage-weighted
/// refinement (P(a read spans >= N PSVs) given read length and error); v1 is
/// count-based, since requiring `min_psv_columns` distinct columns already makes
/// `min_psv_columns` independent error-agreements improbable at typical long-read
/// error rates. Returns false for families below the floor (e.g. DAZ-like: long
/// identical cores with ~1 PSV), so they are skipped wholesale — the primary
/// phantom defense.
///
/// Naming: `min_psv_columns` is the FAMILY's distinguishing-column count, distinct
/// from the per-READ `min_psv` (the agreeing-PSV count a single read needs, used by
/// the later read-assignment task).
pub fn family_identifiability(fg: &FamilyGraph, error_rate: f64, min_psv_columns: usize) -> bool {
    let _ = error_rate; // reserved for the coverage-weighted refinement (documented v1 limitation)
    psv_columns_for_family(fg).len() >= min_psv_columns
}

/// One read's PSV evidence + its junction (intron) chain, aggregated across ALL of
/// the read's alignments overlapping the family (its primary at the well copy AND
/// any secondaries cross-mapped to siblings — they carry the same molecule's bases).
///
/// `psv_votes`: `(copy_id, #PSVs whose allele this read's base matched)`. WHY this
/// is the load-bearing signal: a read physically aligned at copy X's locus reads
/// copy X's genomic coordinates, but the read's *actual base* at a PSV column tells
/// us its true origin — `base == allele_X` ⇒ native copy-X read; `base == allele_Y`
/// (a sibling's allele) ⇒ a copy-Y read CROSS-MAPPED to X. So at each PSV the read
/// covers we vote for whichever copy's allele the base matches. Aggregating these
/// votes across the PSVs the read spans is its copy genotype — the same allele
/// signal that beats the homology-shadow phantom Part A cannot.
///
/// `junctions`: the read's introns (`exons.windows(2)` → `(prev_end, next_start)`),
/// the splice chain a later task links to the assigned copy.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadGenotype {
    pub read_name_hash: u64,
    pub psv_votes: Vec<(usize, usize)>,
    pub junctions: Vec<(u64, u64)>,
}

/// Confidently assign a read to one copy iff that copy has >= `min_psv` supporting
/// PSVs AND strictly dominates every sibling by >= `margin`
/// (`votes[X] - max_other >= margin`). Otherwise `None` — the read is AMBIGUOUS and
/// is never guessed, because a wrong assignment silently fabricates a phantom
/// per-copy isoform downstream. Deterministic: ties and sub-threshold votes → None;
/// when several copies tie for the top vote no copy can beat the runner-up by the
/// margin, so the result is `None`.
///
/// `min_psv` is the per-READ floor (env `RUSTLE_VG_LAYER2_PSV_MIN`, default 3) — the
/// ~3-PSV/0.997-identity sequencing-error floor; it is DISTINCT from the family
/// gate's `min_psv_columns`.
pub fn assign_read_to_copy(g: &ReadGenotype, min_psv: usize, margin: usize) -> Option<usize> {
    // Find the top-voted copy and the best vote among all OTHER copies. We scan in a
    // single deterministic pass; ties for the top are detected via the runner-up.
    let mut best_copy: Option<usize> = None;
    let mut best_votes: usize = 0;
    let mut second_votes: usize = 0;
    for &(copy_id, votes) in &g.psv_votes {
        if votes > best_votes {
            second_votes = best_votes;
            best_votes = votes;
            best_copy = Some(copy_id);
        } else if votes > second_votes {
            second_votes = votes;
        }
    }
    let copy = best_copy?;
    // Beat the error floor AND strictly dominate every sibling by the margin. The
    // margin check folds in tie rejection: a tie makes `best - second == 0 < margin`.
    if best_votes >= min_psv && best_votes.saturating_sub(second_votes) >= margin {
        Some(copy)
    } else {
        None
    }
}

/// Targeted second BAM pass that GENOTYPES each read overlapping the family copies'
/// genomic spans, then aggregates per molecule. This is "Approach A" of the design:
/// re-read the BAM only over `family_loci`, reuse `bam.rs`'s per-position
/// mismatch-vs-FASTA reconstruction, and at every PSV column a read covers, vote for
/// the copy whose allele the read's base matches (see [`ReadGenotype`] for the
/// allele-voting rationale).
///
/// All alignments overlapping a family locus are read — primaries AND secondaries
/// AND supplementaries: a read's primary at the well copy and its secondary
/// cross-mapped to a sibling both carry the SAME molecule's bases, so both should
/// vote; votes are aggregated per `read_name_hash`. (No `RunConfig` is threaded here,
/// so no long-read-class filter is applied — the work is already bounded tightly by
/// `family_loci` × `psv_columns`, and the per-read PSV floor in `assign_read_to_copy`
/// is what actually gates short/uninformative reads out.)
///
/// Coverage rule: a PSV column's per-copy `genomic_pos` is "covered" by an alignment
/// iff it lies within `[ref_start, ref_end)` AND inside one of the read's aligned
/// exons (NOT an intron gap) — an intronic position carries no base. A read aligned
/// at copy X's locus covers only copy X's frame position for that column, so we read
/// the read's base once at the covered position and compare it against EVERY copy's
/// allele in the column, voting for the one it matches (no match → no vote).
///
/// Deterministic: records are processed in BAM order; aggregation lives in
/// `DetHashMap`/`DetHashSet`; the returned Vec is sorted by `read_name_hash`, each
/// `psv_votes` sorted by `copy_id`, and each `junctions` sorted.
///
/// Known v1 limitation (inherited from the eqx reconstruction): a read with a
/// non-canonical base (e.g. `N`) at a PSV position is not recorded as a mismatch by
/// `extract_mismatches_vs_fasta`, so it falls back to the reference base and casts a
/// native-copy vote. This is rare on HiFi reads and only ever inflates the NATIVE
/// vote (never fabricates a cross-copy vote), so it is conservative w.r.t. phantoms;
/// the per-read `min_psv` floor + strict-dominance margin absorb the occasional
/// spurious native vote. Re-deriving the true base would require duplicating the
/// CIGAR walk, which we deliberately avoid.
pub fn genotype_family_reads(
    bam_path: &std::path::Path,
    family_loci: &[(String, u64, u64)],
    psv_columns: &[PsvColumn],
    genome: &crate::genome::GenomeIndex,
) -> anyhow::Result<Vec<ReadGenotype>> {
    use crate::types::{DetHashMap, DetHashSet};

    // votes: read_name_hash -> (copy_id -> #supporting PSVs); junctions per read.
    let mut votes_by_read: DetHashMap<u64, DetHashMap<usize, usize>> = DetHashMap::default();
    let mut junctions_by_read: DetHashMap<u64, DetHashSet<(u64, u64)>> = DetHashMap::default();
    // Track read order of first appearance is unnecessary — we sort by hash at the
    // end — but we must register every read seen even if it casts no votes.
    let mut seen_reads: DetHashSet<u64> = DetHashSet::default();

    let mut reader = crate::bam::open_bam(bam_path, 1)?;
    let header = reader.read_header()?;
    let mut record = noodles_sam::alignment::RecordBuf::default();

    while reader.read_record_buf(&header, &mut record)? > 0 {
        if record.flags().is_unmapped() {
            continue;
        }
        // Chromosome name for this alignment (needed for locus overlap + base lookup).
        let chrom = match record
            .reference_sequence(&header)
            .transpose()
            .ok()
            .flatten()
            .map(|(n, _)| String::from_utf8_lossy(n).into_owned())
        {
            Some(n) => n,
            None => continue,
        };

        // Parse once via the same parser Layer 1 uses, so exons/junctions/hashing
        // never diverge from the rest of the pipeline. Primaries, secondaries, and
        // supplementaries are ALL kept (each carries the molecule's bases).
        let Some(read) = crate::bam::record_to_bundle_read(&record) else {
            continue;
        };
        let ref_start = read.ref_start;
        let ref_end = read.ref_end; // exclusive (last exon end)

        // Only consider alignments that overlap a family locus on the same chrom.
        let overlaps_family = family_loci.iter().any(|(c, s, e)| {
            c == &chrom && ref_start < *e && *s < ref_end
        });
        if !overlaps_family {
            continue;
        }

        // Register the read so it appears in the output even with zero votes (its
        // junctions still matter for later linkage of confidently-assigned reads).
        seen_reads.insert(read.read_name_hash);

        // The read's intron chain = gaps between consecutive aligned exons.
        if read.exons.len() >= 2 {
            let jset = junctions_by_read.entry(read.read_name_hash).or_default();
            for w in read.exons.windows(2) {
                jset.insert((w[0].1, w[1].0));
            }
        }

        // Per-position mismatch set vs the reference FASTA (the `--eqx`-equivalent).
        // A read's base at ref position P = the FASTA base at P UNLESS P is in this
        // set, in which case it is the recorded alt base.
        let mismatches = crate::bam::extract_mismatches_vs_fasta(&record, ref_start, &chrom, genome);
        let mismatch_at = |pos: u64| -> Option<u8> {
            mismatches.iter().find(|(p, _)| *p == pos).map(|(_, b)| *b)
        };

        // Vote at every PSV column this alignment covers.
        for col in psv_columns {
            // A read aligned at one copy's locus covers at most ONE per-copy frame
            // position for this column. Find the covered position (inside an aligned
            // exon), read the base there, then vote for whichever copy's allele it
            // matches.
            for pc in &col.per_copy {
                let pos = pc.genomic_pos;
                if pos < ref_start || pos >= ref_end {
                    continue;
                }
                // Must fall inside an aligned exon, not an intron gap.
                let in_exon = read.exons.iter().any(|(s, e)| pos >= *s && pos < *e);
                if !in_exon {
                    continue;
                }
                // The read's base at this reference position.
                let base = match mismatch_at(pos) {
                    Some(b) => b,
                    None => match genome
                        .fetch_sequence(&chrom, pos, pos + 1)
                        .and_then(|v| v.first().copied())
                    {
                        Some(b) => b.to_ascii_uppercase(),
                        None => continue, // no reference base -> cannot genotype here
                    },
                };
                // Compare the base against EVERY copy's allele in this column and
                // vote for the copy it matches (the cross-mapping detector). At most
                // one position per column is covered by this alignment, so we break
                // after handling it to avoid double-counting if frames ever collide.
                let copy_votes = votes_by_read.entry(read.read_name_hash).or_default();
                for cand in &col.per_copy {
                    if cand.allele == base {
                        *copy_votes.entry(cand.copy_id).or_insert(0) += 1;
                    }
                }
                break;
            }
        }
    }

    // Emit one ReadGenotype per read seen, sorted by read_name_hash; psv_votes
    // sorted by copy_id; junctions sorted. Determinism is load-bearing.
    let mut out: Vec<ReadGenotype> = seen_reads
        .into_iter()
        .map(|h| {
            let mut psv_votes: Vec<(usize, usize)> = votes_by_read
                .get(&h)
                .map(|m| m.iter().map(|(&c, &n)| (c, n)).collect())
                .unwrap_or_default();
            psv_votes.sort_by_key(|(c, _)| *c);
            let mut junctions: Vec<(u64, u64)> = junctions_by_read
                .get(&h)
                .map(|s| s.iter().copied().collect())
                .unwrap_or_default();
            junctions.sort();
            ReadGenotype { read_name_hash: h, psv_votes, junctions }
        })
        .collect();
    out.sort_by_key(|g| g.read_name_hash);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vg_family::family_graph::{ExonClass, FamilyGraph, JunctionEdge, NodeIdx};

    /// Build a one-node family graph by hand (modeled on layer2.rs's
    /// `build_psv_three_copy_fg`). `seqs` = (copy_id, exon sequence, exon span).
    /// The node's union `span` is set to a deliberately WRONG, whole-chromosome
    /// value so any test that accidentally reads `span` instead of `per_copy_spans`
    /// produces obviously-wrong coordinates.
    fn one_node_fg(seqs: &[(usize, &[u8], (u64, u64))]) -> FamilyGraph {
        let node = ExonClass {
            idx: NodeIdx(0),
            chrom: "chrT".into(),
            // Cross-copy union landmine: NEVER the coordinate source.
            span: (0, 1_000_000),
            strand: '+',
            per_copy_sequences: seqs.iter().map(|(c, s, _)| (*c, s.to_vec())).collect(),
            per_copy_spans: seqs.iter().map(|(c, _, sp)| (*c, *sp)).collect(),
            copy_specific: seqs.len() == 1,
            per_copy_cov: Vec::new(),
        };
        FamilyGraph { family_id: 0, nodes: vec![node], edges: Vec::<JunctionEdge>::new() }
    }

    /// Two-copy diagnostic node: copy0 b"ACGTAC"@(100,106), copy1 b"ACCTAC"@(500,506).
    /// They differ ONLY at offset 2 (G vs C). The single PSV column must map that
    /// offset into each copy's OWN frame: copy0 -> 102, copy1 -> 502.
    fn two_copy_one_psv_fg() -> FamilyGraph {
        one_node_fg(&[
            (0, b"ACGTAC", (100, 106)),
            (1, b"ACCTAC", (500, 506)),
        ])
    }

    /// RABL2-like (identifiable): two copies of equal length differing at >= 3
    /// offsets, so `psv_columns_for_family` yields >= 3 distinguishing columns.
    /// copy0 b"ACGTACGT", copy1 b"AGGAAGGT" differ at offsets 1, 3, 4 -> 3 columns.
    fn rabl2_like_fg() -> FamilyGraph {
        one_node_fg(&[
            (0, b"ACGTACGT", (100, 108)),
            (1, b"AGGAAGGT", (500, 508)),
        ])
    }

    /// DAZ-like (unidentifiable): long near-identical cores differing at exactly
    /// ONE offset -> 1 column, below the >= 3 floor.
    fn daz_like_fg() -> FamilyGraph {
        one_node_fg(&[
            (0, b"ACGTACGT", (100, 108)),
            (1, b"ACGTACGA", (500, 508)),
        ])
    }

    /// Exactly-three-columns fixture for the `>=` boundary test: copy0 b"AAAAAA"
    /// vs copy1 b"ACCCAA" differ at offsets 1, 2, 3 -> 3 columns.
    fn three_psv_fg() -> FamilyGraph {
        one_node_fg(&[
            (0, b"AAAAAA", (100, 106)),
            (1, b"ACCCAA", (500, 506)),
        ])
    }

    #[test]
    fn identifiable_family_passes() {
        // >= 3 distinguishable PSV columns on exons -> family is identifiable.
        assert_eq!(psv_columns_for_family(&rabl2_like_fg()).len(), 3);
        assert!(family_identifiability(&rabl2_like_fg(), 0.001, 3));
    }

    #[test]
    fn unidentifiable_family_fails() {
        // 1 PSV column (DAZ-like long identical cores) -> below the floor -> skip.
        assert_eq!(psv_columns_for_family(&daz_like_fg()).len(), 1);
        assert!(!family_identifiability(&daz_like_fg(), 0.001, 3));
    }

    #[test]
    fn identifiability_threshold_boundary() {
        // Exactly 3 columns: passes at min_psv_columns=3 (>=), fails at 4. Locks
        // the inclusive boundary so a future off-by-one cannot slip the gate.
        let fg = three_psv_fg();
        assert_eq!(psv_columns_for_family(&fg).len(), 3);
        assert!(family_identifiability(&fg, 0.001, 3));
        assert!(!family_identifiability(&fg, 0.001, 4));
    }

    #[test]
    fn psv_columns_map_offset_to_each_copy_frame() {
        let fg = two_copy_one_psv_fg();
        let cols = psv_columns_for_family(&fg);
        assert_eq!(cols.len(), 1, "exactly one differing column (offset 2): {cols:?}");
        let c = &cols[0];
        assert_eq!(c.node_idx, 0);
        assert!(
            c.per_copy.iter().any(|p| p.copy_id == 0 && p.genomic_pos == 102 && p.allele == b'G'),
            "copy0 allele G at its own frame 100+2=102: {c:?}"
        );
        assert!(
            c.per_copy.iter().any(|p| p.copy_id == 1 && p.genomic_pos == 502 && p.allele == b'C'),
            "copy1 allele C at its own frame 500+2=502: {c:?}"
        );
    }

    #[test]
    fn psv_columns_skips_non_diagnostic_node() {
        // Both copies carry the SAME sequence -> not diagnostic -> no columns.
        let fg = one_node_fg(&[
            (0, b"ACGTAC", (100, 106)),
            (1, b"ACGTAC", (500, 506)),
        ]);
        let cols = psv_columns_for_family(&fg);
        assert!(cols.is_empty(), "identical copies are not diagnostic: {cols:?}");
    }

    #[test]
    fn psv_columns_skips_unequal_length() {
        // Copies differ in length -> v1 cannot column-align without indel handling.
        // The sequences are also clearly distinct (so the node IS diagnostic), but
        // the length mismatch must still suppress all columns from this node.
        let fg = one_node_fg(&[
            (0, b"ACGTAC", (100, 106)),
            (1, b"ACGTACG", (500, 507)),
        ]);
        let cols = psv_columns_for_family(&fg);
        assert!(cols.is_empty(), "v1 skips unequal-length nodes (no indel align): {cols:?}");
    }

    #[test]
    fn psv_columns_deterministic() {
        // Same family, copies pushed in different order -> identical output.
        let fg_a = one_node_fg(&[
            (0, b"ACGTAC", (100, 106)),
            (1, b"ACCTAC", (500, 506)),
        ]);
        let fg_b = one_node_fg(&[
            (1, b"ACCTAC", (500, 506)),
            (0, b"ACGTAC", (100, 106)),
        ]);
        assert_eq!(
            psv_columns_for_family(&fg_a),
            psv_columns_for_family(&fg_b),
            "column extraction must be order-independent (per_copy sorted by copy_id)"
        );
    }

    #[test]
    fn psv_columns_multiple_offsets_sorted() {
        // Two differing offsets (1 and 4) in one node -> two columns, sorted by
        // genomic_pos. Confirms the (node_idx, first genomic_pos) total order.
        let fg = one_node_fg(&[
            (0, b"AAAAAA", (100, 106)),
            (1, b"ACAACA", (500, 506)),
        ]);
        let cols = psv_columns_for_family(&fg);
        assert_eq!(cols.len(), 2, "offsets 1 and 4 differ: {cols:?}");
        // First column = lower copy0 genomic_pos (101 < 104).
        assert_eq!(cols[0].per_copy.iter().find(|p| p.copy_id == 0).unwrap().genomic_pos, 101);
        assert_eq!(cols[1].per_copy.iter().find(|p| p.copy_id == 0).unwrap().genomic_pos, 104);
    }

    // ── Task 3: per-read copy assignment (pure logic, no BAM) ──────────────
    //
    // assign_read_to_copy is the decision a single read makes about its copy of
    // origin from its aggregated PSV votes. The contract is deliberately strict:
    // it assigns ONLY when the evidence beats the sequencing-error floor (>= min_psv
    // agreeing PSVs) AND one copy strictly dominates every sibling (>= margin). A
    // tie or sub-threshold read is left ambiguous — NEVER guessed — because a wrong
    // assignment silently fabricates a phantom isoform downstream.

    #[test]
    fn read_with_sibling_alleles_assigned_to_sibling() {
        // 3 PSVs match copy1's allele, 0 match copy0 -> copy1 wins (the cross-mapped
        // read case: aligned at copy0's locus but its bases say copy1).
        let g = ReadGenotype {
            read_name_hash: 7,
            psv_votes: vec![(0, 0), (1, 3)],
            junctions: vec![(160, 300)],
        };
        assert_eq!(assign_read_to_copy(&g, 3, 1), Some(1));
    }

    #[test]
    fn ambiguous_tie_unassigned() {
        // Even-split votes, and neither reaches min_psv -> ambiguous -> None.
        let g = ReadGenotype {
            read_name_hash: 8,
            psv_votes: vec![(0, 2), (1, 2)],
            junctions: vec![],
        };
        assert_eq!(assign_read_to_copy(&g, 3, 1), None);
    }

    #[test]
    fn sub_threshold_unassigned() {
        // Single copy with 2 votes but min_psv=3 -> below the error floor -> None.
        let g = ReadGenotype {
            read_name_hash: 9,
            psv_votes: vec![(0, 2)],
            junctions: vec![],
        };
        assert_eq!(assign_read_to_copy(&g, 3, 1), None);
    }

    #[test]
    fn margin_enforced() {
        // Both copies clear min_psv but are tied -> no strict dominance -> None.
        let tied = ReadGenotype {
            read_name_hash: 10,
            psv_votes: vec![(0, 3), (1, 3)],
            junctions: vec![],
        };
        assert_eq!(assign_read_to_copy(&tied, 3, 1), None);
        // One copy clears min_psv and beats the runner-up by the margin -> assigned.
        let dominant = ReadGenotype {
            read_name_hash: 11,
            psv_votes: vec![(0, 4), (1, 3)],
            junctions: vec![],
        };
        assert_eq!(assign_read_to_copy(&dominant, 3, 1), Some(0));
    }

    #[test]
    fn genotype_family_reads_smoke_constructs() {
        // Smoke-level exercise of the pass's signature + empty-input path: with no
        // PSV columns and no loci, an absent BAM path must not panic before it even
        // tries to open the file is not asserted here (it would error); instead we
        // assert the pure pieces compose. The end-to-end BAM I/O assertion (a read
        // carrying sibling alleles genotyping to the sibling) is DEFERRED to Task 7's
        // `sim_psvlink` fixture, per the plan (Task 3 Step 5 permits deferral).
        //
        // This keeps the function compiled + type-checked here and verifies the
        // ReadGenotype it returns flows into assign_read_to_copy.
        let g = ReadGenotype {
            read_name_hash: 42,
            psv_votes: vec![(0, 3), (1, 1)],
            junctions: vec![(100, 200), (250, 400)],
        };
        assert_eq!(assign_read_to_copy(&g, 3, 1), Some(0));
        // Function exists with the documented signature (compile-time check).
        let _f: fn(
            &std::path::Path,
            &[(String, u64, u64)],
            &[PsvColumn],
            &crate::genome::GenomeIndex,
        ) -> anyhow::Result<Vec<ReadGenotype>> = genotype_family_reads;
    }

    #[test]
    fn isoform_source_psvlinked_orders_after_native() {
        // source_rank parity check: PsvLinked must rank strictly after Native and
        // Transferred, and order among PsvLinked by copy_id. Mirrors the closure in
        // emit_family_isoforms so the layer2 regression anchor holds.
        use crate::vg_family::layer2::IsoformSource;
        let rank = |s: &IsoformSource| -> (usize, usize) {
            match s {
                IsoformSource::Native => (0, 0),
                IsoformSource::Transferred { donor_copy } => (1, *donor_copy),
                IsoformSource::PsvLinked { copy_id } => (2, *copy_id),
            }
        };
        assert!(rank(&IsoformSource::Native) < rank(&IsoformSource::PsvLinked { copy_id: 0 }));
        assert!(
            rank(&IsoformSource::Transferred { donor_copy: 9 })
                < rank(&IsoformSource::PsvLinked { copy_id: 0 })
        );
        assert!(
            rank(&IsoformSource::PsvLinked { copy_id: 0 })
                < rank(&IsoformSource::PsvLinked { copy_id: 1 })
        );
    }
}
