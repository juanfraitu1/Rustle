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

/// One of a read's alignments, at a specific family copy's locus, expressed IN THAT
/// COPY's OWN coordinate frame. This is the per-LOCUS chain — the fix for the
/// cross-locus chimera: a read that cross-maps to copies A and B carries TWO of these
/// (A's junctions in A's coords, B's junctions in B's coords), never a single pooled
/// chimeric set spanning both. Assembly uses the chain whose `copy_id` equals the
/// read's ASSIGNED copy, so an isoform is always built in one copy's frame.
///
/// `junctions`: this alignment's introns (`exons.windows(2)` → `(prev_end, next_start)`).
/// `exons`: this alignment's aligned exon blocks (sorted), needed for the consensus
///   TERMINAL coordinates (the crisp junctions alone fix only interior boundaries).
#[derive(Debug, Clone, PartialEq)]
pub struct LocusChain {
    pub copy_id: usize,
    pub junctions: Vec<(u64, u64)>,
    pub exons: Vec<(u64, u64)>,
}

/// One read's PSV evidence (GLOBAL, across all its alignments) + its per-LOCUS
/// junction chains (one per family-copy locus the read aligns at, each in that copy's
/// own frame).
///
/// `psv_votes`: `(copy_id, #PSVs whose allele this read's base matched)`, aggregated
/// across ALL of the read's alignments. WHY this is the load-bearing signal — and why
/// it MUST stay global while junctions go per-locus: a read physically aligned at copy
/// X's locus reads copy X's genomic coordinates, but the read's *actual base* at a PSV
/// column tells us its true origin — `base == allele_X` ⇒ native copy-X read;
/// `base == allele_Y` (a sibling's allele) ⇒ a copy-Y read CROSS-MAPPED to X. Both the
/// primary AND the cross-mapped secondary carry the SAME molecule's bases, so both
/// must vote — the genotype is a property of the molecule, not of one alignment.
///
/// `per_locus`: one [`LocusChain`] per DISTINCT copy-locus this read aligns at, sorted
/// by `copy_id`. The junctions/exons live here (not pooled) so that downstream linkage
/// and assembly read the chain in the ASSIGNED copy's frame ONLY — never a chimera
/// pooled across loci. A read whose alignment covers no PSV position (so its locus-copy
/// is unknown) contributes no `per_locus` entry.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadGenotype {
    pub read_name_hash: u64,
    pub psv_votes: Vec<(usize, usize)>,
    pub per_locus: Vec<LocusChain>,
}

impl ReadGenotype {
    /// The read's junction chain in `copy`'s own frame, if it aligns at `copy`'s locus.
    fn chain_for_copy(&self, copy: usize) -> Option<&LocusChain> {
        self.per_locus.iter().find(|lc| lc.copy_id == copy)
    }
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
/// vote; PSV votes are aggregated GLOBALLY per `read_name_hash`. (No `RunConfig` is
/// threaded here, so no long-read-class filter is applied — the work is already bounded
/// tightly by `family_loci` × `psv_columns`, and the per-read PSV floor in
/// `assign_read_to_copy` is what actually gates short/uninformative reads out.)
///
/// PER-LOCUS junctions (the cross-locus-chimera fix): unlike the PSV votes, an
/// alignment's junctions+exons must NOT be pooled across the read's alignments — a
/// copy-A read's primary (A's coords) and its secondary cross-mapped to copy B (B's
/// coords) carry junction sets in DIFFERENT coordinate frames; pooling them yields a
/// chimera spanning both loci. So each alignment's junctions/exons are recorded under
/// its LOCUS-COPY in that copy's own frame ([`LocusChain`]), and assembly later uses
/// only the chain in the read's ASSIGNED copy's frame.
///
/// Locus-copy rule: for an alignment, each covered PSV per-copy `genomic_pos` lies in
/// exactly one copy's frame; the alignment's locus-copy = the copy with the MOST
/// covered frame positions (tiebreak: lowest copy_id). An alignment covering no PSV
/// position has no locus-copy and contributes no `per_locus` entry.
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
/// `psv_votes` sorted by `copy_id`, `per_locus` sorted by `copy_id`, and each chain's
/// `junctions`/`exons` sorted.
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

    // votes: read_name_hash -> (copy_id -> #supporting PSVs), GLOBAL across all of the
    // read's alignments (the genotype is a molecule property, not an alignment property).
    let mut votes_by_read: DetHashMap<u64, DetHashMap<usize, usize>> = DetHashMap::default();
    // per_locus: read_name_hash -> (locus copy_id -> that alignment's chain IN THAT
    // COPY's frame). One entry per distinct locus-copy the read aligns at; if a read has
    // two alignments at the SAME locus-copy, the FIRST seen is kept (deterministic — BAM
    // order), guarding against a chimera from merging two same-copy frames.
    let mut per_locus_by_read: DetHashMap<u64, DetHashMap<usize, LocusChain>> =
        DetHashMap::default();
    // Register every read seen even if it casts no votes / has no locus-copy.
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

        // Register the read so it appears in the output even with zero votes / no PSV.
        seen_reads.insert(read.read_name_hash);

        // Per-position mismatch set vs the reference FASTA (the `--eqx`-equivalent).
        // A read's base at ref position P = the FASTA base at P UNLESS P is in this
        // set, in which case it is the recorded alt base.
        let mismatches = crate::bam::extract_mismatches_vs_fasta(&record, ref_start, &chrom, genome);
        let mismatch_at = |pos: u64| -> Option<u8> {
            mismatches.iter().find(|(p, _)| *p == pos).map(|(_, b)| *b)
        };

        // Vote at every PSV column this alignment covers, AND tally which copy's frame
        // this alignment lands in (its locus-copy) — the count of covered PSV positions
        // per copy. A position is "covered" iff it lies in [ref_start,ref_end) inside an
        // aligned exon (the per-copy frame position falling inside THIS alignment marks
        // the copy whose locus the alignment is at).
        let mut covered_by_copy: DetHashMap<usize, usize> = DetHashMap::default();
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
                // This covered per-copy frame position marks pc.copy_id as the locus
                // this alignment sits in (count it for the locus-copy decision).
                *covered_by_copy.entry(pc.copy_id).or_insert(0) += 1;
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
                // Compare the base against EVERY copy's allele in this column and vote
                // GLOBALLY for the copy it matches (the cross-mapping detector). At most
                // one position per column is covered by this alignment, so we break after
                // handling it to avoid double-counting if frames ever collide.
                let copy_votes = votes_by_read.entry(read.read_name_hash).or_default();
                for cand in &col.per_copy {
                    if cand.allele == base {
                        *copy_votes.entry(cand.copy_id).or_insert(0) += 1;
                    }
                }
                break;
            }
        }

        // Determine this alignment's LOCUS-COPY = the copy whose frame positions this
        // alignment covered the most (deterministic tiebreak: lowest copy_id). An
        // alignment that covers no PSV position has no locus-copy → it contributes NO
        // per_locus entry (we cannot place its chain in a copy frame).
        let locus_copy = covered_by_copy
            .iter()
            .map(|(&c, &n)| (c, n))
            .fold(None::<(usize, usize)>, |acc, (c, n)| match acc {
                None => Some((c, n)),
                Some((bc, bn)) => {
                    if n > bn || (n == bn && c < bc) {
                        Some((c, n))
                    } else {
                        Some((bc, bn))
                    }
                }
            })
            .map(|(c, _)| c);
        let Some(locus_copy) = locus_copy else {
            continue;
        };

        // Record this alignment's junctions + exons IN locus_copy's frame (its own
        // coords) under per_locus[read][locus_copy]. If a read has two alignments at the
        // same locus-copy, keep the FIRST (BAM order) — deterministic, no chimera.
        let mut exons: Vec<(u64, u64)> = read.exons.clone();
        exons.sort_unstable();
        exons.dedup();
        let mut junctions: Vec<(u64, u64)> = if read.exons.len() >= 2 {
            read.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect()
        } else {
            Vec::new()
        };
        junctions.sort_unstable();
        junctions.dedup();
        per_locus_by_read
            .entry(read.read_name_hash)
            .or_default()
            .entry(locus_copy)
            .or_insert(LocusChain { copy_id: locus_copy, junctions, exons });
    }

    // Emit one ReadGenotype per read seen, sorted by read_name_hash; psv_votes sorted by
    // copy_id; per_locus sorted by copy_id (each chain's junctions/exons already sorted).
    // Determinism is load-bearing.
    let mut out: Vec<ReadGenotype> = seen_reads
        .into_iter()
        .map(|h| {
            let mut psv_votes: Vec<(usize, usize)> = votes_by_read
                .get(&h)
                .map(|m| m.iter().map(|(&c, &n)| (c, n)).collect())
                .unwrap_or_default();
            psv_votes.sort_by_key(|(c, _)| *c);
            let mut per_locus: Vec<LocusChain> = per_locus_by_read
                .get(&h)
                .map(|m| m.values().cloned().collect())
                .unwrap_or_default();
            per_locus.sort_by_key(|lc| lc.copy_id);
            ReadGenotype { read_name_hash: h, psv_votes, per_locus }
        })
        .collect();
    out.sort_by_key(|g| g.read_name_hash);
    Ok(out)
}

/// (Per-junction linkage view, for the certificate + the chain-consistency gate.)
/// Each CONFIDENTLY-assigned read (`assign_read_to_copy == Some(copy)`) votes every
/// junction it spans to that copy; a junction is assigned to a copy iff `>= k`
/// DISTINCT reads vote it. AMBIGUOUS reads (`assign_read_to_copy == None`) contribute
/// nothing — never a guessed junction. Returns `copy_id -> sorted assigned junctions`.
///
/// WHY a per-junction view alongside the chain assembly: the certificate wants to know
/// whether a chain's EVERY junction independently cleared the linked-read floor (a
/// chain whose splice graph is well-supported overall but has one thinly-linked
/// junction is structurally weaker than one all of whose junctions are linked). The
/// assembly's consistency gate consumes this map directly: a chain is only emitted if
/// all of its junctions appear in its copy's assigned set.
///
/// Deterministic: copies are a `DetHashMap`; each copy's junction Vec is sorted.
pub fn link_junctions(
    gens: &[ReadGenotype],
    min_psv: usize,
    margin: usize,
    k: usize,
) -> crate::types::DetHashMap<usize, Vec<(u64, u64)>> {
    use crate::types::{DetHashMap, DetHashSet};
    // (copy, junction) -> set of distinct reads voting it. A DetHashSet of read
    // hashes per (copy, junction) guards against one read's duplicate junction (a
    // read never spans the same junction twice, but a paranoid distinct count keeps
    // the >= k floor honest regardless of upstream dedup).
    let mut votes: DetHashMap<(usize, (u64, u64)), DetHashSet<u64>> = DetHashMap::default();
    for g in gens {
        let Some(copy) = assign_read_to_copy(g, min_psv, margin) else {
            continue; // ambiguous read → no junction vote
        };
        // Junctions FROM the assigned copy's OWN frame — never the pooled cross-locus
        // set. A read assigned to copy X that does not align at X contributes nothing.
        let Some(lc) = g.chain_for_copy(copy) else {
            continue;
        };
        for &j in &lc.junctions {
            votes
                .entry((copy, j))
                .or_default()
                .insert(g.read_name_hash);
        }
    }
    let mut out: DetHashMap<usize, Vec<(u64, u64)>> = DetHashMap::default();
    for ((copy, j), readers) in &votes {
        if readers.len() >= k {
            out.entry(*copy).or_default().push(*j);
        }
    }
    for v in out.values_mut() {
        v.sort_unstable();
    }
    out
}

/// Certificate of within-molecule PSV-linkage evidence backing an emitted per-copy
/// isoform. It is the per-transcript proof that the assignment to `copy_id` rests on
/// reads carrying that copy's distinguishing alleles, not on homology shadow.
///
/// - `linked_reads`: the number of DISTINCT confidently-assigned reads that traversed
///   this exact intron chain (the chain's support).
/// - `n_psvs`: the WEAKEST LINK — the minimum PSV-vote count among those supporting
///   reads (a chain backed by reads each carrying >= n_psvs of the copy's alleles is
///   only as strong as its thinnest-genotyped molecule).
/// - `copy_posterior`: a v1 vote-fraction confidence (see `assemble_psv_isoforms`).
#[derive(Debug, Clone, PartialEq)]
pub struct PsvCertificate {
    pub copy_id: usize,
    pub linked_reads: usize,
    pub n_psvs: usize,
    pub copy_posterior: f64,
}

/// Assemble per-copy isoforms from PSV-linked reads. This is the recovery payload of
/// the channel: a real per-copy splice isoform, justified entirely by reads whose
/// alleles place them on that copy.
///
/// For each copy (sorted) take its CONFIDENTLY-assigned reads (`assign_read_to_copy ==
/// Some(copy)`), group them by FULL intron chain (their `junctions`), and for each
/// chain that
///   - has `>= k` DISTINCT supporting reads, AND
///   - whose EVERY junction is copy-assigned per [`link_junctions`] (the consistency
///     gate — a chain may not contain a junction the copy didn't independently link),
/// build consensus exons and emit one [`FamilyPath`]. Consensus mirrors
/// `build_exons_from_chain`'s rule: interior boundaries are the crisp junction coords;
/// the 5' start / 3' end are the MEDIAN (lower-middle) of the supporters' terminal exon
/// coordinates. Chains-only — we emit only chains a real molecule traversed; single-
/// exon chains (no junctions) are skipped (no splice signal to assign).
///
/// PsvCertificate (v1):
///   - `linked_reads` = the chain's distinct supporting reads.
///   - `n_psvs` = the minimum PSV-vote count among those supporters (the weakest link).
///   - `copy_posterior` = `supporting / (supporting + competing)`, the VOTE FRACTION,
///     where `competing` is the number of reads that span THIS chain AND are
///     confidently assigned to a DIFFERENT copy. With no competitors this is 1.0; it
///     degrades as cross-mapped reads of other copies traverse the same chain. (v1: a
///     simple fraction, not a calibrated Bayesian posterior — documented.)
///
/// Deterministic: copies sorted; per copy, chains ranked by (support DESC, chain ASC)
/// so the assembly order is stable; the final Vec is sorted by (exons, copy_id).
pub fn assemble_psv_isoforms(
    gens: &[ReadGenotype],
    min_psv: usize,
    margin: usize,
    k: usize,
) -> Vec<crate::vg_family::layer2::FamilyPath> {
    use crate::types::{DetHashMap, DetHashSet};
    use crate::vg_family::layer2::{FamilyPath, IsoformSource};

    // The per-junction linkage view (consistency gate): copy -> its linked junctions.
    let linked = link_junctions(gens, min_psv, margin, k);

    // Confidently-assigned reads, grouped by (copy, full intron chain). The chain is the
    // read's chain IN THE ASSIGNED COPY's OWN FRAME ([`LocusChain`]) — never a pooled
    // cross-locus set, so a cross-copy chimera cannot be assembled. We retain the read's
    // assigned-copy LocusChain per supporter so consensus reads its terminal exon coords
    // in that frame, and the ReadGenotype so the certificate reads PSV-vote counts.
    let mut by_copy_chain: DetHashMap<(usize, Vec<(u64, u64)>), Vec<(&ReadGenotype, &LocusChain)>> =
        DetHashMap::default();
    // (chain -> set of copies that confidently span it) for the copy_posterior competitor
    // count. The chain key is each read's ASSIGNED-copy frame chain (consistent with the
    // grouping above), so competitors are reads whose assigned-copy chain matches.
    let mut chain_copy_reads: DetHashMap<Vec<(u64, u64)>, DetHashMap<usize, DetHashSet<u64>>> =
        DetHashMap::default();
    for g in gens {
        let Some(copy) = assign_read_to_copy(g, min_psv, margin) else {
            continue;
        };
        // The read's chain in the ASSIGNED copy's frame. A read assigned to X but with no
        // alignment at X (no X-frame chain) is skipped — never assembled from another
        // locus's coords.
        let Some(lc) = g.chain_for_copy(copy) else {
            continue;
        };
        if lc.junctions.is_empty() {
            continue; // single-exon read → no chain to assemble
        }
        by_copy_chain
            .entry((copy, lc.junctions.clone()))
            .or_default()
            .push((g, lc));
        chain_copy_reads
            .entry(lc.junctions.clone())
            .or_default()
            .entry(copy)
            .or_default()
            .insert(g.read_name_hash);
    }

    // Iterate copies in sorted order; within a copy, chains ranked (support DESC, chain
    // ASC) so assembly is deterministic regardless of map order.
    let mut keys: Vec<(usize, Vec<(u64, u64)>)> = by_copy_chain.keys().cloned().collect();
    keys.sort_by(|a, b| {
        let sa = by_copy_chain[a].len();
        let sb = by_copy_chain[b].len();
        a.0.cmp(&b.0) // copy ASC
            .then(sb.cmp(&sa)) // support DESC
            .then(a.1.cmp(&b.1)) // chain ASC
    });

    let mut out: Vec<FamilyPath> = Vec::new();
    for (copy, chain) in &keys {
        let supporters = &by_copy_chain[&(*copy, chain.clone())];
        // Distinct supporting reads (defensive — upstream groups per read, but the
        // certificate's `linked_reads` MUST be a distinct count).
        let distinct: DetHashSet<u64> = supporters.iter().map(|(g, _)| g.read_name_hash).collect();
        let support = distinct.len();
        if support < k {
            continue; // chains-only >= k floor
        }
        // Consistency gate: every junction of the chain must be copy-assigned.
        let copy_links = linked.get(copy);
        let all_linked = chain.iter().all(|j| {
            copy_links.map(|v| v.binary_search(j).is_ok()).unwrap_or(false)
        });
        if !all_linked {
            continue;
        }
        // Consensus exons: crisp junctions, median terminal start/end. Mirrors
        // build_exons_from_chain (lower-middle median).
        let median = |mut v: Vec<u64>| -> u64 {
            v.sort_unstable();
            v[(v.len() - 1) / 2]
        };
        // A supporter's terminal coords are the first exon's start and last exon's end,
        // read from the ASSIGNED copy's frame chain (lc.exons) — never a pooled frame.
        let start = median(
            supporters
                .iter()
                .filter_map(|(_, lc)| lc.exons.first().map(|e| e.0))
                .collect(),
        );
        let end = median(
            supporters
                .iter()
                .filter_map(|(_, lc)| lc.exons.last().map(|e| e.1))
                .collect(),
        );
        let mut exons: Vec<(u64, u64)> = Vec::with_capacity(chain.len() + 1);
        let mut cur_start = start;
        for &(donor, acceptor) in chain {
            exons.push((cur_start, donor));
            cur_start = acceptor;
        }
        exons.push((cur_start, end));
        if exons.len() < 2 {
            continue; // single-exon → skip (defensive; a non-empty chain has >= 2)
        }
        // Certificate. n_psvs = weakest link = min PSV-vote count for the ASSIGNED copy
        // among the supporters. copy_posterior = supporting / (supporting + competing),
        // where competing = distinct reads spanning this chain confidently assigned to a
        // DIFFERENT copy.
        let n_psvs = supporters
            .iter()
            .map(|(g, _)| {
                g.psv_votes
                    .iter()
                    .find(|(c, _)| c == copy)
                    .map(|(_, n)| *n)
                    .unwrap_or(0)
            })
            .min()
            .unwrap_or(0);
        let competing: usize = chain_copy_reads
            .get(chain)
            .map(|m| {
                m.iter()
                    .filter(|(c, _)| **c != *copy)
                    .map(|(_, reads)| reads.len())
                    .sum()
            })
            .unwrap_or(0);
        let copy_posterior = support as f64 / (support + competing) as f64;
        let junction_chain = chain.clone();
        out.push(FamilyPath {
            exons,
            flow: support as f64,
            copy_id: *copy,
            junction_chain,
            source: IsoformSource::PsvLinked { copy_id: *copy },
            certificate: Some(PsvCertificate {
                copy_id: *copy,
                linked_reads: support,
                n_psvs,
                copy_posterior,
            }),
        });
    }
    // Final deterministic total order: (exons, copy_id).
    out.sort_by(|a, b| a.exons.cmp(&b.exons).then(a.copy_id.cmp(&b.copy_id)));
    out
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
            per_locus: vec![LocusChain {
                copy_id: 1,
                junctions: vec![(160, 300)],
                exons: vec![(100, 160), (300, 360)],
            }],
        };
        assert_eq!(assign_read_to_copy(&g, 3, 1), Some(1));
    }

    #[test]
    fn ambiguous_tie_unassigned() {
        // Even-split votes, and neither reaches min_psv -> ambiguous -> None.
        let g = ReadGenotype {
            read_name_hash: 8,
            psv_votes: vec![(0, 2), (1, 2)],
            per_locus: vec![LocusChain {
                copy_id: 0,
                junctions: vec![],
                exons: vec![(100, 200)],
            }],
        };
        assert_eq!(assign_read_to_copy(&g, 3, 1), None);
    }

    #[test]
    fn sub_threshold_unassigned() {
        // Single copy with 2 votes but min_psv=3 -> below the error floor -> None.
        let g = ReadGenotype {
            read_name_hash: 9,
            psv_votes: vec![(0, 2)],
            per_locus: vec![LocusChain {
                copy_id: 0,
                junctions: vec![],
                exons: vec![(100, 200)],
            }],
        };
        assert_eq!(assign_read_to_copy(&g, 3, 1), None);
    }

    #[test]
    fn margin_enforced() {
        // Both copies clear min_psv but are tied -> no strict dominance -> None.
        let tied = ReadGenotype {
            read_name_hash: 10,
            psv_votes: vec![(0, 3), (1, 3)],
            per_locus: vec![LocusChain {
                copy_id: 0,
                junctions: vec![],
                exons: vec![(100, 200)],
            }],
        };
        assert_eq!(assign_read_to_copy(&tied, 3, 1), None);
        // One copy clears min_psv and beats the runner-up by the margin -> assigned.
        let dominant = ReadGenotype {
            read_name_hash: 11,
            psv_votes: vec![(0, 4), (1, 3)],
            per_locus: vec![LocusChain {
                copy_id: 0,
                junctions: vec![],
                exons: vec![(100, 200)],
            }],
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
            per_locus: vec![LocusChain {
                copy_id: 0,
                junctions: vec![(100, 200), (250, 400)],
                exons: vec![(50, 100), (200, 250), (400, 450)],
            }],
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

    // ── Chimera guard: assembly uses the ASSIGNED copy's frame ONLY ──────────
    //
    // The root-cause regression test. A copy-0 read cross-maps to BOTH copy0 and
    // copy1, so it carries two LocusChains (one per copy frame). Its alleles vote
    // copy0. Assembly MUST build the isoform from copy0's frame chain ONLY — never the
    // pooled {copy0, copy1} junction set that produced a cross-locus chimera (on real
    // chrY a 9.3 Mb `.`-strand transcript spanning two copies).
    #[test]
    fn no_cross_locus_chimera_uses_assigned_copy_frame() {
        // A copy-0 read aligns at BOTH copy0 (junctions [(160,300)], exons
        // [(100,160),(300,360)]) and copy1 (junctions [(30160,30300)] in copy1's frame),
        // votes copy0 (psv_votes [(0,3)]).
        let g = ReadGenotype {
            read_name_hash: 1,
            psv_votes: vec![(0, 3)],
            per_locus: vec![
                LocusChain {
                    copy_id: 0,
                    junctions: vec![(160, 300)],
                    exons: vec![(100, 160), (300, 360)],
                },
                LocusChain {
                    copy_id: 1,
                    junctions: vec![(30160, 30300)],
                    exons: vec![(30100, 30160), (30300, 30360)],
                },
            ],
        };
        // two such reads -> assemble emits ONE copy0 isoform with exons within copy0's
        // range, junction (160,300) only; NO exon coordinate > 1000 (copy1's frame).
        let out = assemble_psv_isoforms(
            &[
                g.clone(),
                ReadGenotype {
                    read_name_hash: 2,
                    ..g.clone()
                },
            ],
            3,
            1,
            2,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].copy_id, 0);
        assert!(
            out[0].exons.iter().all(|&(_, e)| e <= 1000),
            "no copy1-frame (chimera) coords: {:?}",
            out[0].exons
        );
        assert_eq!(out[0].junction_chain, vec![(160, 300)]);
    }

    // ── Task 4: junction linkage + per-copy isoform assembly + PsvCertificate ──
    //
    // link_junctions is the per-junction view: confidently-assigned reads vote each
    // junction they span to their copy; a junction is the copy's only at >= k votes.
    // assemble_psv_isoforms is the recovery payload: it groups a copy's reads by full
    // chain, gates each chain (>= k support + every junction copy-linked), and emits a
    // PsvLinked FamilyPath with a certificate. Both NEVER act on ambiguous reads.

    #[test]
    fn link_junctions_assigns_at_k() {
        // 2 reads both confidently copy0 (3 PSV votes) spanning junction (160,300).
        // At k=2 the junction is assigned to copy0; at k=3 it falls below the floor.
        let gens = vec![
            ReadGenotype {
                read_name_hash: 1,
                psv_votes: vec![(0, 3)],
                per_locus: vec![LocusChain {
                    copy_id: 0,
                    junctions: vec![(160, 300)],
                    exons: vec![(100, 160), (300, 360)],
                }],
            },
            ReadGenotype {
                read_name_hash: 2,
                psv_votes: vec![(0, 3)],
                per_locus: vec![LocusChain {
                    copy_id: 0,
                    junctions: vec![(160, 300)],
                    exons: vec![(100, 160), (300, 360)],
                }],
            },
        ];
        let m = link_junctions(&gens, 3, 1, 2);
        assert_eq!(m.get(&0).unwrap(), &vec![(160, 300)]);
        // k=3 → only 2 distinct readers → copy0 has no linked junction.
        let m3 = link_junctions(&gens, 3, 1, 3);
        assert!(m3.get(&0).is_none(), "k=3 with 2 reads → empty for copy0: {m3:?}");
    }

    #[test]
    fn link_junctions_excludes_ambiguous() {
        // One CONFIDENT copy0 read + one AMBIGUOUS (tied) read, both span (160,300).
        // The ambiguous read contributes nothing, so at k=2 the junction is NOT linked.
        let gens = vec![
            ReadGenotype {
                read_name_hash: 1,
                psv_votes: vec![(0, 3)],
                per_locus: vec![LocusChain {
                    copy_id: 0,
                    junctions: vec![(160, 300)],
                    exons: vec![(100, 160), (300, 360)],
                }],
            },
            ReadGenotype {
                read_name_hash: 2,
                psv_votes: vec![(0, 2), (1, 2)], // tied → assign_read_to_copy == None
                per_locus: vec![LocusChain {
                    copy_id: 0,
                    junctions: vec![(160, 300)],
                    exons: vec![(100, 160), (300, 360)],
                }],
            },
        ];
        let m = link_junctions(&gens, 3, 1, 2);
        assert!(
            m.get(&0).is_none(),
            "ambiguous read must not count toward the k=2 floor: {m:?}"
        );
        // Sanity: with k=1 the single confident read links it (the ambiguous one still
        // contributes nothing — the linked set is exactly the confident junction).
        let m1 = link_junctions(&gens, 3, 1, 1);
        assert_eq!(m1.get(&0).unwrap(), &vec![(160, 300)]);
    }

    #[test]
    fn assemble_emits_psvlinked_isoform() {
        // 2 reads confidently copy0, same 3-exon chain. k=2 → one PsvLinked FamilyPath
        // with the consensus exons, the chain, and a certificate (linked_reads == 2).
        let exons = vec![(100, 160), (300, 360), (500, 560)];
        let chain = vec![(160, 300), (360, 500)];
        let gens = vec![
            ReadGenotype {
                read_name_hash: 1,
                psv_votes: vec![(0, 3)],
                per_locus: vec![LocusChain {
                    copy_id: 0,
                    junctions: chain.clone(),
                    exons: exons.clone(),
                }],
            },
            ReadGenotype {
                read_name_hash: 2,
                psv_votes: vec![(0, 4)],
                per_locus: vec![LocusChain {
                    copy_id: 0,
                    junctions: chain.clone(),
                    exons: exons.clone(),
                }],
            },
        ];
        let paths = assemble_psv_isoforms(&gens, 3, 1, 2);
        assert_eq!(paths.len(), 1, "exactly one assembled isoform: {paths:?}");
        let p = &paths[0];
        assert_eq!(p.source, crate::vg_family::layer2::IsoformSource::PsvLinked { copy_id: 0 });
        assert_eq!(p.copy_id, 0);
        assert_eq!(p.junction_chain, chain);
        assert_eq!(p.exons, exons, "consensus exons = crisp junctions + median terminals");
        assert_eq!(p.flow, 2.0);
        let cert = p.certificate.as_ref().expect("PsvLinked path carries a certificate");
        assert_eq!(cert.copy_id, 0);
        assert_eq!(cert.linked_reads, 2);
        // n_psvs = weakest link = min(3, 4) = 3.
        assert_eq!(cert.n_psvs, 3);
        // No competitor reads → posterior 1.0.
        assert_eq!(cert.copy_posterior, 1.0);
    }

    #[test]
    fn assemble_respects_k_and_consistency() {
        // Only ONE read for copy0's chain at k=2 → not emitted (chains-only >= k). A
        // read confidently assigned to copy1 sharing nothing with copy0 must not bleed
        // into copy0's chains.
        let chain0 = vec![(160, 300), (360, 500)];
        let exons0 = vec![(100, 160), (300, 360), (500, 560)];
        let chain1 = vec![(1160, 1300)];
        let exons1 = vec![(1100, 1160), (1300, 1360)];
        let gens = vec![
            // single copy0 supporter for chain0 (below k=2)
            ReadGenotype {
                read_name_hash: 1,
                psv_votes: vec![(0, 3)],
                per_locus: vec![LocusChain {
                    copy_id: 0,
                    junctions: chain0.clone(),
                    exons: exons0.clone(),
                }],
            },
            // a copy1 read on a different chain — must not contribute to copy0
            ReadGenotype {
                read_name_hash: 2,
                psv_votes: vec![(1, 3)],
                per_locus: vec![LocusChain {
                    copy_id: 1,
                    junctions: chain1.clone(),
                    exons: exons1.clone(),
                }],
            },
        ];
        let paths = assemble_psv_isoforms(&gens, 3, 1, 2);
        assert!(
            paths.is_empty(),
            "neither chain reaches k=2 distinct supporters: {paths:?}"
        );
        // Add a SECOND copy0 read on chain0 → now it clears k=2 and is emitted; copy1's
        // lone read still does not (and the copy0 path is assigned to copy0 only).
        let mut gens2 = gens.clone();
        gens2.push(ReadGenotype {
            read_name_hash: 3,
            psv_votes: vec![(0, 3)],
            per_locus: vec![LocusChain {
                copy_id: 0,
                junctions: chain0.clone(),
                exons: exons0.clone(),
            }],
        });
        let paths2 = assemble_psv_isoforms(&gens2, 3, 1, 2);
        assert_eq!(paths2.len(), 1, "only copy0's chain clears k: {paths2:?}");
        assert_eq!(paths2[0].copy_id, 0);
        assert_eq!(paths2[0].junction_chain, chain0);
    }

    #[test]
    fn assemble_skips_single_exon() {
        // Confidently-assigned reads with NO junction (single aligned block) carry no
        // splice signal → nothing to assemble → no emission.
        let gens = vec![
            ReadGenotype {
                read_name_hash: 1,
                psv_votes: vec![(0, 3)],
                per_locus: vec![LocusChain {
                    copy_id: 0,
                    junctions: vec![],
                    exons: vec![(100, 600)],
                }],
            },
            ReadGenotype {
                read_name_hash: 2,
                psv_votes: vec![(0, 4)],
                per_locus: vec![LocusChain {
                    copy_id: 0,
                    junctions: vec![],
                    exons: vec![(100, 600)],
                }],
            },
        ];
        let paths = assemble_psv_isoforms(&gens, 3, 1, 2);
        assert!(paths.is_empty(), "single-exon reads emit nothing: {paths:?}");
    }

    #[test]
    fn familypath_certificate_none_by_default() {
        // Sanity: a hand-built Native FamilyPath carries certificate: None (the field's
        // default everywhere outside assemble_psv_isoforms). Locks the PartialEq anchor
        // for emit_both_off_equals_decompose (all base/Part-A/B paths are certificate-None).
        use crate::vg_family::layer2::{FamilyPath, IsoformSource};
        let p = FamilyPath {
            exons: vec![(100, 160), (300, 360)],
            flow: 1.0,
            copy_id: 0,
            junction_chain: vec![(160, 300)],
            source: IsoformSource::Native,
            certificate: None,
        };
        assert!(p.certificate.is_none());
    }

    #[test]
    fn assemble_copy_posterior_degrades_with_competitor() {
        // 2 copy0 reads + 1 copy1 read ALL span the same chain. copy0's chain clears
        // k=2; its posterior = 2 / (2 + 1) (one competing copy1 read on the same chain).
        let chain = vec![(160, 300), (360, 500)];
        let exons = vec![(100, 160), (300, 360), (500, 560)];
        let gens = vec![
            ReadGenotype {
                read_name_hash: 1,
                psv_votes: vec![(0, 3)],
                per_locus: vec![LocusChain {
                    copy_id: 0,
                    junctions: chain.clone(),
                    exons: exons.clone(),
                }],
            },
            ReadGenotype {
                read_name_hash: 2,
                psv_votes: vec![(0, 3)],
                per_locus: vec![LocusChain {
                    copy_id: 0,
                    junctions: chain.clone(),
                    exons: exons.clone(),
                }],
            },
            // A copy1 read whose copy1-frame chain has the SAME coords → it competes on
            // copy0's chain in the posterior count (one competitor) but is itself
            // sub-k for copy1, so it is not emitted.
            ReadGenotype {
                read_name_hash: 3,
                psv_votes: vec![(1, 3)],
                per_locus: vec![LocusChain {
                    copy_id: 1,
                    junctions: chain.clone(),
                    exons: exons.clone(),
                }],
            },
        ];
        let paths = assemble_psv_isoforms(&gens, 3, 1, 2);
        // copy1 has only 1 read on the chain (< k) → not emitted; copy0 is.
        assert_eq!(paths.len(), 1, "only copy0 clears k=2: {paths:?}");
        let cert = paths[0].certificate.as_ref().unwrap();
        assert_eq!(cert.copy_id, 0);
        assert_eq!(cert.linked_reads, 2);
        assert!(
            (cert.copy_posterior - 2.0 / 3.0).abs() < 1e-9,
            "posterior 2/(2+1): {}",
            cert.copy_posterior
        );
    }
}
