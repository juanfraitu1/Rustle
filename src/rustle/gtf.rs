//! DEAD CODE -- StringTie-era assembler; slated for removal.
//! NOT part of the multi-copy-family thesis (O1 family-def / O2 copy-assign / O3 ASJ / O4 absent-copy).
//! See docs/RETIREMENT_AND_MIGRATION.md. Do not extend.
//! Write transcripts to GTF.

use std::io::Write;

use crate::path_extract::Transcript;
use crate::types::DetHashMap as HashMap;

fn tx_bounds(tx: &Transcript) -> Option<(u64, u64)> {
    Some((tx.exons.first()?.0, tx.exons.last()?.1))
}

/// Returns true if any exon of `a` overlaps (shares ≥1 base with) any exon of `b`.
/// Mirrors `overlaps.get(i,j)` / `update_overlap()` in print_predcluster.
fn exons_overlap(a: &[(u64, u64)], b: &[(u64, u64)]) -> bool {
    for &(as_, ae) in a {
        for &(bs, be) in b {
            if as_ < be && bs < ae {
                return true;
            }
        }
    }
    false
}

fn uf_find(parent: &mut Vec<usize>, x: usize) -> usize {
    let mut x = x;
    while parent[x] != x {
        parent[x] = parent[parent[x]]; // path compression
        x = parent[x];
    }
    x
}

/// Assign sequential gene/transcript numbering using Union-Find on exon overlap.
/// Two transcripts belong to the same gene iff their exons share at least one base —
/// matching print_predcluster (~line 19230) which uses overlaps.get(i,j).
/// Returns (gene_no, transcript_no_within_gene) per input transcript index.
pub fn assign_gene_tx_numbers(transcripts: &[Transcript]) -> Vec<(usize, usize)> {
    let n = transcripts.len();
    let mut out = vec![(0usize, 0usize); n];
    if n == 0 {
        return out;
    }
    // Sort by chrom, strand, start for stable left-to-right processing.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        let ta = &transcripts[a];
        let tb = &transcripts[b];
        let sa = ta.exons.first().map(|e| e.0).unwrap_or(0);
        let sb = tb.exons.first().map(|e| e.0).unwrap_or(0);
        ta.chrom
            .cmp(&tb.chrom)
            .then_with(|| ta.strand.cmp(&tb.strand))
            .then_with(|| sa.cmp(&sb))
    });

    // Union-Find over sorted positions (parent[i] = parent of order[i]).
    let mut parent: Vec<usize> = (0..n).collect();

    // Precompute bounds for span-overlap quick reject.
    let bounds: Vec<(u64, u64)> = order
        .iter()
        .map(|&idx| tx_bounds(&transcripts[idx]).unwrap_or((0, 0)))
        .collect();

    for i in 0..n {
        let (si, ei) = bounds[i];
        let ti = &transcripts[order[i]];
        // Scan backward while spans can still overlap with i.
        // Since sorted by start, j < i means bounds[j].start <= si.
        // Overlap possible when bounds[j].end > si.
        // We track the prefix-max of end to allow an early-exit:
        // once max_end_so_far[j] <= si we know no earlier j overlaps i.
        // Simple O(n²) scan is fine for typical bundle sizes.
        for j in (0..i).rev() {
            let (sj, ej) = bounds[j];
            let tj = &transcripts[order[j]];
            // Span-level quick reject / early exit.
            if ej <= si || sj >= ei {
                // Since sorted by start and going backward,
                // if ej <= si then all earlier j have starts ≤ sj ≤ si,
                // but their ends might be larger. Only safe to break if we
                // track prefix max end; skip break and just continue.
                continue;
            }
            if tj.chrom != ti.chrom || tj.strand != ti.strand {
                break; // different chrom/strand block — stop
            }
            // Check actual exon-level overlap (overlaps.get).
            if exons_overlap(&ti.exons, &tj.exons) {
                let ri = uf_find(&mut parent, i);
                let rj = uf_find(&mut parent, j);
                if ri != rj {
                    // Merge larger root into smaller (stable numbering).
                    parent[ri.max(rj)] = ri.min(rj);
                }
            }
        }
    }

    // Collect components and assign deterministic gene numbers by first member position.
    let mut root_to_gene: HashMap<usize, usize> = Default::default();
    let mut gene_no = 0usize;
    // Walk in sorted order so gene numbers increase left-to-right.
    for i in 0..n {
        let root = uf_find(&mut parent, i);
        if !root_to_gene.contains_key(&root) {
            gene_no += 1;
            root_to_gene.insert(root, gene_no);
        }
    }

    // Assign transcript numbers within each gene (sorted by original transcript index).
    let mut gene_members: HashMap<usize, Vec<usize>> = Default::default();
    for i in 0..n {
        let root = uf_find(&mut parent, i);
        let gno = root_to_gene[&root];
        gene_members.entry(gno).or_default().push(order[i]);
    }
    for (&gno, members) in &mut gene_members {
        let mut m = members.clone();
        m.sort_unstable();
        for (k, &idx) in m.iter().enumerate() {
            out[idx] = (gno, k + 1);
        }
    }
    out
}

/// Write transcripts to GTF (1-based coordinates) in compatible format.
/// Emits transcript + exon lines with cov/FPKM/TPM attributes.
pub fn write_gtf<W: Write>(
    transcripts: &[Transcript],
    writer: &mut W,
    label: &str,
) -> std::io::Result<()> {
    let gene_tx_no = assign_gene_tx_numbers(transcripts);
    // Abstain floor for the capacity-confidence channel (matches the pipeline
    // RUSTLE_VG_ABSTAIN_FLOOR default). A VG transcript whose capacity_confidence
    // falls below this is an unreliable (phantom / non-anchored) copy: we emit an
    // explicit, filterable low_confidence tag (the copy is kept, not dropped).
    let abstain_floor: f64 = std::env::var("RUSTLE_VG_ABSTAIN_FLOOR")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.05);
    for (i, tx) in transcripts.iter().enumerate() {
        let (gno, tno) = gene_tx_no[i];
        let (gid, tid) = if let (Some(ref g), Some(ref t)) = (&tx.gene_id, &tx.transcript_id) {
            (g.clone(), t.clone())
        } else {
            (
                format!("{}.{}", label, gno.max(1)),
                format!("{}.{}.{}", label, gno.max(1), tno.max(1)),
            )
        };
        let t_start = tx.exons.first().map(|(s, _)| s + 1).unwrap_or(0);
        let t_end = tx.exons.last().map(|(_, e)| *e).unwrap_or(0);
        let source_col = label;
        // Transcript line (the original algorithm format)
        let mut tx_attrs = format!(
            "gene_id \"{}\"; transcript_id \"{}\"; cov \"{:.6}\"; FPKM \"{:.6}\"; TPM \"{:.6}\";",
            gid, tid, tx.coverage, tx.fpkm, tx.tpm
        );
        if let Some(ref s) = tx.source {
            tx_attrs.push_str(&format!(" source \"{}\";", s));
        }
        if tx.is_longread {
            tx_attrs.push_str(&format!(" longcov \"{:.6}\";", tx.longcov));
        }
        if let Some(ref rid) = tx.ref_transcript_id {
            tx_attrs.push_str(&format!(" reference_id \"{}\";", rid));
        }
        if let Some(ref rgid) = tx.ref_gene_id {
            tx_attrs.push_str(&format!(" ref_gene_id \"{}\";", rgid));
        }
        // Variation graph annotations (--vg mode).
        if let Some(fam_id) = tx.vg_family_id {
            tx_attrs.push_str(&format!(" family_id \"FAM_{}\";", fam_id));
        }
        if let Some(copy_id) = tx.vg_copy_id {
            tx_attrs.push_str(&format!(" copy_id \"{}\";", copy_id));
        }
        if let Some(fam_size) = tx.vg_family_size {
            tx_attrs.push_str(&format!(" family_size \"{}\";", fam_size));
        }
        // Phased assembly (--vg-phase): haplotype + phase-set. Only set on transcripts from
        // a phase-split sub-bundle; absent otherwise, so default runs are byte-identical.
        if let Some(hp) = tx.hp_tag {
            tx_attrs.push_str(&format!(" haplotype \"{}\";", hp));
        }
        if let Some(ps) = tx.ps_tag {
            tx_attrs.push_str(&format!(" phase_set \"{}\";", ps));
        }
        // copy_confidence = ATTRIBUTION certainty: fraction of this copy's ambiguous reads
        // that are EVIDENCE-decisive winners (pre-prior margin; vg.rs compute_per_copy_confidence).
        // This is the metric to read for "how sure is the copy assignment". Distinct from
        // capacity_confidence below, which is COVERAGE/flow adequacy, NOT attribution (it can be
        // high on a non-identifiable dominant copy — the RBMY mis-calibration).
        if let Some(conf) = tx.copy_assignment_confidence {
            tx_attrs.push_str(&format!(" copy_confidence \"{:.3}\";", conf));
        }
        // Certified copy-support fraction (--vg copy-support guard). Only set on
        // VG transcripts; non-VG transcripts (vg_copy_id None) never carry it.
        if let Some(supp) = tx.copy_independent_support {
            tx_attrs.push_str(&format!(" copy_independent_support \"{:.3}\";", supp));
        }
        // Capacity-confidence channel (--vg flow-apportionment spec 2026-06-01).
        // Per-transcript anchored/total coverage fraction; only set in VG mode,
        // so non-VG GTF is byte-identical. Mirrors copy_confidence formatting.
        if let Some(cc) = tx.capacity_confidence {
            tx_attrs.push_str(&format!(" capacity_confidence \"{:.3}\";", cc));
            if cc < abstain_floor {
                tx_attrs.push_str(" low_confidence \"true\";");
            }
        }
        // O5/tandem provenance flag (--vg only; spec 2026-06-02). VG-only so
        // non-VG GTF is byte-identical.
        if tx.tandem_copy == Some(true) {
            tx_attrs.push_str(" tandem_copy \"true\";");
        }
        // Jointly-feasible lower bound on abundance (coverage * capacity_confidence).
        if let Some(amin) = tx.abundance_min {
            tx_attrs.push_str(&format!(" abundance_min \"{:.3}\";", amin));
        }
        // Multi-copy family classification verdict (--vg only).
        if let Some(ref v) = tx.family_verdict {
            tx_attrs.push_str(&format!(
                " family_verdict \"{}\"; family_identifiability \"{}\"; family_n_copies \"{}\"; family_n_expressed \"{}\"; family_locus_rel \"{}\";",
                v.class.as_str(), v.identifiability.as_str(), v.n_copies, v.n_expressed, v.locus_rel.as_str()));
            // Audit lever #3: depth-inferred copy number (only when RUSTLE_VG_DEPTH_COPYNUM
            // populated it; default runs leave it None, so the GTF is byte-identical).
            if let Some(dc) = v.depth_copies {
                tx_attrs.push_str(&format!(" family_depth_copies \"{:.2}\";", dc));
            }
            // Audit lever #2: per-copy EM attribution confidence (decisive_frac is the
            // headline; the bucket counts + means expose the identifiability detail).
            if let Some(ref ec) = v.em_confidence {
                tx_attrs.push_str(&format!(
                    " em_decisive_frac \"{:.3}\"; em_n_decisive \"{}\"; em_n_moderate \"{}\"; em_n_uncertain \"{}\"; em_mean_gap \"{:.3}\"; em_mean_sites \"{:.1}\";",
                    ec.decisive_frac(), ec.n_decisive, ec.n_moderate, ec.n_uncertain, ec.mean_gap, ec.mean_sites));
            }
            // Segmental-duplication call (opt-in RUSTLE_VG_SEGDUP_EXTENT; absent otherwise).
            if let Some(ref sd) = v.segdup {
                tx_attrs.push_str(&format!(
                    " family_segdup \"{}\"; family_segdup_extent \"{}\"; family_flank_up \"{}\"; family_flank_down \"{}\";",
                    if sd.is_segdup { "true" } else { "false" }, sd.total_extent, sd.upstream_extent, sd.downstream_extent));
            }
            // Gene-conversion event (opt-in RUSTLE_VG_MOSAIC_ON; absent otherwise). Confirmed =
            // breakpoint recurred across ≥k molecules; otherwise a chimera-suspect candidate.
            if let Some(ref c) = v.conversion {
                tx_attrs.push_str(&format!(
                    " gene_conversion \"{}\"; conversion_copies \"{}>{}\"; conversion_breakpoint \"{}-{}\"; conversion_reads \"{}\";",
                    if c.confirmed { "confirmed" } else { "chimera_suspect" },
                    c.copy_a, c.copy_b, c.breakpoint_ref.0, c.breakpoint_ref.1, c.n_supporting_reads));
                // The emitted recombinant ISOFORM (opt-in RUSTLE_VG_MOSAIC_EMIT) is tagged
                // source="gene_conversion"; mark it so it is distinguishable from the native
                // transcripts that merely carry the locus-level conversion flag.
                if tx.source.as_deref() == Some("gene_conversion") {
                    tx_attrs.push_str(" conversion_isoform \"recombinant\";");
                }
            }
            // Hidden-copy evidence (opt-in RUSTLE_VG_HIDDEN_COPY; absent otherwise). A flag that
            // the reads imply a copy not in the reference — NOT a placed/fabricated copy.
            if let Some(ref h) = v.hidden_copy {
                tx_attrs.push_str(&format!(
                    " hidden_copy_evidence \"true\"; hidden_copy_alt_positions \"{}\"; hidden_copy_read_frac \"{:.2}\";",
                    h.n_alt_positions, h.alt_read_fraction));
            }
        }
        if let Some(rc) = tx.rescue_class {
            tx_attrs.push_str(&format!(" rescue_class \"{}\";", rc));
        }
        // Mark transcripts produced from a synthetic bundle (HMM rescue of
        // unmapped or mask-region reads) so they're identifiable as
        // novel-copy candidates from objective 4 (--vg-discover-novel).
        if tx.synthetic {
            tx_attrs.push_str(" copy_status \"novel\";");
        }
        writeln!(
            writer,
            "{}\t{}\ttranscript\t{}\t{}\t1000\t{}\t.\t{}",
            tx.chrom, source_col, t_start, t_end, tx.strand, tx_attrs
        )?;
        // Parity-decisions: emit one path_emit row per final transcript so
        // the rustle/StringTie outputs can be diffed by (start, end, strand)
        // to surface "rustle-only" predictions and cluster their causes.
        // Coordinates are 1-based inclusive (same as the GTF line above) so
        // they align with StringTie's pd_emit on the C side.
        let exons_csv: String = tx.exons.iter()
            .map(|(s, e)| format!("{}-{}", s + 1, e))
            .collect::<Vec<_>>()
            .join(",");
        let payload = format!(
            r#""nexons":{},"cov":{:.6},"longcov":{:.6},"source":"{}","exons":"{}""#,
            tx.exons.len(),
            tx.coverage,
            tx.longcov,
            tx.source.as_deref().unwrap_or(""),
            exons_csv,
        );
        crate::parity::decisions::emit(
            "path_emit",
            Some(&tx.chrom),
            t_start,
            t_end,
            tx.strand,
            &payload,
        );
        for (j, (start, end)) in tx.exons.iter().enumerate() {
            // Match the original algorithm formatting: report per-exon coverage as stored in `exon_cov` (or fall back to tx.coverage).
            let ecov = tx
                .exon_cov
                .get(j)
                .copied()
                .filter(|v| *v > 0.0)
                .unwrap_or(tx.coverage);
            let exon_attrs = format!(
                "gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\"; cov \"{:.6}\";",
                gid,
                tid,
                j + 1,
                ecov
            );
            writeln!(
                writer,
                "{}\t{}\texon\t{}\t{}\t1000\t{}\t.\t{}",
                tx.chrom,
                source_col,
                start + 1,
                end,
                tx.strand,
                exon_attrs
            )?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::path_extract::Transcript;
    use crate::vg::{FamilyClass, FamilyVerdict, Identifiability, LocusRel};

    fn one_exon_tx() -> Transcript {
        Transcript {
            chrom: "chrT".to_string(),
            strand: '+',
            exons: vec![(100, 200)],
            coverage: 12.0,
            ..Default::default()
        }
    }

    fn render(tx: Transcript) -> String {
        let mut buf: Vec<u8> = Vec::new();
        write_gtf(&[tx], &mut buf, "STRG").unwrap();
        String::from_utf8(buf).unwrap()
    }

    #[test]
    fn gtf_emits_haplotype_attrs_when_phased() {
        let mut tx = one_exon_tx();
        tx.hp_tag = Some(1);
        tx.ps_tag = Some(7);
        let out = render(tx);
        assert!(out.contains("haplotype \"1\";"), "missing haplotype attr in:\n{out}");
        assert!(out.contains("phase_set \"7\";"), "missing phase_set attr in:\n{out}");
    }

    #[test]
    fn gtf_omits_haplotype_attrs_when_unphased() {
        let out = render(one_exon_tx()); // hp_tag/ps_tag default None
        assert!(!out.contains("haplotype "), "unphased tx must not emit haplotype in:\n{out}");
        assert!(!out.contains("phase_set "), "unphased tx must not emit phase_set in:\n{out}");
    }

    #[test]
    fn capacity_attrs_emitted_when_some() {
        let mut tx = one_exon_tx();
        tx.capacity_confidence = Some(0.250);
        tx.abundance_min = Some(3.000);
        let out = render(tx);
        assert!(out.contains("capacity_confidence \"0.250\";"), "missing capacity_confidence in:\n{out}");
        assert!(out.contains("abundance_min \"3.000\";"), "missing abundance_min in:\n{out}");
    }

    #[test]
    fn capacity_attrs_absent_when_none() {
        let tx = one_exon_tx(); // both fields default to None (non-vg)
        let out = render(tx);
        assert!(!out.contains("capacity_confidence"), "capacity_confidence leaked into non-vg GTF:\n{out}");
        assert!(!out.contains("abundance_min"), "abundance_min leaked into non-vg GTF:\n{out}");
        assert!(!out.contains("low_confidence"), "low_confidence must not appear in non-vg GTF:\n{out}");
    }

    #[test]
    fn low_confidence_tag_below_abstain_floor() {
        std::env::remove_var("RUSTLE_VG_ABSTAIN_FLOOR"); // default 0.05
        // capacity_confidence below floor -> low_confidence tag emitted.
        let mut lo = one_exon_tx(); lo.capacity_confidence = Some(0.000);
        let out_lo = render(lo);
        assert!(out_lo.contains("low_confidence \"true\";"),
            "low cc must emit low_confidence tag:\n{out_lo}");
        // capacity_confidence above floor -> no tag.
        let mut hi = one_exon_tx(); hi.capacity_confidence = Some(1.000);
        let out_hi = render(hi);
        assert!(!out_hi.contains("low_confidence"),
            "high cc must NOT emit low_confidence tag:\n{out_hi}");
    }

    #[test]
    fn tandem_copy_attr_emitted_when_some_true() {
        // Emission depends only on tandem_copy == Some(true); vg_copy_id is NOT
        // consulted by the guard (line 205). Verify both with and without
        // vg_copy_id so a stray guard would be caught.
        let mut tx = one_exon_tx();
        tx.tandem_copy = Some(true);
        let out = render(tx);
        assert!(out.contains("tandem_copy \"true\";"), "missing tandem_copy (vg_copy_id=None):\n{out}");

        let mut tx2 = one_exon_tx();
        tx2.vg_copy_id = Some(0);
        tx2.tandem_copy = Some(true);
        let out2 = render(tx2);
        assert!(out2.contains("tandem_copy \"true\";"), "missing tandem_copy (vg_copy_id=Some(0)):\n{out2}");
    }

    #[test]
    fn tandem_copy_attr_absent_otherwise() {
        // None -> absent
        let tx = one_exon_tx();
        let out = render(tx);
        assert!(!out.contains("tandem_copy"), "tandem_copy must not appear when None:\n{out}");
        // Some(false) -> absent
        let mut tx2 = one_exon_tx();
        tx2.tandem_copy = Some(false);
        let out2 = render(tx2);
        assert!(!out2.contains("tandem_copy"), "tandem_copy must not appear when Some(false):\n{out2}");
    }

    #[test]
    fn abstained_tx_retained_and_tagged() {
        // Simulate the pipeline abstain tag: a low-cc tx with no prior verdict
        // gets a synthesized verdict whose identifiability == None. The tx is NOT
        // dropped, so its coverage survives and the GTF carries the tag.
        let mut tx = one_exon_tx();
        tx.capacity_confidence = Some(0.010);
        tx.family_verdict = Some(FamilyVerdict {
            class: FamilyClass::Spillover,
            n_copies: 1,
            n_expressed: 0,
            connectivity: 0.0,
            identifiability: Identifiability::None,
            n_id_classes: 0,
            locus_rel: LocusRel::Single,
            depth_copies: None,
            em_confidence: None, segdup: None, conversion: None, hidden_copy: None,
        });
        let out = render(tx);
        assert!(out.contains("cov \"12.000000\";"), "abstained tx lost its coverage:\n{out}");
        assert!(out.contains("family_identifiability \"none\";"), "abstain tag missing:\n{out}");
    }
}
