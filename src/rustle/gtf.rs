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
        if let Some(rc) = tx.rescue_class {
            tx_attrs.push_str(&format!(" rescue_class \"{}\";", rc));
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
        crate::parity_decisions::emit(
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
