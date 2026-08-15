//! O1 CATALOG → O2 INPUT: parse a `gw_family_catalog` `<out>.copies.tsv` (+ the optional
//! `<out>.copies.fa`) back into the copy set `copy_assign` assigns reads to.
//!
//! # Why this module exists
//!
//! O1 (`gw_family_catalog`) and O2 (`copy_assign`) already share ONE node type
//! (`family_detect::DenovoTranscript` = catalog row = O2 copy), ONE rep-build front end, ONE edge engine
//! and ONE admission primitive — but they shared all of it BY FUNCTION CALL and NOTHING BY FILE. Each
//! binary re-derived its own families from the BAM, so at defaults they built DIFFERENT objects (measured:
//! GSTM catalog 4 copies vs `copy_assign` 0 families on 6031 reads) and their family ids (`GWFAM{i}` vs
//! `CAFAM{i}`) were assigned independently, leaving NO JOIN KEY between the two tables.
//!
//! This module is the file-level contract that closes that gap: `copy_assign --families <copies.tsv>`
//! CONSUMES the O1 roster instead of re-deriving one, so the two objects agree BY CONSTRUCTION rather
//! than by coincidence, and every emitted row carries the catalog's own `family_id`/`tid`.
//!
//! # Contract (fail loudly, never silently drop)
//!
//! Every check below is an ERROR, not a filter. A silently dropped copy is exactly the defect class this
//! audit keeps finding (the `unwrap_or_default()` famCN silent degradation; the subset-BAM traps), and it
//! would be undetectable here — a copy quietly missing from O2's roster looks identical to a copy O2
//! legitimately could not assign reads to.
//!
//! * the TSV must carry the `gw_family_catalog` header and every required column (parsed BY NAME, so a
//!   future appended column cannot shift the meaning of a field);
//! * a family's copies must all sit on ONE chromosome — `copy_assign` is region-scoped, so a cross-chrom
//!   family (RABL2's 5 contigs) is STRUCTURALLY unassignable here and must say so rather than be truncated
//!   to whichever copies happened to fall in the region;
//! * the exon blocks must be well formed and reconstruct the copy's own `start`/`end` and `n_exon`;
//! * with `--copies-fa`, EVERY supplied copy must have a FASTA record whose header coordinates match its
//!   TSV row (the header is `>{family_id}|{copy_idx}|{chrom}:{start}-{end}|{strand}|nexon={n}`);
//! * without `--copies-fa`, the spliced sequence is rebuilt from the genome at the catalog's own exon
//!   coordinates via the SAME `build_spliced_seq` the catalog used, and the strand it derives from the
//!   junction motifs must agree with the strand the catalog recorded (a disagreement means the FASTA is
//!   not the assembly the catalog was built against).

use std::collections::BTreeMap;

use anyhow::{bail, Context, Result};

use super::denovo_pipeline::ColocatedFamily;
use super::family_detect::DenovoTranscript;
use crate::genome::GenomeIndex;

/// One `<out>.copies.tsv` row: a catalog COPY, with its catalog identity (`family_id`, `copy_idx`, `tid`)
/// preserved so the O2 output can be joined back to the O1 table it came from.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CatalogCopy {
    pub family_id: String,
    pub copy_idx: usize,
    pub tid: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub n_reads: u32,
    /// Half-open genomic exon blocks, ascending — the `exons` column, verbatim.
    pub exons: Vec<(u64, u64)>,
}

/// A catalog FAMILY: its rows grouped by `family_id`, in first-seen (file) order.
#[derive(Clone, Debug)]
pub struct CatalogFamily {
    pub family_id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub copies: Vec<CatalogCopy>,
}

/// A `<out>.copies.fa` record: the sequence plus the header fields it is keyed by, so the header can be
/// CHECKED against the TSV row rather than trusted.
#[derive(Clone, Debug)]
pub struct CatalogSeq {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub n_exon: usize,
    pub seq: Vec<u8>,
}

/// `(family_id, copy_idx)` → the catalog's own emitted sequence.
pub type SeqIndex = BTreeMap<(String, usize), CatalogSeq>;

/// Intron `(donor, acceptor)` chain implied by an ascending half-open exon block list. Inverse of
/// `gw_family_catalog::exon_blocks`.
pub fn introns_of(exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
    exons.windows(2).map(|w| (w[0].1, w[1].0)).collect()
}

fn parse_exons(s: &str) -> Result<Vec<(u64, u64)>> {
    let mut out = Vec::new();
    for blk in s.split(',') {
        let blk = blk.trim();
        if blk.is_empty() {
            continue;
        }
        let (a, b) = blk.split_once('-').with_context(|| format!("malformed exon block {blk:?}"))?;
        let a: u64 = a.parse().with_context(|| format!("malformed exon block {blk:?}"))?;
        let b: u64 = b.parse().with_context(|| format!("malformed exon block {blk:?}"))?;
        out.push((a, b));
    }
    Ok(out)
}

/// Parse a `gw_family_catalog` `<out>.copies.tsv`. Columns are located BY HEADER NAME.
pub fn parse_copies_tsv(text: &str) -> Result<Vec<CatalogCopy>> {
    let mut lines = text.lines();
    let header = lines.next().context("--families file is empty (expected a gw_family_catalog copies.tsv)")?;
    let cols: Vec<&str> = header.split('\t').collect();
    let idx = |name: &str| -> Result<usize> {
        cols.iter().position(|c| *c == name).with_context(|| {
            format!("--families: copies.tsv has no `{name}` column (header was: {header:?})")
        })
    };
    let (i_fam, i_ci, i_tid, i_chrom, i_start, i_end, i_nexon, i_strand, i_reads, i_exons) = (
        idx("family_id")?, idx("copy_idx")?, idx("tid")?, idx("chrom")?, idx("start")?, idx("end")?,
        idx("n_exon")?, idx("strand")?, idx("n_reads")?, idx("exons")?,
    );
    let need = cols.len();
    let mut out = Vec::new();
    for (ln, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < need {
            bail!("--families: copies.tsv line {} has {} fields, expected {need}", ln + 2, f.len());
        }
        let at = |i: usize| f[i];
        let strand_s = at(i_strand);
        let strand = match strand_s {
            "+" => '+',
            "-" => '-',
            other => bail!("--families: copies.tsv line {}: strand {other:?} is neither + nor -", ln + 2),
        };
        let exons = parse_exons(at(i_exons))
            .with_context(|| format!("--families: copies.tsv line {}", ln + 2))?;
        let n_exon: usize = at(i_nexon)
            .parse()
            .with_context(|| format!("--families: copies.tsv line {}: bad n_exon", ln + 2))?;
        let c = CatalogCopy {
            family_id: at(i_fam).to_string(),
            copy_idx: at(i_ci).parse().with_context(|| format!("line {}: bad copy_idx", ln + 2))?,
            tid: at(i_tid).to_string(),
            chrom: at(i_chrom).to_string(),
            start: at(i_start).parse().with_context(|| format!("line {}: bad start", ln + 2))?,
            end: at(i_end).parse().with_context(|| format!("line {}: bad end", ln + 2))?,
            strand,
            n_reads: at(i_reads).parse().with_context(|| format!("line {}: bad n_reads", ln + 2))?,
            exons,
        };
        // Structural checks on the row itself: an exon chain that does not reconstruct the copy's own
        // span/exon count means the row is not the object the catalog wrote.
        if c.exons.is_empty() {
            bail!("--families: {} copy {} has no exon blocks", c.family_id, c.copy_idx);
        }
        if c.exons.len() != n_exon {
            bail!(
                "--families: {} copy {}: n_exon={n_exon} but the exons column has {} blocks",
                c.family_id, c.copy_idx, c.exons.len()
            );
        }
        for w in c.exons.windows(2) {
            if w[0].1 > w[1].0 {
                bail!("--families: {} copy {}: exon blocks are not ascending/disjoint", c.family_id, c.copy_idx);
            }
        }
        if c.exons.iter().any(|&(s, e)| s >= e) {
            bail!("--families: {} copy {}: an exon block is empty or reversed", c.family_id, c.copy_idx);
        }
        if c.exons[0].0 != c.start || c.exons[c.exons.len() - 1].1 != c.end {
            bail!(
                "--families: {} copy {}: exon blocks span {}-{} but the row says {}-{}",
                c.family_id, c.copy_idx, c.exons[0].0, c.exons[c.exons.len() - 1].1, c.start, c.end
            );
        }
        out.push(c);
    }
    if out.is_empty() {
        bail!("--families: copies.tsv has a header but no copy rows");
    }
    Ok(out)
}

/// Parse a `gw_family_catalog` `<out>.copies.fa`, keyed by `(family_id, copy_idx)` from the header
/// `>{family_id}|{copy_idx}|{chrom}:{start}-{end}|{strand}|nexon={n}`.
pub fn parse_copies_fa(text: &str) -> Result<SeqIndex> {
    let mut out: SeqIndex = BTreeMap::new();
    let mut cur: Option<((String, usize), CatalogSeq)> = None;
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some((k, v)) = cur.take() {
                out.insert(k, v);
            }
            let parts: Vec<&str> = h.split('|').collect();
            if parts.len() < 5 {
                bail!("--copies-fa: malformed header {h:?} (expected fam|idx|chrom:start-end|strand|nexon=N)");
            }
            let copy_idx: usize =
                parts[1].parse().with_context(|| format!("--copies-fa: bad copy index in {h:?}"))?;
            let (chrom, span) = parts[2]
                .rsplit_once(':')
                .with_context(|| format!("--copies-fa: bad locus field in {h:?}"))?;
            let (s, e) = span
                .split_once('-')
                .with_context(|| format!("--copies-fa: bad span in {h:?}"))?;
            let strand = match parts[3] {
                "+" => '+',
                "-" => '-',
                other => bail!("--copies-fa: strand {other:?} in {h:?} is neither + nor -"),
            };
            let n_exon: usize = parts[4]
                .strip_prefix("nexon=")
                .with_context(|| format!("--copies-fa: missing nexon= in {h:?}"))?
                .parse()
                .with_context(|| format!("--copies-fa: bad nexon in {h:?}"))?;
            cur = Some((
                (parts[0].to_string(), copy_idx),
                CatalogSeq {
                    chrom: chrom.to_string(),
                    start: s.parse().with_context(|| format!("--copies-fa: bad start in {h:?}"))?,
                    end: e.parse().with_context(|| format!("--copies-fa: bad end in {h:?}"))?,
                    strand,
                    n_exon,
                    seq: Vec::new(),
                },
            ));
        } else if let Some((_, v)) = cur.as_mut() {
            v.seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        } else if !line.trim().is_empty() {
            bail!("--copies-fa: sequence line before any `>` header");
        }
    }
    if let Some((k, v)) = cur.take() {
        out.insert(k, v);
    }
    Ok(out)
}

/// Group parsed rows into families, in first-seen order, enforcing the SAME-CHROMOSOME contract.
///
/// A cross-chrom catalog family is not a `copy_assign` object: this binary is region-scoped, so honouring
/// only the copies that fall in the region would silently assign reads against a TRUNCATED roster — which
/// tightens the Bonferroni certificate over the wrong K and mislabels reads whose true copy is on another
/// contig. Refuse instead.
pub fn group_families(copies: Vec<CatalogCopy>) -> Result<Vec<CatalogFamily>> {
    let mut order: Vec<String> = Vec::new();
    let mut by_id: BTreeMap<String, Vec<CatalogCopy>> = BTreeMap::new();
    for c in copies {
        if !by_id.contains_key(&c.family_id) {
            order.push(c.family_id.clone());
        }
        by_id.entry(c.family_id.clone()).or_default().push(c);
    }
    let mut out = Vec::new();
    for fid in order {
        let mut cs = by_id.remove(&fid).expect("family id was recorded in `order`");
        let chroms: std::collections::BTreeSet<&str> = cs.iter().map(|c| c.chrom.as_str()).collect();
        if chroms.len() > 1 {
            bail!(
                "--families: {fid} spans {} chromosomes ({}) — copy_assign is REGION-scoped, so a \
                 cross-chrom family cannot be assigned here without silently truncating its copy set. \
                 Split the catalog or assign this family per-contig.",
                chroms.len(),
                chroms.into_iter().collect::<Vec<_>>().join(",")
            );
        }
        // Same ordering guarantee `colocated_families`/`colocated_from_copies` give the assignment step.
        cs.sort_by_key(|c| c.start);
        let chrom = cs[0].chrom.clone();
        let start = cs.iter().map(|c| c.start).min().unwrap_or(0);
        let end = cs.iter().map(|c| c.end).max().unwrap_or(0);
        out.push(CatalogFamily { family_id: fid, chrom, start, end, copies: cs });
    }
    Ok(out)
}

/// Where a copy's spliced sequence came from — reported so a run is never ambiguous about which substrate
/// its copy set was materialized from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeqSource {
    /// `--copies-fa`: the catalog's OWN emitted bytes. Agreement with O1 then holds by construction.
    CopiesFa,
    /// Rebuilt from `--fasta` at the catalog's exon coordinates via `build_spliced_seq` (the same function
    /// `assemble_gate` used to write the catalog), with the derived strand checked against the TSV.
    Genome,
}

/// Materialize one catalog family as the `ColocatedFamily` the assignment stage consumes, KEEPING the
/// catalog's `family_id` and per-copy `tid` so every emitted row joins back to `copies.tsv`.
///
/// `seqs` = `Some` ⟹ `--copies-fa` (exact catalog bytes, checked against the TSV row); `None` ⟹ rebuild
/// from `genome`. Every failure is an error — nothing is dropped.
pub fn to_colocated(
    fam: &CatalogFamily,
    seqs: Option<&SeqIndex>,
    genome: &GenomeIndex,
) -> Result<(ColocatedFamily, SeqSource)> {
    let mut copies: Vec<DenovoTranscript> = Vec::with_capacity(fam.copies.len());
    let source = if seqs.is_some() { SeqSource::CopiesFa } else { SeqSource::Genome };
    for c in &fam.copies {
        let introns = introns_of(&c.exons);
        let seq: Vec<u8> = match seqs {
            Some(ix) => {
                let rec = ix.get(&(c.family_id.clone(), c.copy_idx)).with_context(|| {
                    format!(
                        "--copies-fa has no record for {} copy {} ({}:{}-{}); the FASTA does not match the \
                         --families table",
                        c.family_id, c.copy_idx, c.chrom, c.start, c.end
                    )
                })?;
                if rec.chrom != c.chrom || rec.start != c.start || rec.end != c.end
                    || rec.strand != c.strand || rec.n_exon != c.exons.len()
                {
                    bail!(
                        "--copies-fa record for {} copy {} says {}:{}-{} {} nexon={} but copies.tsv says \
                         {}:{}-{} {} nexon={}",
                        c.family_id, c.copy_idx, rec.chrom, rec.start, rec.end, rec.strand, rec.n_exon,
                        c.chrom, c.start, c.end, c.strand, c.exons.len()
                    );
                }
                if rec.seq.is_empty() {
                    bail!("--copies-fa record for {} copy {} is empty", c.family_id, c.copy_idx);
                }
                rec.seq.clone()
            }
            None => {
                let (seq, derived) =
                    super::denovo_assemble::build_spliced_seq(genome, &c.chrom, c.start, c.end, &introns)
                        .with_context(|| {
                            format!(
                                "--families: could not rebuild {} copy {} ({}:{}-{}) from --fasta — the \
                                 junction motifs are not canonical in this assembly, or the contig/coords \
                                 are absent. Pass --copies-fa to use the catalog's own sequences.",
                                c.family_id, c.copy_idx, c.chrom, c.start, c.end
                            )
                        })?;
                if derived != c.strand {
                    bail!(
                        "--families: {} copy {} ({}:{}-{}): --fasta gives strand {derived} but copies.tsv \
                         says {} — the FASTA is not the assembly this catalog was built against",
                        c.family_id, c.copy_idx, c.chrom, c.start, c.end, c.strand
                    );
                }
                seq
            }
        };
        copies.push(DenovoTranscript {
            tid: c.tid.clone(),
            chrom: c.chrom.clone(),
            start: c.start,
            end: c.end,
            n_reads: c.n_reads,
            strand: c.strand,
            introns,
            seq,
            ..Default::default()
        });
    }
    Ok((
        ColocatedFamily {
            family_id: fam.family_id.clone(),
            chrom: fam.chrom.clone(),
            start: fam.start,
            end: fam.end,
            copies,
        },
        source,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const HDR: &str = "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons";

    fn row(fam: &str, ci: usize, chrom: &str, s: u64, e: u64, exons: &str, n_exon: usize) -> String {
        format!("{fam}\t{ci}\tDN_{chrom}_{s}_{n_exon}\t{chrom}\t{s}\t{e}\t{n_exon}\t+\t7\t{exons}")
    }

    #[test]
    fn parses_a_catalog_row_and_keeps_the_catalog_identity() {
        let t = format!("{HDR}\n{}\n", row("GWFAM3", 1, "c1", 100, 400, "100-200,300-400", 2));
        let cs = parse_copies_tsv(&t).unwrap();
        assert_eq!(cs.len(), 1);
        assert_eq!(cs[0].family_id, "GWFAM3");
        assert_eq!(cs[0].copy_idx, 1);
        assert_eq!(cs[0].tid, "DN_c1_100_2");
        assert_eq!(cs[0].exons, vec![(100, 200), (300, 400)]);
        assert_eq!(introns_of(&cs[0].exons), vec![(200, 300)]);
    }

    #[test]
    fn columns_are_located_by_name_not_position() {
        // an EXTRA leading column must not shift the meaning of any field
        let hdr = format!("extra\t{HDR}");
        let t = format!("{hdr}\nZ\t{}\n", row("GWFAM0", 0, "c1", 0, 60, "0-60", 1));
        let cs = parse_copies_tsv(&t).unwrap();
        assert_eq!(cs[0].chrom, "c1");
        assert_eq!(cs[0].end, 60);
    }

    #[test]
    fn a_missing_required_column_is_an_error_not_a_default() {
        let t = "family_id\tcopy_idx\ttid\nGWFAM0\t0\tDN_x\n";
        let e = parse_copies_tsv(t).unwrap_err().to_string();
        assert!(e.contains("chrom"), "{e}");
    }

    #[test]
    fn exon_blocks_must_reconstruct_the_rows_own_span() {
        let t = format!("{HDR}\n{}\n", row("GWFAM0", 0, "c1", 100, 400, "100-200,300-390", 2));
        let e = parse_copies_tsv(&t).unwrap_err().to_string();
        assert!(e.contains("exon blocks span"), "{e}");
    }

    #[test]
    fn n_exon_must_agree_with_the_exon_column() {
        let t = format!("{HDR}\n{}\n", row("GWFAM0", 0, "c1", 100, 400, "100-200,300-400", 3));
        let e = parse_copies_tsv(&t).unwrap_err().to_string();
        assert!(e.contains("n_exon=3"), "{e}");
    }

    #[test]
    fn a_header_only_file_is_an_error_rather_than_an_empty_roster() {
        let e = parse_copies_tsv(&format!("{HDR}\n")).unwrap_err().to_string();
        assert!(e.contains("no copy rows"), "{e}");
    }

    #[test]
    fn families_group_in_first_seen_order_and_sort_copies_by_start() {
        let t = format!(
            "{HDR}\n{}\n{}\n{}\n",
            row("GWFAM9", 0, "c1", 500, 600, "500-600", 1),
            row("GWFAM9", 1, "c1", 100, 200, "100-200", 1),
            row("GWFAM2", 0, "c1", 900, 950, "900-950", 1),
        );
        let fams = group_families(parse_copies_tsv(&t).unwrap()).unwrap();
        assert_eq!(fams.iter().map(|f| f.family_id.as_str()).collect::<Vec<_>>(), vec!["GWFAM9", "GWFAM2"]);
        assert_eq!(fams[0].copies.iter().map(|c| c.start).collect::<Vec<_>>(), vec![100, 500]);
        assert_eq!((fams[0].start, fams[0].end), (100, 600));
    }

    #[test]
    fn a_cross_chrom_family_is_refused_loudly_never_truncated() {
        let t = format!(
            "{HDR}\n{}\n{}\n",
            row("GWFAM0", 0, "c1", 0, 60, "0-60", 1),
            row("GWFAM0", 1, "c2", 0, 60, "0-60", 1),
        );
        let e = group_families(parse_copies_tsv(&t).unwrap()).unwrap_err().to_string();
        assert!(e.contains("spans 2 chromosomes"), "{e}");
        assert!(e.contains("truncating"), "{e}");
    }

    #[test]
    fn copies_fa_is_keyed_by_family_and_copy_index_with_checkable_coordinates() {
        let fa = ">GWFAM0|1|c1:100-400|-|nexon=2\nacgt\nACGT\n>GWFAM1|0|c2:0-50|+|nexon=1\nTTTT\n";
        let ix = parse_copies_fa(fa).unwrap();
        assert_eq!(ix.len(), 2);
        let r = &ix[&("GWFAM0".to_string(), 1)];
        assert_eq!(r.seq, b"ACGTACGT".to_vec(), "multi-line records concatenate, uppercased");
        assert_eq!((r.chrom.as_str(), r.start, r.end, r.strand, r.n_exon), ("c1", 100, 400, '-', 2));
    }

    #[test]
    fn a_malformed_copies_fa_header_is_an_error() {
        let e = parse_copies_fa(">GWFAM0|0|c1:0-10\nAC\n").unwrap_err().to_string();
        assert!(e.contains("malformed header"), "{e}");
    }
}
