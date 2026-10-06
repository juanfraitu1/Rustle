//! **STATUS:** SHIPPED-DEFAULT
//!
//! The `--from-gtf` family stage as a library: loci of an assembled GTF, their all-vs-all
//! (`minimap2 -x asm20`), and the copy table (`write_locus_rep_copies`) consumed by
//! `copy_assign --families`. Extracted verbatim from `src/bin/mcl_families.rs` (2026-10-04);
//! `mcl_families` itself re-imports every moved function. ⚠ Every function here is
//! byte-identity-critical: families-stage products are cmp-checked across builds.

use anyhow::{Context, Result};
use crate::genome::GenomeIndex;
use crate::family::annotation_families::{Cluster, CoreRecord, CoreStatus, GeneKey};
use crate::family::denovo_assemble::longest_orf;
use std::collections::BTreeMap;
use crate::types::DetHashSet;
use std::io::Write;

/// Header of `<out>.copies.tsv` (`--from-gtf --emit-units`). Columns 1-11 are `gw_family_catalog`'s `copies.tsv`
/// header, in its order, so `copy_assign --families` (parsed by name), `bench/sim.py copies` (positional 1-9) and
/// every reader of the legacy catalog read it unchanged; the rest are appended.
pub const COPIES_HEADER: &str = "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tmax_family_identity\
     \tsource\tgene_id\tcore_hull\tsd_depth\tcore_bp\trep_frac\tmember_status\tlocus_start\tlocus_end";

/// `exonic_blocks` of a GFF (`--from-gtf`: the `<out>.loci.gff3` this module writes): per gene key
/// (⚠ GFF 1-based, verbatim) the unioned, sorted exon blocks (0-based half-open). With
/// `exonless_span`, a gene with no exon children is one block spanning itself.
pub fn exonic_blocks(gff: &str, exonless_span: bool) -> Result<BTreeMap<GeneKey, Vec<(u64, u64)>>> {
    let text = std::fs::read_to_string(gff).with_context(|| format!("reading {gff}"))?;
    let attr = |a: &str, k: &str| -> Option<String> {
        a.split(';').find_map(|f| f.strip_prefix(k).map(|v| v.to_string()))
    };
    let mut span_of_name: BTreeMap<String, GeneKey> = BTreeMap::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || !matches!(f[2], "gene" | "pseudogene") {
            continue;
        }
        let (Ok(s), Ok(e)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else { continue };
        if let Some(n) = attr(f[8], "Name=") {
            span_of_name.insert(n, (f[0].to_string(), s, e)); // ⚠ GFF 1-based, verbatim
        }
    }
    let mut blocks: BTreeMap<GeneKey, Vec<(u64, u64)>> = BTreeMap::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || f[2] != "exon" {
            continue;
        }
        let Some(name) = attr(f[8], "gene=") else { continue };
        let Some(g) = span_of_name.get(&name) else { continue };
        let (Ok(s), Ok(e)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else { continue };
        if f[0] == g.0 {
            blocks.entry(g.clone()).or_default().push((s - 1, e));
        }
    }
    if exonless_span {
        // a record with no exon children is one exon: its own span (0-based half-open)
        for g in span_of_name.values() {
            blocks.entry(g.clone()).or_insert_with(|| vec![(g.1 - 1, g.2)]);
        }
    }
    let mut out: BTreeMap<GeneKey, Vec<(u64, u64)>> = BTreeMap::new();
    for (g, mut v) in blocks {
        v.sort_unstable();
        let (mut merged, mut cur) = (Vec::new(), v[0]);
        for &(s, e) in &v[1..] {
            if s <= cur.1 {
                cur.1 = cur.1.max(e);
            } else {
                merged.push(cur);
                cur = (s, e);
            }
        }
        merged.push(cur);
        out.insert(g, merged);
    }
    Ok(out)
}

/// Exon-union LENGTH per gene, derived from the same merge as [`exonic_blocks`] so the numerator and the
/// denominator can never disagree about what an exon is.
pub fn lengths_from_blocks(b: &BTreeMap<GeneKey, Vec<(u64, u64)>>) -> BTreeMap<GeneKey, u64> {
    b.iter()
        .map(|(g, v)| (g.clone(), v.iter().map(|&(s, e)| e - s).sum::<u64>().max(1)))
        .collect()
}

/// GFF gene strand per gene key (for `--emit-units` fallbacks and a tie-break when reads carry no strand).
pub fn gene_strands(gff: &str) -> Result<BTreeMap<GeneKey, char>> {
    let text = std::fs::read_to_string(gff).with_context(|| format!("reading {gff}"))?;
    let mut m = BTreeMap::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || !matches!(f[2], "gene" | "pseudogene") {
            continue;
        }
        let (Ok(s), Ok(e)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else { continue };
        m.insert((f[0].to_string(), s, e), f[6].chars().next().unwrap_or('+'));
    }
    Ok(m)
}

/// ⭐ L2: a locus never contains another catalog unit — of ANY family (§6fm). `units[i] = (chain start, chain end,
/// extent)` on ONE contig; the extent is clipped to the nearest chain ends of the other units that lie entirely outside its own
/// span (units overlapping the span — nested, interleaved — do not clip). Without this, in a tandem array a
/// single read-through molecule extends copy 1's locus over copy 2's chain, a copy-2 read then aligns identically
/// inside both targets and ties with zero decisive columns (MCL4: 226 assigned → tied, arm A0 of PREREG L1/L2).
pub fn clip_extents_to_neighbours(units: &[(u64, u64, (u64, u64))]) -> Vec<(u64, u64)> {
    units
        .iter()
        .enumerate()
        .map(|(i, &(s, e, (a, b)))| {
            let (mut lo, mut hi) = (a.min(s), b.max(e));
            for (j, &(s2, e2, _)) in units.iter().enumerate() {
                if j == i {
                    continue;
                }
                if e2 <= s {
                    lo = lo.max(e2);
                } else if s2 >= e {
                    hi = hi.min(s2);
                }
            }
            (lo, hi)
        })
        .collect()
}

/// Interspersed-repeat fraction of an exon chain (`None` when the contig has no RepeatMasker interval at all).
pub fn rep_frac_in(rmsk: &BTreeMap<String, Vec<(u64, u64)>>, chrom: &str, exons: &[(u64, u64)]) -> Option<f64> {
    let v = rmsk.get(chrom)?;
    let (mut tot, mut inter) = (0u64, 0u64);
    for &(s, e) in exons {
        tot += e - s;
        let i = v.partition_point(|x| x.1 <= s);
        for &(a, b) in &v[i..] {
            if a >= e {
                break;
            }
            inter += b.min(e).saturating_sub(a.max(s));
        }
    }
    Some(inter as f64 / tot.max(1) as f64)
}

/// Counts of one `write_locus_rep_copies` call (the params certificate and the log).
#[derive(Debug, Default, Clone, PartialEq)]
pub struct RepCopyStats {
    /// rows written to `copies.tsv`
    pub copies: usize,
    /// families with >= 1 copy / with >= 2 copies
    pub families: usize,
    pub multi_copy_families: usize,
    /// members folded into an exon-overlapping copy of the same family (`copies.merged.tsv`)
    pub merged: usize,
    /// members skipped because the core rule dropped them and `--no-units-include-dropped` was given
    pub skipped_dropped: usize,
    pub dropped_emitted: usize,
    pub noncoding: usize,
    /// copies whose representative has `reads 0` (or no `reads` attribute)
    pub unexpressed: usize,
    /// representatives with strand `.`, written as `+` (the copies contract has no unstranded copy)
    pub unstranded: usize,
    /// loci sharing one `(chrom, start, end)` with another `gene_id` (one graph node; the representative with more
    /// reads is the copy)
    pub key_collisions: usize,
    /// representatives whose exons overlapped or abutted and were coalesced into one block
    pub coalesced: usize,
}

/// Genome bases at `exons` (0-based half-open, ascending), concatenated, reverse-complemented on `-`: a copy's
/// spliced sequence in transcription orientation, as `gw_family_catalog` writes `copies.fa`.
pub fn spliced_exon_sum(genome: &GenomeIndex, chrom: &str, exons: &[(u64, u64)], strand: char) -> Result<Vec<u8>> {
    let mut seq: Vec<u8> = Vec::new();
    for &(s, e) in exons {
        let part = genome
            .fetch_sequence(chrom, s, e)
            .with_context(|| format!("copies: {chrom}:{s}-{e} is not in --fasta"))?;
        seq.extend_from_slice(&part);
    }
    if strand == '-' {
        seq = crate::family::seq_utils::revcomp_keep_case(&seq);
    }
    Ok(seq)
}

/// ⭐ `--from-gtf --emit-units` (user decision 2026-09-25 16:00: ONE default de novo family definition, and copy
/// assignment consumes the SAME families): write the families as a COPY TABLE in the `copy_assign --families /
/// --copies-fa` contract. One copy per member locus of every reported cluster (the rows of `clusters.tsv`,
/// `family_id` = its `cluster_id`); the copy IS the locus representative — its exons (the positional exon sum),
/// its spliced sequence, its `transcript_id` (`tid`, which joins back to the assembled GTF), its `reads`
/// (`n_reads`) and its strand. No BAM is read: the de novo loci already are read-derived.
///
/// `max_family_identity` = the identity of the best families-stage edge (genomic-span `-x asm20` alignment, the
/// one the admission rule scored) from this locus to another member, `NA` when it has no direct edge. ⚠ The legacy
/// catalog's column of that name is an exon-sum alignment identity: same role, different alignment.
/// `locus_start`/`locus_end` = the de novo locus span (all its transcripts), clipped at the exon-chain ends of the
/// neighbouring copies on the contig (the L2 rule, `clip_extents_to_neighbours`).
/// Also applied as for read-chain units: `--no-units-include-dropped`, `--coding-core`, and the §6fb merge of
/// copies of one family that share exon bases (a no-op after `--fold-within-clusters`, kept as the guarantee that
/// no base belongs to two copies of one family). Writes `<out>.copies.tsv/.fa/.regions/.merged.tsv`.
#[allow(clippy::too_many_arguments)]
pub fn write_locus_rep_copies(
    out: &str,
    fasta: &str,
    clusters: &[Cluster],
    g: &crate::family::annotation_families::HomologyGraph,
    core_records: &[Vec<crate::family::annotation_families::CoreRecord>],
    loci: &[GtfLocus],
    rmsk: Option<&BTreeMap<String, Vec<(u64, u64)>>>,
    include_dropped: bool,
    merge_overlapping: bool,
    coding_core: bool,
) -> Result<RepCopyStats> {
    struct RepCopy {
        member: GeneKey,
        gene_id: String,
        tid: String,
        reads: u64,
        strand: char,
        exons: Vec<(u64, u64)>,
        seq: Vec<u8>,
        ident: Option<f64>,
        hull_col: String,
        sd_depth: String,
        core_bp: String,
        status: &'static str,
        orf: usize,
        locus: (u64, u64),
    }
    let mut st = RepCopyStats::default();
    // representative per graph node (`loci.gff3` gene line = the node key, GFF 1-based)
    let mut rep_of: BTreeMap<GeneKey, &GtfLocus> = BTreeMap::new();
    for l in loci {
        let k: GeneKey = (l.chrom.clone(), l.start, l.end);
        match rep_of.get(&k) {
            Some(prev) => {
                st.key_collisions += 1;
                if l.rep_reads > prev.rep_reads {
                    rep_of.insert(k, l);
                }
            }
            None => {
                rep_of.insert(k, l);
            }
        }
    }
    let contigs: DetHashSet<String> =
        clusters.iter().flat_map(|c| c.members.iter().map(|m| m.0.clone())).collect();
    // ⚠ `from_fasta_contigs` with an EMPTY set loads the whole genome: no family, no genome
    let genome = if contigs.is_empty() {
        GenomeIndex::empty()
    } else {
        GenomeIndex::from_fasta_contigs(fasta, &contigs)?
    };
    let node_idx: BTreeMap<&GeneKey, usize> = g.genes.iter().enumerate().map(|(k, gk)| (gk, k)).collect();
    let mut staged: Vec<(String, Vec<RepCopy>, Vec<Option<usize>>)> = Vec::new();
    for (i, c) in clusters.iter().enumerate() {
        let fid = format!("MCL{i}");
        let mut pending: Vec<RepCopy> = Vec::new();
        for (mi, m) in c.members.iter().enumerate() {
            let rec = core_records.get(i).and_then(|v| v.get(mi));
            let status: &'static str = match rec.map(|r| r.status) {
                Some(crate::family::annotation_families::CoreStatus::Dropped) => "dropped",
                Some(crate::family::annotation_families::CoreStatus::KeptTrimmed) => "kept_trimmed",
                Some(_) => "kept_full",
                None => "ungated",
            };
            if status == "dropped" && !include_dropped {
                st.skipped_dropped += 1;
                continue;
            }
            let l = rep_of.get(m).with_context(|| {
                format!("copies: family member {}:{}-{} is not a locus of the --from-gtf GTF (a failed join)", m.0, m.1, m.2)
            })?;
            // the representative's exons, 0-based half-open; overlapping or abutting blocks coalesced
            let mut exons: Vec<(u64, u64)> = Vec::with_capacity(l.rep_exons.len());
            let mut coalesced = false;
            for &(_, a, b) in &l.rep_exons {
                let (s, e) = (a.saturating_sub(1), b);
                match exons.last_mut() {
                    Some(p) if s <= p.1 => {
                        p.1 = p.1.max(e);
                        coalesced = true;
                    }
                    _ => exons.push((s, e)),
                }
            }
            anyhow::ensure!(!exons.is_empty(), "copies: locus {} has a representative without exons", l.gene_id);
            if coalesced {
                st.coalesced += 1;
            }
            let strand = match l.strand.as_str() {
                "+" => '+',
                "-" => '-',
                _ => {
                    st.unstranded += 1;
                    '+'
                }
            };
            let seq = spliced_exon_sum(&genome, &m.0, &exons, strand)?;
            let ident = node_idx.get(m).and_then(|&a| {
                c.members
                    .iter()
                    .filter_map(|o| node_idx.get(o).copied())
                    .filter(|&b| b != a)
                    .filter_map(|b| g.idents.get(&(a.min(b), a.max(b))).copied())
                    .reduce(f64::max)
            });
            let hull_col = match rec.and_then(|r| r.hull) {
                Some((a, b)) => format!("{}-{}", a.saturating_sub(1), b),
                None => "NA".to_string(),
            };
            let (sd_depth, core_bp) =
                rec.map(|r| (r.max_depth.to_string(), r.core_bp.to_string())).unwrap_or_else(|| ("NA".into(), "NA".into()));
            let (us, ue) = (exons[0].0, exons.last().unwrap().1);
            if status == "dropped" {
                st.dropped_emitted += 1;
            }
            pending.push(RepCopy {
                member: m.clone(),
                gene_id: l.gene_id.clone(),
                tid: l.rep.clone(),
                reads: l.rep_reads,
                strand,
                orf: if coding_core { longest_orf(&seq) } else { 0 },
                exons,
                seq,
                ident,
                hull_col,
                sd_depth,
                core_bp,
                status,
                locus: (m.1.saturating_sub(1).min(us), m.2.max(ue)),
            });
        }
        // copies in genomic order within the family (the catalog's order)
        pending.sort_by(|a, b| {
            (a.member.0.as_str(), a.exons[0].0, a.exons.last().unwrap().1, a.tid.as_str())
                .cmp(&(b.member.0.as_str(), b.exons[0].0, b.exons.last().unwrap().1, b.tid.as_str()))
        });
        // CODING CORE: the family's best ORF is the reference (see `--coding-core`)
        if coding_core {
            let best_orf = pending.iter().filter(|u| u.status != "dropped").map(|u| u.orf).max().unwrap_or(0);
            if best_orf > 0 {
                for u in pending.iter_mut() {
                    if u.status != "dropped" && u.orf * 2 < best_orf {
                        u.status = "noncoding";
                        st.noncoding += 1;
                    }
                }
            }
        }
        // §6fb: copies of one family that share exon bases are one copy (kept before dropped, then the longest
        // exon union represents them)
        let n = pending.len();
        let mut merged_into: Vec<Option<usize>> = vec![None; n];
        if merge_overlapping {
            let mut parent: Vec<usize> = (0..n).collect();
            fn find(p: &mut [usize], mut x: usize) -> usize {
                while p[x] != x {
                    p[x] = p[p[x]];
                    x = p[x];
                }
                x
            }
            for a in 0..n {
                for b in (a + 1)..n {
                    if pending[a].member.0 != pending[b].member.0 {
                        continue;
                    }
                    let share = pending[a].exons.iter().any(|&(s1, e1)| pending[b].exons.iter().any(|&(s2, e2)| s1 < e2 && s2 < e1));
                    if share {
                        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                        if ra != rb {
                            parent[ra.max(rb)] = ra.min(rb);
                        }
                    }
                }
            }
            let rank = |u: &RepCopy| (u.status != "dropped", u.exons.iter().map(|(s, e)| e - s).sum::<u64>());
            let mut rep_of_root: BTreeMap<usize, usize> = BTreeMap::new();
            for k in 0..n {
                let r = find(&mut parent, k);
                let e = rep_of_root.entry(r).or_insert(k);
                if rank(&pending[k]) > rank(&pending[*e]) {
                    *e = k;
                }
            }
            for k in 0..n {
                let rep = rep_of_root[&find(&mut parent, k)];
                if rep != k {
                    merged_into[k] = Some(rep);
                }
            }
        }
        staged.push((fid, pending, merged_into));
    }
    // L2: every copy's locus extent stops at the exon-chain ends of the copies (of any family) around it
    let mut clipped: Vec<Vec<(u64, u64)>> = staged.iter().map(|(_, p, _)| p.iter().map(|u| u.locus).collect()).collect();
    {
        let mut by_ctg: BTreeMap<&str, Vec<(usize, usize)>> = BTreeMap::new();
        for (fi, (_, pending, merged_into)) in staged.iter().enumerate() {
            for (k, u) in pending.iter().enumerate() {
                if merged_into[k].is_none() {
                    by_ctg.entry(u.member.0.as_str()).or_default().push((fi, k));
                }
            }
        }
        for ks in by_ctg.values() {
            let spans: Vec<(u64, u64, (u64, u64))> = ks
                .iter()
                .map(|&(fi, k)| {
                    let u = &staged[fi].1[k];
                    (u.exons[0].0, u.exons.last().unwrap().1, u.locus)
                })
                .collect();
            for (&(fi, k), c) in ks.iter().zip(clip_extents_to_neighbours(&spans)) {
                clipped[fi][k] = c;
            }
        }
    }
    let mut ct = std::io::BufWriter::new(std::fs::File::create(format!("{out}.copies.tsv"))?);
    let mut cf = std::io::BufWriter::new(std::fs::File::create(format!("{out}.copies.fa"))?);
    let mut cr = std::io::BufWriter::new(std::fs::File::create(format!("{out}.copies.regions"))?);
    let mut cm = std::io::BufWriter::new(std::fs::File::create(format!("{out}.copies.merged.tsv"))?);
    writeln!(ct, "{COPIES_HEADER}")?;
    writeln!(cm, "family_id\tmerged_member\tmerged_tid\tinto_member\tinto_tid")?;
    for (fi, (fid, pending, merged_into)) in staged.iter().enumerate() {
        let mut idx = 0usize;
        let mut hulls: BTreeMap<&str, (u64, u64)> = BTreeMap::new();
        for (k, u) in pending.iter().enumerate() {
            if let Some(rep) = merged_into[k] {
                let r = &pending[rep];
                writeln!(
                    cm,
                    "{fid}\t{}:{}-{}\t{}\t{}:{}-{}\t{}",
                    u.member.0, u.member.1, u.member.2, u.tid, r.member.0, r.member.1, r.member.2, r.tid
                )?;
                st.merged += 1;
                continue;
            }
            let chrom = u.member.0.as_str();
            let (us, ue) = (u.exons[0].0, u.exons.last().unwrap().1);
            let rep_col = rmsk
                .and_then(|r| rep_frac_in(r, chrom, &u.exons))
                .map(|v| format!("{v:.3}"))
                .unwrap_or_else(|| "NA".into());
            writeln!(
                ct,
                "{fid}\t{idx}\t{}\t{chrom}\t{us}\t{ue}\t{}\t{}\t{}\t{}\t{}\tlocus_rep\t{}\t{}\t{}\t{}\t{rep_col}\t{}\t{}\t{}",
                u.tid,
                u.exons.len(),
                u.strand,
                u.reads,
                u.exons.iter().map(|(s, e)| format!("{s}-{e}")).collect::<Vec<_>>().join(","),
                u.ident.map(|v| format!("{v:.6}")).unwrap_or_else(|| "NA".into()),
                u.gene_id,
                u.hull_col,
                u.sd_depth,
                u.core_bp,
                u.status,
                clipped[fi][k].0,
                clipped[fi][k].1
            )?;
            writeln!(cf, ">{fid}|{idx}|{chrom}:{us}-{ue}|{}|nexon={}", u.strand, u.exons.len())?;
            cf.write_all(&u.seq)?;
            writeln!(cf)?;
            let h = hulls.entry(chrom).or_insert((us, ue));
            h.0 = h.0.min(us);
            h.1 = h.1.max(ue);
            if u.reads == 0 {
                st.unexpressed += 1;
            }
            st.copies += 1;
            idx += 1;
        }
        for (ctg, (a, b)) in hulls {
            writeln!(cr, "{fid}\t{ctg}:{}-{}", a.saturating_sub(5_000).max(1), b + 5_000)?;
        }
        if idx >= 1 {
            st.families += 1;
        }
        if idx >= 2 {
            st.multi_copy_families += 1;
        }
    }
    for w in [&mut ct, &mut cf, &mut cr, &mut cm] {
        w.flush()?;
    }
    Ok(st)
}

/// One de novo locus of an assembled GTF (`--from-gtf`): a `gene_id` group, its span (GFF 1-based, min/max over the
/// exons of all its transcripts) and its REPRESENTATIVE — the transcript with the most `reads`, ties to the longer
/// span, then to the lexicographically last `transcript_id` (the order `max_by_key` has always resolved ties in).
/// The representative's exons are the locus's exons in `loci.gff3` and, with `--emit-units`, the copy's exons in
/// `copies.tsv`: the "positional exon sum" (read-derived coordinates, genome bases).
#[derive(Clone, Debug, PartialEq)]
pub struct GtfLocus {
    pub gene_id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub rep: String,
    pub rep_reads: u64,
    /// The representative's strand column, verbatim (`.` when the transcript line had none).
    pub strand: String,
    /// The representative's exons, GFF 1-based closed, sorted by start (stable).
    pub rep_exons: Vec<(String, u64, u64)>,
}

/// Which transcript of a `--from-gtf` locus is its representative (`--representative`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Representative {
    /// the transcript with the most `reads`; ties to the longer span, then the last `transcript_id` (the shipped rule)
    MostReads,
    /// the transcript with the most junctions; ties to the most `reads`, then the longer span, then the last `transcript_id`
    MostJunctions,
}

const MIN_JUNCTION_GAP: u64 = 50;

/// Junctions of one transcript, for `--representative most-junctions`: the gaps of at least [`MIN_JUNCTION_GAP`] bp between
/// consecutive exons in coordinate order.
pub fn junction_count(exons: &[(String, u64, u64)]) -> usize {
    let mut iv: Vec<(u64, u64)> = exons.iter().map(|x| (x.1, x.2)).collect();
    iv.sort_unstable();
    iv.windows(2).filter(|w| w[1].0.saturating_sub(w[0].1 + 1) >= MIN_JUNCTION_GAP).count()
}

pub fn gtf_loci<R: std::io::BufRead>(reader: R, rule: Representative) -> Result<Vec<GtfLocus>> {
    use std::collections::{BTreeMap, HashMap, HashSet};
    fn attr<'a>(s: &'a str, key: &str) -> Option<&'a str> {
        let pat = format!("{key} \"");
        let i = s.find(&pat)? + pat.len();
        let j = s[i..].find('"')? + i;
        Some(&s[i..j])
    }
    let mut exons: HashMap<String, Vec<(String, u64, u64)>> = HashMap::new();
    let mut gene_of: HashMap<String, String> = HashMap::new();
    let mut strand: HashMap<String, String> = HashMap::new();
    let mut reads: HashMap<String, u64> = HashMap::new();
    let mut gene_order: Vec<String> = Vec::new();
    // genes already in `gene_order` (was a scan of every `gene_of` value per transcript line: O(T^2) on a
    // whole-genome GTF); same first-appearance order for any GTF whose transcript_ids do not switch gene
    let mut seen_genes: HashSet<String> = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < 9 {
            continue;
        }
        let Some(t) = attr(r[8], "transcript_id") else { continue };
        if r[2] == "transcript" {
            let g = attr(r[8], "gene_id").unwrap_or(t).to_string();
            if seen_genes.insert(g.clone()) {
                gene_order.push(g.clone());
            }
            gene_of.insert(t.to_string(), g);
            strand.insert(t.to_string(), r[6].to_string());
            reads.insert(t.to_string(), attr(r[8], "reads").and_then(|v| v.parse().ok()).unwrap_or(0));
        } else if r[2] == "exon" {
            exons.entry(t.to_string()).or_default().push((r[0].to_string(), r[3].parse()?, r[4].parse()?));
        }
    }
    let mut txs_of: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for (t, g) in &gene_of {
        txs_of.entry(g.clone()).or_default().push(t.clone());
    }
    let mut out = Vec::new();
    for g in &gene_order {
        let Some(ts) = txs_of.get(g) else { continue };
        let all: Vec<&(String, u64, u64)> = ts.iter().flat_map(|t| exons.get(t).into_iter().flatten()).collect();
        if all.is_empty() {
            continue;
        }
        let chrom = all[0].0.clone();
        let (s, e) = (all.iter().map(|x| x.1).min().unwrap(), all.iter().map(|x| x.2).max().unwrap());
        let span_of = |t: &String| exons.get(t).map(|v| v.iter().map(|x| x.2).max().unwrap() - v.iter().map(|x| x.1).min().unwrap()).unwrap_or(0);
        let mut sorted_ts = ts.clone();
        sorted_ts.sort();
        // `max_by_key` keeps the LAST maximum: over the sorted ids, a full tie goes to the last `transcript_id` under both
        // rules (`most-junctions` only puts the junction count in front of the `most-reads` key)
        let rep = match rule {
            Representative::MostReads => sorted_ts.iter().max_by_key(|t| (reads.get(*t).copied().unwrap_or(0), span_of(t))),
            Representative::MostJunctions => sorted_ts.iter().max_by_key(|t| {
                (exons.get(*t).map_or(0, |v| junction_count(v)), reads.get(*t).copied().unwrap_or(0), span_of(t))
            }),
        }
        .unwrap()
        .clone();
        let st = strand.get(&rep).cloned().unwrap_or_else(|| ".".into());
        let mut ex = exons.get(&rep).cloned().unwrap_or_default();
        ex.sort_by_key(|x| x.1);
        out.push(GtfLocus {
            gene_id: g.clone(),
            chrom,
            start: s,
            end: e,
            rep_reads: reads.get(&rep).copied().unwrap_or(0),
            rep,
            strand: st,
            rep_exons: ex,
        });
    }
    Ok(out)
}

/// The loci of `--from-gtf` as `--emit-relations` reads them: the graph's own keys and representatives.
pub fn relation_loci(loci: &[GtfLocus]) -> Vec<crate::family::fam_from_gtf::family_relations::LocusIn> {
    loci.iter()
        .map(|l| crate::family::fam_from_gtf::family_relations::LocusIn {
            gene_id: l.gene_id.clone(),
            key: (l.chrom.clone(), l.start as i64, l.end as i64),
            rep: l.rep.clone(),
            rep_reads: l.rep_reads as i64,
        })
        .collect()
}

/// Write `<out>.loci.fa`: one record `>CONTIG:START-END` + the genome's forward strand over that span, per locus span.
/// With `hash`, also the [`ContentHash`](crate::family::run_cache::ContentHash) of every byte written (the PAF
/// cache key), taken as the bytes go out so the file is never read back.
pub fn write_loci_fa(
    genome: &GenomeIndex,
    spans: &[(String, u64, u64)],
    fa_path: &str,
    fasta: &str,
    hash: bool,
) -> Result<Option<crate::family::run_cache::ContentHash>> {
    use crate::family::run_cache::HashingWriter;
    let mut fa = HashingWriter::new(std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(fa_path)?), hash);
    for (c, s, e) in spans {
        let seq = genome.fetch_sequence(c, s - 1, *e).with_context(|| format!("{c}:{s}-{e} not in {fasta}"))?;
        writeln!(fa, ">{c}:{s}-{e}")?;
        fa.write_all(&seq)?;
        writeln!(fa)?;
    }
    fa.flush()?;
    Ok(fa.hash)
}

/// The families PAF cache key (`paf/<fnv(key)>/key.tsv`): the minimap2 command line and build, and the content
/// hash + byte length of the loci FASTA it aligns. v2 (2026-09-28): the hash is the word-wise 128-bit
/// `ContentHash` taken while the FASTA is written; v1 was a byte-wise FNV-1a 64 of the FASTA read back from disk.
pub fn families_paf_key(cmd: &str, minimap2_version: &str, loci_fa: &crate::family::run_cache::ContentHash) -> String {
    format!(
        "rustle families paf v2\ncmd\t{cmd}\nminimap2\t{minimap2_version}\nquery_hash\tcontent128:{}\nquery_bytes\t{}\n",
        loci_fa.hex(),
        loci_fa.len()
    )
}

/// `--from-gtf`: the de novo locus set of an assembled GTF, as the family stage consumes it (see the flag doc).
/// Returns `(loci.gff3, loci.fa, loci.paf)` paths and the loci themselves (for the copy table).
/// `reader` lets the merged `copy_assign` mode parse an in-memory GTF (identical parser, no disk round-trip).
pub fn loci_from_gtf_reader<R: std::io::BufRead>(
    reader: R,
    fasta: &str,
    out: &str,
    threads: usize,
    rule: Representative,
) -> Result<(String, String, String, Vec<GtfLocus>)> {
    use std::collections::HashSet;
    let loci = gtf_loci(reader, rule)?;
    let gff3 = format!("{out}.loci.gff3");
    let fa_path = format!("{out}.loci.fa");
    let paf = format!("{out}.loci.paf");
    let mut g3 = std::fs::File::create(&gff3)?;
    writeln!(g3, "##gff-version 3")?;
    let mut spans: Vec<(String, u64, u64)> = Vec::new();
    for l in &loci {
        let (chrom, s, e, st, g) = (&l.chrom, l.start, l.end, &l.strand, &l.gene_id);
        writeln!(g3, "{chrom}\t.\tgene\t{s}\t{e}\t.\t{st}\t.\tID=gene-{g};Name={g}")?;
        for (_, a, b) in &l.rep_exons {
            writeln!(g3, "{chrom}\t.\texon\t{a}\t{b}\t.\t{st}\t.\tParent=gene-{g};gene={g}")?;
        }
        spans.push((chrom.clone(), s, e));
    }
    use crate::family::run_cache as rc;
    let root = rc::cache_root();
    let contigs: DetHashSet<String> = spans.iter().map(|x| x.0.clone()).collect();
    let genome = GenomeIndex::from_fasta_contigs(fasta, &contigs)?;
    let fa_hash = write_loci_fa(&genome, &spans, &fa_path, fasta, root.is_some())?;
    eprintln!("[mcl_families] --from-gtf: {} loci -> all-vs-all", spans.len());
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let mm_args: Vec<String> = ["-x", "asm20", "-c", "-X", "-N", "50", "-p", "0.1", "--secondary=yes", "-t"]
        .iter()
        .map(|s| s.to_string())
        .chain(std::iter::once(threads.to_string()))
        .collect();
    // PAF cache (`RUSTLE_CACHE_DIR`, see `crate::family::run_cache`): keyed by EVERY byte of the loci FASTA
    // (hashed while it was written above, never re-read), the command line and the minimap2 build; a hit hard-links
    // the cached PAF to `<out>.loci.paf` instead of re-aligning (and instead of copying it). The entry is pinned, so a
    // write through the link is a miss on the next run, never a stale replay.
    let paf_entry = root.zip(fa_hash).map(|(root, h)| {
        let key = families_paf_key(&format!("{mm2} {}", mm_args.join(" ")), &rc::minimap2_version(&mm2), &h);
        rc::Entry::new(&root, "paf", key).pinned()
    });
    let paf_path = std::path::Path::new(&paf);
    if let Some(e) = paf_entry.as_ref().filter(|e| e.is_hit()) {
        // a replay that fails (another run replacing the entry) falls through to running minimap2
        if let Ok(linked) = e.replay("out.paf", paf_path) {
            eprintln!(
                "[cache] all-vs-all PAF replayed from {} ({}; minimap2 skipped)",
                e.dir.display(),
                if linked { "hard link" } else { "copy" }
            );
            return Ok((gff3, fa_path, paf, loci));
        }
    }
    // never truncate in place: `<out>.loci.paf` may be a hard link to a cache entry from an earlier run
    match std::fs::remove_file(paf_path) {
        Ok(()) => {}
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {}
        Err(err) => return Err(err).with_context(|| format!("removing {paf}")),
    }
    let out_paf = std::fs::File::create(&paf)?;
    let status = std::process::Command::new(&mm2)
        .args(&mm_args)
        .arg(&fa_path)
        .arg(&fa_path)
        .stdout(out_paf)
        .stderr(std::process::Stdio::null())
        .status()
        .with_context(|| format!("running {mm2}"))?;
    anyhow::ensure!(status.success(), "minimap2 all-vs-all failed");
    if let Some(e) = paf_entry.as_ref() {
        let stored = e.staging().and_then(|st| {
            e.stage_link(&st, "out.paf", paf_path)?;
            e.commit(&st)
        });
        if let Err(err) = stored {
            eprintln!("[cache] could not store the PAF ({err:#}); continuing");
        }
    }
    Ok((gff3, fa_path, paf, loci))
}

/// File-reading convenience wrapper of [`loci_from_gtf_reader`].
pub fn loci_from_gtf(
    gtf: &str,
    fasta: &str,
    out: &str,
    threads: usize,
    rule: Representative,
) -> Result<(String, String, String, Vec<GtfLocus>)> {
    let f = std::fs::File::open(gtf).with_context(|| format!("opening {gtf}"))?;
    loci_from_gtf_reader(std::io::BufReader::new(f), fasta, out, threads, rule)
}

#[cfg(test)]
mod paf_cache_key_tests {
    use super::*;
    use crate::family::run_cache as rc;

    /// The families PAF key covers every byte of the loci FASTA (hashed as it is written, equal to a re-read of the
    /// file): the same loci FASTA written again (a new mtime, a later run) keeps the key, one changed base inside a
    /// locus span changes it at the same file size, and an entry committed under one key is a miss for the other.
    #[test]
    fn families_paf_key_changes_with_any_loci_fasta_byte_and_only_then() {
        let dir = std::env::temp_dir().join(format!("rustle_fam_paf_key_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let g = dir.join("g.fa");
        let spans = vec![("c1".to_string(), 3u64, 12u64), ("c2".to_string(), 1, 8)];
        let key_of = |genome_text: &str, name: &str| -> (String, Vec<u8>) {
            std::fs::write(&g, genome_text).unwrap();
            let genome = GenomeIndex::from_fasta(g.to_str().unwrap()).unwrap();
            let fa = dir.join(name);
            let h = write_loci_fa(&genome, &spans, fa.to_str().unwrap(), "g.fa", true).unwrap().unwrap();
            let bytes = std::fs::read(&fa).unwrap();
            assert_eq!(rc::ContentHash::of_file(&fa).unwrap().hex(), h.hex(), "in-stream hash == hash of the file");
            assert_eq!(h.len(), bytes.len() as u64);
            // without a cache nothing is hashed and the file is the same
            assert!(write_loci_fa(&genome, &spans, fa.to_str().unwrap(), "g.fa", false).unwrap().is_none());
            assert_eq!(std::fs::read(&fa).unwrap(), bytes);
            (families_paf_key("minimap2 -x asm20 -t 4", "2.28-r1209", &h), bytes)
        };
        let (k1, b1) = key_of(">c1\nACGTACGTACGTAC\n>c2\nGGGGCCCCAA\n", "a.loci.fa");
        assert_eq!(b1, b">c1:3-12\nGTACGTACGT\n>c2:1-8\nGGGGCCCC\n");
        std::thread::sleep(std::time::Duration::from_millis(20));
        let (k1b, b1b) = key_of(">c1\nACGTACGTACGTAC\n>c2\nGGGGCCCCAA\n", "a.loci.fa");
        assert_eq!((&k1, &b1), (&k1b, &b1b), "rewritten, same bytes: same key");
        let (k2, b2) = key_of(">c1\nACGTACGTTCGTAC\n>c2\nGGGGCCCCAA\n", "a.loci.fa");
        assert_eq!(b1.len(), b2.len());
        assert_ne!(k1, k2, "one base inside a span, same size: another key");
        let (k3, _) = key_of(">c1\nACGTACGTACGTAT\n>c2\nGGGGCCCCAA\n", "a.loci.fa");
        assert_eq!(k1, k3, "a base outside every span leaves the loci FASTA, and so the key, unchanged");
        let root = dir.join("cache");
        let e1 = rc::Entry::new(&root, "paf", k1.clone()).pinned();
        let st = e1.staging().unwrap();
        std::fs::write(dir.join("p.paf"), b"c1:3-12\t10\n").unwrap();
        e1.stage_link(&st, "out.paf", &dir.join("p.paf")).unwrap();
        e1.commit(&st).unwrap();
        assert!(rc::Entry::new(&root, "paf", k1).pinned().is_hit_verify(true));
        assert!(!rc::Entry::new(&root, "paf", k2).pinned().is_hit_verify(false), "a changed loci FASTA is a miss");
        let _ = std::fs::remove_dir_all(&dir);
    }
}

// ---- merged 2026-10-05: was `vg_family/family_container.rs`, now the inline module below (one component) ----
#[allow(clippy::all)]
pub mod family_container {
//! The CONTAINER of a family member's extra pieces (its ACCESSORY exon blocks) and their relations to other
//! families: the Rust port of `bench/family_container.py` (frozen sha1 e197ccb3, the binding definition of
//! `docs/PREREG_fusion_container_sim_2026-09-28.md` §1 + Amendment 1), run by `mcl_families --from-gtf
//! --emit-container` after the families are written. It reads the families products and never changes a family.
//!
//! **STATUS:** OPT-IN  (docs/MODULE_STATUS.md; `mcl_families --emit-container`, default off; driver `RUSTLE_FAMILY_CONTAINER=1`)
//!
//! Definition (prereg §1, Amendment 1):
//! * LOCUS m = a member row of `clusters.tsv` (a graph node key `CONTIG:START-END`) plus every annotation record that
//!   `loci.tsv` folds into it; FAMILY F = its `cluster_id`. A record is one assembled `gene_id`: `loci.gff3`'s `gene`
//!   line gives its key and `Name=` its gene_id (two gene_ids with the same span share one key and are both taken).
//! * EXON BLOCKS of m = the union of the exons of ALL transcripts of ALL gene_ids of m's records, merged where they
//!   OVERLAP (share >= 1 base; abutting exons stay separate).
//! * Block b of m is CORE iff some PAF record between a record of m and a record of another member m' of F has an
//!   aligned column (CIGAR M/=/X) whose m-side base lies in b and whose m'-side base is an exon base of m' (m''s own
//!   all-transcript blocks). Every PAF record counts, whatever its identity, length or primary flag. ACCESSORY = not
//!   core.
//! * RELATION: each accessory block gets the same test against the members of every OTHER family F'; the family
//!   relation F -> F' exists iff some member of F carries such a block; `reciprocal` says whether F' -> F exists too.
//!   Core blocks are not tested for relations (their relation columns are `.`).
//! * Unclustered loci get no rows and are never partners. Records between two records of the same locus are skipped.
//!
//! Coordinates: loci.fa / PAF names are record keys `CONTIG:START-END` (GFF 1-based closed) and a record's sequence is
//! the genome's FORWARD strand from START to END, so PAF offset o is genome base START + o (0-based START - 1 + o).
//! PAF `+`: the CIGAR walks target [ts,te) and query [qs,qe) ascending; `-`: target ascending against the reverse
//! complement of query [qs,qe), i.e. forward query offsets from qe-1 downwards. M/=/X consume both sides, I the query,
//! D/N the target; any other op, a CIGAR whose lengths disagree with the PAF columns, or a projected record without
//! `cg:Z` is an error. Outputs are GFF 1-based closed, blocks numbered in GENOMIC order.
//!
//! Outputs (`write`): `<out>.container.tsv` (one row per clustered locus x exon block), `<out>.container_relations.tsv`
//! (one row per directed family relation) and `<out>.container_summary.tsv` (counts): byte for byte the frozen
//! script's `OUT.blocks.tsv`, `OUT.relations.tsv` and `OUT.summary.tsv`, including its input conventions (the first
//! `key "` occurrence of a GTF attribute, an empty `gene_id` read as the transcript id, the fold table's
//! last-write-wins value at the first-write position, Python's `int()` on the decimal fields). Only a lone `\r` line
//! break (which Python's universal newlines would split on) is not reproduced.
    use crate::types::DetHashSet;
use anyhow::{bail, Context, Result};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

/// A record / locus key `(CONTIG, START, END)`, GFF 1-based closed.
pub type Key = (String, i64, i64);

pub const BLOCK_HEADER: [&str; 17] = [
    "family_id",
    "locus",
    "chrom",
    "strand",
    "gene_ids",
    "n_records",
    "block",
    "n_blocks",
    "start",
    "end",
    "bp",
    "class",
    "core_bp",
    "core_partners",
    "rel_families",
    "rel_bp",
    "rel_partners",
];
pub const REL_HEADER: [&str; 7] = ["family_id", "related_family", "n_members", "n_blocks", "bp", "reciprocal", "members"];
pub const SUMMARY_KEYS: [&str; 26] = [
    "families",
    "loci",
    "folded_records",
    "loci_gff3_keys",
    "gene_key_collisions",
    "records_with_key_collision",
    "records_span_mismatch",
    "blocks",
    "core_blocks",
    "accessory_blocks",
    "core_block_bp",
    "accessory_bp",
    "accessory_blocks_related",
    "loci_all_core",
    "loci_with_accessory",
    "families_with_relation",
    "family_relations_directed",
    "family_relations_reciprocal",
    "paf_records",
    "paf_malformed",
    "paf_self",
    "paf_unclustered",
    "paf_same_locus",
    "paf_no_exon_interval",
    "paf_projected",
    "paf_exon_exon",
];

pub fn key_str(k: &Key) -> String {
    format!("{}:{}-{}", k.0, k.1, k.2)
}

/// Python's `int(s)` on a decimal field: surrounding whitespace, an optional sign, ASCII digits with single `_`
/// between them. `None` where Python raises `ValueError` (and beyond i64).
pub fn py_int(s: &str) -> Option<i64> {
    let t = s.trim();
    let (neg, body) = match t.as_bytes().first() {
        Some(b'+') => (false, &t[1..]),
        Some(b'-') => (true, &t[1..]),
        _ => (false, t),
    };
    if body.is_empty() || body.starts_with('_') || body.ends_with('_') || body.contains("__") {
        return None;
    }
    let digits: String = body.chars().filter(|&c| c != '_').collect();
    if !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let v: i64 = digits.parse().ok()?;
    Some(if neg { -v } else { v })
}

fn int_field(s: &str, what: &str) -> Result<i64> {
    py_int(s).with_context(|| format!("invalid literal for int(): {s:?} ({what})"))
}

/// `CONTIG:START-END` -> key (the LAST `:` splits the contig, the FIRST `-` after it the range); `None` if malformed.
pub fn parse_key(name: &str) -> Option<Key> {
    let (c, r) = name.rsplit_once(':')?;
    let (a, b) = r.split_once('-')?;
    Some((c.to_string(), py_int(a)?, py_int(b)?))
}

/// The value of `key "..."` in a GTF attribute column: the FIRST occurrence of `key "` (the same substring rule as
/// `mcl_families::gtf_loci`), up to the next `"`.
pub(crate) fn gtf_attr<'a>(s: &'a str, key: &str) -> Option<&'a str> {
    let pat = format!("{key} \"");
    let i = s.find(&pat)? + pat.len();
    let j = s[i..].find('"')? + i;
    Some(&s[i..j])
}

/// Merge 0-based half-open intervals that OVERLAP (share >= 1 base); abutting intervals stay separate.
pub fn merge_blocks(mut iv: Vec<(i64, i64)>) -> Vec<(i64, i64)> {
    iv.sort();
    let mut out: Vec<(i64, i64)> = Vec::new();
    for (s, e) in iv {
        if let Some(last) = out.last_mut() {
            if s < last.1 {
                if e > last.1 {
                    last.1 = e;
                }
                continue;
            }
        }
        out.push((s, e));
    }
    out
}

/// Bases covered by 0-based half-open intervals (overlapping or abutting).
pub fn union_len<'a, I: IntoIterator<Item = &'a (i64, i64)>>(iv: I) -> i64 {
    let mut v: Vec<(i64, i64)> = iv.into_iter().copied().collect();
    v.sort();
    let mut n = 0;
    let mut cur: Option<(i64, i64)> = None;
    for (s, e) in v {
        cur = match cur {
            None => Some((s, e)),
            Some((cs, ce)) if s > ce => {
                n += ce - cs;
                Some((s, e))
            }
            Some((cs, ce)) => Some((cs, ce.max(e))),
        };
    }
    if let Some((cs, ce)) = cur {
        n += ce - cs;
    }
    n
}

/// `(\d+)([MIDNSHP=X])` tokens whose concatenation is the whole CIGAR, else an error.
fn parse_cigar(cigar: &str) -> Result<Vec<(i64, u8)>> {
    let b = cigar.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let j = i;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
        if i == j || i >= b.len() || !b"MIDNSHP=X".contains(&b[i]) {
            bail!("unparseable CIGAR {:?}", &cigar[..cigar.len().min(60)]);
        }
        let n: i64 = cigar[j..i].parse().with_context(|| format!("CIGAR length {:?}", &cigar[j..i]))?;
        out.push((n, b[i]));
        i += 1;
    }
    Ok(out)
}

/// The aligned runs of one PAF record, in the record's own offset frame: `(t0, q0, n)` = n aligned columns; column i
/// joins target offset t0 + i with query offset q0 + i on `+`, and with q0 + n - 1 - i on `-` (so [q0, q0 + n) is the
/// run's query interval on both strands).
pub fn aligned_runs(cigar: &str, strand: &str, qs: i64, qe: i64, ts: i64, te: i64) -> Result<Vec<(i64, i64, i64)>> {
    let ops = parse_cigar(cigar)?;
    if strand != "+" && strand != "-" {
        bail!("PAF strand {strand:?}");
    }
    let fwd = strand == "+";
    let (mut t, mut qc) = (ts, 0i64);
    let mut out = Vec::new();
    for (n, op) in ops {
        match op {
            b'M' | b'=' | b'X' => {
                out.push((t, if fwd { qs + qc } else { qe - qc - n }, n));
                t += n;
                qc += n;
            }
            b'I' => qc += n,
            b'D' | b'N' => t += n,
            _ => bail!("CIGAR op {} not allowed in a PAF record", op as char),
        }
    }
    if t != te || qc != qe - qs {
        bail!("CIGAR consumes target {} / query {qc} but the record spans {} / {}", t - ts, te - ts, qe - qs);
    }
    Ok(out)
}

/// `(lo, hi, block_index)` for every block of a locus (sorted, disjoint `starts`/`ends`) intersecting [a, b).
fn mask_ranges(starts: &[i64], ends: &[i64], a: i64, b: i64) -> Vec<(i64, i64, usize)> {
    let mut out = Vec::new();
    let mut i = ends.partition_point(|&e| e <= a);
    while i < starts.len() && starts[i] < b {
        let (lo, hi) = (starts[i].max(a), ends[i].min(b));
        if hi > lo {
            out.push((lo, hi, i));
        }
        i += 1;
    }
    out
}

fn overlaps_any(starts: &[i64], ends: &[i64], a: i64, b: i64) -> bool {
    let i = ends.partition_point(|&e| e <= a);
    i < starts.len() && starts[i] < b
}

/// One exon-exon stretch of a record: `(q_lo, q_hi, q_block, t_lo, t_hi, t_block)` in genome coordinates; the two
/// intervals have the same length and are joined column by column (reversed on `-`).
pub type ExonColumns = (i64, i64, usize, i64, i64, usize);

/// Aligned columns of one record joining an exon base of the query locus to an exon base of the target locus.
/// `q_off` / `t_off`: the genome 0-based coordinate of offset 0 of the query / target record; `q_blocks` / `t_blocks`:
/// `(starts, ends)` of the query / target LOCUS blocks, genome 0-based half-open.
#[allow(clippy::too_many_arguments)]
pub fn exon_columns(
    cigar: &str,
    strand: &str,
    qs: i64,
    qe: i64,
    ts: i64,
    te: i64,
    q_off: i64,
    t_off: i64,
    q_blocks: (&[i64], &[i64]),
    t_blocks: (&[i64], &[i64]),
) -> Result<Vec<ExonColumns>> {
    let (qst, qen) = q_blocks;
    let (tst, ten) = t_blocks;
    let fwd = strand == "+";
    let mut out = Vec::new();
    for (t0, q0, n) in aligned_runs(cigar, strand, qs, qe, ts, te)? {
        let (tg, qg) = (t_off + t0, q_off + q0);
        let tr: Vec<(i64, i64, usize)> =
            mask_ranges(tst, ten, tg, tg + n).into_iter().map(|(lo, hi, bi)| (lo - tg, hi - tg, bi)).collect();
        if tr.is_empty() {
            continue;
        }
        let qm = mask_ranges(qst, qen, qg, qg + n);
        if qm.is_empty() {
            continue;
        }
        let qr: Vec<(i64, i64, usize)> = if fwd {
            qm.into_iter().map(|(lo, hi, bi)| (lo - qg, hi - qg, bi)).collect()
        } else {
            // column i <-> query qg + n - 1 - i, so query [lo, hi) <-> i in [qg + n - hi, qg + n - lo)
            qm.into_iter().rev().map(|(lo, hi, bi)| (qg + n - hi, qg + n - lo, bi)).collect()
        };
        let (mut x, mut y) = (0, 0);
        while x < qr.len() && y < tr.len() {
            let (lo, hi) = (qr[x].0.max(tr[y].0), qr[x].1.min(tr[y].1));
            if hi > lo {
                let (ql, qh) = if fwd { (qg + lo, qg + hi) } else { (qg + n - hi, qg + n - lo) };
                out.push((ql, qh, qr[x].2, tg + lo, tg + hi, tr[y].2));
            }
            if qr[x].1 <= tr[y].1 {
                x += 1;
            } else {
                y += 1;
            }
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------------------- inputs

/// A text input, gzip-decoded when the path ends in `.gz` (as the script's `open_text`).
pub fn open_text(path: &str) -> Result<Box<dyn BufRead>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {path}"))?;
    Ok(if path.ends_with(".gz") {
        Box::new(std::io::BufReader::with_capacity(1 << 20, flate2::read::MultiGzDecoder::new(f)))
    } else {
        Box::new(std::io::BufReader::with_capacity(1 << 20, f))
    })
}

/// `clusters.tsv`: the family of each member key, the members per family in file order (families in
/// first-appearance order), and the member keys in first-appearance order.
struct Clusters {
    fam_of: HashMap<Key, usize>,
    fam_ids: Vec<String>,
    members: Vec<Vec<Key>>,
    keys: Vec<Key>,
}

fn read_clusters(r: &mut dyn BufRead, path: &str) -> Result<Clusters> {
    let mut lines = r.lines();
    let header_line = lines.next().transpose()?.unwrap_or_default();
    let header: Vec<&str> = header_line.split('\t').collect();
    let mut col: HashMap<&str, usize> = HashMap::new();
    for (i, h) in header.iter().enumerate() {
        col.insert(h, i);
    }
    for need in ["cluster_id", "chrom", "start", "end"] {
        if !col.contains_key(need) {
            bail!("{path}: no {need} column (header {header:?})");
        }
    }
    let (ci, cc, cs, ce) = (col["cluster_id"], col["chrom"], col["start"], col["end"]);
    let mut out = Clusters { fam_of: HashMap::new(), fam_ids: Vec::new(), members: Vec::new(), keys: Vec::new() };
    let mut order: HashMap<String, usize> = HashMap::new();
    for line in lines {
        let line = line?;
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < header.len() {
            continue;
        }
        let fid = r[ci];
        let key: Key = (r[cc].to_string(), int_field(r[cs], "clusters start")?, int_field(r[ce], "clusters end")?);
        if let Some(&f) = out.fam_of.get(&key) {
            if out.fam_ids[f] != fid {
                bail!("{path}: {} is in {} and {fid} (not a strict partition)", key_str(&key), out.fam_ids[f]);
            }
            continue;
        }
        let f = match order.get(fid) {
            Some(&f) => f,
            None => {
                order.insert(fid.to_string(), out.fam_ids.len());
                out.fam_ids.push(fid.to_string());
                out.members.push(Vec::new());
                out.fam_ids.len() - 1
            }
        };
        out.fam_of.insert(key.clone(), f);
        out.members[f].push(key.clone());
        out.keys.push(key);
    }
    Ok(out)
}

/// `loci.tsv` (`annotation representative`) -> the folds `annotation -> representative` (a != b), in first-insertion
/// order with the last value (a Python dict's semantics).
fn read_folds(r: &mut dyn BufRead, path: &str) -> Result<Vec<(Key, Key)>> {
    let mut lines = r.lines();
    let header_line = lines.next().transpose()?.unwrap_or_default();
    let header: Vec<&str> = header_line.split('\t').collect();
    if header.len() < 2 || header[0] != "annotation" || header[1] != "representative" {
        bail!("{path}: header {header:?} is not `annotation representative`");
    }
    let mut folds: Vec<(Key, Key)> = Vec::new();
    let mut at: HashMap<Key, usize> = HashMap::new();
    for line in lines {
        let line = line?;
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < 2 {
            continue;
        }
        let (Some(a), Some(b)) = (parse_key(r[0]), parse_key(r[1])) else {
            bail!("{path}: malformed keys {:?}", &r[..2]);
        };
        if a != b {
            match at.get(&a) {
                Some(&i) => folds[i].1 = b,
                None => {
                    at.insert(a.clone(), folds.len());
                    folds.push((a, b));
                }
            }
        }
    }
    Ok(folds)
}

/// `loci.gff3` `gene` lines -> ({key: [gene_id, ...]}, {key: strand of its first gene line}).
fn read_loci_gff3(r: &mut dyn BufRead, path: &str) -> Result<(HashMap<Key, Vec<String>>, HashMap<Key, String>)> {
    let mut genes: HashMap<Key, Vec<String>> = HashMap::new();
    let mut strand: HashMap<Key, String> = HashMap::new();
    for line in r.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < 9 || r[2] != "gene" {
            continue;
        }
        let key: Key = (r[0].to_string(), int_field(r[3], "gff3 start")?, int_field(r[4], "gff3 end")?);
        let Some(name) = r[8].split(';').find_map(|kv| kv.strip_prefix("Name=")) else {
            bail!("{path}: gene line without Name=: {}", line.trim());
        };
        genes.entry(key.clone()).or_default().push(name.to_string());
        strand.entry(key).or_insert_with(|| r[6].to_string());
    }
    Ok((genes, strand))
}

/// Assembled GTF -> {gene_id: [(chrom, start1, end1), ...]} over all its transcripts, for the wanted gene_ids.
/// Transcript -> gene comes from `transcript` lines and exons join by `transcript_id`, as `mcl_families::gtf_loci`.
fn read_gtf(r: &mut dyn BufRead, wanted: &HashSet<String>) -> Result<HashMap<String, Vec<(String, i64, i64)>>> {
    let mut gene_of: HashMap<String, String> = HashMap::new();
    let mut exons: HashMap<String, Vec<(String, i64, i64)>> = HashMap::new();
    for line in r.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < 9 {
            continue;
        }
        let Some(t) = gtf_attr(r[8], "transcript_id") else { continue };
        if r[2] == "transcript" {
            let g = gtf_attr(r[8], "gene_id").filter(|g| !g.is_empty()).unwrap_or(t);
            if wanted.contains(g) {
                gene_of.insert(t.to_string(), g.to_string());
            }
        } else if r[2] == "exon" {
            let ex = (r[0].to_string(), int_field(r[3], "gtf exon start")?, int_field(r[4], "gtf exon end")?);
            exons.entry(t.to_string()).or_default().push(ex);
        }
    }
    let mut out: HashMap<String, Vec<(String, i64, i64)>> = HashMap::new();
    for (t, g) in gene_of {
        let v = out.entry(g).or_default();
        if let Some(ex) = exons.get(&t) {
            v.extend(ex.iter().cloned());
        }
    }
    Ok(out)
}

// ------------------------------------------------------------------------------------------------------------ core

/// One clustered locus: its records (the member key first, then the folded ones), gene_ids, exon blocks (genome
/// 0-based half-open, sorted, disjoint) and family.
struct Locus {
    key: Key,
    family: usize,
    records: Vec<Key>,
    gene_ids: Vec<String>,
    starts: Vec<i64>,
    ends: Vec<i64>,
}

#[derive(Default)]
struct Counts(HashMap<&'static str, i64>);
impl Counts {
    fn add(&mut self, k: &'static str, n: i64) {
        *self.0.entry(k).or_insert(0) += n;
    }
    fn inc(&mut self, k: &'static str) {
        self.add(k, 1);
    }
    fn get(&self, k: &str) -> i64 {
        self.0.get(k).copied().unwrap_or(0)
    }
}

fn build_loci(
    cl: &Clusters,
    folds: &[(Key, Key)],
    loci_genes: &HashMap<Key, Vec<String>>,
    gene_exons: &HashMap<String, Vec<(String, i64, i64)>>,
    cnt: &mut Counts,
) -> Result<(Vec<Locus>, HashMap<Key, usize>)> {
    let idx_of: HashMap<&Key, usize> = cl.keys.iter().enumerate().map(|(i, k)| (k, i)).collect();
    let mut records: Vec<Vec<Key>> = cl.keys.iter().map(|m| vec![m.clone()]).collect();
    for (ann, rep) in folds {
        if let Some(&i) = idx_of.get(rep) {
            if cl.fam_of.contains_key(ann) {
                bail!("loci.tsv folds {} into {} but it is itself a cluster member", key_str(ann), key_str(rep));
            }
            records[i].push(ann.clone());
            cnt.inc("folded_records");
        }
    }
    let mut locus_of_record: HashMap<Key, usize> = HashMap::new();
    let mut loci = Vec::with_capacity(records.len());
    for (i, recs) in records.into_iter().enumerate() {
        let m = &cl.keys[i];
        let mut gids: Vec<String> = Vec::new();
        let mut exs: Vec<(i64, i64)> = Vec::new();
        for rec in &recs {
            if locus_of_record.contains_key(rec) {
                bail!("record {} belongs to two loci", key_str(rec));
            }
            locus_of_record.insert(rec.clone(), i);
            let Some(names) = loci_genes.get(rec) else {
                bail!("{} has no gene line in loci.gff3 (was the families stage run --from-gtf?)", key_str(rec));
            };
            let mut rec_ex: Vec<&(String, i64, i64)> = Vec::new();
            for g in names {
                match gene_exons.get(g) {
                    Some(v) if !v.is_empty() => rec_ex.extend(v.iter()),
                    _ => bail!("gene_id {g} of {} has no exons in the GTF", key_str(rec)),
                }
                gids.push(g.clone());
            }
            if names.len() > 1 {
                cnt.inc("records_with_key_collision");
            }
            if rec_ex.iter().any(|(c, _, _)| *c != rec.0) {
                bail!("{} has exons on another contig", key_str(rec));
            }
            let lo = rec_ex.iter().map(|x| x.1).min().unwrap();
            let hi = rec_ex.iter().map(|x| x.2).max().unwrap();
            if (lo, hi) != (rec.1, rec.2) {
                cnt.inc("records_span_mismatch"); // the key is not the gene's all-transcript span
            }
            exs.extend(rec_ex.iter().map(|x| (x.1 - 1, x.2)));
        }
        let blocks = merge_blocks(exs);
        loci.push(Locus {
            key: m.clone(),
            family: cl.fam_of[m],
            records: recs,
            gene_ids: gids,
            starts: blocks.iter().map(|b| b.0).collect(),
            ends: blocks.iter().map(|b| b.1).collect(),
        });
    }
    Ok((loci, locus_of_record))
}

/// hits[m][block][partner] = genome intervals of m's block joined to an exon base of partner (loci by index).
type Hits = Vec<HashMap<usize, HashMap<usize, Vec<(i64, i64)>>>>;

fn add_iv(lst: &mut Vec<(i64, i64)>, lo: i64, hi: i64) {
    // coalesce with the last interval when they touch (the union is what is ever read)
    if let Some(last) = lst.last_mut() {
        if lo <= last.1 && hi >= last.0 {
            *last = (lo.min(last.0), hi.max(last.1));
            return;
        }
    }
    lst.push((lo, hi));
}

fn project_paf(
    r: &mut dyn BufRead,
    path: &str,
    loci: &[Locus],
    locus_of_record: &HashMap<Key, usize>,
    cnt: &mut Counts,
) -> Result<Hits> {
    let mut hits: Hits = (0..loci.len()).map(|_| HashMap::new()).collect();
    // record name -> (its key's contig start, its locus), memoised on the raw name (names repeat on every record)
    let mut memo: HashMap<String, Option<(i64, usize)>> = HashMap::new();
    let mut lookup = |name: &str| -> Option<(i64, usize)> {
        if let Some(v) = memo.get(name) {
            return *v;
        }
        let v = parse_key(name).and_then(|k| locus_of_record.get(&k).map(|&l| (k.1, l)));
        memo.insert(name.to_string(), v);
        v
    };
    for line in r.lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            cnt.inc("paf_malformed");
            continue;
        }
        cnt.inc("paf_records");
        if f[0] == f[5] {
            cnt.inc("paf_self");
            continue;
        }
        let (Some((q_start, lq)), Some((t_start, lt))) = (lookup(f[0]), lookup(f[5])) else {
            cnt.inc("paf_unclustered");
            continue;
        };
        if lq == lt {
            cnt.inc("paf_same_locus");
            continue;
        }
        let (qs, qe) = (int_field(f[2], "PAF qs")?, int_field(f[3], "PAF qe")?);
        let strand = f[4];
        let (ts, te) = (int_field(f[7], "PAF ts")?, int_field(f[8], "PAF te")?);
        let (q_off, t_off) = (q_start - 1, t_start - 1);
        let (lqq, ltt) = (&loci[lq], &loci[lt]);
        if !(overlaps_any(&lqq.starts, &lqq.ends, q_off + qs, q_off + qe)
            && overlaps_any(&ltt.starts, &ltt.ends, t_off + ts, t_off + te))
        {
            cnt.inc("paf_no_exon_interval");
            continue;
        }
        let Some(cg) = f[12..].iter().find_map(|x| x.strip_prefix("cg:Z:")) else {
            bail!("{path}: record {} -> {} has no cg:Z CIGAR (the families PAF is run with -c)", f[0], f[5]);
        };
        cnt.inc("paf_projected");
        let cols = exon_columns(
            cg,
            strand,
            qs,
            qe,
            ts,
            te,
            q_off,
            t_off,
            (&lqq.starts, &lqq.ends),
            (&ltt.starts, &ltt.ends),
        )?;
        for &(ql, qh, qb, tl, th, tb) in &cols {
            add_iv(hits[lq].entry(qb).or_default().entry(lt).or_default(), ql, qh);
            add_iv(hits[lt].entry(tb).or_default().entry(lq).or_default(), tl, th);
        }
        if !cols.is_empty() {
            cnt.inc("paf_exon_exon");
        }
    }
    Ok(hits)
}

/// The container of one families run: the rows of the three tables (fields as written) and the counts.
pub struct Container {
    pub blocks: Vec<Vec<String>>,
    pub relations: Vec<Vec<String>>,
    counts: Counts,
}

fn classify(cl: &Clusters, loci: &[Locus], hits: &Hits, strand_of: &HashMap<Key, String>) -> Container {
    let mut cnt = Counts::default();
    let mut rows: Vec<Vec<String>> = Vec::new();
    // (family, related family) -> (members, n_blocks, bp); a family index IS its first-appearance order
    let mut rel: HashMap<(usize, usize), (Vec<usize>, i64, i64)> = HashMap::new();
    let locus_idx: HashMap<&Key, usize> = loci.iter().enumerate().map(|(i, l)| (&l.key, i)).collect();
    let empty: HashMap<usize, Vec<(i64, i64)>> = HashMap::new();
    for (fid, fam_members) in cl.members.iter().enumerate() {
        for mk in fam_members {
            let mi = locus_idx[mk];
            let l = &loci[mi];
            let nb = l.starts.len();
            let mut has_acc = false;
            for bi in 0..nb {
                let (s, e) = (l.starts[bi], l.ends[bi]);
                let partners = hits[mi].get(&bi).unwrap_or(&empty);
                let mut same: Vec<usize> = partners.keys().copied().filter(|&p| loci[p].family == fid).collect();
                same.sort_by(|&a, &b| loci[a].key.cmp(&loci[b].key));
                let mut row: Vec<String> = vec![
                    cl.fam_ids[fid].clone(),
                    key_str(&l.key),
                    l.key.0.clone(),
                    strand_of.get(&l.key).cloned().unwrap_or_else(|| ".".into()),
                    l.gene_ids.join(","),
                    l.records.len().to_string(),
                    (bi + 1).to_string(),
                    nb.to_string(),
                    (s + 1).to_string(),
                    e.to_string(),
                    (e - s).to_string(),
                ];
                cnt.inc("blocks");
                if !same.is_empty() {
                    let core_bp = union_len(same.iter().flat_map(|p| partners[p].iter()));
                    let core_partners: Vec<String> = same
                        .iter()
                        .map(|p| format!("{}={}", key_str(&loci[*p].key), union_len(partners[p].iter())))
                        .collect();
                    row.extend([
                        "core".to_string(),
                        core_bp.to_string(),
                        core_partners.join(","),
                        ".".into(),
                        ".".into(),
                        ".".into(),
                    ]);
                    cnt.inc("core_blocks");
                    cnt.add("core_block_bp", e - s);
                } else {
                    has_acc = true;
                    let mut other: Vec<usize> = partners.keys().copied().filter(|&p| loci[p].family != fid).collect();
                    other.sort_by(|&a, &b| (loci[a].family, &loci[a].key).cmp(&(loci[b].family, &loci[b].key)));
                    row.extend(["accessory".to_string(), "0".into(), ".".into()]);
                    cnt.inc("accessory_blocks");
                    cnt.add("accessory_bp", e - s);
                    if !other.is_empty() {
                        let mut fams: Vec<usize> = other.iter().map(|&p| loci[p].family).collect();
                        fams.dedup(); // `other` is sorted by family order
                        let rel_bp = union_len(other.iter().flat_map(|p| partners[p].iter()));
                        let rel_partners: Vec<String> = other
                            .iter()
                            .map(|&p| {
                                format!(
                                    "{}|{}={}",
                                    cl.fam_ids[loci[p].family],
                                    key_str(&loci[p].key),
                                    union_len(partners[&p].iter())
                                )
                            })
                            .collect();
                        row.extend([
                            fams.iter().map(|&f| cl.fam_ids[f].as_str()).collect::<Vec<_>>().join(","),
                            rel_bp.to_string(),
                            rel_partners.join(","),
                        ]);
                        cnt.inc("accessory_blocks_related");
                        for &f2 in &fams {
                            let r = rel.entry((fid, f2)).or_default();
                            if r.0.last() != Some(&mi) {
                                r.0.push(mi);
                            }
                            r.1 += 1;
                            r.2 += union_len(
                                other.iter().filter(|&&p| loci[p].family == f2).flat_map(|p| partners[p].iter()),
                            );
                        }
                    } else {
                        row.extend([".".to_string(), "0".into(), ".".into()]);
                    }
                }
                rows.push(row);
            }
            cnt.inc("loci");
            cnt.inc(if has_acc { "loci_with_accessory" } else { "loci_all_core" });
        }
    }
    let mut keys: Vec<(usize, usize)> = rel.keys().copied().collect();
    keys.sort();
    let mut rel_rows = Vec::with_capacity(keys.len());
    for (f1, f2) in keys {
        let r = &rel[&(f1, f2)];
        let reciprocal = rel.contains_key(&(f2, f1));
        rel_rows.push(vec![
            cl.fam_ids[f1].clone(),
            cl.fam_ids[f2].clone(),
            r.0.len().to_string(),
            r.1.to_string(),
            r.2.to_string(),
            if reciprocal { "yes" } else { "no" }.to_string(),
            r.0.iter().map(|&m| key_str(&loci[m].key)).collect::<Vec<_>>().join(","),
        ]);
        if reciprocal {
            cnt.inc("family_relations_reciprocal");
        }
    }
    cnt.add("family_relations_directed", rel_rows.len() as i64);
    cnt.add("families", cl.members.len() as i64);
    cnt.add("families_with_relation", rel.keys().map(|k| k.0).collect::<HashSet<_>>().len() as i64);
    Container { blocks: rows, relations: rel_rows, counts: cnt }
}

/// The whole container on open inputs: the assembled GTF, the families' `clusters.tsv`, `loci.gff3`, the fold table
/// `loci.tsv` (None = no folded records) and the all-vs-all PAF (with `cg:Z`). `paf_name` names the PAF in errors.
pub fn run(
    gtf: &mut dyn BufRead,
    clusters: &mut dyn BufRead,
    loci_gff3: &mut dyn BufRead,
    loci_tsv: Option<&mut dyn BufRead>,
    paf: &mut dyn BufRead,
    paf_name: &str,
) -> Result<Container> {
    let cl = read_clusters(clusters, "clusters.tsv")?;
    let folds = match loci_tsv {
        Some(r) => read_folds(r, "loci.tsv")?,
        None => Vec::new(),
    };
    let (loci_genes, strand_of) = read_loci_gff3(loci_gff3, "loci.gff3")?;
    let mut wanted: HashSet<String> = HashSet::new();
    for m in &cl.keys {
        wanted.extend(loci_genes.get(m).into_iter().flatten().cloned());
    }
    for (ann, rep) in &folds {
        if cl.fam_of.contains_key(rep) {
            wanted.extend(loci_genes.get(ann).into_iter().flatten().cloned());
        }
    }
    let gene_exons = read_gtf(gtf, &wanted)?;
    let mut cnt = Counts::default();
    let (loci, locus_of_record) = build_loci(&cl, &folds, &loci_genes, &gene_exons, &mut cnt)?;
    cnt.add("loci_gff3_keys", loci_genes.len() as i64);
    cnt.add("gene_key_collisions", loci_genes.values().filter(|v| v.len() > 1).count() as i64);
    let hits = project_paf(paf, paf_name, &loci, &locus_of_record, &mut cnt)?;
    let mut c = classify(&cl, &loci, &hits, &strand_of);
    for (k, v) in cnt.0 {
        c.counts.add(k, v);
    }
    Ok(c)
}

/// [`run`] on file paths (`.gz` read through gzip, as the script); `loci_tsv` None or a missing file = no folds.
pub fn run_paths(gtf: &str, clusters: &str, loci_gff3: &str, loci_tsv: Option<&str>, paf: &str) -> Result<Container> {
    let mut lt = match loci_tsv {
        Some(p) if std::path::Path::new(p).exists() => Some(open_text(p)?),
        _ => None,
    };
    run(
        &mut *open_text(gtf)?,
        &mut *open_text(clusters)?,
        &mut *open_text(loci_gff3)?,
        lt.as_deref_mut().map(|r| r as &mut dyn BufRead),
        &mut *open_text(paf)?,
        paf,
    )
}

fn tsv(header: &[&str], rows: &[Vec<String>]) -> String {
    let mut s = header.join("\t");
    s.push('\n');
    for r in rows {
        s.push_str(&r.join("\t"));
        s.push('\n');
    }
    s
}

impl Container {
    pub fn count(&self, k: &str) -> i64 {
        self.counts.get(k)
    }
    /// `<out>.container.tsv` (the script's `OUT.blocks.tsv`).
    pub fn blocks_tsv(&self) -> String {
        tsv(&BLOCK_HEADER, &self.blocks)
    }
    /// `<out>.container_relations.tsv` (the script's `OUT.relations.tsv`).
    pub fn relations_tsv(&self) -> String {
        tsv(&REL_HEADER, &self.relations)
    }
    /// `<out>.container_summary.tsv` (the script's `OUT.summary.tsv`).
    pub fn summary_tsv(&self) -> String {
        let mut s = String::from("key\tvalue\n");
        for k in SUMMARY_KEYS {
            s.push_str(&format!("{k}\t{}\n", self.count(k)));
        }
        s
    }
    /// The script's stderr summary line.
    pub fn summary_line(&self) -> String {
        let kv: Vec<String> = SUMMARY_KEYS.iter().map(|k| format!("{k}={}", self.count(k))).collect();
        format!("[family_container] {}", kv.join(" "))
    }
    /// Write the three tables next to the families products (`<out>.container.tsv`, `<out>.container_relations.tsv`,
    /// `<out>.container_summary.tsv`).
    pub fn write(&self, out: &str) -> Result<()> {
        for (suffix, text) in [
            ("container.tsv", self.blocks_tsv()),
            ("container_relations.tsv", self.relations_tsv()),
            ("container_summary.tsv", self.summary_tsv()),
        ] {
            let path = format!("{out}.{suffix}");
            std::fs::write(&path, text).with_context(|| format!("writing {path}"))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    //! The frozen script's 20 unit tests (`bench/test_family_container.py`, sha1 a2f4a13d), ported one to one on the
    //! same hand-made fixture, plus byte-identity of all three tables against the script's own outputs on that
    //! fixture (`testdata/family_container/`, written by the frozen script).
    //!
    //! The fixture (GFF 1-based closed; PAF offsets are 0-based within each record, record key START = offset 0):
    //!
    //!   family MCL0: A = c1:1001-2000  gA  tA1 1001-1100,1301-1400,1901-2000 ; tA2 1051-1150,1901-1950
    //!                  + folded record A2 = c1:1951-2100 (gA2, one exon 1951-2100)  => blocks 1001-1150, 1301-1400, 1901-2100
    //!                B = c1:5001-6000  gB  tB1 5001-5150, 5301-5400, 5701-6000
    //!   family MCL1: D = c2:1001-2000  gD  tD1 1001-1200, 1601-2000
    //!                E = c2:5001-6000  gE  tE1 5001-5200, 5501-5600, 5901-6000
    //!   unclustered: U = c3:1001-2000  gU  tU1 1001-1100, 1901-2000
    //!
    //!   R1  A->B  +  q[0,150)    t[0,150)   50M2I48M2D50M     A.b1/B.b1 exon-exon, 148 aligned columns each side
    //!   R2  A->B  +  q[399,900)  t[399,900) 501M              touches A.b2/B.b2 by exactly 1 base; ends 0 bases before A.b3
    //!   R3  A->D  -  q[700,1000) t[100,400) 100M5D95M5I100M   reverse strand: only the first 100M is exon-exon
    //!   R4  B->A  +  q[700,760)  t[450,510) 60M               B.b3 onto A's INTRON: not evidence
    //!   R5  B->U  +  q[800,900)  t[0,100)   100M              U unclustered: skipped
    //!   R6  A2->D +  q[100,150)  t[650,700) 50M               the folded record counts for A (A.b3 genome 2051-2100)
    //!   R7  A->A2 +  q[950,1000) t[0,50)    50M               same locus: skipped
    //!   R8  D->E  +  q[600,650)  t[0,50)    50M               D.b2 / E.b1 core
    //!   R9  A->U  +  q[900,1000) t[900,1000) 100M             U unclustered: skipped
    //!   R10 E->B  +  q[500,600)  t[0,100)   100M              E.b2 accessory related to MCL0 (B.b1 is core already)
    use super::*;
    use std::collections::BTreeMap;
    use std::io::Write;

    type Genes = Vec<(&'static str, &'static str, &'static str, Vec<(&'static str, Vec<(i64, i64)>)>)>;
    fn genes() -> Genes {
        vec![
            ("gA", "c1", "+", vec![("tA1", vec![(1001, 1100), (1301, 1400), (1901, 2000)]), ("tA2", vec![(1051, 1150), (1901, 1950)])]),
            ("gA2", "c1", "+", vec![("tA2x", vec![(1951, 2100)])]),
            ("gB", "c1", "+", vec![("tB1", vec![(5001, 5150), (5301, 5400), (5701, 6000)])]),
            ("gD", "c2", "-", vec![("tD1", vec![(1001, 1200), (1601, 2000)])]),
            ("gE", "c2", "-", vec![("tE1", vec![(5001, 5200), (5501, 5600), (5901, 6000)])]),
            ("gU", "c3", "+", vec![("tU1", vec![(1001, 1100), (1901, 2000)])]),
        ]
    }
    fn key(g: &str) -> &'static str {
        match g {
            "gA" => "c1:1001-2000",
            "gA2" => "c1:1951-2100",
            "gB" => "c1:5001-6000",
            "gD" => "c2:1001-2000",
            "gE" => "c2:5001-6000",
            "gU" => "c3:1001-2000",
            _ => panic!("{g}"),
        }
    }
    const CLUSTERS: [(&str, &str); 4] = [("MCL0", "gA"), ("MCL0", "gB"), ("MCL1", "gD"), ("MCL1", "gE")];
    type PafRow = (&'static str, i64, i64, &'static str, &'static str, i64, i64, &'static str);
    const PAF: [PafRow; 10] = [
        ("gA", 0, 150, "+", "gB", 0, 150, "50M2I48M2D50M"),
        ("gA", 399, 900, "+", "gB", 399, 900, "501M"),
        ("gA", 700, 1000, "-", "gD", 100, 400, "100M5D95M5I100M"),
        ("gB", 700, 760, "+", "gA", 450, 510, "60M"),
        ("gB", 800, 900, "+", "gU", 0, 100, "100M"),
        ("gA2", 100, 150, "+", "gD", 650, 700, "50M"),
        ("gA", 950, 1000, "+", "gA2", 0, 50, "50M"),
        ("gD", 600, 650, "+", "gE", 0, 50, "50M"),
        ("gA", 900, 1000, "+", "gU", 900, 1000, "100M"),
        ("gE", 500, 600, "+", "gB", 0, 100, "100M"),
    ];
    fn key_len(g: &str) -> i64 {
        let k = parse_key(key(g)).unwrap();
        k.2 - k.1 + 1
    }

    fn tmpdir(tag: &str) -> std::path::PathBuf {
        static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let n = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let d = std::env::temp_dir().join(format!("rustle_family_container_{tag}_{}_{n}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    /// `write_fixture` of the Python tests, byte for byte. Returns (gtf, fam prefix).
    fn write_fixture(d: &std::path::Path, clusters: &[(&str, &str)], folds: &[(&str, &str)], with_cg: bool) -> (String, String) {
        let gtf = d.join("x.gtf").display().to_string();
        let mut fh = std::fs::File::create(&gtf).unwrap();
        for (g, c, st, txs) in genes() {
            for (t, exs) in &txs {
                writeln!(fh, "{c}\trustle\ttranscript\t{}\t{}\t.\t{st}\t.\tgene_id \"{g}\"; transcript_id \"{t}\"; reads \"5\";", exs[0].0, exs[exs.len() - 1].1).unwrap();
                for (i, (s, e)) in exs.iter().enumerate() {
                    writeln!(fh, "{c}\trustle\texon\t{s}\t{e}\t.\t{st}\t.\tgene_id \"{g}\"; transcript_id \"{t}\"; exon_number \"{}\";", i + 1).unwrap();
                }
            }
        }
        let fam = d.join("x.fam").display().to_string();
        let mut fh = std::fs::File::create(format!("{fam}.loci.gff3")).unwrap();
        writeln!(fh, "##gff-version 3").unwrap();
        for (g, c, st, txs) in genes() {
            let k = parse_key(key(g)).unwrap();
            writeln!(fh, "{c}\t.\tgene\t{}\t{}\t.\t{st}\t.\tID=gene-{g};Name={g}", k.1, k.2).unwrap();
            for (s, e) in &txs[0].1 {
                writeln!(fh, "{c}\t.\texon\t{s}\t{e}\t.\t{st}\t.\tParent=gene-{g};gene={g}").unwrap();
            }
        }
        let mut fh = std::fs::File::create(format!("{fam}.clusters.tsv")).unwrap();
        writeln!(fh, "cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend").unwrap();
        for (fid, g) in clusters {
            let k = parse_key(key(g)).unwrap();
            let n = clusters.iter().filter(|(f, _)| f == fid).count();
            writeln!(fh, "{fid}\t{n}\t1.0000\t1.0000\tNA\t{}\t{}\t{}", k.0, k.1, k.2).unwrap();
        }
        let mut fh = std::fs::File::create(format!("{fam}.loci.tsv")).unwrap();
        writeln!(fh, "annotation\trepresentative").unwrap();
        for (a, r) in folds {
            writeln!(fh, "{}\t{}", key(a), key(r)).unwrap();
        }
        let mut fh = std::fs::File::create(format!("{fam}.loci.paf")).unwrap();
        for (qg, qs, qe, st, tg, ts, te, cg) in PAF {
            let tags = if with_cg { format!("\tNM:i:0\ttp:A:P\tcg:Z:{cg}") } else { "\tNM:i:0\ttp:A:P".to_string() };
            writeln!(
                fh,
                "{}\t{}\t{qs}\t{qe}\t{st}\t{}\t{}\t{ts}\t{te}\t{}\t{}\t60{tags}",
                key(qg),
                key_len(qg),
                key(tg),
                key_len(tg),
                qe - qs,
                (qe - qs).max(te - ts)
            )
            .unwrap();
        }
        (gtf, fam)
    }

    fn run_fam(gtf: &str, fam: &str) -> Result<Container> {
        run_paths(
            gtf,
            &format!("{fam}.clusters.tsv"),
            &format!("{fam}.loci.gff3"),
            Some(&format!("{fam}.loci.tsv")),
            &format!("{fam}.loci.paf"),
        )
    }

    /// The `Container` class fixture: the main fixture, run and written, read back as the Python tests read it.
    struct Out {
        header: Vec<String>,
        rows: BTreeMap<(String, i64), BTreeMap<String, String>>,
        rel: Vec<BTreeMap<String, String>>,
        summary: BTreeMap<String, i64>,
    }
    fn container_fixture() -> Out {
        let d = tmpdir("main");
        let (gtf, fam) = write_fixture(&d, &CLUSTERS, &[("gA2", "gA")], true);
        let out = d.join("out").display().to_string();
        run_fam(&gtf, &fam).unwrap().write(&out).unwrap();
        let read = |p: String| std::fs::read_to_string(p).unwrap();
        let blocks = read(format!("{out}.container.tsv"));
        let mut lines = blocks.lines();
        let header: Vec<String> = lines.next().unwrap().split('\t').map(String::from).collect();
        let mut rows = BTreeMap::new();
        for l in lines {
            let r: BTreeMap<String, String> = header.iter().cloned().zip(l.split('\t').map(String::from)).collect();
            rows.insert((r["locus"].clone(), r["block"].parse().unwrap()), r);
        }
        let rels = read(format!("{out}.container_relations.tsv"));
        let mut lines = rels.lines();
        let rh: Vec<String> = lines.next().unwrap().split('\t').map(String::from).collect();
        let rel = lines.map(|l| rh.iter().cloned().zip(l.split('\t').map(String::from)).collect()).collect();
        let summary = read(format!("{out}.container_summary.tsv"))
            .lines()
            .skip(1)
            .map(|l| {
                let (k, v) = l.split_once('\t').unwrap();
                (k.to_string(), v.parse().unwrap())
            })
            .collect();
        let _ = std::fs::remove_dir_all(&d);
        Out { header, rows, rel, summary }
    }
    impl Out {
        fn row(&self, g: &str, b: i64) -> &BTreeMap<String, String> {
            &self.rows[&(key(g).to_string(), b)]
        }
    }
    fn f<'a>(r: &'a BTreeMap<String, String>, k: &str) -> &'a str {
        r[k].as_str()
    }

    // ---------------------------------------------------------------------------------------- Primitives (6)

    #[test]
    fn merge_blocks_overlap_merges_abutting_stays_separate() {
        assert_eq!(merge_blocks(vec![(10, 20), (0, 10), (15, 30), (40, 50), (45, 46)]), vec![(0, 10), (10, 30), (40, 50)]);
    }

    #[test]
    fn union_len_counts_overlapping_and_abutting_once() {
        assert_eq!(union_len(&[(0, 10), (10, 20), (15, 25), (30, 31)]), 26);
        assert_eq!(union_len(&[]), 0);
    }

    #[test]
    fn aligned_runs_forward_with_eq_x_and_n() {
        let runs = aligned_runs("10=2X3N4I5M1D6M", "+", 100, 127, 50, 77).unwrap();
        assert_eq!(runs, vec![(50, 100, 10), (60, 110, 2), (65, 116, 5), (71, 121, 6)]);
    }

    #[test]
    fn aligned_runs_reverse() {
        // '-': query consumed from qe downwards; run query interval [qe - consumed - n, qe - consumed)
        let runs = aligned_runs("100M5D95M5I100M", "-", 700, 1000, 100, 400).unwrap();
        assert_eq!(runs, vec![(100, 900, 100), (205, 805, 95), (300, 700, 100)]);
    }

    #[test]
    fn cigar_length_mismatch_and_bad_ops_raise() {
        assert!(aligned_runs("10M", "+", 0, 11, 0, 10).is_err());
        assert!(aligned_runs("5S10M", "+", 0, 10, 0, 10).is_err());
        assert!(aligned_runs("10M", ".", 0, 10, 0, 10).is_err());
    }

    #[test]
    fn exon_columns_reverse_maps_each_column_exactly() {
        // query record offset 0 = genome 0; target offset 0 = genome 1000. Query exon [3,5), target exon [1000,1002):
        // '-' 10M on q[0,10) / t[0,10): column i joins t 1000+i with q 9-i, so t [1000,1002) <-> q [8,10) -- not an
        // exon on the query -- and q [3,5) <-> t [1005,1007) -- not an exon on the target. No exon-exon column.
        let q: (&[i64], &[i64]) = (&[3], &[5]);
        assert!(exon_columns("10M", "-", 0, 10, 0, 10, 0, 1000, q, (&[1000], &[1002])).unwrap().is_empty());
        // a target exon at [1005,1007) is exactly the mirror of the query exon: two columns
        let t: (&[i64], &[i64]) = (&[1005], &[1007]);
        assert_eq!(exon_columns("10M", "-", 0, 10, 0, 10, 0, 1000, q, t).unwrap(), vec![(3, 5, 0, 1005, 1007, 0)]);
        // and on '+' the same exons do not meet (q [3,5) <-> t [1003,1005))
        assert!(exon_columns("10M", "+", 0, 10, 0, 10, 0, 1000, q, t).unwrap().is_empty());
    }

    // ----------------------------------------------------------------------------------------- Container (11)

    #[test]
    fn container_header() {
        assert_eq!(container_fixture().header, BLOCK_HEADER.to_vec());
    }

    #[test]
    fn multi_transcript_and_folded_record_union() {
        let o = container_fixture();
        let a: Vec<_> = (1..=3).map(|b| o.row("gA", b)).collect();
        let spans: Vec<(i64, i64)> = a.iter().map(|r| (f(r, "start").parse().unwrap(), f(r, "end").parse().unwrap())).collect();
        assert_eq!(spans, vec![(1001, 1150), (1301, 1400), (1901, 2100)]);
        assert_eq!(f(a[0], "gene_ids"), "gA,gA2");
        assert_eq!(f(a[0], "n_records"), "2");
        assert_eq!(f(a[0], "n_blocks"), "3");
        assert!(!o.rows.contains_key(&(key("gA2").to_string(), 1)), "a folded record is part of its locus, not a row of its own");
    }

    #[test]
    fn forward_cigar_with_insertion_and_deletion() {
        let o = container_fixture();
        let (a1, b1) = (o.row("gA", 1), o.row("gB", 1));
        assert_eq!((f(a1, "class"), f(a1, "core_bp"), f(a1, "core_partners")), ("core", "148", format!("{}=148", key("gB")).as_str()));
        assert_eq!((f(b1, "class"), f(b1, "core_bp")), ("core", "148"));
        assert_eq!(f(a1, "rel_families"), ".");
    }

    #[test]
    fn one_base_touch_is_core_zero_is_accessory() {
        let o = container_fixture();
        assert_eq!((f(o.row("gA", 2), "class"), f(o.row("gA", 2), "core_bp")), ("core", "1"));
        assert_eq!((f(o.row("gB", 2), "class"), f(o.row("gB", 2), "core_bp")), ("core", "1"));
        let a3 = o.row("gA", 3); // R2 ends one base before it; its partner base there IS a B exon base
        assert_eq!(f(a3, "class"), "accessory");
        assert_eq!(f(a3, "core_partners"), ".");
    }

    #[test]
    fn reverse_strand_projection_and_relation() {
        let o = container_fixture();
        let a3 = o.row("gA", 3);
        assert_eq!(f(a3, "rel_families"), "MCL1");
        assert_eq!(f(a3, "rel_bp"), "150"); // R3's 100 (genome 1901-2000) + R6's 50 via the folded record
        assert_eq!(f(a3, "rel_partners"), format!("MCL1|{}=150", key("gD")));
        let d1 = o.row("gD", 1);
        assert_eq!((f(d1, "class"), f(d1, "rel_families"), f(d1, "rel_bp")), ("accessory", "MCL0", "100"));
        assert_eq!(f(d1, "rel_partners"), format!("MCL0|{}=100", key("gA")));
    }

    #[test]
    fn partner_intron_is_not_evidence_and_unclustered_is_ignored() {
        let o = container_fixture();
        let b3 = o.row("gB", 3);
        assert_eq!((f(b3, "class"), f(b3, "rel_families"), f(b3, "rel_bp")), ("accessory", ".", "0"));
    }

    #[test]
    fn core_blocks_are_not_tested_for_relations() {
        let o = container_fixture();
        let (d2, e1, b1) = (o.row("gD", 2), o.row("gE", 1), o.row("gB", 1));
        assert_eq!((f(d2, "class"), f(d2, "core_bp"), f(d2, "rel_families")), ("core", "50", "."));
        assert_eq!((f(e1, "class"), f(e1, "core_bp")), ("core", "50"));
        assert_eq!(f(b1, "rel_families"), "."); // R10 hits B.b1 from MCL1 but B.b1 is core
    }

    #[test]
    fn accessory_without_any_alignment() {
        let o = container_fixture();
        let e3 = o.row("gE", 3);
        assert_eq!((f(e3, "class"), f(e3, "rel_families"), f(e3, "bp")), ("accessory", ".", "100"));
        let e2 = o.row("gE", 2);
        assert_eq!((f(e2, "class"), f(e2, "rel_families"), f(e2, "rel_bp")), ("accessory", "MCL0", "100"));
    }

    #[test]
    fn unclustered_locus_has_no_rows() {
        let o = container_fixture();
        assert!(!o.rows.keys().any(|k| k.0 == key("gU")));
        assert_eq!(o.rows.len(), 3 + 3 + 2 + 3);
    }

    #[test]
    fn family_relation_table() {
        let o = container_fixture();
        let got: Vec<Vec<&str>> = o
            .rel
            .iter()
            .map(|r| ["family_id", "related_family", "n_members", "n_blocks", "bp", "reciprocal", "members"].iter().map(|k| f(r, k)).collect())
            .collect();
        let both = format!("{},{}", key("gD"), key("gE"));
        assert_eq!(
            got,
            vec![vec!["MCL0", "MCL1", "1", "1", "150", "yes", key("gA")], vec!["MCL1", "MCL0", "2", "2", "200", "yes", both.as_str()]]
        );
    }

    #[test]
    fn summary_counts() {
        let s = container_fixture().summary;
        let g = |k: &str| s[k];
        assert_eq!((g("families"), g("loci"), g("folded_records"), g("blocks")), (2, 4, 1, 11));
        assert_eq!((g("core_blocks"), g("accessory_blocks"), g("accessory_blocks_related")), (6, 5, 3));
        assert_eq!((g("paf_records"), g("paf_unclustered"), g("paf_same_locus"), g("paf_no_exon_interval")), (10, 2, 1, 1));
        assert_eq!((g("paf_projected"), g("paf_exon_exon")), (6, 6));
        assert_eq!((g("loci_all_core"), g("loci_with_accessory")), (0, 4));
        assert_eq!(g("records_span_mismatch"), 0);
    }

    // --------------------------------------------------------------------------------------------- Guards (3)

    #[test]
    fn record_without_cigar_raises() {
        let d = tmpdir("nocg");
        let (gtf, fam) = write_fixture(&d, &CLUSTERS, &[("gA2", "gA")], false);
        assert!(run_fam(&gtf, &fam).is_err());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn member_in_two_families_raises() {
        let d = tmpdir("twofam");
        let mut cl = CLUSTERS.to_vec();
        cl.push(("MCL1", "gA"));
        let (gtf, fam) = write_fixture(&d, &cl, &[("gA2", "gA")], true);
        assert!(run_fam(&gtf, &fam).is_err());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn without_the_fold_the_folded_record_is_unclustered() {
        let d = tmpdir("nofold");
        let (gtf, fam) = write_fixture(&d, &CLUSTERS, &[], true);
        let c = run_fam(&gtf, &fam).unwrap();
        let col = |name: &str| BLOCK_HEADER.iter().position(|h| *h == name).unwrap();
        let a: Vec<&Vec<String>> = c.blocks.iter().filter(|r| r[col("locus")] == key("gA")).collect();
        let spans: Vec<(&str, &str)> = a.iter().map(|r| (r[col("start")].as_str(), r[col("end")].as_str())).collect();
        assert_eq!(spans, vec![("1001", "1150"), ("1301", "1400"), ("1901", "2000")]);
        assert_eq!(a[2][col("rel_bp")], "100"); // R6 (from the no-longer-folded record) no longer counts
        assert_eq!(c.count("paf_same_locus"), 0);
        assert_eq!(c.count("paf_unclustered"), 4); // R5, R6, R7, R9
        let _ = std::fs::remove_dir_all(&d);
    }

    // ---------------------------------------------------------------- byte identity with the frozen script

    /// All three tables, byte for byte, against the frozen script's outputs on the same fixture (main and no-fold).
    #[test]
    fn tables_are_byte_identical_to_the_frozen_script() {
        for (tag, folds, blocks, rels, summary) in [
            (
                "main",
                &[("gA2", "gA")][..],
                include_str!("testdata/family_container/main.blocks.tsv"),
                include_str!("testdata/family_container/main.relations.tsv"),
                include_str!("testdata/family_container/main.summary.tsv"),
            ),
            (
                "nofold",
                &[][..],
                include_str!("testdata/family_container/nofold.blocks.tsv"),
                include_str!("testdata/family_container/nofold.relations.tsv"),
                include_str!("testdata/family_container/nofold.summary.tsv"),
            ),
        ] {
            let d = tmpdir(tag);
            let (gtf, fam) = write_fixture(&d, &CLUSTERS, folds, true);
            let c = run_fam(&gtf, &fam).unwrap();
            assert_eq!(c.blocks_tsv(), blocks, "{tag} blocks");
            assert_eq!(c.relations_tsv(), rels, "{tag} relations");
            assert_eq!(c.summary_tsv(), summary, "{tag} summary");
            let _ = std::fs::remove_dir_all(&d);
        }
    }

    #[test]
    fn python_int_and_key_parsing() {
        assert_eq!(py_int(" 12 "), Some(12));
        assert_eq!(py_int("+7"), Some(7));
        assert_eq!(py_int("-3"), Some(-3));
        assert_eq!(py_int("1_000"), Some(1000));
        for bad in ["", "1__0", "_1", "1_", "1.0", "0x1", "+-1", "1e3"] {
            assert_eq!(py_int(bad), None, "{bad:?}");
        }
        assert_eq!(parse_key("chrUn:KI270:1-20"), Some(("chrUn:KI270".to_string(), 1, 20)));
        assert_eq!(parse_key("c1:0100-200"), Some(("c1".to_string(), 100, 200)));
        assert_eq!(parse_key("c1:1-2-3"), None);
        assert_eq!(parse_key("c1"), None);
        assert_eq!(gtf_attr("gene_id \"\"; transcript_id \"t\";", "gene_id"), Some(""));
    }
}
}

// ---- merged 2026-10-05: was `vg_family/family_relations.rs`, now the inline module below (one component) ----
#[allow(clippy::all)]
pub mod family_relations {
//! The RELATION RECORDS of split transcripts and the MEMBERS of each family BY LOCUS: the container output spec v2
//! (`docs/PREREG_container_units_v2_dev_2026-09-30.md` §6.3 / Part C), run by `mcl_families --from-gtf --emit-relations`
//! after the families are written. It reads the families products and never changes a family. The port of the dev
//! prototype `relations.py` (scratch `container_units_v2/lib/`, d90a33da); on the same inputs it writes its tables byte
//! for byte (columns 1-17; column 18 is new, see below).
//!
//! **STATUS:** OPT-IN  (docs/MODULE_STATUS.md; `mcl_families --emit-relations`, default off; driver `RUSTLE_FAMILY_RELATIONS=1`)
//!
//! INPUT. The families-input GTF of `copy_assign --bridge-regroup f1units`: the UNIT transcripts `<T>.U<i>` (attributes
//! `fusion_of` = the original transcript T, `fusion_unit` `i/n`, `fusion_junction` = every cut of T,
//! `fusion_locus` = the pre-split locus key `CONTIG:START-END` (the span of ALL transcripts of T's input gene_id),
//! `fusion_gene` = that gene_id, `fusion_detector`, optional `fusion_evidence`), the clusters (`<out>.clusters.tsv`),
//! the fold table (`<out>.loci.tsv`) and the loci of the graph (one per `gene_id` of the GTF: `mcl_families`' own list,
//! so a locus key here IS a node key there). A GTF without units gives no relation rows and every locus `whole`.
//!
//! `<out>.relations.tsv`: one row per SPLIT transcript T, ranked `REL<k>` by (pre-split locus key, line of T's first
//! unit). Columns: `relation_id transcript fused_gene_id fused_locus strand n_units cut_introns detector unit_ids
//! unit_loci unit_families unit_family_sizes outcome relation separated ref_lenient ref_strict`, then
//! `detector_evidence` when the detector gave any (F1: `reads_TJ;reads_up;reads_down;share` per cut, comma-joined).
//! `cut_introns` = `S-E` per cut in genomic order (the contig is in `fused_locus`, the strand is column 5); `unit_loci` = the units'
//! new locus keys in transcription order (equal keys = one locus); `unit_families` = `MCL<k>` or `-`; `outcome` = SAME
//! (every unit in one family, none unclustered: the split is invisible) / DIFF (>= 2 families) / ONE_UNCL (some unit
//! unclustered, some clustered) / ALL_UNCL; `relation` = `cover` iff DIFF (the fused locus belongs to >= 2 families,
//! related in transcription order), else `.`; `separated` = every unit of T is in a new gene_id of its own (ALL units in
//! pairwise distinct gene_ids, the prototype's rule: not only adjacent pairs, so an A-B-A' fusion whose first and last unit
//! share a gene_id is `false`; and gene_ids, not the printed `unit_loci` keys: two gene_ids of one span count as two);
//! `ref_lenient` / `ref_strict` are BENCHMARK-ONLY (a scorer fills them against a reference partition): always `.`.
//!
//! `<out>.members_by_locus.tsv`: one row per (family, locus), a locus inheriting the families of its units, ONE member
//! per locus per family. Columns: `family locus gene_id kind via_unit_loci n_unit_loci_in_family other_families rep_tid
//! family_size_by_locus`. A clustered new locus whose gene_id holds a unit is the PRE-SPLIT locus (`kind` `fused`, `locus`
//! = its `fusion_locus`, `gene_id` = the `fusion_gene` of the last unit of that gene_id in file order, as the
//! prototype's last-write-wins map); any other is itself (`whole`). Two unit loci of one fused locus in the same family
//! are counted once (`n_unit_loci_in_family` > 1 marks it). Families, loci and keys sort as strings (Python's `sorted`).
//!
//! INVARIANTS (errors, not assertions): every unit line belongs to a complete, consistently numbered `fusion_unit`
//! group; every member of `clusters.tsv` is a locus of the GTF or a key the fold table folds; `unit_families` is
//! `clusters.tsv`'s answer for each unit locus by construction (a test recomputes it).
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::BufRead;

use anyhow::{bail, ensure, Context, Result};

// the families stage's own attribute rule (the first `key "` of the column), shared with the container
use crate::family::fam_from_gtf::family_container::{gtf_attr as attr, key_str, parse_key, py_int, Key};

pub const RELATIONS_HEADER: [&str; 17] = [
    "relation_id",
    "transcript",
    "fused_gene_id",
    "fused_locus",
    "strand",
    "n_units",
    "cut_introns",
    "detector",
    "unit_ids",
    "unit_loci",
    "unit_families",
    "unit_family_sizes",
    "outcome",
    "relation",
    "separated",
    "ref_lenient",
    "ref_strict",
];
pub const EVIDENCE_COLUMN: &str = "detector_evidence";
pub const MEMBERS_HEADER: [&str; 9] = [
    "family",
    "locus",
    "gene_id",
    "kind",
    "via_unit_loci",
    "n_unit_loci_in_family",
    "other_families",
    "rep_tid",
    "family_size_by_locus",
];

/// One locus of the graph: a `gene_id` of the families GTF, its key (`CONTIG:START-END`, GFF 1-based closed, over ALL
/// its transcripts) and its representative transcript (most `reads`, ties to the longer span, then the
/// lexicographically last id): `mcl_families::gtf_loci`'s record, in first-appearance order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocusIn {
    pub gene_id: String,
    pub key: Key,
    pub rep: String,
    pub rep_reads: i64,
}

/// One unit transcript of the GTF.
struct UnitLine {
    /// line index among the GTF's lines (ranks the relations of one fused locus)
    line: usize,
    tid: String,
    gene: String,
    strand: String,
    parent: String,
    unit: usize,
    of: usize,
    cuts: String,
    locus: Key,
    input_gene: Option<String>,
    detector: String,
    evidence: Option<String>,
}

/// Every `transcript` line carrying `fusion_unit`.
fn read_units(r: &mut dyn BufRead, path: &str) -> Result<Vec<UnitLine>> {
    let mut out = Vec::new();
    for (n, line) in r.lines().enumerate() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || f[2] != "transcript" {
            continue;
        }
        let Some(fu) = attr(f[8], "fusion_unit") else { continue };
        let what = |k: &str| format!("{path}:{}: a unit line without {k}: {}", n + 1, line.trim());
        let tid = attr(f[8], "transcript_id").with_context(|| what("transcript_id"))?;
        let (unit, of) = fu
            .split_once('/')
            .and_then(|(a, b)| Some((py_int(a)?, py_int(b)?)))
            .filter(|&(a, b)| a >= 1 && b >= a)
            .with_context(|| format!("{path}:{}: fusion_unit {fu:?} is not `i/n` with 1 <= i <= n", n + 1))?;
        let parent = attr(f[8], "fusion_of").with_context(|| what("fusion_of"))?;
        let locus = attr(f[8], "fusion_locus").with_context(|| what("fusion_locus"))?;
        let locus = parse_key(locus).with_context(|| format!("{path}:{}: fusion_locus {locus:?} is not CONTIG:START-END", n + 1))?;
        let junction = attr(f[8], "fusion_junction").with_context(|| what("fusion_junction"))?;
        let mut cuts: Vec<&str> = Vec::new();
        for tok in junction.split(',').filter(|t| !t.is_empty()) {
            // CONTIG:S-E:STRAND -> S-E
            cuts.push(tok.rsplit(':').nth(1).with_context(|| format!("{path}:{}: fusion_junction {junction:?}", n + 1))?);
        }
        out.push(UnitLine {
            line: n,
            tid: tid.to_string(),
            gene: attr(f[8], "gene_id").unwrap_or(tid).to_string(),
            strand: f[6].to_string(),
            parent: parent.to_string(),
            unit: unit as usize,
            of: of as usize,
            cuts: cuts.join(","),
            locus,
            input_gene: attr(f[8], "fusion_gene").map(str::to_string),
            detector: attr(f[8], "fusion_detector").unwrap_or(".").to_string(),
            evidence: attr(f[8], "fusion_evidence").map(str::to_string),
        });
    }
    Ok(out)
}

/// `clusters.tsv`: member key -> cluster id, and the cluster ids with their member counts, in file order.
struct Clusters {
    of: HashMap<Key, usize>,
    ids: Vec<String>,
    size: Vec<usize>,
    keys: Vec<Key>,
}

fn read_clusters(r: &mut dyn BufRead, path: &str) -> Result<Clusters> {
    let mut lines = r.lines();
    let header_line = lines.next().transpose()?.unwrap_or_default();
    let header: Vec<&str> = header_line.split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).with_context(|| format!("{path}: no {name} column (header {header:?})"));
    let (ci, cc, cs, ce) = (col("cluster_id")?, col("chrom")?, col("start")?, col("end")?);
    let mut out = Clusters { of: HashMap::new(), ids: Vec::new(), size: Vec::new(), keys: Vec::new() };
    let mut index: HashMap<String, usize> = HashMap::new();
    for line in lines {
        let line = line?;
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < header.len() {
            continue;
        }
        let key: Key = (
            r[cc].to_string(),
            py_int(r[cs]).with_context(|| format!("{path}: start {:?}", r[cs]))?,
            py_int(r[ce]).with_context(|| format!("{path}: end {:?}", r[ce]))?,
        );
        let c = *index.entry(r[ci].to_string()).or_insert_with(|| {
            out.ids.push(r[ci].to_string());
            out.size.push(0);
            out.ids.len() - 1
        });
        if let Some(&prev) = out.of.get(&key) {
            ensure!(prev == c, "{path}: {} is in {} and {} (not a strict partition)", key_str(&key), out.ids[prev], r[ci]);
            continue;
        }
        out.of.insert(key.clone(), c);
        out.size[c] += 1;
        out.keys.push(key);
    }
    Ok(out)
}

/// `loci.tsv` (`annotation representative`): the last value per annotation key.
fn read_folds(r: &mut dyn BufRead, path: &str) -> Result<HashMap<Key, Key>> {
    let mut lines = r.lines();
    let header_line = lines.next().transpose()?.unwrap_or_default();
    let header: Vec<&str> = header_line.split('\t').collect();
    ensure!(
        header.len() >= 2 && header[0] == "annotation" && header[1] == "representative",
        "{path}: header {header:?} is not `annotation representative`"
    );
    let mut folds = HashMap::new();
    for line in lines {
        let line = line?;
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < 2 {
            continue;
        }
        let (Some(a), Some(b)) = (parse_key(r[0]), parse_key(r[1])) else { bail!("{path}: malformed keys {:?}", &r[..2]) };
        folds.insert(a, b);
    }
    Ok(folds)
}

/// SAME / DIFF / ONE_UNCL / ALL_UNCL of the units' families (`None` = unclustered).
fn outcome(fams: &[Option<usize>]) -> &'static str {
    let clustered: BTreeSet<usize> = fams.iter().flatten().copied().collect();
    let unclustered = fams.iter().filter(|f| f.is_none()).count();
    if clustered.is_empty() {
        "ALL_UNCL"
    } else if unclustered > 0 {
        "ONE_UNCL"
    } else if clustered.len() == 1 {
        "SAME"
    } else {
        "DIFF"
    }
}

/// One row of `relations.tsv`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RelationRow {
    pub id: String,
    pub transcript: String,
    pub fused_gene_id: String,
    pub fused_locus: String,
    pub strand: String,
    pub n_units: usize,
    pub cut_introns: String,
    pub detector: String,
    pub unit_ids: Vec<String>,
    pub unit_loci: Vec<String>,
    /// `MCL<k>` or `-`
    pub unit_families: Vec<String>,
    pub unit_family_sizes: Vec<usize>,
    pub outcome: &'static str,
    pub separated: bool,
    pub evidence: Option<String>,
}

impl RelationRow {
    pub fn relation(&self) -> &'static str {
        if self.outcome == "DIFF" {
            "cover"
        } else {
            "."
        }
    }
}

/// One row of `members_by_locus.tsv`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MemberRow {
    pub family: String,
    pub locus: String,
    pub gene_id: String,
    pub fused: bool,
    /// the distinct unit-locus keys that put the locus in this family, sorted as strings
    pub via: Vec<String>,
    pub other_families: Vec<String>,
    pub rep_tid: String,
    pub family_size_by_locus: usize,
}

#[derive(Debug)]
pub struct Relations {
    pub relations: Vec<RelationRow>,
    pub members: Vec<MemberRow>,
}

/// The tables. `gtf`: the families-input GTF; `clusters` / `loci_tsv`: `mcl_families`' products of the same run
/// (`loci_tsv` only if that run wrote a fold table); `loci`: the graph's loci (see [`LocusIn`]).
pub fn run(
    gtf: &mut dyn BufRead,
    clusters: &mut dyn BufRead,
    loci_tsv: Option<&mut dyn BufRead>,
    loci: &[LocusIn],
) -> Result<Relations> {
    let units = read_units(gtf, "the families GTF")?;
    let cl = read_clusters(clusters, "clusters.tsv")?;
    let fold = match loci_tsv {
        Some(r) => read_folds(r, "loci.tsv")?,
        None => HashMap::new(),
    };
    let cluster_of = |key: &Key| -> Option<usize> {
        let u = if cl.of.contains_key(key) { Some(key) } else { fold.get(key) };
        u.and_then(|u| cl.of.get(u)).copied()
    };
    let fam_name = |f: Option<usize>| f.map_or("-".to_string(), |c| cl.ids[c].clone());
    let locus_of: HashMap<&str, &LocusIn> = loci.iter().map(|l| (l.gene_id.as_str(), l)).collect();

    // one group per split transcript, in order of first appearance
    let mut groups: Vec<(String, Vec<&UnitLine>)> = Vec::new();
    let mut at: HashMap<&str, usize> = HashMap::new();
    for u in &units {
        let k = *at.entry(u.parent.as_str()).or_insert_with(|| {
            groups.push((u.parent.clone(), Vec::new()));
            groups.len() - 1
        });
        groups[k].1.push(u);
    }
    for (parent, us) in groups.iter_mut() {
        us.sort_by_key(|u| u.unit);
        let n = us[0].of;
        ensure!(
            us.len() == n && us.iter().enumerate().all(|(i, u)| u.unit == i + 1 && u.of == n),
            "the units of {parent} are not exactly 1..{n} of {n} (fusion_unit {:?})",
            us.iter().map(|u| format!("{}/{}", u.unit, u.of)).collect::<Vec<_>>()
        );
        ensure!(
            us.iter().all(|u| u.locus == us[0].locus && u.cuts == us[0].cuts),
            "the units of {parent} disagree on fusion_locus or fusion_junction"
        );
    }
    // ranked by the pre-split locus, then by the line of the transcript's first unit (units replace T in its position)
    let first_line = |us: &[&UnitLine]| us.iter().map(|u| u.line).min().expect("a group has a unit");
    groups.sort_by(|a, b| (&a.1[0].locus, first_line(&a.1)).cmp(&(&b.1[0].locus, first_line(&b.1))));

    let mut relations: Vec<RelationRow> = Vec::new();
    for (k, (parent, us)) in groups.iter().enumerate() {
        let mut keys: Vec<Key> = Vec::new();
        for u in us {
            let l = locus_of.get(u.gene.as_str()).with_context(|| {
                format!("unit {} is in gene_id {}, which is not a locus of the families GTF (are the loci of another GTF?)", u.tid, u.gene)
            })?;
            keys.push(l.key.clone());
        }
        let fams: Vec<Option<usize>> = keys.iter().map(&cluster_of).collect();
        let genes: BTreeSet<&str> = us.iter().map(|u| u.gene.as_str()).collect();
        relations.push(RelationRow {
            id: format!("REL{}", k + 1),
            transcript: parent.clone(),
            fused_gene_id: us[0].input_gene.clone().unwrap_or_else(|| ".".to_string()),
            fused_locus: key_str(&us[0].locus),
            strand: us[0].strand.clone(),
            n_units: us.len(),
            cut_introns: us[0].cuts.clone(),
            detector: us[0].detector.clone(),
            unit_ids: us.iter().map(|u| u.tid.clone()).collect(),
            unit_loci: keys.iter().map(key_str).collect(),
            unit_families: fams.iter().map(|&f| fam_name(f)).collect(),
            unit_family_sizes: fams.iter().map(|f| f.map_or(0, |c| cl.size[c])).collect(),
            outcome: outcome(&fams),
            separated: genes.len() == us.len(),
            evidence: us[0].evidence.clone(),
        });
    }

    // members by locus: a new locus that holds a unit is its pre-split locus, any other is itself
    let mut input_gene: HashMap<&str, &UnitLine> = HashMap::new(); // new gene_id -> its LAST unit line
    let mut pre_split: HashMap<&str, &Key> = HashMap::new(); // input gene_id -> its pre-split locus
    for u in &units {
        input_gene.insert(u.gene.as_str(), u);
        if let Some(g) = u.input_gene.as_deref() {
            pre_split.insert(g, &u.locus);
        }
    }
    struct Acc {
        family: usize,
        locus: String,
        gene_id: String,
        fused: bool,
        via: Vec<String>,
        rep: (i64, String),
    }
    let mut mem: Vec<Acc> = Vec::new();
    let mut slot: HashMap<(usize, String), usize> = HashMap::new();
    for l in loci {
        let Some(c) = cluster_of(&l.key) else { continue };
        let (gene_id, locus, fused) = match input_gene.get(l.gene_id.as_str()) {
            Some(u) => {
                let ig = u.input_gene.as_deref().with_context(|| format!("unit {} has no fusion_gene", u.tid))?;
                let pre = pre_split.get(ig).with_context(|| format!("no fusion_locus for the gene_id {ig}"))?;
                (ig.to_string(), key_str(pre), true)
            }
            None => (l.gene_id.clone(), key_str(&l.key), false),
        };
        let k = *slot.entry((c, locus.clone())).or_insert_with(|| {
            mem.push(Acc { family: c, locus, gene_id, fused, via: Vec::new(), rep: (i64::MIN, String::new()) });
            mem.len() - 1
        });
        let a = &mut mem[k];
        a.via.push(key_str(&l.key));
        // max (reads, tid); a locus's representative ties never repeat a transcript id
        if a.rep.1.is_empty() || (l.rep_reads, l.rep.as_str()) > (a.rep.0, a.rep.1.as_str()) {
            a.rep = (l.rep_reads, l.rep.clone());
        }
    }
    let mut fams_of: HashMap<&str, BTreeSet<usize>> = HashMap::new();
    for a in &mem {
        fams_of.entry(a.locus.as_str()).or_default().insert(a.family);
    }
    let mut by_family: HashMap<usize, usize> = HashMap::new();
    for a in &mem {
        *by_family.entry(a.family).or_insert(0) += 1;
    }
    let members: Vec<MemberRow> = mem
        .iter()
        .map(|a| {
            let mut via: Vec<String> = a.via.clone();
            via.sort();
            via.dedup();
            let mut other: Vec<String> =
                fams_of[a.locus.as_str()].iter().filter(|&&f| f != a.family).map(|&f| cl.ids[f].clone()).collect();
            other.sort();
            MemberRow {
                family: cl.ids[a.family].clone(),
                locus: a.locus.clone(),
                gene_id: a.gene_id.clone(),
                fused: a.fused,
                via,
                other_families: other,
                rep_tid: a.rep.1.clone(),
                family_size_by_locus: by_family[&a.family],
            }
        })
        .collect();

    // invariants
    let known: BTreeSet<&Key> = loci.iter().map(|l| &l.key).chain(fold.keys()).collect();
    if let Some(k) = cl.keys.iter().find(|k| !known.contains(k)) {
        bail!("{} is a member of clusters.tsv but no locus of the families GTF (clusters of another run?)", key_str(k));
    }
    Ok(Relations { relations, members })
}

impl Relations {
    /// `relations.tsv` text: 17 columns, then `detector_evidence` iff any row has evidence.
    pub fn relations_tsv(&self) -> String {
        let with_evidence = self.relations.iter().any(|r| r.evidence.is_some());
        let mut out = RELATIONS_HEADER.join("\t");
        if with_evidence {
            out.push('\t');
            out.push_str(EVIDENCE_COLUMN);
        }
        out.push('\n');
        for r in &self.relations {
            let sizes: Vec<String> = r.unit_family_sizes.iter().map(|n| n.to_string()).collect();
            let cols: [String; 17] = [
                r.id.clone(),
                r.transcript.clone(),
                r.fused_gene_id.clone(),
                r.fused_locus.clone(),
                r.strand.clone(),
                r.n_units.to_string(),
                r.cut_introns.clone(),
                r.detector.clone(),
                r.unit_ids.join(","),
                r.unit_loci.join(","),
                r.unit_families.join(","),
                sizes.join(","),
                r.outcome.to_string(),
                r.relation().to_string(),
                r.separated.to_string(),
                ".".to_string(),
                ".".to_string(),
            ];
            out.push_str(&cols.join("\t"));
            if with_evidence {
                out.push('\t');
                out.push_str(r.evidence.as_deref().unwrap_or("."));
            }
            out.push('\n');
        }
        out
    }

    /// `members_by_locus.tsv` text.
    pub fn members_tsv(&self) -> String {
        let mut out = MEMBERS_HEADER.join("\t");
        out.push('\n');
        for m in &self.members {
            let list = |v: &[String]| if v.is_empty() { ".".to_string() } else { v.join(",") };
            let cols: [String; 9] = [
                m.family.clone(),
                m.locus.clone(),
                m.gene_id.clone(),
                if m.fused { "fused" } else { "whole" }.to_string(),
                if m.fused { list(&m.via) } else { ".".to_string() },
                m.via.len().to_string(),
                list(&m.other_families),
                m.rep_tid.clone(),
                m.family_size_by_locus.to_string(),
            ];
            out.push_str(&cols.join("\t"));
            out.push('\n');
        }
        out
    }

    /// Write `<out>.relations.tsv` and `<out>.members_by_locus.tsv`.
    pub fn write(&self, out: &str) -> Result<()> {
        std::fs::write(format!("{out}.relations.tsv"), self.relations_tsv()).with_context(|| format!("writing {out}.relations.tsv"))?;
        std::fs::write(format!("{out}.members_by_locus.tsv"), self.members_tsv())
            .with_context(|| format!("writing {out}.members_by_locus.tsv"))?;
        Ok(())
    }

    /// Counts for the log and `params.tsv` (the rows' names): relation rows, `cover` rows, separated rows, member rows,
    /// fused members, fused loci in >= 2 families, members counted once for two unit loci.
    pub fn counts(&self) -> [(&'static str, usize); 7] {
        let fused_multi: BTreeSet<&str> =
            self.members.iter().filter(|m| m.fused && !m.other_families.is_empty()).map(|m| m.locus.as_str()).collect();
        [
            ("relations_records", self.relations.len()),
            ("relations_cover", self.relations.iter().filter(|r| r.outcome == "DIFF").count()),
            ("relations_separated", self.relations.iter().filter(|r| r.separated).count()),
            ("relations_members_by_locus", self.members.len()),
            ("relations_members_fused", self.members.iter().filter(|m| m.fused).count()),
            ("relations_fused_loci_in_ge2_families", fused_multi.len()),
            ("relations_members_double_counted", self.members.iter().filter(|m| m.via.len() > 1).count()),
        ]
    }

    /// The outcome classes of the relation rows, for the log (`SAME`, `DIFF`, ... with their counts).
    pub fn outcome_counts(&self) -> BTreeMap<&'static str, usize> {
        let mut m = BTreeMap::new();
        for r in &self.relations {
            *m.entry(r.outcome).or_insert(0) += 1;
        }
        m
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    /// A `transcript` line (and one exon line) of a plain or unit transcript.
    fn line(chrom: &str, gene: &str, tid: &str, reads: i64, exon: (i64, i64), strand: &str, extra: &str) -> String {
        format!(
            "{chrom}\trustle\ttranscript\t{}\t{}\t.\t{strand}\t.\tgene_id \"{gene}\"; transcript_id \"{tid}\"; reads \"{reads}\";{extra}\n\
             {chrom}\trustle\texon\t{}\t{}\t.\t{strand}\t.\tgene_id \"{gene}\"; transcript_id \"{tid}\"; exon_number \"1\";",
            exon.0, exon.1, exon.0, exon.1
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn unit(gene: &str, tid: &str, exon: (i64, i64), strand: &str, i: usize, n: usize, parent: &str, locus: &str, ig: &str) -> String {
        let extra = format!(
            " fusion_of \"{parent}\"; fusion_unit \"{i}/{n}\"; fusion_junction \"c1:401-1999:{strand},c1:2501-2599:{strand}\"; \
             fusion_locus \"{locus}\"; fusion_gene \"{ig}\"; fusion_detector \"f1\"; fusion_evidence \"3;5;4;0.4286\";"
        );
        line("c1", gene, tid, 10, exon, strand, &extra)
    }

    /// The loci the graph would see: gene_id groups in first-appearance order with their span and representative.
    fn loci_of(gtf: &str) -> Vec<LocusIn> {
        let mut order: Vec<String> = Vec::new();
        type Tr = (String, i64, i64, i64, i64); // tid, reads, span, lo, hi
        let mut by: HashMap<String, Vec<Tr>> = HashMap::new();
        let mut chrom: HashMap<String, String> = HashMap::new();
        for l in gtf.lines().filter(|l| l.split('\t').nth(2) == Some("transcript")) {
            let f: Vec<&str> = l.split('\t').collect();
            let (g, t) = (attr(f[8], "gene_id").unwrap().to_string(), attr(f[8], "transcript_id").unwrap().to_string());
            if !by.contains_key(&g) {
                order.push(g.clone());
                chrom.insert(g.clone(), f[0].to_string());
            }
            let (lo, hi): (i64, i64) = (f[3].parse().unwrap(), f[4].parse().unwrap());
            by.entry(g).or_default().push((t, attr(f[8], "reads").unwrap().parse().unwrap(), hi - lo, lo, hi));
        }
        order
            .into_iter()
            .map(|g| {
                let ts = &by[&g];
                let mut sorted = ts.clone();
                sorted.sort();
                let rep = sorted.iter().max_by_key(|t| (t.1, t.2)).unwrap().clone();
                let (lo, hi) = (ts.iter().map(|t| t.3).min().unwrap(), ts.iter().map(|t| t.4).max().unwrap());
                LocusIn { key: (chrom[&g].clone(), lo, hi), gene_id: g, rep: rep.0, rep_reads: rep.1 }
            })
            .collect()
    }

    fn rel(gtf: &str, clusters: &str, folds: Option<&str>) -> Result<Relations> {
        let loci = loci_of(gtf);
        let mut fr = folds.map(|f| std::io::Cursor::new(f.to_string()));
        run(
            &mut std::io::Cursor::new(gtf.to_string()),
            &mut std::io::Cursor::new(clusters.to_string()),
            fr.as_mut().map(|r| r as &mut dyn BufRead),
            &loci,
        )
    }

    fn clusters(rows: &[(&str, i64, i64)]) -> String {
        let mut s = String::from("cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n");
        for (c, a, b) in rows {
            s.push_str(&format!("{c}\t9\t1\t1\tNA\tc1\t{a}\t{b}\n"));
        }
        s
    }

    /// Two loci: the left gene L (A, and the fusion's first unit) and the right gene R (B and its second unit); the
    /// fusion F of the plain GTF was in gene `g` with A and B (locus c1:100-2600) and is replaced by F.U1, F.U2.
    fn two_loci() -> (String, String) {
        let rows = [
            line("c1", "gL", "A", 30, (100, 600), "+", ""),
            unit("gL", "F.U1", (100, 600), "+", 1, 2, "F", "c1:100-2600", "g"),
            line("c1", "gR", "B", 30, (2000, 2600), "+", ""),
            unit("gR", "F.U2", (2000, 2600), "+", 2, 2, "F", "c1:100-2600", "g"),
        ];
        (rows.join("\n") + "\n", clusters(&[("MCL0", 100, 600), ("MCL0", 5000, 5600), ("MCL1", 2000, 2600), ("MCL1", 7000, 7600)]))
    }

    #[test]
    fn a_split_transcript_is_a_relation_and_its_locus_is_a_member_of_both_families() {
        let (gtf, cl) = two_loci();
        // the clusters name loci the GTF does not have (5000-5600, 7000-7600): invalid input, refused
        assert!(rel(&gtf, &cl, None).unwrap_err().to_string().contains("no locus of the families GTF"));
        let gtf = format!(
            "{gtf}{}\n{}\n",
            line("c1", "gX", "X", 5, (5000, 5600), "+", ""),
            line("c1", "gY", "Y", 5, (7000, 7600), "+", "")
        );
        let r = rel(&gtf, &cl, None).unwrap();
        assert_eq!(r.relations.len(), 1);
        let row = &r.relations[0];
        assert_eq!(
            (row.id.as_str(), row.transcript.as_str(), row.fused_gene_id.as_str(), row.fused_locus.as_str()),
            ("REL1", "F", "g", "c1:100-2600")
        );
        assert_eq!((row.n_units, row.cut_introns.as_str(), row.detector.as_str()), (2, "401-1999,2501-2599", "f1"));
        assert_eq!(row.unit_ids, ["F.U1", "F.U2"]);
        assert_eq!(row.unit_loci, ["c1:100-600", "c1:2000-2600"]);
        assert_eq!(row.unit_families, ["MCL0", "MCL1"]);
        assert_eq!(row.unit_family_sizes, [2, 2]);
        assert_eq!((row.outcome, row.relation(), row.separated), ("DIFF", "cover", true));
        assert_eq!(row.evidence.as_deref(), Some("3;5;4;0.4286"));
        let text = r.relations_tsv();
        let header = text.lines().next().unwrap();
        assert!(header.ends_with("\tref_lenient\tref_strict\tdetector_evidence"), "{header}");
        assert_eq!(
            text.lines().nth(1).unwrap(),
            "REL1\tF\tg\tc1:100-2600\t+\t2\t401-1999,2501-2599\tf1\tF.U1,F.U2\tc1:100-600,c1:2000-2600\tMCL0,MCL1\t2,2\tDIFF\tcover\ttrue\t.\t.\t3;5;4;0.4286"
        );
        // both unit loci are `fused` members (the pre-split locus), one per family, each naming the other family
        let m: Vec<(String, &str, &str, bool, String, &str)> = r
            .members
            .iter()
            .map(|m| (m.family.clone(), m.locus.as_str(), m.gene_id.as_str(), m.fused, m.other_families.join(","), m.rep_tid.as_str()))
            .collect();
        let want = |a: &str, b: &'static str, c: &'static str, d: bool, e: &str, f: &'static str| (a.to_string(), b, c, d, e.to_string(), f);
        assert_eq!(
            m,
            [
                want("MCL0", "c1:100-2600", "g", true, "MCL1", "A"),
                want("MCL1", "c1:100-2600", "g", true, "MCL0", "B"),
                want("MCL0", "c1:5000-5600", "gX", false, "", "X"),
                want("MCL1", "c1:7000-7600", "gY", false, "", "Y"),
            ]
        );
        // the gene_ids gL and gR hold a unit each, so A (in gL) and B (in gR) are not members of their own
        assert_eq!(r.members.iter().filter(|m| m.fused).count(), 2);
        let mtext = r.members_tsv();
        assert_eq!(mtext.lines().nth(1).unwrap(), "MCL0\tc1:100-2600\tg\tfused\tc1:100-600\t1\tMCL1\tA\t2");
        assert_eq!(mtext.lines().nth(3).unwrap(), "MCL0\tc1:5000-5600\tgX\twhole\t.\t1\t.\tX\t2");
    }

    #[test]
    fn a_gtf_without_units_has_no_relations_and_every_locus_whole_and_no_evidence_column_without_evidence() {
        let gtf = format!("{}\n{}\n", line("c1", "gX", "X", 5, (5000, 5600), "+", ""), line("c1", "gY", "Y", 5, (7000, 7600), "+", ""));
        let r = rel(&gtf, &clusters(&[("MCL0", 5000, 5600), ("MCL0", 7000, 7600)]), None).unwrap();
        assert!(r.relations.is_empty());
        assert_eq!(r.relations_tsv(), format!("{}\n", RELATIONS_HEADER.join("\t")));
        assert_eq!(r.members.len(), 2);
        assert!(r.members.iter().all(|m| !m.fused && m.family_size_by_locus == 2));
        // a unit without evidence (a list detector): exactly the 17 columns
        let (g, _) = two_loci();
        let g = g.replace(" fusion_evidence \"3;5;4;0.4286\";", "");
        let r = rel(&g, &clusters(&[("MCL0", 100, 600), ("MCL1", 2000, 2600)]), None).unwrap();
        assert_eq!(r.relations_tsv().lines().next().unwrap().split('\t').count(), 17);
    }

    #[test]
    fn outcomes_and_cover_follow_the_units_families() {
        let f = |v: &[Option<usize>]| outcome(v);
        assert_eq!(f(&[Some(1), Some(1)]), "SAME");
        assert_eq!(f(&[Some(1), Some(2)]), "DIFF");
        assert_eq!(f(&[Some(1), None]), "ONE_UNCL");
        assert_eq!(f(&[None, None]), "ALL_UNCL");
        assert_eq!(f(&[Some(1), Some(2), None]), "ONE_UNCL");
        // a split whose units sit in ONE family: the two unit loci of the fused locus are one member, counted once
        let gtf = [
            unit("gL", "F.U1", (100, 600), "+", 1, 2, "F", "c1:100-2600", "g"),
            unit("gR", "F.U2", (2000, 2600), "+", 2, 2, "F", "c1:100-2600", "g"),
            line("c1", "gX", "X", 5, (5000, 5600), "+", ""),
        ]
        .join("\n");
        let r = rel(&format!("{gtf}\n"), &clusters(&[("MCL0", 100, 600), ("MCL0", 2000, 2600), ("MCL0", 5000, 5600)]), None).unwrap();
        assert_eq!((r.relations[0].outcome, r.relations[0].relation()), ("SAME", "."));
        assert_eq!(r.members.len(), 2, "F (fused, counted once) and X");
        let fused = &r.members[0];
        assert_eq!((fused.fused, fused.via.clone(), fused.family_size_by_locus), (true, vec!["c1:100-600".to_string(), "c1:2000-2600".to_string()], 2));
        assert_eq!(r.members_tsv().lines().nth(1).unwrap(), "MCL0\tc1:100-2600\tg\tfused\tc1:100-600,c1:2000-2600\t2\t.\tF.U2\t2");
        assert_eq!(r.counts()[6], ("relations_members_double_counted", 1));
    }

    /// Strings sort as strings: `MCL10` before `MCL2` in `other_families` and in `via_unit_loci`.
    #[test]
    fn families_and_keys_sort_as_python_sorts_strings() {
        let gtf = [
            unit("g1", "F.U1", (100, 600), "+", 1, 3, "F", "c1:100-9000", "g"),
            unit("g2", "F.U2", (2000, 2600), "+", 2, 3, "F", "c1:100-9000", "g"),
            unit("g3", "F.U3", (8000, 9000), "+", 3, 3, "F", "c1:100-9000", "g"),
        ]
        .join("\n");
        let r = rel(&format!("{gtf}\n"), &clusters(&[("MCL2", 100, 600), ("MCL10", 2000, 2600), ("MCL1", 8000, 9000)]), None).unwrap();
        let by_family: HashMap<&str, &MemberRow> = r.members.iter().map(|m| (m.family.as_str(), m)).collect();
        assert_eq!(by_family["MCL2"].other_families, ["MCL1", "MCL10"]);
        assert_eq!(r.relations[0].unit_families, ["MCL2", "MCL10", "MCL1"]);
        assert_eq!(r.relations[0].unit_family_sizes, [1, 1, 1]);
    }

    #[test]
    fn rows_rank_by_pre_split_locus_then_line_and_units_by_number() {
        // two split transcripts of one fused locus, given out of order in the file, and one of an earlier locus
        let gtf = [
            unit("g2", "F2.U2", (2000, 2600), "+", 2, 2, "F2", "c1:100-2600", "g"),
            unit("g1", "F2.U1", (100, 600), "+", 1, 2, "F2", "c1:100-2600", "g"),
            unit("g3", "F1.U1", (100, 600), "+", 1, 2, "F1", "c1:100-2600", "g"),
            unit("g4", "F1.U2", (2000, 2600), "+", 2, 2, "F1", "c1:100-2600", "g"),
            unit("g5", "E.U1", (50, 60), "+", 1, 2, "E", "c1:50-99", "h"),
            unit("g6", "E.U2", (80, 99), "+", 2, 2, "E", "c1:50-99", "h"),
        ]
        .join("\n");
        let r = rel(&format!("{gtf}\n"), &clusters(&[]), None).unwrap();
        let order: Vec<(&str, &str)> = r.relations.iter().map(|x| (x.id.as_str(), x.transcript.as_str())).collect();
        assert_eq!(order, [("REL1", "E"), ("REL2", "F2"), ("REL3", "F1")], "locus key first, then the first unit's line");
        assert_eq!(r.relations[1].unit_ids, ["F2.U1", "F2.U2"]);
        assert!(r.relations.iter().all(|x| x.outcome == "ALL_UNCL" && x.relation() == "."));
    }

    #[test]
    fn a_broken_unit_group_is_an_error() {
        let ok = unit("gL", "F.U1", (100, 600), "+", 1, 2, "F", "c1:100-2600", "g");
        let e = rel(&format!("{ok}\n"), &clusters(&[]), None).unwrap_err().to_string();
        assert!(e.contains("not exactly 1..2 of 2"), "{e}");
        let no_locus = ok.replace(" fusion_locus \"c1:100-2600\";", "");
        let e = rel(&format!("{no_locus}\n"), &clusters(&[]), None).unwrap_err().to_string();
        assert!(e.contains("without fusion_locus"), "{e}");
    }

    #[test]
    fn a_locus_folded_into_a_representative_takes_the_representatives_family() {
        // gX (5000-5600) lies inside gZ (4000-7000): the graph has one node, gZ's, and gX is folded into it
        let gtf = format!("{}\n{}\n", line("c1", "gX", "X", 5, (5000, 5600), "+", ""), line("c1", "gZ", "Z", 9, (4000, 7000), "+", ""));
        let folds = "annotation\trepresentative\nc1:5000-5600\tc1:4000-7000\n";
        let r = rel(&gtf, &clusters(&[("MCL0", 4000, 7000)]), Some(folds)).unwrap();
        let m: Vec<(&str, &str)> = r.members.iter().map(|m| (m.family.as_str(), m.locus.as_str())).collect();
        assert_eq!(m, [("MCL0", "c1:5000-5600"), ("MCL0", "c1:4000-7000")], "both loci are members, in GTF order");
        // without the fold table the clustered node is alone
        let r = rel(&gtf, &clusters(&[("MCL0", 4000, 7000)]), None).unwrap();
        assert_eq!(r.members.len(), 1);
        let bad = "annotation\tnope\n";
        assert!(rel(&gtf, &clusters(&[("MCL0", 4000, 7000)]), Some(bad)).is_err());
    }

    fn fixture(name: &str) -> String {
        let p = format!("{}/tests/fixtures/bridge_units/{name}", env!("CARGO_MANIFEST_DIR"));
        std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("read {p}: {e}"))
    }

    /// The chain on the synthetic fixture of `tests/fixtures/bridge_units`: `run_list` writes the families input, this
    /// module writes the two tables, and they equal the dev prototype's (`relations.py --detector list:cuts.tsv`, d90a33da)
    /// byte for byte: DIFF (F1, FF with a double-counted fused locus, F5 through a fold), ONE_UNCL (F3: its second unit
    /// joined an unclustered locus), ALL_UNCL (FM), a locus of two gene_ids with one key, strings sorted as strings.
    #[test]
    fn the_tables_equal_the_python_prototype_on_the_fixture() {
        use crate::family::bridge_regroup::{run_list, UnitsList};
        let mut lines: Vec<String> = fixture("plain.gtf").lines().map(str::to_string).collect();
        let list = UnitsList::read(&format!("{}/tests/fixtures/bridge_units/cuts.tsv", env!("CARGO_MANIFEST_DIR"))).unwrap();
        let gtf = run_list(&mut lines, &list).unwrap().units.unwrap().families_lines.join("\n") + "\n";
        let r = rel(&gtf, &fixture("clusters.tsv"), Some(&fixture("loci.tsv"))).unwrap();
        assert_eq!(r.relations_tsv(), fixture("expected.relations.tsv"));
        assert_eq!(r.members_tsv(), fixture("expected.members_by_locus.tsv"));
        assert_eq!(r.outcome_counts().into_iter().collect::<Vec<_>>(), [("ALL_UNCL", 1), ("DIFF", 3), ("ONE_UNCL", 1)]);
        let c = r.counts();
        assert_eq!(c.iter().map(|x| x.1).collect::<Vec<_>>(), [5, 3, 5, 13, 7, 3, 1]);
        // the invariants on real-shaped output: every unit is in one row, `unit_families` is the clusters' answer for its locus
        let units: usize = r.relations.iter().map(|x| x.n_units).sum();
        assert_eq!(units, gtf.lines().filter(|l| l.contains("\ttranscript\t") && l.contains("fusion_unit \"")).count());
        let family_of: HashMap<String, String> = fixture("clusters.tsv")
            .lines()
            .skip(1)
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (format!("{}:{}-{}", f[5], f[6], f[7]), f[0].to_string())
            })
            .collect();
        let fold: HashMap<String, String> =
            fixture("loci.tsv").lines().skip(1).map(|l| l.split_once('\t').map(|(a, b)| (a.to_string(), b.to_string())).unwrap()).collect();
        for row in &r.relations {
            let want: Vec<String> = row
                .unit_loci
                .iter()
                .map(|k| family_of.get(k).or_else(|| fold.get(k).and_then(|u| family_of.get(u))).cloned().unwrap_or_else(|| "-".to_string()))
                .collect();
            assert_eq!(row.unit_families, want, "{}", row.id);
        }
        let keys: HashSet<(&str, &str)> = r.members.iter().map(|m| (m.family.as_str(), m.locus.as_str())).collect();
        assert_eq!(keys.len(), r.members.len(), "one member per locus per family");
    }

    /// `separated` = every unit of T is in a gene_id of its own, as the prototype computes it: an A-B-A' fusion (U1 and U3 in
    /// gene gA, U2 in gB) is NOT separated although no ADJACENT pair shares a locus, and two gene_ids of ONE span (gX, gY:
    /// equal `unit_loci` keys) ARE separated, the rule being about gene_ids and not about the printed keys. The F row also
    /// carries two cuts' evidence comma-joined, which the table's last column repeats verbatim.
    #[test]
    fn separated_means_every_unit_in_a_gene_id_of_its_own() {
        let evidence = "2;10;20;0.1667,2;20;10;0.1667";
        let gtf = [
            unit("gA", "F.U1", (100, 600), "+", 1, 3, "F", "c1:100-5600", "g").replace("3;5;4;0.4286", evidence),
            unit("gB", "F.U2", (2000, 2600), "+", 2, 3, "F", "c1:100-5600", "g").replace("3;5;4;0.4286", evidence),
            unit("gA", "F.U3", (5000, 5600), "+", 3, 3, "F", "c1:100-5600", "g").replace("3;5;4;0.4286", evidence),
            unit("gX", "E.U1", (7000, 7600), "+", 1, 2, "E", "c1:7000-7600", "h"),
            unit("gY", "E.U2", (7000, 7600), "+", 2, 2, "E", "c1:7000-7600", "h"),
        ]
        .join("\n");
        let r = rel(&format!("{gtf}\n"), &clusters(&[]), None).unwrap();
        assert_eq!(r.relations.len(), 2);
        let (f, e) = (&r.relations[0], &r.relations[1]);
        assert_eq!((f.transcript.as_str(), e.transcript.as_str()), ("F", "E"));
        assert_eq!(f.unit_loci, ["c1:100-5600", "c1:2000-2600", "c1:100-5600"], "adjacent units never share a locus");
        assert!(!f.separated, "U1 and U3 share gene_id gA");
        assert_eq!(e.unit_loci, ["c1:7000-7600", "c1:7000-7600"], "one printed key");
        assert!(e.separated, "gX and gY are two gene_ids");
        let text = r.relations_tsv();
        let col = |row: usize, k: usize| text.lines().nth(row).unwrap().split('\t').nth(k).unwrap().to_string();
        assert_eq!((col(1, 14).as_str(), col(2, 14).as_str()), ("false", "true"), "column 15: separated");
        // the evidence (m = 2) is column 18 verbatim; the second record has none and a `.`
        assert_eq!(text.lines().next().unwrap().split('\t').nth(17), Some("detector_evidence"));
        assert_eq!(col(1, 17), evidence);
        assert_eq!(col(2, 17), "3;5;4;0.4286");
    }
}
}

pub mod annotation_families {
    //! Multi-copy gene families defined by MCL clustering of the ANNOTATION, **by sequence**, then
    //! corroborated by RNA.
    //!
    //! ⭐ **THE CONTRIBUTION IS NOT THE CLUSTERING — IT IS THE CORROBORATION.** Clustering annotated gene
    //! sequences is what OrthoFinder does, and `docs/NEGATIVE_RESULTS_REGISTER.md:474` already diagnosed the
    //! gap it fills: transitive closure over the same homology graph gives superfamilies of 145 and 114 genes
    //! by chaining subfamilies through domain hubs — *"It is OrthoFinder without MCL, and skipping MCL is
    //! why."* What this module adds is the second stage: **DNA proposes, RNA disposes.** Measured genome-wide
    //! on gorilla (ledger §6de), a repeat clique is a PERFECT clique and a real family is not —
    //! median DNA density **1.000** for the 391/1,101 clusters with zero RNA support versus **0.700** for the
    //! corroborated ones, **at identical median size 4.0**, so the separation is not the size confound.
    //!
    //! ⭐ **NO GENE SYMBOL EVER ENTERS THE DEFINITION.** `RFPL` is 9/9 LOC-named and `GOLGA6` 14/14, and
    //! genome-wide **95.2%** of members across 1,021 product-defined families are `LOC*`-named with **67.9%**
    //! of families entirely so — a `gene=NPIP*` grep recovers **1 of 44** members of the NPIP cluster. Symbols
    //! are attached only when REPORTING (see [`Cluster::members`], which carries coordinates, not names).
    //!
    //! ⚠⚠ **THE COVERAGE DENOMINATOR IS EXONIC, NOT THE GENOMIC SPAN** (ledger §6dc). The median gorilla gene
    //! is **23.15% exon**, so a `cov_longer >= 0.30` floor measured against the SPAN demands more coverage
    //! than the median gene *has exonic sequence at all*: it removed **1,525/2,428 = 62.8%** of genes that
    //! already passed identity and length. Switching the denominator to exonic bases took the pilot graph from
    //! **903 to 1,920 nodes with 0 lost**.
    //!
    //! ⚠ **`cov_longer`, NOT `cov_shorter`.** A ~300 bp Alu covers most of a short fragment and almost none of
    //! a real gene, so a shorter-side weight lets shared repeat drive the clustering — every one of the 22
    //! adjudicated NPIP false merges was Alu-mediated (§6cr). Weighting by the LONGER side asks how much of
    //! BOTH objects is shared, which is what paralogy means and what a shared repeat cannot satisfy.
    //!
    //! **STATUS:** OTHER-BINARY  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    use std::collections::{BTreeMap, BTreeSet};

    /// Admission thresholds for an edge of the homology graph. Defaults are the values every measurement in
    /// §6da–§6de used; changing one changes the graph, so they are recorded in the params certificate.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GraphParams {
        /// Minimum alignment identity (`nmatch / blocklen`).
        pub min_identity: f64,
        /// Minimum fraction of the LONGER sequence covered by the alignment (see the module note).
        pub min_cov_longer: f64,
        /// ⭐ §6x4 CONTAINMENT ESCAPE. `0.0` = OFF, and OFF is byte-identical to every catalog built before
        /// 2026-09-22. When `> 0.0`, a pair ALSO passes if the alignment covers at least this fraction of the
        /// SHORTER sequence, even when `cov_longer` fails.
        ///
        /// ⚠ This is the norm register 913 refuted — *unguarded*. §6x4 measured why: with the exon conjunct
        /// OFF the largest component runs 5.7-16.3x baseline; with it ON and this floor >= 0.70 it holds at
        /// 2.1x while admitting 216 of the 679 loci that `cov_longer` evicts despite them aligning along
        /// ~100% of their own length (§6x3/r1002). ⚠⚠**Never enable this without the exon conjunct**
        /// (`exonic_both_sides` + `min_shared_exon_frac`) — the guard is what separates 2x from 16x.
        pub min_cov_shorter: f64,
        /// Minimum alignment block length in bp.
        pub min_bp: u64,
        /// ⭐ Charge `cov_longer`'s NUMERATOR in EXONIC bases too, instead of aligned genomic span.
        ///
        /// The span numerator and the exonic denominator are in DIFFERENT UNITS: the node's sequence is the
        /// whole gene span (introns included), so `max(qe-qs, te-ts)` counts intronic bases and divides them
        /// by an exon-union length. Measured on the pilot: **6,388/24,286 = 26.3%** of admitted edges have
        /// `aln >= denominator`, so the floor is vacuous for them, and **702/4,276 = 16.4%** of within-cluster
        /// edges are driven by an alignment that is <5% exonic. Default OFF ⟹ byte-identical.
        pub exonic_overlap: bool,
        /// ⭐ Reject a pair whose ANNOTATION INTERVALS overlap on the same contig.
        ///
        /// The only self-comparison guard is `q == t` — string equality of the FASTA header — so two distinct
        /// gene records whose intervals overlap align to each other as if they were paralogs. Measured:
        /// **108/4,276** intra-cluster edges join fully NESTED intervals at identity exactly 1.000, and
        /// 54/139 clusters contain at least one. Default OFF ⟹ byte-identical.
        pub reject_overlapping: bool,
        /// ⭐ §6ey: some single alignment record must map ≥ `min_exonic_bp` exon bases of one record onto exon
        /// bases of the other (EXON-TO-EXON homology), not only touch the longer one's exons. Records are aligned as genomic SPANS, so a pseudogene lying inside another family's gene
        /// carries that host's bases and aligns to the host's paralogs on the host's exons alone (Soto: PMS2P7,
        /// inside SPDYE8's span, joined the SPDYE cluster with 26 false pairs). Homologous copies share exon
        /// bases on both sides by definition; an alignment that touches only one gene's exons is not evidence of
        /// the other gene's homology. No new constant. Default `false` (byte-identical).
        pub exonic_both_sides: bool,
        /// ⭐ Minimum ABSOLUTE exonic bases the pair's alignments must jointly cover on the longer gene.
        ///
        /// An ADDITIVE guard, not a replacement for `cov_longer`. ⚠ Replacing the coverage measure with an
        /// exonic one (`exonic_overlap`) over-restricts: a recent segmental duplication copies INTRONS too,
        /// so a real paralog's alignment is legitimately mostly intronic — measured, it shattered the NPIP
        /// cluster 43 -> 19/14/4/4. This clause instead demands that an edge rest on SOME exonic evidence,
        /// which is what the repeat-driven merges lack (MCL2's 33 unrelated genes: <5% exonic over 1.5-2.5 kb
        /// alignments ⟹ <125 exonic bp). 0 = off ⟹ byte-identical.
        ///
        /// ⭐⭐ **SET IT TO 1, NOT 300** (§6dt). The distribution is a WALL AT ZERO, not a gradient:
        /// **53,305/63,361 = 84.1%** of candidate pairs share **exactly zero** exonic bases, only 97 pairs
        /// fall below 32 bp, and moving the floor 1 -> 300 changes the surviving set by **1.6 points**. So
        /// the threshold is not doing the work — the zero/non-zero boundary is — and the clause is properly
        /// a STRUCTURAL requirement with no free number. At 1 bp it also beats 300 on the adjudicated set
        /// (MCL3 27/27 vs 22; the repeat clique 0/33 vs 1/33). ⚠ 300 was anchored to `min_bp` and that
        /// anchor is REFUTED: 300 sits at the EDGE of the 1-200 plateau and costs MCL3 five members.
        pub min_exonic_bp: u64,
        /// ⭐⭐ §6ks — **`mcl_families` ships this ON at 0.30 (user decision 2026-09-14); the STRUCT default below
        /// stays 0.0**, the same split as `exonic_both_sides` (struct default `false`, CLI default `true`, §6ey):
        /// `GraphParams::default()` is the conservative library/test value, so a test exercising an unrelated
        /// clause via `..GraphParams::default()` is not silently perturbed by this one; `mcl_families.rs`'s own
        /// `#[arg(default_value_t = 0.30)]` is what actually ships. `--min-shared-exon-frac 0.0` reproduces every
        /// catalog built before 2026-09-14 byte-for-byte.
        ///
        /// `min_exonic_bp = 1` (`exonic_both_sides`) is a STRUCTURAL zero/non-zero gate — touch one exonic base of
        /// each gene, on some single record, however small a sliver of either gene that base is. Measured against
        /// Soto et al. 2025's family calls (`docs/o1_ledger.md` §6kr, development chr1/chr15/17): pairs BOTH
        /// definitions keep together share a median 52-89% of the smaller gene's exonic length on their best
        /// record; pairs we joined that Soto kept apart shared a median 11-25%, and >=1,000 of those pairs shared
        /// <5% or none at all — co-duplicated neighbours (NBPF beside NOTCH2NL, a lncRNA beside a GOLGA copy,
        /// GTF2I-adjacent genes) whose single shared exonic base rides on flanking segmental-duplication sequence,
        /// not on paralogy between the two genes' own models. **§6ks: T=0.30 held out on FRESH Soto families never
        /// used to pick it (chr5/7/21) — bipartite F (universe) 0.831 -> 0.881, pairwise precision (universe)
        /// 0.815 -> 1.000, ZERO Soto-verified true pairs lost on any of the 11 held-out families; TBC1D3's 9-copy
        /// family (development, reported not tested) stays whole.**
        ///
        /// This is a FRACTION of `min(exonic_len(gene_a), exonic_len(gene_b))`, not an absolute count, so it scales
        /// with the gene rather than penalising short exons: `shared_exon_bases(best record) / smaller_exonic_len`.
        /// `shared_exon_bases` is the same per-record `min(qx, tx)` `exonic_both_sides` already computes (exonic
        /// bases of each gene falling inside that record's own aligned span), maxed over the pair's records — no new
        /// alignment walk. Implies `exonic_both_sides` (the fraction cannot be computed without it; setting this
        /// without turning that on has no effect other than the wasted comparison — `exonic_both_sides` is ALSO
        /// default-on, so this is inert only behind an explicit `--no-exonic-both-sides`).
        pub min_shared_exon_frac: f64,
    }

    impl Default for GraphParams {
        fn default() -> Self {
            Self {
                min_identity: 0.70,
                min_cov_longer: 0.30,
                min_cov_shorter: 0.0,
                min_bp: 300,
                exonic_overlap: false,
                reject_overlapping: false,
                exonic_both_sides: false,
                min_exonic_bp: 0,
                min_shared_exon_frac: 0.0,
            }
        }
    }

    /// A gene, addressed the way the FASTA headers of an all-vs-all address it.
    ///
    /// ⚠ Coordinates are **GFF 1-based, verbatim** — `NC_073241.2:31346-41669` is exactly GFF `31346 41669`.
    /// A `start - 1` key joins **0/4,477** against the annotation and silently falls back to the span
    /// denominator, producing a byte-identical graph that reads as "the fix is inert" (§6dd).
    /// ⭐ **Always report the join rate.**
    pub type GeneKey = (String, u64, u64);

    /// One weighted, undirected homology graph over annotated genes.
    #[derive(Debug, Default)]
    pub struct HomologyGraph {
        /// Node index -> gene, in the order the nodes were first seen (stable for a given PAF).
        pub genes: Vec<GeneKey>,
        /// `(i, j)` with `i < j` -> weight `identity * cov_longer`, capped at 1.0.
        pub edges: BTreeMap<(usize, usize), f64>,
        /// `(i, j)` -> the best admitted record pair's IDENTITY alone (the abstention forecast, §6er: a copy whose
        /// nearest paralogue is ≥0.99 identical gets no assignment).
        pub idents: BTreeMap<(usize, usize), f64>,
        /// Genes whose exonic length was unknown, so the span was used. Reported, never silent.
        pub missing_exonic: usize,
        /// Edges whose numerator was computed from EXON BLOCKS (`exonic_overlap`). ⚠ **Always report this**:
        /// a zero join is the signature of the §6dd coordinate bug, which produced a byte-identical graph
        /// that read as "the fix is inert".
        pub exonic_overlap_joined: usize,
        /// Edges that reached the numerator but had NO exon blocks, so the span numerator was used.
        pub exonic_overlap_missing: usize,
        /// Pairs dropped by `reject_overlapping`. Reported, never silent.
        pub rejected_overlapping: usize,
        /// Pairs dropped by `min_exonic_bp` — the edge rested on no exonic evidence. Reported, never silent.
        pub rejected_no_exonic: usize,
        /// Pairs dropped by `min_shared_exon_frac` — some exonic evidence existed, but on the best record it
        /// covered too small a fraction of the smaller gene. Reported, never silent.
        pub rejected_low_shared_exon: usize,
        /// ⭐ §6x4: pairs admitted ONLY by the containment escape (`min_cov_shorter`), i.e. `cov_longer` failed
        /// and `cov_shorter` passed. Always 0 when the escape is off. Reported, never silent.
        pub admitted_by_containment: usize,
        /// PAF records between two annotations of ONE locus (`LocusMap`) — a locus aligned to itself —
        /// skipped. Reported, never silent.
        pub same_locus_records: usize,
    }

    /// ⭐ A LOCUS, not an annotation record, is the node (§6ee).
    ///
    /// Two annotation records whose EXON-UNIONS share ≥1 bp on one contig are the SAME transcribed DNA — a
    /// lncRNA model drawn over a gene's exons, two Gnomon models of one transcription unit, an antisense
    /// model over the same exons. Left as separate nodes they align to each other at identity 1.000 and each
    /// collects the same paralogue edges, so a family counts one locus twice and O2 ties by construction
    /// (NPIP: 11 such pairs, 607/1,221 tied reads). Components of that relation are one locus, represented by
    /// the record with the greatest exon-union length (ties: lowest start, then lowest end).
    /// ⚠ GENOMIC overlap is NOT the criterion: a gene inside another gene's intron is a distinct locus, and
    /// merging by genomic overlap folded 7,575 intronic genes genome-wide and dissolved MCL6 (22/23 members
    /// sit in one host's introns). Strand is ignored on purpose: for DNA homology an antisense model over the
    /// same exons is the same locus.
    #[derive(Debug, Default)]
    pub struct LocusMap {
        /// Every annotation in a multi-record locus -> its representative (representatives map to themselves).
        pub rep_of: BTreeMap<GeneKey, GeneKey>,
        /// Number of loci that hold more than one annotation record.
        pub n_multi: usize,
        /// Evidence policy. `false` (REPRESENTATIVE-ONLY, the default since §6eq): a locus's homology is its
        /// representative record's alignments; records of folded-away models are skipped. `true` (ATTRIBUTION):
        /// every model's admitted edge is the locus's edge — measured to reconstruct the duplication BLOCK
        /// (LCR16a + LCR16u merge, §6eg) and kept only as an explicit option.
        pub attribute_edges: bool,
    }

    impl LocusMap {
        pub fn representative<'a>(&'a self, k: &'a GeneKey) -> &'a GeneKey {
            self.rep_of.get(k).unwrap_or(k)
        }
        pub fn is_representative(&self, k: &GeneKey) -> bool {
            self.rep_of.get(k).map_or(true, |r| r == k)
        }
        /// Annotation records folded into another record's locus.
        pub fn n_merged(&self) -> usize {
            self.rep_of.iter().filter(|(k, r)| k != r).count()
        }
    }

    /// Build the locus map from every gene's merged exon blocks (absolute coordinates, half-open, as
    /// `mcl_families::exonic_blocks` emits them). A gene with no blocks is its own locus.
    pub fn loci_from_exon_blocks(exon_blocks: &BTreeMap<GeneKey, Vec<(u64, u64)>>) -> LocusMap {
        let genes: Vec<&GeneKey> = exon_blocks.keys().collect();
        let idx: BTreeMap<&GeneKey, usize> = genes.iter().enumerate().map(|(i, g)| (*g, i)).collect();
        let mut parent: Vec<usize> = (0..genes.len()).collect();
        fn find(p: &mut Vec<usize>, mut x: usize) -> usize {
            while p[x] != x {
                p[x] = p[p[x]];
                x = p[x];
            }
            x
        }
        // Sweep every exon block per contig; a block starting before the running max end overlaps the
        // block that set it, so the two genes are joined.
        let mut by_contig: BTreeMap<&str, Vec<(u64, u64, usize)>> = BTreeMap::new();
        for (g, blocks) in exon_blocks {
            for &(s, e) in blocks {
                by_contig.entry(g.0.as_str()).or_default().push((s, e, idx[g]));
            }
        }
        for (_, mut v) in by_contig {
            v.sort_unstable();
            let mut max_end = 0u64;
            let mut owner = usize::MAX;
            for (s, e, gi) in v {
                if owner != usize::MAX && s < max_end {
                    let (a, b) = (find(&mut parent, gi), find(&mut parent, owner));
                    if a != b {
                        parent[a.max(b)] = a.min(b);
                    }
                }
                if e > max_end || owner == usize::MAX {
                    max_end = max_end.max(e);
                    owner = gi;
                }
            }
        }
        let mut comps: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for i in 0..genes.len() {
            let r = find(&mut parent, i);
            comps.entry(r).or_default().push(i);
        }
        let exlen = |i: usize| exon_blocks[genes[i]].iter().map(|(s, e)| e - s).sum::<u64>();
        let mut m = LocusMap::default();
        for (_, members) in comps {
            if members.len() < 2 {
                continue;
            }
            let rep = *members
                .iter()
                .max_by(|&&a, &&b| exlen(a).cmp(&exlen(b)).then_with(|| (genes[b].1, genes[b].2).cmp(&(genes[a].1, genes[a].2))))
                .unwrap();
            for &i in &members {
                m.rep_of.insert(genes[i].clone(), genes[rep].clone());
            }
            m.n_multi += 1;
        }
        m
    }

    impl HomologyGraph {
        pub fn n_nodes(&self) -> usize {
            self.genes.len()
        }
        pub fn n_edges(&self) -> usize {
            self.edges.len()
        }
    }

    /// Parse one gene key out of a FASTA/PAF name of the form `CONTIG:START-END`.
    pub fn parse_gene_key(name: &str) -> Option<GeneKey> {
        let (c, r) = name.rsplit_once(':')?;
        let (a, b) = r.split_once('-')?;
        Some((c.to_string(), a.parse().ok()?, b.parse().ok()?))
    }

    /// Build the weighted graph from an all-vs-all PAF.
    ///
    /// `exonic_len` maps a gene to its EXON-UNION length; a gene absent from the map falls back to the PAF's
    /// own sequence length (the genomic span) and is counted in [`HomologyGraph::missing_exonic`].
    pub fn graph_from_paf(
        paf: &str,
        exonic_len: &BTreeMap<GeneKey, u64>,
        exon_blocks: &BTreeMap<GeneKey, Vec<(u64, u64)>>,
        p: &GraphParams,
    ) -> HomologyGraph {
        graph_from_paf_loci(paf, exonic_len, exon_blocks, p, None)
    }

    /// [`graph_from_paf`] with an optional [`LocusMap`]: when given, the NODE of every record is its locus
    /// representative, while the record's own lengths and exon blocks still decide admission — so an edge
    /// admitted for ANY annotation of a locus is an edge of the locus (a folded lncRNA model may carry an
    /// alignment its host record lacks: NPIP's 4,743-read copy had only such edges). Records between two
    /// annotations of one locus are skipped (`same_locus_records`). Coordinates are never remapped.
    pub fn graph_from_paf_loci(
        paf: &str,
        exonic_len: &BTreeMap<GeneKey, u64>,
        exon_blocks: &BTreeMap<GeneKey, Vec<(u64, u64)>>,
        p: &GraphParams,
        loci: Option<&LocusMap>,
    ) -> HomologyGraph {
        let mut g = HomologyGraph::default();
        let mut idx: BTreeMap<GeneKey, usize> = BTreeMap::new();
        let mut missing: BTreeSet<GeneKey> = BTreeSet::new();

        let mut node_of = |k: GeneKey, genes: &mut Vec<GeneKey>| -> usize {
            if let Some(&i) = idx.get(&k) {
                return i;
            }
            let i = genes.len();
            genes.push(k.clone());
            idx.insert(k, i);
            i
        };

        // ⭐ Under `exonic_overlap`, coverage is a property of the PAIR, not of one PAF record: minimap2
        // splits one paralogous alignment across many records, so a per-record floor rejects a real family
        // whose exons are jointly covered. Accumulate the pair's intervals, threshold ONCE at the end.
        // (Measured: per-record thresholding shattered the NPIP cluster 43 -> 17/15/3/3/1.)
        type PairAcc = (Vec<(u64, u64)>, Vec<(u64, u64)>, u64, u64, u64);
        let mut acc: BTreeMap<(GeneKey, GeneKey), PairAcc> = BTreeMap::new();

        for line in paf.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 11 {
                continue;
            }
            let (q, t) = (f[0], f[5]);
            if q == t {
                continue;
            }
            let (Some(qk), Some(tk)) = (parse_gene_key(q), parse_gene_key(t)) else { continue };
            if let Some(l) = loci {
                if l.representative(&qk) == l.representative(&tk) {
                    g.same_locus_records += 1;
                    continue;
                }
                if !l.attribute_edges && (!l.is_representative(&qk) || !l.is_representative(&tk)) {
                    g.same_locus_records += 1; // representative-only: a folded model's record carries no evidence
                    continue;
                }
            }
            // Two DISTINCT gene records whose 1-based inclusive intervals intersect are not paralogs; the
            // `q == t` guard above only catches byte-identical headers.
            if p.reject_overlapping && qk.0 == tk.0 && qk.1 <= tk.2 && tk.1 <= qk.2 {
                g.rejected_overlapping += 1;
                continue;
            }
            let (Ok(ql), Ok(qs), Ok(qe)) = (f[1].parse::<u64>(), f[2].parse::<u64>(), f[3].parse::<u64>())
            else {
                continue;
            };
            let (Ok(tl), Ok(ts), Ok(te)) = (f[6].parse::<u64>(), f[7].parse::<u64>(), f[8].parse::<u64>())
            else {
                continue;
            };
            let (Ok(nmatch), Ok(blocklen)) = (f[9].parse::<u64>(), f[10].parse::<u64>()) else { continue };
            if blocklen < p.min_bp {
                continue;
            }
            let identity = nmatch as f64 / blocklen.max(1) as f64;
            if identity < p.min_identity {
                continue;
            }
            let aln = (qe.saturating_sub(qs)).max(te.saturating_sub(ts));
            let dq = match exonic_len.get(&qk) {
                Some(&v) => v,
                None => {
                    missing.insert(qk.clone());
                    ql
                }
            };
            let dt = match exonic_len.get(&tk) {
                Some(&v) => v,
                None => {
                    missing.insert(tk.clone());
                    tl
                }
            };
            if p.exonic_overlap || p.min_exonic_bp > 0 {
                // Defer: bank this record's intervals on each side and decide once the pair is complete.
                let (k, qf) = if qk <= tk {
                    ((qk.clone(), tk.clone()), true)
                } else {
                    ((tk.clone(), qk.clone()), false)
                };
                let e = acc.entry(k).or_insert_with(|| (Vec::new(), Vec::new(), 0, 0, 0));
                if qf {
                    e.0.push((qs, qe));
                    e.1.push((ts, te));
                } else {
                    e.0.push((ts, te));
                    e.1.push((qs, qe));
                }
                e.2 += nmatch;
                e.3 += blocklen;
                if p.exonic_both_sides || p.min_shared_exon_frac > 0.0 {
                    // exon-to-exon: THIS record's interval must touch exon bases on both sides; keep the best record
                    let qx = exon_blocks.get(&qk).map_or(0, |b| exonic_bases_in(b, qk.1, qs, qe));
                    let tx = exon_blocks.get(&tk).map_or(0, |b| exonic_bases_in(b, tk.1, ts, te));
                    e.4 = e.4.max(qx.min(tx));
                }
                continue;
            }
            let cov_longer = (aln as f64 / dq.max(dt).max(1) as f64).min(1.0);
            if cov_longer < p.min_cov_longer {
                continue;
            }
            let w = identity * cov_longer;
            let (qn, tn) = match loci {
                Some(l) => (l.representative(&qk).clone(), l.representative(&tk).clone()),
                None => (qk, tk),
            };
            let (a, b) = (node_of(qn, &mut g.genes), node_of(tn, &mut g.genes));
            let key = if a < b { (a, b) } else { (b, a) };
            let e = g.edges.entry(key).or_insert(0.0);
            if w > *e {
                *e = w;
            }
            let id_e = g.idents.entry(key).or_insert(0.0);
            if identity > *id_e {
                *id_e = identity;
            }
        }
        // Resolve the deferred pairs: merge each side's intervals, charge the LONGER gene's own exonic bases
        // against its own exonic denominator, and threshold once.
        for ((ak, bk), (aiv, biv, nmatch, blocklen, exon_exon)) in acc {
            let da = exonic_len.get(&ak).copied().unwrap_or(0);
            let db = exonic_len.get(&bk).copied().unwrap_or(0);
            let a_is_longer = da >= db;
            let (gk, iv, den) = if a_is_longer { (&ak, aiv.clone(), da) } else { (&bk, biv.clone(), db) };
            // §6x4: the SHORTER side of the same pair, kept for the containment escape below.
            let (sk, siv, sden) = if a_is_longer { (&bk, biv, db) } else { (&ak, aiv, da) };
            let Some(blocks) = exon_blocks.get(gk) else {
                g.exonic_overlap_missing += 1;
                continue;
            };
            g.exonic_overlap_joined += 1;
            let merged = merge_intervals(iv);
            let covered =
                merged.iter().map(|&(s0, e0)| exonic_bases_in(blocks, gk.1, s0, e0)).sum::<u64>();
            if covered < p.min_exonic_bp {
                g.rejected_no_exonic += 1;
                continue;
            }
            if p.exonic_both_sides && exon_exon < p.min_exonic_bp {
                // no single record maps exon bases of one gene onto exon bases of the other: the homology is
                // between spans (a nested pseudogene carrying its host's bases, co-duplicated neighbours), not
                // between the genes
                g.rejected_no_exonic += 1;
                continue;
            }
            if p.min_shared_exon_frac > 0.0 {
                // the same exon-to-exon evidence as above, but as a FRACTION of the smaller gene's exonic length:
                // one shared base passes the structural gate above yet can be a sliver of either gene's model
                // (§6ks — co-duplicated neighbours pass at exactly this point).
                let smaller = da.min(db).max(1);
                if (exon_exon as f64 / smaller as f64) < p.min_shared_exon_frac {
                    g.rejected_low_shared_exon += 1;
                    continue;
                }
            }
            // `exonic_overlap` REPLACES the numerator; otherwise keep the span numerator (unioned over the
            // pair's records, so a split alignment is not penalised for being split).
            let numer = if p.exonic_overlap {
                covered
            } else {
                merged.iter().map(|&(s0, e0)| e0 - s0).sum::<u64>()
            };
            let cov_longer = (numer as f64 / den.max(1) as f64).min(1.0);
            // ⭐ §6x4 CONTAINMENT ESCAPE. `min_cov_shorter == 0.0` skips this block entirely, so `cov_eff`
            // is `cov_longer` and both the gate and the edge weight are byte-identical to every prior catalog.
            let mut cov_eff = cov_longer;
            // ⛔ r1010: a positional-overlap guard was tried here and REFUTED — see the register. Two de novo
            // loci that overlap on the genome are frequently genuine tandem copies, so refusing the escape for
            // them cost more than it saved (chr16 de novo F .230 -> .199, BELOW baseline .214) while doing
            // nothing for the SD-region node set it was designed to protect (.175 either way).
            if p.min_cov_shorter > 0.0 && cov_longer < p.min_cov_longer {
                let merged_s = merge_intervals(siv);
                let numer_s = if p.exonic_overlap {
                    exon_blocks
                        .get(sk)
                        .map_or(0, |b| merged_s.iter().map(|&(s0, e0)| exonic_bases_in(b, sk.1, s0, e0)).sum())
                } else {
                    merged_s.iter().map(|&(s0, e0)| e0 - s0).sum::<u64>()
                };
                let cov_shorter = (numer_s as f64 / sden.max(1) as f64).min(1.0);
                if cov_shorter >= p.min_cov_shorter {
                    g.admitted_by_containment += 1;
                    // ⭐⭐ §6x8: the weight MUST be `cov_shorter`. Charging the pair's mutual coverage instead
                    // (r1019) makes the whole flag a NO-OP — every arm returns to its baseline F to three
                    // decimals — because the escaped edge is then too weak for MCL to route through.
                    // **The gain comes from the WEIGHT, not from the admission**: putting the node in the
                    // graph changes nothing; ranking its edge highly is what changes the partition.
                    cov_eff = cov_shorter;
                }
            }
            if cov_eff < p.min_cov_longer && cov_eff < p.min_cov_shorter {
                continue;
            }
            if p.min_cov_shorter <= 0.0 && cov_eff < p.min_cov_longer {
                continue;
            }
            let identity = nmatch as f64 / blocklen.max(1) as f64;
            let w = identity * cov_eff;
            let (an, bn) = match loci {
                Some(l) => (l.representative(&ak).clone(), l.representative(&bk).clone()),
                None => (ak, bk),
            };
            let (a, b) = (node_of(an, &mut g.genes), node_of(bn, &mut g.genes));
            let key = if a < b { (a, b) } else { (b, a) };
            let e = g.edges.entry(key).or_insert(0.0);
            if w > *e {
                *e = w;
            }
            let id_e = g.idents.entry(key).or_insert(0.0);
            if identity > *id_e {
                *id_e = identity;
            }
        }
        g.missing_exonic = missing.len();
        g
    }

    // ---------------------------------------------------------------------------------------------------
    // ⭐ DUPLICON-FIRST CORE REFINEMENT (§6eh, pre-registered adj/core/PREREG.md)
    //
    // A gene family is the set of loci sharing a CORE — the segment linked by segmental-duplication pairs
    // to most of the family — not the set of annotation models that happen to align. On NPIP the core is the
    // ~23 kb LCR16a duplicon: it sits whole inside every 125–308 kb chimeric model (ABCC1+NPIP, SORL1+NPIP),
    // EIF3C carries 808 bp of it, the ABCC1-region records none, and no SMG1P (LCR16u) record is linked to
    // even half of NPIP. One constant, "half", used three times; no other number.
    // ---------------------------------------------------------------------------------------------------

    /// SEDEF pairs (BED, 0-based half-open, cols 1–3 and 4–6), indexed by the contig of EITHER side.
    #[derive(Debug, Default, Clone)]
    pub struct SdPairs {
        /// contig -> (start, end, other_contig, other_start, other_end), sorted by start.
        by_contig: BTreeMap<String, Vec<(u64, u64, String, u64, u64)>>,
        max_len: BTreeMap<String, u64>,
    }

    impl SdPairs {
        pub fn from_bed_str(text: &str) -> SdPairs {
            let mut sd = SdPairs::default();
            for line in text.lines() {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() < 6 {
                    continue;
                }
                let (Ok(s1), Ok(e1), Ok(s2), Ok(e2)) =
                    (f[1].parse::<u64>(), f[2].parse::<u64>(), f[4].parse::<u64>(), f[5].parse::<u64>())
                else {
                    continue;
                };
                sd.push(f[0], s1, e1, f[3], s2, e2);
                sd.push(f[3], s2, e2, f[0], s1, e1);
            }
            for v in sd.by_contig.values_mut() {
                v.sort_unstable();
            }
            sd
        }
        /// ⭐ Pairs derived from the catalog's OWN alignments (§6fo): every PAF record between two annotated loci
        /// (`chrom:start-end` names, GFF 1-based as written; offsets 0-based on the query/target) becomes one pair of
        /// genomic intervals — the same object SEDEF supplies, restricted to what the annotated loci align to each
        /// other. No external SD caller needed; what is lost is the intergenic SD context (`blocks.tsv`).
        pub fn from_paf_str(text: &str) -> SdPairs {
            let mut sd = SdPairs::default();
            for line in text.lines() {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() < 11 {
                    continue;
                }
                let (Some((qc, qs0, _)), Some((tc, ts0, _))) = (parse_gene_key(f[0]), parse_gene_key(f[5])) else { continue };
                let (Ok(qs), Ok(qe), Ok(ts), Ok(te)) = (f[2].parse::<u64>(), f[3].parse::<u64>(), f[7].parse::<u64>(), f[8].parse::<u64>()) else { continue };
                let (a1, b1) = (qs0.saturating_sub(1) + qs, qs0.saturating_sub(1) + qe);
                let (a2, b2) = (ts0.saturating_sub(1) + ts, ts0.saturating_sub(1) + te);
                if b1 <= a1 || b2 <= a2 {
                    continue;
                }
                sd.push(&qc, a1, b1, &tc, a2, b2);
                sd.push(&tc, a2, b2, &qc, a1, b1);
            }
            for v in sd.by_contig.values_mut() {
                v.sort_unstable();
            }
            sd
        }
        fn push(&mut self, c: &str, s: u64, e: u64, oc: &str, os: u64, oe: u64) {
            self.by_contig.entry(c.to_string()).or_default().push((s, e, oc.to_string(), os, oe));
            let m = self.max_len.entry(c.to_string()).or_insert(0);
            *m = (*m).max(e.saturating_sub(s));
        }
        pub fn n_pairs(&self) -> usize {
            self.by_contig.values().map(|v| v.len()).sum::<usize>() / 2
        }
        /// Pairs whose THIS side overlaps `[s, e)` on `contig`.
        /// ⭐ §6fw guard: are two intervals of ONE contig linked to each other by a single duplication pair?
        /// A cross-copy mis-chain — the artefact the read-through object is accused of being — requires the
        /// donor flank and the acceptor flank to be copies of one another, so an aligner has something to jump
        /// between. When no pair links them, that mechanism cannot produce the junction. Pairs are stored in
        /// both directions (`from_bed_str`/`from_paf_str`), so one query on `a` suffices.
        pub fn links(&self, contig: &str, a: (u64, u64), b: (u64, u64)) -> bool {
            self.overlapping(contig, a.0, a.1).any(|(_, _, oc, os, oe)| oc == contig && *oe > b.0 && *os < b.1)
        }

        /// ⭐ §6iz proposal #1 (o1_ledger §6j4): a family-free generalization of `refine_cluster_cores_with`'s
        /// core rule for a SINGLE genomic interval with no known co-member set — de novo mode's fresh
        /// boundary-discovery case. `refine_cluster_cores_with`'s depth floor is `(n-1)/2` OTHER KNOWN
        /// members, which is meaningless with n=1. Here the floor is `min_partners`, an ABSOLUTE minimum of
        /// DISTINCT partner loci corroborating the same position — the same "one hit proves nothing, two+
        /// agreeing sources do" principle already used by `family_detect::CNT_MIN` (a k-mer must be owned by
        /// `>= CNT_MIN` distinct reps to be family-informative). `window_slop` widens the query around
        /// `[s, e)` before looking for SD evidence, since a truncated/fragmentary de novo node's own span can
        /// sit well inside the true duplicated block without any SD pair touching its exact, too-small
        /// boundaries. Returns the merged hull of the depth-passing segments, or `None` when no segment
        /// clears `min_partners` (including when there is no SD evidence at all near this locus).
        pub fn single_span_core(
            &self,
            contig: &str,
            s: u64,
            e: u64,
            window_slop: u64,
            min_partners: usize,
        ) -> Option<(u64, u64)> {
            let qs = s.saturating_sub(window_slop);
            let qe = e.saturating_add(window_slop);
            let mut hits: Vec<(u64, u64, String, u64, u64)> = self
                .overlapping(contig, qs, qe)
                .map(|(ss, se, oc, os, oe)| (qs.max(*ss), qe.min(*se), oc.clone(), *os, *oe))
                .collect();
            if hits.is_empty() {
                return None;
            }
            // Distinct partner clusters: merge overlapping partner intervals on the SAME partner contig so
            // several SD rows hitting the same sibling copy corroborate as ONE source, not one per row.
            hits.sort_by(|a, b| (a.2.clone(), a.3, a.4).cmp(&(b.2.clone(), b.3, b.4)));
            let mut cluster_id: Vec<usize> = Vec::with_capacity(hits.len());
            let mut next_id = 0usize;
            let mut open: Option<(String, u64, u64, usize)> = None;
            for h in &hits {
                match &mut open {
                    Some((oc, _os, oe, id)) if *oc == h.2 && h.3 <= *oe => {
                        *oe = (*oe).max(h.4);
                        cluster_id.push(*id);
                    }
                    _ => {
                        let id = next_id;
                        next_id += 1;
                        cluster_id.push(id);
                        open = Some((h.2.clone(), h.3, h.4, id));
                    }
                }
            }
            // Depth sweep over the QUERY window: count DISTINCT partner-cluster ids covering each position.
            let mut events: Vec<(u64, i32, usize)> = Vec::with_capacity(hits.len() * 2);
            for (h, &cid) in hits.iter().zip(cluster_id.iter()) {
                events.push((h.0, 1, cid));
                events.push((h.1, -1, cid));
            }
            events.sort_unstable_by_key(|x| x.0);
            let mut count: BTreeMap<usize, i32> = BTreeMap::new();
            let mut depth = 0usize;
            let mut last: Option<u64> = None;
            let mut segs: Vec<(u64, u64)> = Vec::new();
            for (pos, d, cid) in events {
                if let Some(l) = last {
                    if pos > l && depth >= min_partners {
                        if let Some(t) = segs.last_mut().filter(|t| t.1 == l) {
                            t.1 = pos;
                        } else {
                            segs.push((l, pos));
                        }
                    }
                }
                let c = count.entry(cid).or_insert(0);
                if d > 0 {
                    if *c == 0 {
                        depth += 1;
                    }
                    *c += 1;
                } else {
                    *c -= 1;
                    if *c == 0 {
                        depth -= 1;
                    }
                }
                last = Some(pos);
            }
            if segs.is_empty() {
                return None;
            }
            Some((segs[0].0, segs.last().unwrap().1))
        }

        fn overlapping(&self, contig: &str, s: u64, e: u64) -> impl Iterator<Item = &(u64, u64, String, u64, u64)> {
            let v = self.by_contig.get(contig).map(|v| v.as_slice()).unwrap_or(&[]);
            let lo = s.saturating_sub(self.max_len.get(contig).copied().unwrap_or(0));
            let i = v.partition_point(|x| x.0 < lo);
            v[i..].iter().take_while(move |x| x.0 < e).filter(move |x| x.1 > s)
        }
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum CoreStatus {
        /// The cluster failed the depth gate (no SD evidence): every member is left exactly as it was.
        Untouched,
        /// `core_bp >= span/2`: a full or partial copy of the duplicon.
        KeptFull,
        /// `core_bp < span/2` but `>= median_core/2`: a chimeric model carrying a full core; node := core hull.
        KeptTrimmed,
        /// Neither: a record carrying at most a sliver of the family's core.
        Dropped,
    }

    #[derive(Debug, Clone)]
    pub struct CoreRecord {
        pub member: GeneKey,
        pub gate_passed: bool,
        /// Maximum over the member's bases of the number of distinct other members linked there.
        pub max_depth: usize,
        pub core_bp: u64,
        pub span: u64,
        pub median_core: u64,
        pub status: CoreStatus,
        /// 1-based inclusive hull of the core segments, when a core exists.
        pub hull: Option<(u64, u64)>,
    }

    /// ⭐ DUPLICATION BLOCKS (user request 2026-09-05): the family is defined by a shared CORE; the block is the
    /// larger unit of co-duplication that SEDEF sees. Two core hulls belong to one block when some single SEDEF
    /// pair overlaps both (a block-level SD record spans several modules — LCR16a + LCR16u — so the pair links
    /// hulls of DIFFERENT families; a module-level record links hulls of one family). Union-find over every hull
    /// in the catalog; returns one block index per hull (dense, in first-seen order). Families whose hulls share a
    /// block are the same duplication block; a family whose hulls span several blocks is a family that
    /// transposed. `hulls` are (contig, 1-based start, 1-based end). Also returns every DIRECT link (i, j),
    /// i < j, between two hulls joined by one pair — the block is the transitive closure of these, and on 16p the
    /// closure is the whole LCR16 network (30 clusters in one block), so the direct partners are the finer,
    /// quotable relation ("NPIP's hulls are directly SD-linked to LCR16u's").
    pub fn sd_blocks(hulls: &[(String, u64, u64)], sd: &SdPairs) -> (Vec<usize>, Vec<(usize, usize)>) {
        let n = hulls.len();
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(p: &mut Vec<usize>, mut x: usize) -> usize {
            while p[x] != x {
                p[x] = p[p[x]];
                x = p[x];
            }
            x
        }
        // per contig, hull indices sorted by start
        let mut by_contig: BTreeMap<&str, Vec<(u64, u64, usize)>> = BTreeMap::new();
        for (i, (c, s, e)) in hulls.iter().enumerate() {
            by_contig.entry(c.as_str()).or_default().push((*s, *e, i));
        }
        for v in by_contig.values_mut() {
            v.sort_unstable();
        }
        let hulls_over = |c: &str, s: u64, e: u64| -> Vec<usize> {
            by_contig
                .get(c)
                .map(|v| v.iter().filter(|&&(hs, he, _)| hs <= e && s <= he).map(|&(_, _, i)| i).collect())
                .unwrap_or_default()
        };
        let mut links: BTreeSet<(usize, usize)> = BTreeSet::new();
        for (i, (c, s, e)) in hulls.iter().enumerate() {
            for &(ps, pe, ref oc, os, oe) in sd.overlapping(c, *s, *e) {
                let _ = (ps, pe);
                for j in hulls_over(oc, os, oe) {
                    if i != j {
                        links.insert((i.min(j), i.max(j)));
                    }
                    let (a, b) = (find(&mut parent, i), find(&mut parent, j));
                    if a != b {
                        parent[a.max(b)] = a.min(b);
                    }
                }
            }
        }
        let mut id: BTreeMap<usize, usize> = BTreeMap::new();
        let blocks = (0..n)
            .map(|i| {
                let r = find(&mut parent, i);
                let next = id.len();
                *id.entry(r).or_insert(next)
            })
            .collect();
        (blocks, links.into_iter().collect())
    }

    /// Apply the pre-registered core rule to one cluster. Members are GFF 1-based inclusive.
    /// `inclusive` (§6ft polish 2): the majority counts the locus itself — a member's core is the part shared with
    /// at least half of the FAMILY (depth + 1 ≥ n/2) instead of half of the OTHER members ((n − 1)/2). Only the
    /// boundary case moves (NPIP: a 7-kb fragment shared with 15 of 31 others).
    pub fn refine_cluster_cores_with(members: &[GeneKey], sd: &SdPairs, inclusive: bool) -> Vec<CoreRecord> {
        let n = members.len();
        if n < 2 {
            return members
                .iter()
                .map(|m| CoreRecord {
                    member: m.clone(),
                    gate_passed: false,
                    max_depth: 0,
                    core_bp: 0,
                    span: m.2.saturating_sub(m.1) + 1,
                    median_core: 0,
                    status: CoreStatus::Untouched,
                    hull: None,
                })
                .collect();
        }
        let half_others = if inclusive { n as f64 / 2.0 - 1.0 } else { (n - 1) as f64 / 2.0 };
        // Per member: depth profile from every SD pair whose one side overlaps the member and whose other
        // side overlaps a DIFFERENT member; depth counts distinct partner members.
        let mut recs: Vec<(usize, u64, Vec<(u64, u64)>)> = Vec::with_capacity(n); // (max_depth, core_bp, core segments)
        for (ri, r) in members.iter().enumerate() {
            let (rs, re) = (r.1.saturating_sub(1), r.2);
            let mut events: Vec<(u64, i32, usize)> = Vec::new();
            for (s, e, oc, os, oe) in sd.overlapping(&r.0, rs, re) {
                for (oi, o) in members.iter().enumerate() {
                    if oi == ri || o.0 != *oc {
                        continue;
                    }
                    if *oe > o.1.saturating_sub(1) && *os < o.2 {
                        events.push((rs.max(*s), 1, oi));
                        events.push((re.min(*e), -1, oi));
                    }
                }
            }
            events.sort_unstable();
            let mut count = vec![0i32; n];
            let mut depth = 0usize;
            let mut max_depth = 0usize;
            let mut last: Option<u64> = None;
            let mut segs: Vec<(u64, u64)> = Vec::new();
            for (pos, d, oi) in events {
                if let Some(l) = last {
                    if pos > l && depth as f64 >= half_others {
                        if let Some(t) = segs.last_mut().filter(|t| t.1 == l) {
                            t.1 = pos;
                        } else {
                            segs.push((l, pos));
                        }
                    }
                }
                if d > 0 {
                    if count[oi] == 0 {
                        depth += 1;
                    }
                    count[oi] += 1;
                } else {
                    count[oi] -= 1;
                    if count[oi] == 0 {
                        depth -= 1;
                    }
                }
                max_depth = max_depth.max(depth);
                last = Some(pos);
            }
            let core_bp = segs.iter().map(|(a, b)| b - a).sum::<u64>();
            recs.push((max_depth, core_bp, segs));
        }
        let mut ratios: Vec<f64> = recs.iter().map(|(md, _, _)| *md as f64 / (n - 1) as f64).collect();
        ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let gate = median_f(&ratios) >= 0.5;
        let mut cores: Vec<u64> = recs.iter().map(|(_, c, _)| *c).collect();
        cores.sort_unstable();
        let median_core = cores[cores.len() / 2];
        members
            .iter()
            .zip(recs)
            .map(|(m, (max_depth, core_bp, segs))| {
                let span = m.2.saturating_sub(m.1) + 1;
                let hull = if segs.is_empty() { None } else { Some((segs[0].0 + 1, segs.last().unwrap().1)) };
                let status = if !gate {
                    CoreStatus::Untouched
                } else if core_bp as f64 >= span as f64 / 2.0 {
                    CoreStatus::KeptFull
                } else if core_bp as f64 >= median_core as f64 / 2.0 && core_bp > 0 {
                    CoreStatus::KeptTrimmed
                } else {
                    CoreStatus::Dropped
                };
                CoreRecord { member: m.clone(), gate_passed: gate, max_depth, core_bp, span, median_core, status, hull }
            })
            .collect()
    }

    fn median_f(sorted: &[f64]) -> f64 {
        if sorted.is_empty() {
            return 0.0;
        }
        sorted[sorted.len() / 2]
    }

    /// Sort and merge half-open intervals so overlapping alignment records are not double-counted.
    fn merge_intervals(mut v: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
        if v.is_empty() {
            return v;
        }
        v.sort_unstable();
        let mut out = vec![v[0]];
        for &(s, e) in &v[1..] {
            let last = out.last_mut().expect("non-empty");
            if s <= last.1 {
                last.1 = last.1.max(e);
            } else {
                out.push((s, e));
            }
        }
        out
    }

    /// Bases of the local, 0-based half-open interval `[s, e)` that fall inside an exon.
    ///
    /// `blocks` are ABSOLUTE 0-based half-open merged exon intervals; `gene_start` is the GFF **1-based**
    /// gene start, so an absolute base `b` sits at local `b - (gene_start - 1)`. ⚠ Getting that conversion
    /// wrong is the §6dd bug: a `start - 1` key joined 0/4,477 and silently fell back to the span.
    fn exonic_bases_in(blocks: &[(u64, u64)], gene_start: u64, s: u64, e: u64) -> u64 {
        let off = gene_start.saturating_sub(1);
        let mut acc = 0u64;
        for &(bs, be) in blocks {
            let (lo, hi) = (bs.saturating_sub(off).max(s), be.saturating_sub(off).min(e));
            if hi > lo {
                acc += hi - lo;
            }
        }
        acc
    }

    /// Markov clustering over a sparse symmetric graph.
    ///
    /// ⭐ **INFLATION SITS ON A PLATEAU FOR NPIP — AND ONLY FOR NPIP** (§6dd, caveated §6dx). `I = 2.8..3.2`
    /// gives identical NPIP recovery (31/31) — but the real 48-member tandem clique has best-Jaccard **0.00 at
    /// I=3.2**, and low-cohesion SD-embedded subfamilies are stable (1.00) at every inflation. Do not generalise
    /// the plateau beyond NPIP. Measured on the pilot, `I = 2.8..3.2` gives
    /// identical NPIP recovery (31/31), identical principal-cluster size (43) and identical truth
    /// concentration (26), with **zero clusters over 100 members**; the cliff is at 3.6 (NPIP 31 -> 18). A
    /// constant on a plateau is a choice; a constant on a cliff is indefensible. The default is the middle.
    /// ⚠ §6ec (2026-09-04): that cliff is a PRUNE artefact. With the absolute post-inflation `prune` (1e-5) a
    /// near-uniform clique of n nodes empties once n+1 > 10^(5/I) (61 at I=2.8, 37 at 3.2); at prune 1e-9
    /// NPIP is 44/44 from I=2.0 to 4.0 and the 84+22-copy tandem array, dissolved genome-wide at 1e-5,
    /// clusters. The default is kept for byte-identity; `mcl_families --prune` exposes it.
    ///
    /// Determinism: the graph is a `BTreeMap`, columns are normalised in index order, and ties keep the
    /// lowest index — the same PAF and parameters always yield the same clustering.
    pub fn mcl(g: &HomologyGraph, inflation: f64, max_iter: usize, prune: f64) -> Vec<Vec<usize>> {
        let n = g.n_nodes();
        if n == 0 {
            return Vec::new();
        }
        // Column-major sparse matrix with self-loops, as MCL requires.
        let mut col: Vec<BTreeMap<usize, f64>> = vec![BTreeMap::new(); n];
        for (&(a, b), &w) in &g.edges {
            col[a].insert(b, w);
            col[b].insert(a, w);
        }
        for (i, c) in col.iter_mut().enumerate() {
            c.insert(i, 1.0);
        }
        normalize(&mut col);

        for _ in 0..max_iter {
            let next = expand(&col);
            let mut next = next;
            for c in next.iter_mut() {
                for v in c.values_mut() {
                    *v = v.powf(inflation);
                }
                c.retain(|_, v| *v >= prune);
            }
            normalize(&mut next);
            let done = converged(&col, &next);
            col = next;
            if done {
                break;
            }
        }
        attractors(&col, n)
    }

    fn normalize(col: &mut [BTreeMap<usize, f64>]) {
        for c in col.iter_mut() {
            let s: f64 = c.values().sum();
            if s > 0.0 {
                for v in c.values_mut() {
                    *v /= s;
                }
            }
        }
    }

    /// One MCL expansion step: `M <- M * M`, column-major.
    fn expand(col: &[BTreeMap<usize, f64>]) -> Vec<BTreeMap<usize, f64>> {
        col.iter()
            .map(|c| {
                let mut out: BTreeMap<usize, f64> = BTreeMap::new();
                for (&k, &wk) in c {
                    for (&i, &wi) in &col[k] {
                        *out.entry(i).or_insert(0.0) += wi * wk;
                    }
                }
                out
            })
            .collect()
    }

    fn converged(a: &[BTreeMap<usize, f64>], b: &[BTreeMap<usize, f64>]) -> bool {
        a.iter().zip(b.iter()).all(|(x, y)| {
            x.len() == y.len()
                && x.iter().zip(y.iter()).all(|((i, u), (j, v))| i == j && (u - v).abs() < 1e-7)
        })
    }

    /// Read clusters off the converged matrix: every column's heaviest row is its attractor, and columns
    /// sharing an attractor (transitively) are one cluster. Ties keep the lowest index.
    fn attractors(col: &[BTreeMap<usize, f64>], n: usize) -> Vec<Vec<usize>> {
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(p: &mut Vec<usize>, mut x: usize) -> usize {
            while p[x] != x {
                p[x] = p[p[x]];
                x = p[x];
            }
            x
        }
        for (j, c) in col.iter().enumerate() {
            let Some((&r, _)) = c.iter().max_by(|a, b| {
                a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal).then(b.0.cmp(a.0))
            }) else {
                continue;
            };
            let (ra, rb) = (find(&mut parent, j), find(&mut parent, r));
            if ra != rb {
                parent[ra.max(rb)] = ra.min(rb);
            }
        }
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for i in 0..n {
            let r = find(&mut parent, i);
            groups.entry(r).or_default().push(i);
        }
        let mut out: Vec<Vec<usize>> = groups.into_values().collect();
        out.sort_by(|a, b| b.len().cmp(&a.len()).then(a[0].cmp(&b[0])));
        out
    }

    /// A cluster as reported: member genes plus the two statistics that decide whether it is a family.
    #[derive(Debug, Clone)]
    pub struct Cluster {
        pub members: Vec<GeneKey>,
        /// Fraction of the possible member pairs that carry a homology edge. ⭐ A repeat clique is a PERFECT
        /// clique (median **1.000**); a real family is not (**0.700**) — §6de, size-controlled.
        pub density: f64,
        /// ⭐ **COHESION CERTIFICATE**: `e_in / (e_in + e_out)` — the share of this cluster's incident edges
        /// that stay inside it. ⚠ Distinct from `density`, which sees only internal pairs and is therefore
        /// blind to a slice cut out of a hairball.
        ///
        /// Calibrated (§6du, admitted edges at `min_exonic_bp = 1`): adjudicated REAL families measure
        /// **0.970** (NPIP, n=43) and **0.923** (n=27); the adjudicated repeat clique measures **0.283** and
        /// the unadjudicated 101-member cluster **0.557**. Over the pilot, **73/142** clusters sit in a clean
        /// mode at >= 0.9, while **20/46** clusters with n >= 5 fall below 0.75, covering 253/620 members.
        /// ⚠ Calibrated on 3 real + 2 artefact clusters — REPORTED, never used as a filter.
        pub frac_in: f64,
        /// Fraction of members with RNA read support, when a BAM was supplied. `None` means NOT MEASURED —
        /// never conflate it with 0.0, which is the repeat-clique signature.
        pub corroborated: Option<f64>,
    }

    /// Assemble reportable clusters. `min_size` drops singletons (and, at 3, the pairs that carry no
    /// density signal). `corroborated_members` is the caller's RNA verdict per gene, or `None` if no BAM.
    /// ⭐ Fold overlapping annotation records into loci WITHIN each MCL part (§6ey). The fold-first order
    /// (`loci_from_exon_blocks` over the whole annotation, then representative-only edges) loses every record
    /// that overlaps a record of a DIFFERENT family on exon bases — a pseudogene inside a host gene's exon, a
    /// head-to-head pair, an antisense lncRNA, a readthrough model — because the locus keeps ONE representative
    /// (the longest exon union) and discards the others' homology (Soto: ANAPC1P1 → CD8B's locus, PMS2P4 →
    /// SPDYE21, FAM72A → SRGAP2, LRRC37A → ARL17B; 19 of 43 annotated misses). Attribution semantics keep the
    /// evidence but rebuild the duplication BLOCK (§6eg). Folding AFTER clustering dissolves the dilemma: two
    /// records are one locus only if they overlap on exon bases AND sit in the same sequence-homology cluster —
    /// coordinates and sequence are two different criteria, so the fold is not circular. Returns the parts with
    /// each locus reduced to its representative (longest exon union, `loci_from_exon_blocks`' rule) and the map
    /// annotation → representative for every record folded away. Records without exon blocks are never folded.
    pub fn fold_parts_into_loci(
        g: &HomologyGraph,
        parts: &[Vec<usize>],
        exon_blocks: &BTreeMap<GeneKey, Vec<(u64, u64)>>,
    ) -> (Vec<Vec<usize>>, LocusMap) {
        let mut all = LocusMap::default();
        let mut out = Vec::with_capacity(parts.len());
        for part in parts {
            let sub: BTreeMap<GeneKey, Vec<(u64, u64)>> = part
                .iter()
                .filter_map(|&i| exon_blocks.get(&g.genes[i]).map(|b| (g.genes[i].clone(), b.clone())))
                .collect();
            let m = loci_from_exon_blocks(&sub);
            let kept: Vec<usize> = part
                .iter()
                .copied()
                .filter(|&i| m.rep_of.get(&g.genes[i]).map_or(true, |r| *r == g.genes[i]))
                .collect();
            for (k, r) in &m.rep_of {
                if k != r {
                    all.rep_of.insert(k.clone(), r.clone());
                }
            }
            all.n_multi += m.n_multi;
            out.push(kept);
        }
        (out, all)
    }

    pub fn build_clusters(
        g: &HomologyGraph,
        parts: &[Vec<usize>],
        min_size: usize,
        corroborated: Option<&dyn Fn(&GeneKey) -> bool>,
    ) -> Vec<Cluster> {
        parts
            .iter()
            .filter(|p| p.len() >= min_size)
            .map(|p| {
                let n = p.len();
                let possible = n * (n - 1) / 2;
                let mut have = 0usize;
                for (a, x) in p.iter().enumerate() {
                    for y in p.iter().skip(a + 1) {
                        let k = if x < y { (*x, *y) } else { (*y, *x) };
                        if g.edges.contains_key(&k) {
                            have += 1;
                        }
                    }
                }
                // Cohesion: count every edge incident on a member, then split inside/outside. `density`
                // cannot see e_out, which is exactly how a hairball slice passes as a family.
                let inside: BTreeSet<usize> = p.iter().copied().collect();
                let (mut e_in, mut e_out) = (0usize, 0usize);
                for &(x, y) in g.edges.keys() {
                    match (inside.contains(&x), inside.contains(&y)) {
                        (true, true) => e_in += 1,
                        (true, false) | (false, true) => e_out += 1,
                        _ => {}
                    }
                }
                let members: Vec<GeneKey> = p.iter().map(|&i| g.genes[i].clone()).collect();
                let corr = corroborated.map(|f| {
                    members.iter().filter(|m| f(m)).count() as f64 / n.max(1) as f64
                });
                Cluster {
                    members,
                    density: if possible == 0 { 0.0 } else { have as f64 / possible as f64 },
                    frac_in: if e_in + e_out == 0 { 0.0 } else { e_in as f64 / (e_in + e_out) as f64 },
                    corroborated: corr,
                }
            })
            .collect()
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn paf_line(q: &str, ql: u64, qs: u64, qe: u64, t: &str, tl: u64, ts: u64, te: u64, nm: u64, bl: u64) -> String {
            format!("{q}\t{ql}\t{qs}\t{qe}\t+\t{t}\t{tl}\t{ts}\t{te}\t{nm}\t{bl}\t60")
        }

        /// §6fo: a PAF record between two annotated loci becomes one SD-like pair in genomic coordinates (both
        /// directions), offsets applied to the loci's 1-based starts.
        #[test]
        /// The §6fw guard: `links` is true only when one duplication pair holds BOTH intervals, in either
        /// stored direction; unrelated flanks on the same contig, and pairs to another contig, are not links.
        #[test]
        fn links_is_true_only_when_one_pair_holds_both_flanks() {
            // one pair c1:1000-2000 <-> c1:50000-51000, and one c1:1000-2000 <-> c2:7000-8000
            let bed = "c1\t1000\t2000\tc1\t50000\t51000\nc1\t1000\t2000\tc2\t7000\t8000\n";
            let sd = SdPairs::from_bed_str(bed);
            assert!(sd.links("c1", (1500, 1600), (50100, 50200)));
            assert!(sd.links("c1", (50100, 50200), (1500, 1600)), "pairs are stored both ways");
            // the second flank is nowhere near the partner interval
            assert!(!sd.links("c1", (1500, 1600), (90000, 90100)));
            // the only pair reaching c2 does not link two c1 intervals
            assert!(!sd.links("c2", (7100, 7200), (7300, 7400)));
            // neither flank is in a pair at all
            assert!(!sd.links("c1", (300000, 300100), (400000, 400100)));
        }

        /// §6iz proposal #1 / §6j4: two distinct partner clusters overlapping the widened query window agree
        /// on [9900,10600) — the segment where depth>=2 — and disagree (depth 1 only) outside it.
        #[test]
        fn single_span_core_admits_the_segment_two_distinct_partners_agree_on() {
            let bed = "c1\t9800\t10600\tc2\t50000\t50800\nc1\t9900\t10700\tc3\t70000\t70800\n";
            let sd = SdPairs::from_bed_str(bed);
            let core = sd.single_span_core("c1", 10000, 10500, 1000, 2);
            assert_eq!(core, Some((9900, 10600)));
        }

        /// A single corroborating partner never clears an `min_partners=2` floor -- one hit proves nothing.
        #[test]
        fn single_span_core_none_when_only_one_partner_backs_the_locus() {
            let bed = "c1\t9800\t10600\tc2\t50000\t50800\n";
            let sd = SdPairs::from_bed_str(bed);
            assert_eq!(sd.single_span_core("c1", 10000, 10500, 1000, 2), None);
        }

        /// No SD evidence anywhere near the locus at all.
        #[test]
        fn single_span_core_none_when_no_sd_evidence_nearby() {
            let sd = SdPairs::from_bed_str("");
            assert_eq!(sd.single_span_core("c1", 10000, 10500, 1000, 2), None);
        }

        /// Several SD rows hitting the SAME sibling copy (overlapping partner intervals on one partner
        /// contig) corroborate as ONE source, not one per row -- they must NOT by themselves clear a
        /// `min_partners=2` floor.
        #[test]
        fn single_span_core_merges_overlapping_same_partner_rows_into_one_source() {
            let bed = "c1\t9800\t10600\tc2\t50000\t50800\nc1\t9850\t10550\tc2\t50050\t50850\n";
            let sd = SdPairs::from_bed_str(bed);
            assert_eq!(sd.single_span_core("c1", 10000, 10500, 1000, 2), None);
        }

        /// §6j4 real-data check (ignored: reads a real, machine-local file; not part of the CI suite). Verifies
        /// `single_span_core` against the REAL gorilla SD file on the 5 real, currently-wrong-sized NPIP de
        /// novo nodes measured in the ledger, reproducing (not just asserting-in-the-abstract) the real
        /// coverage/ratio improvement reported there.
        #[test]
        #[ignore]
        fn single_span_core_real_gorilla_npip_wrong_sized_nodes() {
            let text = std::fs::read_to_string("/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed")
                .expect("real gorilla SD file must exist for this ad-hoc check");
            let sd = SdPairs::from_bed_str(&text);
            // (truth_start, truth_end, node_start, node_end) on their real contigs, from the real ggo.nodes.tsv
            // dump / npip31.regions oracle (o1_ledger §6j4).
            let cases: &[(&str, u64, u64, u64, u64)] = &[
                ("NC_073242.2", 35184502, 35267752, 35228011, 35230492),
                ("NC_073244.2", 21077140, 21105296, 21075545, 21078820),
                ("NC_073242.2", 28995101, 29024611, 28995202, 28996824),
                ("NC_073242.2", 28300720, 28325984, 28319700, 28337466),
            ];
            for &(contig, ts, te, ns, ne) in cases {
                let node_cov = (ne.min(te) as i64 - ns.max(ts) as i64).max(0) as f64 / (te - ts) as f64;
                let node_ratio = (ne - ns) as f64 / (te - ts) as f64;
                match sd.single_span_core(contig, ns, ne, 15_000, 2) {
                    Some((cs, ce)) => {
                        let cov = (ce.min(te) as i64 - cs.max(ts) as i64).max(0) as f64 / (te - ts) as f64;
                        let ratio = (ce - cs) as f64 / (te - ts) as f64;
                        println!(
                            "{contig}:{ts}-{te}  raw_node cov={node_cov:.3} ratio={node_ratio:.3}  -> sd_core {cs}-{ce} cov={cov:.3} ratio={ratio:.3}"
                        );
                        assert!(cov > node_cov, "SD-core span should improve coverage vs the raw node span");
                    }
                    None => println!(
                        "{contig}:{ts}-{te}  raw_node cov={node_cov:.3} ratio={node_ratio:.3}  -> sd_core: NONE (no SD evidence)"
                    ),
                }
            }
        }

        fn paf_records_become_pairs_in_genomic_coordinates() {
            let paf = "c1:1001-2000\t1000\t100\t600\t+\tc2:5001-7000\t2000\t1500\t2000\t480\t500\t60\n\
                       c1:1001-2000\t1000\t0\t0\t+\tc2:5001-7000\t2000\t0\t0\t0\t0\t0\n";
            let sd = SdPairs::from_paf_str(paf);
            assert_eq!(sd.n_pairs(), 1);
            let hits: Vec<_> = sd.overlapping("c1", 1000, 1200).cloned().collect();
            assert_eq!(hits, vec![(1100, 1600, "c2".to_string(), 6500, 7000)]);
            let back: Vec<_> = sd.overlapping("c2", 6600, 6700).cloned().collect();
            assert_eq!(back, vec![(6500, 7000, "c1".to_string(), 1100, 1600)]);
        }

        #[test]
        fn gene_key_parses_gff_one_based_coordinates_verbatim() {
            // ⚠ The headers carry GFF 1-based coordinates AS WRITTEN. A `start - 1` key joins nothing and
            // silently falls back to the span denominator (§6dd) — pin the convention.
            assert_eq!(
                parse_gene_key("NC_073241.2:31346-41669"),
                Some(("NC_073241.2".to_string(), 31346, 41669))
            );
            // A contig name containing ':' must still split on the LAST one.
            assert_eq!(parse_gene_key("a:b:10-20"), Some(("a:b".to_string(), 10, 20)));
            assert_eq!(parse_gene_key("no-colon"), None);
        }

        #[test]
        fn coverage_denominator_is_exonic_not_span() {
            // Two genes with 10 kb spans but only 1 kb of exon each, sharing a 600 bp alignment.
            // Against the SPAN: 600/10000 = 0.06, rejected. Against EXON: 600/1000 = 0.60, admitted.
            // This is the §6dc defect that excluded 62.8% of genes.
            let paf = paf_line("c:1-10000", 10_000, 0, 600, "c:20001-30000", 10_000, 0, 600, 570, 600);
            let p = GraphParams::default();

            let empty = BTreeMap::new();
            let span_graph = graph_from_paf(&paf, &empty, &BTreeMap::new(), &p);
            assert_eq!(span_graph.n_edges(), 0, "span denominator must reject this pair");
            assert_eq!(span_graph.missing_exonic, 2, "both genes fell back, and that must be counted");

            let mut exonic = BTreeMap::new();
            exonic.insert(("c".to_string(), 1, 10_000), 1_000);
            exonic.insert(("c".to_string(), 20_001, 30_000), 1_000);
            let exon_graph = graph_from_paf(&paf, &exonic, &BTreeMap::new(), &p);
            assert_eq!(exon_graph.n_edges(), 1, "exonic denominator must admit it");
            assert_eq!(exon_graph.missing_exonic, 0);
        }

        #[test]
        fn weight_uses_the_longer_side_so_a_short_fragment_cannot_drive_clustering() {
            // A 300 bp fragment fully covered by a 300 bp alignment to a 10 kb gene. On the SHORTER side that
            // is coverage 1.0; on the LONGER side 0.03, below the floor. Every one of the 22 adjudicated NPIP
            // false merges was exactly this shape — a shared Alu (§6cr).
            let paf = paf_line("c:1-300", 300, 0, 300, "c:1001-11000", 10_000, 0, 300, 290, 300);
            let mut exonic = BTreeMap::new();
            exonic.insert(("c".to_string(), 1, 300), 300);
            exonic.insert(("c".to_string(), 1001, 11_000), 10_000);
            let g = graph_from_paf(&paf, &exonic, &BTreeMap::new(), &GraphParams::default());
            assert_eq!(g.n_edges(), 0, "a fragment-sized alignment must not become an edge");
        }

        #[test]
        fn mcl_splits_two_cliques_joined_by_a_single_bridge() {
            // The superfamily failure this module exists to prevent: transitive closure chains subfamilies
            // through a hub (register:474 — 145- and 114-gene superfamilies). MCL must cut the bridge.
            let mut g = HomologyGraph::default();
            for i in 0..6 {
                g.genes.push(("c".to_string(), i as u64 * 100 + 1, i as u64 * 100 + 50));
            }
            for (a, b) in [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)] {
                g.edges.insert((a, b), 0.95);
            }
            g.edges.insert((2, 3), 0.72); // the single weak bridge

            let parts = mcl(&g, 2.8, 100, 1e-5);
            let big: Vec<&Vec<usize>> = parts.iter().filter(|p| p.len() >= 2).collect();
            assert_eq!(big.len(), 2, "expected the two cliques to separate, got {parts:?}");
            assert!(big.iter().all(|p| p.len() == 3), "each clique keeps its three members: {parts:?}");
        }

        #[test]
        fn absolute_prune_empties_large_uniform_cliques_and_a_size_safe_prune_does_not() {
            // §6ec: after inflation every entry of a near-uniform n-clique is ~(1/(n+1))^I, so an absolute
            // prune p empties the whole column once n+1 > p^(-1/I) — 61 nodes at I=2.8, p=1e-5. That is how
            // the anchored 84+22-copy tandem array dissolved genome-wide. A 100-clique must survive at 1e-9.
            let mut g = HomologyGraph::default();
            for i in 0..100 {
                g.genes.push(("c".to_string(), i as u64 * 100 + 1, i as u64 * 100 + 50));
            }
            // Weights vary deterministically in [0.90, 0.99]: a perfectly uniform clique is a numerical
            // knife-edge (every column's self entry ties its neighbours) and is not what a family looks like.
            for a in 0..100 {
                for b in (a + 1)..100 {
                    g.edges.insert((a, b), 0.90 + 0.09 * (((a * 7 + b * 13) % 17) as f64 / 16.0));
                }
            }
            let old = mcl(&g, 2.8, 100, 1e-5);
            let largest_old = old.iter().map(|p| p.len()).max().unwrap_or(0);
            assert!(largest_old < 50, "1e-5 must shatter a 100-clique at I=2.8 (largest {largest_old}); if it no longer does, the prune semantics changed and §6ec must be re-measured");
            let new = mcl(&g, 2.8, 100, 1e-9);
            let largest_new = new.iter().map(|p| p.len()).max().unwrap_or(0);
            assert!(largest_new >= 95, "1e-9 must keep the 100-clique (largest {largest_new}): {new:?}");
        }

        #[test]
        fn loci_are_exon_overlap_components_not_genomic_overlap() {
            let k = |s: u64, e: u64| ("c".to_string(), s, e);
            let mut bl: BTreeMap<GeneKey, Vec<(u64, u64)>> = BTreeMap::new();
            bl.insert(k(100, 1000), vec![(100, 200), (900, 1000)]); // host: two exons, big intron
            bl.insert(k(300, 400), vec![(300, 400)]); // INTRONIC gene inside the host: a distinct locus
            bl.insert(k(150, 250), vec![(150, 250)]); // overlaps the host's first exon: same locus
            bl.insert(k(950, 1500), vec![(950, 1050), (1400, 1500)]); // overlaps the host's last exon: same locus, longest exon-union? 200 vs host 200 -> tie -> lowest start = host
            bl.insert(k(1450, 1600), vec![(1450, 1600)]); // overlaps the previous only: chained in
            bl.insert(k(2000, 2500), vec![(2000, 2500)]); // separate
            bl.insert(("d".to_string(), 100, 1000), vec![(100, 200)]); // other contig
            let m = loci_from_exon_blocks(&bl);
            assert_eq!(m.n_multi, 1, "{:?}", m.rep_of);
            assert_eq!(m.n_merged(), 3);
            assert!(m.is_representative(&k(300, 400)), "an intronic gene is its own locus");
            assert!(m.is_representative(&k(2000, 2500)));
            assert!(m.is_representative(&("d".to_string(), 100, 1000)));
            // exon-union lengths: host 200, (150,250) 100, (950,1500) 200, (1450,1600) 150 -> tie host vs (950,1500) -> lowest start wins
            for g in [k(150, 250), k(950, 1500), k(1450, 1600)] {
                assert_eq!(m.representative(&g), &k(100, 1000), "{g:?}");
            }
        }

        #[test]
        fn a_nested_annotation_is_not_a_second_node() {
            // Host H (c:1001-3000) contains lncRNA N (c:1201-2200). H aligns to paralogue P, N to paralogue Q
            // (each at full coverage). Without the locus map the family counts H's locus twice: 4 nodes, 2 edges.
            // REPRESENTATIVE-ONLY (default, §6eq): N's record carries no evidence -> 2 nodes {H, P}, 1 edge.
            // ATTRIBUTION (explicit, §6eg): N's edge becomes H's -> 3 nodes {H, P, Q}, 2 edges.
            let mut paf = paf_line("c:1001-3000", 2000, 0, 2000, "c:9001-11000", 2000, 0, 2000, 1900, 2000);
            paf.push('\n');
            paf.push_str(&paf_line("c:1201-2200", 1000, 0, 1000, "c:20001-21000", 1000, 0, 1000, 950, 1000));
            let mut ex = BTreeMap::new();
            ex.insert(("c".to_string(), 1001, 3000), 2000);
            ex.insert(("c".to_string(), 1201, 2200), 1000);
            ex.insert(("c".to_string(), 9001, 11000), 2000);
            ex.insert(("c".to_string(), 20001, 21000), 1000);
            let p = GraphParams { min_bp: 100, ..GraphParams::default() };
            let plain = graph_from_paf(&paf, &ex, &BTreeMap::new(), &p);
            assert_eq!((plain.n_nodes(), plain.n_edges()), (4, 2), "{:?}", plain.genes);
            let mut bl: BTreeMap<GeneKey, Vec<(u64, u64)>> = BTreeMap::new();
            bl.insert(("c".to_string(), 1001, 3000), vec![(1000, 3000)]);
            bl.insert(("c".to_string(), 1201, 2200), vec![(1200, 2200)]); // over the host's exon
            bl.insert(("c".to_string(), 9001, 11000), vec![(9000, 11000)]);
            bl.insert(("c".to_string(), 20001, 21000), vec![(20000, 21000)]);
            let mut m = loci_from_exon_blocks(&bl);
            assert_eq!(m.n_merged(), 1);
            let rep_only = graph_from_paf_loci(&paf, &ex, &BTreeMap::new(), &p, Some(&m));
            assert_eq!((rep_only.n_nodes(), rep_only.n_edges(), rep_only.same_locus_records), (2, 1, 1), "{:?}", rep_only.genes);
            assert!(!rep_only.genes.contains(&("c".to_string(), 1201, 2200)));
            m.attribute_edges = true;
            let attributed = graph_from_paf_loci(&paf, &ex, &BTreeMap::new(), &p, Some(&m));
            assert_eq!((attributed.n_nodes(), attributed.n_edges(), attributed.same_locus_records), (3, 2, 0), "{:?}", attributed.genes);
            assert!(!attributed.genes.contains(&("c".to_string(), 1201, 2200)));
            // a record between H and N (the locus aligned to itself) is skipped under both policies
            paf.push('\n');
            paf.push_str(&paf_line("c:1201-2200", 1000, 0, 1000, "c:1001-3000", 2000, 200, 1200, 1000, 1000));
            let g2 = graph_from_paf_loci(&paf, &ex, &BTreeMap::new(), &p, Some(&m));
            assert_eq!((g2.n_nodes(), g2.n_edges(), g2.same_locus_records), (3, 2, 1));
        }

        #[test]
        fn folding_within_clusters_keeps_overlapping_records_of_different_families_apart() {
            // Host A1 (c:1001-3000) contains B1 (c:1201-2200) on exon bases; A1' (c:1001-2500) is a second record
            // of A1's own gene. A1, A1' ~ A2 (c:9001-11000); B1 ~ B2 (c:20001-21000). Fold-first loses B1 (its
            // locus's representative is A1); fold-within-clusters keeps {A1, A2} and {B1, B2}: four loci.
            let mut g = HomologyGraph::default();
            let keys: Vec<GeneKey> = [
                ("c", 1001u64, 3000u64), ("c", 1001, 2500), ("c", 9001, 11000), ("c", 1201, 2200), ("c", 20001, 21000),
            ]
            .iter()
            .map(|(c, s, e)| (c.to_string(), *s, *e))
            .collect();
            g.genes = keys.clone();
            let mut bl: BTreeMap<GeneKey, Vec<(u64, u64)>> = BTreeMap::new();
            bl.insert(keys[0].clone(), vec![(1000, 3000)]);
            bl.insert(keys[1].clone(), vec![(1000, 2500)]);
            bl.insert(keys[2].clone(), vec![(9000, 11000)]);
            bl.insert(keys[3].clone(), vec![(1200, 2200)]);
            bl.insert(keys[4].clone(), vec![(20000, 21000)]);
            // fold-first: one locus {A1, A1', B1} with A1 as representative -> B1 gone
            let first = loci_from_exon_blocks(&bl);
            assert_eq!(first.rep_of.get(&keys[3]), Some(&keys[0]));
            // fold within the MCL parts {A1, A1', A2}, {B1, B2}
            let parts = vec![vec![0, 1, 2], vec![3, 4]];
            let (folded, m) = fold_parts_into_loci(&g, &parts, &bl);
            assert_eq!(folded, vec![vec![0, 2], vec![3, 4]]);
            assert_eq!(m.rep_of.get(&keys[1]), Some(&keys[0]));
            assert!(!m.rep_of.contains_key(&keys[3]), "B1 must not be folded into A1's locus");
            assert_eq!((m.n_merged(), m.n_multi), (1, 1));
        }

        #[test]
        fn exonic_both_sides_rejects_a_pair_that_touches_only_the_hosts_exons() {
            // Pseudogene N (c:1201-2200, one 200-bp exon at 1400-1600) lies inside host H (c:1001-3000, exons at
            // 1000-1200 and 2800-3000). N's SPAN carries H's bases; it aligns to H's paralogue P (c:9001-11000,
            // exons 9000-9200, 10800-11000) over N's first 200 bp (H's exon 1000-1200 in span coordinates 0..200),
            // i.e. on P's exons but NOT on N's own exon (span offsets 200..400 untouched).
            let paf = paf_line("c:1201-2200", 1000, 0, 200, "c:9001-11000", 2000, 0, 200, 195, 200);
            let mut ex = BTreeMap::new();
            ex.insert(("c".to_string(), 1201, 2200), 200);
            ex.insert(("c".to_string(), 9001, 11000), 400);
            let mut bl: BTreeMap<GeneKey, Vec<(u64, u64)>> = BTreeMap::new();
            bl.insert(("c".to_string(), 1201, 2200), vec![(1400, 1600)]);
            bl.insert(("c".to_string(), 9001, 11000), vec![(9000, 9200), (10800, 11000)]);
            let one_side = GraphParams { min_bp: 100, min_exonic_bp: 1, min_cov_longer: 0.1, ..GraphParams::default() };
            let g1 = graph_from_paf(&paf, &ex, &bl, &one_side);
            assert_eq!(g1.n_edges(), 1, "the longer side (P) is touched on its exons: admitted today");
            let both = GraphParams { exonic_both_sides: true, ..one_side };
            let g2 = graph_from_paf(&paf, &ex, &bl, &both);
            assert_eq!((g2.n_edges(), g2.rejected_no_exonic), (0, 1), "N's own exon is untouched: no edge");
            // a real paralogue pair maps exon onto exon in ONE record: N's exon (span 200..400) vs P's exon 2
            let paf2 = paf_line("c:1201-2200", 1000, 200, 400, "c:9001-11000", 2000, 1800, 2000, 195, 200);
            let g3 = graph_from_paf(&paf2, &ex, &bl, &both);
            assert_eq!(g3.n_edges(), 1);
            // two records that each touch only ONE side's exons do not add up to exon-to-exon homology
            let mut paf3 = paf.clone();
            paf3.push('\n');
            paf3.push_str(&paf_line("c:1201-2200", 1000, 200, 400, "c:9001-11000", 2000, 400, 600, 195, 200));
            let g4 = graph_from_paf(&paf3, &ex, &bl, &both);
            assert_eq!((g4.n_edges(), g4.rejected_no_exonic), (0, 1));
        }

        /// ⭐ §6x4/§6x5 CONTAINMENT ESCAPE (`min_cov_shorter`).
        ///
        /// The population it exists for (§6x3/r1002): a SHORT locus that aligns along ~100% of its own length
        /// to a LONG partner, at passing identity and `alen`, and is rejected only because `cov_longer`
        /// divides by the long partner's exonic length. Short S (400 bp exonic) aligns end to end onto long
        /// L (4,000 bp exonic): `cov_longer` = 400/4000 = 0.10 < 0.30, `cov_shorter` = 400/400 = 1.00.
        #[test]
        fn min_cov_shorter_admits_a_fully_contained_locus_and_is_a_no_op_when_zero() {
            let paf = paf_line("c:1001-1400", 400, 0, 400, "c:5001-9000", 4000, 0, 400, 390, 400);
            let mut ex = BTreeMap::new();
            ex.insert(("c".to_string(), 1001, 1400), 400);
            ex.insert(("c".to_string(), 5001, 9000), 4000);
            let mut bl: BTreeMap<GeneKey, Vec<(u64, u64)>> = BTreeMap::new();
            bl.insert(("c".to_string(), 1001, 1400), vec![(1000, 1400)]);
            bl.insert(("c".to_string(), 5001, 9000), vec![(5000, 9000)]);
            // OFF (the shipped default) — cov_longer 0.10 < 0.30, so no edge, and the counter stays 0.
            let off = GraphParams { min_exonic_bp: 1, exonic_both_sides: true, ..GraphParams::default() };
            let g_off = graph_from_paf(&paf, &ex, &bl, &off);
            assert_eq!((g_off.n_edges(), g_off.admitted_by_containment), (0, 0));
            // ON at 0.90 — cov_shorter 1.00 clears it, and the admission is counted, never silent.
            let on = GraphParams { min_cov_shorter: 0.90, ..off };
            let g_on = graph_from_paf(&paf, &ex, &bl, &on);
            assert_eq!((g_on.n_edges(), g_on.admitted_by_containment), (1, 1));
            // the weight uses the coverage that actually admitted the pair, not the failing cov_longer
            let w = *g_on.edges.values().next().unwrap();
            assert!((w - 0.975).abs() < 1e-6, "identity 0.975 * cov_shorter 1.0, got {w}");
            // ON but above the pair's containment — still rejected, so the floor is a real threshold.
            let strict = GraphParams { min_cov_shorter: 0.90, min_cov_longer: 0.30, ..off };
            let g_strict = graph_from_paf(
                &paf_line("c:1001-1400", 400, 0, 200, "c:5001-9000", 4000, 0, 200, 195, 200),
                &ex, &bl, &strict,
            );
            assert_eq!((g_strict.n_edges(), g_strict.admitted_by_containment), (0, 0),
                "cov_shorter 200/400 = 0.50 < 0.90 and cov_longer 0.05 < 0.30");
        }

        /// ⭐ §6ks: `min_exonic_bp = 1` is a zero/non-zero gate — it is satisfied by ONE shared exonic base,
        /// however small a sliver of either gene's model that base is. `min_shared_exon_frac` asks for a
        /// FRACTION of the smaller gene's own exonic length instead, which is what separated Soto's agreed
        /// pairs (median 52-89% shared) from this project's extra pairs (median 11-25%, ledger §6kr).
        #[test]
        fn min_shared_exon_frac_rejects_a_small_slice_and_keeps_a_majority_overlap() {
            let a: GeneKey = ("g1".to_string(), 1, 2000);
            let b: GeneKey = ("g2".to_string(), 1, 2000);
            let ex: BTreeMap<GeneKey, u64> = [(a.clone(), 1000u64), (b.clone(), 1000u64)].into_iter().collect();
            // both genes' only exon is their first 1000 bp (local offsets 0..1000)
            let bl: BTreeMap<GeneKey, Vec<(u64, u64)>> =
                [(a.clone(), vec![(1u64, 1001u64)]), (b.clone(), vec![(1u64, 1001u64)])].into_iter().collect();

            // one record overlapping both exons in their first 200 bp only: shared/min(1000,1000) = 0.20
            let low = paf_line("g1:1-2000", 2000, 0, 200, "g2:1-2000", 2000, 0, 200, 200, 200);
            let p_off = GraphParams { min_bp: 100, min_exonic_bp: 1, min_cov_longer: 0.1, ..GraphParams::default() };
            let g_off = graph_from_paf(&low, &ex, &bl, &p_off);
            assert_eq!(g_off.n_edges(), 1, "off by default: the structural 1-bp floor alone admits it");
            assert_eq!(GraphParams::default().min_shared_exon_frac, 0.0);

            let p_strict = GraphParams { min_shared_exon_frac: 0.3, ..p_off };
            let g_strict = graph_from_paf(&low, &ex, &bl, &p_strict);
            assert_eq!(
                (g_strict.n_edges(), g_strict.rejected_low_shared_exon),
                (0, 1),
                "20% of the smaller gene's exon is not enough at a 0.3 floor"
            );
            assert_eq!(g_strict.rejected_no_exonic, 0, "counted under its own field, not the 1-bp one");

            // a majority-overlap record instead (600/1000 = 0.60) must survive the same floor
            let high = paf_line("g1:1-2000", 2000, 0, 600, "g2:1-2000", 2000, 0, 600, 600, 600);
            let g_high = graph_from_paf(&high, &ex, &bl, &p_strict);
            assert_eq!(g_high.n_edges(), 1, "60% shared clears a 0.3 floor");
            assert_eq!(g_high.rejected_low_shared_exon, 0);

            // independent of `exonic_both_sides`: setting the fraction alone (that flag OFF) still gates
            let p_alone = GraphParams { min_shared_exon_frac: 0.3, exonic_both_sides: false, ..p_off };
            assert_eq!(graph_from_paf(&low, &ex, &bl, &p_alone).n_edges(), 0);
        }

        #[test]
        fn sd_blocks_link_hulls_through_one_pair_and_keep_unlinked_hulls_apart() {
            // pair 1: c:1000-5000 <-> c:20000-24000 spans two modules (hulls A=c:1000-2000, B=c:3000-4000 on one
            // side; C=c:20000-21000, D=c:23000-24000 on the other) -> one block {A,B,C,D}.
            // hull E=c:50000-51000 has a pair only to c:70000-71000 (no hull) -> its own block.
            let bed = "c\t1000\t5000\tc\t20000\t24000\nc\t50000\t51000\tc\t70000\t71000\n";
            let sd = SdPairs::from_bed_str(bed);
            let hulls: Vec<(String, u64, u64)> = [(1000, 2000), (3000, 4000), (20000, 21000), (23000, 24000), (50000, 51000)]
                .iter()
                .map(|&(s, e)| ("c".to_string(), s, e))
                .collect();
            let (b, links) = sd_blocks(&hulls, &sd);
            assert_eq!(b, vec![0, 0, 0, 0, 1], "{b:?}");
            // direct links: every hull on one side of pair 1 to every hull on the other side; E links to nothing
            assert_eq!(links, vec![(0, 2), (0, 3), (1, 2), (1, 3)]);
        }

        #[test]
        fn mcl_is_deterministic() {
            let mut g = HomologyGraph::default();
            for i in 0..8 {
                g.genes.push(("c".to_string(), i as u64 + 1, i as u64 + 2));
            }
            for (a, b) in [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (5, 6), (6, 7)] {
                g.edges.insert((a, b), 0.9);
            }
            let a = mcl(&g, 2.8, 100, 1e-5);
            let b = mcl(&g, 2.8, 100, 1e-5);
            assert_eq!(a, b, "same graph and parameters must give the same clustering");
        }

        #[test]
        fn density_separates_a_perfect_clique_from_a_real_family() {
            // §6de, genome-wide and size-controlled: zero-corroboration clusters have median density 1.000,
            // corroborated ones 0.700, at identical median size 4.0.
            let mut g = HomologyGraph::default();
            for i in 0..8 {
                g.genes.push(("c".to_string(), i as u64 + 1, i as u64 + 2));
            }
            for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
                g.edges.insert((a, b), 0.9); // a perfect 4-clique
            }
            for (a, b) in [(4, 5), (5, 6), (6, 7), (4, 7)] {
                g.edges.insert((a, b), 0.9); // a 4-cycle: 4 of 6 possible pairs
            }
            let parts = vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]];
            let cs = build_clusters(&g, &parts, 3, None);
            assert_eq!(cs[0].density, 1.0, "the clique is perfect");
            assert!((cs[1].density - 4.0 / 6.0).abs() < 1e-9, "the family is not");
            assert!(cs[0].corroborated.is_none(), "no BAM => NOT MEASURED, never 0.0");
        }

        #[test]
        fn corroboration_is_none_not_zero_when_unmeasured() {
            // ⚠ `None` (no BAM) and `Some(0.0)` (measured, no support — the repeat-clique signature) are
            // different claims and must never be conflated.
            let mut g = HomologyGraph::default();
            for i in 0..3 {
                g.genes.push(("c".to_string(), i as u64 + 1, i as u64 + 2));
            }
            g.edges.insert((0, 1), 0.9);
            g.edges.insert((1, 2), 0.9);
            g.edges.insert((0, 2), 0.9);
            let parts = vec![vec![0, 1, 2]];
            assert!(build_clusters(&g, &parts, 3, None)[0].corroborated.is_none());
            let never = |_: &GeneKey| false;
            assert_eq!(build_clusters(&g, &parts, 3, Some(&never))[0].corroborated, Some(0.0));
        }

        /// ⭐ RED-BEFORE: the MCL0 mechanism. A 794 bp bridge that is 1.1% exonic saturates `cov_longer`
        /// under the span numerator (794 non-exonic bases / a 679 bp exonic denominator) and is REJECTED
        /// once the numerator is charged in exonic bases too.
        #[test]
        fn exonic_overlap_rejects_a_non_exonic_bridge_the_span_numerator_admits() {
            let q: GeneKey = ("NC_1".to_string(), 1, 5000);
            let t: GeneKey = ("NC_2".to_string(), 1, 5000);
            // The alignment sits at local [2400, 3194) — 794 bp, entirely INTRONIC on both genes.
            let paf = "NC_1:1-5000\t5000\t2400\t3194\t+\tNC_2:1-5000\t5000\t2400\t3194\t790\t794\n";
            let exonic: BTreeMap<GeneKey, u64> =
                [(q.clone(), 679u64), (t.clone(), 679u64)].into_iter().collect();
            // Exons live at absolute 0-based [0, 679) — disjoint from the alignment.
            let blocks: BTreeMap<GeneKey, Vec<(u64, u64)>> =
                [(q.clone(), vec![(0u64, 679u64)]), (t.clone(), vec![(0u64, 679u64)])].into_iter().collect();

            let span = graph_from_paf(paf, &exonic, &blocks, &GraphParams::default());
            assert_eq!(span.n_edges(), 1, "span numerator admits the intronic bridge (794/679 -> capped 1.0)");

            let p = GraphParams { exonic_overlap: true, ..GraphParams::default() };
            let exon = graph_from_paf(paf, &exonic, &blocks, &p);
            assert_eq!(exon.n_edges(), 0, "exonic numerator sees 0 shared exonic bases and rejects it");
            assert_eq!(exon.exonic_overlap_joined, 1, "the join MUST fire - a 0 join is the §6dd bug");
            assert_eq!(exon.exonic_overlap_missing, 0);
        }

        /// A genuinely exonic alignment must still be admitted — the fix is a restriction, not a wall.
        #[test]
        fn exonic_overlap_keeps_a_real_exonic_alignment() {
            let q: GeneKey = ("NC_1".to_string(), 1, 5000);
            let t: GeneKey = ("NC_2".to_string(), 1, 5000);
            let paf = "NC_1:1-5000\t5000\t0\t679\t+\tNC_2:1-5000\t5000\t0\t679\t670\t679\n";
            let exonic: BTreeMap<GeneKey, u64> =
                [(q.clone(), 679u64), (t.clone(), 679u64)].into_iter().collect();
            let blocks: BTreeMap<GeneKey, Vec<(u64, u64)>> =
                [(q.clone(), vec![(0u64, 679u64)]), (t.clone(), vec![(0u64, 679u64)])].into_iter().collect();
            let p = GraphParams { exonic_overlap: true, ..GraphParams::default() };
            let g = graph_from_paf(paf, &exonic, &blocks, &p);
            assert_eq!(g.n_edges(), 1, "fully exonic alignment survives");
            assert_eq!(g.exonic_overlap_joined, 1);
        }

        /// ⚠ The §6dd coordinate trap, pinned: blocks are ABSOLUTE, the gene start is GFF 1-based.
        #[test]
        fn exonic_bases_in_converts_absolute_blocks_to_the_local_frame() {
            // Gene at GFF 1-based 1001..2000; exon at absolute 0-based [1000, 1100) == local [0, 100).
            assert_eq!(exonic_bases_in(&[(1000, 1100)], 1001, 0, 100), 100);
            assert_eq!(exonic_bases_in(&[(1000, 1100)], 1001, 50, 150), 50);
            assert_eq!(exonic_bases_in(&[(1000, 1100)], 1001, 200, 300), 0);
        }

        /// The nesting artefact: 108 intra-cluster edges joined fully nested intervals at identity 1.000.
        #[test]
        fn reject_overlapping_drops_a_nested_annotation_pair() {
            let outer: GeneKey = ("NC_1".to_string(), 1000, 9000);
            let inner: GeneKey = ("NC_1".to_string(), 2000, 3000);
            let paf = "NC_1:1000-9000\t8001\t1000\t2001\t+\tNC_1:2000-3000\t1001\t0\t1001\t1001\t1001\n";
            let exonic: BTreeMap<GeneKey, u64> =
                [(outer.clone(), 1001u64), (inner.clone(), 1001u64)].into_iter().collect();
            let b = BTreeMap::new();
            assert_eq!(graph_from_paf(paf, &exonic, &b, &GraphParams::default()).n_edges(), 1);
            let p = GraphParams { reject_overlapping: true, ..GraphParams::default() };
            let g = graph_from_paf(paf, &exonic, &b, &p);
            assert_eq!(g.n_edges(), 0, "a gene nested inside another is not its paralog");
            assert_eq!(g.rejected_overlapping, 1);
        }

        /// Both new clauses are OFF by default, so every prior measurement stays byte-identical.
        #[test]
        fn the_new_clauses_are_off_by_default() {
            let d = GraphParams::default();
            assert!(!d.exonic_overlap, "flipping this is a THESIS EDIT, not a code edit");
            assert!(!d.reject_overlapping);
            assert_eq!(d.min_shared_exon_frac, 0.0, "flipping this is a THESIS EDIT, not a code edit");
        }


        /// ⭐ The shipped remedy: an ADDITIVE exonic floor keeps a real, intron-spanning paralogy edge and
        /// drops one that rests on no exonic sequence. Measured on the pilot: NPIP 43/43 intact, the
        /// 33-gene repeat clique reduced to 1, MCL0 split at its adjudicated cut vertex.
        #[test]
        fn min_exonic_bp_drops_a_non_exonic_edge_and_keeps_an_intron_spanning_one() {
            let q: GeneKey = ("NC_1".to_string(), 1, 20000);
            let t: GeneKey = ("NC_2".to_string(), 1, 20000);
            let exonic: BTreeMap<GeneKey, u64> =
                [(q.clone(), 4000u64), (t.clone(), 4000u64)].into_iter().collect();
            // Exons at absolute [0,2000) and [15000,17000); introns everywhere else.
            let blocks: BTreeMap<GeneKey, Vec<(u64, u64)>> = [
                (q.clone(), vec![(0u64, 2000u64), (15000u64, 17000u64)]),
                (t.clone(), vec![(0u64, 2000u64), (15000u64, 17000u64)]),
            ]
            .into_iter()
            .collect();
            let p = GraphParams { min_exonic_bp: 300, ..GraphParams::default() };

            // A whole-locus duplication: 16 kb spanning introns AND both exons. Must SURVIVE.
            let real = "NC_1:1-20000\t20000\t500\t16500\t+\tNC_2:1-20000\t20000\t500\t16500\t15400\t16000\n";
            let g = graph_from_paf(real, &exonic, &blocks, &p);
            assert_eq!(g.n_edges(), 1, "an intron-spanning real paralogy edge must survive");
            assert_eq!(g.rejected_no_exonic, 0);

            // A purely intronic repeat bridge of the same length. Must DIE.
            let rep = "NC_1:1-20000\t20000\t3000\t13000\t+\tNC_2:1-20000\t20000\t3000\t13000\t8500\t10000\n";
            let g2 = graph_from_paf(rep, &exonic, &blocks, &p);
            assert_eq!(g2.n_edges(), 0, "an edge resting on zero exonic bases must be rejected");
            assert_eq!(g2.rejected_no_exonic, 1, "and the rejection must be COUNTED, never silent");

            // Off by default: the same repeat bridge is admitted on the default path.
            assert_eq!(graph_from_paf(rep, &exonic, &blocks, &GraphParams::default()).n_edges(), 1);
            assert_eq!(GraphParams::default().min_exonic_bp, 0);
        }


        /// ⭐ `frac_in` sees a hairball slice that `density` cannot: a 4-clique cut out of a larger blob has
        /// density 1.000 and cohesion 0.5. Calibrated on the pilot at real 0.923-0.970 vs artefact 0.283.
        #[test]
        fn frac_in_catches_a_slice_that_density_calls_perfect() {
            let mut g = HomologyGraph::default();
            for i in 0..8u64 {
                g.genes.push(("NC_1".to_string(), i * 1000 + 1, i * 1000 + 500));
            }
            // 0..4 is a perfect clique; each of its members also leaks one edge to the outside block.
            for a in 0..4 {
                for b in (a + 1)..4 {
                    g.edges.insert((a, b), 1.0);
                }
                g.edges.insert((a, 4 + a), 1.0);
            }
            let c = build_clusters(&g, &[vec![0, 1, 2, 3]], 3, None);
            assert_eq!(c.len(), 1);
            assert!((c[0].density - 1.0).abs() < 1e-9, "density is blind: a perfect internal clique");
            assert!(
                (c[0].frac_in - 0.6).abs() < 1e-9,
                "cohesion sees the leak: 6 internal / (6 internal + 4 leaving) = 0.6, got {}",
                c[0].frac_in
            );
        }

    }
}
