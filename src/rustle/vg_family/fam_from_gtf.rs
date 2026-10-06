//! **STATUS:** SHIPPED-DEFAULT
//!
//! The `--from-gtf` family stage as a library: loci of an assembled GTF, their all-vs-all
//! (`minimap2 -x asm20`), and the copy table (`write_locus_rep_copies`) consumed by
//! `copy_assign --families`. Extracted verbatim from `src/bin/mcl_families.rs` (2026-10-04);
//! `mcl_families` itself re-imports every moved function. ⚠ Every function here is
//! byte-identity-critical: families-stage products are cmp-checked across builds.

use anyhow::{Context, Result};
use crate::genome::GenomeIndex;
use crate::vg_family::annotation_families::{Cluster, CoreRecord, CoreStatus, GeneKey};
use crate::vg_family::denovo_assemble::longest_orf;
use std::collections::BTreeMap;
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
        seq = crate::vg_family::seq_utils::revcomp_keep_case(&seq);
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
    g: &crate::vg_family::annotation_families::HomologyGraph,
    core_records: &[Vec<crate::vg_family::annotation_families::CoreRecord>],
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
    let contigs: std::collections::HashSet<String> =
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
                Some(crate::vg_family::annotation_families::CoreStatus::Dropped) => "dropped",
                Some(crate::vg_family::annotation_families::CoreStatus::KeptTrimmed) => "kept_trimmed",
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
pub fn relation_loci(loci: &[GtfLocus]) -> Vec<crate::vg_family::family_relations::LocusIn> {
    loci.iter()
        .map(|l| crate::vg_family::family_relations::LocusIn {
            gene_id: l.gene_id.clone(),
            key: (l.chrom.clone(), l.start as i64, l.end as i64),
            rep: l.rep.clone(),
            rep_reads: l.rep_reads as i64,
        })
        .collect()
}

/// Write `<out>.loci.fa`: one record `>CONTIG:START-END` + the genome's forward strand over that span, per locus span.
/// With `hash`, also the [`ContentHash`](crate::vg_family::run_cache::ContentHash) of every byte written (the PAF
/// cache key), taken as the bytes go out so the file is never read back.
pub fn write_loci_fa(
    genome: &GenomeIndex,
    spans: &[(String, u64, u64)],
    fa_path: &str,
    fasta: &str,
    hash: bool,
) -> Result<Option<crate::vg_family::run_cache::ContentHash>> {
    use crate::vg_family::run_cache::HashingWriter;
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
pub fn families_paf_key(cmd: &str, minimap2_version: &str, loci_fa: &crate::vg_family::run_cache::ContentHash) -> String {
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
    use crate::vg_family::run_cache as rc;
    let root = rc::cache_root();
    let contigs: HashSet<String> = spans.iter().map(|x| x.0.clone()).collect();
    let genome = GenomeIndex::from_fasta_contigs(fasta, &contigs)?;
    let fa_hash = write_loci_fa(&genome, &spans, &fa_path, fasta, root.is_some())?;
    eprintln!("[mcl_families] --from-gtf: {} loci -> all-vs-all", spans.len());
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let mm_args: Vec<String> = ["-x", "asm20", "-c", "-X", "-N", "50", "-p", "0.1", "--secondary=yes", "-t"]
        .iter()
        .map(|s| s.to_string())
        .chain(std::iter::once(threads.to_string()))
        .collect();
    // PAF cache (`RUSTLE_CACHE_DIR`, see `crate::vg_family::run_cache`): keyed by EVERY byte of the loci FASTA
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
    use crate::vg_family::run_cache as rc;

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
