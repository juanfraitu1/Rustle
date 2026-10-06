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
pub fn relation_loci(loci: &[GtfLocus]) -> Vec<crate::vg_family::fam_from_gtf::family_relations::LocusIn> {
    loci.iter()
        .map(|l| crate::vg_family::fam_from_gtf::family_relations::LocusIn {
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
use crate::vg_family::fam_from_gtf::family_container::{gtf_attr as attr, key_str, parse_key, py_int, Key};

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
        use crate::vg_family::bridge_regroup::{run_list, UnitsList};
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
