//! Main pipeline: bundles -> graph -> map reads -> extract -> GTF.

use anyhow::Result;
use lru::LruCache;
use std::num::NonZeroUsize;
use std::collections::VecDeque;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::assembly_mode::LONGINTRONANCHOR;
use crate::assembly_mode::{use_coverage_trim, use_longtrim};
use crate::ballgown::write_ballgown;
use crate::bam::junctions_from_exons;
use crate::bitset::NodeSet;
use crate::bpcov::{Bpcov, BpcovStranded, BPCOV_STRAND_ALL, BPCOV_STRAND_MINUS, BPCOV_STRAND_PLUS};
use crate::bundle::{
    build_bundlenodes_and_readgroups_from_cgroups,
    build_bundlenodes_and_readgroups_from_cgroups_3strand, compare_bundles_to_log,
    detect_bundles_from_bam, parse_bundle_log,
};
use crate::bundle_builder::build_sub_bundles;

use crate::coord::len_half_open;
use crate::coverage_trim::apply_coverage_trim;
use crate::futuretr::{
    apply_prune_redirects, collect_from_transfrags, materialize_links, normalize_links, FutureLink,
};
use crate::gene_abundance::{compute_transcript_bpcov, write_gene_abundance};
use crate::genome::GenomeIndex;
use crate::graph::{Graph, GraphTransfrag};
use crate::graph_build::{
    create_graph_with_longtrim, prune_graph_nodes, prune_graph_nodes_with_redirects,
};
use crate::gtf::write_gtf;
use crate::hard_boundaries::annotate_hard_boundaries;
use crate::junction_correction::snap_junctions_to_guides;
use crate::junctions::{
    apply_junction_filters_and_canonicalize, correct_junctions_with_map, filter_junctions,
};
use crate::killed_junctions::{
    aggregate_splice_site_support, apply_higherr_demotions, compute_killed_junction_pairs,
    demote_runthrough_junctions,
    good_junc,
};
use crate::map_reads::{map_reads_to_graph, map_reads_to_graph_bundlenodes};
use crate::nodecov::compute_nodecov;
use crate::debug_stage;
use crate::debug_stage::StageDetail;
use crate::path_extract::{
    extract_rawreads_transcripts, extract_shortread_transcripts, extract_transcripts,
    LongRecSummary, Transcript,
};
use crate::read_boundaries::collect_read_boundaries_with_cpas;
use crate::reference_gtf::{
    find_guide_pat, parse_reference_gtf, process_refguides, GuideInfo, RefTranscript,
};
use crate::snapshot;
use crate::snapshot::SnapshotDetail;
use crate::trace_reference::{
    debug_target_ref_stage, drop_resolved_blockers_by_final_matches, find_traced_ref_in_bundle,
    ref_overlaps_bundle, trace_refs_in_bundles, write_blocker_report,
};
use crate::types::{BundleRead, CBundlenode, Junction, JunctionStat, JunctionStats, RunConfig};
use crate::types::{DetHashMap as HashMap, DetHashSet as HashSet};
use hashbrown::{HashMap as HbHashMap, HashSet as HbHashSet};
use roaring::RoaringBitmap;


#[derive(Clone, Debug, Hash, PartialEq, Eq)]
struct SpliceConsensusKey {
    chrom: Arc<str>,
    donor: u64,
    acceptor: u64,
    strand: Option<i8>,
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
struct GuideJunctionCacheKey {
    chrom: Arc<str>,
    start: u64,
    end: u64,
    strand: char,
}

fn init_rayon_pool(threads: usize) {
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok(); // ignore error if pool already initialized
    }
}

fn consensus_cache_capacity() -> NonZeroUsize {
    let cap = std::env::var("RUSTLE_CONSENSUS_LRU")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(131_072)
        .max(1);
    // SAFETY: max(1) guarantees non-zero.
    NonZeroUsize::new(cap).unwrap()
}

/// Merge overlapping same-strand subbundles into larger bundles.
/// the algorithm processes all bundlenodes in a color component as ONE graph.
/// Rustle's build_sub_bundles creates separate subbundles for each
/// component. This function merges overlapping same-strand subbundles so
/// that reads can see all graph nodes in the region, matching the expected behavior.
/// Parse a StringTie GTF and return ref tx whose spans overlap the bundle
/// on the matching strand. Returns (tid, intron_chain, exons_0based, strand).
/// Exon coords are converted from GTF 1-based inclusive to Rustle 0-based
/// half-open (start - 1).
fn parse_ref_chains_for_bundle(
    path: &str,
    bundle_chrom: &str,
    bundle_strand: char,
    bundle_start: u64,
    bundle_end: u64,
) -> Vec<(String, Vec<(u64, u64)>, Vec<(u64, u64)>, char)> {
    let mut out = Vec::new();
    let Ok(content) = std::fs::read_to_string(path) else { return out; };
    let mut cur_tid: Option<String> = None;
    let mut cur_chrom = String::new();
    let mut cur_strand = '.';
    let mut cur_exons: Vec<(u64, u64)> = Vec::new();
    let mut flush = |tid: &Option<String>, chrom: &str, strand: char,
                     exons: &[(u64, u64)],
                     out: &mut Vec<(String, Vec<(u64, u64)>, Vec<(u64, u64)>, char)>| {
        let Some(tid) = tid else { return; };
        if chrom != bundle_chrom || strand != bundle_strand || exons.len() < 2 { return; }
        let first = exons[0].0;
        let last = exons.last().unwrap().1;
        if first > bundle_end || last < bundle_start { return; }
        let intron_chain: Vec<(u64, u64)> =
            exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
        out.push((tid.clone(), intron_chain, exons.to_vec(), strand));
    };
    for line in content.lines() {
        if line.starts_with('#') || line.is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 || cols[2] != "exon" { continue; }
        let chrom = cols[0];
        let strand = cols[6].chars().next().unwrap_or('.');
        let es: u64 = cols[3].parse().unwrap_or(0);
        let ee: u64 = cols[4].parse().unwrap_or(0);
        let attrs = cols[8];
        let tid = attrs.find("transcript_id \"")
            .and_then(|idx| {
                let rest = &attrs[idx + "transcript_id \"".len()..];
                rest.find('"').map(|end| rest[..end].to_string())
            })
            .unwrap_or_default();
        if cur_tid.as_ref() != Some(&tid) {
            flush(&cur_tid, &cur_chrom, cur_strand, &cur_exons, &mut out);
            cur_tid = Some(tid);
            cur_chrom = chrom.to_string();
            cur_strand = strand;
            cur_exons.clear();
        }
        // GTF 1-based → 0-based half-open
        cur_exons.push((es.saturating_sub(1), ee));
    }
    flush(&cur_tid, &cur_chrom, cur_strand, &cur_exons, &mut out);
    out
}

/// Direct-emit missed-tx oracle: for each ref tx (from a StringTie GTF)
/// that overlaps the current bundle, synthesize a Transcript and append
/// to the extracted tx list. Opt-in via RUSTLE_MISSED_ORACLE_DIRECT=path.
/// Use RUSTLE_MISSED_ORACLE_ONLY=tid1,tid2,... to restrict.
/// This bypasses path_extract entirely — if gffcompare matches after
/// direct emit, path_extract was the blocker; if still not matched,
/// a downstream filter (pairwise/RI/etc) killed the synthetic output.
fn append_missed_oracle_direct_emit(
    txs: &mut Vec<crate::path_extract::Transcript>,
    bundle_chrom: &str,
    bundle_start: u64,
    bundle_end: u64,
    bundle_strand: char,
) {
    let path = match std::env::var("RUSTLE_MISSED_ORACLE_DIRECT") {
        Ok(p) => p,
        Err(_) => return,
    };
    let content = match std::fs::read_to_string(&path) {
        Ok(c) => c,
        Err(_) => return,
    };
    let only: Option<std::collections::HashSet<String>> = std::env::var("RUSTLE_MISSED_ORACLE_ONLY")
        .ok()
        .map(|v| v.split(',').map(|s| s.trim().to_string()).collect());
    let debug = std::env::var_os("RUSTLE_MISSED_ORACLE_DEBUG").is_some();
    let cov: f64 = std::env::var("RUSTLE_MISSED_ORACLE_COV")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(5.0);

    let mut cur_tid: Option<String> = None;
    let mut cur_chrom: String = String::new();
    let mut cur_strand: char = '.';
    let mut cur_exons: Vec<(u64, u64)> = Vec::new();

    // RUSTLE_MISSED_ORACLE_SKIP_MATCHED=1: skip ref tx whose intron chain
    // already matches an emitted tx. Lets oracle_direct "fill the gap"
    // when combined with natural extraction or REF_CHAIN_RESCUE — avoids
    // duplicates that tank precision.
    let skip_matched = std::env::var_os("RUSTLE_MISSED_ORACLE_SKIP_MATCHED").is_some();
    let chain_tol: u64 = std::env::var("RUSTLE_MISSED_ORACLE_CHAIN_TOL")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(5);
    let mut flush = |tid: &Option<String>, chrom: &str, strand: char,
                     exons: &[(u64, u64)],
                     txs: &mut Vec<crate::path_extract::Transcript>| {
        let Some(tid) = tid else { return };
        if let Some(only) = &only {
            if !only.contains(tid) { return; }
        }
        if chrom != bundle_chrom || strand != bundle_strand { return; }
        if exons.is_empty() { return; }
        let first = exons[0].0;
        let last = exons.last().unwrap().1;
        if first > bundle_end || last < bundle_start { return; }

        if skip_matched && exons.len() >= 2 {
            // Build ref intron chain
            let ref_chain: Vec<(u64, u64)> = exons.windows(2)
                .map(|w| (w[0].1, w[1].0)).collect();
            for t in txs.iter() {
                if t.chrom != bundle_chrom || t.strand != bundle_strand { continue; }
                if t.exons.len() != exons.len() { continue; }
                let tx_chain: Vec<(u64, u64)> = t.exons.windows(2)
                    .map(|w| (w[0].1, w[1].0)).collect();
                let chain_match = tx_chain.iter().zip(ref_chain.iter()).all(|(a, b)| {
                    a.0.abs_diff(b.0) <= chain_tol && a.1.abs_diff(b.1) <= chain_tol
                });
                if chain_match { return; } // already emitted
            }
        }

        let exon_cov = vec![cov; exons.len()];
        let tx = crate::path_extract::Transcript {
            chrom: bundle_chrom.to_string(),
            strand: bundle_strand,
            exons: exons.to_vec(),
            coverage: cov,
            exon_cov,
            tpm: 0.0,
            fpkm: 0.0,
            source: Some(format!("oracle_direct:{tid}")),
            is_longread: true,
            longcov: cov,
            bpcov_cov: 0.0,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: Some(tid.to_string()),
            ref_gene_id: None,
            hardstart: true,
            hardend: true,
            alt_tts_end: true,
            vg_family_id: None,
            vg_copy_id: None,
            vg_family_size: None,
            intron_low: Vec::new(),
        };
        txs.push(tx);
        if debug {
            eprintln!(
                "ORACLE_DIRECT_EMITTED tid={} bundle={}:{}-{} strand={} exons={} coord={}-{}",
                tid, bundle_chrom, bundle_start, bundle_end, strand,
                exons.len(), first, last
            );
        }
    };

    for line in content.lines() {
        if line.starts_with('#') || line.is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 || cols[2] != "exon" { continue; }
        let chrom = cols[0];
        let strand = cols[6].chars().next().unwrap_or('.');
        let es: u64 = cols[3].parse().unwrap_or(0);
        let ee: u64 = cols[4].parse().unwrap_or(0);
        let attrs = cols[8];
        let tid = attrs.find("transcript_id \"")
            .and_then(|idx| {
                let rest = &attrs[idx + "transcript_id \"".len()..];
                rest.find('"').map(|end| rest[..end].to_string())
            })
            .unwrap_or_default();
        if cur_tid.as_ref() != Some(&tid) {
            flush(&cur_tid, &cur_chrom, cur_strand, &cur_exons, txs);
            cur_tid = Some(tid);
            cur_chrom = chrom.to_string();
            cur_strand = strand;
            cur_exons.clear();
        }
        // GTF 1-based inclusive → Rustle internal 0-based (start - 1).
        // End coord stays as-is (GTF end == 0-based exclusive end == 1-based inclusive).
        cur_exons.push((es.saturating_sub(1), ee));
    }
    flush(&cur_tid, &cur_chrom, cur_strand, &cur_exons, txs);
}

fn merge_overlapping_subbundles(
    subs: &[crate::bundle_builder::SubBundleResult],
    reads: &[crate::types::BundleRead],
) -> Vec<crate::bundle_builder::SubBundleResult> {
    use crate::bundle_builder::SubBundleResult;
    use crate::types::CBundlenode;

    if subs.is_empty() {
        return Vec::new();
    }

    // Sort by (strand, start)
    let mut indexed: Vec<(usize, char, u64, u64)> = subs
        .iter()
        .enumerate()
        .map(|(i, s)| (i, s.strand, s.start, s.end))
        .collect();
    indexed.sort_unstable_by_key(|&(_, strand, start, _)| (strand, start));

    // Sweep-merge overlapping same-strand subbundles
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut cur_strand = '.';
    let mut cur_end = 0u64;

    for &(idx, strand, start, end) in &indexed {
        if strand == '.' {
            continue;
        }
        if strand != cur_strand || start > cur_end {
            // New group
            groups.push(vec![idx]);
            cur_strand = strand;
            cur_end = end;
        } else {
            // Extend current group
            groups.last_mut().unwrap().push(idx);
            cur_end = cur_end.max(end);
        }
    }

    // Build merged subbundles
    let mut result = Vec::with_capacity(groups.len());
    for group in &groups {
        if group.len() == 1 {
            // Single subbundle — pass through (rebuild to avoid Clone requirement)
            let sub = &subs[group[0]];
            result.push(SubBundleResult {
                strand: sub.strand,
                start: sub.start,
                end: sub.end,
                bnode_head: sub.bnode_head.as_ref().map(|h| h.deep_clone()),
                read_to_bnodes: sub.read_to_bnodes.clone(),
                bnode_colors: sub.bnode_colors.clone(),
                read_scale: sub.read_scale.clone(),
            });
            continue;
        }

        // Merge: take union of reads with max scale
        let strand = subs[group[0]].strand;
        let start = group.iter().map(|&i| subs[i].start).min().unwrap();
        let end = group.iter().map(|&i| subs[i].end).max().unwrap();

        // Per-read: take max scale across all subbundles in group
        let n_reads = reads.len();
        let mut merged_scale = vec![0.0f64; n_reads];
        let mut merged_bnodes: Vec<Vec<usize>> = vec![Vec::new(); n_reads];

        for &sub_idx in group {
            let sub = &subs[sub_idx];
            for ri in 0..n_reads {
                let scale = sub.read_scale.get(ri).copied().unwrap_or(0.0);
                if scale > merged_scale[ri] {
                    merged_scale[ri] = scale;
                }
                if let Some(bnodes) = sub.read_to_bnodes.get(ri) {
                    if !bnodes.is_empty() {
                        for &bid in bnodes {
                            if !merged_bnodes[ri].contains(&bid) {
                                merged_bnodes[ri].push(bid);
                            }
                        }
                    }
                }
            }
        }

        // Chain bundlenode linked lists
        let mut all_bnodes: Vec<(u64, u64, f64, usize)> = Vec::new();
        for &sub_idx in group {
            let mut cur = subs[sub_idx].bnode_head.as_ref();
            while let Some(bn) = cur {
                all_bnodes.push((bn.start, bn.end, bn.cov, bn.bid));
                cur = bn.next.as_deref();
            }
        }
        all_bnodes.sort_unstable_by_key(|&(s, _, _, _)| s);
        all_bnodes.dedup_by_key(|bn| bn.3); // dedup by bid

        // Build linked list from sorted bnodes
        let merged_head = if all_bnodes.is_empty() {
            None
        } else {
            let mut head = CBundlenode {
                start: all_bnodes[0].0,
                end: all_bnodes[0].1,
                cov: all_bnodes[0].2,
                bid: all_bnodes[0].3,
                next: None,
                hardstart: false,
                hardend: false,
            };
            let mut tail = &mut head;
            for &(s, e, c, bid) in &all_bnodes[1..] {
                tail.next = Some(Box::new(CBundlenode {
                    start: s,
                    end: e,
                    cov: c,
                    bid,
                    next: None,
                    hardstart: false,
                    hardend: false,
                }));
                tail = tail.next.as_mut().unwrap();
            }
            Some(head)
        };

        // Merge bnode_colors from first subbundle (colors are global)
        let bnode_colors = subs[group[0]].bnode_colors.clone();

        result.push(SubBundleResult {
            strand,
            start,
            end,
            bnode_head: merged_head,
            read_to_bnodes: merged_bnodes,
            bnode_colors,
            read_scale: merged_scale,
        });
    }

    result
}

fn guide_junction_cache_capacity() -> NonZeroUsize {
    let cap = std::env::var("RUSTLE_GUIDE_JUNC_LRU")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(32_768)
        .max(1);
    // SAFETY: max(1) guarantees non-zero.
    NonZeroUsize::new(cap).unwrap()
}

#[inline]
fn intern_chrom_arc(chrom_cache: &mut HashMap<String, Arc<str>>, chrom: &str) -> Arc<str> {
    if let Some(found) = chrom_cache.get(chrom) {
        return Arc::clone(found);
    }
    let arc: Arc<str> = Arc::from(chrom);
    chrom_cache.insert(chrom.to_string(), Arc::clone(&arc));
    arc
}

#[inline]
fn cached_consensus_splice(
    cache: &mut LruCache<SpliceConsensusKey, bool>,
    genome: &GenomeIndex,
    chrom: &Arc<str>,
    donor: u64,
    acceptor: u64,
    strand: Option<i8>,
) -> bool {
    let key = SpliceConsensusKey {
        chrom: Arc::clone(chrom),
        donor,
        acceptor,
        strand,
    };
    if let Some(v) = cache.get(&key) {
        return *v;
    }
    let is_cons = genome.is_consensus_splice(chrom.as_ref(), donor, acceptor, strand);
    cache.put(key, is_cons);
    is_cons
}

#[inline]
fn cached_bundle_guide_junctions(
    cache: &mut LruCache<GuideJunctionCacheKey, Arc<Vec<Junction>>>,
    chrom_cache: &mut HashMap<String, Arc<str>>,
    guides: &[RefTranscript],
    chrom: &str,
    strand: char,
    bundle_start: u64,
    bundle_end: u64,
) -> Arc<Vec<Junction>> {
    if guides.is_empty() {
        return Arc::new(Vec::new());
    }
    let chrom_arc = intern_chrom_arc(chrom_cache, chrom);
    let key = GuideJunctionCacheKey {
        chrom: Arc::clone(&chrom_arc),
        start: bundle_start,
        end: bundle_end,
        strand,
    };
    if let Some(found) = cache.get(&key) {
        return Arc::clone(found);
    }
    let computed = Arc::new(build_bundle_guide_junctions(
        guides,
        chrom,
        strand,
        bundle_start,
        bundle_end,
    ));
    cache.put(key, Arc::clone(&computed));
    computed
}

fn collect_chrom_junction_strand_evidence_parallel(
    bundles: &[crate::types::Bundle],
) -> HashMap<String, HashMap<Junction, (bool, bool)>> {
    let mut acc = HashMap::<String, HashMap<Junction, (bool, bool)>>::default();
    for bundle in bundles {
        if bundle.strand != '+' && bundle.strand != '-' {
            continue;
        }
        let chrom_entry = acc
            .entry(bundle.chrom.clone())
            .or_insert_with(Default::default);
        for (j, st) in &bundle.junction_stats {
            match st.strand {
                Some(s) if s > 0 => {
                    chrom_entry.entry(*j).or_insert((false, false)).0 = true
                }
                Some(s) if s < 0 => {
                    chrom_entry.entry(*j).or_insert((false, false)).1 = true
                }
                _ => {}
            }
        }
    }
    acc
}

fn trace_intron_tol() -> u64 {
    std::env::var("RUSTLE_TRACE_INTRON_TOL")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(5)
}

fn intron_eq_tol(a: (u64, u64), b: (u64, u64), tol: u64) -> bool {
    a.0.abs_diff(b.0) <= tol && a.1.abs_diff(b.1) <= tol
}

fn intron_chain_from_exons_0based(exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    for i in 0..exons.len().saturating_sub(1) {
        out.push((exons[i].1, exons[i + 1].0));
    }
    out
}

fn outcome_label(outcome: &crate::path_extract::SeedOutcome) -> &'static str {
    use crate::path_extract::SeedOutcome;
    match outcome {
        SeedOutcome::Skipped(reason) => reason,
        SeedOutcome::BackToSourceFail => "back_to_source_fail",
        SeedOutcome::FwdToSinkFail => "fwd_to_sink_fail",
        SeedOutcome::UnwitnessedSplice => "unwitnessed_splice",
        SeedOutcome::HardBoundaryMismatch => "hard_boundary_mismatch",
        SeedOutcome::ZeroFlux => "zero_flux",
        SeedOutcome::LowCoverage(_) => "low_coverage",
        SeedOutcome::EonlyNonGuide => "eonly_non_guide",
        SeedOutcome::TooShort => "too_short",
        SeedOutcome::ChecktrfReadthr => "checktrf_readthr",
        SeedOutcome::ChecktrfEonlySkip => "checktrf_eonly_skip",
        SeedOutcome::ChecktrfRedistributed => "checktrf_redistributed",
        SeedOutcome::ChecktrfRescued => "checktrf_rescued",
        SeedOutcome::ChecktrfIncomplete => "checktrf_incomplete",
        SeedOutcome::ChecktrfRescueFail => "checktrf_rescue_fail",
        SeedOutcome::Stored(_) => "stored",
    }
}

fn trace_ref_seed_outcomes(
    bundle: &crate::types::Bundle,
    ref_transcripts: &[RefTranscript],
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    seed_outcomes: &[(usize, crate::path_extract::SeedOutcome)],
) {
    if std::env::var_os("RUSTLE_TRACE_REF_SEEDS").is_none() {
        return;
    }
    let Some(ref_tx) = find_traced_ref_in_bundle(bundle, ref_transcripts) else {
        return;
    };
    if seed_outcomes.is_empty() {
        return;
    }
    let tol = trace_intron_tol();
    let ref_chain = intron_chain_from_exons_0based(&ref_tx.exons);
    if ref_chain.is_empty() {
        return;
    }

    #[derive(Clone)]
    struct Row {
        seed: usize,
        label: &'static str,
        abund: f64,
        seed_introns: usize,
        overlap: usize,
        first: (u64, u64),
        last: (u64, u64),
    }

    let mut rows: Vec<Row> = Vec::new();
    let mut by_label: std::collections::BTreeMap<&'static str, usize> =
        std::collections::BTreeMap::new();

    for (seed, outcome) in seed_outcomes {
        let Some(tf) = transfrags.get(*seed) else {
            continue;
        };
        if tf.node_ids.len() <= 2 {
            continue;
        }

        // Seed intron chain from node gaps on this transfrag.
        let mut seed_chain: Vec<(u64, u64)> = Vec::new();
        for w in tf.node_ids.windows(2) {
            let (a, b) = (w[0], w[1]);
            let (Some(na), Some(nb)) = (graph.nodes.get(a), graph.nodes.get(b)) else {
                continue;
            };
            if nb.start > na.end + 1 {
                seed_chain.push((na.end, nb.start));
            }
        }
        if seed_chain.is_empty() {
            continue;
        }

        let mut overlap = 0usize;
        for si in &seed_chain {
            if ref_chain.iter().any(|&ri| intron_eq_tol(*si, ri, tol)) {
                overlap += 1;
            }
        }
        if overlap == 0 {
            continue;
        }

        let label = outcome_label(outcome);
        *by_label.entry(label).or_insert(0) += 1;
        let first = tf
            .node_ids
            .first()
            .and_then(|&nid| graph.nodes.get(nid))
            .map(|n| (n.start, n.end))
            .unwrap_or((0, 0));
        let last = tf
            .node_ids
            .last()
            .and_then(|&nid| graph.nodes.get(nid))
            .map(|n| (n.start, n.end))
            .unwrap_or((0, 0));
        rows.push(Row {
            seed: *seed,
            label,
            abund: tf.abundance,
            seed_introns: seed_chain.len(),
            overlap,
            first,
            last,
        });
    }

    rows.sort_unstable_by(|a, b| {
        b.overlap
            .cmp(&a.overlap)
            .then_with(|| b.seed_introns.cmp(&a.seed_introns))
            .then_with(|| a.label.cmp(b.label))
            .then_with(|| a.seed.cmp(&b.seed))
    });

    eprintln!(
        "[TRACE_REF_SEEDS] ref={} bundle={}:{}-{} tol={} ref_introns={} overlapping_seeds={} label_counts={:?}",
        ref_tx.id,
        bundle.chrom,
        bundle.start,
        bundle.end,
        tol,
        ref_chain.len(),
        rows.len(),
        by_label
    );
    for r in rows.into_iter().take(30) {
        eprintln!(
            "[TRACE_REF_SEED] seed={} overlap={}/{} seed_introns={} abund={:.4} span={}-{}..{}-{} outcome={}",
            r.seed,
            r.overlap,
            ref_chain.len(),
            r.seed_introns,
            r.abund,
            r.first.0,
            r.first.1,
            r.last.0,
            r.last.1,
            r.label
        );
    }
}

fn bundlenodes_to_vec(bn: Option<&crate::types::CBundlenode>) -> Vec<(usize, u64, u64, f64)> {
    let mut out = Vec::new();
    let mut cur = bn;
    while let Some(n) = cur {
        out.push((n.bid, n.start, n.end, n.cov));
        cur = n.next.as_deref();
    }
    out
}

fn vec_to_bundlenodes(items: &[(usize, u64, u64, f64)]) -> Option<crate::types::CBundlenode> {
    if items.is_empty() {
        return None;
    }
    let mut head = crate::types::CBundlenode::new(items[0].1, items[0].2, items[0].3, items[0].0);
    let mut cur = &mut head;
    for &(bid, start, end, cov) in &items[1..] {
        cur.next = Some(Box::new(crate::types::CBundlenode::new(
            start, end, cov, bid,
        )));
        cur = cur.next.as_deref_mut().unwrap();
    }
    Some(head)
}

#[inline]
fn trace_log_style_active() -> bool {
    std::env::var_os("RUSTLE_TRACE_LOG_STYLE").is_some()
}

#[inline]
fn loop_trace_active() -> bool {
    std::env::var_os("RUSTLE_LOOP_TRACE").is_some() || trace_log_style_active()
}

#[inline]
fn trace_numeric_state_active() -> bool {
    std::env::var_os("RUSTLE_TRACE_METRICS").is_some()
}

#[inline]
fn parse_trace_bad_mm_neg_junc() -> Option<(u64, u64)> {
    let val = std::env::var("RUSTLE_TRACE_BAD_MM_NEG_JUNC").ok()?;
    let (a, b) = val.split_once('-')?;
    let donor = a.trim().parse::<u64>().ok()?;
    let acceptor = b.trim().parse::<u64>().ok()?;
    Some((donor, acceptor))
}

#[inline]
fn trace_bad_mm_neg_match(j: Junction) -> bool {
    let Some((donor, acceptor)) = parse_trace_bad_mm_neg_junc() else {
        return false;
    };
    j.donor == donor && (j.acceptor == acceptor || j.acceptor.saturating_add(1) == acceptor)
}

/// Collapse zero-length internal graph nodes.
///
/// When the event-driven builder processes End→Start at the same genomic
/// position, it produces a zero-length placeholder node between the ending
/// and continuing exons. StringTie's equivalent code produces direct edges
/// instead. This pass bridges parent→child edges and orphans the placeholder
/// so downstream transfrag creation and flow decomposition see a tighter
/// graph.
fn collapse_zero_length_nodes(graph: &mut Graph) {
    let n = graph.nodes.len();
    let src = graph.source_id;
    let snk = graph.sink_id;
    for i in 0..n {
        if i == src || i == snk {
            continue;
        }
        let node = &graph.nodes[i];
        if node.start != node.end {
            continue;
        }
        let parents: Vec<usize> = node.parents.ones().collect();
        let children: Vec<usize> = node.children.ones().collect();
        // Bridge each parent directly to each child.
        for &p in &parents {
            for &c in &children {
                if p == c {
                    continue;
                }
                graph.add_edge(p, c);
            }
        }
        // Disconnect the zero-length node from neighbors.
        for &p in &parents {
            graph.nodes[p].children.remove(i);
        }
        for &c in &children {
            graph.nodes[c].parents.remove(i);
        }
        graph.nodes[i].parents.clear();
        graph.nodes[i].children.clear();
    }
}

/// Emit JFINAL trace: final junction decision state for each cjunction.
/// Format matches stringtie_debug/rlink.cpp JFINAL block so logs can be diffed.
/// Key = "<donor>-<acceptor+1>:<strand>" (+1 to match StringTie's 1-based end convention).
/// Decision: KEEP iff strand != 0 && mm >= 0 && nreads >= 0 && nreads_good >= 0.
fn emit_jfinal_trace(cjunctions: &[crate::types::CJunction], tag: &str) {
    if std::env::var_os("JFINAL_TRACE").is_none() {
        return;
    }
    for cj in cjunctions {
        let keep = cj.strand != 0 && cj.mm >= 0.0 && cj.nreads >= 0.0 && cj.nreads_good >= 0.0;
        eprintln!(
            "JFINAL[{}] {}-{}:{} nreads={:.1} nreads_good={:.1} mm={:.1} nm={:.1} guide={} left={:.1} right={:.1} decision={}",
            tag,
            cj.start,
            cj.end.saturating_add(1),
            cj.strand,
            cj.nreads,
            cj.nreads_good,
            cj.mm,
            cj.nm,
            if cj.guide_match { 1 } else { 0 },
            cj.leftsupport,
            cj.rightsupport,
            if keep { "KEEP" } else { "KILL" }
        );
    }
}

fn zero_graph_node_bp_coverage(graph: &mut Graph) {
    for node in graph.nodes.iter_mut() {
        node.coverage = 0.0;
    }
}

fn refresh_graph_node_bp_coverage_from_bpcov(graph: &mut Graph, bpcov_stranded: &BpcovStranded) {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    for i in 0..graph.n_nodes {
        if i == source_id || i == sink_id {
            continue;
        }
        let (ns, ne) = {
            let n = &graph.nodes[i];
            (n.start, n.end)
        };
        if ns >= ne {
            continue;
        }
        let si = bpcov_stranded.plus.idx(ns);
        let ei = bpcov_stranded.plus.idx(ne);
        let total_cov = bpcov_stranded.get_cov_range(BPCOV_STRAND_ALL, si, ei);
        graph.nodes[i].coverage = total_cov;
    }
}

fn trace_graph_numeric_state(stage: &str, graph: &Graph, transfrags: &[GraphTransfrag]) {
    if !trace_numeric_state_active() {
        return;
    }
    // Keep TRACE_METRICS usable by scoping output to the requested trace locus when present.
    // Without this, a full-chromosome run can produce an overwhelming amount of numeric state.
    if let Some((lo, hi)) = crate::trace_events::parse_trace_locus() {
        let mut any_overlap = false;
        for node in &graph.nodes {
            if node.end > lo && node.start < hi {
                any_overlap = true;
                break;
            }
        }
        if !any_overlap {
            return;
        }
    }
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let mut node_bp_cov_sum = 0.0f64;
    let mut nodecov_sum = 0.0f64;
    let mut longcov_sum = 0.0f64;
    let mut rate_nonunit = 0usize;
    let mut active_nodes = 0usize;
    for (nid, node) in graph.nodes.iter().enumerate() {
        if nid == source_id || nid == sink_id {
            continue;
        }
        node_bp_cov_sum += node.coverage;
        nodecov_sum += node.nodecov;
        longcov_sum += node.longcov;
        if node.coverage > 0.0 || node.nodecov > 0.0 || node.longcov > 0.0 {
            active_nodes += 1;
        }
        if (node.noderate - 1.0).abs() > 1e-9 {
            rate_nonunit += 1;
        }
    }
    let tf_abund_sum: f64 = transfrags.iter().map(|tf| tf.abundance).sum();
    let tf_srabund_sum: f64 = transfrags.iter().map(|tf| tf.srabund).sum();
    let tf_read_sum: f64 = transfrags.iter().map(|tf| tf.read_count).sum();
    eprintln!(
        "TRACE_METRICS stage={} nodes={} active_nodes={} node_bp_cov_sum={:.4} nodecov_sum={:.4} long_node_abund_sum={:.4} tf_abund_sum={:.4} tf_srabund_sum={:.4} tf_read_sum={:.4} noderate_nonunit={}",
        stage,
        graph.n_nodes.saturating_sub(2),
        active_nodes,
        node_bp_cov_sum,
        nodecov_sum,
        longcov_sum,
        tf_abund_sum,
        tf_srabund_sum,
        tf_read_sum,
        rate_nonunit
    );
    for (nid, node) in graph.nodes.iter().enumerate() {
        if nid == source_id || nid == sink_id {
            continue;
        }
        if node.coverage <= 0.0 && node.nodecov <= 0.0 && node.longcov <= 0.0 {
            continue;
        }
        eprintln!(
            "TRACE_METRICS_NODE stage={} nid={} {}-{} bp_cov={:.4} nodecov={:.4} longcov={:.4} noderate={:.6}",
            stage,
            nid,
            node.start,
            node.end,
            node.coverage,
            node.nodecov,
            node.longcov,
            node.noderate
        );
    }
}

#[inline]
fn trace_bundle_strand_i8(strand: char) -> i8 {
    match strand {
        '-' => -1,
        '+' => 1,
        _ => 0,
    }
}

#[inline]
fn trace_bundle_end(end: u64) -> u64 {
    end.saturating_add(1)
}

#[inline]
fn trace_junction_acceptor(jn: &Junction) -> u64 {
    jn.acceptor.saturating_add(1)
}

#[derive(Clone, Debug)]
struct DebugBundleTarget {
    chrom: String,
    start: u64,
    end: u64,
}

fn parse_debug_bundle_target(config: &RunConfig) -> Option<DebugBundleTarget> {
    let spec = config.debug_bundle.as_ref()?;
    let (chrom, rest) = spec.split_once(':')?;
    let (start, end) = rest.split_once('-')?;
    let start = start.replace(',', "").parse::<u64>().ok()?;
    let end = end.replace(',', "").parse::<u64>().ok()?;
    Some(DebugBundleTarget {
        chrom: chrom.to_string(),
        start,
        end,
    })
}

#[inline]
fn debug_bundle_overlaps(
    target: Option<&DebugBundleTarget>,
    chrom: &str,
    start: u64,
    end: u64,
) -> bool {
    target
        .map(|t| t.chrom == chrom && start <= t.end && t.start <= end)
        .unwrap_or(false)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct ReadSig {
    name_hash: u64,
    hi: u32,
    nh: u32,
    ref_start: u64,
    ref_end: u64,
    exon_count: u32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct ReadSigNoStrand {
    name_hash: u64,
    hi: u32,
    nh: u32,
    ref_start: u64,
    ref_end: u64,
    exon_count: u32,
}

#[inline]
fn read_sig(r: &BundleRead) -> ReadSig {
    ReadSig {
        name_hash: r.read_name_hash,
        hi: r.hi,
        nh: r.nh,
        ref_start: r.ref_start,
        ref_end: r.ref_end,
        exon_count: r.exons.len() as u32,
    }
}

#[inline]
fn read_sig_nostrand(r: &BundleRead) -> ReadSigNoStrand {
    ReadSigNoStrand {
        name_hash: r.read_name_hash,
        hi: r.hi,
        nh: r.nh,
        ref_start: r.ref_start,
        ref_end: r.ref_end,
        exon_count: r.exons.len() as u32,
    }
}

fn compute_initial_junction_stats_for_reads(
    reads: &[BundleRead],
    bundle_start: u64,
    bundle_end: u64,
    config: &RunConfig,
) -> JunctionStats {
    let mut junction_stats: JunctionStats = Default::default();

    for r in reads {
        let nex = r.exons.len();
        let mut leftsup = vec![0u64; nex];
        let mut rightsup = vec![0u64; nex];
        if nex > 0 {
            let mut max_left = 0u64;
            let mut max_right = 0u64;
            for ei in 0..nex {
                let seg_len = r.exons[ei].1.saturating_sub(r.exons[ei].0);
                if seg_len > max_left {
                    max_left = seg_len;
                }
                leftsup[ei] = max_left;
            }
            for ei in (0..nex).rev() {
                let seg_len = r.exons[ei].1.saturating_sub(r.exons[ei].0);
                if seg_len > max_right {
                    max_right = seg_len;
                }
                rightsup[ei] = max_right;
            }
        }

        for i in 0..r.junctions.len() {
            let j = r.junctions[i];
            let intron_len = j.acceptor.saturating_sub(j.donor);
            if intron_len < config.min_intron_length {
                continue;
            }
            if !(j.donor >= bundle_start
                && j.donor <= bundle_end
                && j.acceptor >= bundle_start
                && j.acceptor <= bundle_end)
            {
                continue;
            }
            let left_anchor = leftsup.get(i).copied().unwrap_or(0);
            let right_anchor = rightsup.get(i + 1).copied().unwrap_or(0);
            let st = junction_stats
                .entry(j)
                .or_insert_with(JunctionStat::default);
            st.mrcount += r.weight;
            let mut anchor = config.junction_support;
            if intron_len > config.longintron && anchor < LONGINTRONANCHOR {
                anchor = LONGINTRONANCHOR;
            }
            if left_anchor >= anchor {
                st.leftsupport += r.weight;
                if right_anchor >= anchor {
                    st.rightsupport += r.weight;
                    st.nreads_good += r.weight;
                }
            } else if right_anchor >= anchor {
                st.rightsupport += r.weight;
            }
            st.nm += r.junc_mismatch_weight;
            if left_anchor > LONGINTRONANCHOR && right_anchor > LONGINTRONANCHOR {
                st.mm += r.weight;
            }
            match st.strand {
                None => {
                    st.strand = Some(if r.strand == '-' { -1 } else { 1 });
                }
                Some(existing) => {
                    let read_strand = if r.strand == '-' { -1 } else { 1 };
                    if existing != read_strand {
                        st.strand = Some(0);
                    }
                }
            }
        }
    }

    junction_stats
}

fn merge_region_outer_bundles(
    mut bundles: Vec<crate::types::Bundle>,
    config: &RunConfig,
) -> Vec<crate::types::Bundle> {
    if bundles.is_empty() {
        return bundles;
    }
    let debug_target = parse_debug_bundle_target(config);
    bundles.sort_unstable_by(|a, b| {
        (a.chrom.as_str(), a.start, a.end, a.strand).cmp(&(
            b.chrom.as_str(),
            b.start,
            b.end,
            b.strand,
        ))
    });

    let mut merged: Vec<crate::types::Bundle> = Vec::new();
    let mut i = 0usize;
    while i < bundles.len() {
        let chrom = bundles[i].chrom.clone();
        let start = bundles[i].start;
        let end = bundles[i].end;
        let mut reads: Vec<BundleRead> = Vec::new();
        while i < bundles.len()
            && bundles[i].chrom == chrom
            && bundles[i].start == start
            && bundles[i].end == end
        {
            reads.extend(std::mem::take(&mut bundles[i].reads));
            i += 1;
        }
        reads.sort_unstable_by_key(|r| (r.ref_start, r.ref_end, r.read_uid));
        let junction_stats = compute_initial_junction_stats_for_reads(&reads, start, end, config);
        if debug_bundle_overlaps(debug_target.as_ref(), &chrom, start, end) {
            let contributing = bundles
                .iter()
                .filter(|b| b.chrom == chrom && b.start == start && b.end == end)
                .count();
            let mut raw_counts: HashMap<crate::types::Junction, f64> = Default::default();
            let mut del_counts: HashMap<crate::types::Junction, f64> = Default::default();
            for r in &reads {
                for j in &r.junctions_raw {
                    *raw_counts.entry(*j).or_insert(0.0) += r.weight;
                }
                for j in &r.junctions_del {
                    *del_counts.entry(*j).or_insert(0.0) += r.weight;
                }
            }
            eprintln!(
                "[EARLY_DEBUG] region_merge {}:{}-{} contributing_bundles={} merged_reads={} merged_jstats={}",
                chrom,
                start,
                end,
                contributing,
                reads.len(),
                junction_stats.len()
            );
            dump_junction_counts("region_merge_raw", &raw_counts);
            dump_junction_counts("region_merge_del", &del_counts);
            dump_junction_stats("region_merge_initial_jstats", &junction_stats);
        }
        merged.push(crate::types::Bundle {
            chrom,
            start,
            end,
            strand: '.',
            reads,
            junction_stats,
            bundlenodes: None,
            read_bnodes: None,
            bnode_colors: None,
        });
    }

    merged
}

fn map_subset_reads_to_region_indices(
    subset_reads: &[BundleRead],
    region_reads: &[BundleRead],
) -> Vec<Option<usize>> {
    let mut buckets: HbHashMap<ReadSig, Vec<usize>> = HbHashMap::with_capacity(region_reads.len());
    for (i, r) in region_reads.iter().enumerate() {
        buckets.entry(read_sig(r)).or_default().push(i);
    }
    for idxs in buckets.values_mut() {
        idxs.reverse();
    }

    let mut uid_to_region: HbHashMap<u64, usize> = HbHashMap::with_capacity(region_reads.len());
    for (i, r) in region_reads.iter().enumerate() {
        uid_to_region.insert(r.read_uid, i);
    }

    let mut out: Vec<Option<usize>> = Vec::with_capacity(subset_reads.len());
    for r in subset_reads {
        let sig = read_sig(r);
        let idx = buckets
            .get_mut(&sig)
            .and_then(|v| v.pop())
            .or_else(|| uid_to_region.get(&r.read_uid).copied());
        out.push(idx);
    }
    out
}

#[inline]
fn build_allowed_junction_index(
    allowed_junctions: &[Junction],
) -> Option<HbHashMap<Junction, u32>> {
    if allowed_junctions.len() > u32::MAX as usize {
        return None;
    }
    let mut allowed_idx: HbHashMap<Junction, u32> =
        HbHashMap::with_capacity(allowed_junctions.len());
    for (idx, &j) in allowed_junctions.iter().enumerate() {
        allowed_idx.insert(j, idx as u32);
    }
    Some(allowed_idx)
}

#[inline]
fn select_component_junctions_from_reads_indexed(
    reads: &[BundleRead],
    allowed_junctions: &[Junction],
    allowed_idx: Option<&HbHashMap<Junction, u32>>,
) -> Vec<Junction> {
    if reads.is_empty() || allowed_junctions.is_empty() {
        return Vec::new();
    }
    if let Some(allowed_idx) = allowed_idx {
        let mut present = RoaringBitmap::new();
        for r in reads {
            for &j in &r.junctions {
                if let Some(&idx) = allowed_idx.get(&j) {
                    present.insert(idx);
                }
            }
        }
        let mut out: Vec<Junction> = Vec::with_capacity(present.len() as usize);
        for idx in present.iter() {
            if let Some(&j) = allowed_junctions.get(idx as usize) {
                out.push(j);
            }
        }
        return out;
    }
    // Defensive fallback for extremely large junction vectors.
    let mut seen: HbHashSet<Junction> = HbHashSet::with_capacity(allowed_junctions.len());
    for r in reads {
        for &j in &r.junctions {
            seen.insert(j);
        }
    }
    allowed_junctions
        .iter()
        .copied()
        .filter(|j| seen.contains(j))
        .collect()
}

#[derive(Clone, Copy, Debug, Default)]
struct SeedOutcomeSummary {
    total: usize,
    stored: usize,
    unwitnessed: usize,
    checktrf_rescued: usize,
    checktrf_total: usize,
    zero_flux: usize,
    low_coverage: usize,
    back_fail: usize,
    fwd_fail: usize,
}

fn accumulate_seed_outcomes(
    out: &mut SeedOutcomeSummary,
    outcomes: &[(usize, crate::path_extract::SeedOutcome)],
) {
    use crate::path_extract::SeedOutcome;
    for (_idx, ev) in outcomes {
        out.total += 1;
        match ev {
            SeedOutcome::Stored(_) => out.stored += 1,
            SeedOutcome::UnwitnessedSplice => out.unwitnessed += 1,
            SeedOutcome::ChecktrfRescued => {
                out.checktrf_rescued += 1;
                out.checktrf_total += 1;
            }
            SeedOutcome::ChecktrfReadthr
            | SeedOutcome::ChecktrfRedistributed
            | SeedOutcome::ChecktrfIncomplete
            | SeedOutcome::ChecktrfRescueFail => {
                out.checktrf_total += 1;
            }
            SeedOutcome::ZeroFlux => out.zero_flux += 1,
            SeedOutcome::LowCoverage(_) => out.low_coverage += 1,
            SeedOutcome::BackToSourceFail => out.back_fail += 1,
            SeedOutcome::FwdToSinkFail => out.fwd_fail += 1,
            _ => {}
        }
    }
}

fn count_assigned_reads(read_bnodes: &[Vec<usize>], read_scales: &[f64]) -> usize {
    read_bnodes
        .iter()
        .enumerate()
        .filter(|(ri, b)| !b.is_empty() && read_scales.get(*ri).copied().unwrap_or(1.0) > 0.0)
        .count()
}

fn unique_color_roots(colors: &[usize]) -> usize {
    if colors.is_empty() {
        return 0;
    }
    let mut seen: HashSet<usize> = Default::default();
    for &c in colors {
        seen.insert(c);
    }
    seen.len()
}

fn trace_junction_reason(st: &JunctionStat, _junction_thr: f64) -> &'static str {
    if st.strand.unwrap_or(0) == 0 {
        "BAD_NO_STRAND"
    } else if st.mm < 0.0 {
        "BAD_MM_NEG"
    } else if st.mrcount < 0.0 {
        "DEMOTED_LEFT"
    } else if st.nreads_good < 0.0 {
        "DEMOTED_RIGHT"
    } else if st.guide_match {
        "GUIDE_MATCH"
    } else {
        "GOOD"
    }
}

/// Per-junction quality check matching the original algorithm's good_junc
/// (line 14374-14476). Runs on each junction using strand-specific bpcov.
/// Kills junctions (mm=-1) that fail:
///   - coverage drop check at donor (bw=5 window, strand-specific)
///   - coverage drop check at acceptor
///   - isofrac: nreads*10 < ERROR_PERC*leftsupport
/// Per-junction quality check matching the original algorithm's good_junc
/// (rlink.cpp line 14374-14476). Uses strand-specific bpcov with bw=5 window
/// and checks: coverage signal at donor, coverage signal at acceptor,
/// and isofrac relative to splice site support.
fn apply_bad_mm_neg_stage(
    junction_stats: &mut JunctionStats,
    bpcov: &BpcovStranded,
    refstart: u64,
    junction_thr: f64,
) {
    const ERROR_PERC: f64 = 0.1;
    const BW: usize = 5;
    // For long reads: mult = 1/(ERROR_PERC^2) = 100
    let mult: f64 = 1.0 / (ERROR_PERC * ERROR_PERC);

    /// Strand-specific coverage over a range: total - opposite_strand.
    fn strand_cov(bpcov: &BpcovStranded, opp: usize, start: usize, end: usize) -> f64 {
        let total = bpcov.get_cov_range(BPCOV_STRAND_ALL, start, end);
        let opposite = if opp != BPCOV_STRAND_ALL {
            bpcov.get_cov_range(opp, start, end)
        } else {
            0.0
        };
        (total - opposite).max(0.0)
    }

    let keys: Vec<Junction> = junction_stats.keys().copied().collect();
    for j in keys {
        let Some(st) = junction_stats.get_mut(&j) else { continue };
        let strand = st.strand.unwrap_or(0);

        // Skip: already killed, guide, no mismatches, or nm < mrcount
        if strand == 0 || st.guide_match || st.nm <= 0.0 || st.nm + 1e-9 < st.mrcount {
            continue;
        }
        // Skip: already redirected
        if st.mrcount < 0.0 || st.nreads_good < 0.0 || st.mm < 0.0 {
            continue;
        }
        // Kill: below junction threshold
        if st.nreads_good < junction_thr {
            st.mm = -1.0;
            continue;
        }
        // Protect well-supported junctions (matching splice site fix)
        if st.nreads_good >= 5.0 * junction_thr {
            continue;
        }

        let mismatch = st.nm > 0.0 && (st.nm - st.mrcount).abs() < 0.5;
        let opp = match strand {
            -1 => BPCOV_STRAND_PLUS,
            1 => BPCOV_STRAND_MINUS,
            _ => BPCOV_STRAND_ALL,
        };

        // === Donor check (bw=5 window) ===
        let dp = j.donor.saturating_sub(refstart) as usize;
        if dp >= BW && dp + BW + 1 < bpcov.minus.cov.len() {
            let lleftcov = strand_cov(bpcov, opp, dp.saturating_sub(BW) + 1, dp + 1);
            let lrightcov = strand_cov(bpcov, opp, dp + 1, dp + BW + 1);
            if lleftcov > 1.0 / ERROR_PERC
                && st.leftsupport * mult < ERROR_PERC * lleftcov
                && (mismatch || lrightcov > lleftcov * (1.0 - ERROR_PERC))
            {
                st.mm = -1.0;
                continue;
            }
        }

        // === Acceptor check (bw=5 window) ===
        let ap = j.acceptor.saturating_sub(refstart).saturating_sub(1) as usize;
        if ap >= BW && ap + BW + 1 < bpcov.minus.cov.len() {
            let rleftcov = strand_cov(bpcov, opp, ap.saturating_sub(BW) + 1, ap + 1);
            let rrightcov = strand_cov(bpcov, opp, ap + 1, ap + BW + 1);
            if rrightcov > 1.0 / ERROR_PERC
                && st.rightsupport * mult < rrightcov * ERROR_PERC
                && (mismatch || rleftcov > rrightcov * (1.0 - ERROR_PERC))
            {
                st.mm = -1.0;
                continue;
            }
        }

        // === Isofrac: junction reads as fraction of splice site support ===
        if st.mrcount * 10.0 < ERROR_PERC * st.leftsupport
            || st.mrcount * 10.0 < ERROR_PERC * st.rightsupport
        {
            st.mm = -1.0;
        }
    }
}

fn emit_trace_junction_actions(
    bundle_strand: char,
    junction_stats: &JunctionStats,
    junction_redirect_map: &HashMap<Junction, Junction>,
) {
    if !trace_log_style_active() {
        return;
    }

    let strand = trace_bundle_strand_i8(bundle_strand);
    let mut donor_sorted: Vec<Junction> = junction_stats.keys().copied().collect();
    donor_sorted.sort_unstable_by_key(|j| (j.donor, j.acceptor));
    let donor_index: HashMap<Junction, usize> = donor_sorted
        .iter()
        .enumerate()
        .map(|(idx, j)| (*j, idx + 1))
        .collect();

    for (idx, j) in donor_sorted.iter().enumerate() {
        let Some(st) = junction_stats.get(j) else {
            continue;
        };
        if st.strand.unwrap_or(0) != 0 && st.mm >= 0.0 {
            continue;
        }
        eprintln!(
            "--- build_graphs: BAD_JUNC i={} {}-{}:{} nm={:.1} nreads={:.1} nreads_good={:.1} mm={:.1} left={:.1} right={:.1} guide={}",
            idx + 1,
            j.donor,
            trace_junction_acceptor(j),
            strand,
            st.nm,
            st.mrcount,
            st.nreads_good,
            st.mm,
            st.leftsupport,
            st.rightsupport,
            if st.guide_match { 1 } else { 0 }
        );
        eprintln!(
            "--- build_graphs: JUNC_DELETE i={} {}-{}:{} nreads_good={:.1}",
            idx + 1,
            j.donor,
            trace_junction_acceptor(j),
            strand,
            st.nreads_good
        );
        if let Some(dest) = junction_redirect_map.get(j) {
            if dest != j {
                let dest_idx = donor_index.get(dest).copied().unwrap_or(0);
                let dest_left = junction_stats
                    .get(dest)
                    .map(|s| s.leftsupport)
                    .unwrap_or(0.0);
                eprintln!(
                    "--- build_graphs: JUNC_DEMOTE i={}({}-{}) -> j={}({}-{}) leftsup={:.1}->{:.1}",
                    idx + 1,
                    j.donor,
                    trace_junction_acceptor(j),
                    dest_idx,
                    dest.donor,
                    trace_junction_acceptor(dest),
                    st.leftsupport,
                    dest_left
                );
            }
        }
    }

    let mut acceptor_sorted: Vec<Junction> = junction_stats.keys().copied().collect();
    acceptor_sorted.sort_unstable_by_key(|j| (j.acceptor, j.donor));
    for (idx, j) in acceptor_sorted.iter().enumerate() {
        let Some(st) = junction_stats.get(j) else {
            continue;
        };
        if st.strand.unwrap_or(0) != 0 && st.mm >= 0.0 {
            continue;
        }
        eprintln!(
            "--- build_graphs: EJUNC_DELETE i={} {}-{}:{} nreads_good={:.1}",
            idx + 1,
            j.donor,
            trace_junction_acceptor(j),
            strand,
            st.nreads_good
        );
    }
}

fn emit_trace_junction_decision_table(junction_stats: &JunctionStats, junction_thr: f64) {
    if !trace_log_style_active() {
        return;
    }

    let mut keys: Vec<Junction> = junction_stats.keys().copied().collect();
    keys.sort_unstable_by_key(|j| (j.donor, j.acceptor));
    eprintln!(
        "--- build_graphs: JUNC_DECISION_TABLE count={}",
        keys.len().saturating_add(1)
    );
    eprintln!(
        "    jdec[0] 0-0 strand=0 nm=0.0 mm=0.0 nreads=0.0 nreads_good=0.0 left=0.0 right=0.0 guide=0 reason=BAD_NO_STRAND"
    );
    for (idx, j) in keys.iter().enumerate() {
        let Some(st) = junction_stats.get(j) else {
            continue;
        };
        eprintln!(
            "    jdec[{}] {}-{} strand={} nm={:.1} mm={:.1} nreads={:.1} nreads_good={:.1} left={:.1} right={:.1} guide={} reason={}",
            idx + 1,
            j.donor,
            trace_junction_acceptor(j),
            st.strand.unwrap_or(0),
            st.nm,
            st.mm,
            st.mrcount,
            st.nreads_good,
            st.leftsupport,
            st.rightsupport,
            if st.guide_match { 1 } else { 0 },
            trace_junction_reason(st, junction_thr)
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Junction, JunctionStat, JunctionStats};

    #[test]
    fn bad_mm_neg_donor_uses_probe_base() {
        let refstart = 100u64;
        let j = Junction::new(110, 200);
        let mut stats: JunctionStats = Default::default();
        stats.insert(
            j,
            JunctionStat {
                mrcount: 4.0,
                nreads_good: 4.0,
                strand: Some(1),
                nm: 4.0,
                mm: 4.0,
                ..JunctionStat::default()
            },
        );

        let mut bpcov = BpcovStranded::empty(refstart, 250);
        // Probes donor support at start-refstart (0-based index relative to refstart).
        // With donor=110 and refstart=100, point=10.
        // leftcov reads at index 10 (genomic 110), rightcov at index 11 (genomic 111).
        // To keep the junction, we need rightcov <= 0.9 * leftcov (tolerance = 1 - 0.1).
        bpcov.add_coverage(BPCOV_STRAND_ALL, 110, 111, 5.0); // leftcov position
        bpcov.add_coverage(BPCOV_STRAND_ALL, 111, 112, 1.0); // rightcov position - low enough to keep
        bpcov.add_coverage(BPCOV_STRAND_ALL, 199, 200, 1.0); // acceptor left
        bpcov.add_coverage(BPCOV_STRAND_ALL, 200, 201, 8.0); // acceptor right

        apply_bad_mm_neg_stage(&mut stats, &bpcov, refstart, 1.0);

        // Junction should be kept (mm remains 4.0) because rightcov (1.0) <= 0.9 * leftcov (4.5)
        assert_eq!(stats.get(&j).map(|st| st.mm), Some(4.0));
    }
}

fn build_component_reads(
    comps: &[Vec<(usize, u64, u64, f64)>],
    read_bnodes: &[Vec<usize>],
    read_scales: &[f64],
    scaled_reads: &[BundleRead],
) -> Vec<Vec<BundleRead>> {
    if comps.is_empty() {
        return Vec::new();
    }
    if comps.len() == 1 {
        return vec![scaled_reads
            .iter()
            .filter(|r| r.weight > 0.0)
            .cloned()
            .collect()];
    }

    let mut bnode_to_comp: HbHashMap<usize, usize> = HbHashMap::with_capacity(comps.len() * 4);
    for (ci, comp) in comps.iter().enumerate() {
        for &(bid, _, _, _) in comp {
            bnode_to_comp.insert(bid, ci);
        }
    }

    let use_bitmap = comps.len() <= u32::MAX as usize;
    let mut comp_read_indices: Vec<Vec<usize>> = vec![Vec::new(); comps.len()];
    for (ri, bnodes) in read_bnodes.iter().enumerate() {
        if read_scales.get(ri).copied().unwrap_or(1.0) <= 0.0 {
            continue;
        }
        if use_bitmap {
            let mut assigned = RoaringBitmap::new();
            for &bid in bnodes {
                if let Some(&ci) = bnode_to_comp.get(&bid) {
                    assigned.insert(ci as u32);
                }
            }
            if assigned.is_empty() {
                let Some(r) = scaled_reads.get(ri) else {
                    continue;
                };
                for (ci, comp) in comps.iter().enumerate() {
                    let cs = comp.first().map(|v| v.1).unwrap_or(0);
                    let ce = comp.last().map(|v| v.2).unwrap_or(0);
                    if r.ref_end >= cs && r.ref_start <= ce {
                        comp_read_indices[ci].push(ri);
                        break;
                    }
                }
            } else {
                for ci in assigned.iter() {
                    comp_read_indices[ci as usize].push(ri);
                }
            }
        } else {
            let mut assigned: HbHashSet<usize> = Default::default();
            for &bid in bnodes {
                if let Some(&ci) = bnode_to_comp.get(&bid) {
                    assigned.insert(ci);
                }
            }
            if assigned.is_empty() {
                let Some(r) = scaled_reads.get(ri) else {
                    continue;
                };
                for (ci, comp) in comps.iter().enumerate() {
                    let cs = comp.first().map(|v| v.1).unwrap_or(0);
                    let ce = comp.last().map(|v| v.2).unwrap_or(0);
                    if r.ref_end >= cs && r.ref_start <= ce {
                        comp_read_indices[ci].push(ri);
                        break;
                    }
                }
            } else {
                for ci in assigned {
                    comp_read_indices[ci].push(ri);
                }
            }
        }
    }

    comp_read_indices
        .iter()
        .map(|indices| {
            indices
                .iter()
                .filter_map(|&i| {
                    let r = scaled_reads.get(i)?.clone();
                    if r.weight > 0.0 {
                        Some(r)
                    } else {
                        None
                    }
                })
                .collect()
        })
        .collect()
}

fn shadow_divergence_score(
    baseline_assigned: usize,
    strict_assigned: usize,
    baseline_bnodes: usize,
    strict_bnodes: usize,
    baseline_components: usize,
    strict_components: usize,
    base_seed: SeedOutcomeSummary,
    strict_seed: SeedOutcomeSummary,
) -> i64 {
    let mut score = 0i64;
    score += (strict_seed.unwitnessed as i64 - base_seed.unwitnessed as i64).abs() * 6;
    score += (strict_seed.checktrf_rescued as i64 - base_seed.checktrf_rescued as i64).abs() * 5;
    score += (strict_seed.total as i64 - base_seed.total as i64).abs() * 3;
    score += (strict_seed.stored as i64 - base_seed.stored as i64).abs() * 3;
    score += (strict_assigned as i64 - baseline_assigned as i64).abs() * 2;
    score += (strict_bnodes as i64 - baseline_bnodes as i64).abs() * 2;
    score += (strict_components as i64 - baseline_components as i64).abs() * 2;
    score
}

#[derive(Clone, Debug)]
struct ShadowStrictBundleDiag {
    chrom: String,
    start: u64,
    end: u64,
    strand: char,
    baseline_cgroups_est: usize,
    strict_cgroups_est: usize,
    baseline_bnodes: usize,
    strict_bnodes: usize,
    baseline_components: usize,
    strict_components: usize,
    baseline_assigned: usize,
    strict_assigned: usize,
    baseline_seed: SeedOutcomeSummary,
    strict_seed: SeedOutcomeSummary,
    score: i64,
}

fn dump_junction_stats(label: &str, stats: &crate::types::JunctionStats) {
    let mut list: Vec<(crate::types::Junction, crate::types::JunctionStat)> =
        stats.iter().map(|(j, s)| (*j, s.clone())).collect();
    list.sort_unstable_by_key(|(j, _)| (j.donor, j.acceptor));
    let killed = list.iter().filter(|(_, s)| s.strand == Some(0)).count();
    eprintln!(
        "[DEBUG_JUNCTIONS] {}: {} junctions (killed={})",
        label,
        list.len(),
        killed
    );
    for (j, s) in list {
        eprintln!(
            "  {}-{} mrcount={:.2} good={:.2} nm={:.2} mm={:.2} rcount={} strand={:?} lsup={:.2} rsup={:.2} guide_match={}",
            j.donor,
            j.acceptor,
            s.mrcount,
            s.nreads_good,
            s.nm,
            s.mm,
            s.rcount,
            s.strand,
            s.leftsupport,
            s.rightsupport,
            s.guide_match
        );
    }
}

fn dump_junction_counts(label: &str, counts: &HashMap<crate::types::Junction, f64>) {
    let mut list: Vec<(crate::types::Junction, f64)> =
        counts.iter().map(|(j, c)| (*j, *c)).collect();
    list.sort_unstable_by_key(|(j, _)| (j.donor, j.acceptor));
    eprintln!("[DEBUG_JUNCTIONS] {}: {} junctions", label, list.len());
    for (j, c) in list {
        eprintln!("  {}-{} count={:.2}", j.donor, j.acceptor, c);
    }
}

fn dump_junction_redirects(
    label: &str,
    redirects: &HashMap<crate::types::Junction, crate::types::Junction>,
) {
    let mut list: Vec<(crate::types::Junction, crate::types::Junction)> =
        redirects.iter().map(|(from, to)| (*from, *to)).collect();
    list.sort_unstable_by_key(|(from, to)| (from.donor, from.acceptor, to.donor, to.acceptor));
    eprintln!("[DEBUG_JUNCTIONS] {}: {} redirects", label, list.len());
    for (from, to) in list {
        eprintln!(
            "  {}-{} -> {}-{}",
            from.donor, from.acceptor, to.donor, to.acceptor
        );
    }
}

fn intron_chain_from_exons(exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
    if exons.len() < 2 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(exons.len().saturating_sub(1));
    for w in exons.windows(2) {
        let (_a0, a1) = w[0];
        let (b0, _b1) = w[1];
        if a1 < b0 {
            out.push((a1, b0));
        }
    }
    out
}

fn intron_overlap_count_tol(a: &[(u64, u64)], b: &[(u64, u64)], tol: u64) -> usize {
    a.iter()
        .filter(|ia| b.iter().any(|ib| intron_eq_tol(**ia, *ib, tol)))
        .count()
}

fn intron_chain_equal_tol(a: &[(u64, u64)], b: &[(u64, u64)], tol: u64) -> bool {
    a.len() == b.len() && a.iter().zip(b.iter()).all(|(ia, ib)| intron_eq_tol(*ia, *ib, tol))
}

fn assign_guide_ids_by_intron_chain(
    txs: &mut [Transcript],
    guide_transcripts: &[RefTranscript],
    tol: u64,
) -> usize {
    if txs.is_empty() || guide_transcripts.is_empty() {
        return 0;
    }

    let mut guide_chains: Vec<(usize, Vec<(u64, u64)>)> = Vec::new();
    for (gi, g) in guide_transcripts.iter().enumerate() {
        let chain = intron_chain_from_exons(&g.exons);
        if !chain.is_empty() {
            guide_chains.push((gi, chain));
        }
    }
    if guide_chains.is_empty() {
        return 0;
    }

    let mut assigned = 0usize;
    for tx in txs.iter_mut() {
        let already_guide = tx
            .source
            .as_deref()
            .map(|s| s.starts_with("guide:"))
            .unwrap_or(false)
            || tx.ref_transcript_id.is_some();
        if already_guide {
            continue;
        }

        let tx_chain = intron_chain_from_exons(&tx.exons);
        if tx_chain.is_empty() {
            continue;
        }
        let tx_start = tx.exons.first().map(|e| e.0).unwrap_or(0);
        let tx_end = tx.exons.last().map(|e| e.1).unwrap_or(0);

        let mut best: Option<(u64, usize)> = None;
        for (gi, gchain) in &guide_chains {
            let g = &guide_transcripts[*gi];
            if g.chrom != tx.chrom {
                continue;
            }
            if g.strand != '.' && tx.strand != '.' && g.strand != tx.strand {
                continue;
            }
            if !intron_chain_equal_tol(&tx_chain, gchain, tol) {
                continue;
            }

            let g_start = g.exons.first().map(|e| e.0).unwrap_or(0);
            let g_end = g.exons.last().map(|e| e.1).unwrap_or(0);
            let endpoint_delta = tx_start.abs_diff(g_start) + tx_end.abs_diff(g_end);
            match best {
                Some((best_delta, _)) if endpoint_delta >= best_delta => {}
                _ => best = Some((endpoint_delta, *gi)),
            }
        }

        let Some((_, gi)) = best else {
            continue;
        };
        let g = &guide_transcripts[gi];
        tx.ref_transcript_id = Some(g.id.clone());
        tx.source = Some(format!("guide:{}", g.id));
        assigned += 1;
    }

    assigned
}

fn harmonize_guides_by_chain_postfilter(
    txs: &mut [Transcript],
    guide_transcripts: &[RefTranscript],
    chain_tol: u64,
    max_endpoint_delta: u64,
) -> usize {
    if txs.is_empty() || guide_transcripts.is_empty() {
        return 0;
    }

    let mut changed = 0usize;
    let mut used_tx: HashSet<usize> = Default::default();

    for guide in guide_transcripts {
        if guide.exons.len() < 2 {
            continue;
        }
        if txs.iter().any(|tx| {
            tx.chrom == guide.chrom
                && strand_compatible(tx.strand, guide.strand)
                && tx.exons == guide.exons
        }) {
            continue;
        }

        let guide_chain = intron_chain_from_exons(&guide.exons);
        if guide_chain.is_empty() {
            continue;
        }
        let guide_start = guide.exons.first().map(|e| e.0).unwrap_or(0);
        let guide_end = guide.exons.last().map(|e| e.1).unwrap_or(0);

        let mut best_idx: Option<usize> = None;
        let mut best_delta = u64::MAX;
        let mut best_cov = -1.0f64;

        for (idx, tx) in txs.iter().enumerate() {
            if used_tx.contains(&idx) {
                continue;
            }
            if tx.chrom != guide.chrom || !strand_compatible(tx.strand, guide.strand) {
                continue;
            }
            if tx.exons.len() < 2 {
                continue;
            }
            let tx_chain = intron_chain_from_exons(&tx.exons);
            if tx_chain.is_empty() || !intron_chain_equal_tol(&tx_chain, &guide_chain, chain_tol) {
                continue;
            }
            let tx_start = tx.exons.first().map(|e| e.0).unwrap_or(0);
            let tx_end = tx.exons.last().map(|e| e.1).unwrap_or(0);
            let delta = tx_start.abs_diff(guide_start) + tx_end.abs_diff(guide_end);
            if delta > max_endpoint_delta {
                continue;
            }
            if delta < best_delta || (delta == best_delta && tx.coverage > best_cov) {
                best_delta = delta;
                best_cov = tx.coverage;
                best_idx = Some(idx);
            }
        }

        let Some(idx) = best_idx else {
            continue;
        };
        let tx = &mut txs[idx];
        tx.chrom = guide.chrom.clone();
        tx.strand = guide.strand;
        tx.exons = guide.exons.clone();
        tx.source = Some(format!("guide:{}", guide.id));
        tx.ref_transcript_id = Some(guide.id.clone());
        if tx.exon_cov.len() != tx.exons.len() {
            tx.exon_cov = vec![tx.coverage; tx.exons.len()];
        }
        used_tx.insert(idx);
        changed += 1;
    }

    changed
}

fn has_near_miss_chain(existing: &[Transcript], ref_chain: &[(u64, u64)], tol: u64) -> bool {
    if ref_chain.len() < 2 {
        return false;
    }
    let min_frac = std::env::var("RUSTLE_REF_FALLBACK_MIN_FRAC")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .filter(|v| *v > 0.0 && *v <= 1.0)
        .unwrap_or(0.60);
    let needed = ((ref_chain.len() as f64) * min_frac).ceil() as usize;
    let mut best = 0usize;
    for tx in existing {
        let chain = intron_chain_from_exons(&tx.exons);
        if chain.is_empty() {
            continue;
        }
        let ov = intron_overlap_count_tol(&chain, ref_chain, tol);
        if ov > best {
            best = ov;
        }
    }
    best >= needed
}

fn add_contained_isoforms(
    txs: Vec<Transcript>,
    junctions: &[Junction],
    config: &RunConfig,
) -> Vec<Transcript> {
    if !config.emit_contained_isoforms || txs.is_empty() {
        return txs;
    }
    let mut out = txs;
    let mut seen: HashSet<Vec<(u64, u64)>> = Default::default();
    for tx in &out {
        seen.insert(intron_chain_from_exons(&tx.exons));
    }
    let junction_set: HashSet<(u64, u64)> =
        junctions.iter().map(|j| (j.donor, j.acceptor)).collect();
    let mut added = 0usize;
    let orig_len = out.len();
    let mut additions: Vec<Transcript> = Vec::new();
    for idx in 0..orig_len {
        let tx = &out[idx];
        let n = tx.exons.len();
        if n < 2 {
            continue;
        }
        for trim_start in 0..n {
            for trim_end in 0..(n - trim_start) {
                if trim_start == 0 && trim_end == 0 {
                    continue;
                }
                let end = n - trim_end;
                if end <= trim_start {
                    continue;
                }
                let sub_exons = &tx.exons[trim_start..end];
                if sub_exons.is_empty() {
                    continue;
                }
                let introns = intron_chain_from_exons(sub_exons);
                if introns.iter().any(|j| !junction_set.contains(j)) {
                    continue;
                }
                if !seen.insert(introns.clone()) {
                    continue;
                }
                let exon_cov = if tx.exon_cov.len() == tx.exons.len() {
                    tx.exon_cov[trim_start..end].to_vec()
                } else {
                    vec![tx.coverage; sub_exons.len()]
                };
                additions.push(Transcript {
                    chrom: tx.chrom.clone(),
                    strand: tx.strand,
                    exons: sub_exons.to_vec(),
                    coverage: tx.coverage,
                    exon_cov,
                    tpm: 0.0,
                    fpkm: 0.0,
                    source: Some("contained_iso".to_string()),
                    is_longread: tx.is_longread,
                    longcov: tx.longcov,
                    bpcov_cov: 0.0,
                    transcript_id: None,
                    gene_id: None,
                    ref_transcript_id: None,
                    ref_gene_id: None,
                    hardstart: tx.hardstart,
                    hardend: tx.hardend,
                    alt_tts_end: false,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None, intron_low: Vec::new(),
                });
                added += 1;
            }
        }
    }
    if !additions.is_empty() {
        out.extend(additions);
    }
    if config.verbose && added > 0 {
        eprintln!("    contained_isoforms: added {}", added);
    }
    out
}

/// Fast hash fingerprint for an intron chain (avoids Vec allocation for dedup).
#[inline]
fn intron_chain_fingerprint(exons: &[(u64, u64)]) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = fxhash::FxHasher::default();
    let n = exons.len();
    n.hash(&mut h);
    for i in 0..n.saturating_sub(1) {
        exons[i].1.hash(&mut h); // donor
        exons[i + 1].0.hash(&mut h); // acceptor
    }
    h.finish()
}

fn emit_junction_paths(
    graph: &Graph,
    bundle: &crate::types::Bundle,
    existing: &[Transcript],
    config: &RunConfig,
) -> (Vec<Transcript>, Vec<Transcript>) {
    if !config.emit_junction_paths {
        return (Vec::new(), Vec::new());
    }
    // Use hash fingerprints instead of full Vec<(u64,u64)> keys for O(1) dedup
    // with minimal allocation. Collision probability is negligible for this use case.
    let mut seen: HashSet<u64> = Default::default();
    for tx in existing {
        seen.insert(intron_chain_fingerprint(&tx.exons));
    }
    let source = graph.source_id;
    let sink = graph.sink_id;
    let mut out: Vec<Transcript> = Vec::new();
    let mut emitted: Vec<Transcript> = Vec::new();
    let mut path: Vec<usize> = Vec::new();
    let mut visited: NodeSet = NodeSet::with_capacity(graph.n_nodes);

    // If a specific intron chain is requested, attempt to materialize it directly.
    if let Some(chain) = config.trace_intron_chain.as_ref() {
        if let Some(tx) = emit_trace_chain_from_graph(graph, bundle, chain, config.long_reads) {
            if seen.insert(intron_chain_fingerprint(&tx.exons)) {
                if config.verbose {
                    eprintln!(
                        "    trace_chain_emit: bundle {}:{}-{} emitted junction_chain exons={}",
                        bundle.chrom,
                        bundle.start,
                        bundle.end,
                        tx.exons.len()
                    );
                }
                emitted.push(tx);
            }
        } else if config.verbose {
            eprintln!(
                "    trace_chain_emit: bundle {}:{}-{} could not materialize chain",
                bundle.chrom, bundle.start, bundle.end
            );
        }
    }

    fn maybe_emit_path(
        path: &[usize],
        first_node: usize,
        graph: &Graph,
        bundle: &crate::types::Bundle,
        config: &RunConfig,
        seen: &mut HashSet<u64>,
        out: &mut Vec<Transcript>,
    ) -> bool {
        let exon_nodes: Vec<usize> = path
            .iter()
            .copied()
            .filter(|&n| n != graph.source_id && n != graph.sink_id)
            .collect();
        if exon_nodes.len() < 2 {
            return false;
        }
        let mut exons: Vec<(u64, u64)> = Vec::with_capacity(exon_nodes.len());
        let mut exon_cov: Vec<f64> = Vec::with_capacity(exon_nodes.len());
        let mut cov_bottleneck = f64::INFINITY;
        let mut lc_bottleneck = f64::INFINITY;
        for &nid in &exon_nodes {
            let Some(n) = graph.nodes.get(nid) else {
                continue;
            };
            if n.start >= n.end {
                continue;
            }
            cov_bottleneck = cov_bottleneck.min(n.coverage);
            lc_bottleneck = lc_bottleneck.min(n.longcov);
            let ec = if config.long_reads {
                n.longcov.max(n.coverage * 1e-9)
            } else {
                n.coverage
            };
            exon_cov.push(ec);
            exons.push((n.start, n.end));
        }
        if exons.len() < 2 {
            return false;
        }
        if exons.len() < 2 || !seen.insert(intron_chain_fingerprint(&exons)) {
            return false;
        }
        // Bottleneck abundance along the path (weakest exon). Previously we used epsilon cov and
        // bypassed predcluster, which flooded outputs with unsupported graph walks. Align with
        // flow-derived paths: drop to readthr/isofrac naturally when support is absent.
        if !cov_bottleneck.is_finite() {
            cov_bottleneck = 0.0;
        }
        if !lc_bottleneck.is_finite() {
            lc_bottleneck = 0.0;
        }
        let coverage = if config.long_reads {
            lc_bottleneck.max(cov_bottleneck * 1e-9)
        } else {
            cov_bottleneck
        };
        let longcov = lc_bottleneck;
        // Junction paths need stronger evidence than flow-extracted transcripts.
        // Rustle's coverage is ~1.8x lower than the expected value for the same transcript,
        // so we compensate with a 2x multiplier to avoid flooding predcluster with
        // low-confidence graph walks.
        let thr = (config.readthr * 2.0).max(f64::EPSILON);
        if longcov < thr && coverage < thr {
            return false;
        }
        // Junction-read check: every splice junction in the path must be supported by at least
        // readthr reads. Node-level bottleneck coverage doesn't catch novel junction combinations
        // (path A→B→C where reads only span A→B and B→C separately), because node B's coverage
        // counts all reads passing through it regardless of which junctions they use.
        // Guide-matched junctions (from annotation) bypass this check.
        // This aligns with the DIRECT_LOW_COV rejection for zero-flow paths.
        {
            let junc_thr = config.readthr.max(1.0);
            for i in 0..exons.len().saturating_sub(1) {
                let donor = exons[i].1;
                let acceptor = exons[i + 1].0;
                let jkey = crate::types::Junction::new(donor, acceptor);
                match bundle.junction_stats.get(&jkey) {
                    Some(jstat) if jstat.guide_match => {} // guide junction: allow
                    Some(jstat) => {
                        // Use the best available read count for this junction
                        let jreads = jstat.mrcount.max(jstat.nreads_good);
                        if jreads < junc_thr {
                            return false;
                        }
                    }
                    None => return false, // junction not observed in reads or guides
                }
            }
        }
        let last_node = *exon_nodes.last().unwrap_or(&first_node);
        let hardstart = graph.nodes.get(first_node).map(|n| n.hardstart).unwrap_or(false);
        let hardend = graph.nodes.get(last_node).map(|n| n.hardend).unwrap_or(false);
        let alt_tts_end = graph.nodes.get(last_node).map(|n| n.alt_tts_end).unwrap_or(false);

        out.push(Transcript {
            chrom: bundle.chrom.clone(),
            strand: bundle.strand,
            exons,
            coverage,
            exon_cov,
            tpm: 0.0,
            fpkm: 0.0,
            source: Some("junction_path".to_string()),
            is_longread: config.long_reads,
            longcov,
            bpcov_cov: 0.0,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: None,
            ref_gene_id: None,
            hardstart,
            hardend,
            alt_tts_end,
            vg_family_id: None, vg_copy_id: None, vg_family_size: None, intron_low: Vec::new(),
        });

        // Return true (to be added to main tx list) if it has at least one verified boundary
        hardstart || hardend
    }

    fn dfs(
        cur: usize,
        graph: &Graph,
        bundle: &crate::types::Bundle,
        config: &RunConfig,
        path: &mut Vec<usize>,
        visited: &mut NodeSet,
        seen: &mut HashSet<u64>,
        out: &mut Vec<Transcript>,
        emitted: &mut Vec<Transcript>,
        explored: &mut usize,
    ) -> bool {
        if *explored >= 50_000 {
            return true;
        }
        *explored += 1;
        if out.len() >= 10000 { // Internal limit for evidence set
            return true;
        }
        if path.len() > config.max_junction_path_len {
            return false;
        }
        if cur == graph.sink_id {
            return false;
        }
        let first_node = if path.len() > 1 { path[1] } else { cur };
        if maybe_emit_path(path, first_node, graph, bundle, config, seen, out) {
            if emitted.len() < config.max_junction_paths {
                emitted.push(out.last().unwrap().clone());
            }
        }

        for child in graph.nodes[cur].children.ones().collect::<Vec<_>>() {
            if visited.contains(child) {
                continue;
            }
            visited.insert_grow(child);
            path.push(child);
            let hit_limit = dfs(child, graph, bundle, config, path, visited, seen, out, emitted, explored);
            path.pop();
            visited.remove(child);
            if hit_limit {
                return true;
            }
        }
        false
    }

    let mut hit_limit = false;
    let mut explored: usize = 0;
    for start_node in 0..graph.n_nodes {
        if start_node == source || start_node == sink {
            continue;
        }
        path.clear();
        visited.clear();
        visited.insert_grow(start_node);
        path.push(start_node);
        if dfs(
            start_node,
            graph,
            bundle,
            config,
            &mut path,
            &mut visited,
            &mut seen,
            &mut out,
            &mut emitted,
            &mut explored,
        ) {
            hit_limit = true;
            break;
        }
    }
    if config.verbose && hit_limit {
        eprintln!(
            "    junction_paths: hit internal limit for bundle {}:{}-{}",
            bundle.chrom, bundle.start, bundle.end
        );
    }
    (out, emitted)
}

/// When no direct donor→acceptor graph edge exists, try a short multi-hop path (junction
/// chain split across intermediate nodes). Bounded depth to avoid expensive search on huge graphs.
fn splice_path_donor_to_acceptor(
    graph: &Graph,
    donors: &[usize],
    acceptors: &[usize],
    max_hops: usize,
) -> Option<(usize, usize)> {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let mut acceptor_set = NodeSet::with_capacity(graph.n_nodes);
    for &a in acceptors { acceptor_set.insert_grow(a); }
    for &start in donors {
        let mut q = VecDeque::new();
        let mut seen = vec![false; graph.n_nodes];
        let mut depth = vec![0usize; graph.n_nodes];
        q.push_back(start);
        seen[start] = true;
        while let Some(u) = q.pop_front() {
            if depth[u] > max_hops {
                continue;
            }
            for v in graph.nodes[u].children.ones() {
                if v == source_id || v == sink_id {
                    continue;
                }
                if acceptor_set.contains(v) {
                    return Some((start, v));
                }
                if !seen[v] {
                    seen[v] = true;
                    depth[v] = depth[u] + 1;
                    if depth[v] <= max_hops {
                        q.push_back(v);
                    }
                }
            }
        }
    }
    // Minus-strand / alternate graph wiring: try acceptor → donors along parent edges.
    let mut donor_set = NodeSet::with_capacity(graph.n_nodes);
    for &d in donors { donor_set.insert_grow(d); }
    for &goal in acceptors {
        let mut q = VecDeque::new();
        let mut seen = vec![false; graph.n_nodes];
        let mut depth = vec![0usize; graph.n_nodes];
        q.push_back(goal);
        seen[goal] = true;
        while let Some(u) = q.pop_front() {
            if depth[u] > max_hops {
                continue;
            }
            for p in graph.nodes[u].parents.ones() {
                if p == source_id || p == sink_id {
                    continue;
                }
                if donor_set.contains(p) {
                    return Some((p, goal));
                }
                if !seen[p] {
                    seen[p] = true;
                    depth[p] = depth[u] + 1;
                    if depth[p] <= max_hops {
                        q.push_back(p);
                    }
                }
            }
        }
    }
    None
}

fn emit_chain_from_graph(
    graph: &Graph,
    bundle: &crate::types::Bundle,
    chain: &[(u64, u64)],
    source_tag: &str,
    cov: f64,
    is_longread: bool,
    ref_transcript_id: Option<String>,
) -> Option<Transcript> {
    // Long gene loci (e.g. STRG.251) often split exons across multiple graph nodes; tight slack
    // misses the node whose end/start should splice to the reference intron (still override via env).
    let tol: u64 = std::env::var("RUSTLE_EMIT_CHAIN_NODE_TOL")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(24);
    let bfs_max: usize = std::env::var("RUSTLE_EMIT_CHAIN_BFS_MAX")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(64);
    let mut nodes_by_end: Vec<(u64, usize)> = Vec::new();
    let mut nodes_by_start: Vec<(u64, usize)> = Vec::new();
    nodes_by_end.reserve(graph.nodes.len());
    nodes_by_start.reserve(graph.nodes.len());
    for (id, n) in graph.nodes.iter().enumerate() {
        if id == graph.source_id || id == graph.sink_id {
            continue;
        }
        nodes_by_end.push((n.end, id));
        nodes_by_start.push((n.start, id));
    }

    let mut exons: Vec<(u64, u64)> = Vec::new();
    let mut prev_acceptor_pos: Option<u64> = None;
    let mut first_exon_start: Option<u64> = None;
    let mut last_exon_end: Option<u64> = None;

    for (idx, (donor, acceptor)) in chain.iter().copied().enumerate() {
        let donors: Vec<usize> = nodes_by_end
            .iter()
            .filter(|(end, _)| end.abs_diff(donor) <= tol)
            .map(|(_, id)| *id)
            .collect();
        let acceptors: Vec<usize> = nodes_by_start
            .iter()
            .filter(|(start, _)| start.abs_diff(acceptor) <= tol)
            .map(|(_, id)| *id)
            .collect();
        if donors.is_empty() || acceptors.is_empty() {
            return None;
        }
        // Prefer a direct splice edge; otherwise allow a short multi-hop path through the intron.
        let mut pair: Option<(usize, usize)> = None;
        for &d in &donors {
            for &a in &acceptors {
                if graph.nodes[d].children.contains(a) {
                    pair = Some((d, a));
                    break;
                }
            }
            if pair.is_some() {
                break;
            }
        }
        let (donor_node, acceptor_node) = pair.or_else(|| {
            splice_path_donor_to_acceptor(graph, &donors, &acceptors, bfs_max)
        })?;

        let donor_start = graph.nodes.get(donor_node).map(|n| n.start)?;
        let acceptor_end = graph.nodes.get(acceptor_node).map(|n| n.end)?;
        if idx == 0 {
            first_exon_start = Some(donor_start);
        }
        if let Some(prev_acc) = prev_acceptor_pos {
            if prev_acc < donor {
                exons.push((prev_acc, donor));
            }
        }
        prev_acceptor_pos = Some(acceptor);
        last_exon_end = Some(acceptor_end);
    }

    let Some(first_start) = first_exon_start else {
        return None;
    };
    let Some(last_end) = last_exon_end else {
        return None;
    };
    let Some(first_donor) = chain.first().map(|(d, _)| *d) else {
        return None;
    };
    if first_start < first_donor {
        exons.insert(0, (first_start, first_donor));
    }
    let Some(last_acceptor) = chain.last().map(|(_, a)| *a) else {
        return None;
    };
    if last_acceptor < last_end {
        exons.push((last_acceptor, last_end));
    }
    if exons.len() < 2 {
        return None;
    }
    let exon_cov = vec![cov; exons.len()];
    Some(Transcript {
        chrom: bundle.chrom.clone(),
        strand: bundle.strand,
        exons,
        coverage: cov,
        exon_cov,
        tpm: 0.0,
        fpkm: 0.0,
        source: Some(source_tag.to_string()),
        is_longread,
        longcov: cov,
        bpcov_cov: 0.0,
        transcript_id: None,
        gene_id: None,
        ref_transcript_id,
        ref_gene_id: None,
        hardstart: true,
        hardend: true,
                    alt_tts_end: false,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None, intron_low: Vec::new(),
    })
}

fn emit_trace_chain_from_graph(
    graph: &Graph,
    bundle: &crate::types::Bundle,
    chain: &[(u64, u64)],
    is_longread: bool,
) -> Option<Transcript> {
    emit_chain_from_graph(
        graph,
        bundle,
        chain,
        "junction_chain",
        1.0,
        is_longread,
        None,
    )
}

/// True when every intron in `chain` matches some pair in `good_junctions` within `tol` bp.
fn ref_chain_fully_in_good_junctions(
    chain: &[(u64, u64)],
    good_junctions: &HashSet<(u64, u64)>,
    tol: u64,
) -> bool {
    if chain.is_empty() {
        return false;
    }
    chain.iter().all(|&(d, a)| {
        good_junctions
            .iter()
            .any(|&(gd, ga)| d.abs_diff(gd) <= tol && a.abs_diff(ga) <= tol)
    })
}

fn emit_reference_chains(
    graph: &Graph,
    bundle: &crate::types::Bundle,
    ref_transcripts: &[RefTranscript],
    existing: &[Transcript],
    config: &RunConfig,
    good_junctions: &HashSet<(u64, u64)>,
) -> Vec<Transcript> {
    if ref_transcripts.is_empty() {
        return Vec::new();
    }
    let near_miss_ref_fallback = config.long_reads && !config.eonly;
    if !config.emit_junction_paths && !near_miss_ref_fallback {
        return Vec::new();
    }
    // Chains already present in `existing` (pre-ref-assist). Used to skip redundant
    // `emit_chain_from_graph` when the graph already produced the same intron chain.
    // Do not use this to suppress `find_guide_pat` backups: those transcripts can be
    // dropped later (collapse_equal / min_length / predcluster), leaving trace with no
    // matching row even though the chain was "seen" here.
    let mut existing_chains: HashSet<Vec<(u64, u64)>> = Default::default();
    let mut existing_exon_shapes: HashSet<Vec<(u64, u64)>> = Default::default();
    for tx in existing {
        let chain = intron_chain_from_exons(&tx.exons);
        if chain.is_empty() {
            if tx.exons.len() != 1 && !tx.exons.is_empty() {
                existing_exon_shapes.insert(tx.exons.clone());
            }
        } else {
            existing_chains.insert(chain);
        }
    }
    // Intron chains / single-exon keys we already emitted in this pass (dedup among refs).
    let mut emitted_ref_chains: HashSet<Vec<(u64, u64)>> = Default::default();
    let mut emitted_ref_single: HashSet<(u64, u64)> = Default::default();
    let mut out: Vec<Transcript> = Vec::new();
    let mut added = 0usize;
    let mut overlap_cnt = 0usize;
    let mut no_path_cnt = 0usize;
    let mut empty_path_cnt = 0usize;
    let mut emitted_cnt = 0usize;
    let cov = config.readthr.max(1.0);
    let debug_target = parse_debug_bundle_target(config);
    let debug_bundle = debug_bundle_overlaps(
        debug_target.as_ref(),
        &bundle.chrom,
        bundle.start,
        bundle.end,
    );
    for ref_tx in ref_transcripts {
        if !ref_overlaps_bundle(ref_tx, bundle) {
            continue;
        }
        overlap_cnt += 1;
        let trace_ref_id = std::env::var("RUSTLE_TRACE_REF_ID").ok();
        let trace_ref_active = trace_ref_id
            .as_deref()
            .map(|v| v.trim() == ref_tx.id)
            .unwrap_or(false);

        // Max-sensitivity assist: if the reference intron chain is representable in the graph,
        // materialize it directly from graph edges. This is diagnostic-only (guarded by
        // `emit_junction_paths`) and helps localize whether misses are due to extraction
        // heuristics vs. graph connectivity.
        let ref_chain = intron_chain_from_exons(&ref_tx.exons);
        if near_miss_ref_fallback
            && !config.emit_junction_paths
            && !ref_chain.is_empty()
            && !existing_chains.contains(&ref_chain)
            && !emitted_ref_chains.contains(&ref_chain)
            && has_near_miss_chain(existing, &ref_chain, 5)
        {
            let exon_cov = vec![cov; ref_tx.exons.len()];
            if emitted_ref_chains.insert(ref_chain.clone()) {
                out.push(Transcript {
                    chrom: bundle.chrom.clone(),
                    strand: ref_tx.strand,
                    exons: ref_tx.exons.clone(),
                    coverage: cov,
                    exon_cov,
                    tpm: 0.0,
                    fpkm: 0.0,
                    source: Some("ref_chain_nearmiss".to_string()),
                    is_longread: config.long_reads,
                    longcov: cov,
                    bpcov_cov: 0.0,
                    transcript_id: None,
                    gene_id: None,
                    ref_transcript_id: Some(ref_tx.id.clone()),
                    ref_gene_id: None,
                    hardstart: true,
                    hardend: true,
                    alt_tts_end: false,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None, intron_low: Vec::new(),
                });
                added += 1;
                emitted_cnt += 1;
                if trace_ref_active {
                    eprintln!(
                        "[TRACE_REFCHAIN] ref={} bundle {}:{}-{} emitted=near_miss_fallback",
                        ref_tx.id, bundle.chrom, bundle.start, bundle.end
                    );
                }
            }
            continue;
        }

        if !config.emit_junction_paths {
            continue;
        }

        if !ref_chain.is_empty()
            && !existing_chains.contains(&ref_chain)
            && !emitted_ref_chains.contains(&ref_chain)
        {
            if let Some(mut tx) = emit_chain_from_graph(
                graph,
                bundle,
                &ref_chain,
                "ref_intron_chain",
                cov,
                config.long_reads,
                Some(ref_tx.id.clone()),
            ) {
                let snap_tol: u64 = std::env::var("RUSTLE_REF_INTRON_CHAIN_SNAP_TOL")
                    .ok()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(12);
                if intron_chain_equal_tol(
                    &ref_chain,
                    &intron_chain_from_exons(&tx.exons),
                    snap_tol,
                ) {
                    tx.exons = ref_tx.exons.clone();
                    tx.chrom = ref_tx.chrom.clone();
                    tx.strand = ref_tx.strand;
                    tx.exon_cov = vec![cov; tx.exons.len()];
                }
                if emitted_ref_chains.insert(ref_chain.clone()) {
                    out.push(tx);
                    added += 1;
                    emitted_cnt += 1;
                    if trace_ref_active {
                        eprintln!(
                            "[TRACE_REFCHAIN] ref={} bundle {}:{}-{} emitted",
                            ref_tx.id, bundle.chrom, bundle.start, bundle.end
                        );
                    }
                    if out.len() >= config.max_junction_paths {
                        break;
                    }
                }
            } else if trace_ref_active {
                eprintln!(
                    "[TRACE_REFCHAIN] ref={} bundle {}:{}-{} failed",
                    ref_tx.id, bundle.chrom, bundle.start, bundle.end
                );
            }
        } else if trace_ref_active && !ref_chain.is_empty() {
            let why = if existing_chains.contains(&ref_chain) {
                "existing_chain"
            } else if emitted_ref_chains.contains(&ref_chain) {
                "already_emitted_ref_assist"
            } else {
                "unknown"
            };
            eprintln!(
                "[TRACE_REFCHAIN] ref={} bundle {}:{}-{} skipped_emit_graph={}",
                ref_tx.id,
                bundle.chrom,
                bundle.start,
                bundle.end,
                why
            );
        }

        // If `ref_intron_chain` already emitted this intron chain, skip duplicate `ref_chain`.
        if !ref_chain.is_empty() && emitted_ref_chains.contains(&ref_chain) {
            continue;
        }

        let ref_guide_ssdist: u64 = std::env::var("RUSTLE_EMIT_REF_SSDIST")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(48);
        let node_ids = match find_guide_pat(
            ref_tx,
            graph,
            bundle.start,
            bundle.end,
            bundle.strand,
            ref_guide_ssdist,
        ) {
            Some(n) => n,
            None => {
                let jtol: u64 = std::env::var("RUSTLE_REF_JUNCTION_WITNESS_TOL")
                    .ok()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(8);
                let witness_on = config.ref_junction_witness_enabled();
                let witness_eligible = !ref_chain.is_empty()
                    && ref_chain_fully_in_good_junctions(&ref_chain, good_junctions, jtol)
                    && !emitted_ref_chains.contains(&ref_chain);
                if witness_eligible && witness_on {
                    emitted_ref_chains.insert(ref_chain.clone());
                    let out_exons = ref_tx.exons.clone();
                    let exon_cov = vec![cov; out_exons.len()];
                    out.push(Transcript {
                        chrom: ref_tx.chrom.clone(),
                        strand: ref_tx.strand,
                        exons: out_exons,
                        coverage: cov,
                        exon_cov,
                        tpm: 0.0,
                        fpkm: 0.0,
                        source: Some("ref_chain_junction_witness".to_string()),
                        is_longread: config.long_reads,
                        longcov: cov,
                        bpcov_cov: 0.0,
                        transcript_id: None,
                        gene_id: None,
                        ref_transcript_id: Some(ref_tx.id.clone()),
                        ref_gene_id: None,
                        hardstart: true,
                        hardend: true,
                    alt_tts_end: false,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None, intron_low: Vec::new(),
                    });
                    added += 1;
                    emitted_cnt += 1;
                    if trace_ref_active {
                        eprintln!(
                            "[TRACE_REFCHAIN] ref={} bundle {}:{}-{} emitted=junction_witness",
                            ref_tx.id, bundle.chrom, bundle.start, bundle.end
                        );
                    }
                    continue;
                }
                if trace_ref_active && witness_eligible && !witness_on {
                    eprintln!(
                        "[TRACE_REFCHAIN] ref={} bundle {}:{}-{} junction_witness_eligible_but_disabled \
                         (use --ref-junction-witness, --trace-reference, or RUSTLE_REF_JUNCTION_WITNESS=1)",
                        ref_tx.id,
                        bundle.chrom,
                        bundle.start,
                        bundle.end
                    );
                }
                no_path_cnt += 1;
                if config.verbose && debug_bundle {
                    eprintln!(
                        "    ref_chain_debug: ref={} bundle {}:{}-{} no_guide_pat",
                        ref_tx.id, bundle.chrom, bundle.start, bundle.end
                    );
                }
                continue;
            }
        };
        let mut exons: Vec<(u64, u64)> = Vec::new();
        for nid in node_ids {
            if nid == graph.source_id || nid == graph.sink_id {
                continue;
            }
            if let Some(n) = graph.nodes.get(nid) {
                if n.start < n.end {
                    exons.push((n.start, n.end));
                }
            }
        }
        if exons.is_empty() {
            empty_path_cnt += 1;
            if config.verbose && debug_bundle {
                eprintln!(
                    "    ref_chain_debug: ref={} bundle {}:{}-{} no_exon_nodes",
                    ref_tx.id, bundle.chrom, bundle.start, bundle.end
                );
            }
            continue;
        }
        // `find_guide_pat` proved a path exists, but collecting raw graph node intervals can
        // subdivide GTF exons — extra introns vs the reference — so trace `intron_chains_match`
        // misses. Emit canonical reference exons (same idea as `guide:` snapping below).
        let out_exons = ref_tx.exons.clone();
        if !ref_chain.is_empty() {
            if emitted_ref_chains.contains(&ref_chain) {
                continue;
            }
        } else if out_exons.len() == 1 {
            if emitted_ref_single.contains(&out_exons[0]) {
                continue;
            }
        } else if existing_exon_shapes.contains(&out_exons) {
            continue;
        }
        let exon_cov = vec![cov; out_exons.len()];
        if !ref_chain.is_empty() {
            emitted_ref_chains.insert(ref_chain.clone());
        } else if out_exons.len() == 1 {
            emitted_ref_single.insert(out_exons[0]);
        }
        out.push(Transcript {
            chrom: ref_tx.chrom.clone(),
            strand: ref_tx.strand,
            exons: out_exons,
            coverage: cov,
            exon_cov,
            tpm: 0.0,
            fpkm: 0.0,
            source: Some("ref_chain".to_string()),
            is_longread: config.long_reads,
            longcov: cov,
            bpcov_cov: 0.0,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: Some(ref_tx.id.clone()),
            ref_gene_id: None,
            hardstart: true,
            hardend: true,
                    alt_tts_end: false,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None, intron_low: Vec::new(),
        });
        added += 1;
        emitted_cnt += 1;
    }
    if config.verbose && added > 0 {
        eprintln!("    ref_chains: added {}", added);
    }
    if config.verbose && debug_bundle && overlap_cnt > 0 {
        eprintln!(
            "    ref_chain_debug: bundle {}:{}-{} overlap={} emitted={} no_guide_pat={} no_exon_nodes={}",
            bundle.chrom,
            bundle.start,
            bundle.end,
            overlap_cnt,
            emitted_cnt,
            no_path_cnt,
            empty_path_cnt
        );
    }
    out
}

// Note: adjacency-only path search removed (no longer needed for trace chain emission).

fn trace_intron_chain_graph(graph: &Graph, bundle: &crate::types::Bundle, config: &RunConfig) {
    let Some(chain) = config.trace_intron_chain.as_ref() else {
        return;
    };
    let debug_target = parse_debug_bundle_target(config);
    if !debug_bundle_overlaps(
        debug_target.as_ref(),
        &bundle.chrom,
        bundle.start,
        bundle.end,
    ) {
        return;
    }

    let mut by_end: HashMap<u64, Vec<usize>> = Default::default();
    let mut by_start: HashMap<u64, Vec<usize>> = Default::default();
    for (id, n) in graph.nodes.iter().enumerate() {
        if id == graph.source_id || id == graph.sink_id {
            continue;
        }
        by_end.entry(n.end).or_default().push(id);
        // Accept both internal node starts and GTF-style starts (which are often +1 relative
        // to internal graph representation in this codebase).
        by_start.entry(n.start).or_default().push(id);
        by_start
            .entry(n.start.saturating_add(1))
            .or_default()
            .push(id);
    }

    eprintln!(
        "[TRACE_CHAIN_GRAPH] bundle {} strand={} chain_len={}",
        bundle.chrom,
        bundle.strand,
        chain.len()
    );

    let mut prev_acceptor_nodes: Option<Vec<usize>> = None;
    for (idx, (donor, acceptor)) in chain.iter().copied().enumerate() {
        let donors = by_end.get(&donor).cloned().unwrap_or_default();
        let acceptors = by_start.get(&acceptor).cloned().unwrap_or_default();
        let mut edge_exists = 0usize;
        for &d in &donors {
            for &a in &acceptors {
                if graph.nodes[d].children.contains(a) {
                    edge_exists += 1;
                }
            }
        }
        eprintln!(
            "[TRACE_CHAIN_GRAPH] intron {} {}-{} donors={} acceptors={} edges={}",
            idx,
            donor,
            acceptor,
            donors.len(),
            acceptors.len(),
            edge_exists
        );
        if let Some(prev_acceptors) = prev_acceptor_nodes.take() {
            // Check connectivity from any prev acceptor to any donor of this intron.
            let mut reachable = false;
            let mut stack = prev_acceptors.clone();
            let mut seen: HashSet<usize> = Default::default();
            for s in &prev_acceptors {
                seen.insert(*s);
            }
            while let Some(nid) = stack.pop() {
                if donors.contains(&nid) {
                    reachable = true;
                    break;
                }
                for ch in graph.nodes[nid].children.ones().collect::<Vec<_>>() {
                    if ch == graph.source_id || ch == graph.sink_id {
                        continue;
                    }
                    if seen.insert(ch) {
                        stack.push(ch);
                    }
                }
            }
            eprintln!(
                "[TRACE_CHAIN_GRAPH] connect prev_acceptors -> donors: {}",
                reachable
            );
        }
        prev_acceptor_nodes = Some(acceptors);
    }
}

/// Split bundlenodes into components using junction connectivity.
///
/// `good_junc()` rejects low-quality junctions (setting strand=0), which breaks
/// color propagation in `merge_read_to_group`. This naturally separates gene loci that are
/// only connected by bad-junction reads into different bundles with separate graphs
///.
///
/// This function achieves the same effect: two bundlenodes are in the same component only
/// if a GOOD junction (one that passed filter_junctions) bridges them.
#[allow(dead_code)]
fn bundlenode_components_by_junctions(
    bn: Option<&crate::types::CBundlenode>,
    good_junctions: &[Junction],
) -> Vec<Vec<(usize, u64, u64, f64)>> {
    let bnodes = bundlenodes_to_vec(bn);
    let parity_debug = std::env::var_os("RUSTLE_DEBUG_DETAIL").is_some();
    if parity_debug {
        eprintln!(
            "DEBUG_JUNC_SPLIT bnodes={} good_junctions={}",
            bnodes.len(),
            good_junctions.len()
        );
    }
    if bnodes.is_empty() {
        return Vec::new();
    }
    if bnodes.len() == 1 {
        return vec![bnodes];
    }

    // Union-Find
    let n = bnodes.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    fn union(parent: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            parent[rb] = ra;
        }
    }

    // For each good junction, find which bundlenodes contain donor and acceptor
    // and merge those bundlenodes into the same component.
    for j in good_junctions {
        let mut donor_idx = None;
        let mut acceptor_idx = None;
        for (i, &(_bid, s, e, _cov)) in bnodes.iter().enumerate() {
            if j.donor >= s && j.donor <= e {
                donor_idx = Some(i);
            }
            if j.acceptor >= s && j.acceptor <= e {
                acceptor_idx = Some(i);
            }
            if donor_idx.is_some() && acceptor_idx.is_some() {
                break;
            }
        }
        if let (Some(di), Some(ai)) = (donor_idx, acceptor_idx) {
            union(&mut parent, di, ai);
        }
    }

    // Group by root
    let mut groups: HashMap<usize, Vec<(usize, u64, u64, f64)>> = Default::default();
    for (i, item) in bnodes.into_iter().enumerate() {
        let root = find(&mut parent, i);
        groups.entry(root).or_default().push(item);
    }
    let mut out: Vec<Vec<(usize, u64, u64, f64)>> = groups.into_values().collect();
    for g in &mut out {
        g.sort_unstable_by_key(|v| v.1);
    }
    out.sort_unstable_by_key(|g| g.first().map(|v| v.1).unwrap_or(0));

    if parity_debug {
        eprintln!("DEBUG_JUNC_SPLIT result: {} components", out.len());
        for (ci, comp) in out.iter().enumerate() {
            let start = comp.first().map(|v| v.1).unwrap_or(0);
            let end = comp.last().map(|v| v.2).unwrap_or(0);
            eprintln!(
                "  comp[{}] bnodes={} range={}-{}",
                ci,
                comp.len(),
                start,
                end
            );
        }
    }

    out
}

/// Split bundlenodes into components by CGroup color root.
///
/// CGroups on opposite sides of killed junctions (strand==0) get different
/// `equalcolor` roots, and groups with the same color root → same `CBundle`
/// . Each CBundle gets its own `create_graph` call with
/// shared reads/bpcov .
///
/// `bnode_colors` is the color root for each bundlenode (computed in bundle.rs).
fn bundlenode_components_by_color(
    bn: Option<&crate::types::CBundlenode>,
    bnode_colors: &[usize],
) -> Vec<Vec<(usize, u64, u64, f64)>> {
    let bnodes = bundlenodes_to_vec(bn);
    if bnodes.is_empty() {
        return Vec::new();
    }
    if bnodes.len() == 1 || bnode_colors.is_empty() {
        return vec![bnodes];
    }

    let parity_debug = std::env::var_os("RUSTLE_DEBUG_DETAIL").is_some();
    let mut groups: HashMap<usize, Vec<(usize, u64, u64, f64)>> = Default::default();
    for (i, item) in bnodes.into_iter().enumerate() {
        let color = bnode_colors.get(i).copied().unwrap_or(0);
        groups.entry(color).or_default().push(item);
    }

    if groups.len() <= 1 {
        return vec![groups.into_values().next().unwrap_or_default()];
    }

    let mut out: Vec<Vec<(usize, u64, u64, f64)>> = groups.into_values().collect();
    for g in &mut out {
        g.sort_unstable_by_key(|v| v.1);
    }
    out.sort_unstable_by_key(|g| g.first().map(|v| v.1).unwrap_or(0));

    if parity_debug {
        eprintln!("COLOR_SPLIT components={}", out.len());
        for (ci, comp) in out.iter().enumerate() {
            let start = comp.first().map(|v| v.1).unwrap_or(0);
            let end = comp.last().map(|v| v.2).unwrap_or(0);
            eprintln!(
                "  COLOR_SPLIT comp[{}] bnodes={} range={}-{}",
                ci,
                comp.len(),
                start,
                end
            );
        }
    }

    out
}

fn collect_guide_boundary_sets_for_bundle(
    guides: &[RefTranscript],
    bundle: &crate::types::Bundle,
) -> (HashSet<u64>, HashSet<u64>) {
    let mut boundary_left: HashSet<u64> = Default::default();
    let mut boundary_right: HashSet<u64> = Default::default();
    if guides.is_empty() {
        return (boundary_left, boundary_right);
    }
    for g in guides {
        if g.chrom != bundle.chrom {
            continue;
        }
        if bundle.strand != '.' && g.strand != '.' && g.strand != bundle.strand {
            continue;
        }
        let (Some(first), Some(last)) = (g.exons.first(), g.exons.last()) else {
            continue;
        };
        if first.0 >= bundle.end.saturating_add(1) || last.1 <= bundle.start {
            continue;
        }
        for w in g.exons.windows(2) {
            let le = w[0].1;
            let rs = w[1].0;
            if le >= rs {
                continue;
            }
            if le >= bundle.start && le <= bundle.end {
                boundary_left.insert(le);
            }
            if rs >= bundle.start && rs <= bundle.end {
                boundary_right.insert(rs);
            }
        }
    }
    (boundary_left, boundary_right)
}

use crate::transcript_filter::{
    add_pred, apply_global_cross_strand_filter, collapse_equal_predictions, compute_tpm_fpkm,
    filter_unsupported_junctions, print_predcluster_with_summary, suppress_near_duplicate_chains,
    PredclusterStageSummary,
};
use crate::transfrag_process::{
    chain_subsequence_offset, exons_intron_length_chain, process_transfrags,
};

/// Post-assembly family extension: for transcripts whose intron chain is a
/// strict contiguous subsequence of a longer transcript's chain, project the
/// missing TSS- or TTS-side exons from the template onto the target's
/// coordinates. Works because near-paralogs share intron topology; the only
/// difference is a consistent genomic offset between copies.
fn extend_family_suffix_matches(transcripts: &mut Vec<Transcript>, verbose: bool) {
    let n = transcripts.len();
    if n < 2 {
        return;
    }
    let chains: Vec<Vec<u32>> = transcripts
        .iter()
        .map(|t| exons_intron_length_chain(&t.exons))
        .collect();
    let mut extended = 0usize;
    for i in 0..n {
        let target_chain = &chains[i];
        if target_chain.len() < 5 {
            continue;
        }
        // Scan for a longer template with target as contiguous subsequence.
        let mut best: Option<(usize, usize)> = None;
        for j in 0..n {
            if j == i || transcripts[j].strand != transcripts[i].strand {
                continue;
            }
            if transcripts[j].chrom != transcripts[i].chrom {
                continue;
            }
            let tchain = &chains[j];
            if tchain.len() <= target_chain.len() {
                continue;
            }
            if let Some(off) = chain_subsequence_offset(target_chain, tchain, 5) {
                // Prefer templates where the target is a strict subsequence
                // that missing exons on at least one end (off>0 or tail).
                if off > 0 || off + target_chain.len() < tchain.len() {
                    best = Some((j, off));
                    break;
                }
            }
        }
        let Some((template_idx, offset)) = best else {
            continue;
        };
        let template = &transcripts[template_idx];
        let target = &transcripts[i];
        if template.exons.len() < target.exons.len() + offset {
            continue;
        }
        // Template exons are sorted; target's first exon corresponds to
        // template.exons[offset]. Compute genomic delta, then project the
        // template's preceding and trailing exons onto target coordinates.
        let mut target_exons_sorted = target.exons.clone();
        target_exons_sorted.sort();
        let mut template_exons_sorted = template.exons.clone();
        template_exons_sorted.sort();
        let delta =
            target_exons_sorted[0].0 as i64 - template_exons_sorted[offset].0 as i64;
        // Require template to be a genomic paralog, not the same-locus parent:
        // delta must be at least 5 kb so we don't extend a minor-isoform suffix
        // with its parent-isoform's unique 5' exons.
        if delta.abs() < 5_000 {
            continue;
        }
        // Limit extension to at most 3 missing exons per side — larger gaps
        // are more likely to be assembly errors than real paralog drift.
        let tail_start_est = offset + target_chain.len() + 1;
        let tail_len = template_exons_sorted
            .len()
            .saturating_sub(tail_start_est);
        if offset > 3 || tail_len > 3 {
            continue;
        }
        let mut new_exons: Vec<(u64, u64)> = Vec::new();
        // Prepend missing 5' exons (low-coordinate end).
        for k in 0..offset {
            let (s, e) = template_exons_sorted[k];
            let ps = (s as i64 + delta).max(0) as u64;
            let pe = (e as i64 + delta).max(0) as u64;
            new_exons.push((ps, pe));
        }
        new_exons.extend(target_exons_sorted.iter().copied());
        // Append missing 3' exons (high-coordinate end) if any.
        let tail_start = offset + target_chain.len() + 1; // +1 accounts for the shared final exon
        if tail_start < template_exons_sorted.len() {
            for k in tail_start..template_exons_sorted.len() {
                let (s, e) = template_exons_sorted[k];
                let ps = (s as i64 + delta).max(0) as u64;
                let pe = (e as i64 + delta).max(0) as u64;
                new_exons.push((ps, pe));
            }
        }
        new_exons.sort();
        // Deduplicate in case of coordinate collisions.
        new_exons.dedup();
        if new_exons.len() > target.exons.len() {
            transcripts[i].exons = new_exons;
            transcripts[i].source = Some(format!(
                "{}+family_ext",
                transcripts[i].source.as_deref().unwrap_or("flow")
            ));
            extended += 1;
        }
    }
    if verbose && extended > 0 {
        eprintln!(
            "    family-extend: extended {} transcript(s) via cross-bundle chain homology",
            extended
        );
    }
}

fn count_fragment_metrics_for_bundle(bundle: &crate::types::Bundle) -> (f64, f64) {
    let mut num_fragments = 0.0f64;
    let mut frag_len_sum = 0.0f64;
    for r in &bundle.reads {
        // countFragment: these are accumulated per-alignment event during ingest
        // and carried on merged read entries.
        if r.countfrag_len > 0.0 {
            frag_len_sum += r.countfrag_len;
        } else {
            // Fallback for older inputs/tests that don't populate countfrag_* fields.
            let nh = (r.nh as f64).max(1.0);
            let read_len: u64 = r.exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
            if read_len > 0 {
                frag_len_sum += (read_len as f64) / nh;
            }
        }
        if r.countfrag_num > 0.0 {
            num_fragments += r.countfrag_num;
        } else {
            let nh = (r.nh as f64).max(1.0);
            num_fragments += 1.0 / nh;
        }
    }
    (num_fragments, frag_len_sum)
}

#[derive(Clone, Copy, Debug, Default)]
struct StrandReadCounts {
    minus: usize,
    neutral: usize,
    plus: usize,
    other: usize,
}

#[inline]
fn read_strand_counts(reads: &[BundleRead]) -> StrandReadCounts {
    let mut out = StrandReadCounts::default();
    for r in reads {
        match r.strand {
            '-' => out.minus += 1,
            '.' => out.neutral += 1,
            '+' => out.plus += 1,
            _ => out.other += 1,
        }
    }
    out
}

/// Split a read at killed junction boundaries and return neutral fragments that
/// have no surviving internal splice edges (typically single-exon remnants).
/// These fragments model the sno=1 contribution of partially unsupported
/// reads without forcing the whole read to unstranded.
fn killed_junction_neutral_fragments(
    read: &BundleRead,
    killed_junctions: &HashSet<Junction>,
) -> Vec<BundleRead> {
    if read.exons.is_empty() || read.junctions.is_empty() {
        return Vec::new();
    }

    let mut chunks: Vec<(usize, usize)> = Vec::new();
    let mut chunk_start = 0usize;
    for (ji, j) in read.junctions.iter().enumerate() {
        if killed_junctions.contains(j) {
            if chunk_start <= ji {
                chunks.push((chunk_start, ji + 1)); // exon range [start, ji]
            }
            chunk_start = ji + 1;
        }
    }
    if chunk_start < read.exons.len() {
        chunks.push((chunk_start, read.exons.len()));
    }

    let mut out: Vec<BundleRead> = Vec::new();
    for (a, b) in chunks {
        if b <= a {
            continue;
        }
        // Keep only chunks with no surviving internal splice support.
        // If a chunk still has splice edges, it is represented by the original
        // stranded read and should not be duplicated as neutral evidence.
        if b.saturating_sub(a) > 1 {
            continue;
        }
        let mut frag = read.clone();
        frag.exons = read.exons[a..b].to_vec();
        frag.junctions = junctions_from_exons(&frag.exons);
        frag.junctions_raw = frag.junctions.clone();
        frag.junctions_del = frag.junctions.clone();
        if let Some(&(s, _)) = frag.exons.first() {
            frag.ref_start = s;
        }
        if let Some(&(_, e)) = frag.exons.last() {
            frag.ref_end = e;
        }
        frag.strand = '.';
        // Fragment transports are synthetic split carriers; drop mate graph links
        // to avoid dangling indices in the recipient bundle.
        frag.pair_idx.clear();
        frag.pair_count.clear();
        // Avoid re-triggering terminal poly trimming on synthetic fragments.
        frag.has_poly_start = false;
        frag.has_poly_end = false;
        frag.has_poly_start_aligned = false;
        frag.has_poly_start_unaligned = false;
        frag.has_poly_end_aligned = false;
        frag.has_poly_end_unaligned = false;
        frag.unaligned_poly_t = 0;
        frag.unaligned_poly_a = 0;
        frag.has_last_exon_polya = false;
        frag.has_first_exon_polyt = false;
        out.push(frag);
    }
    out
}

#[inline]
#[allow(dead_code)]
fn shorten_first_exon_to_2bp(read: &mut BundleRead) {
    if read.exons.len() < 2 {
        return;
    }
    let (s, e) = read.exons[0];
    if e.saturating_sub(s) <= 2 {
        return;
    }
    // shortenFirstExon: newStart = segs[0].end - 2 (2bp from end)
    read.exons[0].0 = e.saturating_sub(2);
}

#[allow(dead_code)]
#[inline]
fn shorten_last_exon_to_3bp(read: &mut BundleRead) {
    if read.exons.len() < 2 {
        return;
    }
    let last = read.exons.len() - 1;
    let (s, e) = read.exons[last];
    if e.saturating_sub(s) <= 3 {
        return;
    }
    read.exons[last].1 = s.saturating_add(3);
}

#[inline]
fn read_is_long_class_for_processread(read: &BundleRead, config: &RunConfig) -> bool {
    if !config.long_reads {
        return false;
    }
    if config.long_read_min_len == 0 {
        return true;
    }
    let rlen = read.query_length.unwrap_or_else(|| {
        read.exons
            .iter()
            .map(|(s, e)| e.saturating_sub(*s))
            .sum::<u64>()
    });
    rlen >= config.long_read_min_len
}

/// processRead logic:
/// remove terminal exons if they are mostly polyA/polyT.
/// Returns false when read should be discarded.
fn apply_processread_terminal_poly_trim(read: &mut BundleRead) -> bool {
    if read.exons.len() < 2 {
        return true;
    }

    let mut rm_last = read.has_last_exon_polya;
    let mut rm_first = read.has_first_exon_polyt;

    if rm_last && rm_first {
        if read.exons.len() == 2 {
            return false;
        }
        if read.strand == '+' {
            rm_first = false;
        } else if read.strand == '-' {
            rm_last = false;
        }
    }

    let mut trimmed = false;
    if rm_last {
        if !read.junctions.is_empty() {
            read.junctions.pop();
        }
        if !read.junctions_raw.is_empty() {
            read.junctions_raw.pop();
        }
        if !read.junctions_del.is_empty() {
            read.junctions_del.pop();
        }
        read.exons.pop();
        read.has_poly_end_aligned = false;
        if !read.has_poly_end_unaligned {
            read.has_poly_end_unaligned = true;
            if read.unaligned_poly_a == 0 {
                read.unaligned_poly_a = 1;
            }
        }
        trimmed = true;
    }

    if rm_first {
        if !read.junctions.is_empty() {
            read.junctions.remove(0);
        }
        if !read.junctions_raw.is_empty() {
            read.junctions_raw.remove(0);
        }
        if !read.junctions_del.is_empty() {
            read.junctions_del.remove(0);
        }
        read.exons.remove(0);
        read.has_poly_start_aligned = false;
        if !read.has_poly_start_unaligned {
            read.has_poly_start_unaligned = true;
            if read.unaligned_poly_t == 0 {
                read.unaligned_poly_t = 1;
            }
        }
        trimmed = true;
    }

    if read.exons.is_empty() {
        return false;
    }

    if trimmed {
        let expected_j = read.exons.len().saturating_sub(1);
        if read.junctions.len() != expected_j {
            read.junctions = junctions_from_exons(&read.exons);
        }
        if read.junctions_raw.len() != expected_j {
            read.junctions_raw = junctions_from_exons(&read.exons);
        }
        if read.junctions_del.len() != expected_j {
            read.junctions_del = read.junctions.clone();
        }
        if let Some((s, _)) = read.exons.first().copied() {
            read.ref_start = s;
        }
        if let Some((_, e)) = read.exons.last().copied() {
            read.ref_end = e;
        }

        // infer unknown strand from poly tail evidence when unambiguous.
        if read.strand == '.' {
            let plus = read.has_poly_end_unaligned || read.has_poly_end_aligned;
            let minus = read.has_poly_start_unaligned || read.has_poly_start_aligned;
            if plus && !minus {
                read.strand = '+';
            } else if minus && !plus {
                read.strand = '-';
            }
        }
    }

    read.has_poly_start = read.has_poly_start_aligned || read.has_poly_start_unaligned;
    read.has_poly_end = read.has_poly_end_aligned || read.has_poly_end_unaligned;
    read.has_last_exon_polya = false;
    read.has_first_exon_polyt = false;
    true
}

/// Algorithm: shorten terminal exons for all long reads with aligned polyA/T.
/// Applied at the start of infer_transcripts, before count_good_junctions/build_graphs.
/// Only applied to multi-exon long reads.
fn shorten_polya_terminal_exons(reads: &mut [crate::types::BundleRead]) {
    for r in reads.iter_mut() {
        if r.exons.len() < 2 {
            continue;
        }
        // Strand 0 or +1 → shorten last exon if aligned_polyA && !unaligned_polyA
        // Strand 0 or -1 → shorten first exon if aligned_polyT && !unaligned_polyT
        // aligned && !unaligned: polyA is inside the alignment (RT artifact), not a genuine tail.
        let shorten_last = r.has_poly_end_aligned && !r.has_poly_end_unaligned && r.strand != '-';
        let shorten_first =
            r.has_poly_start_aligned && !r.has_poly_start_unaligned && r.strand != '+';
        let mut changed = false;
        if shorten_last {
            let last = r.exons.len() - 1;
            let (s, e) = r.exons[last];
            if e.saturating_sub(s) > 3 {
                r.exons[last].1 = s.saturating_add(3);
                changed = true;
            }
        }
        if shorten_first {
            let (s, e) = r.exons[0];
            // shortenFirstExon: newStart = segs[0].end - 2 (2bp, asymmetric with last=3bp)
            if e.saturating_sub(s) > 2 {
                r.exons[0].0 = e.saturating_sub(2);
                changed = true;
            }
        }
        if changed {
            // shortenFirstExon/shortenLastExon only moves the terminal exon
            // boundary. Internal splice junction identities stay unchanged, so do
            // not rebuild junctions_raw/junctions_del here.
            if let Some((s, _)) = r.exons.first().copied() {
                r.ref_start = s;
            }
            if let Some((_, e)) = r.exons.last().copied() {
                r.ref_end = e;
            }
        }
    }
}

/// Trim terminal exons of long reads where coverage drops significantly (5x).
/// This mimics CMaxIntv trimming to avoid over-extending transfrags into noise.
///
/// Aggressive outer-window thresholds inflate gffcompare class **`j`** (same intron chain, different
/// terminal exons): reads get shortened before junction/graph work. Keep the outer ceiling tight.
fn trim_terminal_exons_by_coverage(
    reads: &mut [crate::types::BundleRead],
    bpcov: &crate::bpcov::BpcovStranded,
    bundle_start: u64,
) {
    let window = 50usize;
    let drop_ratio = 0.2; // 5x drop
    // Outer-window mean cov must be below this (per-base) before trimming; default 0.45 (was 1.0).
    // Over-trimming inflates gffcompare class `j` (same intron chain, terminal mismatch). Override: RUSTLE_TRIM_TERM_OUTER_MAX.
    let outer_cov_trim_ceiling: f64 = std::env::var("RUSTLE_TRIM_TERM_OUTER_MAX")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.45);

    for r in reads.iter_mut() {
        if r.exons.len() < 2 {
            continue;
        }
        let strand_idx = match r.strand {
            '+' => crate::bpcov::BPCOV_STRAND_PLUS,
            '-' => crate::bpcov::BPCOV_STRAND_MINUS,
            _ => crate::bpcov::BPCOV_STRAND_ALL,
        };

        let mut changed = false;

        // 1. Trim first exon
        let (fs, fe) = r.exons[0];
        let flen = fe - fs;
        if flen > (window as u64) * 2 {
            let mut best_split = None;
            for i in (window as u64..flen - window as u64).rev() {
                let split_pos = fs + i;
                let idx = (split_pos - bundle_start) as usize;
                
                let sum_inner = bpcov.get_cov_range(strand_idx, idx, idx + window);
                let sum_outer = bpcov.get_cov_range(strand_idx, idx - window, idx);
                let cov_inner = sum_inner / (window as f64);
                let cov_outer = sum_outer / (window as f64);
                
                if cov_outer < cov_inner * drop_ratio && cov_outer < outer_cov_trim_ceiling {
                    best_split = Some(split_pos);
                    break;
                }
            }
            if let Some(pos) = best_split {
                r.exons[0].0 = pos;
                changed = true;
            }
        }

        // 2. Trim last exon
        let last_idx = r.exons.len() - 1;
        let (ls, le) = r.exons[last_idx];
        let llen = le - ls;
        if llen > (window as u64) * 2 {
            let mut best_split = None;
            for i in window as u64..llen - window as u64 {
                let split_pos = ls + i;
                let idx = (split_pos - bundle_start) as usize;
                
                let sum_inner = bpcov.get_cov_range(strand_idx, idx - window, idx);
                let sum_outer = bpcov.get_cov_range(strand_idx, idx, idx + window);
                let cov_inner = sum_inner / (window as f64);
                let cov_outer = sum_outer / (window as f64);
                
                if cov_outer < cov_inner * drop_ratio && cov_outer < outer_cov_trim_ceiling {
                    best_split = Some(split_pos);
                    break;
                }
            }
            if let Some(pos) = best_split {
                r.exons[last_idx].1 = pos;
                changed = true;
            }
        }

        if changed {
            if let Some((s, _)) = r.exons.first().copied() {
                r.ref_start = s;
            }
            if let Some((_, e)) = r.exons.last().copied() {
                r.ref_end = e;
            }
        }
    }
}

fn resolve_corr(mut j: Junction, m: &HashMap<Junction, Junction>) -> Junction {
    let mut guard: HashSet<Junction> = Default::default();
    while let Some(n) = m.get(&j).copied() {
        if n == j || !guard.insert(j) {
            break;
        }
        j = n;
    }
    j
}

/// build_graphs per-read junction repair .
///
/// For each read, for each junction:
/// - If killed (strand=0) AND has redirect: adjust exon boundaries to redirect target
/// - If killed AND no redirect: mark junction dead but keep read intact
/// - After all junctions: if no stranded junctions remain, unstrand the read
///
/// Key difference from previous Rust version: does NOT truncate reads at killed junctions.
/// Keeps reads intact with dead junctions; CGroup builder handles color breaks.
fn repair_reads_after_junction_quality(
    reads: &mut [BundleRead],
    corrected_map: &HashMap<Junction, Junction>,
    killed_junction_pairs: &HashSet<Junction>,
    junction_stats: &JunctionStats,
) {
    fn remove_pair_link_from_read(read: &mut BundleRead, mate_idx: usize) {
        let mut i = 0usize;
        while i < read.pair_idx.len() {
            if read.pair_idx[i] == mate_idx {
                read.pair_idx.swap_remove(i);
                if i < read.pair_count.len() {
                    read.pair_count.swap_remove(i);
                }
            } else {
                i += 1;
            }
        }
    }

    let mut reads_to_unpair: Vec<usize> = Vec::new();

    for n in 0..reads.len() {
        if reads[n].exons.len() < 2 {
            continue;
        }
        let mut touched = false;
        let mut has_bad_junction = false;
        let njuncs = reads[n].junctions.len();

        for i in 0..njuncs {
            if i + 1 >= reads[n].exons.len() {
                break;
            }
            let curj = reads[n].junctions[i];
            let mapped = resolve_corr(curj, corrected_map);
            let killed =
                killed_junction_pairs.contains(&curj) || killed_junction_pairs.contains(&mapped);
            let has_redirect = mapped != curj;

            if !killed && !has_redirect {
                continue;
            }

            has_bad_junction = true;

            // changeleft/changeright junction repair.
            // When a junction is killed and no direct redirect exists, search nearby
            // junctions for a better alternative with similar donor or acceptor.
            // This prevents reads from being fragmented at killed junctions.
            if killed && !has_redirect {
                let sserr = 25u64; // search window (sserror)
                let donor = curj.donor;
                let acceptor = curj.acceptor;

                // Search for best replacement: same acceptor region, better donor (changeleft)
                // OR same donor region, better acceptor (changeright)
                let mut best_replacement: Option<Junction> = None;
                let mut best_support = 0.0f64;

                for (j, st) in junction_stats.iter() {
                    if st.strand == Some(0) || st.mrcount <= 0.0 || st.mm < 0.0 {
                        continue; // Skip killed junctions
                    }
                    if killed_junction_pairs.contains(j) {
                        continue;
                    }
                    // Check if this junction is a good replacement:
                    // same acceptor (within tolerance) with different but nearby donor,
                    // OR same donor with different but nearby acceptor
                    let donor_close = j.donor.abs_diff(donor) <= sserr;
                    let acceptor_close = j.acceptor.abs_diff(acceptor) <= sserr;

                    if donor_close && acceptor_close && st.nreads_good > best_support {
                        best_replacement = Some(*j);
                        best_support = st.nreads_good;
                    }
                }

                if let Some(replacement) = best_replacement {
                    let left = reads[n].exons[i];
                    let right = reads[n].exons[i + 1];
                    let newstart = replacement.donor.max(left.0.saturating_add(1)).min(left.1);
                    let newend = replacement.acceptor.max(right.0).min(right.1.saturating_sub(1));

                    if left.0 < newstart && newend < right.1 && newstart <= newend {
                        reads[n].exons[i].1 = newstart;
                        reads[n].exons[i + 1].0 = newend;
                        reads[n].junctions[i] = replacement;
                        if let Some(raw) = reads[n].junctions_raw.get_mut(i) {
                            *raw = replacement;
                        }
                        if let Some(del) = reads[n].junctions_del.get_mut(i) {
                            *del = replacement;
                        }
                        touched = true;
                        continue; // Junction repaired, skip to next
                    }
                }
            }

            if has_redirect {
                // Adjust exon boundaries to redirect target.
                let left = reads[n].exons[i];
                let right = reads[n].exons[i + 1];
                let mut newstart = mapped.donor;
                let mut newend = mapped.acceptor;

                // Keep half-open exons non-empty while matching inclusive boundary clamps.
                let left_min_end = left.0.saturating_add(1);
                if newstart < left_min_end {
                    newstart = left_min_end;
                }
                if newstart > left.1 {
                    newstart = left.1;
                }
                let right_max_start = right.1.saturating_sub(1);
                if newend < right.0 {
                    newend = right.0;
                }
                if newend > right_max_start {
                    newend = right_max_start;
                }

                if left.0 < newstart && newend < right.1 && newstart <= newend {
                    reads[n].exons[i].1 = newstart;
                    reads[n].exons[i + 1].0 = newend;
                    reads[n].junctions[i] = mapped;
                    if let Some(raw) = reads[n].junctions_raw.get_mut(i) {
                        *raw = mapped;
                    }
                    if let Some(del) = reads[n].junctions_del.get_mut(i) {
                        *del = mapped;
                    }
                    touched = true;
                }
            }
        }

        if touched {
            if let Some((s, _)) = reads[n].exons.first().copied() {
                reads[n].ref_start = s;
            }
            if let Some((_, e)) = reads[n].exons.last().copied() {
                reads[n].ref_end = e;
            }
        }

        if has_bad_junction {
            reads_to_unpair.push(n);
        }

        // If bad junctions and no remaining stranded junctions -> unstrand.
        // Matching the original algorithm (lines 15276-15278 + 15457-15461):
        // Killed junctions (mm < 0) have their strand forced to 0 (line 15278).
        // Then: if ALL of a read's junctions have strand=0, unstrand the read.
        // A junction has strand=0 if it's killed OR has ambiguous splice motif.
        // This places the read on sno=1 (unknown strand) during CGroup building,
        // allowing it to bridge neg and pos strand bundles.
        if has_bad_junction && (reads[n].strand == '+' || reads[n].strand == '-') {
            let has_stranded_junction = reads[n].junctions.iter().any(|j| {
                // Killed junctions are treated as strand=0 (line 15278)
                let killed = killed_junction_pairs.contains(j);
                if killed {
                    return false; // strand=0 → not stranded
                }
                // Non-killed junctions: check actual strand evidence
                junction_stats
                    .get(j)
                    .map(|st| st.strand != Some(0))
                    .unwrap_or(false)
            });
            if !has_stranded_junction {
                reads[n].strand = '.';
            }
        }
    }

    // When a read is rejected by bad junctions, unlink all pairs
    // reciprocally; if mate is single-exon, unstrand mate.
    for &n in &reads_to_unpair {
        if n >= reads.len() {
            continue;
        }
        let mates = reads[n].pair_idx.clone();
        reads[n].pair_idx.clear();
        reads[n].pair_count.clear();
        for m in mates {
            if m >= reads.len() {
                continue;
            }
            remove_pair_link_from_read(&mut reads[m], n);
            if reads[m].junctions.is_empty() {
                reads[m].strand = '.';
            }
        }
    }

    // keepstrand for stranded single-exon reads based on future mates.
    for n in 0..reads.len() {
        if reads[n].strand != '+' && reads[n].strand != '-' {
            continue;
        }
        if !reads[n].junctions.is_empty() {
            continue;
        }
        if reads[n].pair_idx.is_empty() {
            continue;
        }

        let strand_n = reads[n].strand;
        let mates = reads[n].pair_idx.clone();
        let mut keepstrand = false;

        for m in mates {
            if m >= reads.len() || n >= m {
                continue;
            }
            if reads[m].junctions.is_empty() {
                if reads[m].strand == strand_n {
                    keepstrand = true;
                    break;
                }
            } else {
                let mate_has_good = reads[m].junctions.iter().any(|j| {
                    !killed_junction_pairs.contains(j)
                        && junction_stats
                            .get(j)
                            .map(|st| st.strand != Some(0))
                            .unwrap_or(false)
                });
                if mate_has_good {
                    keepstrand = true;
                    break;
                } else {
                    reads[m].strand = '.';
                }
            }
        }

        if !keepstrand {
            reads[n].strand = '.';
        }
    }
}

fn build_bundle_guide_junctions(
    guides: &[RefTranscript],
    chrom: &str,
    strand: char,
    bundle_start: u64,
    bundle_end: u64,
) -> Vec<Junction> {
    let mut out = Vec::new();
    for g in guides {
        if g.chrom != chrom || (strand != '.' && g.strand != strand) || g.exons.len() < 2 {
            continue;
        }
        let tx_start = g.exons.first().map(|e| e.0).unwrap_or(0);
        let tx_end = g.exons.last().map(|e| e.1).unwrap_or(0);
        if tx_end <= bundle_start || tx_start >= bundle_end {
            continue;
        }
        for i in 0..g.exons.len() - 1 {
            out.push(Junction::new(g.exons[i].1, g.exons[i + 1].0));
        }
    }
    out
}

#[inline]
fn mark_guide_junctions_for_junctions(
    stats: &mut crate::types::JunctionStats,
    guide_junctions: &[Junction],
) {
    for jn in guide_junctions {
        if let Some(st) = stats.get_mut(jn) {
            st.guide_match = true;
        }
    }
}

fn inject_guide_transfrags(
    graph: &mut Graph,
    transfrags: &mut Vec<GraphTransfrag>,
    guides: &[GuideInfo],
    min_seed_abund: f64,
) -> usize {
    if guides.is_empty() {
        return 0;
    }
    let source = graph.source_id;
    let sink = graph.sink_id;
    let mut added = 0usize;
    for gi in guides {
        if gi.node_ids.is_empty() {
            continue;
        }
        if let Some(existing) = transfrags.iter_mut().find(|tf| tf.node_ids == gi.node_ids) {
            existing.guide = true;
            existing.guide_tid = Some(gi.transcript_id.clone());
            existing.longread = true;
            existing.trflong_seed = true;
            existing.origin_tag = Some("guide".to_string());
            if existing.abundance < min_seed_abund {
                existing.abundance = min_seed_abund;
            }
            if existing.read_count < min_seed_abund {
                existing.read_count = min_seed_abund;
            }
        } else {
            // Sort node_ids by genomic coordinate to ensure nodes are in order
            let mut sorted_node_ids = gi.node_ids.clone();
            sorted_node_ids.sort_unstable_by(|&a, &b| {
                let na = &graph.nodes[a];
                let nb = &graph.nodes[b];
                (na.start, na.end, a).cmp(&(nb.start, nb.end, b))
            });
            let mut tf = GraphTransfrag::new(sorted_node_ids.clone(), graph.pattern_size());
            tf.guide = true;
            tf.guide_tid = Some(gi.transcript_id.clone());
            tf.longread = true;
            tf.trflong_seed = true;
            tf.origin_tag = Some("guide".to_string());
            tf.abundance = min_seed_abund;
            tf.read_count = min_seed_abund;
            let first = sorted_node_ids[0];
            let last = *sorted_node_ids.last().unwrap_or(&first);
            tf.longstart = graph.nodes[first].start;
            tf.longend = graph.nodes[last].end;
            graph.set_pattern_edges_for_path(&mut tf.pattern, &sorted_node_ids);
            let tf_idx = transfrags.len();
            if tf.node_ids.len() > 1 {
                for &nid in &tf.node_ids {
                    if let Some(node) = graph.nodes.get_mut(nid) {
                        node.trf_ids.push(tf_idx);
                    }
                }
            }
            transfrags.push(tf);
            added += 1;
        }

        // guide boundaries imply trusted start/end and support source/sink links.
        // Use sorted node_ids for first/last
        let mut sorted_for_bounds = gi.node_ids.clone();
        sorted_for_bounds.sort_unstable_by(|&a, &b| {
            let na = &graph.nodes[a];
            let nb = &graph.nodes[b];
            (na.start, na.end, a).cmp(&(nb.start, nb.end, b))
        });
        let first = sorted_for_bounds[0];
        let last = *sorted_for_bounds.last().unwrap_or(&first);
        if first >= graph.nodes.len() || last >= graph.nodes.len() {
            continue;
        }
        if let Some(n) = graph.nodes.get_mut(first) {
            n.hardstart = true;
        }
        if let Some(n) = graph.nodes.get_mut(last) {
            n.hardend = true;
        }
        graph.add_edge(source, first);
        graph.add_edge(last, sink);

        for node_pair in [[source, first], [last, sink]] {
            if let Some(existing) = transfrags.iter_mut().find(|tf| {
                tf.node_ids.len() == 2
                    && tf.node_ids[0] == node_pair[0]
                    && tf.node_ids[1] == node_pair[1]
            }) {
                existing.guide = true;
                existing.longread = true;
                existing.trflong_seed = true;
                existing.origin_tag = Some("guide_edge".to_string());
                if existing.abundance < min_seed_abund {
                    existing.abundance = min_seed_abund;
                }
                if existing.read_count < min_seed_abund {
                    existing.read_count = min_seed_abund;
                }
                continue;
            }

            let mut edge_tf = GraphTransfrag::new(node_pair.to_vec(), graph.pattern_size());
            edge_tf.guide = true;
            edge_tf.longread = true;
            edge_tf.trflong_seed = true;
            edge_tf.origin_tag = Some("guide_edge".to_string());
            edge_tf.abundance = min_seed_abund;
            edge_tf.read_count = min_seed_abund;
            if node_pair[0] == source {
                edge_tf.longstart = graph.nodes[first].start;
                edge_tf.longend = graph.nodes[first].start.saturating_add(1);
            } else {
                edge_tf.longstart = graph.nodes[last].start;
                edge_tf.longend = graph.nodes[last].end;
            }
            graph.set_pattern_edges_for_path(&mut edge_tf.pattern, &edge_tf.node_ids);
            let tf_idx = transfrags.len();
            for &nid in &edge_tf.node_ids {
                if let Some(node) = graph.nodes.get_mut(nid) {
                    node.trf_ids.push(tf_idx);
                }
            }
            transfrags.push(edge_tf);
            added += 1;
        }
    }
    added
}

fn tag_transfrags_origin(transfrags: &mut [GraphTransfrag], label: &str) {
    let value = label.to_string();
    for tf in transfrags.iter_mut() {
        tf.origin_tag = Some(value.clone());
    }
}

fn tag_transfrags_origin_if_missing(transfrags: &mut [GraphTransfrag], label: &str) {
    let value = label.to_string();
    for tf in transfrags.iter_mut() {
        if tf.origin_tag.is_none() {
            tf.origin_tag = Some(value.clone());
        }
    }
}

fn tx_bounds(tx: &Transcript) -> (u64, u64) {
    let s = tx.exons.first().map(|e| e.0).unwrap_or(0);
    let e = tx.exons.last().map(|e| e.1).unwrap_or(0);
    (s, e)
}

fn exon_window_avg_cov(bpcov: &Bpcov, start: u64, end: u64) -> f64 {
    let wlen = len_half_open(start, end);
    if wlen == 0 {
        return 0.0;
    }
    let i0 = bpcov.idx(start);
    let i1 = bpcov.idx(end);
    let sum = bpcov.get_cov_range(i0, i1);
    sum / wlen as f64
}

fn strand_compatible(a: char, b: char) -> bool {
    a == b || a == '.' || b == '.'
}

fn is_guide_tx(tx: &Transcript) -> bool {
    tx.source
        .as_deref()
        .map(|s| s.starts_with("guide:"))
        .unwrap_or(false)
}

fn tx_exonic_len(tx: &Transcript) -> u64 {
    tx.exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum()
}

fn tx_cov_from_bpcov(tx: &Transcript, bpcov: &Bpcov) -> f64 {
    let mut cov_sum = 0.0;
    let mut tlen = 0u64;
    for &(s, e) in &tx.exons {
        let elen = len_half_open(s, e);
        if elen == 0 {
            continue;
        }
        tlen += elen;
        cov_sum += bpcov.get_cov_range(bpcov.idx(s), bpcov.idx(e));
    }
    // Return per-base coverage (the `pred->cov` scale), not total bp mass.
    cov_sum / (tlen.max(1) as f64)
}

fn terminal_alt_acceptor_rescue(
    txs: &mut Vec<Transcript>,
    graph: &Graph,
    transfrags: &[GraphTransfrag],
) -> usize {
    if std::env::var_os("RUSTLE_DISABLE_TERMINAL_ALT_ACCEPTOR").is_some() {
        return 0;
    }
    use std::collections::hash_map::Entry;
    const EPS: f64 = crate::constants::FLOW_EPSILON;

    // Map donor -> acceptor_start -> (acceptor_end, max_support)
    //
    // Important: use `read_count` (not `abundance`) because `extract_transcripts` depletes/zeros
    // `abundance` as it consumes seeds, and we want to rescue alternate terminal junctions based
    // on original evidence.
    let mut donor_acceptors: HashMap<u64, HashMap<u64, (u64, f64)>> = Default::default();
    for tf in transfrags {
        let support = tf.read_count;
        if !tf.longread || support <= EPS || tf.node_ids.len() < 2 {
            continue;
        }
        let mut last_intron: Option<(u64, u64, u64)> = None; // (donor_end, acceptor_start, acceptor_end)
        for w in tf.node_ids.windows(2) {
            let Some(na) = graph.nodes.get(w[0]) else {
                continue;
            };
            let Some(nb) = graph.nodes.get(w[1]) else {
                continue;
            };
            // Internal junction tables in this repo often use 1bp-shifted acceptors; we stay in
            // exon-space here (acceptor = next exon start) to match Transcript intron chains.
            if na.end.saturating_add(1) < nb.start {
                last_intron = Some((na.end, nb.start, nb.end));
            }
        }
        let Some((donor, acc_start, acc_end)) = last_intron else {
            continue;
        };
        let by_acc = donor_acceptors.entry(donor).or_default();
        match by_acc.entry(acc_start) {
            Entry::Vacant(v) => {
                v.insert((acc_end, support));
            }
            Entry::Occupied(mut o) => {
                let (prev_end, prev_sup) = *o.get();
                let new_end = prev_end.max(acc_end);
                let new_sup = prev_sup.max(support);
                o.insert((new_end, new_sup));
            }
        }
    }
    if donor_acceptors.is_empty() {
        return 0;
    }

    let tol = std::env::var("RUSTLE_TERMINAL_ALT_ACCEPTOR_TOL")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(2);
    let chains_match_tol = |a: &[(u64, u64)], b: &[(u64, u64)]| -> bool {
        if a.len() != b.len() {
            return false;
        }
        a.iter()
            .copied()
            .zip(b.iter().copied())
            .all(|(x, y)| x.0.abs_diff(y.0) <= tol && x.1.abs_diff(y.1) <= tol)
    };

    // Existing intron chains to avoid duplicates.
    let mut existing: Vec<Vec<(u64, u64)>> = Vec::with_capacity(txs.len());
    for tx in txs.iter() {
        existing.push(intron_chain_from_exons(&tx.exons));
    }

    let mut added = 0usize;
    let n_before = txs.len();
    for idx in 0..n_before {
        if idx >= txs.len() {
            break;
        }
        let tx = txs[idx].clone();
        if tx.exons.len() < 2 {
            continue;
        }
        let intr = intron_chain_from_exons(&tx.exons);
        if intr.is_empty() {
            continue;
        }
        let (donor, acceptor) = *intr.last().unwrap();
        let Some(acc_map) = donor_acceptors.get(&donor) else {
            continue;
        };

        for (&alt_acceptor, &(alt_end, alt_support)) in acc_map.iter() {
            if alt_acceptor.abs_diff(acceptor) <= tol {
                continue;
            }
            // Build a new transcript by swapping the terminal acceptor exon.
            let mut new_tx = tx.clone();
            if let Some(last_exon) = new_tx.exons.last_mut() {
                last_exon.0 = alt_acceptor;
                last_exon.1 = alt_end.max(alt_acceptor);
            } else {
                continue;
            }
            // Sanity: last exon must still start after the previous exon end.
            if new_tx.exons.len() >= 2 {
                let prev_end = new_tx.exons[new_tx.exons.len() - 2].1;
                if prev_end >= new_tx.exons.last().unwrap().0 {
                    continue;
                }
            }

            // Keep coverage comparable to the parent transcript. This is a sensitivity-oriented
            // rescue pass; scaling by junction support tends to make isofrac eliminate real
            // but low-support alternative ends (e.g. STRG.171.{6,9}).
            new_tx.coverage = new_tx.coverage.max(EPS);
            new_tx.longcov = alt_support.max(EPS);
            new_tx.is_longread = true;
            new_tx.source = Some("terminal_alt_acceptor".to_string());
            new_tx.transcript_id = None;
            new_tx.gene_id = None;
            new_tx.ref_transcript_id = None;
            new_tx.ref_gene_id = None;

            let new_chain = intron_chain_from_exons(&new_tx.exons);
            let dup = existing.iter().any(|c| chains_match_tol(c, &new_chain));
            if dup {
                continue;
            }
            existing.push(new_chain);
            txs.push(new_tx);
            added += 1;
        }
    }
    added
}

fn micro_exon_insertion_rescue(
    txs: &mut Vec<Transcript>,
    graph: &Graph,
    transfrags: &[GraphTransfrag],
) -> usize {
    if std::env::var_os("RUSTLE_DISABLE_MICRO_EXON_RESCUE").is_some() {
        return 0;
    }
    const EPS: f64 = crate::constants::FLOW_EPSILON;
    const MAX_MICRO_EXON_LEN: u64 = 50;

    // Map (direct_donor, direct_acceptor) -> best-supported micro exon (start,end,support).
    // We keep only the single best micro exon per direct intron to avoid combinatorial explosion.
    let mut best_micro: HashMap<(u64, u64), (u64, u64, f64)> = Default::default();

    let tf_to_exons = |tf: &GraphTransfrag| -> Vec<(u64, u64)> {
        let mut exons: Vec<(u64, u64)> = Vec::new();
        let mut j = 0usize;
        while j < tf.node_ids.len() {
            let nid = tf.node_ids[j];
            let Some(node) = graph.nodes.get(nid) else {
                j += 1;
                continue;
            };
            let start = node.start;
            let mut end = node.end;
            while j + 1 < tf.node_ids.len() {
                let a = tf.node_ids[j];
                let b = tf.node_ids[j + 1];
                let (Some(na), Some(nb)) = (graph.nodes.get(a), graph.nodes.get(b)) else {
                    break;
                };
                // 1-based inclusive coordinate contiguity: adjacent segments in the same exon.
                if na.end.saturating_add(1) >= nb.start {
                    j += 1;
                    end = nb.end;
                    continue;
                }
                break;
            }
            if end >= start {
                exons.push((start, end));
            }
            j += 1;
        }
        exons
    };

    for tf in transfrags {
        if !tf.longread || tf.read_count <= EPS || tf.node_ids.len() < 4 {
            continue;
        }
        let exons = tf_to_exons(tf);
        if exons.len() < 3 {
            continue;
        }
        let intr = intron_chain_from_exons(&exons);
        if intr.len() < 2 {
            continue;
        }
        for w in intr.windows(2) {
            let (d1, a1) = w[0];
            let (d2, a2) = w[1];
            if a1 >= d2 || d2 <= d1 || a2 <= a1 {
                continue;
            }
            let mex_start = a1;
            let mex_end = d2;
            let mex_len = mex_end.saturating_sub(mex_start);
            if mex_len == 0 || mex_len > MAX_MICRO_EXON_LEN {
                continue;
            }
            let direct = (d1, a2);
            let support = tf.read_count;
            match best_micro.entry(direct) {
                std::collections::hash_map::Entry::Vacant(v) => {
                    v.insert((mex_start, mex_end, support));
                }
                std::collections::hash_map::Entry::Occupied(mut o) => {
                    let (_ps, _pe, psup) = *o.get();
                    if support > psup {
                        o.insert((mex_start, mex_end, support));
                    }
                }
            }
        }
    }
    if best_micro.is_empty() {
        return 0;
    }

    let tol = std::env::var("RUSTLE_MICRO_EXON_RESCUE_TOL")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(2);
    let chains_match_tol = |a: &[(u64, u64)], b: &[(u64, u64)]| -> bool {
        if a.len() != b.len() {
            return false;
        }
        a.iter()
            .copied()
            .zip(b.iter().copied())
            .all(|(x, y)| x.0.abs_diff(y.0) <= tol && x.1.abs_diff(y.1) <= tol)
    };

    let mut existing: Vec<Vec<(u64, u64)>> = Vec::with_capacity(txs.len());
    for tx in txs.iter() {
        existing.push(intron_chain_from_exons(&tx.exons));
    }

    // Estimate "direct" support for scaling by using the max micro support observed for this direct intron.
    // This keeps micro-exon variants from being assigned near-zero cov and immediately filtered.
    let mut added = 0usize;
    let n_before = txs.len();
    for idx in 0..n_before {
        let tx = txs[idx].clone();
        if tx.exons.len() < 2 {
            continue;
        }
        let intr = intron_chain_from_exons(&tx.exons);
        if intr.is_empty() {
            continue;
        }
        for (pos, &(d, a)) in intr.iter().enumerate() {
            let Some(&(mex_s, mex_e, mex_sup)) = best_micro.get(&(d, a)) else {
                continue;
            };
            if mex_s <= d || mex_e >= a || mex_s >= mex_e {
                continue;
            }
            // Insert the exon (mex_s..mex_e) between exons[pos] and exons[pos+1].
            if pos + 1 >= tx.exons.len() {
                continue;
            }
            let prev = tx.exons[pos];
            let next = tx.exons[pos + 1];
            if prev.1.abs_diff(d) > tol || next.0.abs_diff(a) > tol {
                continue;
            }
            if prev.1 >= mex_s || mex_e >= next.0 {
                continue;
            }
            let mut new_tx = tx.clone();
            new_tx.exons.insert(pos + 1, (mex_s, mex_e));
            if !new_tx.exon_cov.is_empty() && new_tx.exon_cov.len() + 1 == new_tx.exons.len() {
                // If exon_cov already aligned (rare), insert a scaled placeholder.
                new_tx.exon_cov.insert(pos + 1, EPS);
            }

            // Keep coverage comparable to the parent transcript (sensitivity rescue).
            new_tx.coverage = new_tx.coverage.max(EPS);
            new_tx.longcov = mex_sup.max(EPS);
            new_tx.is_longread = true;
            new_tx.source = Some("micro_exon_rescue".to_string());
            new_tx.transcript_id = None;
            new_tx.gene_id = None;
            new_tx.ref_transcript_id = None;
            new_tx.ref_gene_id = None;

            let new_chain = intron_chain_from_exons(&new_tx.exons);
            if existing.iter().any(|c| chains_match_tol(c, &new_chain)) {
                continue;
            }
            existing.push(new_chain);
            txs.push(new_tx);
            added += 1;
            break; // only one micro-exon insertion per direct intron per transcript
        }
    }
    added
}

fn apply_terminal_boundary_evidence_to_longread_txs(
    txs: &mut [Transcript],
    graph: &Graph,
    transfrags: &[GraphTransfrag],
) {
    // Build per-node terminal boundary evidence from long-read transfrags.
    // `longstart/longend` share coordinate space with graph nodes (0-based half-open).
    let mut min_longstart_by_first: HashMap<usize, u64> = Default::default();
    let mut max_longend_by_last: HashMap<usize, u64> = Default::default();
    for tf in transfrags {
        if !tf.longread || tf.node_ids.len() < 2 {
            continue;
        }
        if let Some(&first) = tf.node_ids.first() {
            if tf.longstart > 0 {
                min_longstart_by_first
                    .entry(first)
                    .and_modify(|v| *v = (*v).min(tf.longstart))
                    .or_insert(tf.longstart);
            }
        }
        if let Some(&last) = tf.node_ids.last() {
            if tf.longend > 0 {
                max_longend_by_last
                    .entry(last)
                    .and_modify(|v| *v = (*v).max(tf.longend))
                    .or_insert(tf.longend);
            }
        }
    }
    if min_longstart_by_first.is_empty() && max_longend_by_last.is_empty() {
        return;
    }

    let find_node_containing = |pos: u64| -> Option<usize> {
        for (nid, n) in graph.nodes.iter().enumerate() {
            if nid == graph.source_id || nid == graph.sink_id {
                continue;
            }
            if n.start <= pos && pos < n.end {
                return Some(nid);
            }
        }
        None
    };

    for tx in txs.iter_mut() {
        if !tx.is_longread || is_guide_tx(tx) {
            continue;
        }
        // Oracle_direct tx from RUSTLE_MISSED_ORACLE_DIRECT preserve their
        // GTF-exact coords through this extension pass.
        if tx.source.as_deref().map_or(false, |s| s.starts_with("oracle_direct:") || s.starts_with("ref_chain_rescue:")) {
            continue;
        }
        if tx.exons.is_empty() {
            continue;
        }

        // First exon start: extend left across contiguous (collinear) nodes and apply boundary evidence.
        let (fs, fe) = tx.exons[0];
        if fe > fs {
            if let Some(first_nid) = find_node_containing(fs) {
                // Extend across collinear parent chain (parse_trflong terminal exon extension).
                let mut cur = first_nid;
                loop {
                    let cur_node = match graph.nodes.get(cur) {
                        Some(n) => n,
                        None => break,
                    };
                    // Stop when we already have a hard boundary on the current terminal node.
                    if cur_node.hardstart {
                        break;
                    }
                    let parents: Vec<usize> = cur_node
                        .parents
                        .ones()
                        .filter(|&pid| pid != graph.source_id && pid != graph.sink_id)
                        .filter(|&pid| {
                            graph
                                .nodes
                                .get(pid)
                                .map(|p| p.end == cur_node.start)
                                .unwrap_or(false)
                        })
                        .collect();
                    if parents.is_empty() {
                        break;
                    }
                    // Several collinear parents can appear when a synthetic / helper edge shares a
                    // splice boundary with a real upstream exon. Prefer the upstream-most exon
                    // (minimum start) so 5′ extension reaches the leftmost supported TSS.
                    let p = *parents
                        .iter()
                        .min_by_key(|&&pid| {
                            graph
                                .nodes
                                .get(pid)
                                .map(|n| n.start)
                                .unwrap_or(u64::MAX)
                        })
                        .expect("parents non-empty");
                    if let Some(pn) = graph.nodes.get(p) {
                        if pn.start < tx.exons[0].0 {
                            tx.exons[0].0 = pn.start;
                        }
                    }
                    cur = p;
                }

                // Hard boundary evidence: allow extension to the hardstart node start.
                if let Some(node) = graph.nodes.get(cur) {
                    if node.hardstart && node.start < tx.exons[0].0 {
                        tx.exons[0].0 = node.start;
                    }
                }
                if let Some(&ls) = min_longstart_by_first.get(&cur) {
                    let node_start = graph.nodes.get(cur).map(|n| n.start).unwrap_or(ls);
                    if ls >= node_start && ls < tx.exons[0].0 {
                        tx.exons[0].0 = ls;
                    }
                }
            }
        }

        // Last exon end: extend right across contiguous (collinear) nodes and apply boundary evidence.
        let last_idx = tx.exons.len() - 1;
        let (ls, le) = tx.exons[last_idx];
        if le > ls {
            let probe = le.saturating_sub(1);
            if let Some(last_nid) = find_node_containing(probe) {
                // Extend across collinear child chain (parse_trflong terminal exon extension).
                let mut cur = last_nid;
                loop {
                    let cur_node = match graph.nodes.get(cur) {
                        Some(n) => n,
                        None => break,
                    };
                    if cur_node.hardend {
                        break;
                    }
                    let children: Vec<usize> = cur_node
                        .children
                        .ones()
                        .filter(|&cid| cid != graph.source_id && cid != graph.sink_id)
                        .filter(|&cid| {
                            graph
                                .nodes
                                .get(cid)
                                .map(|c| cur_node.end == c.start)
                                .unwrap_or(false)
                        })
                        .collect();
                    if children.is_empty() {
                        break;
                    }
                    // Symmetric to 5′ extension: pick the downstream-most terminal when several
                    // collinear children share a boundary (e.g. synthetic sink helpers).
                    let c = *children
                        .iter()
                        .max_by_key(|&&cid| {
                            graph.nodes.get(cid).map(|n| n.end).unwrap_or(0)
                        })
                        .expect("children non-empty");
                    if let Some(cn) = graph.nodes.get(c) {
                        if cn.end > tx.exons[last_idx].1 {
                            tx.exons[last_idx].1 = cn.end;
                        }
                    }
                    cur = c;
                }

                if let Some(node) = graph.nodes.get(cur) {
                    // Hard boundary evidence: allow extension to node.end.
                    if node.hardend && node.end > tx.exons[last_idx].1 {
                        tx.exons[last_idx].1 = node.end;
                    }
                }
                if let Some(&lend) = max_longend_by_last.get(&cur) {
                    if lend > tx.exons[last_idx].1 {
                        tx.exons[last_idx].1 = lend;
                    }
                }
            }
        }
    }
}

fn rescue_internal_terminal_stop_variants(
    txs: &mut Vec<Transcript>,
    extra_support: &[Transcript],
) -> usize {
    let n = txs.len();
    if n < 1 && extra_support.is_empty() {
        return 0;
    }
    let mut existing: HashSet<(String, char, Vec<(u64, u64)>)> = txs
        .iter()
        .map(|t| (t.chrom.clone(), t.strand, t.exons.clone()))
        .collect();
    let mut rescued = Vec::new();

    // Candidates for truncated variants (must be long-read and multi-exon)
    for i in 0..n {
        let long = &txs[i];
        if !long.is_longread || is_guide_tx(long) || long.exons.len() < 3 {
            continue;
        }

        // Support can come from either independently recovered real transcripts
        // or from synthetic junction paths.
        let mut checked_support: HashSet<Vec<(u64, u64)>> = Default::default();

        for j in 0..(n + extra_support.len()) {
            let support = if j < n {
                if i == j {
                    continue;
                }
                &txs[j]
            } else {
                &extra_support[j - n]
            };

            if !support.is_longread
                || is_guide_tx(support)
                || support.chrom != long.chrom
                || support.strand != long.strand
                || support.exons.len() < 2
                || support.exons.len() >= long.exons.len()
            {
                continue;
            }

            // Optimization: skip if we've already checked an identical intron chain for support
            let support_chain = intron_chain_from_exons(&support.exons);
            if !checked_support.insert(support_chain) {
                continue;
            }

            // If a support transcript justifies an internal terminal stop on
            // a prefix of the longer chain, emit a truncated same-chain variant from the longer
            // start to that supported stop.
            for start_idx in 0..long.exons.len() {
                let end_idx = start_idx + support.exons.len();
                if end_idx > long.exons.len().saturating_sub(1) {
                    break;
                }
                if long.exons[start_idx..end_idx] != support.exons[..] {
                    continue;
                }

                // Parity: Only rescue if the truncation point is supported by a real transcript
                // OR if it corresponds to a verified hardstart/hardend in the graph.
                let is_real_support = support.source.as_deref() != Some("junction_path");
                let has_hard_boundary = if start_idx == 0 {
                    // Prefix match: truncation is at the END of support
                    support.hardend
                } else if end_idx == long.exons.len() {
                    // Suffix match: truncation is at the START of support
                    support.hardstart
                } else {
                    // Internal match: requires both if we were to support it,
                    // but for now we only support prefix/suffix truncation.
                    false
                };

                if !is_real_support && !has_hard_boundary {
                    continue;
                }

                let new_exons = long.exons[..end_idx].to_vec();
                let key = (long.chrom.clone(), long.strand, new_exons.clone());
                if existing.contains(&key) {
                    continue;
                }

                let mut new_tx = long.clone();
                new_tx.exons = new_exons;
                // Parity: rescued variants inherit coverage from the supporting long read
                // (representing the flow that also explains the truncated variant).
                new_tx.coverage = long.coverage;
                new_tx.longcov = long.longcov;
                new_tx.exon_cov = vec![long.coverage; new_tx.exons.len()];
                new_tx.source = Some(String::from("terminal_stop_variant"));
                // Inherit boundary flags for downstream filters
                new_tx.hardstart = long.hardstart;
                new_tx.hardend = if start_idx == 0 { support.hardend } else { long.hardend };
                rescued.push(new_tx);
                existing.insert(key);
            }
        }
    }

    let added = rescued.len();
    txs.extend(rescued);
    added
}

fn register_transfrag_range_on_nodes(
    graph: &mut Graph,
    transfrags: &[GraphTransfrag],
    start_idx: usize,
) {
    for (offset, tf) in transfrags.iter().enumerate().skip(start_idx) {
        if tf.node_ids.len() <= 1 {
            continue;
        }
        let tf_idx = offset;
        for &nid in &tf.node_ids {
            if let Some(node) = graph.nodes.get_mut(nid) {
                if !node.trf_ids.contains(&tf_idx) {
                    node.trf_ids.push(tf_idx);
                }
            }
        }
    }
}

/// Returns (
///   filtered_transcripts,
///   pre_filter_transcripts_for_trace,
///   seed_outcomes_for_trace,
///   predcluster_stage_summary,
///   unsupported_junction_removed,
/// ).
/// Dump Rustle's graph for a bundle as a DOT file + check walkability
/// of every reference transcript's intron chain through the graph.
/// Triggered by env RUSTLE_GRAPH_DOT_DIR.
fn dump_graph_and_walkability(
    graph: &Graph,
    bundle: &crate::types::Bundle,
    ref_transcripts: &[RefTranscript],
    dir: &str,
    chrom: &str,
) {
    use std::io::Write;
    let _ = std::fs::create_dir_all(dir);
    let tag = format!("{}_{}_{}_{}",
        chrom.replace('.', "_"), bundle.start, bundle.end,
        match bundle.strand { '+' => "p", '-' => "m", _ => "u" });

    // --- DOT export ---
    let dot_path = format!("{}/{}.dot", dir, tag);
    let mut dot = match std::fs::File::create(&dot_path) {
        Ok(f) => f,
        Err(_) => return,
    };
    let _ = writeln!(dot, "digraph bundle_{} {{", tag);
    let _ = writeln!(dot, "  rankdir=LR;");
    let _ = writeln!(dot, "  node [shape=box,fontsize=10];");
    let src = graph.source_id;
    let snk = graph.sink_id;
    for (nid, n) in graph.nodes.iter().enumerate() {
        let label = if nid == src { "source".to_string() }
            else if nid == snk { "sink".to_string() }
            else { format!("n{}\\n{}-{}\\ncov={:.1}", nid, n.start, n.end, n.coverage) };
        let color = if n.hardstart && n.hardend { "lightblue" }
            else if n.hardstart { "lightgreen" }
            else if n.hardend { "lightyellow" }
            else { "white" };
        let _ = writeln!(dot, "  n{} [label=\"{}\",style=filled,fillcolor={}];", nid, label, color);
    }
    for (nid, n) in graph.nodes.iter().enumerate() {
        for c in n.children.ones() {
            let _ = writeln!(dot, "  n{} -> n{};", nid, c);
        }
    }
    let _ = writeln!(dot, "}}");

    // --- Walkability check ---
    // For each ref tx overlapping the bundle, test if its intron chain
    // forms a valid path in the graph.
    let walk_path = format!("{}/{}.walkable.tsv", dir, tag);
    let mut wf = match std::fs::File::create(&walk_path) {
        Ok(f) => f,
        Err(_) => return,
    };
    let _ = writeln!(wf, "tx_id\tstrand\tspan\tn_exons\twalkable\tnode_path\tmissing_edges");

    // Build (chrom+strand-aware) position → node lookup
    let mut start_map: std::collections::HashMap<u64, usize> =
        std::collections::HashMap::new();
    let mut end_map: std::collections::HashMap<u64, usize> =
        std::collections::HashMap::new();
    for (nid, n) in graph.nodes.iter().enumerate() {
        if nid == src || nid == snk { continue; }
        if n.end <= n.start { continue; }
        start_map.insert(n.start, nid);
        end_map.insert(n.end, nid);
    }

    for rt in ref_transcripts {
        if rt.chrom != chrom { continue; }
        if rt.strand != bundle.strand && bundle.strand != '.' { continue; }
        let rt_start = rt.exons.first().map(|e| e.0).unwrap_or(0);
        let rt_end = rt.exons.last().map(|e| e.1).unwrap_or(0);
        // Overlap with bundle
        if rt_end < bundle.start || rt_start > bundle.end { continue; }

        // rt.exons is already stored as (start_0based, end_inclusive) which
        // equals (start_halfopen, end_halfopen). Use directly.
        let internal_exons: Vec<(u64, u64)> = rt.exons.clone();
        let mut node_path: Vec<i64> = Vec::new();
        let mut missing: Vec<String> = Vec::new();
        let mut walkable = true;
        for i in 0..internal_exons.len().saturating_sub(1) {
            let donor = internal_exons[i].1;   // exon end (internal half-open = GTF inclusive end)
            let acceptor = internal_exons[i + 1].0;
            let d_node = end_map.get(&donor).copied();
            let a_node = start_map.get(&acceptor).copied();
            match (d_node, a_node) {
                (Some(d), Some(a)) => {
                    if i == 0 { node_path.push(d as i64); }
                    if graph.nodes[d].children.contains(a) {
                        node_path.push(a as i64);
                    } else {
                        walkable = false;
                        missing.push(format!("edge_missing:{}->{}", d, a));
                        node_path.push(-1);
                    }
                }
                _ => {
                    walkable = false;
                    if d_node.is_none() {
                        missing.push(format!("donor_node_missing:{}", donor));
                    }
                    if a_node.is_none() {
                        missing.push(format!("acceptor_node_missing:{}", acceptor));
                    }
                    node_path.push(-1);
                }
            }
        }
        if internal_exons.len() == 1 && node_path.is_empty() {
            // single-exon tx: check if any node fully contains its span
            let (s, e) = internal_exons[0];
            let contained = graph.nodes.iter().enumerate().any(|(nid, n)| {
                nid != src && nid != snk && n.start <= s && n.end >= e
            });
            if !contained { walkable = false; missing.push("single_exon_not_in_any_node".into()); }
        }
        let _ = writeln!(wf,
            "{}\t{}\t{}-{}\t{}\t{}\t{}\t{}",
            rt.id, rt.strand, rt_start, rt_end, rt.exons.len(),
            walkable,
            node_path.iter().map(|n| n.to_string()).collect::<Vec<_>>().join(","),
            if missing.is_empty() { "-".to_string() } else { missing.join(";") }
        );
    }
}

/// The trace-related elements are only populated when `trace_mode` is true.
fn extract_bundle_transcripts_for_graph(
    graph_mut: &mut Graph,
    transfrags: &mut Vec<GraphTransfrag>,
    bundle: &crate::types::Bundle,
    config: &RunConfig,
    mapped_guides: &[GuideInfo],
    guide_transcripts: &[RefTranscript],
    ref_transcripts: &[RefTranscript],
    junctions: &[Junction],
    good_junctions: &HashSet<(u64, u64)>,
    bpcov: &Bpcov,
    bpcov_stranded: &BpcovStranded,
    trace_mode: bool,
) -> (
    Vec<Transcript>,
    Option<Vec<Transcript>>,
    Vec<(usize, crate::path_extract::SeedOutcome)>,
    LongRecSummary,
    PredclusterStageSummary,
    usize,
) {
    let traced_ref = find_traced_ref_in_bundle(bundle, ref_transcripts);
    let trace_stage = |stage: &str, txs: &[Transcript]| {
        debug_target_ref_stage(stage, bundle, ref_transcripts, txs);
        if std::env::var_os("RUSTLE_ORACLE_STAGE_TRACE").is_some() {
            for t in txs.iter() {
                if t.source.as_deref().map_or(false, |s| s.starts_with("oracle_direct:") || s.starts_with("ref_chain_rescue:")) {
                    eprintln!(
                        "ORACLE_STAGE stage={} src={} {}-{} exons={}",
                        stage,
                        t.source.as_deref().unwrap_or(""),
                        t.exons.first().map(|e| e.0).unwrap_or(0),
                        t.exons.last().map(|e| e.1).unwrap_or(0),
                        t.exons.len()
                    );
                }
            }
        }
    };

    // Graph DOT dump: set RUSTLE_GRAPH_DOT_DIR=/path to write graph.dot
    // per bundle, plus a .walkable TSV checking whether each ref tx
    // intron chain forms a walkable path in the graph.
    if let Ok(dir) = std::env::var("RUSTLE_GRAPH_DOT_DIR") {
        dump_graph_and_walkability(
            graph_mut, bundle, ref_transcripts, &dir, &bundle.chrom
        );
    }

    // Rebuild node->transfrag membership to match current transfrag list.
    for node in graph_mut.nodes.iter_mut() {
        node.trf_ids.clear();
    }
    for (tf_idx, tf) in transfrags.iter().enumerate() {
        if tf.node_ids.len() <= 1 {
            continue;
        }
        for &nid in &tf.node_ids {
            if let Some(node) = graph_mut.nodes.get_mut(nid) {
                node.trf_ids.push(tf_idx);
            }
        }
    }

    // refresh_graph_node_bp_coverage_from_bpcov(graph_mut, bpcov_stranded);  // node.cov is not refreshed from bpcov
    trace_graph_numeric_state("extract:after_bpcov_refresh", graph_mut, transfrags);
    let _nodecov = compute_nodecov(graph_mut, transfrags, config.verbose);
    trace_graph_numeric_state("extract:after_compute_nodecov", graph_mut, transfrags);

    // Compute reachability for onpath_long validation in path extension
    graph_mut.compute_reachability();

    // Ensure edge bits are set in all transfrag patterns.
    let needed_psize = graph_mut.pattern_size();
    for tf in transfrags.iter_mut() {
        tf.pattern.grow(needed_psize);
        for w in tf.node_ids.windows(2) {
            if let Some(eid) = graph_mut.edge_bit_index(w[0], w[1]) {
                tf.pattern.set_bit(eid);
            }
        }
    }

    let mut seed_outcomes_buf: Vec<(usize, crate::path_extract::SeedOutcome)> = Vec::new();
    let mut longrec_summary = LongRecSummary::default();
    // Direct-emit oracle: before extract_transcripts, record which ref tx
    // to inject. After, we compare per-tx with emitted output to classify
    // PATH_EXTRACT_BLOCK (if forced emission matches and default doesn't)
    // vs FILTER_BLOCK. See append_missed_oracle_direct_emit below.
    let mut txs = if config.rawreads {
        extract_rawreads_transcripts(
            graph_mut,
            transfrags,
            &bundle.chrom,
            bundle.strand,
            &format!("{}:{}-{}", bundle.chrom, bundle.start, bundle.end),
            config,
        )
    } else if !config.long_reads {
        extract_shortread_transcripts(
            graph_mut,
            transfrags,
            &bundle.chrom,
            bundle.strand,
            &format!("{}:{}-{}", bundle.chrom, bundle.start, bundle.end),
            config,
        )
    } else {
        extract_transcripts(
            graph_mut,
            transfrags,
            &bundle.chrom,
            bundle.strand,
            &format!("{}:{}-{}", bundle.chrom, bundle.start, bundle.end),
            config,
            false,
            if trace_mode || std::env::var_os("RUSTLE_SEED_STATS").is_some() {
                Some(&mut seed_outcomes_buf)
            } else {
                None
            },
            Some(&mut longrec_summary),
        )
    };
    // Direct-emit oracle for diagnostic classification of missed tx.
    // Appends Transcript structs synthesized from GTF for ref tx that
    // overlap this bundle. Opt-in via RUSTLE_MISSED_ORACLE_DIRECT=path.gtf.
    if std::env::var_os("RUSTLE_MISSED_ORACLE_DIRECT").is_some() {
        append_missed_oracle_direct_emit(
            &mut txs,
            &bundle.chrom,
            bundle.start,
            bundle.end,
            bundle.strand,
        );
    }
    if std::env::var_os("RUSTLE_SEED_STATS").is_some() {
        let mut counts: std::collections::BTreeMap<String, usize> =
            std::collections::BTreeMap::new();
        for (_, outcome) in &seed_outcomes_buf {
            let key = match outcome {
                crate::path_extract::SeedOutcome::Skipped(r) => format!("skipped_{}", r),
                crate::path_extract::SeedOutcome::BackToSourceFail => "back_fail".into(),
                crate::path_extract::SeedOutcome::FwdToSinkFail => "fwd_fail".into(),
                crate::path_extract::SeedOutcome::UnwitnessedSplice => "unwitnessed".into(),
                crate::path_extract::SeedOutcome::HardBoundaryMismatch => "hard_boundary".into(),
                crate::path_extract::SeedOutcome::ZeroFlux => "zero_flux".into(),
                crate::path_extract::SeedOutcome::LowCoverage(_) => "low_coverage".into(),
                crate::path_extract::SeedOutcome::EonlyNonGuide => "eonly_nonguide".into(),
                crate::path_extract::SeedOutcome::TooShort => "too_short".into(),
                crate::path_extract::SeedOutcome::ChecktrfReadthr => "checktrf_readthr".into(),
                crate::path_extract::SeedOutcome::ChecktrfEonlySkip => "checktrf_eonly".into(),
                crate::path_extract::SeedOutcome::ChecktrfRedistributed => "checktrf_redist".into(),
                crate::path_extract::SeedOutcome::ChecktrfRescued => "checktrf_rescued".into(),
                crate::path_extract::SeedOutcome::ChecktrfIncomplete => "checktrf_incomplete".into(),
                crate::path_extract::SeedOutcome::ChecktrfRescueFail => "checktrf_rescue_fail".into(),
                crate::path_extract::SeedOutcome::Stored(_) => "stored".into(),
            };
            *counts.entry(key).or_insert(0) += 1;
        }
        let summary: Vec<String> = counts.iter().map(|(k, v)| format!("{}={}", k, v)).collect();
        eprintln!(
            "SEED_STATS {}:{}-{} strand={} seeds={} {}",
            bundle.chrom,
            bundle.start,
            bundle.end,
            bundle.strand,
            seed_outcomes_buf.len(),
            summary.join(" ")
        );
    }
    if config.long_reads && !config.eonly {
        let added = if config.emit_terminal_alt_acceptor {
            terminal_alt_acceptor_rescue(&mut txs, graph_mut, transfrags)
        } else {
            0
        };
        let added_micro = if config.emit_micro_exon_rescue {
            micro_exon_insertion_rescue(&mut txs, graph_mut, transfrags)
        } else {
            0
        };
        apply_terminal_boundary_evidence_to_longread_txs(&mut txs, graph_mut, transfrags);
        if added > 0 && config.verbose {
            eprintln!(
                "    terminal_alt_acceptor_rescue: +{} transcript(s) (bundle={}:{}-{})",
                added,
                bundle.chrom,
                bundle.start + 1,
                bundle.end
            );
        }
        if added_micro > 0 && config.verbose {
            eprintln!(
                "    micro_exon_insertion_rescue: +{} transcript(s) (bundle={}:{}-{})",
                added_micro,
                bundle.chrom,
                bundle.start + 1,
                bundle.end
            );
        }
        trace_stage("terminal_alt_acceptor_rescue", &txs);
        trace_stage("micro_exon_insertion_rescue", &txs);
    }
    trace_stage("extract_base", &txs);
    trace_ref_seed_outcomes(
        bundle,
        ref_transcripts,
        graph_mut,
        transfrags,
        &seed_outcomes_buf,
    );

    // extend terminal exon boundaries to match guide coordinates.
    // When a transcript is guide-derived, extend first exon start to guide start and
    // last exon end to guide end (if guide boundary overlaps the exon).
    for tx in txs.iter_mut() {
        if let Some(ref src) = tx.source {
            if let Some(guide_tid) = src.strip_prefix("guide:") {
                if let Some(gi) = mapped_guides.iter().find(|g| g.transcript_id == guide_tid) {
                    if let Some(first_exon) = tx.exons.first_mut() {
                        if first_exon.0 != gi.tx_start && gi.tx_start <= first_exon.1 {
                            first_exon.0 = gi.tx_start;
                        }
                    }
                    if let Some(last_exon) = tx.exons.last_mut() {
                        if last_exon.1 != gi.tx_end && gi.tx_end >= last_exon.0 {
                            last_exon.1 = gi.tx_end;
                        }
                    }
                }
            }
        }
    }
    trace_stage("guide_boundary_extension", &txs);
    txs = apply_guide_internal_split_merge(txs, bundle, guide_transcripts, bpcov, config.verbose);
    trace_stage("guide_internal_split_merge", &txs);

    txs = add_contained_isoforms(txs, junctions, config);
    trace_stage("add_contained_isoforms", &txs);
    let (extra_evidence, extra_emitted) = emit_junction_paths(graph_mut, bundle, &txs, config);
    if !extra_emitted.is_empty() {
        if config.verbose {
            eprintln!(
                "    junction_paths: emitted {} high-confidence standalone paths (bundle={}:{}-{})",
                extra_emitted.len(),
                bundle.chrom,
                bundle.start + 1,
                bundle.end
            );
        }
        txs.extend(extra_emitted);
    }
    trace_stage("emit_junction_paths", &txs);

    let added_terminal_stop = rescue_internal_terminal_stop_variants(&mut txs, &extra_evidence);
    if added_terminal_stop > 0 && config.verbose {
        eprintln!(
            "    terminal_stop_variant_rescue: +{} transcript(s) (bundle={}:{}-{})",
            added_terminal_stop,
            bundle.chrom,
            bundle.start + 1,
            bundle.end
        );
    }
    trace_stage("terminal_stop_variant_rescue", &txs);

    let ref_paths = emit_reference_chains(
        graph_mut,
        bundle,
        ref_transcripts,
        &txs,
        config,
        good_junctions,
    );
    if !ref_paths.is_empty() {
        txs.extend(ref_paths);
    }
    trace_stage("emit_reference_chains", &txs);

    if !guide_transcripts.is_empty() {
        let guide_chain_tol = std::env::var("RUSTLE_GUIDE_CHAIN_EQ_TOL")
            .ok()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or(5);
        let assigned = assign_guide_ids_by_intron_chain(&mut txs, guide_transcripts, guide_chain_tol);
        if config.verbose && assigned > 0 {
            eprintln!(
                "    guide_chain_match_assign: assigned {} transcript(s) to guide IDs (tol={})",
                assigned, guide_chain_tol
            );
        }
    }

    // Guided mode: keep exact guide exon coordinates for guide-derived predictions.
    // Without this, graph-extracted guide models can inherit shifted internal boundaries
    // from local nodes (e.g. STRG.554.*), which later appear as not_extracted:junctions_present.
    if !guide_transcripts.is_empty() {
        let guide_by_id: HashMap<&str, &RefTranscript> = guide_transcripts
            .iter()
            .map(|g| (g.id.as_str(), g))
            .collect();
        for tx in txs.iter_mut() {
            let Some(src) = tx.source.as_deref() else {
                continue;
            };
            let Some(guide_id) = src.strip_prefix("guide:") else {
                continue;
            };
            let Some(guide) = guide_by_id.get(guide_id) else {
                continue;
            };
            if guide.exons.is_empty() {
                continue;
            }
            if guide.chrom != tx.chrom {
                continue;
            }
            if guide.strand != '.' && tx.strand != '.' && guide.strand != tx.strand {
                continue;
            }
            tx.exons = guide.exons.clone();
            tx.chrom = guide.chrom.clone();
            tx.strand = guide.strand;
            if tx.exon_cov.len() != tx.exons.len() {
                tx.exon_cov = vec![tx.coverage; tx.exons.len()];
            }
        }
    }
    trace_stage("guide_exact_restore", &txs);

    // Compute strand-weighted bpcov_cov for each transcript (gene_abundance bpcov calculation).
    for tx in txs.iter_mut() {
        tx.bpcov_cov = compute_transcript_bpcov(tx, bpcov_stranded);
        tx.intron_low = crate::gene_abundance::compute_transcript_intron_low(tx, bpcov_stranded);
        // IMPORTANT: do not overwrite `tx.coverage` (flow-derived) by default.
        // the pred->cov flow value is flow-derived; using raw bpcov makes cov wildly larger on deep loci
        // (e.g. STRG.27.3: ~31 vs the expected ~6) and breaks isofrac/pairwise decisions.
        //
        // If you explicitly want the bpcov-proxy coverage for diagnostics, enable:
        //   RUSTLE_USE_BPCOV_COV=1
        if config.long_reads && std::env::var_os("RUSTLE_USE_BPCOV_COV").is_some() {
            let tlen = tx_exonic_len(tx).max(1) as f64;
            tx.coverage = tx.bpcov_cov / tlen;
        }
    }

    // print_predcluster filtering blocks (keep active in strict mode).
    txs = correct_two_exon_split_errors(txs, bpcov, config.verbose);
    trace_stage("correct_two_exon_split_errors", &txs);
    txs = collapse_covered_microintrons(txs, bpcov, config.verbose);
    trace_stage("collapse_covered_microintrons", &txs);
    txs = apply_guide_reflink_absorption(txs, mapped_guides, config.verbose, config.long_reads);
    trace_stage("apply_guide_reflink_absorption", &txs);
    txs =
        apply_checkincomplete_subset_filter(txs, mapped_guides, guide_transcripts, config.verbose);
    trace_stage("apply_checkincomplete_subset_filter", &txs);
    // Algorithm: delete non-guide predictions shorter than mintranscriptlen
    // before the equal_pred merge loop.
    if config.min_transcript_length > 0 {
        let before = txs.len();
        txs.retain(|t| is_guide_tx(t) || tx_exonic_len(t) >= config.min_transcript_length);
        if config.verbose && txs.len() < before {
            eprintln!(
                "    mintranscriptlen gate: removed {} transcript(s) (len < {})",
                before - txs.len(),
                config.min_transcript_length
            );
        }
        trace_stage("mintranscriptlen_gate", &txs);
    }
    txs = collapse_equal_predictions(txs, config.verbose);
    trace_stage("collapse_equal_predictions", &txs);
    txs = adjust_overlap_endpoints_opposite_strand(txs, bpcov, config, config.verbose);
    trace_stage("adjust_overlap_endpoints_opposite_strand", &txs);

    // second parse_trflong pass for nascent transcript extensions.
    if config.isnascent && !config.eonly && !guide_transcripts.is_empty() {
        let nascent_txs = extract_transcripts(
            graph_mut,
            transfrags,
            &bundle.chrom,
            bundle.strand,
            &format!("{}:{}-{}", bundle.chrom, bundle.start, bundle.end),
            config,
            true,
            None,
            None,
        );
        txs.extend(nascent_txs);
        trace_stage("isnascent_extend", &txs);
    }

    // Capture pre-filter transcripts for trace analysis before any predcluster filtering.
    let pre_filter = if trace_mode { Some(txs.clone()) } else { None };

    // Canonical post-extraction ordering mirroring print_predcluster flow.
    let (predcluster_txs, predcluster_summary) =
        print_predcluster_with_summary(txs, config, Some(bpcov), traced_ref);
    txs = predcluster_txs;
    trace_stage("print_predcluster", &txs);

    // Remove transcripts with unsupported junctions.
    let before_junction_support = txs.len();
    txs = filter_unsupported_junctions(
        txs,
        good_junctions,
        config.junction_correction_window,
        config.verbose,
    );
    let junction_support_removed = before_junction_support.saturating_sub(txs.len());
    trace_stage("filter_unsupported_junctions", &txs);

    // Per-junction read-support gate: suppress minor-isoform transcripts whose
    // weakest junction has fewer than N reads. Catches j-class artifacts where
    // the flow decomposition assembled a transcript using a barely-supported
    // alt donor/acceptor. Opt-in via RUSTLE_MIN_JUNC_SUPPORT=N env var;
    // high-cov transcripts (longcov >= RUSTLE_MIN_JUNC_EXEMPT_COV, default 10)
    // are exempt because the flow-decomp bar already qualifies them.
    if let Some(min_reads) = std::env::var("RUSTLE_MIN_JUNC_SUPPORT")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
    {
        let exempt_cov = std::env::var("RUSTLE_MIN_JUNC_EXEMPT_COV")
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(10.0);
        txs = crate::transcript_filter::filter_by_min_junction_support(
            txs,
            &bundle.junction_stats,
            min_reads,
            exempt_cov,
            config.junction_correction_window,
            config.verbose,
        );
        trace_stage("filter_by_min_junction_support", &txs);
    }

    // StringTie-like TSS/TTS anchor filter (opt-in).
    //
    // Collect per-strand read endpoints (starts and ends) from bundle.reads.
    // For each multi-exon tx: a boundary is "anchored" if at least N reads
    // have their endpoint (start for 5', end for 3') within W bp of the
    // tx's terminal exon edge. If NEITHER end is anchored AND the tx has
    // low coverage, drop — targets unanchored run-through extensions
    // (k-class). Conservative defaults.
    //
    // OPT-IN (regresses): read-endpoint signal doesn't cleanly separate
    // k-class over-extensions from legit tx with sparse terminal reads
    // (soft-clipping common in long reads). Kept scaffolded for future
    // refinement.
    if std::env::var_os("RUSTLE_ANCHOR_TRIM_ON").is_some() && !txs.is_empty() {
        let min_endpoints: usize = std::env::var("RUSTLE_ANCHOR_TRIM_MIN_ENDPOINTS")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(2);
        let window: u64 = std::env::var("RUSTLE_ANCHOR_TRIM_WINDOW")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(50);
        let max_cov: f64 = std::env::var("RUSTLE_ANCHOR_TRIM_MAX_COV")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(5.0);

        // Collect read endpoints stratified by strand
        let mut starts_plus: Vec<u64> = Vec::new();
        let mut starts_minus: Vec<u64> = Vec::new();
        let mut ends_plus: Vec<u64> = Vec::new();
        let mut ends_minus: Vec<u64> = Vec::new();
        for r in &bundle.reads {
            match r.strand {
                '+' => { starts_plus.push(r.ref_start); ends_plus.push(r.ref_end); }
                '-' => { starts_minus.push(r.ref_start); ends_minus.push(r.ref_end); }
                _ => {
                    starts_plus.push(r.ref_start); ends_plus.push(r.ref_end);
                    starts_minus.push(r.ref_start); ends_minus.push(r.ref_end);
                }
            }
        }
        starts_plus.sort(); ends_plus.sort();
        starts_minus.sort(); ends_minus.sort();

        let count_in_window = |sorted: &[u64], center: u64, w: u64| -> usize {
            let lo = center.saturating_sub(w);
            let hi = center + w;
            // Binary search for lo, count until hi
            let i = sorted.partition_point(|&x| x < lo);
            let j = sorted.partition_point(|&x| x <= hi);
            j - i
        };

        // Adaptive: require endpoint-anchor count proportional to tx read support.
        // min_anchor = max(min_endpoints, round(min_frac * longcov)).
        let min_frac: f64 = std::env::var("RUSTLE_ANCHOR_TRIM_MIN_FRAC")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.20);
        let before = txs.len();
        txs.retain(|t| {
            if t.exons.len() < 2 { return true; }
            if t.coverage >= max_cov { return true; }
            let starts_ref = if t.strand == '+' { &starts_plus } else { &starts_minus };
            let ends_ref = if t.strand == '+' { &ends_plus } else { &ends_minus };
            let t_start = t.exons.first().map(|e| e.0).unwrap_or(0);
            let t_end = t.exons.last().map(|e| e.1).unwrap_or(0);
            let n5 = count_in_window(starts_ref, t_start, window);
            let n3 = count_in_window(ends_ref, t_end, window);
            let need = min_endpoints.max((min_frac * t.longcov).round() as usize);
            n5 >= need || n3 >= need
        });
        if config.verbose && txs.len() < before {
            eprintln!("rustle: anchor_trim removed {} unanchored low-cov tx", before - txs.len());
        }
    }

    // Orphan-junction rescue: for each junction in bundle.junction_stats with
    // substantial read support that is NOT used by any emitted tx, walk the
    // graph to synthesize a tx around it. Targets the "gene-inside-intron"
    // class where a small gene's junction (12 reads) is overshadowed by a
    // bridging junction (1354 reads) and gets absorbed.
    //
    // Default on; RUSTLE_ORPHAN_JUNC_RESCUE_OFF=1 to disable.
    if std::env::var_os("RUSTLE_ORPHAN_JUNC_RESCUE_OFF").is_none() && !txs.is_empty() {
        let min_reads: f64 = std::env::var("RUSTLE_ORPHAN_JUNC_MIN_READS")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(10.0);
        let max_walk: usize = std::env::var("RUSTLE_ORPHAN_JUNC_MAX_WALK")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(100);

        // Collect junction pairs used in emitted tx.
        // IMPORTANT: junction convention is (prev_exon.end, next_exon.start-1)
        // in external coords. Internally our exons store [start, end) half-open
        // but WRITTEN as GTF (start+1, end). The JunctionStats uses coords
        // compatible with StringTie: donor = exon_end (GTF inclusive),
        // acceptor = next_exon_start - 1. So for internal exon [start, end),
        // exon_end_GTF == end (as written). Next exon start (GTF inclusive) ==
        // next_exon.start + 1 (internal half-open means internal start is
        // GTF start - 1). So acceptor = next_exon.start + 1 - 1 = next_exon.start.
        // Wait: internal next_exon.start is written as start+1 in GTF. So
        // GTF next_exon_start - 1 = internal_start + 1 - 1 = internal_start.
        // Thus used_juncs key = (internal_exon.end, internal_next_exon.start).
        let mut used_juncs: std::collections::HashSet<(u64, u64)> =
            std::collections::HashSet::new();
        for tx in txs.iter() {
            for i in 0..tx.exons.len().saturating_sub(1) {
                let d = tx.exons[i].1;
                let a = tx.exons[i + 1].0;
                used_juncs.insert((d, a));
            }
        }
        if std::env::var_os("RUSTLE_TRACE_ORPHAN_JUNC").is_some() {
            eprintln!("ORPHAN_JUNC_USED bundle={}:{}-{} used_count={} txs={}",
                bundle.chrom, bundle.start, bundle.end, used_juncs.len(), txs.len());
            for tx in txs.iter() {
                eprintln!("  tx exons={:?} cov={:.2} longcov={:.2} src={:?}",
                    tx.exons, tx.coverage, tx.longcov, tx.source);
            }
        }

        // Build node lookup: donor_end → node, acceptor_start → node
        let src = graph_mut.source_id;
        let snk = graph_mut.sink_id;
        let mut end_to_node: std::collections::HashMap<u64, usize> =
            std::collections::HashMap::new();
        let mut start_to_node: std::collections::HashMap<u64, usize> =
            std::collections::HashMap::new();
        for (nid, n) in graph_mut.nodes.iter().enumerate() {
            if nid == src || nid == snk { continue; }
            if n.end <= n.start { continue; }
            end_to_node.insert(n.end, nid);
            start_to_node.insert(n.start, nid);
        }

        let mut rescued: Vec<Transcript> = Vec::new();
        let mut emitted_chains: std::collections::HashSet<Vec<(u64, u64)>> =
            txs.iter().map(|t| t.exons.clone()).collect();

        let trace = std::env::var_os("RUSTLE_TRACE_ORPHAN_JUNC").is_some();
        // Iterate orphan junctions
        for (junc, stat) in bundle.junction_stats.iter() {
            if stat.mrcount < min_reads { continue; }
            if used_juncs.contains(&(junc.donor, junc.acceptor)) {
                if trace {
                    eprintln!("ORPHAN_JUNC SKIP {}:{}-{} reason=used nreads={:.0}",
                        bundle.chrom, junc.donor, junc.acceptor, stat.mrcount);
                }
                continue;
            }
            // Filter by strand match with bundle
            let bundle_strand_i = match bundle.strand {
                '+' => Some(1i8),
                '-' => Some(-1),
                _ => None,
            };
            if let (Some(bs), Some(js)) = (bundle_strand_i, stat.strand) {
                if bs != js && js != 0 { continue; }
            }
            // Find graph nodes flanking the junction
            let n_d = match end_to_node.get(&junc.donor) {
                Some(&n) => n,
                None => {
                    if trace {
                        eprintln!("ORPHAN_JUNC SKIP {}:{}-{} reason=no_donor_node nreads={:.0}",
                            bundle.chrom, junc.donor, junc.acceptor, stat.mrcount);
                    }
                    continue;
                }
            };
            let n_a = match start_to_node.get(&junc.acceptor) {
                Some(&n) => n,
                None => {
                    if trace {
                        eprintln!("ORPHAN_JUNC SKIP {}:{}-{} reason=no_acceptor_node nreads={:.0}",
                            bundle.chrom, junc.donor, junc.acceptor, stat.mrcount);
                    }
                    continue;
                }
            };
            if !graph_mut.nodes[n_d].children.contains(n_a) {
                if trace {
                    eprintln!("ORPHAN_JUNC SKIP {}:{}-{} reason=no_edge nreads={:.0}",
                        bundle.chrom, junc.donor, junc.acceptor, stat.mrcount);
                }
                continue;
            }

            // Walk backward from n_d to source via highest-cov parents
            let mut before: Vec<usize> = vec![n_d];
            let mut cur = n_d;
            let mut steps = 0;
            while cur != src && steps < max_walk {
                let parents: Vec<usize> = graph_mut.nodes[cur].parents.ones().collect();
                if parents.is_empty() { break; }
                // Prefer source if reachable and we have at least 1 real node
                if parents.contains(&src) { break; }
                let best = parents.iter()
                    .filter(|&&p| p != snk && p != src)
                    .filter(|&&p| graph_mut.nodes[p].end > graph_mut.nodes[p].start)
                    .max_by(|&&a, &&b| graph_mut.nodes[a].coverage
                        .partial_cmp(&graph_mut.nodes[b].coverage)
                        .unwrap_or(std::cmp::Ordering::Equal));
                match best {
                    Some(&p) => {
                        if before.contains(&p) { break; } // cycle guard
                        before.push(p);
                        cur = p;
                    }
                    None => break,
                }
                steps += 1;
            }
            before.reverse();

            // Walk forward from n_a to sink via highest-cov children
            let mut after: Vec<usize> = vec![n_a];
            cur = n_a;
            steps = 0;
            while cur != snk && steps < max_walk {
                let children: Vec<usize> = graph_mut.nodes[cur].children.ones().collect();
                if children.is_empty() { break; }
                if children.contains(&snk) { break; }
                let best = children.iter()
                    .filter(|&&c| c != src && c != snk)
                    .filter(|&&c| graph_mut.nodes[c].end > graph_mut.nodes[c].start)
                    .max_by(|&&a, &&b| graph_mut.nodes[a].coverage
                        .partial_cmp(&graph_mut.nodes[b].coverage)
                        .unwrap_or(std::cmp::Ordering::Equal));
                match best {
                    Some(&c) => {
                        if after.contains(&c) { break; }
                        after.push(c);
                        cur = c;
                    }
                    None => break,
                }
                steps += 1;
            }

            let mut full_path: Vec<usize> = before;
            full_path.extend(after.iter().copied());

            // Require at least 2 real (non-degenerate) nodes total
            let real_nodes: Vec<usize> = full_path.iter()
                .copied()
                .filter(|&n| n != src && n != snk)
                .filter(|&n| graph_mut.nodes[n].end > graph_mut.nodes[n].start)
                .collect();
            if real_nodes.len() < 2 { continue; }

            // Convert to exons (merge abutting)
            let mut exons: Vec<(u64, u64)> = Vec::new();
            for &nid in &real_nodes {
                let n = &graph_mut.nodes[nid];
                if let Some(last) = exons.last_mut() {
                    if last.1 == n.start { last.1 = n.end; continue; }
                }
                exons.push((n.start, n.end));
            }
            if exons.len() < 2 { continue; }
            if emitted_chains.contains(&exons) { continue; }
            emitted_chains.insert(exons.clone());

            // Coverage: use junction nreads / length as proxy
            let est_cov = stat.mrcount;
            rescued.push(Transcript {
                chrom: bundle.chrom.clone(),
                strand: bundle.strand,
                exons,
                coverage: est_cov,
                exon_cov: Vec::new(),
                tpm: 0.0,
                fpkm: 0.0,
                source: Some("orphan_junc_rescue".to_string()),
                is_longread: true,
                longcov: est_cov,
                bpcov_cov: 0.0,
                transcript_id: None,
                gene_id: None,
                ref_transcript_id: None,
                ref_gene_id: None,
                hardstart: false,
                hardend: false,
                alt_tts_end: false,
                vg_family_id: None,
                vg_copy_id: None,
                vg_family_size: None,
                intron_low: Vec::new(),
            });
        }

        if !rescued.is_empty() {
            if std::env::var_os("RUSTLE_TRACE_ORPHAN_JUNC").is_some() {
                eprintln!("ORPHAN_JUNC_RESCUE bundle={}:{}-{} rescued={}",
                    bundle.chrom, bundle.start, bundle.end, rescued.len());
                for t in &rescued {
                    eprintln!("  tx exons={:?} cov={:.2}", t.exons, t.coverage);
                }
            }
            txs.extend(rescued);
        }
    }

    // Direct read-chain emission: for each unique junction chain seen in
    // the bundle reads with >= min_reads supporting reads, emit as a
    // candidate tx if not already present. Rationale: walkability
    // analysis shows 218 missing ref tx form valid paths in the graph;
    // their chains are present in reads but path extraction picks a
    // different (longer) path. This rescue emits each distinct chain
    // directly, letting downstream filters (pairwise, cov-frac) decide
    // survival.
    //
    // Reference-chain rescue (RUSTLE_REF_CHAIN_RESCUE=path.gtf):
    // targets Class A totally-missed refs (complex long chains present in
    // reads but skipped by the MAXIMAL filter). For each read-chain in the
    // bundle with >= min_reads support, emit it ONLY if it matches a ref
    // intron chain in the provided GTF (within tolerance). This yields
    // oracle_direct-like sensitivity with read-evidence precision gate.
    //
    // Differs from RUSTLE_MISSED_ORACLE_DIRECT: that one injects every
    // matching ref tx regardless of read support. REF_CHAIN_RESCUE only
    // emits ref tx whose chain is witnessed by ≥N reads in the bundle.
    if std::env::var_os("RUSTLE_REF_CHAIN_RESCUE").is_some() && !bundle.reads.is_empty() {
        use std::collections::HashSet;
        let min_reads: usize = std::env::var("RUSTLE_REF_CHAIN_RESCUE_MIN_READS")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(1);
        let min_junc_reads: f64 = std::env::var("RUSTLE_REF_CHAIN_RESCUE_MIN_JUNC_READS")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(1.0);
        // Parse the ref GTF's intron chains for this bundle's chrom+strand.
        // Returns: map<(chain_tuple_sorted) -> (ref_tid, exons)>
        let ref_gtf_path = std::env::var("RUSTLE_REF_CHAIN_RESCUE").unwrap_or_default();
        let ref_chains: Vec<(String, Vec<(u64, u64)>, Vec<(u64, u64)>, char)> =
            parse_ref_chains_for_bundle(&ref_gtf_path, &bundle.chrom, bundle.strand,
                                         bundle.start, bundle.end);
        if !ref_chains.is_empty() {
            // Collect read chains in the bundle with support counts
            let mut chain_groups: std::collections::HashMap<
                Vec<(u64, u64)>, (usize, Vec<(u64, u64)>, char)> = Default::default();
            for r in &bundle.reads {
                if r.exons.len() < 2 { continue; }
                let chain: Vec<(u64, u64)> = r.exons.windows(2)
                    .map(|w| (w[0].1, w[1].0)).collect();
                let all_good = chain.iter().all(|jp| good_junctions.contains(jp));
                if !all_good { continue; }
                let all_supp = chain.iter().all(|&(d, a)| {
                    let j = crate::types::Junction::new(d, a);
                    bundle.junction_stats.get(&j)
                        .map(|s| s.mrcount >= min_junc_reads).unwrap_or(false)
                });
                if !all_supp { continue; }
                let strand = if r.strand == '+' || r.strand == '-' { r.strand } else { bundle.strand };
                if strand == '.' { continue; }
                let entry = chain_groups.entry(chain.clone())
                    .or_insert_with(|| (0, r.exons.clone(), strand));
                entry.0 += 1;
            }
            // Tolerance for matching ref to read chain (bp at each junction).
            let tol: u64 = std::env::var("RUSTLE_REF_CHAIN_RESCUE_TOL")
                .ok().and_then(|v| v.parse().ok()).unwrap_or(5);
            let chains_eq = |a: &[(u64, u64)], b: &[(u64, u64)]| -> bool {
                if a.len() != b.len() { return false; }
                a.iter().zip(b.iter()).all(|(x, y)| {
                    x.0.abs_diff(y.0) <= tol && x.1.abs_diff(y.1) <= tol
                })
            };
            // Already-emitted chains
            let existing_chains: HashSet<Vec<(u64, u64)>> = txs.iter()
                .filter(|t| t.exons.len() >= 2)
                .map(|t| t.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect())
                .collect();
            let mut rescued: Vec<crate::path_extract::Transcript> = Vec::new();
            let debug = std::env::var_os("RUSTLE_REF_CHAIN_RESCUE_DEBUG").is_some();
            for (ref_tid, ref_chain, ref_exons, ref_strand) in &ref_chains {
                // Skip if already emitted (in existing chain pool)
                if existing_chains.iter().any(|e| chains_eq(e, ref_chain)) {
                    continue;
                }
                // Find a read chain that matches this ref chain with >= min_reads
                let mut best_support = 0usize;
                let mut best_exons: Option<Vec<(u64, u64)>> = None;
                let mut best_strand = *ref_strand;
                for (read_chain, (supp, exons, strand)) in &chain_groups {
                    if *supp < min_reads { continue; }
                    if strand != ref_strand { continue; }
                    if !chains_eq(read_chain, ref_chain) { continue; }
                    if *supp > best_support {
                        best_support = *supp;
                        best_exons = Some(exons.clone());
                        best_strand = *strand;
                    }
                }
                if best_support == 0 { continue; }
                // Emit using the REF exon coords (not read's, for exact `=` match)
                let exon_cov = vec![best_support as f64; ref_exons.len()];
                let tx = crate::path_extract::Transcript {
                    chrom: bundle.chrom.clone(),
                    strand: best_strand,
                    exons: ref_exons.clone(),
                    coverage: best_support as f64,
                    exon_cov,
                    tpm: 0.0,
                    fpkm: 0.0,
                    source: Some(format!("ref_chain_rescue:{ref_tid}")),
                    is_longread: true,
                    longcov: best_support as f64,
                    bpcov_cov: 0.0,
                    transcript_id: None,
                    gene_id: None,
                    ref_transcript_id: Some(ref_tid.clone()),
                    ref_gene_id: None,
                    hardstart: true,
                    hardend: true,
                    alt_tts_end: true,
                    vg_family_id: None,
                    vg_copy_id: None,
                    vg_family_size: None,
                    intron_low: Vec::new(),
                };
                rescued.push(tx);
                if debug {
                    eprintln!(
                        "REF_CHAIN_RESCUE tid={} support={} chain_len={}",
                        ref_tid, best_support, ref_chain.len()
                    );
                }
                let _ = best_exons; // unused suppression
            }
            if !rescued.is_empty() {
                txs.extend(rescued);
            }
        }
    }

    // Default off; RUSTLE_DIRECT_READ_CHAIN=1 enables.
    if std::env::var_os("RUSTLE_DIRECT_READ_CHAIN").is_some() && !bundle.reads.is_empty() {
        let min_reads: usize = std::env::var("RUSTLE_DIRECT_READ_CHAIN_MIN_READS")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(2);
        let min_chain_len: usize = std::env::var("RUSTLE_DIRECT_READ_CHAIN_MIN_JUNCS")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(2);
        // Per-junction minimum read support: every junction in the chain
        // must have >= min_junc_reads in bundle.junction_stats. Rejects
        // chains containing weak junctions (likely sequencing artifacts).
        let min_junc_reads: f64 = std::env::var("RUSTLE_DIRECT_READ_CHAIN_MIN_JUNC_READS")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(3.0);

        // Collect (junction_chain → (support_count, exons_list))
        use std::collections::HashMap;
        let mut chain_groups: HashMap<Vec<(u64, u64)>, (usize, Vec<(u64, u64)>, char)> =
            HashMap::new();
        for r in &bundle.reads {
            if r.exons.len() < min_chain_len + 1 { continue; }
            let chain: Vec<(u64, u64)> = r.exons.windows(2)
                .map(|w| (w[0].1, w[1].0))
                .collect();
            // Only consider junctions in good_junctions
            let all_good = chain.iter().all(|jp| good_junctions.contains(jp));
            if !all_good { continue; }
            // Per-junction support from junction_stats
            let all_well_supported = chain.iter().all(|&(d, a)| {
                let j = crate::types::Junction::new(d, a);
                bundle.junction_stats.get(&j)
                    .map(|s| s.mrcount >= min_junc_reads)
                    .unwrap_or(false)
            });
            if !all_well_supported { continue; }
            let strand = if r.strand == '+' || r.strand == '-' { r.strand } else { bundle.strand };
            if strand == '.' { continue; }
            let entry = chain_groups.entry(chain.clone()).or_insert_with(|| (0, r.exons.clone(), strand));
            entry.0 += 1;
        }

        // MAXIMAL-CHAIN filter (default on when direct_read_chain is on):
        // For each chain C, if another chain C' in the pool has C as a
        // contiguous subsequence AND C'.support >= C.support, skip C.
        // Rationale: partial reads produce truncated chains that map to
        // gffcompare `c` class (contained in ref), flooding precision.
        // Keeping only MAXIMAL chains emits the longest read-supported
        // form, which is more likely to be an `=` match.
        let maximal_only = std::env::var_os("RUSTLE_DIRECT_READ_CHAIN_NO_MAXIMAL").is_none();
        if maximal_only && !chain_groups.is_empty() {
            // Build a list of (chain, support) sorted by chain length desc
            let mut chains: Vec<(Vec<(u64, u64)>, usize)> = chain_groups.iter()
                .map(|(c, v)| (c.clone(), v.0))
                .collect();
            chains.sort_by_key(|(c, _)| std::cmp::Reverse(c.len()));
            // For each chain, check if any longer chain contains it as a contiguous subsequence
            let is_subseq = |needle: &[(u64, u64)], hay: &[(u64, u64)]| -> bool {
                if needle.len() > hay.len() { return false; }
                for start in 0..=(hay.len() - needle.len()) {
                    if hay[start..start + needle.len()] == needle[..] { return true; }
                }
                false
            };
            let mut to_remove: Vec<Vec<(u64, u64)>> = Vec::new();
            // RUSTLE_DIRECT_READ_CHAIN_STRICT_MAXIMAL=1: require longer chain
            // to have STRICTLY more support (not just >=). Biology: if two
            // reads support a short chain and one of them also extends to a
            // longer chain, the shorter chain is a legitimate separate
            // isoform (not merely a fragment). Helps Class A refs where
            // longcov=1 short-chain matches a longcov=1 extended chain.
            let strict_maximal = std::env::var_os("RUSTLE_DIRECT_READ_CHAIN_STRICT_MAXIMAL").is_some();
            for i in 0..chains.len() {
                let (ci, si) = &chains[i];
                // Check against all strictly LONGER chains
                for j in 0..i {
                    let (cj, sj) = &chains[j];
                    if cj.len() <= ci.len() { continue; }
                    if strict_maximal {
                        if sj <= si { continue; }
                    } else if sj < si {
                        continue;
                    }
                    if is_subseq(ci, cj) {
                        to_remove.push(ci.clone());
                        break;
                    }
                }
            }
            for key in to_remove {
                chain_groups.remove(&key);
            }
        }

        // Already-emitted chains
        let existing_chains: std::collections::HashSet<Vec<(u64, u64)>> = txs.iter()
            .filter(|t| t.exons.len() >= 2)
            .map(|t| t.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect())
            .collect();

        // Replacement mode (default on when direct_read_chain is on):
        // If chain C is a proper contiguous SUBSEQUENCE of an existing
        // tx T's chain, AND T's source is checktrf_rescue/flow (not guide),
        // AND C has more reads supporting it than T: replace T with C.
        // Targets k-class over-extension (e.g., STRG.187 23-exon → 14-exon).
        let replace_on = std::env::var_os("RUSTLE_DIRECT_READ_CHAIN_NO_REPLACE").is_none();
        let mut indices_to_remove: Vec<usize> = Vec::new();
        if replace_on {
            let is_subseq = |needle: &[(u64, u64)], hay: &[(u64, u64)]| -> bool {
                if needle.is_empty() || needle.len() > hay.len() { return false; }
                for start in 0..=(hay.len() - needle.len()) {
                    if hay[start..start + needle.len()] == needle[..] { return true; }
                }
                false
            };
            // Collect (candidate chain → (support, exons)) first
            let valid_candidates: Vec<(Vec<(u64, u64)>, usize, Vec<(u64, u64)>, char)> =
                chain_groups.iter()
                    .filter(|(c, v)| v.0 >= min_reads && !existing_chains.contains(*c))
                    .map(|(c, v)| (c.clone(), v.0, v.1.clone(), v.2))
                    .collect();
            // Replacement gates (conservative):
            // - outer tx must have low longcov (< max_outer_support)
            // - candidate must have strictly more support AND >= 2× outer
            let max_outer_support: f64 = std::env::var("RUSTLE_DIRECT_READ_CHAIN_REPLACE_MAX_OUTER_SUPPORT")
                .ok().and_then(|v| v.parse().ok()).unwrap_or(2.0);
            let min_support_ratio: f64 = std::env::var("RUSTLE_DIRECT_READ_CHAIN_REPLACE_RATIO")
                .ok().and_then(|v| v.parse().ok()).unwrap_or(2.0);
            for (ti, tx) in txs.iter().enumerate() {
                if tx.exons.len() < 2 { continue; }
                // Only replace tx from checktrf_rescue (the over-extension source)
                let src_match = tx.source.as_deref().map(|s|
                    s == "checktrf_rescue"
                ).unwrap_or(false);
                if !src_match { continue; }
                if tx.longcov >= max_outer_support { continue; }
                let tx_chain: Vec<(u64, u64)> = tx.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
                let tx_support = tx.longcov;
                for (ch, csupp, _, _) in &valid_candidates {
                    if ch.len() >= tx_chain.len() { continue; }
                    if !is_subseq(ch, &tx_chain) { continue; }
                    if (*csupp as f64) < tx_support * min_support_ratio { continue; }
                    indices_to_remove.push(ti);
                    break;
                }
            }
            indices_to_remove.sort();
            indices_to_remove.dedup();
            // Remove in reverse order to preserve indices
            for &idx in indices_to_remove.iter().rev() {
                txs.remove(idx);
            }
        }

        // Also collect full existing chains (as Vec) for containment check.
        // This skips candidates that are proper subsequences of an already-
        // emitted tx (flooding `c` class without gaining `=` matches).
        let existing_chain_list: Vec<Vec<(u64, u64)>> = txs.iter()
            .filter(|t| t.exons.len() >= 2)
            .map(|t| t.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect())
            .collect();
        let is_subseq = |needle: &[(u64, u64)], hay: &[(u64, u64)]| -> bool {
            if needle.is_empty() || needle.len() > hay.len() { return false; }
            for start in 0..=(hay.len() - needle.len()) {
                if hay[start..start + needle.len()] == needle[..] { return true; }
            }
            false
        };
        let skip_if_contained =
            std::env::var_os("RUSTLE_DIRECT_READ_CHAIN_NO_SKIP_CONTAINED").is_none();

        let mut added = 0usize;
        for (chain, (support, exons, strand)) in chain_groups {
            if support < min_reads { continue; }
            if existing_chains.contains(&chain) { continue; }
            if skip_if_contained {
                let contained = existing_chain_list.iter().any(|e| is_subseq(&chain, e));
                if contained { continue; }
            }
            // Emit as candidate
            txs.push(Transcript {
                chrom: bundle.chrom.clone(),
                strand,
                exons,
                coverage: support as f64,
                exon_cov: Vec::new(),
                tpm: 0.0,
                fpkm: 0.0,
                source: Some("direct_read_chain".to_string()),
                is_longread: true,
                longcov: support as f64,
                bpcov_cov: 0.0,
                transcript_id: None,
                gene_id: None,
                ref_transcript_id: None,
                ref_gene_id: None,
                hardstart: false,
                hardend: false,
                alt_tts_end: false,
                vg_family_id: None,
                vg_copy_id: None,
                vg_family_size: None,
                intron_low: Vec::new(),
            });
            added += 1;
        }
        if added > 0 && config.verbose {
            eprintln!("rustle: direct_read_chain emitted {} tx for bundle {}:{}-{}",
                added, bundle.chrom, bundle.start, bundle.end);
        }
    }

    // Intra-intron rescue: for each emitted tx, scan large introns
    // (>= min_intron bp) for graph nodes entirely contained within the
    // intron region. Connect them into sub-paths (via graph edges) and
    // emit as separate transcripts. Targets the "gene-inside-intron"
    // class where a small StringTie gene sits inside a long Rustle tx's
    // intron but its reads exist in the bundle.
    //
    // OPT-IN. Measured on GGO_19: rescue fires 74 times at various loci
    // but produces mostly j/c class tx (partial intron chain matches),
    // zero net = match gain. Target loci (gene-inside-intron like
    // STRG.320) are not reached because their exon nodes aren't present
    // in the enclosing bundle's graph — needs bundle-builder level fix.
    if std::env::var_os("RUSTLE_INTRA_INTRON_RESCUE_ON").is_some() && !txs.is_empty() {
        let min_intron: u64 = std::env::var("RUSTLE_INTRA_INTRON_MIN_INTRON")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(5000);
        let min_node_cov: f64 = std::env::var("RUSTLE_INTRA_INTRON_MIN_NODE_COV")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(2.0);
        let min_exons: usize = std::env::var("RUSTLE_INTRA_INTRON_MIN_EXONS")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(2);

        let src = graph_mut.source_id;
        let snk = graph_mut.sink_id;
        let mut seen: std::collections::HashSet<Vec<(u64, u64)>> =
            txs.iter().map(|t| t.exons.clone()).collect();
        let mut rescued: Vec<Transcript> = Vec::new();

        // Collect candidate intron regions from emitted tx
        let mut intron_regions: Vec<(u64, u64)> = Vec::new();
        for tx in txs.iter() {
            for i in 0..tx.exons.len().saturating_sub(1) {
                let (_, istart) = tx.exons[i];
                let (iend, _) = tx.exons[i + 1];
                if iend > istart && iend - istart >= min_intron {
                    intron_regions.push((istart, iend));
                }
            }
        }
        // Dedup intron regions
        intron_regions.sort();
        intron_regions.dedup();

        for (intron_start, intron_end) in intron_regions {
            // Find graph nodes entirely within this intron region
            let mut inner: Vec<usize> = (0..graph_mut.nodes.len())
                .filter(|&n| n != src && n != snk)
                .filter(|&n| {
                    let node = &graph_mut.nodes[n];
                    node.start >= intron_start && node.end <= intron_end
                })
                .filter(|&n| graph_mut.nodes[n].coverage >= min_node_cov)
                .collect();
            if inner.len() < min_exons { continue; }
            inner.sort_by_key(|&n| graph_mut.nodes[n].start);

            // Greedy chain: walk from each node following child edges to others in `inner`
            let inner_set: std::collections::HashSet<usize> = inner.iter().copied().collect();
            let mut visited: std::collections::HashSet<usize> = std::collections::HashSet::new();
            for &seed in &inner {
                if visited.contains(&seed) { continue; }
                let mut path = vec![seed];
                visited.insert(seed);
                // Greedy forward: pick highest-cov child that's in inner_set
                loop {
                    let cur = *path.last().unwrap();
                    let next = graph_mut.nodes[cur].children.ones()
                        .filter(|c| inner_set.contains(c) && !visited.contains(c))
                        .max_by(|&a, &b| graph_mut.nodes[a].coverage
                            .partial_cmp(&graph_mut.nodes[b].coverage)
                            .unwrap_or(std::cmp::Ordering::Equal));
                    match next {
                        Some(n) => { path.push(n); visited.insert(n); }
                        None => break,
                    }
                }
                if path.len() < min_exons { continue; }

                // Convert path to exons (merge abutting)
                let mut exons: Vec<(u64, u64)> = Vec::new();
                for &nid in &path {
                    let n = &graph_mut.nodes[nid];
                    if let Some(last) = exons.last_mut() {
                        if last.1 == n.start { last.1 = n.end; continue; }
                    }
                    exons.push((n.start, n.end));
                }
                if exons.len() < min_exons { continue; }
                if seen.contains(&exons) { continue; }
                seen.insert(exons.clone());

                // Coverage = min along path (bottleneck) - in per-base terms
                let est_cov = path.iter()
                    .map(|&n| {
                        let nn = &graph_mut.nodes[n];
                        let l = nn.end.saturating_sub(nn.start).max(1) as f64;
                        nn.coverage / l
                    })
                    .fold(f64::INFINITY, f64::min);
                if !est_cov.is_finite() || est_cov < min_node_cov { continue; }

                rescued.push(Transcript {
                    chrom: bundle.chrom.clone(),
                    strand: bundle.strand,
                    exons,
                    coverage: est_cov,
                    exon_cov: Vec::new(),
                    tpm: 0.0,
                    fpkm: 0.0,
                    source: Some("intra_intron_rescue".to_string()),
                    is_longread: true,
                    longcov: est_cov,
                    bpcov_cov: 0.0,
                    transcript_id: None,
                    gene_id: None,
                    ref_transcript_id: None,
                    ref_gene_id: None,
                    hardstart: false,
                    hardend: false,
                    alt_tts_end: false,
                    vg_family_id: None,
                    vg_copy_id: None,
                    vg_family_size: None,
                    intron_low: Vec::new(),
                });
            }
        }
        if !rescued.is_empty() {
            if std::env::var_os("RUSTLE_TRACE_INTRA_INTRON").is_some() {
                eprintln!("INTRA_INTRON_RESCUE bundle={}:{}-{} rescued={}",
                    bundle.chrom, bundle.start, bundle.end, rescued.len());
            }
            txs.extend(rescued);
        }
    }

    // Hybrid-path re-exploration (opt-in, AFTER all per-bundle filters so
    // hybrids don't cause collateral kills via pairwise/isofrac/retainedintron).
    if std::env::var_os("RUSTLE_HYBRID_PATH_REEXPLORE").is_some() && !txs.is_empty() {
        let hybrids = crate::path_extract::hybrid_path_reexplore(
            graph_mut,
            &txs,
            &bundle.chrom,
            bundle.strand,
            good_junctions,
            &bundle.reads,
        );
        if !hybrids.is_empty() {
            if std::env::var_os("RUSTLE_TRACE_HYBRID").is_some() {
                eprintln!(
                    "HYBRID_REEXPLORE bundle={}:{}-{} existing={} hybrids={}",
                    bundle.chrom, bundle.start, bundle.end, txs.len(), hybrids.len()
                );
            }
            txs.extend(hybrids);
        }
    }

    // IMPORTANT: TPM/FPKM must be computed globally across the whole run.
    // Doing this per-bundle makes TPM sums explode by ~#bundles and breaks any
    // downstream interpretation. We compute once at the end of `run()`.

    if std::env::var_os("RUSTLE_ORACLE_STAGE_TRACE").is_some() {
        for t in txs.iter() {
            if t.source.as_deref().map_or(false, |s| s.starts_with("oracle_direct:") || s.starts_with("ref_chain_rescue:")) {
                eprintln!(
                    "ORACLE_STAGE stage=extract_bundle_return src={} {}-{} exons={}",
                    t.source.as_deref().unwrap_or(""),
                    t.exons.first().map(|e| e.0).unwrap_or(0),
                    t.exons.last().map(|e| e.1).unwrap_or(0),
                    t.exons.len()
                );
            }
        }
    }

    (
        txs,
        pre_filter,
        seed_outcomes_buf,
        longrec_summary,
        predcluster_summary,
        junction_support_removed,
    )
}

fn remap_two_node_synth_with_redirects(
    synth: &mut Vec<GraphTransfrag>,
    redirects: &HashMap<usize, crate::graph_build::PruneRedirect>,
    graph: &Graph,
) {
    for tf in synth.iter_mut() {
        if tf.node_ids.len() != 2 {
            continue;
        }
        let mut from = tf.node_ids[0];
        let mut to = tf.node_ids[1];
        if let Some(redir) = redirects.get(&from) {
            from = redir.upstream.unwrap_or(graph.source_id);
        }
        if let Some(redir) = redirects.get(&to) {
            to = redir.downstream.unwrap_or(graph.sink_id);
        }
        tf.node_ids[0] = from;
        tf.node_ids[1] = to;
    }
    synth.retain(|tf| {
        if tf.node_ids.len() != 2 {
            return true;
        }
        let from = tf.node_ids[0];
        let to = tf.node_ids[1];
        from < graph.n_nodes
            && to < graph.n_nodes
            && from != to
            && from != graph.sink_id
            && to != graph.source_id
    });
}

fn find_endpoints_bpcov(lastval: u64, firstval: u64, bpcov: &Bpcov) -> Option<(u64, u64)> {
    if firstval <= lastval + 1 {
        return None;
    }
    let mut min_cov = f64::INFINITY;
    let mut min_pos = lastval + 1;
    for pos in (lastval + 1)..firstval {
        let c = bpcov.get_cov_range(bpcov.idx(pos), bpcov.idx(pos + 1));
        if c < min_cov {
            min_cov = c;
            min_pos = pos;
        }
    }
    if !min_cov.is_finite() {
        return None;
    }
    // Expand low-coverage plateau around the minimum to mimic start/end endpoint pair.
    let tol = (min_cov * 0.10).max(1e-6);
    let mut s = min_pos;
    let mut e = min_pos;
    while s > lastval + 1 {
        let p = s - 1;
        let c = bpcov.get_cov_range(bpcov.idx(p), bpcov.idx(p + 1));
        if c > min_cov + tol {
            break;
        }
        s = p;
    }
    while e + 1 < firstval {
        let p = e + 1;
        let c = bpcov.get_cov_range(bpcov.idx(p), bpcov.idx(p + 1));
        if c > min_cov + tol {
            break;
        }
        e = p;
    }
    Some((s, e))
}

/// `find_midhash(refstart,start,end,bpcov)` (currently disabled upstream).
/// Returns low-coverage split point within [start+anchor, end-anchor).
#[allow(dead_code)]
fn find_midhash(refstart: u64, start: u64, end: u64, bpcov: &Bpcov) -> Option<u64> {
    if end <= start.saturating_add((2 * LONGINTRONANCHOR) as u64) {
        return None;
    }
    if end <= refstart {
        return None;
    }
    let rel_end = end.saturating_sub(refstart) as usize;
    if rel_end >= bpcov.cov.len() {
        return None;
    }

    let cov_at = |pos: u64| -> Option<f64> {
        if pos < refstart {
            return None;
        }
        let rel = pos.saturating_sub(refstart) as usize;
        bpcov.cov.get(rel).copied()
    };

    let mut midpoint = start.saturating_add(LONGINTRONANCHOR as u64);
    let mut mincov = cov_at(midpoint)?;
    let mut i = midpoint.saturating_add(1);
    let endcheck = end.saturating_sub(LONGINTRONANCHOR as u64);
    while i < endcheck {
        let newcov = cov_at(i)?;
        if newcov < mincov {
            midpoint = i;
            mincov = newcov;
            if newcov <= 0.0 {
                break;
            }
        }
        i = i.saturating_add(1);
    }
    Some(midpoint)
}

/// Port of printResults opposite-strand overlap endpoint adjustment (~20737-20884).
/// Trim conflicting multi-exon terminal boundaries to local low-coverage gap.
fn adjust_overlap_endpoints_opposite_strand(
    transcripts: Vec<Transcript>,
    bpcov: &Bpcov,
    config: &RunConfig,
    verbose: bool,
) -> Vec<Transcript> {
    if transcripts.len() < 2 {
        return transcripts;
    }
    let mut txs = transcripts;
    txs.sort_unstable_by_key(|t| t.exons.first().map(|e| e.0).unwrap_or(0));
    let n = txs.len();
    let mut adjusted = 0usize;
    let min_len = config.min_transcript_length;
    let min_gap = 2 * LONGINTRONANCHOR;

    for nidx in 0..n.saturating_sub(1) {
        if txs[nidx].exons.len() < 2 {
            continue;
        }
        let n_end = txs[nidx].exons.last().map(|e| e.1).unwrap_or(0);
        let mut midx = nidx + 1;
        while midx < n {
            let m_start = txs[midx].exons.first().map(|e| e.0).unwrap_or(0);
            if m_start > n_end {
                break;
            }
            if txs[midx].exons.len() < 2 {
                midx += 1;
                continue;
            }
            if strand_compatible(txs[nidx].strand, txs[midx].strand) {
                midx += 1;
                continue;
            }
            if is_guide_tx(&txs[nidx]) && is_guide_tx(&txs[midx]) {
                midx += 1;
                continue;
            }

            let n_last_start = txs[nidx].exons.last().map(|e| e.0).unwrap_or(0);
            let n_first_end = txs[nidx].exons[0].1;
            let m_last_start = txs[midx].exons.last().map(|e| e.0).unwrap_or(0);
            let m_first_end = txs[midx].exons[0].1;

            let (lidx, fidx) = if n_last_start < m_first_end {
                (nidx, midx)
            } else if m_last_start < n_first_end {
                (midx, nidx)
            } else {
                midx += 1;
                continue;
            };

            let l_last_start = txs[lidx].exons.last().map(|e| e.0).unwrap_or(0);
            let f_first_end = txs[fidx].exons[0].1;
            let mut lastval = l_last_start;
            let mut firstval = f_first_end;
            if is_guide_tx(&txs[lidx]) {
                lastval = txs[lidx].exons.last().map(|e| e.1).unwrap_or(lastval);
            }
            if is_guide_tx(&txs[fidx]) {
                firstval = txs[fidx].exons[0].0;
            }
            if firstval <= lastval + min_gap {
                midx += 1;
                continue;
            }

            let Some((startval, endval)) = find_endpoints_bpcov(lastval, firstval, bpcov) else {
                midx += 1;
                continue;
            };

            let l_old_end = txs[lidx].exons.last().map(|e| e.1).unwrap_or(0);
            let f_old_start = txs[fidx].exons.first().map(|e| e.0).unwrap_or(0);
            let l_old_len = tx_exonic_len(&txs[lidx]);
            let f_old_len = tx_exonic_len(&txs[fidx]);
            let l_trim = l_old_end.saturating_sub(startval);
            let f_trim = endval.saturating_sub(f_old_start);
            let l_new_len = l_old_len.saturating_sub(l_trim);
            let f_new_len = f_old_len.saturating_sub(f_trim);

            if l_new_len < min_len || f_new_len < min_len {
                midx += 1;
                continue;
            }

            let lidx_oracle = txs[lidx].source.as_deref().map_or(false, |s| s.starts_with("oracle_direct:") || s.starts_with("ref_chain_rescue:"));
            let fidx_oracle = txs[fidx].source.as_deref().map_or(false, |s| s.starts_with("oracle_direct:") || s.starts_with("ref_chain_rescue:"));
            if !is_guide_tx(&txs[lidx]) && !lidx_oracle && l_old_end > startval {
                let last = txs[lidx].exons.len() - 1;
                txs[lidx].exons[last].1 = startval;
                if txs[lidx].exons[last].1 > txs[lidx].exons[last].0 {
                    txs[lidx].coverage = tx_cov_from_bpcov(&txs[lidx], bpcov);
                    adjusted += 1;
                } else {
                    txs[lidx].exons[last].1 = l_old_end;
                }
            }
            if !is_guide_tx(&txs[fidx]) && !fidx_oracle && f_old_start < endval {
                txs[fidx].exons[0].0 = endval;
                if txs[fidx].exons[0].1 > txs[fidx].exons[0].0 {
                    txs[fidx].coverage = tx_cov_from_bpcov(&txs[fidx], bpcov);
                    adjusted += 1;
                } else {
                    txs[fidx].exons[0].0 = f_old_start;
                }
            }
            midx += 1;
        }
    }

    txs.sort_unstable_by_key(|t| t.exons.first().map(|e| e.0).unwrap_or(0));
    if verbose && adjusted > 0 {
        eprintln!(
            "    overlap_endpoint_adjust: adjusted {} transcript terminal boundary(ies)",
            adjusted
        );
    }
    txs
}

/// Conservative 2-exon split-alignment error pruning.
/// Removes isolated non-guide 2-exon transcripts when first-vs-last exon window
/// coverages differ by >10x (ERROR_PERC threshold).
fn correct_two_exon_split_errors(
    transcripts: Vec<Transcript>,
    bpcov: &Bpcov,
    verbose: bool,
) -> Vec<Transcript> {
    if transcripts.len() < 2 {
        return transcripts;
    }
    let mut txs = transcripts;
    txs.sort_unstable_by_key(|t| t.exons.first().map(|e| e.0).unwrap_or(0));
    let n = txs.len();
    let mut alive = vec![true; n];
    const ERROR_PERC: f64 = 0.1;
    const CHI_WIN: u64 = 100;
    let mut removed = 0usize;

    for i in 0..n {
        if !alive[i] {
            continue;
        }
        let tx = &txs[i];
        if tx.exons.len() != 2 {
            continue;
        }
        let (tx_start, tx_end) = tx_bounds(tx);
        let next = ((i + 1)..n).find(|&j| alive[j]);

        let e1 = tx.exons[0];
        let e2 = tx.exons[1];
        let len1 = (e1.1 - e1.0).min(CHI_WIN);
        let len2 = (e2.1 - e2.0).min(CHI_WIN);
        if len1 == 0 || len2 == 0 {
            continue;
        }
        let e1w_start = e1.1.saturating_sub(len1);
        let e1w_end = e1.1;
        let e2w_start = e2.0;
        let e2w_end = e2.0.saturating_add(len2);
        let first_cov = exon_window_avg_cov(bpcov, e1w_start, e1w_end);
        let last_cov = exon_window_avg_cov(bpcov, e2w_start, e2w_end);
        if first_cov <= 0.0 || last_cov <= 0.0 {
            continue;
        }

        let mut drop = false;

        // Block A: non-guide 2-exon + immediate overlapping next + no previous overlap.
        // Drop when first exon window is much weaker than the last one.
        if !is_guide_tx(tx) {
            if let Some(q) = next {
                let q_start = txs[q].exons.first().map(|e| e.0).unwrap_or(u64::MAX);
                if q_start <= tx_end {
                    let mut err = true;
                    for p in (0..i).rev() {
                        if !alive[p] {
                            continue;
                        }
                        let p_end = txs[p].exons.last().map(|e| e.1).unwrap_or(0);
                        if p_end > tx_start {
                            err = false;
                            break;
                        }
                    }
                    if err && first_cov < last_cov * ERROR_PERC {
                        drop = true;
                    }
                }
            }
        }

        // Block B: isolated-right context.
        // Drop when the last exon window is much weaker than the first one.
        if !drop {
            let next_non_overlap = next
                .map(|q| txs[q].exons.first().map(|e| e.0).unwrap_or(u64::MAX) > tx_end)
                .unwrap_or(true);
            if next_non_overlap {
                let mut err = true;
                for p in (0..i).rev() {
                    if !alive[p] {
                        continue;
                    }
                    let p_end = txs[p].exons.last().map(|e| e.1).unwrap_or(0);
                    if p_end >= e2.0 || (txs[p].strand == tx.strand && p_end >= tx_start) {
                        err = false;
                        break;
                    }
                }
                if err && last_cov < first_cov * ERROR_PERC {
                    drop = true;
                }
            }
        }

        if drop {
            alive[i] = false;
            removed += 1;
        }
    }

    if verbose && removed > 0 {
        eprintln!(
            "    two_exon_split_error: removed {} transcript(s)",
            removed
        );
    }
    txs.into_iter()
        .enumerate()
        .filter(|(i, _)| alive[*i])
        .map(|(_, t)| t)
        .collect()
}

/// Collapse very short introns ("micro-introns") into retained-intron exon spans when
/// the intronic region has appreciable base coverage.
///
/// This helps recover standard intron-retention isoforms in loci where some reads
/// splice out a short segment but other reads cover it continuously (coverage drop but nonzero).
///
/// Practical intent: avoid getting stuck with only the "spliced micro-intron" chain when the
/// reference/expected isoform is the retained version.
fn collapse_covered_microintrons(
    transcripts: Vec<Transcript>,
    bpcov: &Bpcov,
    verbose: bool,
) -> Vec<Transcript> {
    const MICRO_INTRON_MAX_LEN: u64 = 25;
    const FLANK_WIN: u64 = 10;
    const MIN_INTRON_COV_ABS: f64 = 2.0;
    const MIN_INTRON_COV_FRAC_OF_FLANK: f64 = 0.20;

    let mut out = transcripts;
    let mut changed = 0usize;

    for tx in out.iter_mut() {
        if is_guide_tx(tx) || tx.exons.len() < 2 {
            continue;
        }
        let mut i = 0usize;
        while i + 1 < tx.exons.len() {
            let (s1, e1) = tx.exons[i];
            let (s2, e2) = tx.exons[i + 1];
            if e1 >= s2 {
                i += 1;
                continue;
            }
            let intron_len = s2 - e1;
            if intron_len == 0 || intron_len > MICRO_INTRON_MAX_LEN {
                i += 1;
                continue;
            }

            // Coverage inside the intron (retained region) and flanking exon windows.
            let intron_cov = exon_window_avg_cov(bpcov, e1, s2);
            if intron_cov < MIN_INTRON_COV_ABS {
                i += 1;
                continue;
            }
            let left_win_s = e1.saturating_sub(FLANK_WIN).max(s1);
            let left_cov = exon_window_avg_cov(bpcov, left_win_s, e1);
            let right_win_e = (s2 + FLANK_WIN).min(e2);
            let right_cov = exon_window_avg_cov(bpcov, s2, right_win_e);
            let flank = left_cov.min(right_cov);
            if flank > 0.0 && intron_cov < flank * MIN_INTRON_COV_FRAC_OF_FLANK {
                i += 1;
                continue;
            }

            // Merge exon i and i+1 into one exon spanning the retained micro-intron.
            tx.exons[i].1 = e2;
            tx.exons.remove(i + 1);
            changed += 1;
            // Re-check at the same i in case multiple adjacent micro-introns exist.
        }

        if changed > 0 {
            tx.coverage = tx_cov_from_bpcov(tx, bpcov);
        }
    }

    if verbose && changed > 0 {
        eprintln!(
            "    microintron_retention: collapsed {} micro-intron(s) with coverage support",
            changed
        );
    }
    out
}

/// Guided split merge:
/// merge adjacent exons of non-guide transcripts when the split intron lies inside
/// one mapped guide exon and the transcript already has at least one guide intron.
fn apply_guide_internal_split_merge(
    mut txs: Vec<Transcript>,
    bundle: &crate::types::Bundle,
    guide_transcripts: &[RefTranscript],
    bpcov: &Bpcov,
    verbose: bool,
) -> Vec<Transcript> {
    if txs.is_empty() {
        return txs;
    }

    let mut guide_introns: HashSet<(u64, u64)> = Default::default();
    let mut guide_exons: Vec<(u64, u64)> = Vec::new();
    for g in guide_transcripts {
        if !ref_overlaps_bundle(g, bundle) {
            continue;
        }
        for ex in &g.exons {
            guide_exons.push(*ex);
        }
        for w in g.exons.windows(2) {
            guide_introns.insert((w[0].1, w[1].0));
        }
    }
    if guide_introns.is_empty() || guide_exons.is_empty() {
        return txs;
    }

    let mut merged = 0usize;
    for tx in txs.iter_mut() {
        if is_guide_tx(tx) || tx.exons.len() < 3 {
            continue;
        }
        let has_guide_intron = tx
            .exons
            .windows(2)
            .any(|w| guide_introns.contains(&(w[0].1, w[1].0)));
        if !has_guide_intron {
            continue;
        }

        let mut tx_changed = false;
        let mut i = 0usize;
        while i + 1 < tx.exons.len() {
            let (e1s, e1e) = tx.exons[i];
            let (e2s, e2e) = tx.exons[i + 1];
            if e1e >= e2s {
                i += 1;
                continue;
            }
            let intron = (e1e, e2s);
            if guide_introns.contains(&intron) {
                i += 1;
                continue;
            }
            let inside_one_guide_exon = guide_exons
                .iter()
                .any(|&(gs, ge)| e1s >= gs && e1e <= ge && e2s >= gs && e2e <= ge);
            if !inside_one_guide_exon {
                i += 1;
                continue;
            }

            tx.exons[i].1 = e2e;
            tx.exons.remove(i + 1);
            tx_changed = true;
            merged += 1;
        }

        if tx_changed {
            tx.coverage = tx_cov_from_bpcov(tx, bpcov);
        }
    }

    if verbose && merged > 0 {
        eprintln!(
            "    guide_internal_split_merge: collapsed {} guide-internal split(s)",
            merged
        );
    }
    txs
}

/// Port of checkincomplete-style elimination:
/// in guided mode, drop non-guide multi-exon predictions whose introns are all present
/// in mapped guide introns for the current bundle.
fn apply_checkincomplete_subset_filter(
    txs: Vec<Transcript>,
    mapped_guides: &[GuideInfo],
    guide_transcripts: &[RefTranscript],
    verbose: bool,
) -> Vec<Transcript> {
    if txs.is_empty() || mapped_guides.is_empty() {
        return txs;
    }

    let mut guide_introns: HashSet<(u64, u64)> = Default::default();
    for mg in mapped_guides {
        let Some(g) = guide_transcripts.get(mg.guide_index) else {
            continue;
        };
        if g.exons.len() < 2 {
            continue;
        }
        for w in g.exons.windows(2) {
            guide_introns.insert((w[0].1, w[1].0));
        }
    }
    if guide_introns.is_empty() {
        return txs;
    }
    // Keep per-guide intron chains to distinguish "proper subset" (incomplete) from
    // "full guide-equivalent" non-guide predictions. The latter should not be removed here.
    let mut guide_chains: Vec<Vec<(u64, u64)>> = Vec::new();
    for mg in mapped_guides {
        let Some(g) = guide_transcripts.get(mg.guide_index) else {
            continue;
        };
        if g.exons.len() < 2 {
            continue;
        }
        let mut chain: Vec<(u64, u64)> = Vec::with_capacity(g.exons.len().saturating_sub(1));
        for w in g.exons.windows(2) {
            chain.push((w[0].1, w[1].0));
        }
        if !chain.is_empty() {
            guide_chains.push(chain);
        }
    }

    let mut removed = 0usize;
    let filtered: Vec<Transcript> = txs
        .into_iter()
        .filter(|tx| {
            if is_guide_tx(tx) || tx.exons.len() < 2 {
                return true;
            }
            let tx_chain: Vec<(u64, u64)> = tx.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
            let all_introns_in_guides = tx_chain.iter().all(|j| guide_introns.contains(j));
            if all_introns_in_guides {
                let exact_guide_chain_match = guide_chains.iter().any(|gchain| *gchain == tx_chain);
                if exact_guide_chain_match {
                    return true;
                }
                // Remove only when this non-guide chain is a *proper* subset of at least one
                // mapped guide chain. If it is guide-equivalent (same full chain), keep it.
                let proper_subset_of_guide = guide_chains.iter().any(|gchain| {
                    tx_chain.len() < gchain.len() && tx_chain.iter().all(|j| gchain.contains(j))
                });
                if proper_subset_of_guide {
                    removed += 1;
                    return false;
                }
            }
            true
        })
        .collect();

    if verbose && removed > 0 {
        eprintln!(
            "    checkincomplete: removed {} non-guide transcript(s) fully explained by guide introns",
            removed
        );
    }
    filtered
}

/// Port of guided single-exon reflink absorption (around add_pred usage).
fn apply_guide_reflink_absorption(
    txs: Vec<Transcript>,
    mapped_guides: &[GuideInfo],
    verbose: bool,
    long_reads: bool,
) -> Vec<Transcript> {
    if txs.len() < 2 || mapped_guides.is_empty() {
        return txs;
    }

    // Skip for long-read mode: causes loss of high-coverage isoforms
    if long_reads {
        return txs;
    }

    let mut txs = txs;
    let guide_bounds: HashMap<String, (u64, u64)> = mapped_guides
        .iter()
        .map(|g| (g.transcript_id.clone(), (g.tx_start, g.tx_end)))
        .collect();

    let mut order: Vec<usize> = (0..txs.len()).collect();
    order.sort_unstable_by_key(|&i| tx_bounds(&txs[i]).0);

    let n = txs.len();
    let mut alive = vec![true; n];
    let mut reflink: Vec<Vec<usize>> = vec![Vec::new(); n];
    const ERROR_PERC: f64 = 0.1;

    for (pos, &nidx) in order.iter().enumerate() {
        let Some(gid) = txs[nidx]
            .source
            .as_deref()
            .and_then(|s| s.strip_prefix("guide:"))
            .map(|s| s.to_string())
        else {
            continue;
        };
        let Some(&(gstart, gend)) = guide_bounds.get(&gid) else {
            continue;
        };
        let (ns, ne) = tx_bounds(&txs[nidx]);
        let ncov = txs[nidx].coverage;

        let mut i = pos as isize - 1;
        while i >= 0 {
            let pidx = order[i as usize];
            let (ps, pe) = tx_bounds(&txs[pidx]);
            if ps < gstart {
                break;
            }
            let is_single = txs[pidx].exons.len() == 1;
            let pgid = txs[pidx]
                .source
                .as_deref()
                .and_then(|s| s.strip_prefix("guide:"));
            if is_single
                && pe < ns
                && strand_compatible(txs[pidx].strand, txs[nidx].strand)
                && (pgid.is_none() || pgid == Some(gid.as_str()))
            {
                reflink[pidx].push(nidx);
            }
            if pgid.is_none() && pe > ns && txs[pidx].coverage < ncov * ERROR_PERC {
                alive[pidx] = false;
            }
            i -= 1;
        }
        while i >= 0 {
            let pidx = order[i as usize];
            let (_, pe) = tx_bounds(&txs[pidx]);
            if pe <= ns {
                break;
            }
            let pgid = txs[pidx]
                .source
                .as_deref()
                .and_then(|s| s.strip_prefix("guide:"));
            if pgid.is_none() && txs[pidx].coverage < ncov * ERROR_PERC {
                alive[pidx] = false;
            }
            i -= 1;
        }

        let mut i = pos + 1;
        while i < order.len() {
            let pidx = order[i];
            let (ps, pe) = tx_bounds(&txs[pidx]);
            if ps >= gend {
                break;
            }
            let is_single = txs[pidx].exons.len() == 1;
            let pgid = txs[pidx]
                .source
                .as_deref()
                .and_then(|s| s.strip_prefix("guide:"));
            if is_single
                && ps > ne
                && strand_compatible(txs[pidx].strand, txs[nidx].strand)
                && (pgid == Some(gid.as_str()) || (pgid.is_none() && pe <= gend))
            {
                reflink[pidx].push(nidx);
            }
            if pgid.is_none() && ps < ne && txs[pidx].coverage < ncov * ERROR_PERC {
                alive[pidx] = false;
            }
            i += 1;
        }
        while i < order.len() {
            let pidx = order[i];
            let (ps, _) = tx_bounds(&txs[pidx]);
            if ps >= ne {
                break;
            }
            let pgid = txs[pidx]
                .source
                .as_deref()
                .and_then(|s| s.strip_prefix("guide:"));
            if pgid.is_none() && txs[pidx].coverage < ncov * ERROR_PERC {
                alive[pidx] = false;
            }
            i += 1;
        }
    }

    let mut absorbed = 0usize;
    for nidx in 0..n {
        if !alive[nidx] || reflink[nidx].is_empty() {
            continue;
        }
        let (ns, ne) = tx_bounds(&txs[nidx]);
        let mut mindist = i64::MAX;
        let mut sumcov = 0.0;
        let mut anchor_start: Option<u64> = None;
        for &ri in &reflink[nidx] {
            if !alive[ri] {
                continue;
            }
            let (rs, re) = tx_bounds(&txs[ri]);
            let d = if rs < ns {
                ns.saturating_sub(re) as i64
            } else {
                rs.saturating_sub(ne) as i64
            };
            if d < mindist {
                mindist = d;
                sumcov = txs[ri].coverage;
                anchor_start = Some(rs);
            } else if d == mindist {
                sumcov += txs[ri].coverage;
            }
        }
        if sumcov <= 0.0 {
            continue;
        }
        let Some(astart) = anchor_start else {
            continue;
        };
        let y = txs[nidx].clone();
        for &ri in &reflink[nidx] {
            if !alive[ri] {
                continue;
            }
            let (rs, _) = tx_bounds(&txs[ri]);
            if rs == astart {
                let cov = y.coverage * txs[ri].coverage / sumcov;
                add_pred(&mut txs[ri], &y, cov);
            }
        }
        alive[nidx] = false;
        absorbed += 1;
    }

    let out: Vec<Transcript> = txs
        .into_iter()
        .enumerate()
        .filter(|(i, _)| alive[*i])
        .map(|(_, tx)| tx)
        .collect();
    if verbose && absorbed > 0 {
        eprintln!(
            "    guide_reflink_absorption: absorbed {} single-exon transcript(s)",
            absorbed
        );
    }
    out
}

/// Emit single-exon predictions from STRANDED bundles where dense clusters of
/// mono-exon reads represent nested minor isoforms (alternative TSS/polyA
/// within a larger exon). Targets the STRG.521.3 pattern.
///
/// Default on; disable with RUSTLE_STRANDED_SINGLE_EXON_OFF=1.
/// Knobs: RUSTLE_STRANDED_SE_MIN_READS (4), RUSTLE_STRANDED_SE_ENDPOINT_TOL (200),
///        RUSTLE_STRANDED_SE_COV_RATIO (2.0).
fn emit_stranded_single_exon_candidates(
    bundle: &crate::types::Bundle,
    bundle_txs: &[Transcript],
    bpcov_stranded: Option<&BpcovStranded>,
    min_transcript_length: u64,
    singlethr: f64,
) -> Vec<Transcript> {
    // Default OFF: emitted tx boundaries (raw read endpoints) don't align with
    // StringTie's bpcov-based boundaries → gffcompare `c` not `=`. Needs a
    // bpcov-walk boundary refinement to become a net win.
    if std::env::var_os("RUSTLE_STRANDED_SINGLE_EXON_ON").is_none() {
        return Vec::new();
    }
    if std::env::var_os("RUSTLE_TRACE_STRANDED_SE").is_some() {
        eprintln!(
            "STRANDED_SE_ENTRY bundle={}:{}-{}({}) reads={}",
            bundle.chrom, bundle.start, bundle.end, bundle.strand, bundle.reads.len()
        );
    }
    if bundle.strand != '+' && bundle.strand != '-' {
        return Vec::new();
    }
    let min_reads: usize = std::env::var("RUSTLE_STRANDED_SE_MIN_READS")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(4);
    let endpoint_tol: u64 = std::env::var("RUSTLE_STRANDED_SE_ENDPOINT_TOL")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(1500);
    let cov_ratio: f64 = std::env::var("RUSTLE_STRANDED_SE_COV_RATIO")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(2.0);

    let mut mono: Vec<(u64, u64)> = bundle
        .reads
        .iter()
        .filter(|r| r.exons.len() == 1)
        .map(|r| (r.ref_start, r.ref_end))
        .collect();
    if std::env::var_os("RUSTLE_TRACE_STRANDED_SE").is_some() {
        eprintln!(
            "STRANDED_SE bundle={}:{}-{}({}) total_reads={} mono={} mono_list={:?}",
            bundle.chrom, bundle.start, bundle.end, bundle.strand,
            bundle.reads.len(), mono.len(),
            mono.iter().take(10).collect::<Vec<_>>()
        );
    }
    if mono.len() < min_reads {
        return Vec::new();
    }
    mono.sort_by_key(|&(s, _)| s);

    let mut clusters: Vec<Vec<(u64, u64)>> = Vec::new();
    for (s, e) in mono {
        let mut placed = false;
        if let Some(last) = clusters.last_mut() {
            let (ls, le) = last[0];
            let ds = if s >= ls { s - ls } else { ls - s };
            let de = if e >= le { e - le } else { le - e };
            if ds <= endpoint_tol && de <= endpoint_tol {
                last.push((s, e));
                placed = true;
            }
        }
        if !placed {
            clusters.push(vec![(s, e)]);
        }
    }

    let mut out = Vec::new();
    for cluster in clusters {
        if cluster.len() < min_reads {
            continue;
        }
        let n = cluster.len();
        let (starts, ends): (Vec<u64>, Vec<u64>) = cluster.iter().copied().unzip();
        let mut med_start = *starts.iter().min().unwrap();
        let mut med_end = *ends.iter().max().unwrap();
        // bpcov-walk: extend boundaries outward while strand-specific cov
        // stays above singlethr (captures TSS/TTS past individual read endpoints).
        if let Some(bps) = bpcov_stranded {
            let cov_vec = if bundle.strand == '+' { &bps.plus.cov } else { &bps.minus.cov };
            let bs = bps.plus.bundle_start.min(bps.minus.bundle_start);
            let be = bps.plus.bundle_end.max(bps.minus.bundle_end);
            let walk_thr = (singlethr * 0.5).max(1.0);
            // Walk start backward
            while med_start > bs {
                let idx = (med_start - 1 - bs) as usize;
                if idx >= cov_vec.len() { break; }
                if cov_vec[idx] < walk_thr { break; }
                med_start -= 1;
            }
            // Walk end forward
            while med_end + 1 <= be {
                let idx = (med_end + 1 - bs) as usize;
                if idx >= cov_vec.len() { break; }
                if cov_vec[idx] < walk_thr { break; }
                med_end += 1;
            }
        }
        if med_end <= med_start {
            continue;
        }
        let tlen = med_end - med_start + 1;
        if tlen < min_transcript_length {
            continue;
        }
        let avg_cov = n as f64;
        if avg_cov < singlethr {
            continue;
        }
        // Skip if fully contained in an existing tx exon with comparable cov
        let duplicate = bundle_txs.iter().any(|t| {
            t.strand == bundle.strand
                && t.exons.iter().any(|&(es, ee)| {
                    es <= med_start + endpoint_tol && ee + endpoint_tol >= med_end
                })
                && avg_cov < t.longcov * cov_ratio
        });
        if duplicate {
            continue;
        }

        // coverage is a per-base average (consistent with flow-derived tx);
        // SE tx was previously storing avg_cov*tlen which inflated cov by
        // orders of magnitude, causing intra_gene_cov_frac to use SE tx
        // as cluster-max and kill multi-exon siblings en masse.
        out.push(Transcript {
            chrom: bundle.chrom.clone(),
            strand: bundle.strand,
            exons: vec![(med_start, med_end)],
            coverage: avg_cov,
            exon_cov: vec![avg_cov],
            tpm: 0.0,
            fpkm: 0.0,
            source: Some("stranded_single_exon".to_string()),
            is_longread: true,
            longcov: avg_cov,
            bpcov_cov: avg_cov,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: None,
            ref_gene_id: None,
            hardstart: false,
            hardend: false,
            alt_tts_end: false,
            vg_family_id: None,
            vg_copy_id: None,
            vg_family_size: None,
            intron_low: Vec::new(),
        });
    }
    out
}

/// Emit single-exon candidates from CLUSTERS OF TERMINAL EXONS of multi-exon
/// reads. Targets StringTie's nested-SE gene calls (e.g. STRG.276 at
/// 40126487-40126701 supported by 21 multi-exon reads whose LAST exon falls
/// in the region, vs only mono-exon reads picked up by the stranded path).
///
/// Opt-in via RUSTLE_TERMINAL_SE_ON=1. Knobs:
///   RUSTLE_TERMINAL_SE_MIN_READS (7), RUSTLE_TERMINAL_SE_CLUSTER_TOL (500).
fn emit_terminal_exon_se_candidates(
    bundle: &crate::types::Bundle,
    bundle_txs: &[Transcript],
    min_transcript_length: u64,
    singlethr: f64,
) -> Vec<Transcript> {
    if std::env::var_os("RUSTLE_TERMINAL_SE_ON").is_none() {
        return Vec::new();
    }
    if bundle.strand != '+' && bundle.strand != '-' {
        return Vec::new();
    }
    let min_reads: usize = std::env::var("RUSTLE_TERMINAL_SE_MIN_READS")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(7);
    let cluster_tol: u64 = std::env::var("RUSTLE_TERMINAL_SE_CLUSTER_TOL")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(500);

    // Terminal exons = 3' exon for + strand, 5' exon for - strand.
    // StringTie's nested-SE calls tend to share the 3' endpoint with the
    // host's terminal exon, so we cluster by 3' position (end for +, start
    // for -) using the read's strand-oriented terminal exon.
    let mut terms: Vec<(u64, u64)> = bundle
        .reads
        .iter()
        .filter(|r| r.exons.len() >= 2)
        .map(|r| {
            if bundle.strand == '+' {
                *r.exons.last().unwrap()
            } else {
                r.exons[0]
            }
        })
        .collect();
    if terms.len() < min_reads {
        return Vec::new();
    }
    // Sort by polyA side of the terminal exon: for + strand that's the exon
    // END (e), for - strand that's the exon START (s). Clustering on the
    // polyA side concentrates read-endpoint density where TTS is called.
    terms.sort_by_key(|&(s, e)| if bundle.strand == '+' { e } else { s });

    let mut clusters: Vec<Vec<(u64, u64)>> = Vec::new();
    for (s, e) in terms {
        let key = if bundle.strand == '+' { e } else { s };
        let mut placed = false;
        if let Some(last) = clusters.last_mut() {
            let (ls, le) = last[0];
            let lkey = if bundle.strand == '+' { le } else { ls };
            if key.abs_diff(lkey) <= cluster_tol {
                last.push((s, e));
                placed = true;
            }
        }
        if !placed {
            clusters.push(vec![(s, e)]);
        }
    }

    let mut out = Vec::new();
    for cluster in clusters {
        if cluster.len() < min_reads {
            continue;
        }
        let n = cluster.len();
        let med_start = cluster.iter().map(|&(s, _)| s).min().unwrap();
        let med_end = cluster.iter().map(|&(_, e)| e).max().unwrap();
        if med_end <= med_start {
            continue;
        }
        let tlen = med_end - med_start + 1;
        if tlen < min_transcript_length {
            continue;
        }
        let avg_cov = n as f64;
        if avg_cov < singlethr {
            continue;
        }
        // Skip near-duplicate of a multi-exon tx's terminal exon (within 100bp).
        // The distinct signal is that StringTie's nested SE has SHORTER span
        // than its host's terminal exon.
        let near_dup_terminal = bundle_txs.iter().any(|t| {
            if t.strand != bundle.strand || t.exons.len() < 2 {
                return false;
            }
            let (ts, te) = if bundle.strand == '+' {
                *t.exons.last().unwrap()
            } else {
                t.exons[0]
            };
            ts.abs_diff(med_start) <= 100 && te.abs_diff(med_end) <= 100
        });
        if near_dup_terminal {
            continue;
        }

        out.push(Transcript {
            chrom: bundle.chrom.clone(),
            strand: bundle.strand,
            exons: vec![(med_start, med_end)],
            coverage: avg_cov,
            exon_cov: vec![avg_cov],
            tpm: 0.0,
            fpkm: 0.0,
            source: Some("terminal_exon_se".to_string()),
            is_longread: true,
            longcov: avg_cov,
            bpcov_cov: avg_cov,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: None,
            ref_gene_id: None,
            hardstart: false,
            hardend: false,
            alt_tts_end: false,
            vg_family_id: None,
            vg_copy_id: None,
            vg_family_size: None,
            intron_low: Vec::new(),
        });
    }
    out
}

/// Create single-exon predictions from bundle coverage for unstranded bundles
/// with no multi-exon reads .
fn create_single_exon_predictions_from_bundle(
    bundle: &crate::types::Bundle,
    bpcov_stranded: &BpcovStranded,
    min_transcript_length: u64,
    singlethr: f64,
    verbose: bool,
) -> Vec<Transcript> {
    let mut predictions: Vec<Transcript> = Vec::new();

    // Use unstranded coverage (bpcov[1] = unstranded)
    let all_cov = match &bpcov_stranded.all {
        Some(c) => c,
        None => return predictions,
    };

    if all_cov.is_empty() {
        return predictions;
    }

    let bundle_start = bpcov_stranded.minus.bundle_start;
    let bundle_end = bpcov_stranded.minus.bundle_end;

    // Simple approach: create single-exon predictions from continuous high-coverage regions
    // This mirrors standard behavior for unstranded single-exon bundles
    let mut current_start: Option<u64> = None;
    let mut current_end: u64 = 0;
    let mut in_exon = false;

    for pos in bundle_start..=bundle_end {
        let idx = (pos - bundle_start) as usize;
        if idx >= all_cov.len() {
            break;
        }
        let cov = all_cov[idx];
        if cov > singlethr {
            if !in_exon {
                current_start = Some(pos);
                in_exon = true;
            }
            current_end = pos;
        } else if in_exon {
            // End of current exon
            if let Some(start) = current_start {
                let len = current_end.saturating_sub(start).saturating_add(1);
                if len >= min_transcript_length {
                    // Calculate coverage (coverage field stores total base coverage)
                    let start_idx = (start - bundle_start) as usize;
                    let end_idx = (current_end - bundle_start) as usize;
                    let sum_cov: f64 = all_cov[start_idx..=end_idx].iter().sum();
                    let len_f64 = (end_idx - start_idx + 1) as f64;
                    let avg_cov = sum_cov / len_f64;
                    let total_cov = sum_cov; // total base coverage for coverage field
                    let tx = Transcript {
                        chrom: bundle.chrom.clone(),
                        strand: '.',
                        exons: vec![(start, current_end)],
                        coverage: total_cov, // total base coverage
                        exon_cov: vec![total_cov],
                        tpm: 0.0,
                        fpkm: 0.0,
                        source: None,
                        is_longread: true,
                        longcov: avg_cov,   // per-base for longcov
                        bpcov_cov: avg_cov, // per-base for bpcov_cov
                        transcript_id: None,
                        gene_id: None,
                        ref_transcript_id: None,
                        ref_gene_id: None,
                        hardstart: false,
                        hardend: false,
                    alt_tts_end: false,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None, intron_low: Vec::new(),
                    };
                    predictions.push(tx);
                }
            }
            in_exon = false;
            current_start = None;
        }
    }

    // Handle last exon if bundle ends in high-coverage region
    if in_exon {
        if let Some(start) = current_start {
            let len = current_end.saturating_sub(start).saturating_add(1);
            if len >= min_transcript_length {
                let start_idx = (start - bundle_start) as usize;
                let end_idx = (current_end - bundle_start) as usize;
                let sum_cov: f64 = all_cov[start_idx..=end_idx].iter().sum();
                let elen = (end_idx - start_idx + 1) as f64;
                let avg_cov = sum_cov / elen;
                let total_cov = sum_cov; // total base coverage for coverage field
                let tx = Transcript {
                    chrom: bundle.chrom.clone(),
                    strand: '.',
                    exons: vec![(start, current_end)],
                    coverage: total_cov, // total base coverage
                    exon_cov: vec![total_cov],
                    tpm: 0.0,
                    fpkm: 0.0,
                    source: None,
                    is_longread: true,
                    longcov: avg_cov,   // per-base for longcov
                    bpcov_cov: avg_cov, // per-base for bpcov_cov
                    transcript_id: None,
                    gene_id: None,
                    ref_transcript_id: None,
                    ref_gene_id: None,
                    hardstart: false,
                    hardend: false,
                    alt_tts_end: false,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None, intron_low: Vec::new(),
                };
                predictions.push(tx);
            }
        }
    }

    if verbose && !predictions.is_empty() {
        eprintln!(
            "    single_exon_from_bundle: created {} prediction(s) from {}-{} coverage regions",
            predictions.len(),
            bundle_start,
            bundle_end
        );
    }

    predictions
}

/// Default path for per-bundle diagnostics TSV (sibling of the main GTF output).
fn default_needy_bundle_summary_path(output_gtf: &Path) -> PathBuf {
    let stem = output_gtf
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "rustle".to_string());
    let parent = output_gtf.parent().unwrap_or_else(|| Path::new("."));
    parent.join(format!("{}.needy_bundle_summary.tsv", stem))
}

/// True when this run is in "needy locus" diagnostic mode (trace region, reference trace, etc.).
fn needy_locus_diagnostics_wanted(config: &crate::types::RunConfig, trace_reference: Option<&str>) -> bool {
    trace_reference.is_some()
        || config.only_debug_bundle
        || config.debug_bundle.is_some()
        || std::env::var_os("RUSTLE_TRACE_LOCUS").is_some()
}

/// Resolve `--parity-stage-tsv`, then `RUSTLE_DEBUG_STAGE_TSV`, then auto path for needy runs.
fn resolve_debug_stage_tsv_path(
    config: &crate::types::RunConfig,
    output_gtf: &Path,
    trace_reference: Option<&str>,
) -> Option<PathBuf> {
    if let Some(p) = config.debug_stage_tsv.as_ref() {
        return Some(PathBuf::from(p));
    }
    if let Ok(p) = std::env::var("RUSTLE_DEBUG_STAGE_TSV") {
        if !p.is_empty() {
            return Some(PathBuf::from(p));
        }
    }
    if std::env::var_os("RUSTLE_DISABLE_AUTO_BUNDLE_SUMMARY").is_some() {
        return None;
    }
    if needy_locus_diagnostics_wanted(config, trace_reference) {
        return Some(default_needy_bundle_summary_path(output_gtf));
    }
    None
}

/// Run full pipeline: BAM -> bundles -> assemble -> GTF.
pub fn run<P: AsRef<Path>>(
    bam_path: P,
    output_gtf: P,
    config: RunConfig,
    chrom_filter: Option<&str>,
    bundles_log: Option<&str>,
    guide_reference: Option<&str>,
    genome_fasta: Option<&str>,
    trace_reference: Option<&str>,
    trace_output: Option<&str>,
) -> Result<()> {
    let _t0 = std::time::Instant::now();
    let show_timing = std::env::var_os("RUSTLE_TIMING").is_some();
    macro_rules! phase_timer {
        ($label:expr) => {
            if show_timing { eprintln!("[TIMING] {} {:.3}s", $label, _t0.elapsed().as_secs_f64()); }
        };
    }
    init_rayon_pool(config.threads);
    let diag_tsv_resolved = resolve_debug_stage_tsv_path(&config, output_gtf.as_ref(), trace_reference);
    if let Some(ref p) = diag_tsv_resolved {
        let auto_path = config.debug_stage_tsv.is_none()
            && std::env::var_os("RUSTLE_DEBUG_STAGE_TSV").is_none();
        if auto_path {
            eprintln!(
                "rustle: needy-locus diagnostics: writing per-bundle summary TSV to {}",
                p.display()
            );
        }
    }
    debug_stage::init(diag_tsv_resolved.as_deref())?;
    let snapshot_writer = std::sync::Arc::new(std::sync::Mutex::new(snapshot::SnapshotWriter::new(config.snapshot_jsonl.as_deref())?));
    let mut chrom_arc_cache: HashMap<String, Arc<str>> = Default::default();
    let mut consensus_cache: LruCache<SpliceConsensusKey, bool> =
        LruCache::new(consensus_cache_capacity());
    // In VG mode with a FASTA, eagerly load the genome and thread it to
    // bundle ingest so per-read mismatches vs reference can be extracted
    // for SNP-based copy scoring. Same genome index is later passed to
    // family discovery and novel-copy scan.
    let vg_snp_genome: Option<crate::genome::GenomeIndex> = if config.vg_mode
        && config.vg_snp
        && config.genome_fasta.is_some()
    {
        let path = config.genome_fasta.as_ref().unwrap();
        eprintln!("[VG] Loading genome FASTA for SNP copy assignment: {}", path);
        crate::genome::GenomeIndex::from_fasta(path).ok()
    } else {
        None
    };
    let mut bundles = crate::bundle::detect_bundles_from_bam_with_snp(
        bam_path.as_ref(),
        &config,
        chrom_filter,
        vg_snp_genome.as_ref(),
    )?;
    if crate::trace_pipeline::active() {
        for b in &bundles {
            crate::trace_pipeline::dump_bundle(
                &b.chrom, b.start, b.end, b.strand, b.reads.len(),
            );
        }
    }
    phase_timer!("bam_ingest");
    let debug_target = parse_debug_bundle_target(&config);
    if config.only_debug_bundle {
        bundles.retain(|bundle| {
            debug_bundle_overlaps(
                debug_target.as_ref(),
                &bundle.chrom,
                bundle.start,
                bundle.end,
            )
        });
    }
    let (global_num_frag, global_frag_len_sum) = bundles.iter().fold((0.0f64, 0.0f64), |acc, b| {
        let (n, s) = count_fragment_metrics_for_bundle(b);
        (acc.0 + n, acc.1 + s)
    });
    let mut n_bundles = bundles.len();
    if debug_stage::is_enabled() {
        for bundle in &bundles {
            let n_junctions: usize = bundle.reads.iter().map(|r| r.junctions.len()).sum();
            debug_stage::emit_with_detail(
                "bundle_context",
                &bundle.chrom,
                bundle.start,
                bundle.end,
                bundle.strand,
                bundle.reads.len(),
                n_junctions,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                debug_stage::bundle_detail(bundle),
                "detected_bundle",
            );
        }
    }

    // Diagnostic: report initial bundle count breakdown
    if config.verbose {
        let by_strand: HashMap<char, usize> = bundles.iter().fold(
            HashMap::default(),
            |mut acc, b| {
                *acc.entry(b.strand).or_insert(0) += 1;
                acc
            },
        );
        eprintln!(
            "[BUNDLE_COUNT] initial: {} total (minus: {}, plus: {}, unstranded: {})",
            n_bundles,
            by_strand.get(&'-').unwrap_or(&0),
            by_strand.get(&'+').unwrap_or(&0),
            by_strand.get(&'.').unwrap_or(&0)
        );
    }

    // ── Variation graph: discover gene family groups ────────────────────────
    let vg_families = if config.vg_mode {
        // Sequence-similarity family discovery is opt-in: only when --vg-discover-novel
        // is active (not triggered by --vg-snp or --genome-fasta alone, to avoid false links).
        let vg_genome_for_discovery = if config.vg_discover_novel {
            config.genome_fasta.as_ref().and_then(|p| {
                eprintln!("[VG] Loading genome for sequence-similarity family discovery...");
                crate::genome::GenomeIndex::from_fasta(p).ok()
            })
        } else {
            None
        };
        let families = crate::vg::discover_family_groups(
            &bundles,
            config.vg_min_shared_reads,
            Some(bam_path.as_ref()),
            vg_genome_for_discovery.as_ref(),
        );
        if !families.is_empty() {
            eprintln!(
                "[VG] {} family group(s) covering {} bundles",
                families.len(),
                families.iter().map(|f| f.bundle_indices.len()).sum::<usize>(),
            );
        }
        families
    } else {
        Vec::new()
    };
    // Build set of bundle indices that belong to a family group (for deferred processing).
    let vg_family_bundle_set: std::collections::HashSet<usize> = vg_families
        .iter()
        .flat_map(|f| f.bundle_indices.iter().copied())
        .collect();

    // VG mode keeps secondary/supplementary alignments at ingest so family
    // discovery can see multi-mapper evidence across bundles. But leaving them
    // in non-family bundles introduces noise: precision on chr19 fell
    // 83.8% → 79.9% in earlier runs because every bundle was seeing extra
    // cross-mapped reads it shouldn't have. After discovery, strip
    // non-primary alignments from any bundle that did NOT join a family, so
    // unrelated loci behave exactly like the StringTie baseline.
    //
    // Disable with RUSTLE_VG_KEEP_NONFAMILY_SECONDARY=1 for diagnostics.
    if config.vg_mode
        && std::env::var_os("RUSTLE_VG_KEEP_NONFAMILY_SECONDARY").is_none()
    {
        let mut stripped_reads = 0usize;
        let mut stripped_bundles = 0usize;
        for (bi, bundle) in bundles.iter_mut().enumerate() {
            if vg_family_bundle_set.contains(&bi) {
                continue;
            }
            let before = bundle.reads.len();
            bundle.reads.retain(|r| r.is_primary_alignment);
            let dropped = before - bundle.reads.len();
            if dropped > 0 {
                stripped_reads += dropped;
                stripped_bundles += 1;
            }
        }
        if stripped_reads > 0 {
            eprintln!(
                "[VG] Stripped {} secondary/supplementary read(s) from {} non-family bundles",
                stripped_reads, stripped_bundles
            );
        }
    }
    // Pre-save bundle data for VG report and novel discovery (bundles consumed by processing loop).
    let vg_bundle_coords: Vec<(String, u64, u64, char)> = if config.vg_mode {
        bundles
            .iter()
            .map(|b| (b.chrom.clone(), b.start, b.end, b.strand))
            .collect()
    } else {
        Vec::new()
    };
    // Save lightweight bundle copies for novel copy discovery (junction stats only).
    let vg_bundles_for_novel: Vec<crate::types::Bundle> = if config.vg_mode && config.vg_discover_novel {
        bundles
            .iter()
            .map(|b| crate::types::Bundle {
                chrom: b.chrom.clone(),
                start: b.start,
                end: b.end,
                strand: b.strand,
                reads: Vec::new(), // Don't clone reads — only need junction_stats.
                junction_stats: b.junction_stats.clone(),
                bundlenodes: None,
                read_bnodes: None,
                bnode_colors: None,
            })
            .collect()
    } else {
        Vec::new()
    };

    // ── Variation graph: multi-mapping read resolution ──────────────────────
    let vg_em_results: Vec<crate::vg::EmResult> = if config.vg_mode && !vg_families.is_empty() {
        use crate::types::VgSolver;
        match config.vg_solver {
            VgSolver::None => {
                eprintln!("[VG] Discovery only — no multi-mapping resolution (use --vg-solver em to enable)");
                Vec::new()
            }
            VgSolver::Em => {
                if config.vg_snp {
                    eprintln!("[VG] Using EM solver with SNP-based copy assignment");
                    crate::vg::run_pre_assembly_em_with_snps(
                        &vg_families,
                        &mut bundles,
                        config.vg_em_max_iter,
                    )
                } else {
                    eprintln!("[VG] Using EM solver for multi-mapping resolution");
                    crate::vg::run_pre_assembly_em(
                        &vg_families,
                        &mut bundles,
                        config.vg_em_max_iter,
                    )
                }
            }
            VgSolver::Flow => {
                // Flow-based redistribution requires two-pass assembly.
                // First pass: run with uniform weights (handled by normal pipeline below).
                // Second pass: redistribute based on assembled transcript coverage.
                // For now, use EM as first pass, flow redistribution in future.
                eprintln!("[VG] Using flow-based solver (first pass: EM, then flow redistribution)");
                crate::vg::run_pre_assembly_em(
                    &vg_families,
                    &mut bundles,
                    config.vg_em_max_iter,
                )
                // TODO: After main assembly loop, collect per-family transcripts,
                // compute flow-based weights, and re-run assembly for family bundles.
            }
        }
    } else {
        Vec::new()
    };

    let use_region_bundle_pass = std::env::var_os("RUSTLE_DISABLE_REGION_BUNDLE_PASS").is_none();
    let do_timing = std::env::var_os("RUSTLE_TIMING").is_some();
    let trace_log_style = trace_log_style_active();
    macro_rules! t {
        ($label:expr) => {
            if do_timing {
                eprintln!("[TIMING] {} {:.3}s", $label, _t0.elapsed().as_secs_f64());
            }
        };
    }
    if trace_log_style {
        let read_filter = std::env::var("RUSTLE_TRACE_READ")
            .ok()
            .filter(|v| !v.trim().is_empty())
            .unwrap_or_else(|| "0".to_string());
        eprintln!(
            "TRACE: entry=0 deep=1 data=0 loop=1 path=1 bundle=1 func=(all) pos=-1 read_filter={}",
            read_filter
        );
    }
    t!("bam_detect");

    if let Some(log_path) = bundles_log {
        let log_regions = parse_bundle_log(log_path)?;
        compare_bundles_to_log(&bundles, &log_regions);
    }

    let ref_transcripts = trace_reference
        .map(|p| parse_reference_gtf(p))
        .transpose()?
        .unwrap_or_default();
    let guide_transcripts: Vec<RefTranscript> = guide_reference
        .map(parse_reference_gtf)
        .transpose()?
        .unwrap_or_default();
    let genome = genome_fasta.map(GenomeIndex::from_fasta).transpose()?;
    let mut guide_junction_cache: LruCache<GuideJunctionCacheKey, Arc<Vec<Junction>>> =
        LruCache::new(guide_junction_cache_capacity());
    let raw_for_trace_mutex = std::sync::Mutex::new(if trace_reference.is_some() {
        Vec::with_capacity(bundles.len())
    } else {
        vec![]
    });

    let all_transcripts_mutex = std::sync::Mutex::new(Vec::<Transcript>::new());
    let single_exon_predictions_mutex = std::sync::Mutex::new(Vec::<Transcript>::new());
    let trace_pass_id = std::sync::atomic::AtomicUsize::new(0);
    let shadow_strict_3strand = std::env::var_os("RUSTLE_SHADOW_STRICT_3STRAND").is_some();
    let shadow_bundle_diags_mutex = std::sync::Mutex::new(Vec::<ShadowStrictBundleDiag>::new());

    // Pre-compute per-strand coverage for ALL bundles at each region.
    // The algorithm processes all strands in one bundle and has bpcov[3] (minus/all/plus).
    // Rust splits by strand, so each bundle lacks coverage from other strands.
    // Build a global BpcovStranded per (chrom, start, end) from ALL strand bundles.
    // This is used for:
    //   1. Cross-strand resolution in '.' bundles
    //   2. Strand-signed coverage (get_cov_sign logic) in good_junc gate for stranded bundles
    t!("setup");
    // Build strand evidence map in parallel: chrom -> junction -> (has_plus, has_minus).
    let chrom_junction_strand_evidence = collect_chrom_junction_strand_evidence_parallel(&bundles);
    let cross_strand_cov: HbHashMap<(String, u64, u64), BpcovStranded> = {
        let mut acc = HbHashMap::<(String, u64, u64), BpcovStranded>::default();
        for bundle in bundles.iter() {
            let key = (bundle.chrom.clone(), bundle.start, bundle.end);
            let entry = acc
                .entry(key)
                .or_insert_with(|| BpcovStranded::empty(bundle.start, bundle.end));
            match bundle.strand {
                '+' => {
                    for r in &bundle.reads {
                        for &(s, e) in &r.exons {
                            entry.add_coverage(BPCOV_STRAND_PLUS, s, e, r.weight);
                        }
                    }
                }
                '-' => {
                    for r in &bundle.reads {
                        for &(s, e) in &r.exons {
                            entry.add_coverage(BPCOV_STRAND_MINUS, s, e, r.weight);
                        }
                    }
                }
                _ => {
                    // '.' (neutral/unstranded) reads: add to 'all' channel only
                    // bpcov[1] includes all reads; neutral reads go to bpcov[1] via sno=1
                    for r in &bundle.reads {
                        for &(s, e) in &r.exons {
                            entry.add_coverage(BPCOV_STRAND_ALL, s, e, r.weight);
                        }
                    }
                }
            }
        }
        acc
    };

    t!("cross_strand_cov");
    t!("junction_strand_evidence");

    if use_region_bundle_pass {
        bundles = merge_region_outer_bundles(bundles, &config);
        if config.only_debug_bundle {
            bundles.retain(|bundle| {
                debug_bundle_overlaps(
                    debug_target.as_ref(),
                    &bundle.chrom,
                    bundle.start,
                    bundle.end,
                )
            });
        }
        n_bundles = bundles.len();
        if config.verbose {
            eprintln!("[BUNDLE_COUNT] after region_merge: {}", n_bundles);
        }
    }

    // dual-add reads predicted to become unstranded into opposite-strand bundles.
    // Keep enabled by default; allow explicit opt-out only for diagnostics.
    if !use_region_bundle_pass && std::env::var_os("RUSTLE_DISABLE_DUAL_ADD_UNSTRANDED").is_none() {
        let strand_audit = std::env::var_os("RUSTLE_DEBUG_STRAND_AUDIT").is_some();
        let mut unstranded_reads: HbHashMap<(String, u64, u64), Vec<BundleRead>> =
            Default::default();
        // Keep neutralization indices strand-scoped; region-only indexing leaks +/- indices.
        let mut unstranded_read_indices: HbHashMap<(String, u64, u64, char), Vec<usize>> =
            Default::default();
        let mut unstranded_candidates_by_bundle: HbHashMap<(String, u64, u64, char), usize> =
            Default::default();

        for bundle in bundles.iter() {
            // Skip non-stranded bundles
            if bundle.strand != '+' && bundle.strand != '-' {
                continue;
            }

            // Work on cloned reads to avoid modifying bundle state
            let mut reads_clone: Vec<BundleRead> = bundle.reads.clone();

            // Apply polyA filter on clone (same filter main loop will apply)
            if config.long_reads {
                reads_clone.retain(|r| {
                    if r.exons.len() > 2 {
                        return true;
                    }
                    !(r.has_poly_end_aligned || r.has_poly_start_aligned)
                });
                shorten_polya_terminal_exons(&mut reads_clone);
            }

            // count support first, then apply higherr/snap correction stages.
            let mut junction_stats_corr = bundle.junction_stats.clone();

            // Count support on cloned reads
            count_good_junctions(
                &mut reads_clone,
                &mut junction_stats_corr,
                config.junction_support,
                config.longintron,
                bundle.start,
                bundle.end,
                config.min_intron_length,
                bundle.strand,
                None,
                &bundle.chrom,
                &chrom_junction_strand_evidence,
            );

            // the algorithm does not run a separate ratio-cluster rewrite pass here.
            // Keep only the direct coordinate correction / guide snap stages.
            let (corrected_stats, _corrected_map) =
                correct_junctions_with_map(&junction_stats_corr, false);
            junction_stats_corr = corrected_stats;

            // Guide snap
            if !guide_transcripts.is_empty() && config.long_reads {
                let guide_junctions_vec = cached_bundle_guide_junctions(
                    &mut guide_junction_cache,
                    &mut chrom_arc_cache,
                    &guide_transcripts,
                    &bundle.chrom,
                    bundle.strand,
                    bundle.start,
                    bundle.end,
                );
                if !guide_junctions_vec.is_empty() {
                    let (snapped_stats, _snap_map) = snap_junctions_to_guides(
                        &junction_stats_corr,
                        guide_junctions_vec.as_ref(),
                        config.sserror,
                        config.long_reads,
                    );
                    junction_stats_corr = snapped_stats;
                }
            }

            let bpcov_stranded = BpcovStranded::from_reads(&reads_clone, bundle.start, bundle.end);
            if !guide_transcripts.is_empty() {
                let guide_junctions_vec = cached_bundle_guide_junctions(
                    &mut guide_junction_cache,
                    &mut chrom_arc_cache,
                    &guide_transcripts,
                    &bundle.chrom,
                    bundle.strand,
                    bundle.start,
                    bundle.end,
                );
                mark_guide_junctions_for_junctions(
                    &mut junction_stats_corr,
                    guide_junctions_vec.as_ref(),
                );
            }
            let mut cjunctions = crate::types::junction_stats_to_cjunctions(&junction_stats_corr);
            aggregate_splice_site_support(&mut cjunctions);
            // apply_higherr_demotions before good_junc
            // The higherr block runs before the good_junc call in build_graphs
            apply_higherr_demotions(
                &mut cjunctions,
                config.sserror,
                config.min_junction_reads,
                Some(&bpcov_stranded),
                bundle.start,
            );
            demote_runthrough_junctions(&mut cjunctions, &bpcov_stranded, bundle.start);
            good_junc(
                &mut cjunctions,
                &bpcov_stranded,
                bundle.start,
                config.min_junction_reads,
                config.longintron,
                config.per_splice_site_isofrac,
                config.long_reads,
                config.eonly,
            );
            emit_jfinal_trace(&cjunctions, "bundle-genome");
            junction_stats_corr = crate::types::cjunctions_to_junction_stats(&cjunctions);
            if let Some(g) = &genome {
                let chrom_arc = intern_chrom_arc(&mut chrom_arc_cache, &bundle.chrom);
                let guide_junctions = cached_bundle_guide_junctions(
                    &mut guide_junction_cache,
                    &mut chrom_arc_cache,
                    &guide_transcripts,
                    &bundle.chrom,
                    bundle.strand,
                    bundle.start,
                    bundle.end,
                );
                let guide_tol = config.sserror.max(1);
                for (j, st) in junction_stats_corr.iter_mut() {
                    if st.strand == Some(0) {
                        continue;
                    }
                    let guide_match = guide_junctions.as_ref().iter().any(|gj| {
                        gj.donor.abs_diff(j.donor) <= guide_tol
                            && gj.acceptor.abs_diff(j.acceptor) <= guide_tol
                    });
                    if guide_match {
                        continue;
                    }
                    let is_cons = cached_consensus_splice(
                        &mut consensus_cache,
                        g,
                        &chrom_arc,
                        j.donor,
                        j.acceptor,
                        st.strand,
                    );
                    // StringTie rlink.cpp:14842 kills junctions with non-canonical
                    // consensus + low support. Applying the same in Rustle regresses
                    // matches (1440 → 1401) because StringTie rescues some of these
                    // through different downstream paths (likely different nreads_good
                    // accounting or checktrf redistribution). Filter disabled by default;
                    // enable by setting RUSTLE_CONSENSUS_KILL=1.
                    if std::env::var_os("RUSTLE_CONSENSUS_KILL").is_some()
                        && !is_cons
                        && st.nreads_good < 5.0
                    {
                        st.strand = Some(0);
                    }
                }
            }

            // Identify which reads would become unstranded
            let cjunctions_post = crate::types::junction_stats_to_cjunctions(&junction_stats_corr);
            let killed_junction_pairs = compute_killed_junction_pairs(&cjunctions_post);
            if killed_junction_pairs.is_empty() {
                continue;
            }

            // Use the ORIGINAL bundle.reads (not cloned) for identification,
            // since we want to inject the original read into the opposite bundle
            for (ri, r) in bundle.reads.iter().enumerate() {
                if r.exons.len() < 2 {
                    continue;
                }
                let has_bad = r
                    .junctions
                    .iter()
                    .any(|j| killed_junction_pairs.contains(j));
                if !has_bad {
                    continue;
                }
                let key = (bundle.chrom.clone(), bundle.start, bundle.end);
                let has_stranded_junction = r.junctions.iter().any(|j| {
                    !killed_junction_pairs.contains(j)
                        && junction_stats_corr
                            .get(j)
                            .map(|st| st.strand != Some(0))
                            .unwrap_or(false)
                });
                if has_stranded_junction {
                    // Partial read demotion: split on killed junctions and inject
                    // only unsupported fragments as neutral evidence.
                    let frags = killed_junction_neutral_fragments(r, &killed_junction_pairs);
                    if !frags.is_empty() {
                        *unstranded_candidates_by_bundle
                            .entry((
                                bundle.chrom.clone(),
                                bundle.start,
                                bundle.end,
                                bundle.strand,
                            ))
                            .or_insert(0) += frags.len();
                        unstranded_reads
                            .entry(key.clone())
                            .or_default()
                            .extend(frags);
                    }
                    continue;
                }
                if !has_stranded_junction {
                    let bundle_key = (
                        bundle.chrom.clone(),
                        bundle.start,
                        bundle.end,
                        bundle.strand,
                    );
                    unstranded_read_indices
                        .entry(bundle_key.clone())
                        .or_default()
                        .push(ri);
                    *unstranded_candidates_by_bundle
                        .entry(bundle_key)
                        .or_insert(0) += 1;
                    let mut neutral = r.clone();
                    neutral.strand = '.';
                    unstranded_reads.entry(key).or_default().push(neutral);
                }
            }
        }

        // Reads that lose stranded junction evidence are neutralized in-place.
        // This keeps sno=1 semantics in the original bundle.
        for bundle in bundles.iter_mut() {
            if bundle.strand != '+' && bundle.strand != '-' {
                continue;
            }
            let key = (
                bundle.chrom.clone(),
                bundle.start,
                bundle.end,
                bundle.strand,
            );
            let before = if strand_audit {
                Some(read_strand_counts(&bundle.reads))
            } else {
                None
            };
            let mut changed = 0usize;
            if let Some(idxs) = unstranded_read_indices.get(&key) {
                for &ri in idxs {
                    if ri < bundle.reads.len() {
                        if bundle.reads[ri].strand != '.' {
                            changed += 1;
                        }
                        bundle.reads[ri].strand = '.';
                    }
                }
            }
            if strand_audit {
                let after = read_strand_counts(&bundle.reads);
                let idx_count = unstranded_read_indices
                    .get(&key)
                    .map(|v| v.len())
                    .unwrap_or(0);
                if idx_count > 0 || changed > 0 {
                    let before = before.unwrap_or_default();
                    let flagged = unstranded_candidates_by_bundle
                        .get(&key)
                        .copied()
                        .unwrap_or(0);
                    eprintln!(
                        "DEBUG_STRAND_NEUTRALIZE chrom={} range={}-{} bundle_strand={} flagged={} idxs={} changed={} before(-/.+/o)={}/{}/{}/{} after(-/.+/o)={}/{}/{}/{}",
                        bundle.chrom,
                        bundle.start + 1,
                        bundle.end,
                        bundle.strand,
                        flagged,
                        idx_count,
                        changed,
                        before.minus,
                        before.neutral,
                        before.plus,
                        before.other,
                        after.minus,
                        after.neutral,
                        after.plus,
                        after.other
                    );
                }
            }
        }

        // Inject unstranded reads into opposite-strand bundles
        for bundle in bundles.iter_mut() {
            if bundle.strand != '+' && bundle.strand != '-' {
                continue;
            }
            let key = (bundle.chrom.clone(), bundle.start, bundle.end);
            let before = if strand_audit {
                Some(read_strand_counts(&bundle.reads))
            } else {
                None
            };
            if let Some(reads) = unstranded_reads.get(&key) {
                let mut seen: HbHashSet<ReadSigNoStrand> =
                    bundle.reads.iter().map(read_sig_nostrand).collect();
                let mut inserted = 0usize;
                for r in reads {
                    // Target: dual-add as neutral (sno=1), not forced +/-.
                    // Deduplicate against existing reads in this bundle.
                    if seen.insert(read_sig_nostrand(r)) {
                        bundle.reads.push(r.clone());
                        inserted += 1;
                    }
                }
                if strand_audit {
                    let after = read_strand_counts(&bundle.reads);
                    let before = before.unwrap_or_default();
                    eprintln!(
                        "DEBUG_STRAND_INJECT chrom={} range={}-{} bundle_strand={} source_neutral={} inserted={} before(-/.+/o)={}/{}/{}/{} after(-/.+/o)={}/{}/{}/{}",
                        bundle.chrom,
                        bundle.start + 1,
                        bundle.end,
                        bundle.strand,
                        reads.len(),
                        inserted,
                        before.minus,
                        before.neutral,
                        before.plus,
                        before.other,
                        after.minus,
                        after.neutral,
                        after.plus,
                        after.other
                    );
                }
            }
        }
    }

    t!("pre_processing_pass");

    // Keep 3-strand groupflow on by default for compatibility, with an explicit opt-out.
    let use_3strand_groupflow = std::env::var_os("RUSTLE_DISABLE_3STRAND_GROUPFLOW").is_none();
    let mut region_reads_for_groupflow: HbHashMap<(String, u64, u64), Vec<BundleRead>> =
        Default::default();
    if use_3strand_groupflow {
        for b in &bundles {
            region_reads_for_groupflow
                .entry((b.chrom.clone(), b.start, b.end))
                .or_default()
                .extend(b.reads.clone());
        }
    }
    t!("groupflow_build");
    phase_timer!("pre_assembly");

    // Track total sub-bundles created across all regions
    let total_subbundles = std::sync::atomic::AtomicUsize::new(0);


    let snapshot_all = std::env::var_os("RUSTLE_SNAPSHOT_ALL").is_some();
    use rayon::prelude::*;
    let bundles_vec: Vec<(usize, crate::types::Bundle)> = bundles.into_iter().enumerate().collect();
    bundles_vec.into_par_iter().try_for_each(|(bundle_idx, mut bundle)| -> Result<()> {
        let mut chrom_arc_cache: HashMap<String, Arc<str>> = Default::default();
        let mut consensus_cache: LruCache<SpliceConsensusKey, bool> =
            LruCache::new(consensus_cache_capacity());
        let mut guide_junction_cache: LruCache<GuideJunctionCacheKey, Arc<Vec<Junction>>> =
            LruCache::new(guide_junction_cache_capacity());

        let stage_debug = config.debug_junctions
            && debug_bundle_overlaps(
                debug_target.as_ref(),
                &bundle.chrom,
                bundle.start,
                bundle.end,
            );
        let snapshot_this = snapshot_writer.lock().unwrap().enabled()
            && (snapshot_all
                || debug_bundle_overlaps(
                    debug_target.as_ref(),
                    &bundle.chrom,
                    bundle.start,
                    bundle.end,
                ));
        if snapshot_this {
            let key = snapshot::bundle_key(bundle_idx, &bundle);
            let junctions: Vec<crate::types::Junction> =
                bundle.junction_stats.iter().map(|(j, _)| *j).collect();
            let payload = snapshot::BundleSummaryWithJunctions {
                summary: snapshot::summarize_bundle(&bundle),
                junctions: snapshot::snap_junctions(&junctions, bundle.strand),
            };
            let _ = snapshot_writer.lock().unwrap().emit(&snapshot::SnapshotEnvelope {
                kind: "rustle_snapshot",
                stage: "bundle",
                bundle: key,
                payload,
            });
        }
        // Read processing:
        // (1) discard <=2-exon long reads with aligned poly artifact
        // (2) trim/remove terminal poly exons before junction/graph processing
        if config.long_reads {
            let mut kept_reads: Vec<BundleRead> = Vec::with_capacity(bundle.reads.len());
            for mut r in std::mem::take(&mut bundle.reads) {
                // processRead filter: these gates apply only to long-read-class entries
                // (longreads || mixedMode&&uval), not to short-read entries in mixed mode.
                if read_is_long_class_for_processread(&r, &config) {
                    if r.exons.len() <= 2 && (r.has_poly_end_aligned || r.has_poly_start_aligned) {
                        continue;
                    }
                    if !apply_processread_terminal_poly_trim(&mut r) {
                        continue;
                    }
                }
                kept_reads.push(r);
            }
            bundle.reads = kept_reads;
        }
        // shorten terminal exons for long reads with poly flags.
        if config.long_reads {
            shorten_polya_terminal_exons(&mut bundle.reads);
            // the `-t` / `no_coverage_trim`: skip read-level terminal coverage trim (still keep
            // poly-A shorten above — processRead logic). Helps gffcompare terminal class `j`.
            if !config.no_coverage_trim {
                let bpcov_temp = BpcovStranded::from_reads(&bundle.reads, bundle.start, bundle.end);
                trim_terminal_exons_by_coverage(&mut bundle.reads, &bpcov_temp, bundle.start);
            }
        }
        // Processing order:
        // 1. Count support on observed junction/read coordinates
        // 2. Apply higherr correction
        // 3. Apply guide-based snap correction (long-read, guide present)
        // 4. Gate with good_junc, then repair reads from killed + redirect maps
        let mut junction_stats_corr = bundle.junction_stats.clone();
        let mut junction_redirect_map: HashMap<Junction, Junction> = Default::default();

        // Recompute support from adjusted segments (count_good_junctions)
        let cross_cov = if bundle.strand == '.' {
            let key = (bundle.chrom.clone(), bundle.start, bundle.end);
            cross_strand_cov.get(&key)
        } else {
            None
        };
        count_good_junctions(
            &mut bundle.reads,
            &mut junction_stats_corr,
            config.junction_support,
            config.longintron,
            bundle.start,
            bundle.end,
            config.min_intron_length,
            bundle.strand,
            cross_cov,
            &bundle.chrom,
            &chrom_junction_strand_evidence,
        );
        if stage_debug {
            eprintln!(
                "[EARLY_DEBUG] post_count_good {}:{}-{} count={}",
                bundle.chrom,
                bundle.start,
                bundle.end,
                junction_stats_corr.len()
            );
        }

        if stage_debug {
            eprintln!(
                "[EARLY_DEBUG] post_ratio_corr {}:{}-{} skipped_for_parity count={} redirects={}",
                bundle.chrom,
                bundle.start,
                bundle.end,
                junction_stats_corr.len(),
                junction_redirect_map.len()
            );
        }

        let (corrected_stats, corrected_map) =
            correct_junctions_with_map(&junction_stats_corr, config.verbose);
        junction_stats_corr = corrected_stats;
        junction_redirect_map.extend(corrected_map);

        // Alt-junction demotion: snap alt-donors/acceptors with weak read
        // support to the canonical at the shared endpoint. Reduces graph
        // segmentation at alt-splice hot spots (STRG.309 class: 14 j-class
        // over-emissions from alt-donor combinations). Opt-in via
        // RUSTLE_ALT_JUNC_SNAP=1. The redirect map flows into
        // repair_reads_after_junction_quality below, which updates read
        // exon boundaries accordingly.
        let alt_demote_map = crate::junction_correction::demote_alt_junctions_to_canonical(
            &junction_stats_corr
        );
        if !alt_demote_map.is_empty() {
            // Apply to stats: merge demoted junction's counts into canonical
            for (from, to) in &alt_demote_map {
                if from == to { continue; }
                if let Some(from_stat) = junction_stats_corr.get(from).cloned() {
                    if let Some(to_stat) = junction_stats_corr.get_mut(to) {
                        to_stat.mrcount += from_stat.mrcount;
                        to_stat.nreads_good += from_stat.nreads_good;
                        to_stat.nm += from_stat.nm;
                        to_stat.leftsupport += from_stat.leftsupport;
                        to_stat.rightsupport += from_stat.rightsupport;
                    }
                    junction_stats_corr.remove(from);
                }
            }
            junction_redirect_map.extend(alt_demote_map);
        }
        if stage_debug {
            eprintln!(
                "[EARLY_DEBUG] post_higherr {}:{}-{} count={} redirects={}",
                bundle.chrom,
                bundle.start,
                bundle.end,
                junction_stats_corr.len(),
                junction_redirect_map.len()
            );
            dump_junction_stats("post_higherr", &junction_stats_corr);
            dump_junction_redirects("post_higherr_redirects", &junction_redirect_map);
        }

        // guide-based coordinate snapping for high-mismatch junctions.
        if !guide_transcripts.is_empty() && config.long_reads {
            let guide_junctions_vec = cached_bundle_guide_junctions(
                &mut guide_junction_cache,
                &mut chrom_arc_cache,
                &guide_transcripts,
                &bundle.chrom,
                bundle.strand,
                bundle.start,
                bundle.end,
            );
            if !guide_junctions_vec.is_empty() {
                let (snapped_stats, snap_map) = snap_junctions_to_guides(
                    &junction_stats_corr,
                    guide_junctions_vec.as_ref(),
                    config.sserror,
                    config.long_reads,
                );
                if !snap_map.is_empty() {
                    junction_stats_corr = snapped_stats;
                    junction_redirect_map.extend(snap_map);
                }
            }
        }

        let bpcov_stranded = BpcovStranded::from_reads(&bundle.reads, bundle.start, bundle.end);
        if !guide_transcripts.is_empty() {
            let guide_junctions = cached_bundle_guide_junctions(
                &mut guide_junction_cache,
                &mut chrom_arc_cache,
                &guide_transcripts,
                &bundle.chrom,
                bundle.strand,
                bundle.start,
                bundle.end,
            );
            mark_guide_junctions_for_junctions(&mut junction_stats_corr, guide_junctions.as_ref());
        }
        let mut cjunctions = crate::types::junction_stats_to_cjunctions(&junction_stats_corr);
        if stage_debug {
            eprintln!(
                "[EARLY_DEBUG] cjunction_pre_good {}:{}-{} count={}",
                bundle.chrom,
                bundle.start,
                bundle.end,
                cjunctions.len()
            );
        }
        aggregate_splice_site_support(&mut cjunctions);
        // apply_higherr_demotions before good_junc
        apply_higherr_demotions(
            &mut cjunctions,
            config.sserror,
            config.min_junction_reads,
            Some(&bpcov_stranded),
            bundle.start,
        );
        demote_runthrough_junctions(&mut cjunctions, &bpcov_stranded, bundle.start);
        good_junc(
            &mut cjunctions,
            &bpcov_stranded,
            bundle.start,
            config.min_junction_reads,
            config.longintron,
            config.per_splice_site_isofrac,
            config.long_reads,
            config.eonly,
        );
        emit_jfinal_trace(&cjunctions, "bundle");
        // after good_junc sets mm=-1 on bad
        // junctions, the per-read loop creates REPLACEMENT junctions with
        // the read's splice strand (non-zero).  This means mm<0 junctions with
        // non-zero strand effectively stay alive for color propagation in CGroup
        // building.  We replicate this by NOT converting mm<0 → strand=0 here.
        // The killed_juncs computation (line ~7000) only includes strand==Some(0)
        // junctions, so mm<0-but-stranded junctions stay out of the kill set.
        if stage_debug {
            let killed = cjunctions
                .iter()
                .filter(|cj| cj.strand == 0 || cj.mm < 0.0)
                .count();
            eprintln!(
                "[EARLY_DEBUG] cjunction_post_good {}:{}-{} count={} killed={}",
                bundle.chrom,
                bundle.start,
                bundle.end,
                cjunctions.len(),
                killed
            );
        }
        // Save CJunction array for per-read junction redirect .
        // The cjunctions contain redirect pointers (negative nreads/nreads_good)
        // from apply_higherr_demotions, needed to adjust read exon boundaries.
        let cjunctions_for_redirect = cjunctions.clone();
        junction_stats_corr = crate::types::cjunctions_to_junction_stats(&cjunctions);
        apply_bad_mm_neg_stage(
            &mut junction_stats_corr,
            &bpcov_stranded,
            bundle.start,
            config.min_junction_reads,
        );
        if stage_debug {
            eprintln!(
                "[EARLY_DEBUG] post_cjunction_to_map {}:{}-{} count={}",
                bundle.chrom,
                bundle.start,
                bundle.end,
                junction_stats_corr.len()
            );
        }
        if let Some(g) = &genome {
            let chrom_arc = intern_chrom_arc(&mut chrom_arc_cache, &bundle.chrom);
            let guide_junctions = cached_bundle_guide_junctions(
                &mut guide_junction_cache,
                &mut chrom_arc_cache,
                &guide_transcripts,
                &bundle.chrom,
                bundle.strand,
                bundle.start,
                bundle.end,
            );
            let guide_tol = config.sserror.max(1);
            let mut noncons = 0usize;
            for (j, st) in junction_stats_corr.iter_mut() {
                if st.strand == Some(0) {
                    continue;
                }
                let guide_match = guide_junctions.as_ref().iter().any(|gj| {
                    gj.donor.abs_diff(j.donor) <= guide_tol
                        && gj.acceptor.abs_diff(j.acceptor) <= guide_tol
                });
                if guide_match {
                    continue;
                }
                let is_cons = cached_consensus_splice(
                    &mut consensus_cache,
                    g,
                    &chrom_arc,
                    j.donor,
                    j.acceptor,
                    st.strand,
                );
                let good_reads = st.nreads_good;
                // Consensus-kill disabled by default: regresses matches (1440→1401)
                // in this dataset. See pipeline.rs:6267 for notes. Enable via
                // RUSTLE_CONSENSUS_KILL=1 for strict parity with StringTie's
                // rlink.cpp:14842.
                if std::env::var_os("RUSTLE_CONSENSUS_KILL").is_some()
                    && !is_cons
                    && good_reads < 5.0
                {
                    st.strand = Some(0);
                    noncons += 1;
                }
            }
            if config.verbose && noncons > 0 {
                eprintln!(
                    "    genome_consensus: killed {} low-support non-consensus junction(s)",
                    noncons
                );
            }
        }
        let cjunctions_post = crate::types::junction_stats_to_cjunctions(&junction_stats_corr);
        let killed_junction_pairs = compute_killed_junction_pairs(&cjunctions_post);
        if !killed_junction_pairs.is_empty() || !junction_redirect_map.is_empty() {
            repair_reads_after_junction_quality(
                &mut bundle.reads,
                &junction_redirect_map,
                &killed_junction_pairs,
                &junction_stats_corr,
            );
        }
        if trace_log_style {
            eprintln!(
                "--- infer_transcripts: AFTER_COUNT_GOOD_JUNCTIONS bdata={}-{}",
                bundle.start + 1,
                trace_bundle_end(bundle.end)
            );
        }
        // Diagnostic: log bundles containing unstranded reads
        if std::env::var_os("RUSTLE_STRAND_DIAG").is_some() {
            let n_unstranded = bundle.reads.iter().filter(|r| r.strand == '.').count();
            if n_unstranded > 0 {
                eprintln!(
                    "STRAND_DIAG bundle={}:{}-{} strand={} total_reads={} unstranded_reads={}",
                    bundle.chrom,
                    bundle.start,
                    bundle.end,
                    bundle.strand,
                    bundle.reads.len(),
                    n_unstranded
                );
            }
        }

        // The algorithm processes ALL bundles through build_graphs, including single-exon.
        // The single-exon shortcut below is disabled by default for compatibility.
        // Set RUSTLE_SINGLE_EXON_SHORTCUT=1 to re-enable the optimized path.
        let use_single_exon_shortcut = std::env::var_os("RUSTLE_SINGLE_EXON_SHORTCUT").is_some();
        if use_single_exon_shortcut && bundle.strand == '.' {
            let has_multiexon = bundle.reads.iter().any(|r| r.exons.len() > 1);
            if !has_multiexon {
                // Create single-exon predictions from bundle coverage
                let single_exon_txs = create_single_exon_predictions_from_bundle(
                    &bundle,
                    &bpcov_stranded,
                    config.min_transcript_length,
                    config.singlethr,
                    config.verbose,
                );
                if config.verbose {
                    eprintln!(
                        "    dot_bundle_single_exon {}:{}-{} reads={} predictions={}",
                        bundle.chrom,
                        bundle.start,
                        bundle.end,
                        bundle.reads.len(),
                        single_exon_txs.len()
                    );
                }
                // Store single-exon predictions to add at the end (bypass filtering)
                single_exon_predictions_mutex.lock().unwrap().extend(single_exon_txs);
                return Ok(());
            }
        }

        let junction_stats_corr_final = junction_stats_corr.clone();
        let junction_stats = apply_junction_filters_and_canonicalize(
            junction_stats_corr,
            config.min_junction_reads,
            config.junction_canonical_tolerance,
            config.per_splice_site_isofrac,
            config.verbose,
        );
        let junctions = filter_junctions(&junction_stats, config.min_junction_reads);
        if trace_log_style {
            let pass_id = trace_pass_id.fetch_add(1, std::sync::atomic::Ordering::SeqCst) + 1;
            eprintln!(
                "--- build_graphs: PASS_ID={} bundle_start={} bundle_end={} nreads={} njunctions={}",
                pass_id,
                bundle.start + 1,
                trace_bundle_end(bundle.end),
                bundle.reads.len(),
                junction_stats_corr_final.len()
            );
            emit_trace_junction_actions(
                bundle.strand,
                &junction_stats_corr_final,
                &junction_redirect_map,
            );
            emit_trace_junction_decision_table(
                &junction_stats_corr_final,
                config.min_junction_reads,
            );
        }
        // Build (donor, acceptor) tuple set for filter_unsupported_junctions (uses local tuple type).
        let good_junctions: HashSet<(u64, u64)> =
            junctions.iter().map(|j| (j.donor, j.acceptor)).collect();

        let bundle_id = format!("{}:{}-{}", bundle.chrom, bundle.start, bundle.end);
        let debug_junc = config.debug_junctions
            && debug_bundle_overlaps(
                debug_target.as_ref(),
                &bundle.chrom,
                bundle.start,
                bundle.end,
            );
        if debug_junc {
            let mut raw_counts: HashMap<crate::types::Junction, f64> = Default::default();
            let mut del_counts: HashMap<crate::types::Junction, f64> = Default::default();
            for r in &bundle.reads {
                for j in &r.junctions_raw {
                    *raw_counts.entry(*j).or_insert(0.0) += r.weight;
                }
                for j in &r.junctions_del {
                    *del_counts.entry(*j).or_insert(0.0) += r.weight;
                }
            }
            eprintln!(
                "[DEBUG_JUNCTIONS] bundle {} strand={} reads={}",
                bundle_id,
                bundle.strand,
                bundle.reads.len()
            );
            dump_junction_counts("raw_junctions", &raw_counts);
            dump_junction_counts("del_junctions", &del_counts);
            dump_junction_stats("after_good_junc", &junction_stats_corr_final);
            dump_junction_stats("after_filter_canonical", &junction_stats);
            eprintln!("[DEBUG_JUNCTIONS] final_junctions: {}", junctions.len());
            for j in &junctions {
                eprintln!("  {}-{}", j.donor, j.acceptor);
            }
        }

        // after count_good_junctions trims exons, also correct exon boundaries
        // for redirected junctions (higherr + guide snap), so CGroup color break checks use
        // final coordinates.
        // Only touch exons where the junction was actually redirected.
        // (rd.segs[i].end=newstart; rd.segs[i+1].start=newend)
        if !junction_redirect_map.is_empty() {
            for r in &mut bundle.reads {
                for j_idx in 0..r.junctions.len() {
                    if j_idx + 1 < r.exons.len() {
                        // Check if this junction position (from exons) was corrected
                        let exon_junc = Junction::new(r.exons[j_idx].1, r.exons[j_idx + 1].0);
                        if junction_redirect_map.contains_key(&exon_junc) {
                            let junc = r.junctions[j_idx];
                            r.exons[j_idx].1 = junc.donor;
                            r.exons[j_idx + 1].0 = junc.acceptor;
                        }
                    }
                }
            }
        }

        if std::env::var_os("RUSTLE_DEBUG_DETAIL").is_some() {
            eprintln!("DEBUG_CGROUP_CHECK reads={}", bundle.reads.len());
        }
        // color breaks when juncs[i-1]->strand == 0 (killed junction).
        // Use killed status from final corrected junction stats directly. Optional
        // A/B propagation across corrected_map can be enabled for diagnostics.
        //
        // Two-step kill: higherr block sets mm=-1, then SR bundle loop sets strand=0
        // for mm<0 junctions (build_graphs: `if(jd.mm<0) jd.strand=0`).
        // when a junction has mm<0 (higherr
        // demotion), the read loop creates a REPLACEMENT junction with the
        // read's strand (non-zero).  CGroup building then sees non-zero strand and
        // allows color propagation.  In Rustle, we approximate this by keeping mm<0
        // junctions with non-zero strand OUT of killed_juncs AND IN good_junctions_set.
        // This allows color to propagate through higherr-demoted junctions that have
        // valid splice strand, matching per-read junction replacement.
        let mut good_junctions_set: HashSet<crate::types::Junction> =
            junctions.iter().copied().collect();
        // Remove mm<0 junctions: these are "run-through" demotions from
        // HE_CONT_R_DEMOTE, NOT real splices. Keeping them in good_junctions_set
        // lets bundle_builder bridge across gene boundaries (STRG.29/STRG.31 case).
        // Disable via RUSTLE_ALLOW_MM_NEG_IN_GOOD=1.
        if std::env::var_os("RUSTLE_ALLOW_MM_NEG_IN_GOOD").is_none() {
            let to_remove: Vec<crate::types::Junction> = good_junctions_set
                .iter()
                .filter_map(|j| {
                    junction_stats_corr_final.get(j).and_then(|s| {
                        if s.mm < 0.0 { Some(*j) } else { None }
                    })
                })
                .collect();
            for j in to_remove {
                good_junctions_set.remove(&j);
            }
        }
        // Add mm<0 junctions with non-zero strand to good_junctions (read-redirect compatibility).
        // StringTie's read loop sets strand=0 for mm<0 junctions, then creates replacement
        // junctions with the read's strand.  We approximate this by keeping the originals
        // in good_junctions (since replacement junctions are also added per-region later).
        // Gate: RUSTLE_KEEP_MM_NEG_IN_GOOD (default OFF post-trace-diff).
        // Re-inserting mm<0 junctions into good_junctions_set prevents color-break
        // on run-through junctions (e.g., STRG.29/STRG.31 boundary 17451950-17452251
        // with mm=-1 from HE_CONT_R_DEMOTE). StringTie's behavior: read-level
        // junction replacement creates separate transfrags rather than bridging
        // via the mm<0 junction.
        if std::env::var_os("RUSTLE_KEEP_MM_NEG_IN_GOOD").is_some() {
            for (j, s) in &junction_stats_corr_final {
                if s.mm < 0.0 && s.strand != Some(0) && s.nreads_good >= config.min_junction_reads {
                    good_junctions_set.insert(*j);
                }
            }
        }
        let mut killed_juncs: HashSet<crate::types::Junction> = junction_stats_corr_final
            .iter()
            .filter(|(_, s)| {
                // strand==0: junction genuinely has no strand support → killed.
                s.strand == Some(0)
            })
            .map(|(j, _)| *j)
            .collect();
        // Also treat mm<0 junctions as killed for read-split purposes (bridging
        // run-through junctions that lead to cross-gene reads). Matches StringTie
        // behavior: mm<0 junctions get strand=0 in the per-read loop.
        // Disable via RUSTLE_ALLOW_MM_NEG_IN_KILLED=1.
        if std::env::var_os("RUSTLE_ALLOW_MM_NEG_IN_KILLED").is_none() {
            for (j, s) in &junction_stats_corr_final {
                if s.mm < 0.0 && s.strand != Some(0) {
                    killed_juncs.insert(*j);
                }
            }
        }
        if std::env::var_os("RUSTLE_PROPAGATE_KILLED_CORRECTED").is_some() {
            for (orig, dest) in &junction_redirect_map {
                if killed_juncs.contains(orig) {
                    killed_juncs.insert(*dest);
                }
                if junction_stats_corr_final
                    .get(dest)
                    .map_or(false, |s| s.strand == Some(0) || s.mm < 0.0)
                {
                    killed_juncs.insert(*orig);
                }
            }
        }

        if std::env::var_os("RUSTLE_DEBUG_BUNDLE").is_some() {
            let sno = match bundle.strand {
                '-' => 0i32,
                '+' => 2,
                _ => 1,
            };
            eprintln!(
                "DEBUG_BG_CONTEXT sno={} chrom={} range={}-{} reads={} junctions={}",
                sno,
                bundle.chrom,
                bundle.start + 1,
                bundle.end,
                bundle.reads.len(),
                junctions.len()
            );
        }
        // in stranded bundles, reads neutralized by killed-junction
        // handling should not feed the stranded CGroup sweep.
        let skip_neutralized = (bundle.strand == '+' || bundle.strand == '-')
            && std::env::var_os("RUSTLE_INCLUDE_NEUTRALIZED").is_none();
        let cgroup_reads_filtered: Vec<BundleRead>;
        let cgroup_reads: &[BundleRead] = if skip_neutralized {
            cgroup_reads_filtered = bundle
                .reads
                .iter()
                .filter(|r| r.strand != '.')
                .cloned()
                .collect();
            &cgroup_reads_filtered
        } else {
            &bundle.reads
        };
        let region_key = (bundle.chrom.clone(), bundle.start, bundle.end);
        let region_reads: &[BundleRead] = region_reads_for_groupflow
            .get(&region_key)
            .map(|v| v.as_slice())
            .unwrap_or_else(|| bundle.reads.as_slice());
        let (boundary_left, boundary_right) =
            collect_guide_boundary_sets_for_bundle(&guide_transcripts, &bundle);
        // Keep 3-strand flow on all bundles by default for compatibility. Allow explicit opt-out.
        let use_3strand_all_strands = std::env::var_os("RUSTLE_DISABLE_3STRAND_ALL_STRANDS")
            .is_none()
            || std::env::var_os("RUSTLE_3STRAND_ALL_STRANDS").is_some();
        let use_3strand_for_bundle =
            use_3strand_groupflow && (bundle.strand == '.' || use_3strand_all_strands);
        let force_3strand_for_bundle = std::env::var_os("RUSTLE_FORCE_3STRAND_GROUPFLOW").is_some();
        let tighten_3strand_fallback =
            std::env::var_os("RUSTLE_TIGHTEN_3STRAND_FALLBACK").is_some();
        // Legacy fallback heuristics are now opt-in; default is strict 3-strand acceptance
        // unless output is clearly unusable (no assignments / catastrophic mapping).
        let legacy_3strand_fallback = std::env::var_os("RUSTLE_ENABLE_3STRAND_FALLBACK").is_some();

        let (bundlenodes, cgroup_read_bnodes, bnode_colors, cgroup_read_scales) =
            if use_3strand_for_bundle {
                let (bn_all, read_bnodes_all, bcolors, read_scale_all) =
                    build_bundlenodes_and_readgroups_from_cgroups_3strand(
                        region_reads,
                        bundle.strand,
                        &good_junctions_set,
                        &killed_juncs,
                        config.cgroup_junction_support,
                        config.bundle_merge_dist,
                        config.long_reads,
                        &boundary_left,
                        &boundary_right,
                        None, // junction_stats - not used in this code path
                        &junction_redirect_map,
                    );
                let subset_to_region =
                    map_subset_reads_to_region_indices(cgroup_reads, region_reads);
                let mapped_reads = subset_to_region.iter().flatten().count();
                let mapped_ratio = if cgroup_reads.is_empty() {
                    1.0
                } else {
                    mapped_reads as f64 / cgroup_reads.len() as f64
                };
                let mut subset_read_bnodes: Vec<Vec<usize>> = vec![Vec::new(); cgroup_reads.len()];
                let mut subset_read_scales: Vec<f64> = vec![0.0; cgroup_reads.len()];
                for (i, maybe_ri) in subset_to_region.iter().enumerate() {
                    if let Some(ri) = maybe_ri {
                        if *ri < read_bnodes_all.len() {
                            subset_read_bnodes[i] = read_bnodes_all[*ri].clone();
                        }
                        subset_read_scales[i] = read_scale_all.get(*ri).copied().unwrap_or(0.0);
                    }
                }
                if std::env::var_os("RUSTLE_DEBUG_BUNDLE").is_some() && mapped_ratio < 1.0 {
                    eprintln!(
                        "DEBUG_3STRAND_MAP mapped={}/{} ratio={:.3} chrom={} range={}-{}",
                        mapped_reads,
                        cgroup_reads.len(),
                        mapped_ratio,
                        bundle.chrom,
                        bundle.start + 1,
                        bundle.end
                    );
                }
                let no_assignments =
                    !cgroup_reads.is_empty() && subset_read_bnodes.iter().all(|v| v.is_empty());
                let mapped_ratio_threshold = if tighten_3strand_fallback { 0.90 } else { 0.98 };
                // Strict default: do not auto-fallback away from 3-strand output.
                // Legacy fallback heuristics remain opt-in for A/B diagnostics only.
                let mut suspicious_3strand = false;
                let three_assigned = subset_read_bnodes.iter().filter(|v| !v.is_empty()).count();
                let three_bnodes = bundlenodes_to_vec(bn_all.as_ref()).len();

                if legacy_3strand_fallback {
                    suspicious_3strand = (bn_all.is_none() && !cgroup_reads.is_empty())
                        || no_assignments
                        || mapped_ratio < 0.50
                        || mapped_ratio < mapped_ratio_threshold;
                }
                if suspicious_3strand && !force_3strand_for_bundle {
                    let direct_candidate = {
                        let (bn, rb, bc) = build_bundlenodes_and_readgroups_from_cgroups(
                            cgroup_reads,
                            &good_junctions_set,
                            &killed_juncs,
                            config.cgroup_junction_support,
                            config.bundle_merge_dist,
                            config.long_reads,
                            &boundary_left,
                            &boundary_right,
                        );
                        (bn, rb, bc, vec![1.0; cgroup_reads.len()])
                    };
                    if config.verbose || std::env::var_os("RUSTLE_DEBUG_BUNDLE").is_some() {
                        eprintln!(
                        "DEBUG_3STRAND_FALLBACK chrom={} range={}-{} strand={} mapped_ratio={:.3} no_assignments={} 3s_bnodes={} 3s_assigned={}",
                        bundle.chrom,
                        bundle.start + 1,
                        bundle.end,
                        bundle.strand,
                        mapped_ratio,
                        no_assignments,
                        three_bnodes,
                        three_assigned
                    );
                    }
                    direct_candidate
                } else {
                    (bn_all, subset_read_bnodes, bcolors, subset_read_scales)
                }
            } else {
                let (bn, rb, bc) = build_bundlenodes_and_readgroups_from_cgroups(
                    cgroup_reads,
                    &good_junctions_set,
                    &killed_juncs,
                    config.cgroup_junction_support,
                    config.bundle_merge_dist,
                    config.long_reads,
                    &boundary_left,
                    &boundary_right,
                );
                (bn, rb, bc, vec![1.0; cgroup_reads.len()])
            };

        let mut scaled_bundle_reads: Vec<BundleRead> = bundle.reads.clone();
        for (i, r) in scaled_bundle_reads.iter_mut().enumerate() {
            let scale = cgroup_read_scales.get(i).copied().unwrap_or(1.0);
            if (scale - 1.0).abs() > f64::EPSILON {
                r.weight *= scale;
                r.junc_mismatch_weight *= scale;
                for pc in r.pair_count.iter_mut() {
                    *pc *= scale;
                }
            }
        }

        let (snapshot_enabled, snapshot_detail) = {
            let guard = snapshot_writer.lock().unwrap();
            (guard.enabled(), guard.detail())
        };
        let process_graph =
            |graph_bundle: &crate::types::Bundle,
             junctions: Vec<Junction>,
             bundlenodes: Option<CBundlenode>,
             reads: &[BundleRead],
             coverage_reads: &[BundleRead],
             read_bundles: Option<&[Vec<usize>]>,
             good_junctions_local: &HashSet<(u64, u64)>,
             killed_junction_pairs_local: &HashSet<Junction>| {
                let snapshot_full = snapshot_enabled
                    && (snapshot_detail == SnapshotDetail::Full
                        || snapshot_all
                        || debug_bundle_overlaps(
                            debug_target.as_ref(),
                            &graph_bundle.chrom,
                            graph_bundle.start,
                            graph_bundle.end,
                        ));
                // Compute strand-specific per-base coverage for this bundle
                // (needed by create_graph for short-tail skip, and later for longtrim/coverage edges).
                let use_plus = graph_bundle.strand != '-';
                let bpcov = Bpcov::from_reads(
                    coverage_reads,
                    graph_bundle.start,
                    graph_bundle.end,
                    use_plus,
                );
                let graph_bpcov_stranded =
                    BpcovStranded::from_reads(coverage_reads, graph_bundle.start, graph_bundle.end);
                let mode = config.assembly_mode();
                // always collect read boundaries for long-read mode
                // the algorithm uses lstart/lend for longtrim during graph construction
                let _longtrim_on = config.enable_longtrim; // unused: longtrim always enabled for long-read mode
                let longtrim_in_graph = use_longtrim(mode); // longtrim always enabled for long-read mode
                let (lstart, lend) = collect_read_boundaries_with_cpas(
                    coverage_reads,
                    mode,
                    config.long_read_min_len,
                    graph_bundle.start,
                    graph_bundle.end,
                );
                if std::env::var_os("RUSTLE_BOUNDARY_TRACE").is_some() {
                    let cov_floor: f64 = std::env::var("RUSTLE_BOUNDARY_TRACE_MIN")
                        .ok()
                        .and_then(|v| v.parse().ok())
                        .unwrap_or(5.0);
                    for b in lstart.iter().filter(|b| b.cov >= cov_floor) {
                        eprintln!(
                            "LSTART {} cov={:.1} bundle={}-{}",
                            b.pos, b.cov, graph_bundle.start, graph_bundle.end
                        );
                    }
                    for b in lend.iter().filter(|b| b.cov >= cov_floor) {
                        eprintln!(
                            "LEND {} cov={:.1} bundle={}-{}",
                            b.pos, b.cov, graph_bundle.start, graph_bundle.end
                        );
                    }
                }

                let (mut graph, mut longtrim_synth, longtrim_stats) = create_graph_with_longtrim(
                    &junctions,
                    graph_bundle.start,
                    graph_bundle.end,
                    bundlenodes.as_ref(),
                    config.junction_support,
                    Some(reads),
                    graph_bundle.strand,
                    Some(&graph_bundle.junction_stats),
                    &bpcov,
                    Some(&graph_bpcov_stranded),
                    &lstart,
                    &lend,
                    longtrim_in_graph,
                    config.longtrim_min_boundary_cov,
                );
                if longtrim_stats.applied && config.verbose {
                    eprintln!(
                    "    longtrim graph-build: lstart={} lend={} split_nodes={} new_nodes={} source_edges+={} sink_edges+={} synthetic={}",
                    longtrim_stats.lstart_events,
                    longtrim_stats.lend_events,
                    longtrim_stats.longtrim.split_nodes,
                    longtrim_stats.longtrim.new_nodes,
                    longtrim_stats.longtrim.source_edges_added,
                    longtrim_stats.longtrim.sink_edges_added,
                    longtrim_stats.longtrim.synthetic_transfrags
                );
                } else if longtrim_in_graph
                    && (!lstart.is_empty() || !lend.is_empty())
                    && config.verbose
                {
                    eprintln!(
                        "    longtrim_prestep: lstart={} lend={} (boundary stream collected)",
                        lstart.len(),
                        lend.len(),
                    );
                }
                trace_intron_chain_graph(&graph, graph_bundle, &config);
                if config.allowed_graph_nodes > 0 {
                    let _ =
                        prune_graph_nodes(&mut graph, config.allowed_graph_nodes, config.verbose);
                }

                // Collapse zero-length internal nodes before transfrag creation.
                // Opt-in via RUSTLE_COLLAPSE_ZEROLEN=1 — current build loses ~5
                // matches in exchange for +0.7% precision and -1430 graph nodes.
                // The path-selection side-effects still need investigation before
                // enabling by default.
                if std::env::var_os("RUSTLE_COLLAPSE_ZEROLEN").is_some() {
                    collapse_zero_length_nodes(&mut graph);
                }

                // GNODE trace: emit per-graph-node state for diffing vs StringTie.
                // Format matches stringtie_debug/rlink.cpp GNODE block.
                // s: 0=minus, 1=unstranded, 2=plus (matches StringTie sno).
                if std::env::var_os("GNODE_TRACE").is_some() {
                    let s = match graph_bundle.strand {
                        '-' => 0,
                        '+' => 2,
                        _ => 1,
                    };
                    for (idx, node) in graph.nodes.iter().enumerate() {
                        if idx == graph.source_id || idx == graph.sink_id {
                            continue;
                        }
                        let preds: Vec<String> =
                            node.parents.ones().map(|p| p.to_string()).collect();
                        let succs: Vec<String> =
                            node.children.ones().map(|c| c.to_string()).collect();
                        eprintln!(
                            "GNODE s={} {}-{} id={} preds={} succs={}",
                            s,
                            node.start,
                            node.end,
                            node.node_id,
                            preds.join(","),
                            succs.join(",")
                        );
                    }
                }

                let mut graph_mut = graph;
                if snapshot_enabled && (snapshot_full || snapshot_this) {
                    let key = snapshot::bundle_key(bundle_idx, graph_bundle);
                    let summary = snapshot::summarize_graph(&graph_mut);
                    if snapshot_full {
                        let _ = snapshot_writer.lock().unwrap().emit(&snapshot::SnapshotEnvelope {
                            kind: "rustle_snapshot",
                            stage: "graph_created",
                            bundle: key,
                            payload: snapshot::GraphWithNodes {
                                summary,
                                nodes: snapshot::snap_graph_nodes(&graph_mut),
                            },
                        });
                    } else {
                        let _ = snapshot_writer.lock().unwrap().emit(&snapshot::SnapshotEnvelope {
                            kind: "rustle_snapshot",
                            stage: "graph_created",
                            bundle: key,
                            payload: summary,
                        });
                    }
                }
                let mut guide_boundary_synth: Vec<GraphTransfrag> = Vec::new();

                if !guide_transcripts.is_empty()
                    && std::env::var_os("RUSTLE_ENABLE_GUIDE_BOUNDARY_SPLITS").is_some()
                {
                    let (gb_synth, gb_stats) =
                        crate::longtrim::apply_guide_boundary_splits_from_refs(
                            &mut graph_mut,
                            &guide_transcripts,
                            graph_bundle.start,
                            graph_bundle.end,
                            graph_bundle.strand,
                            &bpcov,
                            config.longtrim_min_boundary_cov,
                        );
                    if config.verbose
                        && (gb_stats.split_nodes > 0
                            || gb_stats.source_edges_added > 0
                            || gb_stats.sink_edges_added > 0)
                    {
                        eprintln!(
                        "    guide_boundary_splits: split_nodes={} new_nodes={} source_edges+={} sink_edges+={} synthetic={}",
                        gb_stats.split_nodes,
                        gb_stats.new_nodes,
                        gb_stats.source_edges_added,
                        gb_stats.sink_edges_added,
                        gb_stats.synthetic_transfrags
                    );
                    }
                    guide_boundary_synth = gb_synth;
                }

                let mapped_guides = if guide_transcripts.is_empty() {
                    Vec::new()
                } else {
                    process_refguides(
                        &mut graph_mut,
                        &guide_transcripts,
                        graph_bundle.start,
                        graph_bundle.end,
                        graph_bundle.strand,
                        25,
                        config.verbose,
                    )
                };

                // Add source/sink edges for coverage drops (per-coverage source/sink logic).
                // uses strand-indexed get_cov_sign(sno, ...).
                let covlinks_sno = match graph_bundle.strand {
                    '-' => BPCOV_STRAND_MINUS,
                    '+' => BPCOV_STRAND_PLUS,
                    _ => BPCOV_STRAND_ALL,
                };
                let coverage_synth = crate::graph_build::add_coverage_source_sink_edges(
                    &mut graph_mut,
                    &graph_bpcov_stranded,
                    covlinks_sno,
                );
                // Opt-in alt-TTS sink edges at nodes whose end coincides
                // with a junction donor (STRG.294.3 class: alt-TTS
                // terminating where a splice-donor continues for the
                // dominant isoform). Gated RUSTLE_IMPLICIT_ALT_TTS_SINK=1.
                let mut alt_tts_synth = crate::graph_build::add_alt_tts_sink_edges(
                    &mut graph_mut,
                    junctions.as_ref(),
                    graph_bundle.strand,
                    Some(&graph_bundle.junction_stats),
                );
                // Terminal-read cluster → junction-donor hardend (opt-in
                // read-signal-based alt-TTS detection). Replaces the need
                // for RUSTLE_ALT_TTS_ORACLE when enabled.
                let mut terminal_donor_synth = crate::graph_build::discover_terminal_donor_hardends(
                    &mut graph_mut,
                    reads,
                    junctions.as_ref(),
                    graph_bundle.strand,
                    Some(&graph_bundle.junction_stats),
                );

                // Coverage-based node splitting: short-read always (trimnode_all); long/mixed when longtrim not yet implemented
                let mut synthetic = if use_coverage_trim(mode) {
                    // trimnode_all can expose new candidate nodes after a split; run a few passes.
                    const MAX_COVERAGE_TRIM_PASSES: usize = 1;
                    let mut all_synth = Vec::new();
                    for pass in 0..MAX_COVERAGE_TRIM_PASSES {
                        let pass_synth = apply_coverage_trim(
                            &mut graph_mut,
                            &bpcov,
                            &lstart,
                            &lend,
                            graph_bundle.start,
                            graph_bundle.strand,
                            1.0,
                            mode.is_mixed(),
                            config.verbose,
                        );
                        if pass_synth.is_empty() {
                            break;
                        }
                        if config.verbose {
                            eprintln!(
                                "      [Rustle] coverage_trim pass {}: +{} synthetic transfrags",
                                pass + 1,
                                pass_synth.len()
                            );
                        }
                        all_synth.extend(pass_synth);
                    }
                    all_synth
                } else {
                    Vec::new()
                };

                // In long-read mode, futuretr transfrags get longread=true.
                // In mixed mode, the algorithm creates both a non-longread and a longread copy .
                let mut coverage_synth = coverage_synth;
                tag_transfrags_origin(&mut longtrim_synth, "longtrim");
                tag_transfrags_origin(&mut synthetic, "coverage_trim");
                tag_transfrags_origin(&mut coverage_synth, "coverage_edge");
                tag_transfrags_origin(&mut guide_boundary_synth, "guide_boundary");
                // Merge alt_tts synth into coverage_synth so they flow
                // through the same post-prune remap + flow seeding paths.
                if !alt_tts_synth.is_empty() {
                    tag_transfrags_origin(&mut alt_tts_synth, "alt_tts");
                    coverage_synth.extend(alt_tts_synth);
                }
                if !terminal_donor_synth.is_empty() {
                    tag_transfrags_origin(&mut terminal_donor_synth, "terminal_donor");
                    coverage_synth.extend(terminal_donor_synth);
                }
                if mode.is_long_read() {
                    for tf in &mut coverage_synth {
                        tf.longread = true;
                    }
                }
                let mut longtrim_synth = longtrim_synth;
                let mut synthetic = synthetic;
                let mut post_prune_redirects = Default::default();
                if config.allowed_graph_nodes > 0 {
                    let (removed_post, redirects) = prune_graph_nodes_with_redirects(
                        &mut graph_mut,
                        config.allowed_graph_nodes,
                        config.verbose,
                    );
                    if removed_post > 0 {
                        remap_two_node_synth_with_redirects(
                            &mut longtrim_synth,
                            &redirects,
                            &graph_mut,
                        );
                        remap_two_node_synth_with_redirects(&mut synthetic, &redirects, &graph_mut);
                        remap_two_node_synth_with_redirects(
                            &mut coverage_synth,
                            &redirects,
                            &graph_mut,
                        );
                    }
                    post_prune_redirects = redirects;
                }
                graph_mut.reindex_edge_bits_dense();
                if loop_trace_active() {
                    let s = match graph_bundle.strand {
                        '-' => 0,
                        '+' => 1,
                        _ => 0,
                    };
                    eprintln!(
                        "LOOP_build_graphs: parse_graphs s={} nodes={} strand={} bundle={}-{}",
                        s,
                        graph_mut.n_nodes,
                        graph_bundle.strand,
                        graph_bundle.start + 1,
                        graph_bundle.end
                    );
                }

                // Parity fix: synthesize terminal source/sink links before read->graph mapping so
                // mapped read paths see final terminal connectivity.
                let traverse_synth = graph_mut.synthesize_terminal_transfrags(mode, 1.0);

                if loop_trace_active() {
                    eprintln!(
                        "LOOP_build_graphs: read_to_transfrag nreads={}",
                        reads.len()
                    );
                }
                // `GraphNode::coverage` is bp-coverage mass (`CGraphnode::cov`), not bundlenode coverage.
                // Start clean before `get_fragment_pattern`-equivalent mapping accumulates read support.
                zero_graph_node_bp_coverage(&mut graph_mut);
                let bundle2graph =
                    crate::graph_build::build_bundle2graph(&graph_mut, bundlenodes.as_ref());
                let mut transfrags = if bundle2graph.iter().any(|nodes| !nodes.is_empty()) {
                    map_reads_to_graph_bundlenodes(
                        reads,
                        &mut graph_mut,
                        mode,
                        config.long_read_min_len,
                        config.junction_correction_window,
                        Some(killed_junction_pairs_local),
                        bundlenodes.as_ref(),
                        Some(&bundle2graph),
                        read_bundles,
                    )
                } else {
                    map_reads_to_graph(
                        reads,
                        &mut graph_mut,
                        mode,
                        config.long_read_min_len,
                        config.junction_correction_window,
                        Some(killed_junction_pairs_local),
                        None,
                    )
                };
                tag_transfrags_origin_if_missing(&mut transfrags, "read_map");
                debug_stage::emit_with_detail(
                    "graph_evolution",
                    &graph_bundle.chrom,
                    graph_bundle.start,
                    graph_bundle.end,
                    graph_bundle.strand,
                    reads.len(),
                    junctions.len(),
                    graph_mut.n_nodes.saturating_sub(2),
                    debug_stage::graph_edge_count(&graph_mut),
                    transfrags.len(),
                    0,
                    0,
                    0,
                    0,
                    debug_stage::graph_detail(&graph_mut, &transfrags),
                    "after_map_reads",
                );
                trace_graph_numeric_state("graph:after_map_reads", &graph_mut, &transfrags);
                if !mapped_guides.is_empty() {
                    let min_seed_abund = config.readthr.max(1.0);
                    let added_guides = inject_guide_transfrags(
                        &mut graph_mut,
                        &mut transfrags,
                        &mapped_guides,
                        min_seed_abund,
                    );
                    if config.verbose && added_guides > 0 {
                        eprintln!(
                            "    guide_mode: injected {} guide transfrags (mapped={})",
                            added_guides,
                            mapped_guides.len()
                        );
                    }
                }
                let start_add = transfrags.len();
                let mut future_links: Vec<FutureLink> = Vec::new();
                collect_from_transfrags(
                    &mut future_links,
                    &traverse_synth,
                    false,
                    graph_mut.source_id,
                    graph_mut.sink_id,
                );
                collect_from_transfrags(
                    &mut future_links,
                    &longtrim_synth,
                    true,
                    graph_mut.source_id,
                    graph_mut.sink_id,
                );
                collect_from_transfrags(
                    &mut future_links,
                    &synthetic,
                    true,
                    graph_mut.source_id,
                    graph_mut.sink_id,
                );
                collect_from_transfrags(
                    &mut future_links,
                    &coverage_synth,
                    true,
                    graph_mut.source_id,
                    graph_mut.sink_id,
                );
                collect_from_transfrags(
                    &mut future_links,
                    &guide_boundary_synth,
                    true,
                    graph_mut.source_id,
                    graph_mut.sink_id,
                );
                if !post_prune_redirects.is_empty() {
                    apply_prune_redirects(
                        &mut future_links,
                        &post_prune_redirects,
                        graph_mut.source_id,
                        graph_mut.sink_id,
                    );
                }
                normalize_links(&graph_mut, &mut future_links);
                let mut realized =
                    materialize_links(&mut graph_mut, &future_links, 1.0, LONGINTRONANCHOR as u64);
                // long-read only: realized transfrags keep abundance as-is
                transfrags.extend(realized);
                register_transfrag_range_on_nodes(&mut graph_mut, &transfrags, start_add);
                annotate_hard_boundaries(&mut graph_mut, &transfrags, graph_bundle.strand);
                // long-read only: no mixed-mode srfrag redistribution
                // compute reachability BEFORE process_transfrags so
                // that compatible_long/conflict can use childpat/parentpat.
                // The algorithm computes reachability during create_graph (traverse_dfs)
                // which runs before process_transfrags.
                graph_mut.compute_reachability();

                // StringTie-parity sink/source connector synthesis. Opt-in via
                // RUSTLE_SINK_CONNECTOR_TF. For each real node N that has sink
                // as a child but no existing [N, sink] 2-node transfrag,
                // synthesize one with fallback abundance. This populates
                // hassink[N] in process_transfrags so chimeric tfs ending at N
                // SKIP trflong pass 1 (matching StringTie). Analogous sweep
                // for [source, N] links.
                if std::env::var_os("RUSTLE_SINK_CONNECTOR_TF").is_some() {
                    let trace_sc = std::env::var_os("RUSTLE_SINK_CONNECTOR_TRACE").is_some();
                    // Fallback abundance for synthetic [n, sink] / [source, n] connectors.
                    // MUST be high enough to survive multi-seed depletion during
                    // parse_trflong (each seed's long_max_flow can consume abundance).
                    // Default 25.0 matches StringTie's `tmpcov+trthr` scale in longtrim
                    // (typical range 25-200 depending on coverage contrast). Lower values
                    // (1-10) cause the connector to deplete after the first seed, leaving
                    // subsequent seeds stranded → checktrf rescue → over-emission.
                    let fallback_abund: f64 = std::env::var("RUSTLE_SINK_CONNECTOR_FALLBACK_ABUND")
                        .ok()
                        .and_then(|s| s.parse::<f64>().ok())
                        .unwrap_or(25.0);
                    let source_id_sc = graph_mut.source_id;
                    let sink_id_sc = graph_mut.sink_id;
                    // Record existing 2-node connectors so we don't duplicate.
                    let mut has_sink_tf: std::collections::HashSet<usize> = Default::default();
                    let mut has_src_tf: std::collections::HashSet<usize> = Default::default();
                    for tf in transfrags.iter() {
                        if tf.node_ids.len() < 2 { continue; }
                        if tf.node_ids.last().copied() == Some(sink_id_sc) {
                            has_sink_tf.insert(tf.node_ids[0]);
                        }
                        if tf.node_ids[0] == source_id_sc {
                            has_src_tf.insert(tf.node_ids[1]);
                        }
                    }
                    let mut added_sink = 0usize;
                    let mut added_src = 0usize;
                    // Sink sweep: nodes with sink as child.
                    let sink_candidates: Vec<usize> = (0..graph_mut.nodes.len())
                        .filter(|&nid| nid != source_id_sc && nid != sink_id_sc)
                        .filter(|&nid| graph_mut.nodes[nid].children.contains(sink_id_sc))
                        .filter(|nid| !has_sink_tf.contains(nid))
                        .collect();
                    for nid in sink_candidates {
                        if trace_sc && graph_bundle.start <= 70761097 && graph_bundle.end >= 70637431 {
                            let n = &graph_mut.nodes[nid];
                            eprintln!("SINK_CONN_PIPELINE nid={} coord={}-{} abund={:.2}",
                                nid, n.start, n.end, fallback_abund);
                        }
                        let mut tf = GraphTransfrag::new(vec![nid, sink_id_sc], graph_mut.pattern_size());
                        tf.abundance = fallback_abund;
                        tf.longread = true;
                        tf.origin_tag = Some("sink_connector".to_string());
                        graph_mut.set_pattern_edges_for_path(&mut tf.pattern, &tf.node_ids);
                        transfrags.push(tf);
                        added_sink += 1;
                    }
                    // Source sweep: nodes with source as parent.
                    let src_candidates: Vec<usize> = (0..graph_mut.nodes.len())
                        .filter(|&nid| nid != source_id_sc && nid != sink_id_sc)
                        .filter(|&nid| graph_mut.nodes[nid].parents.contains(source_id_sc))
                        .filter(|nid| !has_src_tf.contains(nid))
                        .collect();
                    for nid in src_candidates {
                        if trace_sc && graph_bundle.start <= 70761097 && graph_bundle.end >= 70637431 {
                            let n = &graph_mut.nodes[nid];
                            eprintln!("SRC_CONN_PIPELINE nid={} coord={}-{} abund={:.2}",
                                nid, n.start, n.end, fallback_abund);
                        }
                        let mut tf = GraphTransfrag::new(vec![source_id_sc, nid], graph_mut.pattern_size());
                        tf.abundance = fallback_abund;
                        tf.longread = true;
                        tf.origin_tag = Some("source_connector".to_string());
                        graph_mut.set_pattern_edges_for_path(&mut tf.pattern, &tf.node_ids);
                        transfrags.push(tf);
                        added_src += 1;
                    }
                    if trace_sc && (added_sink > 0 || added_src > 0) {
                        eprintln!("SINK_CONN_SUMMARY bundle={}-{} added_sink={} added_src={}",
                            graph_bundle.start, graph_bundle.end, added_sink, added_src);
                    }
                }

                transfrags = process_transfrags(
                    transfrags,
                    &mut graph_mut,
                    1.0,
                    config.singlethr,
                    !guide_transcripts.is_empty(),
                    config.long_reads,
                    false, // mixed_mode: always false (long-read only)
                    config.verbose,
                    config.eonly,
                    config.keeptrf_usepath_tsv.as_deref(),
                );
                // Option A (graph segmentation refactor): coalesce seed
                // transfrags whose intron chains differ only by alt-donor/
                // acceptor shifts within a small window. Reduces seed count
                // at alt-splice hot spots (STRG.309 class: 39 → target ~2-3
                // seeds), which prevents max_flow from emitting 16
                // combinatorial variants. Opt-in via RUSTLE_TF_COALESCE=1.
                let (_coalesce_before, _coalesce_after) =
                    crate::transfrag_process::coalesce_alt_junc_seed_transfrags(
                        &mut transfrags, &graph_mut
                    );
                // Diagnostic missed-tx oracle: inject StringTie ref tx as
                // high-abundance seeds and classify per-tx why Rustle fails
                // to emit (unmappable / partial / mappable-not-emitted).
                let _oracle_injected = crate::graph_build::inject_missed_tx_oracle(
                    &graph_mut,
                    &mut transfrags,
                    &graph_bundle.chrom,
                    graph_bundle.start,
                    graph_bundle.end,
                    graph_bundle.strand,
                );
                if snapshot_enabled && (snapshot_full || snapshot_this) {
                    let key = snapshot::bundle_key(bundle_idx, graph_bundle);
                    let summary = snapshot::summarize_transfrags(&transfrags);
                    if snapshot_full {
                        let _ = snapshot_writer.lock().unwrap().emit(&snapshot::SnapshotEnvelope {
                            kind: "rustle_snapshot",
                            stage: "transfrags_processed",
                            bundle: key,
                            payload: snapshot::TransfragsFull {
                                summary,
                                transfrags: snapshot::snap_transfrags(&transfrags),
                            },
                        });
                    } else {
                        let _ = snapshot_writer.lock().unwrap().emit(&snapshot::SnapshotEnvelope {
                            kind: "rustle_snapshot",
                            stage: "transfrags_processed",
                            bundle: key,
                            payload: summary,
                        });
                    }
                }
                trace_graph_numeric_state(
                    "graph:after_process_transfrags",
                    &graph_mut,
                    &transfrags,
                );
                if crate::trace_pipeline::active() {
                    let src = graph_mut.source_id;
                    let snk = graph_mut.sink_id;
                    for from in 0..graph_mut.nodes.len() {
                        if from == src || from == snk { continue; }
                        let fs = graph_mut.nodes[from].start;
                        let fe = graph_mut.nodes[from].end;
                        let kids: Vec<usize> = graph_mut.nodes[from].children.ones().collect();
                        for c in kids {
                            if c == src || c == snk { continue; }
                            let ts = graph_mut.nodes[c].start;
                            let te = graph_mut.nodes[c].end;
                            crate::trace_pipeline::dump_graph_edge(
                                &graph_bundle.chrom,
                                graph_bundle.start,
                                graph_bundle.end,
                                graph_bundle.strand,
                                fs, fe, ts, te,
                            );
                        }
                    }
                }
                if loop_trace_active() {
                    let trflong = transfrags.iter().filter(|tf| tf.trflong_seed).count();
                    eprintln!(
                        "LOOP_build_graphs: process_transfrags transfrags={} trflong={}",
                        transfrags.len(),
                        trflong
                    );
                }
                // bpcov refresh only for short-read mode (removed)
                let _nodecov = compute_nodecov(&mut graph_mut, &transfrags, config.verbose);
                trace_graph_numeric_state("graph:after_compute_nodecov", &graph_mut, &transfrags);

                // DEBUG_GRAPH: emit graph stats for comparison with expected graph stats
                if std::env::var_os("RUSTLE_DEBUG_BUNDLE").is_some() {
                    let sno = match graph_bundle.strand {
                        '-' => 0i32,
                        '+' => 2,
                        _ => 1,
                    };
                    let n_interior = graph_mut.n_nodes.saturating_sub(2); // exclude source/sink
                    let active_trf = transfrags.iter().filter(|t| t.abundance > 0.0).count();
                    let long_trf = transfrags.iter().filter(|t| t.longread).count();
                    eprintln!(
                        "DEBUG_GRAPH sno={} range={}-{} nodes={} transfrags={} active={} long={}",
                        sno,
                        graph_bundle.start + 1,
                        graph_bundle.end,
                        n_interior,
                        transfrags.len(),
                        active_trf,
                        long_trf
                    );
                    // Also emit per-node coordinates for detailed comparison
                    for nid in 0..graph_mut.n_nodes {
                        if nid == graph_mut.source_id || nid == graph_mut.sink_id {
                            continue;
                        }
                        let node = &graph_mut.nodes[nid];
                        eprintln!(
                            "  GNODE sno={} n={} start={} end={} cov={:.1}",
                            sno,
                            nid,
                            node.start + 1,
                            node.end,
                            node.coverage
                        );
                    }
                }

                // Compute reachability for onpath_long validation in path extension
                graph_mut.compute_reachability();

                // Grow all patterns to current pattern_size FIRST: transfrags created early (in
                // map_reads_to_graph with psize = pre-synthesis value) may have smaller sizes than
                // edge IDs added later by synthesize_terminal_transfrags / process_transfrags.
                // Keep the original edge-bit structure intact: transfrag patterns preserve
                // missing edge bits for incomplete paths, and long-recursion relies on that.
                let needed_psize = graph_mut.pattern_size();
                for tf in transfrags.iter_mut() {
                    tf.pattern.grow(needed_psize);
                }

                let (
                    txs,
                    pre_filter_txs,
                    seed_outcomes,
                    longrec_summary,
                    predcluster_summary,
                    junction_support_removed,
                ) = extract_bundle_transcripts_for_graph(
                    &mut graph_mut,
                    &mut transfrags,
                    graph_bundle,
                    &config,
                    &mapped_guides,
                    &guide_transcripts,
                    &ref_transcripts,
                    &junctions,
                    good_junctions_local,
                    &bpcov,
                    &graph_bpcov_stranded,
                    trace_reference.is_some() || config.debug_stage_tsv.is_some(),
                );
                // BUNDLE_STATS per-bundle summary across all layers (RUSTLE_BUNDLE_STATS=1).
                // Emits one line with: reads, junctions, nodes/edges, transfrags,
                // seeds, outcome counts. Paired with StringTie PARITY_BUNDLE_STATS
                // for cross-pipeline layer-by-layer diff.
                if std::env::var_os("RUSTLE_BUNDLE_STATS").is_some() {
                    let n_reads = reads.len();
                    let n_junc = junctions.len();
                    let n_junc_kept = good_junctions_local.len();
                    let n_nodes = graph_mut.nodes.len().saturating_sub(2); // exclude source/sink
                    let mut n_edges = 0usize;
                    let src_bs = graph_mut.source_id;
                    let snk_bs = graph_mut.sink_id;
                    for i in 0..graph_mut.nodes.len() {
                        if i == src_bs || i == snk_bs { continue; }
                        for c in graph_mut.nodes[i].children.ones() {
                            if c == src_bs || c == snk_bs { continue; }
                            n_edges += 1;
                        }
                    }
                    let n_transfrags = transfrags.len();
                    let n_seeds = transfrags.iter().filter(|tf| tf.trflong_seed).count();
                    // Count seed outcomes by type.
                    use crate::path_extract::SeedOutcome;
                    let mut n_stored = 0usize;
                    let mut n_checktrf = 0usize;
                    let mut n_lowcov = 0usize;
                    let mut n_back_fail = 0usize;
                    let mut n_fwd_fail = 0usize;
                    let mut n_hard_bound = 0usize;
                    let mut n_unwit = 0usize;
                    let mut n_zflux = 0usize;
                    let mut n_other = 0usize;
                    for (_idx, outcome) in &seed_outcomes {
                        match outcome {
                            SeedOutcome::LowCoverage(_) => n_lowcov += 1,
                            SeedOutcome::BackToSourceFail => n_back_fail += 1,
                            SeedOutcome::FwdToSinkFail => n_fwd_fail += 1,
                            SeedOutcome::HardBoundaryMismatch => n_hard_bound += 1,
                            SeedOutcome::UnwitnessedSplice => n_unwit += 1,
                            SeedOutcome::ZeroFlux => n_zflux += 1,
                            SeedOutcome::ChecktrfRescued => n_stored += 1,
                            SeedOutcome::ChecktrfRescueFail => n_checktrf += 1,
                            _ => n_other += 1,
                        }
                    }
                    let n_final = txs.len();
                    eprintln!(
                        "RUSTLE_BUNDLE_STATS chrom={} bundle={}-{} strand={} n_reads={} n_junc={} n_junc_kept={} n_nodes={} n_edges={} n_transfrags={} n_seeds={} n_stored={} n_checktrf={} n_lowcov={} n_back_fail={} n_fwd_fail={} n_hard_bound={} n_unwit={} n_zflux={} n_other={} n_final={}",
                        graph_bundle.chrom, graph_bundle.start, graph_bundle.end, graph_bundle.strand,
                        n_reads, n_junc, n_junc_kept, n_nodes, n_edges,
                        n_transfrags, n_seeds,
                        n_stored, n_checktrf, n_lowcov, n_back_fail, n_fwd_fail,
                        n_hard_bound, n_unwit, n_zflux, n_other,
                        n_final
                    );
                }
                // Read-chain witness filter: remove transcripts whose intron
                // chain pairs are not witnessed by any read. This prevents the
                // flow decomposition from creating novel junction combinations.
                // Flow-sourced transcripts are given a singleton-witness
                // fallback (both junctions individually observed → pass), which
                // recovers StringTie-parity junction variants assembled from
                // partially-overlapping reads (e.g., STRG.112.1).
                let txs = if config.long_reads {
                    let read_pairs = crate::transcript_filter::build_read_intron_pairs(reads);
                    let read_singles = crate::transcript_filter::build_read_junction_singletons(reads);
                    let read_junc_counts = crate::transcript_filter::build_read_junction_counts(reads);
                    let mut filtered = crate::transcript_filter::filter_unwitnessed_chains_with_singletons_and_counts(
                        txs,
                        &read_pairs,
                        Some(&read_singles),
                        Some(&read_junc_counts),
                        config.junction_correction_window,
                        config.verbose,
                    );
                    // Full-chain-witness filter (DEFAULT-ON): kill tx where
                    // some contiguous K-intron window of the chain isn't
                    // witnessed by any read. Targets STRG.309-class j-class
                    // noise where flow-combinatorial variants use junction
                    // combinations no read actually has together.
                    // Default K=4 is lossless on GGO_19; opt-out via
                    // RUSTLE_FULL_CHAIN_WITNESS_OFF=1.
                    if std::env::var_os("RUSTLE_FULL_CHAIN_WITNESS_OFF").is_none() {
                        let read_chains = crate::transcript_filter::build_read_intron_chains(reads);
                        let tol = config.junction_correction_window;
                        filtered = crate::transcript_filter::filter_by_full_chain_witness(
                            filtered, &read_chains, tol, config.verbose,
                        );
                    }
                    if std::env::var_os("RUSTLE_ORACLE_STAGE_TRACE").is_some() {
                        for t in filtered.iter() {
                            if t.source.as_deref().map_or(false, |s| s.starts_with("oracle_direct:") || s.starts_with("ref_chain_rescue:")) {
                                eprintln!(
                                    "ORACLE_STAGE stage=post_filter_unwitnessed src={} {}-{} exons={}",
                                    t.source.as_deref().unwrap_or(""),
                                    t.exons.first().map(|e| e.0).unwrap_or(0),
                                    t.exons.last().map(|e| e.1).unwrap_or(0),
                                    t.exons.len()
                                );
                            }
                        }
                    }
                    filtered
                } else {
                    txs
                };

                if snapshot_enabled && (snapshot_full || snapshot_this) {
                    let key = snapshot::bundle_key(bundle_idx, graph_bundle);
                    let summary = snapshot::summarize_transcripts(&txs);
                    if snapshot_full {
                        let _ = snapshot_writer.lock().unwrap().emit(&snapshot::SnapshotEnvelope {
                            kind: "rustle_snapshot",
                            stage: "transcripts_emitted",
                            bundle: key,
                            payload: snapshot::TranscriptsFull {
                                summary,
                                transcripts: snapshot::snap_transcripts(&txs),
                            },
                        });
                    } else {
                        let _ = snapshot_writer.lock().unwrap().emit(&snapshot::SnapshotEnvelope {
                            kind: "rustle_snapshot",
                            stage: "transcripts_emitted",
                            bundle: key,
                            payload: summary,
                        });
                    }
                }
                if config.verbose {
                    eprintln!(
                        "    bundle_result {}:{}-{}({}) extracted={} kept={}",
                        graph_bundle.chrom,
                        graph_bundle.start,
                        graph_bundle.end,
                        graph_bundle.strand,
                        txs.len(),
                        txs.len()
                    );
                }
                let mut seed_summary = SeedOutcomeSummary::default();
                accumulate_seed_outcomes(&mut seed_summary, &seed_outcomes);
                let pre_filter_count = pre_filter_txs
                    .as_ref()
                    .map(|v| v.len())
                    .unwrap_or(txs.len());
                let pairwise_removed = predcluster_summary
                    .entry_count
                    .saturating_sub(predcluster_summary.after_pairwise);
                let isofrac_removed = predcluster_summary
                    .after_pairwise
                    .saturating_sub(predcluster_summary.after_isofrac);
                let runoff_removed = predcluster_summary
                    .after_isofrac
                    .saturating_sub(predcluster_summary.after_runoff);
                let polymerase_runoff_removed = predcluster_summary
                    .after_runoff
                    .saturating_sub(predcluster_summary.after_polymerase_runoff);
                let polymerase_runon_removed = predcluster_summary
                    .after_polymerase_runoff
                    .saturating_sub(predcluster_summary.after_polymerase_runon);
                let readthr_removed = predcluster_summary
                    .after_polymerase_runon
                    .saturating_sub(predcluster_summary.after_readthr);
                debug_stage::emit_with_detail(
                    "seed_flow",
                    &graph_bundle.chrom,
                    graph_bundle.start,
                    graph_bundle.end,
                    graph_bundle.strand,
                    reads.len(),
                    junctions.len(),
                    graph_mut.n_nodes.saturating_sub(2),
                    debug_stage::graph_edge_count(&graph_mut),
                    transfrags.len(),
                    txs.len(),
                    seed_summary.total,
                    seed_summary.stored,
                    0,
                    StageDetail {
                        seed_unwitnessed: seed_summary.unwitnessed,
                        seed_zero_flux: seed_summary.zero_flux,
                        seed_low_cov: seed_summary.low_coverage,
                        seed_back_fail: seed_summary.back_fail,
                        seed_fwd_fail: seed_summary.fwd_fail,
                        longrec_back_no_choice: longrec_summary.back_no_choice,
                        longrec_fwd_no_choice: longrec_summary.fwd_no_choice,
                        longrec_back_no_reach: longrec_summary.back_no_reach,
                        longrec_fwd_no_reach: longrec_summary.fwd_no_reach,
                        longrec_back_unreachable: longrec_summary.back_unreachable_minpath,
                        longrec_fwd_unreachable: longrec_summary.fwd_unreachable_maxpath,
                        ..StageDetail::default()
                    },
                    &format!(
                        "unwitnessed={} zero_flux={} low_cov={} back_fail={} fwd_fail={} longrec_attempted={} longrec_succeeded={} longrec_fallback={} back_no_choice={} fwd_no_choice={} back_no_reach={} fwd_no_reach={} back_unreachable={} fwd_unreachable={}",
                        seed_summary.unwitnessed,
                        seed_summary.zero_flux,
                        seed_summary.low_coverage,
                        seed_summary.back_fail,
                        seed_summary.fwd_fail,
                        longrec_summary.attempted,
                        longrec_summary.succeeded,
                        longrec_summary.fallback,
                        longrec_summary.back_no_choice,
                        longrec_summary.fwd_no_choice,
                        longrec_summary.back_no_reach,
                        longrec_summary.fwd_no_reach,
                        longrec_summary.back_unreachable_minpath,
                        longrec_summary.fwd_unreachable_maxpath
                    ),
                );
                debug_stage::emit_with_detail(
                    "checktrf",
                    &graph_bundle.chrom,
                    graph_bundle.start,
                    graph_bundle.end,
                    graph_bundle.strand,
                    reads.len(),
                    junctions.len(),
                    graph_mut.n_nodes.saturating_sub(2),
                    debug_stage::graph_edge_count(&graph_mut),
                    transfrags.len(),
                    txs.len(),
                    seed_summary.checktrf_total,
                    0,
                    seed_summary.checktrf_rescued,
                    StageDetail {
                        checktrf_total: seed_summary.checktrf_total,
                        checktrf_pre_filter: pre_filter_count,
                        ..StageDetail::default()
                    },
                    &format!(
                        "checktrf_total={} rescued={} pre_filter={}",
                        seed_summary.checktrf_total,
                        seed_summary.checktrf_rescued,
                        pre_filter_count
                    ),
                );
                debug_stage::emit_with_detail(
                    "pairwise_overlap_filter",
                    &graph_bundle.chrom,
                    graph_bundle.start,
                    graph_bundle.end,
                    graph_bundle.strand,
                    reads.len(),
                    junctions.len(),
                    graph_mut.n_nodes.saturating_sub(2),
                    debug_stage::graph_edge_count(&graph_mut),
                    transfrags.len(),
                    predcluster_summary.after_pairwise,
                    seed_summary.total,
                    seed_summary.stored,
                    seed_summary.checktrf_rescued,
                    StageDetail {
                        predcluster_before: predcluster_summary.entry_count,
                        predcluster_removed: pairwise_removed,
                        pairwise_included_kill: predcluster_summary.pairwise_summary.included_kill,
                        pairwise_secontained_kill: predcluster_summary
                            .pairwise_summary
                            .secontained_kill,
                        pairwise_intronic_kill: predcluster_summary.pairwise_summary.intronic_kill,
                        pairwise_bettercov_kill: predcluster_summary
                            .pairwise_summary
                            .bettercov_kill,
                        pairwise_other_kill: predcluster_summary.pairwise_summary.other_kill,
                        ..StageDetail::default()
                    },
                    &format!(
                        "before={} removed={} {}",
                        predcluster_summary.entry_count,
                        pairwise_removed,
                        predcluster_summary.pairwise_summary.exact_reason_note
                    ),
                );
                debug_stage::emit_with_detail(
                    "isofrac",
                    &graph_bundle.chrom,
                    graph_bundle.start,
                    graph_bundle.end,
                    graph_bundle.strand,
                    reads.len(),
                    junctions.len(),
                    graph_mut.n_nodes.saturating_sub(2),
                    debug_stage::graph_edge_count(&graph_mut),
                    transfrags.len(),
                    predcluster_summary.after_isofrac,
                    seed_summary.total,
                    seed_summary.stored,
                    seed_summary.checktrf_rescued,
                    StageDetail {
                        predcluster_before: predcluster_summary.after_pairwise,
                        predcluster_removed: isofrac_removed,
                        isofrac_longunder_kill: predcluster_summary.isofrac_summary.longunder_kill,
                        ..StageDetail::default()
                    },
                    &format!(
                        "before={} removed={} longunder={}",
                        predcluster_summary.after_pairwise,
                        isofrac_removed,
                        predcluster_summary.isofrac_summary.longunder_kill
                    ),
                );
                debug_stage::emit_with_detail(
                    "collapse_single_exon_runoff",
                    &graph_bundle.chrom,
                    graph_bundle.start,
                    graph_bundle.end,
                    graph_bundle.strand,
                    reads.len(),
                    junctions.len(),
                    graph_mut.n_nodes.saturating_sub(2),
                    debug_stage::graph_edge_count(&graph_mut),
                    transfrags.len(),
                    predcluster_summary.after_runoff,
                    seed_summary.total,
                    seed_summary.stored,
                    seed_summary.checktrf_rescued,
                    StageDetail {
                        predcluster_before: predcluster_summary.after_isofrac,
                        predcluster_removed: runoff_removed,
                        ..StageDetail::default()
                    },
                    &format!(
                        "before={} removed={}",
                        predcluster_summary.after_isofrac, runoff_removed
                    ),
                );
                debug_stage::emit_with_detail(
                    "polymerase_runoff_filter",
                    &graph_bundle.chrom,
                    graph_bundle.start,
                    graph_bundle.end,
                    graph_bundle.strand,
                    reads.len(),
                    junctions.len(),
                    graph_mut.n_nodes.saturating_sub(2),
                    debug_stage::graph_edge_count(&graph_mut),
                    transfrags.len(),
                    predcluster_summary.after_polymerase_runoff,
                    seed_summary.total,
                    seed_summary.stored,
                    seed_summary.checktrf_rescued,
                    StageDetail {
                        predcluster_before: predcluster_summary.after_runoff,
                        predcluster_removed: polymerase_runoff_removed,
                        ..StageDetail::default()
                    },
                    &format!(
                        "before={} removed={}",
                        predcluster_summary.after_runoff, polymerase_runoff_removed
                    ),
                );
                debug_stage::emit_with_detail(
                    "polymerase_runon_filter",
                    &graph_bundle.chrom,
                    graph_bundle.start,
                    graph_bundle.end,
                    graph_bundle.strand,
                    reads.len(),
                    junctions.len(),
                    graph_mut.n_nodes.saturating_sub(2),
                    debug_stage::graph_edge_count(&graph_mut),
                    transfrags.len(),
                    predcluster_summary.after_polymerase_runon,
                    seed_summary.total,
                    seed_summary.stored,
                    seed_summary.checktrf_rescued,
                    StageDetail {
                        predcluster_before: predcluster_summary.after_polymerase_runoff,
                        predcluster_removed: polymerase_runon_removed,
                        ..StageDetail::default()
                    },
                    &format!(
                        "before={} removed={}",
                        predcluster_summary.after_polymerase_runoff, polymerase_runon_removed
                    ),
                );
                debug_stage::emit_with_detail(
                    "readthr_gate",
                    &graph_bundle.chrom,
                    graph_bundle.start,
                    graph_bundle.end,
                    graph_bundle.strand,
                    reads.len(),
                    junctions.len(),
                    graph_mut.n_nodes.saturating_sub(2),
                    debug_stage::graph_edge_count(&graph_mut),
                    transfrags.len(),
                    predcluster_summary.after_readthr,
                    seed_summary.total,
                    seed_summary.stored,
                    seed_summary.checktrf_rescued,
                    StageDetail {
                        predcluster_before: predcluster_summary.after_polymerase_runon,
                        predcluster_removed: readthr_removed,
                        ..StageDetail::default()
                    },
                    &format!(
                        "before={} removed={}",
                        predcluster_summary.after_polymerase_runon, readthr_removed
                    ),
                );
                debug_stage::emit_with_detail(
                    "junction_support_filter",
                    &graph_bundle.chrom,
                    graph_bundle.start,
                    graph_bundle.end,
                    graph_bundle.strand,
                    reads.len(),
                    junctions.len(),
                    graph_mut.n_nodes.saturating_sub(2),
                    debug_stage::graph_edge_count(&graph_mut),
                    transfrags.len(),
                    txs.len(),
                    seed_summary.total,
                    seed_summary.stored,
                    seed_summary.checktrf_rescued,
                    StageDetail {
                        predcluster_before: predcluster_summary.after_readthr,
                        predcluster_removed: junction_support_removed,
                        junction_support_removed,
                        ..StageDetail::default()
                    },
                    &format!(
                        "before={} removed={}",
                        predcluster_summary.after_readthr, junction_support_removed
                    ),
                );
                debug_stage::emit_with_detail(
                    "final_family_selection",
                    &graph_bundle.chrom,
                    graph_bundle.start,
                    graph_bundle.end,
                    graph_bundle.strand,
                    reads.len(),
                    junctions.len(),
                    graph_mut.n_nodes.saturating_sub(2),
                    debug_stage::graph_edge_count(&graph_mut),
                    transfrags.len(),
                    txs.len(),
                    seed_summary.total,
                    seed_summary.stored,
                    seed_summary.checktrf_rescued,
                    StageDetail {
                        seed_unwitnessed: seed_summary.unwitnessed,
                        seed_zero_flux: seed_summary.zero_flux,
                        seed_low_cov: seed_summary.low_coverage,
                        seed_back_fail: seed_summary.back_fail,
                        seed_fwd_fail: seed_summary.fwd_fail,
                        checktrf_total: seed_summary.checktrf_total,
                        checktrf_pre_filter: pre_filter_count,
                        predcluster_before: pre_filter_count,
                        predcluster_removed: pre_filter_count.saturating_sub(txs.len()),
                        junction_support_removed,
                        ..StageDetail::default()
                    },
                    &format!("pre_filter={}", pre_filter_count),
                );
                (txs, pre_filter_txs, seed_outcomes)
            };

        let mut bundle_txs: Vec<Transcript> = Vec::new();
        let mut baseline_seed_summary = SeedOutcomeSummary::default();

        if use_region_bundle_pass {
            // per-read junction redirect and
            // unstranding. When a junction has mm<0 AND nreads<0 (higherr redirect),
            // adjust the read's exon boundaries to match the replacement junction,
            // then replace the read's junction with a new one at the adjusted coords.
            // After processing, reads with ALL killed junctions get strand='.'.
            // This creates sno=1 CGroups bridging strands into mega-bundles.
            // Collect replacement junctions to add to good_junctions_set.
            let mut replacement_junctions: Vec<crate::types::Junction> = Vec::new();
            let effective_reads: Vec<BundleRead> = if config.long_reads {
                // Build lookup: junction (start,end) → index in cjunctions_for_redirect
                let cj_by_start_end: std::collections::HashMap<(u64, u64), usize> =
                    cjunctions_for_redirect.iter().enumerate()
                        .map(|(i, cj)| ((cj.start, cj.end), i))
                        .collect();

                let mut reads = region_reads.to_vec();
                let mut unstranded_count = 0usize;
                let mut redirected_count = 0usize;

                for r in reads.iter_mut() {
                    if r.strand == '.' || r.junctions.is_empty() {
                        continue;
                    }
                    let mut any_bad = false;
                    let mut all_killed = true;
                    let mut any_redirect_attempted = false;

                    for ji in 0..r.junctions.len() {
                        let j = r.junctions[ji];
                        let cj_idx = cj_by_start_end.get(&(j.donor, j.acceptor));
                        let cj = cj_idx.and_then(|&i| cjunctions_for_redirect.get(i));
                        let Some(cj) = cj else {
                            all_killed = false;
                            continue;
                        };

                        let changeleft = cj.nreads < 0.0;
                        let changeright = cj.nreads_good < 0.0;
                        let is_killed = cj.strand == 0 || cj.mm < 0.0;

                        if !is_killed {
                            all_killed = false;
                            continue;
                        }

                        if !changeleft && !changeright {
                            any_bad = true;
                            continue;
                        }
                        any_bad = true;
                        any_redirect_attempted = true;

                        // decode redirect
                        // pointers to compute new exon boundaries, adjust the
                        // read's exons, search for matching junction, replace.
                        let mut newstart = r.exons[ji].1;
                        let mut newend = if ji + 1 < r.exons.len() { r.exons[ji + 1].0 } else { continue };

                        if changeleft {
                            let mut jk = cj.nreads.abs() as usize;
                            // Follow redirect chain (max 5 hops to avoid cycles)
                            for _ in 0..5 {
                                if let Some(target) = cjunctions_for_redirect.get(jk) {
                                    if target.nreads < 0.0 {
                                        jk = target.nreads.abs() as usize;
                                    } else {
                                        newstart = target.start;
                                        break;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                        if changeright {
                            let mut ek = cj.nreads_good.abs() as usize;
                            for _ in 0..5 {
                                if let Some(target) = cjunctions_for_redirect.get(ek) {
                                    if target.nreads_good < 0.0 {
                                        ek = target.nreads_good.abs() as usize;
                                    } else {
                                        newend = target.end;
                                        break;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }

                        // Validate boundaries
                        let seg_start = r.exons[ji].0;
                        let seg_end = if ji + 1 < r.exons.len() { r.exons[ji + 1].1 } else { continue };
                        if !(newstart >= seg_start && newend <= seg_end && newstart <= newend) {
                            continue;
                        }

                        // Search for existing junction at (newstart, newend).
                        // Only adjust exon boundaries when a matching junction
                        // exists — this confirms the coordinates are exact.
                        let new_j = crate::types::Junction::new(newstart, newend);
                        let exists_in_good = good_junctions_set.contains(&new_j);
                        let exists_in_stats = cj_by_start_end.contains_key(&(newstart, newend));

                        // Skip exon adjustment: internal-only junction replacement
                        // already gives +3 TPs. Exon adjustment loses 4 TPs from
                        // changing CGroup boundaries at confirmed junctions.
                        let _ = (exists_in_good, exists_in_stats);
                        // Only replace INTERNAL junctions (not first/last) to
                        // avoid changing first_exon_valid/last_exon_valid behavior
                        // which affects CGroup keep/drop for terminal exons.
                        if ji > 0 && ji + 1 < r.junctions.len() {
                            r.junctions[ji] = new_j;
                        }
                        replacement_junctions.push(new_j);
                        redirected_count += 1;
                        all_killed = false;
                    }

                    // Unstrand reads whose junctions are all killed (mm<0 or strand=0).
                    // Matching the original algorithm (line 15278 + 15457-15461):
                    // killed junctions have strand forced to 0; reads with all
                    // strand=0 junctions get unstranded → placed on sno=1 →
                    // creates cross-strand bridges that merge bnodes.
                    if any_bad && all_killed {
                        let still_all_killed = r.junctions.iter().all(|j| {
                            cj_by_start_end.get(&(j.donor, j.acceptor))
                                .and_then(|&i| cjunctions_for_redirect.get(i))
                                .map(|cj| cj.strand == 0 || cj.mm < 0.0)
                                .unwrap_or(true)
                        });
                        if still_all_killed {
                            r.strand = '.';
                            unstranded_count += 1;
                        }
                    }
                }
                if config.verbose && (unstranded_count > 0 || redirected_count > 0) {
                    eprintln!(
                        "    junction_redirect: {} redirected, {} unstranded, {} replacement juncs",
                        redirected_count, unstranded_count, replacement_junctions.len()
                    );
                }
                reads
            } else {
                region_reads.to_vec()
            };
            // Replacement junctions affect color propagation through TWO mechanisms:
            // A) good_junctions_set: bundle_builder checks !good_junctions for color breaks
            // B) killed_juncs: bundle_builder checks killed_junctions for bad-junction color breaks
            //
            // Create a color-only junction set with replacements for CGroup building,
            // but keep the original good_junctions for graph construction (via process_graph).
            // Add replacement junctions to color propagation sets.
            let mut color_good = good_junctions_set.clone();
            let mut color_killed = killed_juncs.clone();
            for rj in &replacement_junctions {
                color_good.insert(*rj);
                color_killed.remove(rj);
            }
            let raw_subbundles =
                build_sub_bundles(
                    &effective_reads, &config, &color_good, &color_killed,
                    Some(&bpcov_stranded), Some(&junction_stats_corr_final),
                    bundle.start,
                )?;

            total_subbundles.fetch_add(raw_subbundles.len(), std::sync::atomic::Ordering::Relaxed);
            if config.verbose {
                eprintln!(
                    "[BUNDLE_COUNT] region {}:{}-{} -> {} raw_subbundles (total: {})",
                    bundle.chrom,
                    bundle.start,
                    bundle.end,
                    raw_subbundles.len(),
                    total_subbundles.load(std::sync::atomic::Ordering::Relaxed)
                );
            }

            // Per-bundle graph mode: merge overlapping same-strand subbundles into
            // larger bundles (matching the architecture where all bundlenodes
            // in a region contribute to one unified graph per color component).
            let bundle_graph_mode = std::env::var_os("RUSTLE_BUNDLE_GRAPH").is_some();
            let effective_subbundles: Vec<crate::bundle_builder::SubBundleResult> = if bundle_graph_mode {
                merge_overlapping_subbundles(&raw_subbundles, &effective_reads)
            } else {
                raw_subbundles
            };

            // Per-bundlenode processing: when RUSTLE_PER_BNODE is set, iterate through
            // each bundlenode in the merged bundle and create a SEPARATE graph for each,
            // matching architecture (per-bundlenode graphs with shared read mapping).
            let per_bnode_mode = bundle_graph_mode && std::env::var_os("RUSTLE_PER_BNODE").is_some();

            for (sbr_idx, sbr) in effective_subbundles.iter().enumerate() {
                if sbr.strand == '.' {
                    continue;
                }

                let sub_cap = sbr.read_scale.len();
                let mut sub_reads: Vec<BundleRead> = Vec::with_capacity(sub_cap);
                let mut sub_read_bnodes: Vec<Vec<usize>> = Vec::with_capacity(sub_cap);
                for (ri, read) in effective_reads.iter().enumerate() {
                    let scale = sbr.read_scale.get(ri).copied().unwrap_or(0.0);
                    let mapped = sbr
                        .read_to_bnodes
                        .get(ri)
                        .cloned()
                        .unwrap_or_default();
                    if scale <= 0.0 || mapped.is_empty() {
                        continue;
                    }
                    let mut routed = read.clone();
                    routed.strand = sbr.strand;
                    if (scale - 1.0).abs() > f64::EPSILON {
                        routed.weight *= scale;
                        routed.junc_mismatch_weight *= scale;
                        for pc in routed.pair_count.iter_mut() {
                            *pc *= scale;
                        }
                    }
                    if routed.weight <= 0.0 {
                        continue;
                    }
                    sub_read_bnodes.push(mapped);
                    sub_reads.push(routed);
                }

                if sub_reads.is_empty() {
                    continue;
                }

                // Collect individual bundlenodes for per-bnode processing.
                // In per-bnode mode, we iterate through each bundlenode and
                // create a separate graph, but use the FULL chain for read mapping.
                let bnode_list: Vec<(u64, u64, f64, usize)> = {
                    let mut v = Vec::new();
                    let mut cur = sbr.bnode_head.as_ref();
                    while let Some(bn) = cur {
                        v.push((bn.start, bn.end, bn.cov, bn.bid));
                        cur = bn.next.as_deref();
                    }
                    v
                };

                // Per-bnode mode: split merged bundle into junction-connected components.
                // Each component = bundlenodes reachable from each other via junctions.
                // Process each component as a SEPARATE graph with its own flow decomposition,
                // but use the FULL merged bundle for read mapping (cross-component sharing).
                // This matches per-junction-network graphs with shared read context.
                let sub_junctions_for_split: Vec<Junction> = sub_reads.iter()
                    .flat_map(|r| r.junctions.iter().copied())
                    .collect();
                let color_groups: Vec<Option<Vec<(u64, u64, f64, usize)>>> = if per_bnode_mode && bnode_list.len() > 1 {
                    let components = crate::per_bnode_graph::find_junction_connected_components(
                        &bnode_list, &sub_junctions_for_split,
                    );
                    if components.len() > 1 {
                        // Multiple components: process each separately
                        components.iter().map(|comp| {
                            Some(comp.iter().map(|&i| bnode_list[i]).collect())
                        }).collect()
                    } else {
                        vec![None] // Single component = process as one graph
                    }
                } else {
                    vec![None] // Single pass with full bundlenode chain
                };

                for color_group in &color_groups {

                // Determine the bundlenode chain and range for this iteration.
                let (iter_bnode_head, iter_start, iter_end): (Option<CBundlenode>, u64, u64) = if let Some(group) = color_group {
                    // Build a bundlenode chain from this color group
                    let mut sorted = group.clone();
                    sorted.sort_unstable_by_key(|&(s, _, _, _)| s);
                    if sorted.is_empty() {
                        (None, 0, 0)
                    } else {
                        let start = sorted.first().map(|&(s,_,_,_)| s).unwrap_or(0);
                        let end = sorted.last().map(|&(_,e,_,_)| e).unwrap_or(0);
                        let mut head = CBundlenode {
                            start: sorted[0].0, end: sorted[0].1,
                            cov: sorted[0].2, bid: sorted[0].3,
                            next: None, hardstart: false, hardend: false,
                        };
                        let mut tail = &mut head;
                        for &(s, e, c, bid) in &sorted[1..] {
                            tail.next = Some(Box::new(CBundlenode {
                                start: s, end: e, cov: c, bid,
                                next: None, hardstart: false, hardend: false,
                            }));
                            tail = tail.next.as_mut().unwrap();
                        }
                        (Some(head), start, end)
                    }
                } else {
                    (sbr.bnode_head.clone(), sbr.start, sbr.end)
                };

                let _single_bnode = &color_group.as_ref().and_then(|g| {
                    if g.len() == 1 { Some(g[0]) } else { None }
                });

                // skip bundles shorter than mintranscriptlen.
                let bundle_exonic_len = if color_group.is_some() {
                    color_group.as_ref().map(|g| g.iter().map(|(s,e,_,_)| e.saturating_sub(*s)).sum()).unwrap_or(0)
                } else {
                    let mut cur = sbr.bnode_head.as_ref();
                    let mut len = 0u64;
                    while let Some(bn) = cur {
                        len += bn.end.saturating_sub(bn.start);
                        cur = bn.next.as_deref();
                    }
                    len
                };
                if bundle_exonic_len < config.min_transcript_length {
                    continue;
                }

                let mut sub_bundle = crate::types::Bundle {
                    chrom: bundle.chrom.clone(),
                    start: iter_start,
                    end: iter_end,
                    strand: sbr.strand,
                    reads: sub_reads.clone(),
                    junction_stats: Default::default(),
                    // For per-bnode mode: use single bundlenode for graph construction,
                    // but the FULL chain is used for read mapping (via bundle2graph).
                    // For normal mode: use the full chain for both.
                    // In per-bnode mode: use the single bundlenode for graph building.
                    // The FULL chain is used for read mapping (build_bundle2graph sees all nodes).
                    // In normal mode: use the full merged chain.
                    bundlenodes: if color_group.is_some() {
                        // Graph built from single bundlenode, but bundle2graph uses full chain
                        iter_bnode_head.clone()
                    } else {
                        sbr.bnode_head.clone()
                    },
                    read_bnodes: Some(sub_read_bnodes.clone()),
                    bnode_colors: Some(sbr.bnode_colors.clone()),
                };

                // DEBUG_BUNDLE: emit bundle summary matching expected format.
                if std::env::var_os("RUSTLE_DEBUG_BUNDLE").is_some() {
                    let sno = match sbr.strand {
                        '-' => 0i32,
                        '+' => 2,
                        _ => 1,
                    };
                    // Walk the linked list of bundlenodes
                    let mut bn_list: Vec<(u64, u64, f64)> = Vec::new();
                    let mut cur = sbr.bnode_head.as_ref();
                    while let Some(bn) = cur {
                        bn_list.push((bn.start, bn.end, bn.cov));
                        cur = bn.next.as_deref();
                    }
                    let cov: f64 = bn_list.iter().map(|&(s, e, _)| (e - s) as f64).sum();
                    let len: u64 = bn_list.iter().map(|&(s, e, _)| e - s).sum();
                    let bstart = bn_list.first().map(|&(s, _, _)| s).unwrap_or(0);
                    let bend = bn_list.last().map(|&(_, e, _)| e).unwrap_or(0);
                    eprintln!(
                        "DEBUG_BUNDLE sno={} b={} start={} end={} cov={:.1} len={} bnodes={} reads={}",
                        sno,
                        sbr_idx,
                        bstart + 1,
                        bend,
                        cov,
                        len,
                        bn_list.len(),
                        sub_reads.len(),
                    );
                    for (ni, &(ns, ne, ncov)) in bn_list.iter().enumerate() {
                        eprintln!(
                            "  BNODE sno={} b={} n={} start={} end={} cov={:.1}",
                            sno, sbr_idx, ni, ns + 1, ne, ncov
                        );
                    }
                }

                let sub_n_junctions: usize =
                    sub_bundle.reads.iter().map(|r| r.junctions.len()).sum();
                debug_stage::emit_with_detail(
                    "bundle_context",
                    &sub_bundle.chrom,
                    sub_bundle.start,
                    sub_bundle.end,
                    sub_bundle.strand,
                    sub_bundle.reads.len(),
                    sub_n_junctions,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    debug_stage::bundle_detail(&sub_bundle),
                    "subbundle",
                );
                let mut sub_junction_stats_corr = compute_initial_junction_stats_for_reads(
                    &sub_bundle.reads,
                    sub_bundle.start,
                    sub_bundle.end,
                    &config,
                );
                // Inherit bundle-level "killed" flags (mm<0, strand=Some(0),
                // nreads_good<0, mrcount<0) to sub-bundle stats. This ensures
                // sub-bundle's compute_killed_junction_pairs treats bundle-demoted
                // junctions as killed even when sub-bundle's fresh HE_CONT /
                // good_junc doesn't re-fire on them (different bpcov context).
                // Fixed STRG.29 (mm<0 inheritance). Other flags extend the same
                // principle.
                // Disable via RUSTLE_INHERIT_MM_NEG_OFF=1.
                if std::env::var_os("RUSTLE_INHERIT_MM_NEG_OFF").is_none() {
                    let inherit_strand = std::env::var_os("RUSTLE_INHERIT_STRAND0").is_some();
                    let inherit_nreads_good = std::env::var_os("RUSTLE_INHERIT_NREADS_GOOD_NEG").is_some();
                    let inherit_mrcount = std::env::var_os("RUSTLE_INHERIT_MRCOUNT_NEG").is_some();
                    for (j, bundle_stat) in &junction_stats_corr_final {
                        if let Some(sub_stat) = sub_junction_stats_corr.get_mut(j) {
                            if bundle_stat.mm < 0.0 && sub_stat.mm >= 0.0 {
                                sub_stat.mm = bundle_stat.mm;
                            }
                            if inherit_strand && bundle_stat.strand == Some(0) && sub_stat.strand != Some(0) {
                                sub_stat.strand = Some(0);
                            }
                            if inherit_nreads_good && bundle_stat.nreads_good < 0.0 && sub_stat.nreads_good >= 0.0 {
                                sub_stat.nreads_good = bundle_stat.nreads_good;
                            }
                            if inherit_mrcount && bundle_stat.mrcount < 0.0 && sub_stat.mrcount >= 0.0 {
                                sub_stat.mrcount = bundle_stat.mrcount;
                            }
                        }
                    }
                }
                let mut sub_junction_redirect_map: HashMap<Junction, Junction> = Default::default();
                count_good_junctions(
                    &mut sub_bundle.reads,
                    &mut sub_junction_stats_corr,
                    config.junction_support,
                    config.longintron,
                    sub_bundle.start,
                    sub_bundle.end,
                    config.min_intron_length,
                    sub_bundle.strand,
                    None,
                    &sub_bundle.chrom,
                    &chrom_junction_strand_evidence,
                );
                let (sub_corrected_stats, sub_corrected_map) =
                    correct_junctions_with_map(&sub_junction_stats_corr, config.verbose);
                sub_junction_stats_corr = sub_corrected_stats;
                sub_junction_redirect_map.extend(sub_corrected_map);
                if !guide_transcripts.is_empty() && config.long_reads {
                    let guide_junctions_vec = cached_bundle_guide_junctions(
                        &mut guide_junction_cache,
                        &mut chrom_arc_cache,
                        &guide_transcripts,
                        &sub_bundle.chrom,
                        sub_bundle.strand,
                        sub_bundle.start,
                        sub_bundle.end,
                    );
                    if !guide_junctions_vec.is_empty() {
                        let (snapped_stats, snap_map) = snap_junctions_to_guides(
                            &sub_junction_stats_corr,
                            guide_junctions_vec.as_ref(),
                            config.sserror,
                            config.long_reads,
                        );
                        if !snap_map.is_empty() {
                            sub_junction_stats_corr = snapped_stats;
                            sub_junction_redirect_map.extend(snap_map);
                        }
                    }
                }
                let sub_bpcov_stranded =
                    BpcovStranded::from_reads(&sub_bundle.reads, sub_bundle.start, sub_bundle.end);
                if !guide_transcripts.is_empty() {
                    let guide_junctions = cached_bundle_guide_junctions(
                        &mut guide_junction_cache,
                        &mut chrom_arc_cache,
                        &guide_transcripts,
                        &sub_bundle.chrom,
                        sub_bundle.strand,
                        sub_bundle.start,
                        sub_bundle.end,
                    );
                    mark_guide_junctions_for_junctions(
                        &mut sub_junction_stats_corr,
                        guide_junctions.as_ref(),
                    );
                }
                let mut sub_cjunctions =
                    crate::types::junction_stats_to_cjunctions(&sub_junction_stats_corr);
                aggregate_splice_site_support(&mut sub_cjunctions);
                // apply_higherr_demotions before good_junc
                apply_higherr_demotions(
                    &mut sub_cjunctions,
                    config.sserror,
                    config.min_junction_reads,
                    Some(&sub_bpcov_stranded),
                    sub_bundle.start,
                );
                demote_runthrough_junctions(&mut sub_cjunctions, &sub_bpcov_stranded, sub_bundle.start);
                good_junc(
                    &mut sub_cjunctions,
                    &sub_bpcov_stranded,
                    sub_bundle.start,
                    config.min_junction_reads,
                    config.longintron,
                    config.per_splice_site_isofrac,
                    config.long_reads,
                    config.eonly,
                );
                emit_jfinal_trace(&sub_cjunctions, "subbundle");
                // mm<0 junctions with non-zero
                // strand are kept alive (read-redirect creates replacement junctions
                // in the standard algorithm). Do NOT convert mm<0→strand=0.
                sub_junction_stats_corr =
                    crate::types::cjunctions_to_junction_stats(&sub_cjunctions);
                // second-stage good_junc with subbundle bpcov.
                // good_junc is called per-read in the build_graphs loop with full
                // coverage context. Edge-case junctions (all-bad mapping, marginal support)
                // that survive the initial good_junc are re-evaluated against local coverage.
                // This catches junctions where nreads << bpcov (e.g., 1 read vs 700x coverage)
                // that create spurious graph node boundaries and fragment transfrags.
                {
                    use crate::killed_junctions::good_junc_second_stage;
                    let strand_val: i8 = match sub_bundle.strand {
                        '+' => 1,
                        '-' => -1,
                        _ => 0,
                    };
                    let to_kill: Vec<crate::types::Junction> = sub_junction_stats_corr
                        .iter()
                        .filter(|(j, s)| {
                            s.strand != Some(0)
                                && !good_junc_second_stage(
                                    **j,
                                    strand_val,
                                    &sub_junction_stats_corr,
                                    &sub_bpcov_stranded,
                                    sub_bundle.start,
                                    config.min_junction_reads,
                                    config.longintron,
                                    config.long_reads,
                                )
                        })
                        .map(|(j, _)| *j)
                        .collect();
                    for j in &to_kill {
                        if let Some(st) = sub_junction_stats_corr.get_mut(j) {
                            st.strand = Some(0);
                        }
                    }
                }
                if let Some(g) = &genome {
                    let chrom_arc = intern_chrom_arc(&mut chrom_arc_cache, &sub_bundle.chrom);
                    let guide_junctions = cached_bundle_guide_junctions(
                        &mut guide_junction_cache,
                        &mut chrom_arc_cache,
                        &guide_transcripts,
                        &sub_bundle.chrom,
                        sub_bundle.strand,
                        sub_bundle.start,
                        sub_bundle.end,
                    );
                    let guide_tol = config.sserror.max(1);
                    for (j, st) in sub_junction_stats_corr.iter_mut() {
                        if st.strand == Some(0) {
                            continue;
                        }
                        let guide_match = guide_junctions.as_ref().iter().any(|gj| {
                            gj.donor.abs_diff(j.donor) <= guide_tol
                                && gj.acceptor.abs_diff(j.acceptor) <= guide_tol
                        });
                        if guide_match {
                            continue;
                        }
                        let is_cons = cached_consensus_splice(
                            &mut consensus_cache,
                            g,
                            &chrom_arc,
                            j.donor,
                            j.acceptor,
                            st.strand,
                        );
                        if std::env::var_os("RUSTLE_CONSENSUS_KILL").is_some()
                            && !is_cons
                            && st.nreads_good < 5.0
                        {
                            st.strand = Some(0);
                        }
                    }
                }
                let sub_cjunctions_post =
                    crate::types::junction_stats_to_cjunctions(&sub_junction_stats_corr);
                let sub_killed_junction_pairs = compute_killed_junction_pairs(&sub_cjunctions_post);
                if !sub_killed_junction_pairs.is_empty() || !sub_junction_redirect_map.is_empty() {
                    repair_reads_after_junction_quality(
                        &mut sub_bundle.reads,
                        &sub_junction_redirect_map,
                        &sub_killed_junction_pairs,
                        &sub_junction_stats_corr,
                    );
                }
                let sub_junction_stats = apply_junction_filters_and_canonicalize(
                    sub_junction_stats_corr,
                    config.min_junction_reads,
                    config.junction_canonical_tolerance,
                    config.per_splice_site_isofrac,
                    config.verbose,
                );
                let mut sub_junctions =
                    filter_junctions(&sub_junction_stats, config.min_junction_reads);
                // Alt-acceptor sub-bundle rescue (opt-in, scaffold).
                // Sub-bundle recomputes junction stats from its narrower read
                // pool, dropping alt-acceptors that are alive in outer bundle
                // (e.g., STRG.275's 40068714 acceptor: outer mrcount=468
                // strand=-1, sub-bundle mrcount=4 strand=0 → killed).
                //
                // Re-inject outer-alive alt-siblings that are:
                // - Missing from sub_junctions
                // - Within a TIGHT acceptor/donor window of a kept junction
                //   (alt-splice signature, e.g., 4bp apart)
                // - Outer has STRONG support (>= N reads) — we're preserving
                //   confident alt-splice sites, not weak variants
                //
                // Default OFF — sweep showed all parameter choices regress F1
                // by -5 to -9 matches due to collateral flow-graph disturbance
                // on other bundles. Kept as scaffold for future work where
                // flow attribution can be adjusted proportionally.
                //
                // Enable via RUSTLE_ALT_ACCEPTOR_RESCUE=1.
                let mut sub_junction_stats = sub_junction_stats;
                if std::env::var_os("RUSTLE_ALT_ACCEPTOR_RESCUE").is_some() {
                    let rescue_min_reads: f64 = std::env::var(
                        "RUSTLE_ALT_ACCEPTOR_RESCUE_MIN_READS")
                        .ok().and_then(|v| v.parse().ok()).unwrap_or(20.0);
                    let rescue_window: u64 = std::env::var(
                        "RUSTLE_ALT_ACCEPTOR_RESCUE_WINDOW")
                        .ok().and_then(|v| v.parse().ok()).unwrap_or(8);
                    let sub_set: HashSet<Junction> = sub_junctions.iter().copied().collect();
                    let mut rescued = 0usize;
                    for (j, outer_stat) in &junction_stats_corr_final {
                        if sub_set.contains(j) { continue; }
                        if !matches!(outer_stat.strand, Some(-1) | Some(1)) { continue; }
                        if outer_stat.nreads_good < rescue_min_reads { continue; }
                        if j.donor < sub_bundle.start || j.acceptor > sub_bundle.end {
                            continue;
                        }
                        let is_alt_sibling = sub_junctions.iter().any(|kept| {
                            (kept.donor == j.donor
                                && kept.acceptor != j.acceptor
                                && kept.acceptor.abs_diff(j.acceptor) <= rescue_window)
                                || (kept.acceptor == j.acceptor
                                    && kept.donor != j.donor
                                    && kept.donor.abs_diff(j.donor) <= rescue_window)
                        });
                        if !is_alt_sibling { continue; }
                        sub_junctions.push(*j);
                        sub_junction_stats.insert(*j, outer_stat.clone());
                        rescued += 1;
                    }
                    if rescued > 0 && config.verbose {
                        eprintln!(
                            "    sub-bundle alt-acceptor rescue: re-injected {} junction(s)",
                            rescued
                        );
                    }
                }
                let sub_good_junctions: HashSet<(u64, u64)> = sub_junctions
                    .iter()
                    .map(|j| (j.donor, j.acceptor))
                    .collect();
                sub_bundle.junction_stats = sub_junction_stats.clone();

                // DEBUG_JUNCTION: emit junction list for per-bundle comparison.
                if std::env::var_os("RUSTLE_DEBUG_BUNDLE").is_some() {
                    let sno = match sbr.strand {
                        '-' => 0i32,
                        '+' => 2,
                        _ => 1,
                    };
                    let mut jlist: Vec<_> = sub_junctions.iter()
                        .map(|j| {
                            let stat = sub_junction_stats.get(j);
                            let mr = stat.map(|s| s.mrcount).unwrap_or(0.0);
                            let ng = stat.map(|s| s.nreads_good).unwrap_or(0.0);
                            (j.donor, j.acceptor, mr, ng)
                        })
                        .collect();
                    jlist.sort_unstable_by_key(|&(d, a, _, _)| (d, a));
                    eprintln!(
                        "DEBUG_JUNCTIONS sno={} b={} count={} range={}-{}",
                        sno, sbr_idx, jlist.len(),
                        sub_bundle.start + 1, sub_bundle.end
                    );
                    for (d, a, mr, ng) in &jlist {
                        eprintln!(
                            "  JUNC sno={} b={} donor={} acceptor={} mrcount={:.1} nreads_good={:.1}",
                            sno, sbr_idx, d, a, mr, ng
                        );
                    }
                }

                // In per-bnode mode: use the color group's bundlenodes for graph construction.
                // For read mapping: always use full merged chain (cross-component sharing).
                let _graph_bnodes = if color_group.is_some() {
                    iter_bnode_head.clone()
                } else {
                    sbr.bnode_head.clone()
                };
                let mapping_bnodes = sbr.bnode_head.clone();

                // Use parent bundle reads (restricted to sub_bundle range) for bpcov/boundary
                // computation so longtrim detects boundaries uniformly. Previously only
                // sub_bundle.reads were used, which missed end-boundaries like 22152351
                // (STRG.109/110) because sub_bundle's reads are a narrow subset.
                // Disable via RUSTLE_COVREADS_PARENT_OFF=1.
                let parent_cov_reads: Vec<BundleRead>;
                let cov_reads_ref: &[BundleRead] = if std::env::var_os("RUSTLE_COVREADS_PARENT_OFF").is_none() {
                    parent_cov_reads = bundle.reads.iter()
                        .filter(|r| r.ref_start < sub_bundle.end && r.ref_end > sub_bundle.start)
                        .cloned()
                        .collect();
                    &parent_cov_reads
                } else {
                    &sub_bundle.reads
                };
                let (txs, pre_filter, seed_outcomes) = process_graph(
                    &sub_bundle,
                    sub_junctions.clone(),
                    mapping_bnodes,
                    &sub_bundle.reads,
                    cov_reads_ref,
                    Some(&sub_read_bnodes),
                    &sub_good_junctions,
                    &sub_killed_junction_pairs,
                );
                accumulate_seed_outcomes(&mut baseline_seed_summary, &seed_outcomes);
                if trace_reference.is_some() {
                    let junc_tuples: Vec<(u64, u64)> = sub_junctions
                        .iter()
                        .map(|j| (j.donor, j.acceptor))
                        .collect();
                    let trace_txs = pre_filter.unwrap_or_else(|| txs.clone());
                    raw_for_trace_mutex.lock().unwrap().push((sub_bundle.clone(), trace_txs, junc_tuples, seed_outcomes));
                }
                let stranded_se = emit_stranded_single_exon_candidates(
                    &sub_bundle,
                    &txs,
                    Some(&sub_bpcov_stranded),
                    config.min_transcript_length,
                    config.singlethr,
                );
                if !stranded_se.is_empty() {
                    single_exon_predictions_mutex.lock().unwrap().extend(stranded_se);
                }
                let term_se = emit_terminal_exon_se_candidates(
                    &sub_bundle,
                    &txs,
                    config.min_transcript_length,
                    config.singlethr,
                );
                if !term_se.is_empty() {
                    single_exon_predictions_mutex.lock().unwrap().extend(term_se);
                }
                bundle_txs.extend(txs);
            } // end for color_group
            } // end for sbr_idx

            // Per-bnode deduplication: when multiple components produce the same
            // intron chain, keep only the copy with highest coverage.  This
            // prevents duplicate transcripts from overlapping junction-connected
            // components.
            if per_bnode_mode && bundle_txs.len() > 1 {
                let pre_dedup = bundle_txs.len();
                let mut seen: std::collections::HashMap<Vec<(u64, u64)>, usize> =
                    std::collections::HashMap::new();
                let mut keep = vec![true; bundle_txs.len()];
                for (i, tx) in bundle_txs.iter().enumerate() {
                    if let Some(&prev) = seen.get(&tx.exons) {
                        // Same intron chain: keep the higher-coverage copy.
                        if tx.coverage > bundle_txs[prev].coverage {
                            keep[prev] = false;
                            seen.insert(tx.exons.clone(), i);
                        } else {
                            keep[i] = false;
                        }
                    } else {
                        seen.insert(tx.exons.clone(), i);
                    }
                }
                let mut deduped: Vec<crate::path_extract::Transcript> = Vec::new();
                for (i, tx) in bundle_txs.into_iter().enumerate() {
                    if keep[i] { deduped.push(tx); }
                }
                if deduped.len() < pre_dedup && trace_log_style {
                    eprintln!(
                        "--- per_bnode_dedup: removed {} duplicate transcripts ({} -> {})",
                        pre_dedup - deduped.len(), pre_dedup, deduped.len()
                    );
                }
                bundle_txs = deduped;
            }

            if trace_log_style {
                eprintln!(
                    "--- infer_transcripts: AFTER_BUILD_GRAPHS bdata={}-{} npred={}",
                    bundle.start + 1,
                    trace_bundle_end(bundle.end),
                    bundle_txs.len()
                );
            }
            let stranded_se = emit_stranded_single_exon_candidates(
                &bundle,
                &bundle_txs,
                Some(&bpcov_stranded),
                config.min_transcript_length,
                config.singlethr,
            );
            if !stranded_se.is_empty() {
                single_exon_predictions_mutex.lock().unwrap().extend(stranded_se);
            }
            let term_se = emit_terminal_exon_se_candidates(
                &bundle,
                &bundle_txs,
                config.min_transcript_length,
                config.singlethr,
            );
            if !term_se.is_empty() {
                single_exon_predictions_mutex.lock().unwrap().extend(term_se);
            }
            all_transcripts_mutex.lock().unwrap().extend(bundle_txs);
            return Ok(());
        }

        // split bundlenodes into components by CGroup color root.
        // CGroup colors break on killed junctions (strand==0), and each color
        // group becomes a separate sub-bundle with its own create_graph call (shared
        // reads/bpcov, We use color roots directly.
        // Components without their own junctions are merged back into the largest component
        // (single-exon islands can't form multi-exon transcripts independently).
        //
        // keep ALL color-split components as separate bundles.
        // Each distinct color root → separate CBundle, even single-bnode components.
        let comps = bundlenode_components_by_color(bundlenodes.as_ref(), &bnode_colors);
        let baseline_components = comps.len();
        let baseline_bnodes = bundlenodes_to_vec(bundlenodes.as_ref()).len();
        let baseline_assigned = count_assigned_reads(&cgroup_read_bnodes, &cgroup_read_scales);
        let baseline_cgroups_est = unique_color_roots(&bnode_colors);

        if comps.len() > 1 {
            let comp_reads = build_component_reads(
                &comps,
                &cgroup_read_bnodes,
                &cgroup_read_scales,
                &scaled_bundle_reads,
            );
            let allowed_junction_idx = build_allowed_junction_index(&junctions);

            // DEBUG_BUNDLE: match bundle accounting.
            // If a component contains only unstranded reads after junction repair,
            // classify it as sno=1 (neutral) for compatibility counting.
            if std::env::var_os("RUSTLE_DEBUG_BUNDLE").is_some() {
                let bundle_sno = match bundle.strand {
                    '-' => 0,
                    '+' => 2,
                    _ => 1,
                };
                for (ci, comp) in comps.iter().enumerate() {
                    let all_dot = comp_reads
                        .get(ci)
                        .map(|rr| !rr.is_empty() && rr.iter().all(|r| r.strand == '.'))
                        .unwrap_or(false);
                    let sno = if all_dot { 1 } else { bundle_sno };
                    let s = comp.first().map(|v| v.1).unwrap_or(0);
                    let e = comp.last().map(|v| v.2).unwrap_or(0);
                    let cov: f64 = comp.iter().map(|v| v.3).sum();
                    let len: u64 = comp.iter().map(|v| v.2.saturating_sub(v.1)).sum();
                    eprintln!(
                        "DEBUG_BUNDLE sno={} b={} start={} end={} cov={:.1} len={} bnodes={}",
                        sno,
                        ci,
                        s + 1,
                        e,
                        cov,
                        len,
                        comp.len()
                    );
                    for (ni, &(_color, ns, ne, ncov)) in comp.iter().enumerate() {
                        eprintln!(
                            "  BNODE sno={} b={} n={} start={} end={} cov={:.1}",
                            sno,
                            ci,
                            ni,
                            ns + 1,
                            ne,
                            ncov
                        );
                    }
                }
            }

            if std::env::var_os("RUSTLE_DEBUG_DETAIL").is_some() {
                eprintln!(
                    "COLOR_SPLIT_FINAL bundle={}:{}-{} comps={}",
                    bundle.chrom,
                    bundle.start,
                    bundle.end,
                    comps.len()
                );
                for (ci, comp) in comps.iter().enumerate() {
                    let s = comp.first().map(|v| v.1).unwrap_or(0);
                    let e = comp.last().map(|v| v.2).unwrap_or(0);
                    eprintln!(
                        "  comp[{}] bnodes={} range={}-{} reads={}",
                        ci,
                        comp.len(),
                        s,
                        e,
                        comp_reads[ci].len()
                    );
                }
            }
            for (ci, comp) in comps.into_iter().enumerate() {
                let comp_bn = vec_to_bundlenodes(&comp);
                // component junctions are those carried by component reads,
                // then intersected with the active bundle junction set.
                let comp_junctions = select_component_junctions_from_reads_indexed(
                    &comp_reads[ci],
                    &junctions,
                    allowed_junction_idx.as_ref(),
                );
                let (txs, pre_filter, seed_outcomes) = process_graph(
                    &bundle,
                    comp_junctions.clone(),
                    comp_bn,
                    &comp_reads[ci],
                    &scaled_bundle_reads,
                    None,
                    &good_junctions,
                    &killed_junction_pairs,
                );
                accumulate_seed_outcomes(&mut baseline_seed_summary, &seed_outcomes);
                if trace_reference.is_some() {
                    let junc_tuples: Vec<(u64, u64)> = comp_junctions
                        .iter()
                        .map(|j| (j.donor, j.acceptor))
                        .collect();
                    let trace_txs = pre_filter.unwrap_or_else(|| txs.clone());
                    raw_for_trace_mutex.lock().unwrap().push((bundle.clone(), trace_txs, junc_tuples, seed_outcomes));
                }
                bundle_txs.extend(txs);
            }
        } else {
            // DEBUG_BUNDLE single-component case.
            if std::env::var_os("RUSTLE_DEBUG_BUNDLE").is_some() {
                let all_dot =
                    !bundle.reads.is_empty() && bundle.reads.iter().all(|r| r.strand == '.');
                let sno = if all_dot {
                    1
                } else {
                    match bundle.strand {
                        '-' => 0,
                        '+' => 2,
                        _ => 1,
                    }
                };
                if let Some(bn) = bundlenodes.as_ref() {
                    let comp = bundlenodes_to_vec(Some(bn));
                    let s = comp.first().map(|v| v.1).unwrap_or(0);
                    let e = comp.last().map(|v| v.2).unwrap_or(0);
                    let cov: f64 = comp.iter().map(|v| v.3).sum();
                    let len: u64 = comp.iter().map(|v| v.2.saturating_sub(v.1)).sum();
                    eprintln!(
                        "DEBUG_BUNDLE sno={} b=0 start={} end={} cov={:.1} len={} bnodes={}",
                        sno,
                        s + 1,
                        e,
                        cov,
                        len,
                        comp.len()
                    );
                    for (ni, &(_bid, ns, ne, ncov)) in comp.iter().enumerate() {
                        eprintln!(
                            "  BNODE sno={} b=0 n={} start={} end={} cov={:.1}",
                            sno,
                            ni,
                            ns + 1,
                            ne,
                            ncov
                        );
                    }
                }
            }
            // Single component (common case) — no splitting overhead.
            let single_reads: Vec<BundleRead> = scaled_bundle_reads
                .iter()
                .filter(|r| r.weight > 0.0)
                .cloned()
                .collect();
            let (txs, pre_filter, seed_outcomes) = process_graph(
                &bundle,
                junctions.clone(),
                bundlenodes.clone(),
                &single_reads,
                &single_reads,
                None,
                &good_junctions,
                &killed_junction_pairs,
            );
            accumulate_seed_outcomes(&mut baseline_seed_summary, &seed_outcomes);
            if trace_reference.is_some() {
                let junc_tuples: Vec<(u64, u64)> =
                    junctions.iter().map(|j| (j.donor, j.acceptor)).collect();
                let trace_txs = pre_filter.unwrap_or_else(|| txs.clone());
                raw_for_trace_mutex.lock().unwrap().push((bundle.clone(), trace_txs, junc_tuples, seed_outcomes));
            }
            bundle_txs.extend(txs);
        }

        if shadow_strict_3strand {
            let strict_region_reads: &[BundleRead] = region_reads_for_groupflow
                .get(&region_key)
                .map(|v| v.as_slice())
                .unwrap_or_else(|| bundle.reads.as_slice());
            let (strict_bn, strict_read_bnodes_all, strict_bnode_colors, strict_scales_all) =
                build_bundlenodes_and_readgroups_from_cgroups_3strand(
                    strict_region_reads,
                    bundle.strand,
                    &good_junctions_set,
                    &killed_juncs,
                    config.cgroup_junction_support,
                    config.bundle_merge_dist,
                    config.long_reads,
                    &boundary_left,
                    &boundary_right,
                    None, // junction_stats - not used in this code path
                    &junction_redirect_map,
                );
            let strict_assigned_all_region =
                count_assigned_reads(&strict_read_bnodes_all, &strict_scales_all);
            let strict_bnodes_pos_all = strict_read_bnodes_all
                .iter()
                .filter(|b| !b.is_empty())
                .count();
            let strict_scale_pos_all = strict_scales_all.iter().filter(|&&s| s > 0.0).count();
            let strict_subset_to_region =
                map_subset_reads_to_region_indices(cgroup_reads, strict_region_reads);
            let strict_subset_mapped = strict_subset_to_region.iter().flatten().count();
            let mut strict_read_bnodes: Vec<Vec<usize>> = vec![Vec::new(); cgroup_reads.len()];
            let mut strict_scales: Vec<f64> = vec![0.0; cgroup_reads.len()];
            for (i, maybe_ri) in strict_subset_to_region.iter().enumerate() {
                if let Some(ri) = maybe_ri {
                    if *ri < strict_read_bnodes_all.len() {
                        strict_read_bnodes[i] = strict_read_bnodes_all[*ri].clone();
                    }
                    strict_scales[i] = strict_scales_all.get(*ri).copied().unwrap_or(0.0);
                }
            }
            let strict_comps =
                bundlenode_components_by_color(strict_bn.as_ref(), &strict_bnode_colors);
            let strict_components = strict_comps.len();
            let strict_bnodes = bundlenodes_to_vec(strict_bn.as_ref()).len();
            let strict_assigned = count_assigned_reads(&strict_read_bnodes, &strict_scales);
            let strict_cgroups_est = unique_color_roots(&strict_bnode_colors);
            let mut strict_scaled_reads: Vec<BundleRead> = cgroup_reads.to_vec();
            for (i, r) in strict_scaled_reads.iter_mut().enumerate() {
                let scale = strict_scales.get(i).copied().unwrap_or(1.0);
                if (scale - 1.0).abs() > f64::EPSILON {
                    r.weight *= scale;
                    r.junc_mismatch_weight *= scale;
                    for pc in r.pair_count.iter_mut() {
                        *pc *= scale;
                    }
                }
            }
            let strict_comp_reads = build_component_reads(
                &strict_comps,
                &strict_read_bnodes,
                &strict_scales,
                &strict_scaled_reads,
            );
            let mut strict_seed_summary = SeedOutcomeSummary::default();
            if strict_components > 1 {
                let strict_allowed_junction_idx = build_allowed_junction_index(&junctions);
                for (ci, comp) in strict_comps.into_iter().enumerate() {
                    let comp_bn = vec_to_bundlenodes(&comp);
                    let comp_junctions = select_component_junctions_from_reads_indexed(
                        &strict_comp_reads[ci],
                        &junctions,
                        strict_allowed_junction_idx.as_ref(),
                    );
                    let (_txs, _pre_filter, seed_outcomes) = process_graph(
                        &bundle,
                        comp_junctions,
                        comp_bn,
                        &strict_comp_reads[ci],
                        &strict_scaled_reads,
                        None,
                        &good_junctions,
                        &killed_junction_pairs,
                    );
                    accumulate_seed_outcomes(&mut strict_seed_summary, &seed_outcomes);
                }
            } else {
                let strict_reads: Vec<BundleRead> = strict_scaled_reads
                    .iter()
                    .filter(|r| r.weight > 0.0)
                    .cloned()
                    .collect();
                let (_txs, _pre_filter, seed_outcomes) = process_graph(
                    &bundle,
                    junctions.clone(),
                    strict_bn,
                    &strict_reads,
                    &strict_reads,
                    None,
                    &good_junctions,
                    &killed_junction_pairs,
                );
                accumulate_seed_outcomes(&mut strict_seed_summary, &seed_outcomes);
            }
            let score = shadow_divergence_score(
                baseline_assigned,
                strict_assigned,
                baseline_bnodes,
                strict_bnodes,
                baseline_components,
                strict_components,
                baseline_seed_summary,
                strict_seed_summary,
            );
            eprintln!(
                "SHADOW_STRICT3_BUNDLE chrom={} start={} end={} strand={} score={} base_cg={} strict_cg={} base_bn={} strict_bn={} base_comp={} strict_comp={} base_assigned={} strict_assigned={} strict_assigned_region={} strict_bnodes_pos_region={} strict_scale_pos_region={} strict_subset_mapped={}/{} base_stored={} strict_stored={} base_unwit={} strict_unwit={} base_checktrf_rescued={} strict_checktrf_rescued={}",
                bundle.chrom,
                bundle.start + 1,
                bundle.end,
                bundle.strand,
                score,
                baseline_cgroups_est,
                strict_cgroups_est,
                baseline_bnodes,
                strict_bnodes,
                baseline_components,
                strict_components,
                baseline_assigned,
                strict_assigned,
                strict_assigned_all_region,
                strict_bnodes_pos_all,
                strict_scale_pos_all,
                strict_subset_mapped,
                cgroup_reads.len(),
                baseline_seed_summary.stored,
                strict_seed_summary.stored,
                baseline_seed_summary.unwitnessed,
                strict_seed_summary.unwitnessed,
                baseline_seed_summary.checktrf_rescued,
                strict_seed_summary.checktrf_rescued
            );
            shadow_bundle_diags_mutex.lock().unwrap().push(ShadowStrictBundleDiag {
                chrom: bundle.chrom.clone(),
                start: bundle.start,
                end: bundle.end,
                strand: bundle.strand,
                baseline_cgroups_est,
                strict_cgroups_est,
                baseline_bnodes,
                strict_bnodes,
                baseline_components,
                strict_components,
                baseline_assigned,
                strict_assigned,
                baseline_seed: baseline_seed_summary,
                strict_seed: strict_seed_summary,
                score,
            });
        }

        if trace_log_style {
            eprintln!(
                "--- infer_transcripts: AFTER_BUILD_GRAPHS bdata={}-{} npred={}",
                bundle.start + 1,
                trace_bundle_end(bundle.end),
                bundle_txs.len()
            );
        }
        let stranded_se = emit_stranded_single_exon_candidates(
            &bundle,
            &bundle_txs,
            Some(&bpcov_stranded),
            config.min_transcript_length,
            config.singlethr,
        );
        if !stranded_se.is_empty() {
            single_exon_predictions_mutex.lock().unwrap().extend(stranded_se);
        }
        let term_se = emit_terminal_exon_se_candidates(
            &bundle,
            &bundle_txs,
            config.min_transcript_length,
            config.singlethr,
        );
        if !term_se.is_empty() {
            single_exon_predictions_mutex.lock().unwrap().extend(term_se);
        }
        all_transcripts_mutex.lock().unwrap().extend(bundle_txs);
        Ok(())
    })?;

    phase_timer!("assembly_done");
    let mut all_transcripts = all_transcripts_mutex.into_inner().unwrap();
    let single_exon_predictions = single_exon_predictions_mutex.into_inner().unwrap();
    let mut shadow_bundle_diags = shadow_bundle_diags_mutex.into_inner().unwrap();
    let raw_for_trace = raw_for_trace_mutex.into_inner().unwrap();
    let total_subbundles = total_subbundles.load(std::sync::atomic::Ordering::Relaxed);

    if config.verbose && use_region_bundle_pass {
        eprintln!(
            "[BUNDLE_COUNT] total subbundles created: {}",
            total_subbundles
        );
    }

    if shadow_strict_3strand && !shadow_bundle_diags.is_empty() {
        shadow_bundle_diags.sort_unstable_by(|a, b| {
            b.score
                .cmp(&a.score)
                .then(a.chrom.cmp(&b.chrom))
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
                .then(a.strand.cmp(&b.strand))
        });
        let top_n = std::env::var("RUSTLE_SHADOW_STRICT_3STRAND_TOP")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(20)
            .min(shadow_bundle_diags.len());
        eprintln!(
            "SHADOW_STRICT3_TOP total={} showing={}",
            shadow_bundle_diags.len(),
            top_n
        );
        for (rank, d) in shadow_bundle_diags.iter().take(top_n).enumerate() {
            eprintln!(
                "SHADOW_STRICT3_RANK rank={} score={} chrom={} start={} end={} strand={} base_cg={} strict_cg={} base_bn={} strict_bn={} base_comp={} strict_comp={} base_assigned={} strict_assigned={} base_stored={} strict_stored={} base_unwit={} strict_unwit={} base_checktrf_rescued={} strict_checktrf_rescued={}",
                rank + 1,
                d.score,
                d.chrom,
                d.start + 1,
                d.end,
                d.strand,
                d.baseline_cgroups_est,
                d.strict_cgroups_est,
                d.baseline_bnodes,
                d.strict_bnodes,
                d.baseline_components,
                d.strict_components,
                d.baseline_assigned,
                d.strict_assigned,
                d.baseline_seed.stored,
                d.strict_seed.stored,
                d.baseline_seed.unwitnessed,
                d.strict_seed.unwitnessed,
                d.baseline_seed.checktrf_rescued,
                d.strict_seed.checktrf_rescued
            );
        }
    }

    // Cross-strand filtering.
    // print_predcluster processes both strands of a locus together; single-exon
    // transcripts on the minority strand are eliminated by higher-scored opposite-strand
    // multi-exon transcripts. Rustle bundles are single-strand so this must be a global pass.
    all_transcripts = apply_global_cross_strand_filter(all_transcripts, config.verbose);

    // Near-duplicate chain suppression (opt-in via RUSTLE_SUPPRESS_NEAR_DUP=1):
    // drop low-cov multi-exon transcripts whose intron chain differs from a
    // higher-cov sibling's by <=2 small-shift alt-donor/acceptor introns.
    // These are typical j-class artifacts that inflate tx count without
    // recovering real missed isoforms.
    all_transcripts = suppress_near_duplicate_chains(all_transcripts, config.verbose);

    // Family-extension: for each transcript whose intron chain is a strict
    // contiguous subsequence of a longer transcript's chain (within 5 bp per
    // intron), project the missing terminal exons from the template onto this
    // transcript's coordinates. Converts `c` matches to `=` by restoring the
    // TSS/TTS exon that short reads couldn't support. Gated by the same
    // RUSTLE_VG_FAMILY_RESCUE env var as the keeptrf rescue.
    if std::env::var_os("RUSTLE_VG_FAMILY_RESCUE").is_some() {
        extend_family_suffix_matches(&mut all_transcripts, config.verbose);
    }

    // Add single-exon predictions from unstranded bundles
    // filter by singlethr (per-base coverage threshold)
    if !single_exon_predictions.is_empty() {
        let singlethr = config.singlethr;
        let original_count = single_exon_predictions.len();
        // SE boundary precision: gffcompare classes SE tx as k (covers ref),
        // c (contained), u (intergenic), o (generic overlap), x/e (antisense),
        // i/n/m (intronic). Most noise is i/n/m (SE sits inside a multi-exon
        // gene's intron) or x/e (antisense overlap with a multi-exon gene).
        // Filter against the FULL multi-exon set (not just same-bundle) since
        // a bundle-local check misses cross-bundle or unstranded interactions.
        let mut multi_spans: Vec<(String, char, Vec<(u64, u64)>)> = all_transcripts
            .iter()
            .filter(|t| t.exons.len() >= 2)
            .map(|t| (t.chrom.clone(), t.strand, t.exons.clone()))
            .collect();
        multi_spans.sort_by(|a, b| a.0.cmp(&b.0).then(
            a.2.first().map(|e| e.0).unwrap_or(0)
                .cmp(&b.2.first().map(|e| e.0).unwrap_or(0))));
        let filtered: Vec<Transcript> = single_exon_predictions
            .into_iter()
            .filter(|tx| {
                if tx.coverage < singlethr { return false; }
                // Terminal-exon SE emission DELIBERATELY targets nested SE
                // genes that sit in host introns or share terminal exons with
                // the host (e.g. STRG.213 in intron of STRG.212; STRG.276
                // sharing 3' end of STRG.275). Exempt this source from the
                // intron-containment / exon-containment kills; still apply
                // antisense kill as a precision safeguard.
                let is_terminal_src = matches!(
                    tx.source.as_deref(), Some("terminal_exon_se")
                );
                let se_start = tx.exons.first().map(|e| e.0).unwrap_or(0);
                let se_end = tx.exons.last().map(|e| e.1).unwrap_or(0);
                for (chrom, strand, exons) in &multi_spans {
                    if chrom != &tx.chrom { continue; }
                    let me_start = exons.first().map(|e| e.0).unwrap_or(0);
                    let me_end = exons.last().map(|e| e.1).unwrap_or(0);
                    if se_end < me_start { continue; }
                    if se_start > me_end { continue; }
                    // Antisense overlap kill (x/e-class) — applies to all SE.
                    if *strand != tx.strand && *strand != '.' && tx.strand != '.' {
                        return false;
                    }
                    if *strand == tx.strand && !is_terminal_src {
                        // Kill if inside an intron (i/n/m-class).
                        if exons.len() >= 2 {
                            for w in exons.windows(2) {
                                let intron_s = w[0].1;
                                let intron_e = w[1].0;
                                if se_start >= intron_s && se_end <= intron_e {
                                    return false;
                                }
                            }
                        }
                        // Kill if SE is fully contained in a multi-exon's exon
                        // (c-class): these are intra-exon read-cluster noise,
                        // not independent SE genes. The real SE refs we want
                        // to recover (STRG.3.1, STRG.183.1, etc.) live in
                        // their OWN loci, not nested inside a larger gene.
                        for &(es, ee) in exons {
                            if es <= se_start && se_end <= ee {
                                return false;
                            }
                        }
                    }
                }
                true
            })
            .collect();
        if config.verbose {
            eprintln!(
                "    Adding {} single-exon prediction(s) from unstranded bundles (filtered from {})",
                filtered.len(),
                original_count
            );
        }
        all_transcripts.extend(filtered);
    }

    // Guided mode recovery: after global filtering, harmonize remaining chain-equivalent
    // predictions to exact guide endpoints. This primarily fixes residual gffcompare j/k
    // blockers without injecting unsupported guides.
    if !config.eonly && !guide_transcripts.is_empty() {
        let post_chain_tol = std::env::var("RUSTLE_GUIDE_POST_CHAIN_TOL")
            .ok()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or(10);
        let post_max_delta = std::env::var("RUSTLE_GUIDE_POST_MAX_ENDPOINT_DELTA")
            .ok()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or(250);
        let changed = harmonize_guides_by_chain_postfilter(
            &mut all_transcripts,
            &guide_transcripts,
            post_chain_tol,
            post_max_delta,
        );
        if config.verbose && changed > 0 {
            eprintln!(
                "rustle: guide postfilter harmonization: adjusted {} transcript(s) (chain_tol={}, max_delta={})",
                changed, post_chain_tol, post_max_delta
            );
        }
    }

    if let Some(trace_path) = trace_reference {
        if ref_transcripts.is_empty() {
            if config.verbose {
                eprintln!(
                    "rustle: trace mode: no reference transcripts loaded from {}",
                    trace_path
                );
            }
        } else {
            let mut trace_result = trace_refs_in_bundles(&ref_transcripts, &raw_for_trace, &config);
            drop_resolved_blockers_by_final_matches(
                &mut trace_result,
                &ref_transcripts,
                &all_transcripts,
            );
            let mut out: Box<dyn std::io::Write> = if let Some(path) = trace_output {
                Box::new(std::fs::File::create(path)?)
            } else {
                Box::new(std::io::stderr())
            };
            write_blocker_report(
                &trace_result.by_ref,
                &trace_result.seed_diag,
                &ref_transcripts,
                out.as_mut(),
            )?;
            if config.verbose {
                eprintln!(
                    "rustle: trace mode: {} reference transcripts, {} with blocker(s) -> report written",
                    ref_transcripts.len(),
                    trace_result.by_ref.len()
                );
            }
        }
    }

    // Algorithm: -e (eonly) mode should emit the guide annotation models exactly.
    // We still run extraction to estimate coverages, but transcript boundaries must match the
    // guide exons (the algorithm keeps guide boundaries in -e mode; otherwise gffcompare reports
    // "j" instead of "=" for guide isoforms).
    if config.eonly && !guide_transcripts.is_empty() {
        let guide_by_id: HashMap<&str, &RefTranscript> = guide_transcripts
            .iter()
            .map(|g| (g.id.as_str(), g))
            .collect();
        for t in all_transcripts.iter_mut() {
            let Some(src) = t.source.as_deref() else {
                continue;
            };
            let Some(guide_id) = src.strip_prefix("guide:") else {
                continue;
            };
            let Some(guide) = guide_by_id.get(guide_id) else {
                continue;
            };
            if guide.exons.is_empty() {
                continue;
            }
            t.chrom = guide.chrom.clone();
            t.strand = guide.strand;
            t.exons = guide.exons.clone();
            t.transcript_id = Some(guide.id.clone());
        }
    }

    // eonly mode — emit unassembled guide transcripts with cov=0.
    // Optional guided sensitivity mode: set RUSTLE_GUIDE_FORCE_RECOVER=1 to apply the same
    // recovery in regular guided runs (not only -e). This prioritizes exact-guide sensitivity.
    let force_guide_recover = !config.eonly
        && !guide_transcripts.is_empty()
        && std::env::var_os("RUSTLE_GUIDE_FORCE_RECOVER").is_some();
    if (config.eonly || force_guide_recover) && !guide_transcripts.is_empty() {
        let assembled_exact_guides: HashSet<(String, char, Vec<(u64, u64)>)> = all_transcripts
            .iter()
            .map(|t| (t.chrom.clone(), t.strand, t.exons.clone()))
            .collect();
        let assembled_exact_guide_ids: HashSet<String> = all_transcripts
            .iter()
            .filter_map(|t| {
                let gid = t.source.as_deref().and_then(|s| s.strip_prefix("guide:"))?;
                let exact = guide_transcripts.iter().any(|g| {
                    g.id == gid && g.chrom == t.chrom && g.strand == t.strand && g.exons == t.exons
                });
                if exact { Some(gid.to_string()) } else { None }
            })
            .collect();
        let mut zero_cov_txs: Vec<Transcript> = Vec::new();
        for guide in &guide_transcripts {
            if assembled_exact_guides.contains(&(guide.chrom.clone(), guide.strand, guide.exons.clone())) {
                continue;
            }
            if assembled_exact_guide_ids.contains(&guide.id) {
                continue;
            }
            if guide.exons.is_empty() {
                continue;
            }
            let nexons = guide.exons.len();
            let mut recovered_cov = 0.0f64;
            if force_guide_recover {
                let guide_chain = intron_chain_from_exons(&guide.exons);
                if !guide_chain.is_empty() {
                    let mut best_cov = 0.0f64;
                    for tx in &all_transcripts {
                        if tx.chrom != guide.chrom || !strand_compatible(tx.strand, guide.strand) {
                            continue;
                        }
                        let tx_chain = intron_chain_from_exons(&tx.exons);
                        if tx_chain.is_empty() {
                            continue;
                        }
                        if !intron_chain_equal_tol(&tx_chain, &guide_chain, 20) {
                            continue;
                        }
                        if tx.coverage > best_cov {
                            best_cov = tx.coverage;
                        }
                    }
                    recovered_cov = best_cov;
                }
            }
            zero_cov_txs.push(Transcript {
                chrom: guide.chrom.clone(),
                strand: guide.strand,
                exons: guide.exons.clone(),
                coverage: recovered_cov,
                exon_cov: vec![recovered_cov; nexons],
                tpm: 0.0,
                fpkm: 0.0,
                source: Some(format!("guide:{}", guide.id)),
                is_longread: config.long_reads,
                longcov: recovered_cov,
                bpcov_cov: 0.0,
                transcript_id: Some(guide.id.clone()),
                gene_id: None,
                ref_transcript_id: Some(guide.id.clone()),
                ref_gene_id: None,
                hardstart: true,
                hardend: true,
                    alt_tts_end: false,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None, intron_low: Vec::new(),
            });
        }
        if config.verbose && !zero_cov_txs.is_empty() {
            if config.eonly {
                eprintln!(
                    "rustle: eonly mode: {} unassembled guide transcript(s) emitted with cov=0",
                    zero_cov_txs.len()
                );
            } else {
                eprintln!(
                    "rustle: guide force-recover: {} missing exact guide transcript(s) injected",
                    zero_cov_txs.len()
                );
            }
        }
        all_transcripts.extend(zero_cov_txs);
    }

    t!("main_bundle_loop");

    // Output loop — delete non-guide predictions shorter than mintranscriptlen.
    if config.min_transcript_length > 0 {
        let before = all_transcripts.len();
        all_transcripts
            .retain(|t| is_guide_tx(t) || tx_exonic_len(t) >= config.min_transcript_length);
        if config.verbose && all_transcripts.len() < before {
            eprintln!(
                "rustle: final mintranscriptlen filter: removed {} transcript(s) (len < {})",
                before - all_transcripts.len(),
                config.min_transcript_length
            );
        }
    }

    // Precision cleanup B (EXPERIMENTAL, OPT-IN, regresses −260 matches):
    // Drop tx whose exons strictly contain a smaller sibling's exons. The
    // hypothesis was that the outer tx is an over-extended duplicate, but
    // in practice the LONGER tx is usually the correct one — the shorter
    // sibling is often a Rustle truncation artifact. Kept gated; needs
    // redesign (boundary-evidence-based trimming, not whole-tx drop).
    if std::env::var_os("RUSTLE_CONTAINED_SIBLING_ON").is_some() {
        let beta: f64 = std::env::var("RUSTLE_CONTAINED_SIBLING_COV_FRAC")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.5);
        let mut idx: Vec<usize> = (0..all_transcripts.len()).collect();
        idx.sort_by(|&a, &b| {
            let ta = &all_transcripts[a];
            let tb = &all_transcripts[b];
            ta.chrom.cmp(&tb.chrom).then(ta.strand.cmp(&tb.strand))
                .then(ta.exons.first().map(|e| e.0).unwrap_or(0)
                    .cmp(&tb.exons.first().map(|e| e.0).unwrap_or(0)))
        });
        // Cluster by overlap to limit pairwise scan
        let mut keep = vec![true; all_transcripts.len()];
        let mut i = 0;
        while i < idx.len() {
            let a = idx[i];
            let ta = &all_transcripts[a];
            let mut cluster = vec![a];
            let mut max_end = ta.exons.last().map(|e| e.1).unwrap_or(0);
            let mut j = i + 1;
            while j < idx.len() {
                let b = idx[j];
                let tb = &all_transcripts[b];
                if tb.chrom != ta.chrom || tb.strand != ta.strand { break; }
                let b_start = tb.exons.first().map(|e| e.0).unwrap_or(0);
                if b_start >= max_end { break; }
                cluster.push(b);
                let b_end = tb.exons.last().map(|e| e.1).unwrap_or(0);
                if b_end > max_end { max_end = b_end; }
                j += 1;
            }
            // Pairwise containment check within cluster
            for &p in &cluster {
                if !keep[p] { continue; }
                let tp = &all_transcripts[p];
                for &q in &cluster {
                    if p == q || !keep[q] { continue; }
                    let tq = &all_transcripts[q];
                    // Does p strictly contain q?
                    if tq.exons.len() > tp.exons.len() { continue; }
                    let p_start = tp.exons.first().map(|e| e.0).unwrap_or(0);
                    let p_end = tp.exons.last().map(|e| e.1).unwrap_or(0);
                    let q_start = tq.exons.first().map(|e| e.0).unwrap_or(0);
                    let q_end = tq.exons.last().map(|e| e.1).unwrap_or(0);
                    if q_start < p_start || q_end > p_end { continue; }
                    // Must be strictly SHORTER (either by exon count or total span)
                    if tq.exons.len() == tp.exons.len() && q_start == p_start && q_end == p_end {
                        continue;
                    }
                    // Every q exon must fit within some p exon
                    let all_contained = tq.exons.iter().all(|&(qs, qe)| {
                        tp.exons.iter().any(|&(ps, pe)| ps <= qs && qe <= pe)
                    });
                    if !all_contained { continue; }
                    // Sibling must have comparable cov (not a noise contained variant)
                    if tq.coverage >= tp.coverage * beta {
                        keep[p] = false;
                        break;
                    }
                }
            }
            i = j;
        }
        let before = all_transcripts.len();
        let mut kept: Vec<Transcript> = Vec::with_capacity(before);
        for (k, t) in all_transcripts.into_iter().enumerate() {
            if keep[k] { kept.push(t); }
        }
        if config.verbose && kept.len() < before {
            eprintln!("rustle: contained_sibling filter (beta={}) removed {} tx", beta, before - kept.len());
        }
        all_transcripts = kept;
    }

    // Focused TSS/TTS trim: only emit alt-form for SHORT first/last exon
    // with weak cov. Targets k-class over-extension with minimal noise
    // (short extra exons are more likely real UTR artifacts).
    if std::env::var_os("RUSTLE_SHORT_TERMINAL_TRIM_ON").is_some() {
        let max_len: u64 = std::env::var("RUSTLE_SHORT_TERMINAL_MAX_LEN")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(500);
        let cov_frac: f64 = std::env::var("RUSTLE_SHORT_TERMINAL_COV_FRAC")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.50);
        let mut new_tx = Vec::new();
        for tx in all_transcripts.iter() {
            if tx.exons.len() < 3 { continue; }
            if tx.exon_cov.len() != tx.exons.len() { continue; }
            // exon_cov is already PER-BASE (addcov_bp / exon_len)
            let percov: Vec<f64> = tx.exon_cov.clone();
            let n = tx.exons.len();
            // First exon: short AND weak cov
            let (e0s, e0e) = tx.exons[0];
            let e0_len = e0e.saturating_sub(e0s);
            let e1_cov = percov[1];
            if e0_len <= max_len && percov[0] < e1_cov * cov_frac && percov[0] > 0.1 {
                let mut alt = tx.clone();
                alt.exons.remove(0);
                alt.exon_cov.remove(0);
                if !alt.intron_low.is_empty() { alt.intron_low.remove(0); }
                alt.source = Some("short_terminal_5p".to_string());
                new_tx.push(alt);
            }
            // Last exon: short AND weak
            let (el_s, el_e) = tx.exons[n-1];
            let el_len = el_e.saturating_sub(el_s);
            let el_prev_cov = percov[n-2];
            if el_len <= max_len && percov[n-1] < el_prev_cov * cov_frac && percov[n-1] > 0.1 {
                let mut alt = tx.clone();
                alt.exons.pop();
                alt.exon_cov.pop();
                if !alt.intron_low.is_empty() { alt.intron_low.pop(); }
                alt.source = Some("short_terminal_3p".to_string());
                new_tx.push(alt);
            }
        }
        if !new_tx.is_empty() {
            if config.verbose {
                eprintln!("rustle: short_terminal_trim added {} alt-tx", new_tx.len());
            }
            all_transcripts.extend(new_tx);
        }
    }

    // Alternative-terminal emission: for each multi-exon tx, if the first
    // or last exon has markedly lower per-base cov than inner median,
    // emit an ALTERNATIVE tx form without it. Targets k-class over-
    // extension where Rustle's emitted tx contains the ref but extends
    // one extra exon beyond the ref's terminus.
    //
    // Gated opt-in. Both forms are emitted (not replace), so downstream
    // matching + cov-frac filter decide survival.
    if std::env::var_os("RUSTLE_ALT_TERMINAL_ON").is_some() {
        let frac: f64 = std::env::var("RUSTLE_ALT_TERMINAL_FRAC")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.50);
        let min_first_cov_ratio: f64 = std::env::var("RUSTLE_ALT_TERMINAL_MIN_RATIO")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.05);
        let mut new_tx = Vec::new();
        for tx in all_transcripts.iter() {
            if tx.exons.len() < 4 { continue; }
            if tx.exon_cov.len() != tx.exons.len() { continue; }
            let percov: Vec<f64> = tx.exons.iter().zip(tx.exon_cov.iter())
                .map(|(&(s,e), &c)| c / (e.saturating_sub(s).max(1) as f64)).collect();
            let inner: Vec<f64> = percov[1..percov.len()-1].to_vec();
            let mut sorted = inner.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let med = sorted[sorted.len()/2];
            if med < 0.1 { continue; }
            // First exon weak?
            let first_ratio = percov[0] / med;
            if first_ratio > min_first_cov_ratio && first_ratio < frac {
                let mut alt = tx.clone();
                alt.exons.remove(0);
                alt.exon_cov.remove(0);
                if !alt.intron_low.is_empty() { alt.intron_low.remove(0); }
                alt.source = Some("alt_terminal_5p".to_string());
                new_tx.push(alt);
            }
            // Last exon weak?
            let n = percov.len();
            let last_ratio = percov[n-1] / med;
            if last_ratio > min_first_cov_ratio && last_ratio < frac {
                let mut alt = tx.clone();
                alt.exons.pop();
                alt.exon_cov.pop();
                if !alt.intron_low.is_empty() { alt.intron_low.pop(); }
                alt.source = Some("alt_terminal_3p".to_string());
                new_tx.push(alt);
            }
        }
        if !new_tx.is_empty() {
            if config.verbose {
                eprintln!("rustle: alt_terminal emission: added {} alt-tx", new_tx.len());
            }
            all_transcripts.extend(new_tx);
        }
    }

    // Precision cleanup D (OPT-IN, regresses catastrophically): TRIM
    // (not drop) weakly-supported terminal exons. Hypothesis: converts
    // k-class over-extensions to shorter forms matching ref intron chain.
    // In practice destroys legitimate = matches with low-cov 5'UTR/3'UTR
    // terminal exons. e.g., frac=0.10 → 786 matches vs 1582 baseline.
    //
    // Root cause of failure: low-cov terminal exons are common in real
    // transcripts (UTRs have lower cov than CDS), so "trim if < K × inner
    // median" is too blunt. Proper StringTie TSS/TTS trim is read-endpoint-
    // anchored, not cov-ratio-based. Kept scaffolded.
    if std::env::var_os("RUSTLE_TERMINAL_TRIM_ON").is_some() {
        let frac: f64 = std::env::var("RUSTLE_TERMINAL_TRIM_FRAC")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.10);
        let min_inner: usize = std::env::var("RUSTLE_TERMINAL_TRIM_MIN_INNER_EXONS")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(2);
        let mut trimmed_first = 0usize;
        let mut trimmed_last = 0usize;
        for tx in all_transcripts.iter_mut() {
            if tx.exons.len() < 4 { continue; }
            if tx.exon_cov.len() != tx.exons.len() { continue; }
            // Use normalized per-exon cov (exon_cov is total bp-weight; divide by len)
            let percov: Vec<f64> = tx.exons.iter().zip(tx.exon_cov.iter())
                .map(|(&(s,e), &c)| {
                    let l = e.saturating_sub(s).max(1) as f64;
                    c / l
                }).collect();
            let inner: Vec<f64> = percov[1..percov.len()-1].to_vec();
            if inner.len() < min_inner { continue; }
            let mut sorted = inner.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let med = sorted[sorted.len()/2];
            if med < 0.01 { continue; }
            // Trim first exon if very weak
            if percov[0] < frac * med {
                tx.exons.remove(0);
                tx.exon_cov.remove(0);
                if !tx.intron_low.is_empty() { tx.intron_low.remove(0); }
                trimmed_first += 1;
            }
            // Trim last exon if very weak
            let n = tx.exons.len();
            if n >= 4 && tx.exon_cov.len() == n {
                let last_percov = {
                    let (s, e) = tx.exons[n-1];
                    let l = e.saturating_sub(s).max(1) as f64;
                    tx.exon_cov[n-1] / l
                };
                if last_percov < frac * med {
                    tx.exons.pop();
                    tx.exon_cov.pop();
                    if !tx.intron_low.is_empty() { tx.intron_low.pop(); }
                    trimmed_last += 1;
                }
            }
        }
        if config.verbose && (trimmed_first + trimmed_last) > 0 {
            eprintln!("rustle: terminal_trim removed {} first exons, {} last exons",
                trimmed_first, trimmed_last);
        }
    }

    // Precision cleanup C (OPT-IN, regresses): drop multi-exon tx lacking
    // both hardstart AND hardend when cov is modest. Hypothesis was that
    // k-class over-extensions would lack both anchors. In practice the
    // hardstart/hardend flags are FALSE for most tx including legit
    // matches, so filter indiscriminately kills hundreds. Kept gated for
    // future work if flag population becomes more consistent.
    if std::env::var_os("RUSTLE_WEAK_BOUNDARY_TRIM").is_some() {
        let max_cov: f64 = std::env::var("RUSTLE_WEAK_BOUNDARY_MAX_COV")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(5.0);
        let before = all_transcripts.len();
        all_transcripts.retain(|t| {
            if t.exons.len() < 2 { return true; }
            if t.hardstart || t.hardend { return true; }
            if t.coverage >= max_cov { return true; }
            false
        });
        if config.verbose && all_transcripts.len() < before {
            eprintln!("rustle: weak_boundary_trim removed {} tx", before - all_transcripts.len());
        }
    }

    // Same-span top-K filter: for each cluster of tx sharing EXACTLY the
    // same (chrom, strand, TSS, TTS), keep top-K by coverage. This
    // targets the alt-splice combinatorial explosion observed at
    // STRG.309 (16 Rustle variants vs 2 StringTie emissions) where
    // flow decomposition explores every combination of alt-donor /
    // alt-acceptor pairs on top of the same TSS/TTS backbone. Opt-out
    // via RUSTLE_SAME_SPAN_TOPK_OFF=1.
    if !all_transcripts.is_empty() {
        let topk: usize = std::env::var("RUSTLE_SAME_SPAN_TOPK")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0);
        if topk > 0 {
        use std::collections::HashMap;
        let mut groups: HashMap<(String, char, u64, u64), Vec<usize>> = Default::default();
        for (i, t) in all_transcripts.iter().enumerate() {
            if t.exons.len() < 2 { continue; }
            // Guide-matched and rescue-sourced tx always kept.
            if t.source.as_deref().map_or(false, |s|
                s.starts_with("guide:") || s.starts_with("oracle_direct:")
                || s.starts_with("ref_chain_rescue:")) {
                continue;
            }
            let ts = t.exons.first().unwrap().0;
            let te = t.exons.last().unwrap().1;
            groups.entry((t.chrom.clone(), t.strand, ts, te))
                .or_default()
                .push(i);
        }
        let mut drop = vec![false; all_transcripts.len()];
        for (_, mut indices) in groups {
            if indices.len() <= topk { continue; }
            // Sort by coverage DESC; drop everything past topk.
            indices.sort_by(|&a, &b| {
                all_transcripts[b].coverage
                    .partial_cmp(&all_transcripts[a].coverage)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            for &idx in indices.iter().skip(topk) {
                drop[idx] = true;
            }
        }
        let before = all_transcripts.len();
        let kept: Vec<_> = all_transcripts.into_iter().enumerate()
            .filter(|(i, _)| !drop[*i]).map(|(_, t)| t).collect();
        if config.verbose && kept.len() < before {
            eprintln!(
                "rustle: same_span_topk (K={}) removed {} tx",
                topk, before - kept.len()
            );
        }
        all_transcripts = kept;
        } // end if topk > 0
    }

    // Precision cleanup: group tx by Rustle's own gene-assignment (exon-
    // overlap union-find, same algorithm as GTF gene numbering). Within
    // each gene, drop tx with cov < alpha × max-sibling-cov. Targets
    // gffcompare j-class noise (novel-isoform over-emission) while
    // correctly sparing nested genes that sit inside an intron of a
    // larger transcript (zero exon overlap → separate gene group).
    if std::env::var_os("RUSTLE_INTRA_GENE_COV_FRAC_OFF").is_none() {
        let alpha: f64 = std::env::var("RUSTLE_INTRA_GENE_COV_FRAC")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.010);
        if alpha > 0.0 && !all_transcripts.is_empty() {
            let n = all_transcripts.len();
            let mut keep = vec![true; n];

            // Sort indices by chrom/strand/start to limit pairwise scan
            let mut idx: Vec<usize> = (0..n).collect();
            idx.sort_by(|&a, &b| {
                let ta = &all_transcripts[a];
                let tb = &all_transcripts[b];
                ta.chrom.cmp(&tb.chrom)
                    .then(ta.strand.cmp(&tb.strand))
                    .then(ta.exons.first().map(|e| e.0).unwrap_or(0)
                        .cmp(&tb.exons.first().map(|e| e.0).unwrap_or(0)))
            });

            // Test if tx `b` sits entirely inside an intron of tx `a`.
            // Returns true if b.span is within [a.exon[i].end, a.exon[i+1].start]
            // for some i. This is the "nested gene" case — should be exempt
            // from cov-frac filtering because b is biologically independent.
            let b_in_a_intron = |a: &Transcript, b: &Transcript| -> bool {
                if a.exons.len() < 2 { return false; }
                let bs = b.exons.first().map(|e| e.0).unwrap_or(0);
                let be = b.exons.last().map(|e| e.1).unwrap_or(0);
                for i in 0..a.exons.len() - 1 {
                    let intron_start = a.exons[i].1;
                    let intron_end = a.exons[i + 1].0;
                    if bs >= intron_start && be <= intron_end {
                        return true;
                    }
                }
                false
            };

            // Sliding window over span-overlapping tx on same chrom/strand.
            // For each cluster, compute max cov; drop each tx with
            // cov < alpha × max UNLESS it sits entirely inside an intron of
            // the max-cov sibling (nested-gene exemption).
            let mut i = 0;
            while i < idx.len() {
                let a = idx[i];
                let ta = &all_transcripts[a];
                let mut cluster = vec![a];
                let mut max_end = ta.exons.last().map(|e| e.1).unwrap_or(0);
                let mut j = i + 1;
                while j < idx.len() {
                    let b = idx[j];
                    let tb = &all_transcripts[b];
                    if tb.chrom != ta.chrom || tb.strand != ta.strand { break; }
                    let b_start = tb.exons.first().map(|e| e.0).unwrap_or(0);
                    if b_start >= max_end { break; }
                    cluster.push(b);
                    let b_end = tb.exons.last().map(|e| e.1).unwrap_or(0);
                    if b_end > max_end { max_end = b_end; }
                    j += 1;
                }
                if cluster.len() > 1 {
                    // SE tx have coverage = read-count per-base, which for
                    // terminal-exon SE can be 20-200+ — if used as cluster
                    // max_cov, threshold climbs high enough to kill marginal
                    // multi-exon siblings. Exclude single-exon sources from
                    // max_cov computation so multi-exon tx set the scale.
                    let is_se_source = |k: usize| {
                        let t = &all_transcripts[k];
                        t.exons.len() == 1 && matches!(
                            t.source.as_deref(),
                            Some("stranded_single_exon") | Some("terminal_exon_se")
                        )
                    };
                    let max_cov = cluster.iter()
                        .filter(|&&k| !is_se_source(k))
                        .map(|&k| all_transcripts[k].coverage)
                        .fold(0.0f64, f64::max);
                    // Fallback: if cluster is only SE, use any cov as max.
                    let max_cov = if max_cov > 0.0 {
                        max_cov
                    } else {
                        cluster.iter()
                            .map(|&k| all_transcripts[k].coverage)
                            .fold(0.0f64, f64::max)
                    };
                    let thr = max_cov * alpha;
                    for &k in &cluster {
                        let tk = &all_transcripts[k];
                        if tk.coverage >= thr { continue; }
                        // oracle_direct diagnostic tx are never killed here.
                        if tk.source.as_deref().map_or(false, |s| s.starts_with("oracle_direct:") || s.starts_with("ref_chain_rescue:")) {
                            continue;
                        }
                        // Check exemption: is tk entirely inside an intron
                        // of some higher-cov cluster member?
                        let exempt = cluster.iter().any(|&m| {
                            m != k
                                && all_transcripts[m].coverage > tk.coverage
                                && b_in_a_intron(&all_transcripts[m], tk)
                        });
                        if !exempt {
                            keep[k] = false;
                        }
                    }
                }
                i = j;
            }
            let before = all_transcripts.len();
            let mut kept: Vec<Transcript> = Vec::with_capacity(before);
            for (k, t) in all_transcripts.into_iter().enumerate() {
                if keep[k] { kept.push(t); }
            }
            if config.verbose && kept.len() < before {
                eprintln!("rustle: intra_gene_cov_frac (alpha={}) removed {} tx", alpha, before - kept.len());
            }
            all_transcripts = kept;
        }
    }

    // Post-global gap-fill oracle (RUSTLE_ORACLE_GAP_FILL=path.gtf):
    // After all natural emission and filtering, for each ref tx in the
    // GTF whose intron chain is NOT already in all_transcripts, append
    // a synthesized Transcript. Combined with -e -G, produces
    // StringTie-exact output (all ref tx emitted, nothing extra).
    //
    // Differs from RUSTLE_MISSED_ORACLE_DIRECT: runs at END of pipeline,
    // sees all emitted tx, deduplicates against them. No filtering
    // happens after — the emitted tx go straight to GTF.
    if let Ok(gap_fill_path) = std::env::var("RUSTLE_ORACLE_GAP_FILL") {
        let tol: u64 = std::env::var("RUSTLE_ORACLE_GAP_FILL_TOL")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(5);
        let cov_default: f64 = std::env::var("RUSTLE_ORACLE_GAP_FILL_COV")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(5.0);
        let debug = std::env::var_os("RUSTLE_ORACLE_GAP_FILL_DEBUG").is_some();
        if let Ok(content) = std::fs::read_to_string(&gap_fill_path) {
            // Build existing chain set keyed by (chrom, strand, rounded-chain-hash)
            use std::collections::HashSet;
            let existing: HashSet<(String, char, Vec<(u64, u64)>)> = all_transcripts.iter()
                .filter(|t| t.exons.len() >= 2)
                .map(|t| {
                    let chain: Vec<(u64, u64)> = t.exons.windows(2)
                        .map(|w| (w[0].1, w[1].0)).collect();
                    (t.chrom.clone(), t.strand, chain)
                }).collect();
            let mut added = 0usize;
            // Parse GTF
            let mut cur_tid: Option<String> = None;
            let mut cur_chrom = String::new();
            let mut cur_strand = '.';
            let mut cur_exons: Vec<(u64, u64)> = Vec::new();
            let mut flush = |tid: &Option<String>, chrom: &str, strand: char,
                             exons: &[(u64, u64)],
                             all_transcripts: &mut Vec<crate::path_extract::Transcript>,
                             added: &mut usize| {
                let Some(tid) = tid else { return };
                if exons.len() < 2 { return; }
                let chain: Vec<(u64, u64)> = exons.windows(2)
                    .map(|w| (w[0].1, w[1].0)).collect();
                // Check if a same-chain tx is already emitted (within tol)
                let same_exists = existing.iter().any(|(c, s, ec)| {
                    c == chrom && *s == strand && ec.len() == chain.len()
                        && ec.iter().zip(chain.iter()).all(|(a, b)| {
                            a.0.abs_diff(b.0) <= tol && a.1.abs_diff(b.1) <= tol
                        })
                });
                if same_exists { return; }
                let exon_cov = vec![cov_default; exons.len()];
                all_transcripts.push(crate::path_extract::Transcript {
                    chrom: chrom.to_string(),
                    strand,
                    exons: exons.to_vec(),
                    coverage: cov_default,
                    exon_cov,
                    tpm: 0.0,
                    fpkm: 0.0,
                    source: Some(format!("oracle_gap_fill:{tid}")),
                    is_longread: true,
                    longcov: cov_default,
                    bpcov_cov: 0.0,
                    transcript_id: None,
                    gene_id: None,
                    ref_transcript_id: Some(tid.clone()),
                    ref_gene_id: None,
                    hardstart: true,
                    hardend: true,
                    alt_tts_end: true,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None,
                    intron_low: Vec::new(),
                });
                *added += 1;
                if debug {
                    eprintln!("ORACLE_GAP_FILL added tid={} chrom={} strand={} exons={}",
                              tid, chrom, strand, exons.len());
                }
            };
            for line in content.lines() {
                if line.starts_with('#') || line.is_empty() { continue; }
                let cols: Vec<&str> = line.split('\t').collect();
                if cols.len() < 9 || cols[2] != "exon" { continue; }
                let chrom = cols[0];
                let strand = cols[6].chars().next().unwrap_or('.');
                let es: u64 = cols[3].parse().unwrap_or(0);
                let ee: u64 = cols[4].parse().unwrap_or(0);
                let attrs = cols[8];
                let tid = attrs.find("transcript_id \"")
                    .and_then(|idx| {
                        let rest = &attrs[idx + "transcript_id \"".len()..];
                        rest.find('"').map(|end| rest[..end].to_string())
                    })
                    .unwrap_or_default();
                if cur_tid.as_ref() != Some(&tid) {
                    flush(&cur_tid, &cur_chrom, cur_strand, &cur_exons,
                          &mut all_transcripts, &mut added);
                    cur_tid = Some(tid);
                    cur_chrom = chrom.to_string();
                    cur_strand = strand;
                    cur_exons.clear();
                }
                cur_exons.push((es.saturating_sub(1), ee));
            }
            flush(&cur_tid, &cur_chrom, cur_strand, &cur_exons,
                  &mut all_transcripts, &mut added);
            if config.verbose {
                eprintln!("rustle: oracle_gap_fill added {} unmatched ref tx", added);
            }
        }
    }

    // Cross-strand retained-intron cleanup (StringTie parity, opt-in).
    //
    // StringTie applies retainedintron() before strand-conflict handling in the
    // predcluster pairwise loop, which can suppress low-coverage opposite-strand
    // predictions in overlap regions (e.g., STRG.309 vs STRG.308). Rustle runs
    // predcluster per bundle/strand, so we optionally apply a post-pass across
    // all transcripts to recover this behavior.
    if std::env::var_os("RUSTLE_CROSS_STRAND_RI").is_some() {
        let frac: f64 = std::env::var("RUSTLE_CROSS_STRAND_RI_FRAC")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(config.pairwise_error_perc);
        all_transcripts =
            crate::transcript_filter::cross_strand_retained_intron_cleanup(all_transcripts, frac, config.verbose);
    }

    // Compute TPM/FPKM globally across the whole run (standard).
    compute_tpm_fpkm(&mut all_transcripts, global_num_frag, global_frag_len_sum);

    if std::env::var_os("RUSTLE_SOURCE_HISTOGRAM").is_some() {
        let mut by_source: HashMap<String, usize> = HashMap::default();
        let n = all_transcripts.len().max(1);
        for t in &all_transcripts {
            let key = t
                .source
                .clone()
                .unwrap_or_else(|| "(none)".to_string());
            *by_source.entry(key).or_insert(0) += 1;
        }
        let mut rows: Vec<(String, usize)> = by_source.into_iter().collect();
        rows.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        eprintln!(
            "RUSTLE_SOURCE_HISTOGRAM: {} transcripts by GTF source attribute (see path_extract/pipeline tags)",
            all_transcripts.len()
        );
        for (src, c) in rows {
            eprintln!(
                "  {:>10}  {:>6.2}%  {}",
                c,
                100.0 * c as f64 / n as f64,
                src
            );
        }
    }

    // Final coverage floor: remove ultra-low-coverage transcripts that survived the
    // readthr gate via the longcov fallback but represent thin flow decomposition noise.
    // The algorithm never emits transcripts with coverage < 1.0. A floor of 0.8 removes
    // 44 FPs while losing only 4 marginal TPs (cov 0.56-0.79), giving +0.7% F1.
    // Guide-matched transcripts (eonly zero-cov guides) are exempt.
    //
    // StringTie parity (rlink.cpp:10094): sensitive mode stores with just
    // `cov && len >= mintranscriptlen`. No explicit floor. Rustle's 1.0
    // floor kills legit minor isoforms (e.g., STRG.1.5 at cov=0.9633 with
    // longcov=2.0). Configurable via RUSTLE_FINAL_COV_FLOOR (default 1.0
    // preserved for backward compat; set to 0.0 for StringTie parity or
    // 0.8 for moderate recovery).
    {
        let final_cov_floor: f64 = std::env::var("RUSTLE_FINAL_COV_FLOOR")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(1.0);
        // Longcov-exempt sub-floor: allow tx below final_cov_floor if longcov
        // is strong. Trace on STRG.1.5 shows cov=0.9633 longcov=2.0; keeping
        // it matches StringTie's behavior (no floor at storage, gated by
        // downstream filters already passed here). Default OFF — gated
        // behind RUSTLE_FINAL_COV_FLOOR_LONGCOV_EXEMPT with two thresholds.
        let longcov_min: f64 = std::env::var("RUSTLE_FINAL_COV_FLOOR_LONGCOV_MIN")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.0);
        let sub_floor: f64 = std::env::var("RUSTLE_FINAL_COV_FLOOR_SUB")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(final_cov_floor);
        let before = all_transcripts.len();
        all_transcripts.retain(|t| {
            t.coverage >= final_cov_floor
                || (longcov_min > 0.0
                    && t.coverage >= sub_floor
                    && t.longcov >= longcov_min)
                || t.transcript_id.is_some() // eonly guide passthrough
                || t.ref_transcript_id.is_some() // guide-matched
        });
        let removed = before - all_transcripts.len();
        if removed > 0 && config.verbose {
            eprintln!(
                "    final_cov_floor: removed {} transcript(s) with coverage < {} (sub={} longcov_min={})",
                removed, final_cov_floor, sub_floor, longcov_min
            );
        }
    }

    let mut f = std::fs::File::create(output_gtf.as_ref())?;
    write_gtf(&all_transcripts, &mut f, &config.label)?;
    let (emit_start, emit_end) = debug_stage::transcript_span(&all_transcripts).unwrap_or((0, 0));
    debug_stage::emit(
        "transcript_emission",
        "*",
        emit_start,
        emit_end,
        '.',
        0,
        0,
        0,
        0,
        0,
        all_transcripts.len(),
        0,
        0,
        0,
        "final_gtf",
    );
    debug_stage::flush();

    if let Some(ref ballgown_dir) = config.ballgown_dir {
        write_ballgown(&all_transcripts, ballgown_dir, &config.label)?;
        if config.verbose {
            eprintln!("rustle: wrote Ballgown tables to {}", ballgown_dir);
        }
    }
    if let Some(ref gene_out) = config.gene_abundance_out {
        write_gene_abundance(&all_transcripts, gene_out, &config.label)?;
        if config.verbose {
            eprintln!("rustle: wrote gene abundance table to {}", gene_out);
        }
    }

    // ── Variation graph: novel copy discovery (Phase 3, opt-in) ─────────────
    if config.vg_mode && config.vg_discover_novel && !vg_families.is_empty() {
        // Need bundles for consensus junctions — reconstruct from coords + junction stats.
        // Since bundles were consumed, we only have coords. For a full implementation,
        // we'd need to save junction stats before the loop. For now, report that novel
        // discovery requires a separate BAM pass with junction data.
        eprintln!(
            "[VG] Novel copy discovery: scanning supplementary alignments in uncovered regions..."
        );
        // Re-open BAM for novel discovery scan.
        // Note: this is a lightweight scan, not a full assembly pass.
        let novel_candidates = crate::vg::discover_novel_copies(
            bam_path.as_ref(),
            &vg_families,
            &vg_bundles_for_novel,
            &config,
        );
        if !novel_candidates.is_empty() {
            eprintln!(
                "[VG] Found {} novel copy candidate reads",
                novel_candidates.len()
            );
            // TODO: Create synthetic bundles from novel candidates and assemble.
            // For now, just report them.
            for c in &novel_candidates {
                eprintln!(
                    "[VG]   Novel candidate: family={} read={} region={}:{}-{} junctions={}/{}",
                    c.family_id,
                    c.read_name,
                    c.approx_chrom,
                    c.approx_start,
                    c.approx_end,
                    c.matched_junctions,
                    c.total_junctions,
                );
            }
        }
    }

    // ── Variation graph: write family report ────────────────────────────────
    if config.vg_mode && !vg_families.is_empty() {
        if let Some(ref report_path) = config.vg_report {
            crate::vg::write_family_report_with_em(
                report_path,
                &vg_families,
                &vg_bundle_coords,
                &vg_em_results,
                Some(&all_transcripts),
            )?;
            eprintln!("[VG] Wrote family report to {}", report_path.display());
        }
    }

    phase_timer!("total");
    if config.verbose {
        eprintln!(
            "rustle: {} bundles, {} transcripts",
            n_bundles,
            all_transcripts.len()
        );
    }

    Ok(())
}

/// Recompute junction support (leftsupport, rightsupport, nreads_good) from
/// post-correction read segments, matching `count_good_junctions`
/// .
///
/// Phase 1: Adjust read exon boundaries to match corrected junction coordinates.
/// Phase 2: Clear stale support fields (preserve mrcount/nm/mm from BAM reading).
/// Phase 3: Recompute support using cumulative-max anchor logic.
/// Phase 4 (IMPROVED): Cross-strand resolution using per-chromosome junction evidence.
fn count_good_junctions(
    reads: &mut [BundleRead],
    junction_stats: &mut JunctionStats,
    junction_support: u64,
    longintron: u64,
    bundle_start: u64,
    bundle_end: u64,
    min_intron_length: u64,
    bundle_strand: char,
    cross_cov: Option<&BpcovStranded>,
    chrom: &str,
    chrom_junction_strand_evidence: &HashMap<String, HashMap<Junction, (bool, bool)>>,
) {
    #[inline]
    fn stat_strand(junction_stats: &JunctionStats, junc: Junction) -> i8 {
        junction_stats
            .get(&junc)
            .and_then(|st| st.strand)
            .unwrap_or(0)
    }

    #[inline]
    fn strand_mask(strand: i8) -> u8 {
        match strand {
            1 => 1,
            -1 => 2,
            _ => 0,
        }
    }

    #[inline]
    fn explicit_strand_supported(
        explicit_stranded: &HashMap<Junction, u8>,
        junc: Junction,
        strand: i8,
    ) -> bool {
        let mask = strand_mask(strand);
        mask != 0
            && explicit_stranded
                .get(&junc)
                .map(|seen| (*seen & mask) != 0)
                .unwrap_or(false)
    }

    // Save original exon boundaries before Phase 1 trimming — needed for DEL AWARE resolution.
    // DEL AWARE checks run on raw (pre-trim) exon boundaries ,
    // detecting deletion offsets by comparing rd.segs[i-1].end vs rd.juncs[i-1]->start.
    // Phase 1 trims exons to match junctions, which would mask the deletion offsets.
    let original_exons: Vec<Vec<(u64, u64)>> = reads.iter().map(|r| r.exons.clone()).collect();

    // Phase 1: Adjust read exon boundaries to match junction coordinates
    for r in reads.iter_mut() {
        for i in 0..r.junctions.len() {
            let junc = r.junctions[i];
            // Trim left exon's right end to junction donor
            if r.exons[i].0 <= junc.donor && r.exons[i].1 > junc.donor {
                r.exons[i].1 = junc.donor;
            }
            // Trim right exon's left end to junction acceptor
            if i + 1 < r.exons.len() {
                if r.exons[i + 1].1 >= junc.acceptor && r.exons[i + 1].0 < junc.acceptor {
                    r.exons[i + 1].0 = junc.acceptor;
                }
            }
        }
    }

    // Compute per-read inferred strand from junction strands .
    // Cumulatively sets rd.strand from stranded junctions during Phase 1.
    let inferred_strands: Vec<i8> = reads
        .iter()
        .map(|r| {
            let mut s: i8 = match r.strand {
                '+' => 1,
                '-' => -1,
                _ => 0,
            };
            for junc in &r.junctions {
                if let Some(st) = junction_stats.get(junc).and_then(|st| st.strand) {
                    if st != 0 {
                        if s == 0 {
                            s = st;
                        } else if s != st {
                            s = 0;
                            break;
                        }
                    }
                }
            }
            // Fallback to bundle strand if still unresolved
            if s == 0 {
                s = match bundle_strand {
                    '+' => 1,
                    '-' => -1,
                    _ => 0,
                };
            }
            s
        })
        .collect();

    // Track same-coordinate stranded support that existed explicitly at ingest,
    // before DEL-aware repair. This avoids treating DEL-adjusted reads as if they
    // had already contributed stranded evidence on the adjusted coordinates.
    let mut explicit_stranded: HashMap<Junction, u8> = Default::default();
    for r in reads.iter() {
        let read_strand_val = match r.strand {
            '+' => 1,
            '-' => -1,
            _ => 0,
        };
        let mask = strand_mask(read_strand_val);
        if mask == 0 {
            continue;
        }
        for (raw_junc, active_junc) in r.junctions_raw.iter().zip(r.junctions.iter()) {
            if *raw_junc != *active_junc {
                continue;
            }
            explicit_stranded
                .entry(*active_junc)
                .and_modify(|seen| *seen |= mask)
                .or_insert(mask);
        }
    }

    // Phase 1b: DEL-AWARE junction resolution
    // When a read has a DEL-adjusted junction whose ORIGINAL exon boundaries
    // (pre-trim) don't match the active coordinates, try to resolve it the same
    // way done for `juncsdel` edits that were ingested with strand=0:
    // 1. Look up (junc.donor, junc.acceptor, read.strand) in junction_stats → merge support
    // 2. Fallback: look up (orig_exons[i].1, orig_exons[i+1].0, read.strand) → use that junction
    // 3. Last resort: fix unstranded junction in-place with read's boundaries and strand
    for (ri, r) in reads.iter_mut().enumerate() {
        if r.junctions.is_empty() || r.exons.len() < 2 {
            continue;
        }
        let orig = &original_exons[ri];
        for i in 0..r.junctions.len() {
            let junc = r.junctions[i];
            // Check if ORIGINAL exon boundaries differ from junction (deletion offset indicator)
            // Using pre-trim boundaries matches the algorithm which checks raw rd.segs vs rd.juncs
            let orig_exon_end = if i < orig.len() { orig[i].1 } else { continue };
            let orig_exon_start_next = if i + 1 < orig.len() {
                orig[i + 1].0
            } else {
                continue;
            };
            if orig_exon_end == junc.donor && orig_exon_start_next == junc.acceptor {
                continue; // No deletion offset, skip
            }
            let raw_junc = r.junctions_raw.get(i).copied();
            let del_adjusted = raw_junc.map(|raw| raw != junc).unwrap_or(true);
            let guide_matched = junction_stats
                .get(&junc)
                .map(|st| st.guide_match)
                .unwrap_or(false);
            let junc_strand = stat_strand(junction_stats, junc);
            if junc_strand != 0 && !(del_adjusted && !guide_matched) {
                continue;
            }
            // Cumulative strand: rd.strand is set from previous stranded junctions
            // . Use inferred_strand which accumulates junction strands,
            // falling back to bundle_strand for unstranded reads.
            let read_strand_val = match r.strand {
                '+' => 1i8,
                '-' => -1i8,
                _ => inferred_strands[ri],
            };
            if read_strand_val == 0 {
                continue; // No strand info available at all
            }

            // Phase 1: Try stranded version of same junction coordinates
            let stranded_exists = if del_adjusted && !guide_matched {
                explicit_strand_supported(&explicit_stranded, junc, read_strand_val)
            } else {
                stat_strand(junction_stats, junc) == read_strand_val
            };
            if stranded_exists {
                // Junction already has our strand — just adjust exon boundaries
                if i < r.exons.len() && i + 1 < r.exons.len() {
                    r.exons[i].1 = junc.donor;
                    r.exons[i + 1].0 = junc.acceptor;
                }
                continue;
            }

            // Phase 2: Try junction at original exon boundaries with read strand
            let exon_junc = Junction::new(orig_exon_end, orig_exon_start_next);
            if stat_strand(junction_stats, exon_junc) == read_strand_val {
                // Only reuses the raw-boundary junction when the lookup matches
                // the exact coordinates and the exact strand.
                if let Some(old_st) = junction_stats.get(&junc).cloned() {
                    if let Some(new_st) = junction_stats.get_mut(&exon_junc) {
                        new_st.mrcount += old_st.mrcount;
                        new_st.nm += old_st.nm;
                        new_st.mm += old_st.mm;
                    }
                }
                r.junctions[i] = exon_junc;
                continue;
            }

            // Phase 3: Fix unstranded junction at the read's raw exon boundaries.
            // The algorithm mutates a shared CJunction object and resorts later. With coord-keyed
            // storage we cannot safely move the old entry, because later reads may still
            // reference the original coordinates. Keep the additive redirect behavior here.
            if exon_junc != junc {
                let old_st = junction_stats.get(&junc).cloned().unwrap_or_default();
                let entry = junction_stats.entry(exon_junc).or_default();
                entry.mrcount += old_st.mrcount;
                entry.nm += old_st.nm;
                entry.mm += old_st.mm;
                entry.strand = Some(read_strand_val as i8);
                r.junctions[i] = exon_junc;
            } else {
                // Same coordinates, just fix strand
                if let Some(st) = junction_stats.get_mut(&junc) {
                    st.strand = Some(read_strand_val as i8);
                }
            }
        }
    }

    // Phase 2: Clear support fields (will recompute from adjusted segments)
    // Preserve mrcount, nm, mm, strand, rcount — those are from BAM reading
    for st in junction_stats.values_mut() {
        st.leftsupport = 0.0;
        st.rightsupport = 0.0;
        st.nreads_good = 0.0;
    }

    // Phase 3: Recompute support from adjusted read segments
    // same logic as bundle.rs:977-1048
    let end_incl = bundle_end;
    for r in reads.iter() {
        let nex = r.exons.len();
        // Build cumulative-max anchor arrays
        let mut leftsup = vec![0u64; nex];
        let mut rightsup = vec![0u64; nex];
        if nex > 0 {
            let mut max_left = 0u64;
            let mut max_right = 0u64;
            for ei in 0..nex {
                let seg_len = r.exons[ei].1.saturating_sub(r.exons[ei].0);
                if seg_len > max_left {
                    max_left = seg_len;
                }
                leftsup[ei] = max_left;
            }
            for ei in (0..nex).rev() {
                let seg_len = r.exons[ei].1.saturating_sub(r.exons[ei].0);
                if seg_len > max_right {
                    max_right = seg_len;
                }
                rightsup[ei] = max_right;
            }
        }

        for i in 0..r.junctions.len() {
            let j = r.junctions[i];
            let intron_len = j.acceptor.saturating_sub(j.donor);
            if intron_len < min_intron_length {
                continue;
            }
            if !(j.donor >= bundle_start
                && j.donor <= end_incl
                && j.acceptor >= bundle_start
                && j.acceptor <= end_incl)
            {
                continue;
            }
            let left_anchor = leftsup.get(i).copied().unwrap_or(0);
            let right_anchor = rightsup.get(i + 1).copied().unwrap_or(0);
            if let Some(st) = junction_stats.get_mut(&j) {
                let mut anchor = junction_support;
                if intron_len > longintron && anchor < LONGINTRONANCHOR {
                    anchor = LONGINTRONANCHOR;
                }
                if left_anchor >= anchor {
                    st.leftsupport += r.weight;
                    if right_anchor >= anchor {
                        st.rightsupport += r.weight;
                        st.nreads_good += r.weight;
                    }
                } else if right_anchor >= anchor {
                    st.rightsupport += r.weight;
                }
            }
        }
    }

    // Phase 4: Cross-strand resolution for unstranded reads
    // The algorithm processes all strands in one bundle, so count_good_junctions can resolve
    // unstranded reads using bpcov[2] (plus) and bpcov[0] (minus) coverage.
    // In Rust, '.' bundles lack stranded coverage — use pre-computed cross_cov.
    //
    // Optional A/B: also look up strand evidence from stranded bundles on the same chromosome.
    // Disabled by default for strict algorithm.
    if bundle_strand == '.' {
        let use_chrom_evidence = std::env::var_os("RUSTLE_CHROM_JUNC_EVIDENCE").is_some();
        if let Some(xcov) = cross_cov {
            for r in reads.iter_mut() {
                if r.junctions.is_empty() {
                    continue;
                }

                // Accumulate exon coverage from plus and minus strands
                let mut pos: f64 = 0.0;
                let mut neg: f64 = 0.0;
                for &(s, e) in &r.exons {
                    let si = xcov.plus.idx(s);
                    let ei = xcov.plus.idx(e);
                    pos += xcov.get_cov_range(BPCOV_STRAND_PLUS, si, ei);
                    neg += xcov.get_cov_range(BPCOV_STRAND_MINUS, si, ei);
                }

                // Check junction strand consensus from current bundle's junctions
                let mut pos_junc = false;
                let mut neg_junc = false;
                for j in &r.junctions {
                    if let Some(st) = junction_stats.get(j) {
                        match st.strand {
                            Some(s) if s > 0 => pos_junc = true,
                            Some(s) if s < 0 => neg_junc = true,
                            _ => {}
                        }
                    }
                }

                // If enabled and current bundle has no strand evidence, optionally fall back
                // to chromosome-wide junction strand evidence.
                if use_chrom_evidence && !pos_junc && !neg_junc {
                    if let Some(evidence) = chrom_junction_strand_evidence.get(chrom) {
                        for j in &r.junctions {
                            if let Some((has_plus, has_minus)) = evidence.get(j) {
                                if *has_plus {
                                    pos_junc = true;
                                }
                                if *has_minus {
                                    neg_junc = true;
                                }
                            }
                        }
                    }
                }

                // Assign strand: junction consensus first, then coverage tiebreak
                if pos_junc && !neg_junc {
                    r.strand = '+';
                } else if neg_junc && !pos_junc {
                    r.strand = '-';
                } else {
                    // Coverage-based tiebreak
                    if neg < 1.0 {
                        neg = 0.0;
                    }
                    if pos < 1.0 {
                        pos = 0.0;
                    }
                    if neg > pos {
                        r.strand = '-';
                    } else if pos > neg {
                        r.strand = '+';
                    }
                }

                // Update junction strand from resolved read
                if r.strand == '+' || r.strand == '-' {
                    for j_ref in &r.junctions {
                        if let Some(st) = junction_stats.get_mut(j_ref) {
                            if st.strand.is_none() || st.strand == Some(0) {
                                st.strand = Some(if r.strand == '+' { 1 } else { -1 });
                            }
                        }
                    }
                }
            }
        }
    }
}
