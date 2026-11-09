//! r1600 (`docs/PREREG_assign_router_2026-11-09.md`): read every window's records in one pass over the BAM.
//!
//! The BAI chunks of all windows are merged into non-overlapping compressed ranges, which are split into contiguous
//! groups and decoded by one worker each. Each record overlapping at least one window is converted once, and the
//! result goes to every window it overlaps. Groups are joined in group order, so each window's list is in file
//! order — the order its own indexed query returns — at any worker count. The overlap test is noodles' query test:
//! same reference, and the record's [alignment_start, alignment_end] intersects the window's region interval
//! (`"{chrom}:{lo+1}-{hi}"`); a record without a start or an end matches no window.
use crate::denovo_assemble::BamIndexCache;
use anyhow::{anyhow, Context, Result};
use noodles_bgzf::VirtualPosition;
use noodles_core::region::Interval;
use noodles_csi::BinningIndex;
use noodles_sam::alignment::Record as _;
use std::collections::HashMap;

/// A window's search key on its reference: (interval start, interval end, window index).
type Key = (u64, u64, usize);

/// Every window's records, in input order, each list in file order. `convert` sees each routed record once, with the
/// sorted indices of the windows it overlaps; its `Some(t)` is given to each of them, `None` to none.
pub fn route<T: Clone + Send>(
    bam_path: &str,
    cache: &BamIndexCache,
    windows: &[(String, u64, u64)],
    workers: usize,
    convert: &(dyn Fn(&noodles_bam::Record, &[usize]) -> Result<Option<T>> + Sync),
) -> Result<Vec<Vec<T>>> {
    let header = cache.header();
    let mut intervals: Vec<Interval> = Vec::with_capacity(windows.len());
    let mut by_ref: HashMap<usize, Vec<Key>> = HashMap::new();
    let mut chunks = Vec::new();
    for (w, (chrom, lo, hi)) in windows.iter().enumerate() {
        let region: noodles_core::Region =
            format!("{chrom}:{}-{hi}", lo + 1).parse().with_context(|| format!("routing window {chrom}:{lo}-{hi}"))?;
        let id = header
            .reference_sequences()
            .get_index_of(region.name())
            .ok_or_else(|| anyhow!("routing window {chrom}:{lo}-{hi}: invalid reference sequence name"))?;
        let iv = region.interval();
        chunks.extend(cache.index().query(id, iv).with_context(|| format!("routing window {chrom}:{lo}-{hi}"))?);
        let s = iv.start().map_or(1, |p| p.get() as u64);
        let e = iv.end().map_or(u64::MAX, |p| p.get() as u64);
        by_ref.entry(id).or_default().push((s, e, w));
        intervals.push(iv);
    }
    // per reference: windows by start, and the running max of their ends (bounds the backward scan)
    let mut prefix_max: HashMap<usize, Vec<u64>> = HashMap::new();
    for (id, v) in by_ref.iter_mut() {
        v.sort_unstable();
        let mut m = 0;
        prefix_max.insert(*id, v.iter().map(|k| { m = m.max(k.1); m }).collect());
    }
    chunks.sort_by_key(|c| (c.start(), c.end()));
    let mut ranges: Vec<(VirtualPosition, VirtualPosition)> = Vec::new();
    for c in chunks {
        match ranges.last_mut() {
            Some(r) if c.start() <= r.1 => r.1 = r.1.max(c.end()),
            _ => ranges.push((c.start(), c.end())),
        }
    }
    // contiguous groups of about equal compressed bytes, cut at range boundaries
    let bytes = |r: &(VirtualPosition, VirtualPosition)| r.1.compressed().saturating_sub(r.0.compressed()) + 1;
    let total: u64 = ranges.iter().map(bytes).sum();
    let n = workers.max(1) as u64;
    let mut groups: Vec<&[(VirtualPosition, VirtualPosition)]> = Vec::new();
    let (mut start, mut acc, mut k) = (0usize, 0u64, 1u64);
    for (i, r) in ranges.iter().enumerate() {
        acc += bytes(r);
        if acc * n >= total * k && k < n {
            groups.push(&ranges[start..=i]);
            start = i + 1;
            k += 1;
        }
    }
    groups.push(&ranges[start..]);

    let matches = |rec: &noodles_bam::Record| -> std::io::Result<Vec<usize>> {
        let (Some(id), Some(s), Some(e)) =
            (rec.reference_sequence_id().transpose()?, rec.alignment_start().transpose()?, rec.alignment_end().transpose()?)
        else {
            return Ok(Vec::new());
        };
        let Some(v) = by_ref.get(&id) else { return Ok(Vec::new()) };
        let pm = &prefix_max[&id];
        let aln: Interval = (s..=e).into();
        let (s, e) = (s.get() as u64, e.get() as u64);
        let mut out = Vec::new();
        let mut j = v.partition_point(|k| k.0 <= e);
        while j > 0 {
            j -= 1;
            if pm[j] < s {
                break;
            }
            if intervals[v[j].2].intersects(aln) {
                out.push(v[j].2);
            }
        }
        out.sort_unstable();
        Ok(out)
    };
    let work = |group: &[(VirtualPosition, VirtualPosition)]| -> Result<Vec<(Vec<usize>, T)>> {
        let mut reader = noodles_bam::io::Reader::new(std::fs::File::open(bam_path).with_context(|| format!("opening {bam_path}"))?);
        let mut rec = noodles_bam::Record::default();
        let mut out = Vec::new();
        for &(s, e) in group {
            reader.get_mut().seek(s).with_context(|| format!("seeking to compressed offset {}", s.compressed()))?;
            loop {
                let at = reader.get_ref().virtual_position();
                if at >= e {
                    break;
                }
                let ctx = || format!("reading the record at compressed offset {}", at.compressed());
                if reader.read_record(&mut rec).with_context(ctx)? == 0 {
                    break;
                }
                let ws = matches(&rec).with_context(ctx)?;
                if ws.is_empty() {
                    continue;
                }
                if let Some(t) = convert(&rec, &ws).with_context(|| {
                    let (c, lo, hi) = &windows[ws[0]];
                    format!("routing window {c}:{lo}-{hi} at compressed offset {}", at.compressed())
                })? {
                    out.push((ws, t));
                }
            }
        }
        Ok(out)
    };
    let parts: Vec<Result<Vec<(Vec<usize>, T)>>> = std::thread::scope(|sc| {
        let hs: Vec<_> = groups.iter().map(|g| sc.spawn(|| work(g))).collect();
        hs.into_iter().map(|h| h.join().expect("router worker panicked")).collect()
    });
    let mut lists: Vec<Vec<T>> = vec![Vec::new(); windows.len()];
    for part in parts {
        for (ws, t) in part? {
            let (last, rest) = ws.split_last().expect("a routed record has a window");
            for &w in rest {
                lists[w].push(t.clone());
            }
            lists[*last].push(t);
        }
    }
    Ok(lists)
}
