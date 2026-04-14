//! TSS/TTS clustering and trim points (C++ reference cluster_positions_with_counts, add_cpas_trimpoint).

/// C++ reference CPAS_POS_BIN (merge positions within this many bp).
pub const CPAS_POS_BIN: u64 = 5;
/// C++ reference CPAS_MIN_SUPPORT (minimum weight in a cluster to emit).
pub const CPAS_MIN_SUPPORT: u64 = 20;

/// Cluster (pos, count) into representative (center, total_count) (C++ reference cluster_positions_with_counts).
/// Positions within `window` bp are merged; center = round(weighted average); clusters with total count < min_support are dropped.
pub fn cluster_positions_with_counts(
    positions: &[(u64, f64)],
    window: u64,
    min_support: u64,
) -> Vec<(u64, f64)> {
    if positions.is_empty() {
        return Vec::new();
    }
    let mut raw: Vec<(u64, f64)> = positions.to_vec();
    raw.sort_by_key(|&(p, _)| p);

    let mut i = 0usize;
    let (mut pos, mut w) = (raw[0].0, raw[0].1);
    while i + 1 < raw.len() && raw[i + 1].0 == pos {
        i += 1;
        w += raw[i].1;
    }
    let mut b = pos;
    let mut sumw = w;
    let mut sumwx = w * (pos as f64);
    let mut out = Vec::new();

    let bin = if window > 0 { window } else { CPAS_POS_BIN };

    i += 1;
    while i < raw.len() {
        pos = raw[i].0;
        w = raw[i].1;
        while i + 1 < raw.len() && raw[i + 1].0 == pos {
            i += 1;
            w += raw[i].1;
        }
        if pos.saturating_sub(b) <= bin {
            b = pos;
            sumw += w;
            sumwx += w * (pos as f64);
        } else {
            if sumw >= min_support as f64 {
                let center = (sumwx / sumw).round() as u64;
                out.push((center, sumw));
            }
            b = pos;
            sumw = w;
            sumwx = w * (pos as f64);
        }
        i += 1;
    }
    if sumw >= min_support as f64 {
        let center = (sumwx / sumw).round() as u64;
        out.push((center, sumw));
    }
    out
}

/// Merge or append (pos, cov) into list; if existing point within window, merge (average pos, add cov) (C++ reference add_cpas_trimpoint).
pub fn add_cpas_trimpoint(list: &mut Vec<(u64, f64)>, pos: u64, cov: f64, window: u64) {
    let win = if window > 0 { window } else { CPAS_POS_BIN };
    for (p, c) in list.iter_mut() {
        if p.abs_diff(pos) <= win {
            *p = (*p + pos) / 2;
            *c += cov;
            return;
        }
    }
    list.push((pos, cov));
}
