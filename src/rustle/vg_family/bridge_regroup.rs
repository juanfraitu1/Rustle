//! ⭐ BRIDGE-AWARE REGROUPING of the assembled GTF: `copy_assign --assemble-only --bridge-regroup f1|f1v2`.
//!
//! **STATUS:** OPT-IN  (docs/MODULE_STATUS.md; `copy_assign --bridge-regroup`, default `off`, `--assemble-only` only; driver `RUSTLE_BRIDGE_REGROUP=f1|f1v2`)
//!
//! The port of two frozen post-processors of the emitted GTF: **F1** = `bench/f1_bridge.py --mode full` (sha1
//! 37ee8e77, `docs/PREREG_f1_bridge_locus_2026-09-28.md` §1) and **F1v2** = F1 filtered by `f1v2.py --rule min` (sha1
//! b4e788ad, `docs/PREREG_f1v2_readshare_2026-09-29.md` §1). Held out, F1 was EFFECTIVE on both gorilla samples but cut
//! ~40 annotated genes per sample, and on human A119b its bridge splits were 219 SEP / 451 FRAG; F1v2 was EFFECTIVE on
//! both human libraries. The fused-locus gain is fusions moved into explicit `fusion_of` relation records, not
//! removed: counting each bridge as its own locus, F1v2 is not better than `--gtf-regroup` (the F1v2 Outcome's
//! independent verification). Read both Outcomes before quoting either. On the held-out samples' BAMs the port
//! writes the scripts' GTFs, families GTFs and side tables byte for byte (`bench/ASSEMBLY_POLISH.md`, 2026-09-29
//! addendum 3).
//!
//! WHY. One readthrough transcript joins a gene and its downstream neighbour into one `gene_id` through shared
//! junctions; the most-read transcript then decides the family of the whole locus. A BRIDGE is such a linking
//! transcript whose two sides each carry read-proven ends of their own.
//!
//! RULE (F1). Per input `gene_id` g and strand, for every intron J = [s, e] (1-based closed) of g's transcripts:
//! * T_J = g's transcripts using J; R_J = g's other transcripts on J's strand, in components of same-strand exon
//!   overlap (>= 1 shared base, `--gtf-regroup`'s adjacency). A component is UP when it has bases left of the intron
//!   and none right of it, DOWN when the reverse (transcript orientation: mirrored on `-`), STRADDLE when both,
//!   INSIDE when neither.
//! * STRUCTURAL(J): no STRADDLE, >= 1 UP and >= 1 DOWN (T_J is then the only link); impossible with < 3 transcripts in
//!   g or < 2 in R_J.
//! * UP-proof(J): the U population (the DEDUPLICATED spliced primaries of J's strand whose 5' end is upstream of the
//!   donor and whose 3' end lies inside the intron) holds >= 1 PAS-proven 3' cluster, by `--polish-tes`'s own cluster
//!   rule, which the caller passes in ([`EndCluster`]).
//! * DOWN-proof(J): V1(J) >= 1, the readthrough filter's V1 over ALL spliced primaries ([`rt_v1`]).
//! * BRIDGE(J) = STRUCTURAL and UP-proof and DOWN-proof, every J judged on the original locus.
//!
//! F1v2 keeps a bridge junction only when also MINORITY(J): reads(T_J) < reads(UP_J) and reads(T_J) < reads(DOWN_J),
//! with reads(X) the `reads` attribute summed over X's transcripts and UP_J / DOWN_J the union of the UP / DOWN
//! components; equivalently share = reads(T_J) / (reads(T_J) + min(reads(UP_J), reads(DOWN_J))) < 1/2. A tie abstains.
//!
//! REGROUP. The bridge transcripts of g are the union of T_J over its (kept) bridge junctions. The others split into
//! same-contig, same-strand exon-overlap pieces named as `--gtf-regroup` (RG3) names them: the piece whose
//! representative max(reads, span, -line) is best keeps g, the others become `<g>.rg<k>`. Without a bridge the names
//! are therefore RG3's, by construction: [`rg3_pieces`] is the one function both flags call (`--gtf-regroup` on every
//! transcript of g, this pass on the non-bridge ones), over the one [`parse`]. Bridges group by exon overlap into
//! `<g>.fus<k>` (k = 1.. in line order). Their `transcript`
//! lines get `fusion_of "<piece>,..."` (the pieces they overlap, 5' to 3') and `fusion_junction
//! "<chrom>:<s>-<e>:<strand>,..."` (their own introns that are bridge junctions of the contig). No line is added,
//! removed or reordered and no coordinate moves, so intron chains are unchanged. The bridges are RELATIONS, not loci:
//! the families input (`<out>.families.gtf`, [`Outcome::in_families`]) is the GTF without them, so a bridge is never
//! a locus representative or a family node.
//!
//! EVIDENCE ([`BridgeEvidence`]): every PRIMARY record (not unmapped, secondary or supplementary; QC-fail and
//! duplicate records count) with >= 1 `N`. Its introns come from every `N` op, its end is one past its last M/D/N/=/X
//! base, and its strand is `ts == '+'` XOR reverse (absent `ts` = `+`). These are the script's reads, not
//! `--polish-tes`'s evidence, whose strand is the alignment's and which skips QC-fail records (on human testis 1,230
//! of 1.10 M spliced primaries carry `ts:A:-`). The pass-1 reader collects them, so a `--genome-wide` run reads every
//! record of each contig, as the script's `bam.fetch(contig)` did; a partial `--region` reads that region's records.
//!
//! Two deviations from the scripts, neither reachable on their held-out inputs:
//! * the scripts keyed the set of bridge junctions without the contig, so `fusion_junction` could name an intron that
//!   is a bridge junction on ANOTHER contig only. Here the set is per contig. The held-out batch merges counted 0 such
//!   coordinate coincidences (`bj_cross_contig_collisions`).
//! * a transcript without a `gene_id` is never regrouped (as in `--gtf-regroup`); the scripts grouped them under
//!   `None`. The assembler always writes a `gene_id`.
//!
//! Row order of the side tables is the scripts' single-call order (contigs by name, then gene_ids by first line, then
//! junctions by (start, end, strand)). The held-out runs merged contig batches, so their tables hold the same rows
//! with the contig blocks in batch order.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use anyhow::{Context, Result};
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::alignment::record::data::field::Value;

use crate::genome::GenomeIndex;
use crate::vg_family::denovo_assemble::{rt_real_starts, rt_v1};

/// `--bridge-regroup`'s arms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    /// Every bridge junction (`f1_bridge.py --mode full`).
    F1,
    /// The bridge junctions whose link carries fewer reads than each side it separates (`f1v2.py --rule min`).
    F1v2,
}

impl Mode {
    /// `off` = `Ok(None)`; `f1` / `f1v2`; anything else is an error.
    pub fn parse(v: &str) -> Result<Option<Mode>> {
        match v {
            "off" => Ok(None),
            "f1" => Ok(Some(Mode::F1)),
            "f1v2" => Ok(Some(Mode::F1v2)),
            other => anyhow::bail!("--bridge-regroup must be off, f1 or f1v2, got `{other}`"),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Mode::F1 => "f1",
            Mode::F1v2 => "f1v2",
        }
    }
}

// ============================================================================================ read evidence

/// One contig's read evidence, in 1-based genomic coordinates (the script's).
#[derive(Debug, Default)]
struct ContigEvidence {
    /// the U population's deduplication keys (a 128-bit hash of strand, start, end and intron chain); emptied once the
    /// contig's region has been read
    seen: HashSet<u128>,
    /// per strand (`+`, `-`): the deduplicated records' (3' end, 5' end)
    ends: [Vec<(u64, u64)>; 2],
    /// per strand: every record's (5' end, 3' end, first donor), [`rt_v1`]'s rows
    rows: [Vec<(u64, u64, u64)>; 2],
}

/// F1's read evidence per contig (see the module header for the population). Fed one record at a time by the
/// streaming pass-1 reader (a field of `Pass1Acc`) or by [`bridge_evidence_region`] on the buffered path.
#[derive(Debug, Default)]
pub struct BridgeEvidence {
    contigs: BTreeMap<String, ContigEvidence>,
    /// [`BridgeEvidence::push`]'s intron buffer, reused across records
    introns: Vec<(u64, u64)>,
}

impl BridgeEvidence {
    /// One PRIMARY alignment (the caller excludes unmapped, secondary and supplementary records): its 0-based start,
    /// its CIGAR ops and its `ts == '+'` (asked only when the record has an intron). The script's CIGAR walk: every
    /// `N` is an intron `(pos + 1, pos + len)`, and M/D/N/=/X advance the reference position.
    pub fn push(&mut self, chrom: &str, reverse: bool, ref_start: u64, ops: &[Op], ts_plus: impl FnOnce() -> bool) {
        let BridgeEvidence { contigs, introns } = self;
        introns.clear();
        let mut pos = ref_start;
        for op in ops {
            let len = op.len() as u64;
            match op.kind() {
                Kind::Skip => {
                    introns.push((pos + 1, pos + len));
                    pos += len;
                }
                Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => pos += len,
                _ => {}
            }
        }
        if introns.is_empty() {
            return;
        }
        let minus = ts_plus() == reverse;
        let ref_end = pos;
        if !contigs.contains_key(chrom) {
            contigs.insert(chrom.to_string(), ContigEvidence::default());
        }
        let c = contigs.get_mut(chrom).expect("inserted above");
        let k = usize::from(minus);
        let (e5, e3) = if minus { (ref_end, ref_start + 1) } else { (ref_start + 1, ref_end) };
        // the read's own first donor in transcript orientation (the script's `first_donor`)
        let donor = if minus { introns[introns.len() - 1].1 + 1 } else { introns[0].0 - 1 };
        c.rows[k].push((e5, e3, donor));
        if c.seen.insert(record_key(minus, ref_start, ref_end, introns)) {
            c.ends[k].push((e3, e5));
        }
    }

    /// Drop the deduplication keys and the spare capacity once a region has been read, so that a genome-wide run,
    /// which computes every region before draining any, never holds every contig's keys at once.
    pub fn seal(&mut self) {
        for c in self.contigs.values_mut() {
            c.seen = HashSet::new();
            for k in 0..2 {
                c.ends[k].shrink_to_fit();
                c.rows[k].shrink_to_fit();
            }
        }
        self.introns = Vec::new();
    }

    /// Fold one (sealed) region's evidence in. Each contig must come from ONE region (`copy_assign` checks one region
    /// per contig before any read is touched): a second region of a contig would count its boundary records twice.
    pub fn absorb(&mut self, other: BridgeEvidence) -> Result<()> {
        for (chrom, c) in other.contigs {
            anyhow::ensure!(
                !self.contigs.contains_key(&chrom),
                "--bridge-regroup: contig {chrom} was read by two regions (its V1 rows would count twice)"
            );
            self.contigs.insert(chrom, c);
        }
        Ok(())
    }

    /// Records held (spliced primaries, every contig).
    pub fn len(&self) -> usize {
        self.contigs.values().map(|c| c.rows[0].len() + c.rows[1].len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// The U population's deduplication key, the script's `(strand, reference_start, reference_end, introns)`, as a
/// 128-bit hash (two independently salted 64-bit halves).
fn record_key(minus: bool, ref_start: u64, ref_end: u64, introns: &[(u64, u64)]) -> u128 {
    use std::hash::{Hash, Hasher};
    let half = |salt: u8| {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        (salt, minus, ref_start, ref_end, introns).hash(&mut h);
        h.finish()
    };
    (u128::from(half(0)) << 64) | u128::from(half(1))
}

/// The script's `ts == '+'` for a record's FIRST `ts` field (`get_tag('ts') if has_tag('ts') else '+'`): absent = `+`;
/// `ts:A` and `ts:Z` compare by value; any other type is not `+`. A malformed field before it reads as absent, as
/// [`crate::bam::record_ts`] does.
pub fn ts_is_plus<R: noodles_sam::alignment::Record + ?Sized>(record: &R) -> bool {
    for entry in record.data().iter() {
        let Ok((tag, value)) = entry else { return true };
        if tag == crate::bam::TS_TAG {
            return match value {
                Value::Character(c) => c == b'+',
                Value::String(s) => {
                    let s: &[u8] = s;
                    s == b"+"
                }
                _ => false,
            };
        }
    }
    true
}

/// The buffered path's evidence (the `--assemble-only` runs that do not stream, e.g. `--materialize-reads`): one
/// indexed pass over `[lo, hi)` of `chrom`, feeding `ev` exactly as the streaming reader does.
pub fn bridge_evidence_region(bam_path: &str, chrom: &str, lo: u64, hi: u64, ev: &mut BridgeEvidence) -> Result<()> {
    let bai_path = format!("{bam_path}.bai");
    anyhow::ensure!(std::path::Path::new(&bai_path).exists(), "--bridge-regroup needs a .bai index");
    let file = std::fs::File::open(bam_path)?;
    let buf = std::io::BufReader::with_capacity(1 << 20, file);
    let bgzf = noodles_bgzf::MultithreadedReader::with_worker_count(std::num::NonZeroUsize::MIN, buf);
    let mut reader = noodles_bam::io::Reader::from(bgzf);
    let header = reader.read_header()?;
    let index = noodles_bam::bai::read(&bai_path)?;
    let region: noodles_core::Region = format!("{chrom}:{}-{}", lo + 1, hi).parse()?;
    let mut ops: Vec<Op> = Vec::with_capacity(256);
    for result in reader.query(&header, &index, &region)? {
        let record = result?;
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        let Some(start) = record.alignment_start() else { continue };
        let ref_start = (usize::from(start?) as u64).saturating_sub(1);
        ops.clear();
        for op in record.cigar().iter() {
            ops.push(op?);
        }
        ev.push(chrom, flags.is_reverse_complemented(), ref_start, &ops, || ts_is_plus(&record));
    }
    Ok(())
}

// ============================================================================================ the rule

/// One 3' cluster of U ends as `--polish-tes` builds it: oriented mode (5'->3' coordinate: genomic on `+`, negated on
/// `-`), ends in the cluster, PAS-proven (canonical PAS and not internally primed), not internally primed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EndCluster {
    pub mode: i64,
    pub n: usize,
    pub proven: bool,
    pub unprimed: bool,
}

/// One transcript of the GTF: its `transcript` line and `exon` lines (1-based closed, sorted, non-overlapping). The
/// scripts' and `rg3.py`'s one transcript record, shared with `--gtf-regroup`.
pub struct Tx {
    pub tid: String,
    /// its `gene_id`; a transcript without one is never regrouped
    pub gene: Option<String>,
    pub chrom: String,
    pub strand: String,
    /// the `reads` attribute (absent or non-integer = 0)
    pub reads: i64,
    /// end - start + 1 of the `transcript` line
    pub span: i64,
    pub exons: Vec<(i64, i64)>,
}

impl Tx {
    fn introns(&self) -> impl Iterator<Item = (i64, i64)> + '_ {
        self.exons.windows(2).map(|w| (w[0].1 + 1, w[1].0 - 1))
    }
}

/// The value of attribute `key` (the first `key "` in `attrs`, up to the next `"`), as the scripts' `attr`.
fn attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    let mut from = 0;
    while let Some(i) = attrs[from..].find(key) {
        let k = from + i;
        if let Some(v) = attrs[k + key.len()..].strip_prefix(" \"") {
            return v.find('"').map(|j| &v[..j]);
        }
        from = k + 1;
    }
    None
}

/// The scripts' (and `rg3.py`'s) `parse`: every `transcript` line (index = its order among them) and its `exon`
/// lines. Errors, prefixed by `flag` (the caller's flag name): a duplicate `transcript` line, an `exon` line before
/// its `transcript` line, unparsable coordinates, a transcript without an exon line, and exons that are not strictly
/// ordered and non-overlapping (the scripts assumed the last two; [`ContigProof::up`]'s slice query needs them).
pub fn parse(lines: &[String], flag: &str) -> Result<Vec<Tx>> {
    let mut txs: Vec<Tx> = Vec::new();
    let mut by_id: HashMap<String, usize> = HashMap::new();
    for line in lines {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            continue;
        }
        let Some(tid) = attr(f[8], "transcript_id") else { continue };
        if f[2] == "transcript" {
            anyhow::ensure!(!by_id.contains_key(tid), "{flag}: duplicate transcript line for {tid}");
            let start: i64 = f[3].parse().with_context(|| format!("{flag}: bad start on transcript {tid}"))?;
            let end: i64 = f[4].parse().with_context(|| format!("{flag}: bad end on transcript {tid}"))?;
            // the scripts: `int(rv) if rv.lstrip('-').isdigit() else 0`
            let reads = attr(f[8], "reads")
                .filter(|v| v.strip_prefix('-').unwrap_or(v).chars().all(|c| c.is_ascii_digit()))
                .and_then(|v| v.parse::<i64>().ok())
                .unwrap_or(0);
            by_id.insert(tid.to_string(), txs.len());
            txs.push(Tx {
                tid: tid.to_string(),
                gene: attr(f[8], "gene_id").map(str::to_string),
                chrom: f[0].to_string(),
                strand: f[6].to_string(),
                reads,
                span: end - start + 1,
                exons: Vec::new(),
            });
        } else if f[2] == "exon" {
            let Some(&i) = by_id.get(tid) else {
                anyhow::bail!("{flag}: exon line before the transcript line of {tid}");
            };
            let a: i64 = f[3].parse().with_context(|| format!("{flag}: bad exon start on {tid}"))?;
            let b: i64 = f[4].parse().with_context(|| format!("{flag}: bad exon end on {tid}"))?;
            txs[i].exons.push((a, b));
        }
    }
    for t in txs.iter_mut() {
        anyhow::ensure!(!t.exons.is_empty(), "{flag}: transcript {} has no exon line", t.tid);
        t.exons.sort_unstable();
        for &(a, b) in &t.exons {
            anyhow::ensure!(a <= b, "{flag}: transcript {} has an exon {a}-{b} whose end precedes its start", t.tid);
        }
        // strictly ordered and non-overlapping, so every intron of `introns()` has s <= e + 1 (a duplicated or
        // overlapping exon pair would give an intron with e < s - 1, and `ContigProof::up`'s ends[lo..hi] a lo > hi)
        for w in t.exons.windows(2) {
            anyhow::ensure!(
                w[1].0 > w[0].1,
                "{flag}: transcript {} has overlapping or duplicated exons {}-{} and {}-{}",
                t.tid,
                w[0].0,
                w[0].1,
                w[1].0,
                w[1].1
            );
        }
    }
    Ok(txs)
}

fn find(p: &mut [usize], mut x: usize) -> usize {
    while p[x] != x {
        p[x] = p[p[x]];
        x = p[x];
    }
    x
}

/// The scripts' `components`: the same-strand exon-overlap components of `ts` (indices into `txs`), by a per-strand
/// start-sorted sweep; components in the order of their first member, members in `ts` order. The contig is NOT part
/// of the key (the callers pass one contig's transcripts, or one gene_id's).
fn components(txs: &[Tx], ts: &[usize]) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..ts.len()).collect();
    let mut by_strand: BTreeMap<&str, Vec<(i64, i64, usize)>> = BTreeMap::new();
    for (i, &t) in ts.iter().enumerate() {
        let v = by_strand.entry(txs[t].strand.as_str()).or_default();
        v.extend(txs[t].exons.iter().map(|&(a, b)| (a, b, i)));
    }
    for ex in by_strand.values_mut() {
        ex.sort_unstable();
        let mut max_end: Option<(i64, usize)> = None;
        for &(a, b, i) in ex.iter() {
            if let Some((me, owner)) = max_end {
                if a <= me {
                    let (ra, rb) = (find(&mut parent, i), find(&mut parent, owner));
                    if ra != rb {
                        parent[ra] = rb;
                    }
                }
            }
            if max_end.map_or(true, |(me, _)| b > me) {
                max_end = Some((b, i));
            }
        }
    }
    let mut slot: HashMap<usize, usize> = HashMap::new();
    let mut out: Vec<Vec<usize>> = Vec::new();
    for (i, &t) in ts.iter().enumerate() {
        let r = find(&mut parent, i);
        let k = *slot.entry(r).or_insert_with(|| {
            out.push(Vec::new());
            out.len() - 1
        });
        out[k].push(t);
    }
    out
}

/// `--gtf-regroup`'s (RG3's) pieces of one `gene_id` and their names, the one rule both flags apply (`rg3.py`
/// ec17e540 `pieces` + `regroup`; `f1_bridge.py`'s `regroup` on the non-bridge transcripts). `ts` (the gene's
/// transcripts to split, in line order) fall into same-contig, same-strand exon-overlap components, contig by contig
/// in first-line order; a piece's representative is its max (reads, span, -index) member; the piece whose
/// representative is best keeps `gene` and the others become `<gene>.rg<k>`, k = 2.. in the order of their
/// representatives' indices. Returns the pieces and their names, parallel; a single piece keeps `gene`.
pub fn rg3_pieces(txs: &[Tx], gene: &str, ts: &[usize]) -> (Vec<Vec<usize>>, Vec<String>) {
    let key = |t: usize| (txs[t].reads, txs[t].span, std::cmp::Reverse(t));
    let mut chroms: Vec<&str> = Vec::new();
    let mut by_chrom: HashMap<&str, Vec<usize>> = HashMap::new();
    for &t in ts {
        by_chrom
            .entry(txs[t].chrom.as_str())
            .or_insert_with(|| {
                chroms.push(txs[t].chrom.as_str());
                Vec::new()
            })
            .push(t);
    }
    let comps: Vec<Vec<usize>> = chroms.iter().flat_map(|c| components(txs, &by_chrom[c])).collect();
    let reps: Vec<usize> =
        comps.iter().map(|c| *c.iter().max_by_key(|&&t| key(t)).expect("a component has a member")).collect();
    let mut name: Vec<String> = vec![String::new(); comps.len()];
    if let Some(keep) = (0..comps.len()).max_by_key(|&i| key(reps[i])) {
        let mut order: Vec<usize> = (0..comps.len()).collect();
        order.sort_unstable_by_key(|&i| reps[i]);
        let mut k = 2;
        for i in order {
            if i == keep {
                name[i] = gene.to_string();
            } else {
                name[i] = format!("{gene}.rg{k}");
                k += 1;
            }
        }
    }
    (comps, name)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Side {
    Up,
    Down,
    Straddle,
    Inside,
}

/// The scripts' `side`: where a component lies against the intron [s, e], in transcript orientation.
fn side(txs: &[Tx], comp: &[usize], s: i64, e: i64, strand: &str) -> Side {
    let lo = comp.iter().map(|&t| txs[t].exons[0].0).min().expect("a component has a member");
    let hi = comp.iter().map(|&t| txs[t].exons[txs[t].exons.len() - 1].1).max().expect("a component has a member");
    match (lo < s, hi > e) {
        (true, true) => Side::Straddle,
        (false, false) => Side::Inside,
        (before, _) if strand == "+" => {
            if before {
                Side::Up
            } else {
                Side::Down
            }
        }
        (before, _) => {
            if before {
                Side::Down
            } else {
                Side::Up
            }
        }
    }
}

/// One STRUCTURAL junction of one gene_id (a row of the junctions table).
struct Junction {
    gene: String,
    chrom: String,
    s: i64,
    e: i64,
    strand: String,
    /// T_J, one entry per transcript line using J (the script's list)
    tj: Vec<usize>,
    /// the transcripts of the UP / DOWN components (F1v2's sides)
    up_tx: Vec<usize>,
    down_tx: Vec<usize>,
    n_up: usize,
    n_down: usize,
    n_inside: usize,
    u: usize,
    /// the U population's 3' clusters: (genomic mode, ends, proven, unprimed), in 5'->3' order
    clusters: Vec<(i64, usize, bool, bool)>,
    up_proof: bool,
    v1: u64,
    down_proof: bool,
    /// F1's BRIDGE(J)
    bridge: bool,
}

/// One contig's evidence sorted for queries, with its start clusters and sequence.
struct ContigProof {
    ev: ContigEvidence,
    real: [Vec<bool>; 2],
    genome: Arc<GenomeIndex>,
}

impl ContigProof {
    fn new(mut ev: ContigEvidence, genome: Arc<GenomeIndex>) -> Self {
        for k in 0..2 {
            ev.ends[k].sort_unstable();
            ev.rows[k].sort_unstable();
        }
        let real = [rt_real_starts(&ev.rows[0]), rt_real_starts(&ev.rows[1])];
        ContigProof { ev, real, genome }
    }

    /// The script's `Evidence.up`: the U population of [s, e] (its size) and its 3' clusters.
    fn up(
        &self,
        chrom: &str,
        s: i64,
        e: i64,
        strand: &str,
        end_clusters: &dyn Fn(&[i64], &[u8], bool) -> Vec<EndCluster>,
    ) -> (usize, Vec<(i64, usize, bool, bool)>) {
        let (k, minus) = match strand {
            "+" => (0, false),
            "-" => (1, true),
            _ => return (0, Vec::new()),
        };
        let ends = &self.ev.ends[k];
        let (lo, hi) = (ends.partition_point(|x| (x.0 as i64) < s), ends.partition_point(|x| (x.0 as i64) <= e));
        let mut up: Vec<i64> = ends[lo..hi]
            .iter()
            .filter(|&&(_, e5)| if minus { (e5 as i64) > e } else { (e5 as i64) < s })
            .map(|&(e3, _)| if minus { -(e3 as i64) } else { e3 as i64 })
            .collect();
        up.sort_unstable();
        let seq: &[u8] = self.genome.chroms().find(|(n, _)| *n == chrom).map(|(_, s)| s).unwrap_or(&[]);
        let cl = end_clusters(&up, seq, minus)
            .into_iter()
            .map(|c| (if minus { -c.mode } else { c.mode }, c.n, c.proven, c.unprimed))
            .collect();
        (up.len(), cl)
    }

    /// The script's `Evidence.v1_count` of [s, e]: [`rt_v1`] on the 1-based rows, queried with `[s, e + 1)`.
    fn v1(&self, s: i64, e: i64, strand: &str) -> u64 {
        let (k, minus) = match strand {
            "+" => (0, false),
            "-" => (1, true),
            _ => return 0,
        };
        let q = |x: i64| u64::try_from(x).unwrap_or(0);
        rt_v1(&self.ev.rows[k], &self.real[k], q(s), q(e + 1), minus)
    }
}

/// The scripts' `decide`: every STRUCTURAL junction with its proofs, contigs by name, gene_ids by first line.
fn decide(
    txs: &[Tx],
    evidence: &mut BridgeEvidence,
    genome_of: &dyn Fn(&str) -> Result<Arc<GenomeIndex>>,
    end_clusters: &dyn Fn(&[i64], &[u8], bool) -> Vec<EndCluster>,
) -> Result<Vec<Junction>> {
    let mut order: Vec<usize> = (0..txs.len()).filter(|&i| txs[i].gene.is_some()).collect();
    order.sort_by(|&a, &b| txs[a].chrom.cmp(&txs[b].chrom)); // stable: line order within a contig
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut slot: HashMap<(&str, &str), usize> = HashMap::new();
    for &i in &order {
        let key = (txs[i].gene.as_deref().expect("filtered above"), txs[i].chrom.as_str());
        let k = *slot.entry(key).or_insert_with(|| {
            groups.push(Vec::new());
            groups.len() - 1
        });
        groups[k].push(i);
    }
    let mut rows: Vec<Junction> = Vec::new();
    // one contig's evidence at a time (groups of one contig are consecutive), built at its first STRUCTURAL junction
    let mut proof: Option<(String, ContigProof)> = None;
    for ts in &groups {
        if ts.len() < 3 {
            continue;
        }
        let (gene, chrom) = (txs[ts[0]].gene.as_deref().expect("filtered above"), txs[ts[0]].chrom.as_str());
        let mut by_junction: BTreeMap<((i64, i64), &str), Vec<usize>> = BTreeMap::new();
        for &t in ts {
            for iv in txs[t].introns() {
                by_junction.entry((iv, txs[t].strand.as_str())).or_default().push(t);
            }
        }
        for (((s, e), strand), tj) in by_junction {
            let tj_set: HashSet<usize> = tj.iter().copied().collect();
            let r: Vec<usize> = ts.iter().copied().filter(|t| txs[*t].strand == strand && !tj_set.contains(t)).collect();
            if r.len() < 2 {
                continue;
            }
            let comps = components(txs, &r);
            let sides: Vec<Side> = comps.iter().map(|c| side(txs, c, s, e, strand)).collect();
            let count = |x: Side| sides.iter().filter(|&&y| y == x).count();
            let (n_up, n_down) = (count(Side::Up), count(Side::Down));
            if count(Side::Straddle) > 0 || n_up == 0 || n_down == 0 {
                continue;
            }
            if proof.as_ref().map_or(true, |(c, _)| c != chrom) {
                drop(proof.take()); // the previous contig's evidence is freed before the next one is sorted
                let ev = evidence.contigs.remove(chrom).unwrap_or_default();
                proof = Some((chrom.to_string(), ContigProof::new(ev, genome_of(chrom)?)));
            }
            let p = &proof.as_ref().expect("built above").1;
            let (u, clusters) = p.up(chrom, s, e, strand, end_clusters);
            let up_proof = clusters.iter().any(|c| c.2);
            let v1 = p.v1(s, e, strand);
            let down_proof = v1 >= 1;
            let members = |x: Side| -> Vec<usize> {
                comps.iter().zip(&sides).filter(|(_, &y)| y == x).flat_map(|(c, _)| c.iter().copied()).collect()
            };
            rows.push(Junction {
                gene: gene.to_string(),
                chrom: chrom.to_string(),
                s,
                e,
                strand: strand.to_string(),
                tj,
                up_tx: members(Side::Up),
                down_tx: members(Side::Down),
                n_up,
                n_down,
                n_inside: count(Side::Inside),
                u,
                clusters,
                up_proof,
                v1,
                down_proof,
                bridge: up_proof && down_proof,
            });
        }
    }
    Ok(rows)
}

/// F1v2's MINORITY(J) and its table row (`f1v2.py`'s `decide_v2`). T_J is the set of transcripts using J.
fn minority_row(txs: &[Tx], j: &Junction) -> (bool, String) {
    let mut tj: Vec<usize> = j.tj.clone();
    tj.sort_unstable();
    tj.dedup();
    let sum = |v: &[usize]| v.iter().map(|&t| txs[t].reads).sum::<i64>();
    let max = |v: &[usize]| v.iter().map(|&t| txs[t].reads).max().expect("a side has a transcript");
    let (rb, ru, rd) = (sum(&tj), sum(&j.up_tx), sum(&j.down_tx));
    let den = rb + ru.min(rd);
    let share = if den > 0 { rb as f64 / den as f64 } else { 1.0 };
    let keep = rb < ru && rb < rd;
    let row = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{rb}\t{}\t{}\t{ru}\t{}\t{}\t{rd}\t{}\t{}\t{}\t{}",
        j.gene,
        j.chrom,
        j.s,
        j.e,
        j.strand,
        tj.len(),
        max(&tj),
        j.up_tx.len(),
        max(&j.up_tx),
        j.down_tx.len(),
        max(&j.down_tx),
        py_round4(share),
        py_bool(keep),
        tid_list(txs, &tj),
    );
    (keep, row)
}

/// Python's `str(round(x, 4))` for a float in [0, 1]: `{:.4}` rounds the exact binary value half to even as `round`
/// does, and the shortest round-trip form of the result is `repr`'s (checked against Python on every rb / (rb + m),
/// rb < 400, m < 2500).
fn py_round4(x: f64) -> String {
    let r: f64 = format!("{x:.4}").parse().expect("a formatted float parses");
    let s = format!("{r}");
    if s.contains('.') {
        s
    } else {
        format!("{s}.0")
    }
}

fn py_bool(b: bool) -> &'static str {
    if b {
        "True"
    } else {
        "False"
    }
}

/// T_J's transcript ids, sorted and comma-joined (the tables' `TJ`).
fn tid_list(txs: &[Tx], ts: &[usize]) -> String {
    let mut v: Vec<&str> = ts.iter().map(|&t| txs[t].tid.as_str()).collect();
    v.sort_unstable();
    v.dedup();
    v.join(",")
}

/// The scripts' `regroup`: every transcript's new gene_id (`None` without a gene_id) and each bridge transcript's
/// (`fusion_of`, `fusion_junction`). `bj` = the bridge junctions as (contig, s, e, strand).
#[allow(clippy::type_complexity)]
fn regroup(
    txs: &[Tx],
    bridges: &HashSet<usize>,
    bj: &HashSet<(&str, i64, i64, &str)>,
) -> (Vec<Option<String>>, HashMap<usize, (String, String)>) {
    let mut new: Vec<Option<String>> = vec![None; txs.len()];
    let mut rel: HashMap<usize, (String, String)> = HashMap::new();
    let mut genes: Vec<&str> = Vec::new();
    let mut members: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, t) in txs.iter().enumerate() {
        if let Some(g) = t.gene.as_deref() {
            members
                .entry(g)
                .or_insert_with(|| {
                    genes.push(g);
                    Vec::new()
                })
                .push(i);
        }
    }
    for g in genes {
        let ts = &members[g];
        // the pieces: RG3's exon-overlap components of the non-bridge transcripts, named as RG3 names them
        let non_bridge: Vec<usize> = ts.iter().copied().filter(|t| !bridges.contains(t)).collect();
        let (comps, name) = rg3_pieces(txs, g, &non_bridge);
        for (i, c) in comps.iter().enumerate() {
            for &t in c {
                new[t] = Some(name[i].clone());
            }
        }
        // the bridges: exon-overlap groups among themselves, `<g>.fus<k>` in line order, with their relation
        let bt: Vec<usize> = ts.iter().copied().filter(|t| bridges.contains(t)).collect();
        if bt.is_empty() {
            continue;
        }
        let mut bcomps = components(txs, &bt);
        bcomps.sort_by_key(|c| *c.iter().min().expect("a component has a member"));
        for (k, c) in bcomps.iter().enumerate() {
            let bname = format!("{g}.fus{}", k + 1);
            for &t in c {
                new[t] = Some(bname.clone());
                let tx = &txs[t];
                let mut pieces: Vec<(i64, &str)> = Vec::new();
                for (i, pc) in comps.iter().enumerate() {
                    let touches = pc.iter().any(|&x| {
                        txs[x].strand == tx.strand
                            && txs[x].exons.iter().any(|&(a, b)| tx.exons.iter().any(|&(c0, d)| a <= d && c0 <= b))
                    });
                    if touches {
                        let lo = pc.iter().map(|&x| txs[x].exons[0].0).min().expect("a component has a member");
                        pieces.push((lo, name[i].as_str()));
                    }
                }
                pieces.sort_unstable();
                if tx.strand == "-" {
                    pieces.reverse();
                }
                let fusion_of: Vec<&str> = pieces.iter().map(|p| p.1).collect();
                let junctions: Vec<String> = tx
                    .introns()
                    .filter(|&(s, e)| bj.contains(&(tx.chrom.as_str(), s, e, tx.strand.as_str())))
                    .map(|(s, e)| format!("{}:{s}-{e}:{}", tx.chrom, tx.strand))
                    .collect();
                rel.insert(t, (fusion_of.join(","), junctions.join(",")));
            }
        }
    }
    (new, rel)
}

/// What one `--bridge-regroup` pass decided (the log line and `params.tsv`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Stats {
    pub transcripts: usize,
    pub gene_ids: usize,
    /// junctions passing STRUCTURAL
    pub structural_junctions: usize,
    pub up_proof: usize,
    pub down_proof: usize,
    /// F1's bridge junctions
    pub f1_bridge_junctions: usize,
    /// the bridge junctions acted on (F1: all of F1's; F1v2: those passing MINORITY)
    pub bridge_junctions: usize,
    pub bridge_transcripts: usize,
    pub gene_ids_with_bridge: usize,
    /// input gene_ids with a bridge or split into >= 2 pieces
    pub gene_ids_split: usize,
    pub gene_ids_after: usize,
    /// gene_ids of the families input (bridges excluded)
    pub families_gene_ids: usize,
    /// lines whose gene_id changed
    pub lines_changed: usize,
    /// lines of bridge transcripts (absent from the families input)
    pub family_lines_dropped: usize,
}

/// The pass's products besides the rewritten lines.
pub struct Outcome {
    pub stats: Stats,
    /// transcript ids of the bridge transcripts
    pub bridges: HashSet<String>,
    /// every STRUCTURAL junction with its evidence (`f1_bridge.py`'s `junctions.tsv`)
    pub junctions_tsv: String,
    /// every F1 bridge junction with its reads and MINORITY decision (`f1v2.py`'s `bridges.tsv`); F1v2 only
    pub bridges_tsv: Option<String>,
}

impl Outcome {
    /// Is `line` part of the families input? Every line except those of a bridge transcript.
    pub fn in_families(&self, line: &str) -> bool {
        if line.is_empty() || line.starts_with('#') {
            return true;
        }
        let f: Vec<&str> = line.split('\t').collect();
        !(f.len() >= 9 && attr(f[8], "transcript_id").is_some_and(|t| self.bridges.contains(t)))
    }
}

const JUNCTIONS_HEADER: &str =
    "gene\tchrom\ts\te\tstrand\tn_TJ\treads_TJ\tup\tdown\tinside\tU\tclusters\tup_proof\tV1\tdown_proof\tbridge\tTJ";
const BRIDGES_HEADER: &str = "gene\tchrom\ts\te\tstrand\tn_TJ\treads_TJ\tmax_TJ\tn_up\treads_up\tmax_up\tn_down\treads_down\tmax_down\tshare\tkeep\tTJ";

/// Run `--bridge-regroup` on the final lines of the emitted GTF: rewrite the `gene_id`s in place, append the relation
/// attributes to the bridge transcripts' `transcript` lines, and return the tables. `evidence` is consumed contig by
/// contig; `genome_of` gives a contig's sequence (uppercase); `end_clusters` is `--polish-tes`'s 3' cluster rule on
/// sorted oriented ends.
pub fn run(
    lines: &mut [String],
    mode: Mode,
    evidence: &mut BridgeEvidence,
    genome_of: &dyn Fn(&str) -> Result<Arc<GenomeIndex>>,
    end_clusters: &dyn Fn(&[i64], &[u8], bool) -> Vec<EndCluster>,
) -> Result<Outcome> {
    let txs = parse(lines, "--bridge-regroup")?;
    let junctions = decide(&txs, evidence, genome_of, end_clusters)?;
    let mut junctions_tsv = String::from(JUNCTIONS_HEADER);
    for j in &junctions {
        let clusters: Vec<String> = j
            .clusters
            .iter()
            .map(|&(m, n, p, u)| format!("{m}:{n}:{}", if p { "P" } else if u { "u" } else { "-" }))
            .collect();
        junctions_tsv.push_str(&format!(
            "\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            j.gene,
            j.chrom,
            j.s,
            j.e,
            j.strand,
            j.tj.len(),
            j.tj.iter().map(|&t| txs[t].reads).sum::<i64>(),
            j.n_up,
            j.n_down,
            j.n_inside,
            j.u,
            if clusters.is_empty() { ".".to_string() } else { clusters.join(";") },
            py_bool(j.up_proof),
            j.v1,
            py_bool(j.down_proof),
            py_bool(j.bridge),
            tid_list(&txs, &j.tj),
        ));
    }
    junctions_tsv.push('\n');
    // the junctions acted on: F1's bridges, or those of them passing MINORITY
    let mut kept: Vec<bool> = junctions.iter().map(|j| j.bridge).collect();
    let bridges_tsv = if mode == Mode::F1v2 {
        let rows: Vec<(usize, bool, String)> = junctions
            .iter()
            .enumerate()
            .filter(|(_, j)| j.bridge)
            .map(|(i, j)| {
                let (keep, row) = minority_row(&txs, j);
                (i, keep, row)
            })
            .collect();
        // f1v2.py takes its header from its first row: a table without rows is the lone column `gene`
        let mut out = String::from(if rows.is_empty() { "gene" } else { BRIDGES_HEADER });
        for (i, keep, row) in rows {
            kept[i] = keep;
            out.push('\n');
            out.push_str(&row);
        }
        out.push('\n');
        Some(out)
    } else {
        None
    };
    let mut bridges: HashSet<usize> = HashSet::new();
    let mut bj: HashSet<(&str, i64, i64, &str)> = HashSet::new();
    for (j, _) in junctions.iter().zip(&kept).filter(|(_, &k)| k) {
        bridges.extend(j.tj.iter().copied());
        bj.insert((j.chrom.as_str(), j.s, j.e, j.strand.as_str()));
    }
    let (new, rel) = regroup(&txs, &bridges, &bj);

    // the scripts' `rewrite`: gene_id "<old>" -> gene_id "<new>" (first occurrence) on every line of a relabelled
    // transcript, then the relation appended to a bridge's `transcript` line
    let id_of: HashMap<&str, usize> = txs.iter().enumerate().map(|(i, t)| (t.tid.as_str(), i)).collect();
    let (mut lines_changed, mut family_lines_dropped) = (0usize, 0usize);
    for line in lines.iter_mut() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            continue;
        }
        let Some(&i) = attr(f[8], "transcript_id").and_then(|t| id_of.get(t)) else { continue };
        let mut attrs: Option<String> = None;
        if let (Some(nv), Some(old)) = (new[i].as_deref(), attr(f[8], "gene_id")) {
            if nv != old {
                attrs = Some(f[8].replacen(&format!("gene_id \"{old}\""), &format!("gene_id \"{nv}\""), 1));
                lines_changed += 1;
            }
        }
        if let Some((fusion_of, junction)) = rel.get(&i) {
            family_lines_dropped += 1;
            if f[2] == "transcript" {
                let base = attrs.take().unwrap_or_else(|| f[8].to_string());
                attrs =
                    Some(format!("{} fusion_of \"{fusion_of}\"; fusion_junction \"{junction}\";", base.trim_end()));
            }
        }
        if let Some(a) = attrs {
            let mut out = f.clone();
            out[8] = &a;
            *line = out.join("\t");
        }
    }

    let bgenes: HashSet<&str> = bridges.iter().filter_map(|&t| txs[t].gene.as_deref()).collect();
    let mut split: HashSet<&str> = txs
        .iter()
        .enumerate()
        .filter(|(i, t)| !bridges.contains(i) && t.gene.is_some() && new[*i].as_deref() != t.gene.as_deref())
        .filter_map(|(_, t)| t.gene.as_deref())
        .collect();
    split.extend(bgenes.iter().copied());
    let stats = Stats {
        transcripts: txs.len(),
        gene_ids: txs.iter().filter_map(|t| t.gene.as_deref()).collect::<HashSet<_>>().len(),
        structural_junctions: junctions.len(),
        up_proof: junctions.iter().filter(|j| j.up_proof).count(),
        down_proof: junctions.iter().filter(|j| j.down_proof).count(),
        f1_bridge_junctions: junctions.iter().filter(|j| j.bridge).count(),
        bridge_junctions: kept.iter().filter(|&&k| k).count(),
        bridge_transcripts: bridges.len(),
        gene_ids_with_bridge: bgenes.len(),
        gene_ids_split: split.len(),
        gene_ids_after: new.iter().flatten().collect::<HashSet<_>>().len(),
        families_gene_ids: new
            .iter()
            .enumerate()
            .filter(|(i, _)| !bridges.contains(i))
            .filter_map(|(_, n)| n.as_deref())
            .collect::<HashSet<_>>()
            .len(),
        lines_changed,
        family_lines_dropped,
    };
    Ok(Outcome {
        stats,
        bridges: bridges.iter().map(|&t| txs[t].tid.clone()).collect(),
        junctions_tsv,
        bridges_tsv,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A transcript's GTF lines (`transcript` + `exon`s; 1-based closed exons).
    fn tx(chrom: &str, gene: &str, tid: &str, reads: i64, exons: &[(i64, i64)], strand: &str) -> Vec<String> {
        let a = format!("gene_id \"{gene}\"; transcript_id \"{tid}\"; reads \"{reads}\"; TPM \"1.0\";");
        let (lo, hi) = (exons[0].0, exons[exons.len() - 1].1);
        let mut v = vec![format!("{chrom}\trustle\ttranscript\t{lo}\t{hi}\t.\t{strand}\t.\t{a}")];
        for (k, &(s, e)) in exons.iter().enumerate() {
            v.push(format!(
                "{chrom}\trustle\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"{gene}\"; transcript_id \"{tid}\"; exon_number \"{}\";",
                k + 1
            ));
        }
        v
    }

    /// An alignment over 1-based closed `exons` as (0-based start, M/N ops).
    fn aln(exons: &[(u64, u64)]) -> (u64, Vec<Op>) {
        let mut ops = vec![Op::new(Kind::Match, (exons[0].1 - exons[0].0 + 1) as usize)];
        for w in exons.windows(2) {
            ops.push(Op::new(Kind::Skip, (w[1].0 - w[0].1 - 1) as usize));
            ops.push(Op::new(Kind::Match, (w[1].1 - w[1].0 + 1) as usize));
        }
        (exons[0].0 - 1, ops)
    }

    fn push(ev: &mut BridgeEvidence, exons: &[(u64, u64)], minus: bool) {
        let (start, ops) = aln(exons);
        // a reverse alignment with `ts:A:+` is a `-` transcript
        ev.push("c1", minus, start, &ops, || true);
    }

    /// The test's 3' cluster rule: every end in one cluster when >= 2 (mode = the 3'-most), PAS-proven iff `pas`.
    fn fake_clusters(pas: bool) -> impl Fn(&[i64], &[u8], bool) -> Vec<EndCluster> {
        move |ends: &[i64], _: &[u8], _: bool| {
            if ends.len() < 2 {
                return Vec::new();
            }
            vec![EndCluster { mode: *ends.last().unwrap(), n: ends.len(), proven: pas, unprimed: true }]
        }
    }

    fn genome(_: &str) -> Result<Arc<GenomeIndex>> {
        Ok(Arc::new(GenomeIndex::from_seqs(&[("c1", &[b'C'; 3000])])))
    }

    /// The `+` locus: X (reads `rx`) ends inside J = [351, 1300]; the bridge B (reads `rb`) joins X's first exon to
    /// Y's last two exons through J; Y (reads 8) starts inside J at its own promoter.
    fn plus_locus(rx: i64, rb: i64) -> Vec<String> {
        let mut v = tx("c1", "G", "X", rx, &[(101, 200), (301, 400)], "+");
        v.extend(tx("c1", "G", "B", rb, &[(101, 200), (301, 350), (1301, 1400), (1501, 1600)], "+"));
        v.extend(tx("c1", "G", "Y", 8, &[(1101, 1200), (1301, 1400), (1501, 1600)], "+"));
        v
    }

    /// X's reads (3' ends 398-400 inside J, deduplicated to three) and Y's (three at one start: a real cluster, first
    /// exon inside J), mirrored onto the `-` strand when `minus` (x -> 2001 - x).
    fn locus_evidence(minus: bool, with_y: bool) -> BridgeEvidence {
        let m = |ex: &[(u64, u64)]| -> Vec<(u64, u64)> {
            if minus {
                ex.iter().rev().map(|&(a, b)| (2001 - b, 2001 - a)).collect()
            } else {
                ex.to_vec()
            }
        };
        let mut ev = BridgeEvidence::default();
        for end in [398, 399, 400, 400] {
            push(&mut ev, &m(&[(101, 200), (301, end)]), minus);
        }
        if with_y {
            for _ in 0..3 {
                push(&mut ev, &m(&[(1101, 1200), (1301, 1400), (1501, 1600)]), minus);
            }
        }
        ev
    }

    fn run_on(lines: &mut [String], mode: Mode, mut ev: BridgeEvidence, pas: bool) -> Outcome {
        run(lines, mode, &mut ev, &genome, &fake_clusters(pas)).expect("run")
    }

    fn gene_of(lines: &[String], tid: &str) -> String {
        let l = lines.iter().find(|l| l.contains("\ttranscript\t") && attr(l, "transcript_id") == Some(tid)).unwrap();
        attr(l, "gene_id").unwrap().to_string()
    }

    #[test]
    fn evidence_takes_every_n_and_ends_one_past_the_last_reference_base() {
        let mut ev = BridgeEvidence::default();
        let (m, d, n, i, s) = (Kind::Match, Kind::Deletion, Kind::Skip, Kind::Insertion, Kind::SoftClip);
        // 10M 5D 100N 3I 20M 2S from 1-based 1000: intron [1015, 1114], last base 1134, first donor 1014
        let ops = [Op::new(m, 10), Op::new(d, 5), Op::new(n, 100), Op::new(i, 3), Op::new(m, 20), Op::new(s, 2)];
        ev.push("c1", false, 999, &ops, || true);
        assert_eq!(ev.contigs["c1"].rows[0], vec![(1000, 1134, 1014)]);
        assert_eq!(ev.contigs["c1"].ends[0], vec![(1134, 1000)]);
        // a leading N and N-I-N are introns of their own (the script's walk, not the assembler's exon blocks)
        let ops = [Op::new(n, 50), Op::new(m, 10), Op::new(n, 20), Op::new(i, 2), Op::new(n, 30), Op::new(m, 10)];
        ev.push("c2", false, 0, &ops, || true);
        assert_eq!(ev.contigs["c2"].rows[0], vec![(1, 120, 0)]);
        // `=` / `X` advance the reference; an unspliced record is never kept and its `ts` never read
        ev.push("c3", false, 0, &[Op::new(Kind::SequenceMatch, 5), Op::new(Kind::SequenceMismatch, 1)], || panic!("ts read"));
        assert!(!ev.contigs.contains_key("c3"));
        assert_eq!(ev.len(), 2);
    }

    #[test]
    fn strand_is_ts_xor_reverse() {
        let mut ev = BridgeEvidence::default();
        let (start, ops) = aln(&[(101, 200), (301, 400)]);
        ev.push("c1", false, start, &ops, || true); // forward, ts + : +
        ev.push("c1", true, start, &ops, || false); // reverse, ts - : +
        ev.push("c1", true, start, &ops, || true); // reverse, ts + : -
        ev.push("c1", false, start, &ops, || false); // forward, ts - : -
        let c = &ev.contigs["c1"];
        assert_eq!(c.rows[0], vec![(101, 400, 200), (101, 400, 200)]);
        // `-`: 5' end = the last base, 3' end = the first, first donor = the base after the last intron
        assert_eq!(c.rows[1], vec![(400, 101, 301), (400, 101, 301)]);
        assert_eq!((c.ends[0].len(), c.ends[1].len()), (1, 1), "one U record per strand after deduplication");
    }

    #[test]
    fn u_ends_deduplicate_on_strand_start_end_and_chain_but_v1_rows_do_not() {
        let mut ev = BridgeEvidence::default();
        push(&mut ev, &[(101, 200), (301, 400)], false);
        push(&mut ev, &[(101, 200), (301, 400)], false);
        push(&mut ev, &[(101, 200), (301, 401)], false); // another end
        push(&mut ev, &[(101, 200), (311, 400)], false); // another chain
        let c = &ev.contigs["c1"];
        assert_eq!(c.rows[0].len(), 4);
        assert_eq!(c.ends[0], vec![(400, 101), (401, 101), (400, 101)]);
    }

    #[test]
    fn seal_drops_the_keys_and_absorb_refuses_a_contig_read_twice() {
        let (mut all, mut region, mut other) = (BridgeEvidence::default(), locus_evidence(false, true), BridgeEvidence::default());
        region.seal();
        assert!(region.contigs["c1"].seen.is_empty());
        assert_eq!((region.contigs["c1"].ends[0].len(), region.len()), (4, 7), "the evidence itself is kept");
        all.absorb(region).unwrap();
        push(&mut other, &[(101, 200), (301, 400)], false);
        assert!(all.absorb(other).is_err());
    }

    #[test]
    fn ts_is_plus_reads_the_first_ts_field() {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value as BufValue;
        let rec = |v: Option<BufValue>| {
            let data: noodles_sam::alignment::record_buf::Data =
                v.map(|v| (Tag::new(b't', b's'), v)).into_iter().collect();
            noodles_sam::alignment::RecordBuf::builder().set_data(data).build()
        };
        assert!(ts_is_plus(&rec(None)), "absent = +");
        assert!(ts_is_plus(&rec(Some(BufValue::Character(b'+')))));
        assert!(!ts_is_plus(&rec(Some(BufValue::Character(b'-')))));
        assert!(ts_is_plus(&rec(Some(BufValue::String("+".into())))), "ts:Z:+ compares by value");
        assert!(!ts_is_plus(&rec(Some(BufValue::Int32(43)))), "another type is not '+'");
    }

    #[test]
    fn attr_is_the_first_key_quote_match() {
        let a = "gene_id \"G\"; matched_reads \"0\"; reads \"7\";";
        assert_eq!(attr(a, "reads"), Some("0"), "as the scripts' attr: the first `reads \"` wins");
        assert_eq!(attr("reads \"7\"; matched_reads \"0\";", "reads"), Some("7"));
        assert_eq!(attr("transcript_idx \"a\"; transcript_id \"b\";", "transcript_id"), Some("b"));
        assert_eq!(attr("gene_id \"G;", "gene_id"), None);
    }

    #[test]
    fn sides_are_in_transcript_orientation() {
        let mut lines = tx("c1", "G", "a", 1, &[(100, 200)], "+");
        lines.extend(tx("c1", "G", "b", 1, &[(350, 500)], "+"));
        lines.extend(tx("c1", "G", "c", 1, &[(100, 120), (450, 500)], "+"));
        lines.extend(tx("c1", "G", "d", 1, &[(320, 380)], "+"));
        let txs = parse(&lines, "--bridge-regroup").unwrap();
        let s = |t: usize, strand: &str| side(&txs, &[t], 300, 400, strand);
        assert_eq!((s(0, "+"), s(1, "+"), s(2, "+"), s(3, "+")), (Side::Up, Side::Down, Side::Straddle, Side::Inside));
        assert_eq!((s(0, "-"), s(1, "-"), s(2, "-"), s(3, "-")), (Side::Down, Side::Up, Side::Straddle, Side::Inside));
    }

    #[test]
    fn components_join_same_strand_exon_overlap_only() {
        let mut lines = tx("c1", "G", "a", 1, &[(100, 200), (500, 600)], "+");
        lines.extend(tx("c1", "G", "b", 1, &[(300, 400)], "+")); // inside a's intron: not overlap
        lines.extend(tx("c1", "G", "c", 1, &[(600, 700)], "+")); // one shared base with a
        lines.extend(tx("c1", "G", "d", 1, &[(650, 800)], "-")); // other strand
        lines.extend(tx("c1", "G", "e", 1, &[(790, 900)], "+")); // chains through nothing on its strand
        let txs = parse(&lines, "--bridge-regroup").unwrap();
        assert_eq!(components(&txs, &[0, 1, 2, 3, 4]), vec![vec![0, 2], vec![1], vec![3], vec![4]]);
    }

    /// RG3's naming on three pieces over two contigs: the best representative (reads, then span, then the earliest
    /// line) keeps the name, the others are `.rg2`, `.rg3` in representative-index order; one piece keeps the name.
    #[test]
    fn rg3_pieces_name_by_representative() {
        let mut lines = tx("c1", "G", "p", 1, &[(1, 10), (20, 30)], "+");
        lines.extend(tx("c1", "G", "q", 7, &[(100, 110), (120, 130)], "+"));
        lines.extend(tx("c2", "G", "r", 7, &[(100, 110), (120, 131)], "+")); // same reads, longer: the best
        lines.extend(tx("c1", "G", "m", 1, &[(25, 28)], "+")); // overlaps p, shorter: p's piece, p its representative
        let txs = parse(&lines, "--bridge-regroup").unwrap();
        let (comps, names) = rg3_pieces(&txs, "G", &[0, 1, 2, 3]);
        assert_eq!(comps, vec![vec![0, 3], vec![1], vec![2]]);
        assert_eq!(names, vec!["G.rg2".to_string(), "G.rg3".to_string(), "G".to_string()]);
        assert_eq!(rg3_pieces(&txs, "G", &[0, 3]), (vec![vec![0, 3]], vec!["G".to_string()]));
        assert_eq!(rg3_pieces(&txs, "G", &[]), (Vec::new(), Vec::new()));
    }

    /// The bridge: X keeps the name (10 reads > 8), Y becomes `.rg2`, B `.fus1` with its relation; the junction table
    /// holds the one STRUCTURAL junction; nothing but gene_id and the two appended attributes changes.
    #[test]
    fn f1_splits_a_proven_bridge_and_keeps_it_as_a_relation() {
        let src = plus_locus(10, 2);
        let mut lines = src.clone();
        let out = run_on(&mut lines, Mode::F1, locus_evidence(false, true), true);
        assert_eq!((gene_of(&lines, "X"), gene_of(&lines, "Y"), gene_of(&lines, "B")), ("G".into(), "G.rg2".into(), "G.fus1".into()));
        assert_eq!(
            lines[3],
            "c1\trustle\ttranscript\t101\t1600\t.\t+\t.\tgene_id \"G.fus1\"; transcript_id \"B\"; reads \"2\"; TPM \"1.0\"; \
             fusion_of \"G,G.rg2\"; fusion_junction \"c1:351-1300:+\";"
        );
        assert_eq!(lines.len(), src.len());
        for (a, b) in src.iter().zip(&lines) {
            assert_eq!(a.split('\t').take(8).collect::<Vec<_>>(), b.split('\t').take(8).collect::<Vec<_>>(), "coordinates never move");
        }
        assert_eq!(
            out.junctions_tsv,
            format!("{JUNCTIONS_HEADER}\nG\tc1\t351\t1300\t+\t1\t2\t1\t1\t0\t3\t400:3:P\tTrue\t3\tTrue\tTrue\tB\n")
        );
        assert!(out.bridges_tsv.is_none());
        let fam: Vec<&String> = lines.iter().filter(|l| out.in_families(l)).collect();
        assert_eq!(fam.len(), src.len() - 5, "the bridge's transcript and 4 exon lines leave the families input");
        assert!(fam.iter().all(|l| attr(l, "transcript_id") != Some("B")));
        let st = &out.stats;
        assert_eq!((st.structural_junctions, st.f1_bridge_junctions, st.bridge_junctions, st.bridge_transcripts), (1, 1, 1, 1));
        assert_eq!((st.gene_ids, st.gene_ids_split, st.gene_ids_after, st.families_gene_ids), (1, 1, 3, 2));
        assert_eq!((st.lines_changed, st.family_lines_dropped), (9, 5), "Y: 4 lines relabelled, B: 5");
    }

    /// Without either proof there is no bridge, and without a bridge the output is the input (one exon-overlap piece).
    #[test]
    fn a_bridge_needs_both_proofs() {
        for (pas, with_y) in [(false, true), (true, false)] {
            let src = plus_locus(10, 2);
            let mut lines = src.clone();
            let out = run_on(&mut lines, Mode::F1, locus_evidence(false, with_y), pas);
            assert_eq!(lines, src, "pas {pas}, own starts {with_y}");
            assert!(out.bridges.is_empty());
            let row = out.junctions_tsv.lines().nth(1).unwrap().to_string();
            assert!(row.ends_with("\tFalse\tB"), "{row}");
            assert!(!row.contains("\tTrue\t3\tTrue\t"));
        }
    }

    /// The `-` strand mirrors every test: UP is the 3'-side genomic RIGHT, `fusion_of` runs 5' to 3'.
    #[test]
    fn the_minus_strand_mirrors() {
        let m = |ex: &[(i64, i64)]| -> Vec<(i64, i64)> { ex.iter().rev().map(|&(a, b)| (2001 - b, 2001 - a)).collect() };
        let mut lines = tx("c1", "G", "X", 10, &m(&[(101, 200), (301, 400)]), "-");
        lines.extend(tx("c1", "G", "B", 2, &m(&[(101, 200), (301, 350), (1301, 1400), (1501, 1600)]), "-"));
        lines.extend(tx("c1", "G", "Y", 8, &m(&[(1101, 1200), (1301, 1400), (1501, 1600)]), "-"));
        let out = run_on(&mut lines, Mode::F1, locus_evidence(true, true), true);
        assert!(lines[3].ends_with("fusion_of \"G,G.rg2\"; fusion_junction \"c1:701-1650:-\";"), "{}", lines[3]);
        assert_eq!(gene_of(&lines, "Y"), "G.rg2");
        assert_eq!(out.junctions_tsv.lines().nth(1).unwrap(), "G\tc1\t701\t1650\t-\t1\t2\t1\t1\t0\t3\t1601:3:P\tTrue\t3\tTrue\tTrue\tB");
    }

    /// F1v2: the link must carry fewer reads than EACH side (share < 1/2); a tie with a side abstains.
    #[test]
    fn f1v2_keeps_minority_bridges_and_a_tie_abstains() {
        let row = |rb: i64| {
            let mut lines = plus_locus(10, rb);
            let out = run_on(&mut lines, Mode::F1v2, locus_evidence(false, true), true);
            (out.bridges_tsv.unwrap().lines().nth(1).unwrap().to_string(), gene_of(&lines, "B"), out.stats)
        };
        let (r, g, st) = row(2);
        assert_eq!(r, "G\tc1\t351\t1300\t+\t1\t2\t2\t1\t10\t10\t1\t8\t8\t0.2\tTrue\tB");
        assert_eq!((g.as_str(), st.f1_bridge_junctions, st.bridge_junctions), ("G.fus1", 1, 1));
        let (r, g, st) = row(7);
        assert!(r.ends_with("\t0.4667\tTrue\tB"), "{r}");
        assert_eq!((g.as_str(), st.bridge_junctions), ("G.fus1", 1));
        // 8 = Y's 8: fewer than X's 10 but not fewer than Y's, share exactly 1/2 -> not a bridge; nothing changes
        let (r, g, st) = row(8);
        assert!(r.ends_with("\t0.5\tFalse\tB"), "{r}");
        assert_eq!((g.as_str(), st.f1_bridge_junctions, st.bridge_junctions, st.bridge_transcripts), ("G", 1, 0, 0));
        let (r, _, st) = row(12);
        assert!(r.ends_with("\t0.6\tFalse\tB"), "{r}");
        assert_eq!(st.bridge_junctions, 0);
    }

    #[test]
    fn f1v2_without_an_f1_bridge_writes_the_lone_gene_header() {
        let mut lines = plus_locus(10, 2);
        let out = run_on(&mut lines, Mode::F1v2, locus_evidence(false, false), true);
        assert_eq!(out.bridges_tsv.as_deref(), Some("gene\n"));
    }

    #[test]
    fn py_round4_is_python_round_then_repr() {
        let cases = [(1.0, "1.0"), (0.0, "0.0"), (0.5, "0.5"), (1.0 / 3.0, "0.3333"), (2.0 / 3.0, "0.6667")];
        for (x, want) in cases {
            assert_eq!(py_round4(x), want);
        }
        // an exact binary tie rounds half to even; an inexact one by its binary value (both as Python's round)
        assert_eq!(py_round4(1.0 / 32.0), "0.0312");
        assert_eq!(py_round4(3.0 / 32.0), "0.0938");
        assert_eq!(py_round4(1.0 / 160.0), "0.0063");
    }

    #[test]
    fn gene_ids_under_three_transcripts_and_gene_less_lines_are_left_alone() {
        let mut lines = vec!["# a comment".to_string()];
        lines.extend(tx("c1", "G", "X", 10, &[(101, 200), (301, 400)], "+"));
        lines.extend(tx("c1", "G", "B", 2, &[(101, 200), (301, 350), (1301, 1400), (1501, 1600)], "+"));
        lines.push("c1\trustle\ttranscript\t5\t9\t.\t+\t.\ttranscript_id \"N\"; reads \"3\";".to_string());
        lines.push("c1\trustle\texon\t5\t9\t.\t+\t.\ttranscript_id \"N\";".to_string());
        let src = lines.clone();
        let out = run_on(&mut lines, Mode::F1, locus_evidence(false, true), true);
        assert_eq!(lines, src);
        assert_eq!(out.junctions_tsv, format!("{JUNCTIONS_HEADER}\n"));
        assert!(src.iter().all(|l| out.in_families(l)));
        assert_eq!((out.stats.transcripts, out.stats.gene_ids), (3, 1));
    }

    #[test]
    fn malformed_gtfs_are_errors() {
        let p = |l: &[String]| parse(l, "--bridge-regroup");
        let mut dup = tx("c1", "G", "X", 1, &[(1, 10)], "+");
        dup.extend(tx("c1", "G", "X", 1, &[(20, 30)], "+"));
        assert!(p(&dup).is_err());
        let orphan = vec!["c1\trustle\texon\t1\t10\t.\t+\t.\tgene_id \"G\"; transcript_id \"X\";".to_string()];
        assert!(p(&orphan).is_err());
        let bare = vec!["c1\trustle\ttranscript\t1\t10\t.\t+\t.\tgene_id \"G\"; transcript_id \"X\";".to_string()];
        assert!(p(&bare).is_err(), "a transcript without exon lines");
        for (exons, what) in [
            (vec![(100, 200), (100, 200)], "a duplicated exon line"),
            (vec![(100, 200), (150, 300)], "overlapping exons"),
            (vec![(100, 200), (200, 300)], "exons sharing a base"),
            (vec![(100, 200), (300, 250)], "an exon whose end precedes its start"),
        ] {
            let e = p(&tx("c1", "G", "X", 1, &exons, "+")).err().unwrap_or_else(|| panic!("{what} must be an error"));
            assert!(e.to_string().contains("--bridge-regroup: transcript X has"), "{what}: {e}");
        }
        assert!(p(&tx("c1", "G", "X", 1, &[(100, 200), (201, 300)], "+")).is_ok(), "adjacent exons are ordered");
    }

    /// The review's scenario: a duplicated exon in a transcript whose junction is STRUCTURAL (a transcript wholly
    /// left of it and one wholly right), with a read ending on the duplicated exon's last base, so `up()`'s slice
    /// query would have had lo > hi. The pass returns the parse error instead of panicking.
    #[test]
    fn a_duplicated_exon_is_an_error_not_a_panic() {
        let mut lines = tx("c1", "G", "L", 1, &[(10, 50)], "+");
        lines.extend(tx("c1", "G", "D", 1, &[(100, 200), (100, 200)], "+"));
        lines.extend(tx("c1", "G", "R", 1, &[(300, 400)], "+"));
        let mut ev = BridgeEvidence::default();
        push(&mut ev, &[(100, 150), (160, 200)], false);
        let src = lines.clone();
        let err = run(&mut lines, Mode::F1, &mut ev, &genome, &fake_clusters(true)).err().expect("an error");
        assert!(err.to_string().contains("transcript D has overlapping or duplicated exons 100-200 and 100-200"), "{err}");
        assert_eq!(lines, src, "nothing was rewritten");
    }
}
