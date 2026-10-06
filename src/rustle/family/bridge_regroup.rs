//! ⭐ BRIDGE-AWARE REGROUPING of the assembled GTF: `copy_assign --assemble-only --bridge-regroup f1|f1v2|f1units`.
//!
//! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; `copy_assign --bridge-regroup`, default `f1v2` under `--assemble-only` since 2026-09-29, `off` = the 2026-09-25 products; driver `RUSTLE_BRIDGE_REGROUP=off|f1|f1v2|f1units`, unset = f1v2; `f1units` and `--bridge-units-list` are OPT-IN)
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
//! UNITS (`--bridge-regroup f1units`, OPT-IN; `docs/PREREG_container_units_v2_dev_2026-09-30.md` Part C, frozen rule U1). The
//! bridges above are dropped from the families input, so their copy and partner halves are never related. `f1units` keeps
//! every read and exon instead: F1's read evidence WITHOUT the share rule (a bridge junction is cut whatever its read
//! share: the rule existed because F1 removed the bridge) nominates the cut junctions, and each bridge transcript T with
//! cuts J_1 < .. < J_m is replaced IN THE FAMILIES INPUT ONLY (`<out>.families.gtf`, at T's line position) by m + 1 UNITS
//! `<T>.U1 ..` in transcription order (`-` strand: U1 is the rightmost), one per exon run between cuts, with T's `reads`
//! and the attributes `fusion_of "<T>"`, `fusion_unit "i/n"`, `fusion_junction` (every cut of T, genomic order),
//! `fusion_locus` (the pre-split locus key: the span of ALL transcripts of T's input gene_id), `fusion_gene` (that
//! gene_id), `fusion_detector` and, for F1, `fusion_evidence` (`reads_TJ;reads_up;reads_down;share` per cut, the cuts of
//! one transcript comma-joined). `<out>.gtf` under `f1units` is exactly what `--bridge-regroup f1` writes: the transcripts
//! that were cut stay whole there and are the only ones tagged `fusion_of` / `fusion_junction` (in `<g>.fus<k>`), and every
//! other gene_id gets RG3's names, so only `<out>.families.gtf` differs from f1's.
//! GENE_IDS of the families input. The assembler's own locus rule ([`native_components`] and [`best_rep`], property-tested
//! equal to `family_detect::collapse_loci_groups`: junction-sharing union-find keyed (contig, donor, acceptor),
//! strand-blind, representative max (reads, span, -line)) runs over the whole families input, and its names apply inside
//! every input gene_id that holds a unit, with each single-exon unit attached to the junction-bearing same-strand
//! transcript it shares most exonic bases with (ties: the lower line; overlap index over the whole families input, so a
//! unit can join a neighbouring gene's piece, which then takes the piece's name); names as the assembler names them with
//! the representative taken after attachment (the best component of each input gene keeps it, the others `<g>.nat<k>`).
//! Every other gene_id keeps today's RG3 names (the SCOPED form: the native rule alone, without a cut, regroups loci no cut
//! touched and loses members). No component spans two input gene_ids on the shipped GTFs; one that spanned a touched and
//! an untouched gene_id would give the touched gene's transcripts a name derived from its representative's gene_id.
//! `--bridge-units-list FILE` names the (transcript, junctions) to cut instead (any detector, e.g. annotation overlap in
//! GUIDED mode); F1's evidence is then not read at all, and a list with no transcript row is an error (run without it).
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

use crate::family::denovo_assemble::{rt_real_starts, rt_v1};
use crate::genome::GenomeIndex;

/// `--bridge-regroup`'s arms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    /// Every bridge junction (`f1_bridge.py --mode full`).
    F1,
    /// The bridge junctions whose link carries fewer reads than each side it separates (`f1v2.py --rule min`).
    F1v2,
    /// F1's bridge junctions with NO read-share rule, each bridge transcript cut into UNITS in the families input only
    /// (see the module header); `<out>.gtf` is F1's.
    F1Units,
}

impl Mode {
    /// `off` = `Ok(None)`; `f1` / `f1v2` / `f1units`; anything else is an error.
    pub fn parse(v: &str) -> Result<Option<Mode>> {
        match v {
            "off" => Ok(None),
            "f1" => Ok(Some(Mode::F1)),
            "f1v2" => Ok(Some(Mode::F1v2)),
            "f1units" => Ok(Some(Mode::F1Units)),
            other => {
                anyhow::bail!("--bridge-regroup must be off, f1, f1v2 or f1units, got `{other}`")
            }
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Mode::F1 => "f1",
            Mode::F1v2 => "f1v2",
            Mode::F1Units => "f1units",
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
    pub fn push(
        &mut self,
        chrom: &str,
        reverse: bool,
        ref_start: u64,
        ops: &[Op],
        ts_plus: impl FnOnce() -> bool,
    ) {
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
                Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    pos += len
                }
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
        let (e5, e3) = if minus {
            (ref_end, ref_start + 1)
        } else {
            (ref_start + 1, ref_end)
        };
        // the read's own first donor in transcript orientation (the script's `first_donor`)
        let donor = if minus {
            introns[introns.len() - 1].1 + 1
        } else {
            introns[0].0 - 1
        };
        c.rows[k].push((e5, e3, donor));
        if c.seen
            .insert(record_key(minus, ref_start, ref_end, introns))
        {
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
        self.contigs
            .values()
            .map(|c| c.rows[0].len() + c.rows[1].len())
            .sum()
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
pub fn bridge_evidence_region(
    bam_path: &str,
    chrom: &str,
    lo: u64,
    hi: u64,
    ev: &mut BridgeEvidence,
) -> Result<()> {
    let bai_path = format!("{bam_path}.bai");
    anyhow::ensure!(
        std::path::Path::new(&bai_path).exists(),
        "--bridge-regroup needs a .bai index"
    );
    let file = std::fs::File::open(bam_path)?;
    let buf = std::io::BufReader::with_capacity(1 << 20, file);
    let bgzf =
        noodles_bgzf::MultithreadedReader::with_worker_count(std::num::NonZeroUsize::MIN, buf);
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
        let Some(start) = record.alignment_start() else {
            continue;
        };
        let ref_start = (usize::from(start?) as u64).saturating_sub(1);
        ops.clear();
        for op in record.cigar().iter() {
            ops.push(op?);
        }
        ev.push(
            chrom,
            flags.is_reverse_complemented(),
            ref_start,
            &ops,
            || ts_is_plus(&record),
        );
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
#[derive(Clone)]
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
        let Some(tid) = attr(f[8], "transcript_id") else {
            continue;
        };
        if f[2] == "transcript" {
            anyhow::ensure!(
                !by_id.contains_key(tid),
                "{flag}: duplicate transcript line for {tid}"
            );
            let start: i64 = f[3]
                .parse()
                .with_context(|| format!("{flag}: bad start on transcript {tid}"))?;
            let end: i64 = f[4]
                .parse()
                .with_context(|| format!("{flag}: bad end on transcript {tid}"))?;
            // the scripts: `int(rv) if rv.lstrip('-').isdigit() else 0`
            let reads = attr(f[8], "reads")
                .filter(|v| {
                    v.strip_prefix('-')
                        .unwrap_or(v)
                        .chars()
                        .all(|c| c.is_ascii_digit())
                })
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
            let a: i64 = f[3]
                .parse()
                .with_context(|| format!("{flag}: bad exon start on {tid}"))?;
            let b: i64 = f[4]
                .parse()
                .with_context(|| format!("{flag}: bad exon end on {tid}"))?;
            txs[i].exons.push((a, b));
        }
    }
    for t in txs.iter_mut() {
        anyhow::ensure!(
            !t.exons.is_empty(),
            "{flag}: transcript {} has no exon line",
            t.tid
        );
        t.exons.sort_unstable();
        for &(a, b) in &t.exons {
            anyhow::ensure!(
                a <= b,
                "{flag}: transcript {} has an exon {a}-{b} whose end precedes its start",
                t.tid
            );
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
    let comps: Vec<Vec<usize>> = chroms
        .iter()
        .flat_map(|c| components(txs, &by_chrom[c]))
        .collect();
    let reps: Vec<usize> = comps
        .iter()
        .map(|c| {
            *c.iter()
                .max_by_key(|&&t| key(t))
                .expect("a component has a member")
        })
        .collect();
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
    let lo = comp
        .iter()
        .map(|&t| txs[t].exons[0].0)
        .min()
        .expect("a component has a member");
    let hi = comp
        .iter()
        .map(|&t| txs[t].exons[txs[t].exons.len() - 1].1)
        .max()
        .expect("a component has a member");
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
        let (lo, hi) = (
            ends.partition_point(|x| (x.0 as i64) < s),
            ends.partition_point(|x| (x.0 as i64) <= e),
        );
        let mut up: Vec<i64> = ends[lo..hi]
            .iter()
            .filter(|&&(_, e5)| {
                if minus {
                    (e5 as i64) > e
                } else {
                    (e5 as i64) < s
                }
            })
            .map(|&(e3, _)| if minus { -(e3 as i64) } else { e3 as i64 })
            .collect();
        up.sort_unstable();
        let seq: &[u8] = self
            .genome
            .chroms()
            .find(|(n, _)| *n == chrom)
            .map(|(_, s)| s)
            .unwrap_or(&[]);
        let cl = end_clusters(&up, seq, minus)
            .into_iter()
            .map(|c| {
                (
                    if minus { -c.mode } else { c.mode },
                    c.n,
                    c.proven,
                    c.unprimed,
                )
            })
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
        let key = (
            txs[i].gene.as_deref().expect("filtered above"),
            txs[i].chrom.as_str(),
        );
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
        let (gene, chrom) = (
            txs[ts[0]].gene.as_deref().expect("filtered above"),
            txs[ts[0]].chrom.as_str(),
        );
        let mut by_junction: BTreeMap<((i64, i64), &str), Vec<usize>> = BTreeMap::new();
        for &t in ts {
            for iv in txs[t].introns() {
                by_junction
                    .entry((iv, txs[t].strand.as_str()))
                    .or_default()
                    .push(t);
            }
        }
        for (((s, e), strand), tj) in by_junction {
            let tj_set: HashSet<usize> = tj.iter().copied().collect();
            let r: Vec<usize> = ts
                .iter()
                .copied()
                .filter(|t| txs[*t].strand == strand && !tj_set.contains(t))
                .collect();
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
                comps
                    .iter()
                    .zip(&sides)
                    .filter(|(_, &y)| y == x)
                    .flat_map(|(c, _)| c.iter().copied())
                    .collect()
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

/// The reads behind F1v2's MINORITY(J): T_J (the distinct transcripts using J) and the `reads` summed over T_J, over
/// the UP side and over the DOWN side.
fn link_reads(txs: &[Tx], j: &Junction) -> (Vec<usize>, i64, i64, i64) {
    let mut tj: Vec<usize> = j.tj.clone();
    tj.sort_unstable();
    tj.dedup();
    let sum = |v: &[usize]| v.iter().map(|&t| txs[t].reads).sum::<i64>();
    let (rb, ru, rd) = (sum(&tj), sum(&j.up_tx), sum(&j.down_tx));
    (tj, rb, ru, rd)
}

/// share = reads(T_J) / (reads(T_J) + min(reads(UP), reads(DOWN))); no reads anywhere = 1 (abstain).
fn link_share(rb: i64, ru: i64, rd: i64) -> f64 {
    let den = rb + ru.min(rd);
    if den > 0 {
        rb as f64 / den as f64
    } else {
        1.0
    }
}

/// `reads_TJ;reads_up;reads_down;share`: the `bridges.tsv` numbers of one bridge junction, which `f1units` carries
/// on its units as `fusion_evidence` (the share is reported, never applied).
fn evidence_string(txs: &[Tx], j: &Junction) -> String {
    let (_, rb, ru, rd) = link_reads(txs, j);
    format!("{rb};{ru};{rd};{}", py_round4(link_share(rb, ru, rd)))
}

/// F1v2's MINORITY(J) and its table row (`f1v2.py`'s `decide_v2`). T_J is the set of transcripts using J.
fn minority_row(txs: &[Tx], j: &Junction) -> (bool, String) {
    let (tj, rb, ru, rd) = link_reads(txs, j);
    let max = |v: &[usize]| {
        v.iter()
            .map(|&t| txs[t].reads)
            .max()
            .expect("a side has a transcript")
    };
    let share = link_share(rb, ru, rd);
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
        let non_bridge: Vec<usize> = ts
            .iter()
            .copied()
            .filter(|t| !bridges.contains(t))
            .collect();
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
                            && txs[x]
                                .exons
                                .iter()
                                .any(|&(a, b)| tx.exons.iter().any(|&(c0, d)| a <= d && c0 <= b))
                    });
                    if touches {
                        let lo = pc
                            .iter()
                            .map(|&x| txs[x].exons[0].0)
                            .min()
                            .expect("a component has a member");
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

// ============================================================================================ units (f1units)

/// One cut: the intron `(s, e)` (1-based closed, genomic) of a transcript at which it is split into units, and the
/// detector's evidence for it (F1: `reads_TJ;reads_up;reads_down;share`; a list: none).
#[derive(Clone, Debug, PartialEq, Eq)]
struct Cut {
    s: i64,
    e: i64,
    evidence: Option<String>,
}

/// What nominated the cuts of one pass, and the cuts: transcript index -> its cuts, sorted by position.
struct CutPlan {
    /// `f1` or `list:<file name>`: the units' `fusion_detector`, the units table's `label`, the relations' `detector`
    detector: String,
    cuts: BTreeMap<usize, Vec<Cut>>,
}

/// F1's cuts without the share rule: every transcript of T_J is cut at every kept bridge junction J.
fn f1_cut_plan(txs: &[Tx], junctions: &[Junction], kept: &[bool]) -> CutPlan {
    let mut cuts: BTreeMap<usize, Vec<Cut>> = BTreeMap::new();
    for (j, _) in junctions.iter().zip(kept).filter(|(_, &k)| k) {
        let evidence = evidence_string(txs, j);
        for &t in &j.tj {
            cuts.entry(t).or_default().push(Cut {
                s: j.s,
                e: j.e,
                evidence: Some(evidence.clone()),
            });
        }
    }
    for v in cuts.values_mut() {
        v.sort_by_key(|c| (c.s, c.e));
        v.dedup_by_key(|c| (c.s, c.e));
    }
    CutPlan {
        detector: "f1".to_string(),
        cuts,
    }
}

/// One cut of `--bridge-units-list`: the intron and, when the list wrote them, its contig and strand.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct ListedCut {
    pub s: i64,
    pub e: i64,
    pub contig: Option<String>,
    pub strand: Option<String>,
}

/// `--bridge-units-list FILE`: the transcripts to cut into units and where, for any detector (annotation overlap in
/// GUIDED mode: `bench/units_from_annotation.py`). A tab-separated file in `RUSTLE_READTHROUGH_JUNCTIONS=list:`'s
/// style: `#` lines and blank lines are skipped, the first other line is the header and must name a `tid` (or
/// `transcript_id`) column and a `junctions` (or `junction`) column; other columns are ignored. A row is a transcript id
/// and its introns to cut, as `S-E` (1-based closed, genomic) or `CONTIG:S-E:STRAND` (the `fusion_junction` form,
/// contig and strand then checked against the transcript), separated by `,` or `;`. A transcript listed twice gets the
/// union. Applied to the assembled GTF, every listed intron must be an intron of its transcript, and the list must
/// name at least one transcript of it.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct UnitsList {
    pub path: String,
    /// listed transcripts in order of first appearance, each with its distinct cuts in position order
    pub rows: Vec<(String, Vec<ListedCut>)>,
}

impl UnitsList {
    pub fn read(path: &str) -> Result<UnitsList> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("--bridge-units-list {path}: cannot read the list"))?;
        Self::parse(path, &text)
    }

    /// See [`UnitsList`]; `path` only names the file in errors.
    pub fn parse(path: &str, text: &str) -> Result<UnitsList> {
        let mut cols: Option<(usize, usize)> = None;
        let mut out = UnitsList {
            path: path.to_string(),
            rows: Vec::new(),
        };
        let mut at: HashMap<String, usize> = HashMap::new();
        for (i, line) in text.lines().enumerate() {
            let line = line.trim_end_matches('\r');
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            let Some((ti, ji)) = cols else {
                let col = |names: &[&str]| f.iter().position(|h| names.contains(&h.trim()));
                match (col(&["tid", "transcript_id"]), col(&["junctions", "junction"])) {
                    (Some(t), Some(j)) => cols = Some((t, j)),
                    _ => anyhow::bail!(
                        "{path}:{}: the units list's header must name tid|transcript_id and junctions|junction \
                         (tab-separated), got {line:?}",
                        i + 1
                    ),
                }
                continue;
            };
            let field = |k: usize| -> Result<&str> {
                f.get(k).map(|x| x.trim()).ok_or_else(|| {
                    anyhow::anyhow!(
                        "{path}:{}: {} field(s), column {} missing",
                        i + 1,
                        f.len(),
                        k + 1
                    )
                })
            };
            let tid = field(ti)?;
            anyhow::ensure!(!tid.is_empty(), "{path}:{}: empty transcript id", i + 1);
            let mut cuts: Vec<ListedCut> = Vec::new();
            for tok in field(ji)?
                .split([',', ';'])
                .map(str::trim)
                .filter(|t| !t.is_empty())
            {
                cuts.push(
                    parse_listed_cut(tok)
                        .map_err(|e| anyhow::anyhow!("{path}:{}: {tid}: {e}", i + 1))?,
                );
            }
            anyhow::ensure!(
                !cuts.is_empty(),
                "{path}:{}: {tid} names no junction",
                i + 1
            );
            let k = *at.entry(tid.to_string()).or_insert_with(|| {
                out.rows.push((tid.to_string(), Vec::new()));
                out.rows.len() - 1
            });
            out.rows[k].1.extend(cuts);
        }
        anyhow::ensure!(
            cols.is_some(),
            "{path}: no header line (empty list: run without --bridge-units-list)"
        );
        anyhow::ensure!(
            !out.rows.is_empty(),
            "{path}: the units list holds no transcript row (empty list: run without --bridge-units-list)"
        );
        for (_, cuts) in out.rows.iter_mut() {
            cuts.sort();
            cuts.dedup();
        }
        Ok(out)
    }

    /// `list:<file name>`: how the units and the relations name this detector.
    pub fn detector(&self) -> String {
        let name = std::path::Path::new(&self.path)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or(&self.path);
        format!("list:{name}")
    }
}

/// `S-E` or `CONTIG:S-E:STRAND` (`1 <= S <= E`, strand `+` or `-`).
fn parse_listed_cut(tok: &str) -> std::result::Result<ListedCut, String> {
    let (qual, range) = match tok.rsplit_once(':') {
        None => (None, tok),
        Some((rest, strand)) => {
            let (contig, range) = rest
                .rsplit_once(':')
                .ok_or_else(|| format!("{tok:?} is not S-E or CONTIG:S-E:STRAND"))?;
            if strand != "+" && strand != "-" {
                return Err(format!("{tok:?}: strand must be + or -"));
            }
            (Some((contig.to_string(), strand.to_string())), range)
        }
    };
    let (a, b) = range
        .split_once('-')
        .ok_or_else(|| format!("{tok:?} is not S-E or CONTIG:S-E:STRAND"))?;
    let (s, e) = match (a.parse::<i64>(), b.parse::<i64>()) {
        (Ok(s), Ok(e)) => (s, e),
        _ => return Err(format!("{tok:?}: S and E must be integers")),
    };
    if s < 1 || e < s {
        return Err(format!(
            "{tok:?}: the intron must be 1-based closed with 1 <= S <= E"
        ));
    }
    let (contig, strand) = qual.unzip();
    Ok(ListedCut {
        s,
        e,
        contig,
        strand,
    })
}

/// The cuts a list names, per transcript of the GTF, and how many listed transcripts the GTF lacks.
fn list_cut_plan(txs: &[Tx], list: &UnitsList) -> Result<(CutPlan, usize)> {
    let by_id: HashMap<&str, usize> = txs
        .iter()
        .enumerate()
        .map(|(i, t)| (t.tid.as_str(), i))
        .collect();
    let mut cuts: BTreeMap<usize, Vec<Cut>> = BTreeMap::new();
    let mut unmatched = 0usize;
    for (tid, listed) in &list.rows {
        let Some(&i) = by_id.get(tid.as_str()) else {
            unmatched += 1;
            continue;
        };
        let t = &txs[i];
        let introns: HashSet<(i64, i64)> = t.introns().collect();
        let mut v: Vec<Cut> = Vec::new();
        for c in listed {
            anyhow::ensure!(
                introns.contains(&(c.s, c.e)),
                "{}: {tid}: {}-{} is not an intron of that transcript",
                list.path,
                c.s,
                c.e
            );
            anyhow::ensure!(
                c.contig.as_deref().is_none_or(|x| x == t.chrom)
                    && c.strand.as_deref().is_none_or(|x| x == t.strand),
                "{}: {tid}: {}-{} names another contig or strand than the transcript ({}, {})",
                list.path,
                c.s,
                c.e,
                t.chrom,
                t.strand
            );
            v.push(Cut {
                s: c.s,
                e: c.e,
                evidence: None,
            });
        }
        v.sort_by_key(|c| (c.s, c.e));
        v.dedup_by_key(|c| (c.s, c.e));
        cuts.insert(i, v);
    }
    anyhow::ensure!(
        !cuts.is_empty(),
        "{}: none of its {} transcripts is in this GTF (a list written for another assembly?)",
        list.path,
        list.rows.len()
    );
    Ok((
        CutPlan {
            detector: list.detector(),
            cuts,
        },
        unmatched,
    ))
}

/// `t`'s exons cut at the introns `cuts`: one exon run per unit, in TRANSCRIPTION order (the genomic order reversed on
/// `-`), each in genomic order. m distinct cuts give m + 1 runs.
pub fn split_units(t: &Tx, cuts: &[(i64, i64)]) -> Vec<Vec<(i64, i64)>> {
    let cut: HashSet<(i64, i64)> = cuts.iter().copied().collect();
    let mut runs: Vec<Vec<(i64, i64)>> = Vec::new();
    let mut cur: Vec<(i64, i64)> = vec![t.exons[0]];
    for w in t.exons.windows(2) {
        if cut.contains(&(w[0].1 + 1, w[1].0 - 1)) {
            runs.push(std::mem::take(&mut cur));
        }
        cur.push(w[1]);
    }
    runs.push(cur);
    if t.strand == "-" {
        runs.reverse();
    }
    runs
}

/// The assembler's locus rule on GTF transcripts (`family_detect::collapse_loci_groups`, which `copy_assign` names
/// `gene_id`s with): the transcripts sharing one intron are one component, strand-blind. Components by their smallest
/// member, members ascending. A single-exon transcript is always a component of its own. The keys are `(contig,
/// exon_j.end, exon_{j+1}.start)` where the assembler's are `(contig, donor, acceptor)`: the same partition. Property-tested,
/// with [`best_rep`], against `collapse_loci_groups` (`native_components_equal_collapse_loci_groups`).
fn native_components(txs: &[Tx]) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..txs.len()).collect();
    let mut owner: HashMap<(&str, i64, i64), usize> = HashMap::new();
    for (i, t) in txs.iter().enumerate() {
        for w in t.exons.windows(2) {
            match owner.entry((t.chrom.as_str(), w[0].1, w[1].0)) {
                std::collections::hash_map::Entry::Occupied(o) => {
                    let (ra, rb) = (find(&mut parent, i), find(&mut parent, *o.get()));
                    if ra != rb {
                        parent[ra] = rb;
                    }
                }
                std::collections::hash_map::Entry::Vacant(v) => {
                    v.insert(i);
                }
            }
        }
    }
    let mut slot: HashMap<usize, usize> = HashMap::new();
    let mut out: Vec<Vec<usize>> = Vec::new();
    for i in 0..txs.len() {
        let r = find(&mut parent, i);
        let k = *slot.entry(r).or_insert_with(|| {
            out.push(Vec::new());
            out.len() - 1
        });
        out[k].push(i);
    }
    out
}

/// A component's representative: max (reads, span, -index), the assembler's own (most reads, then longest span, then
/// the earliest transcript; `collapse_loci_groups`' rule, property-tested with [`native_components`]).
fn best_rep(txs: &[Tx], members: &[usize]) -> usize {
    *members
        .iter()
        .max_by_key(|&&i| (txs[i].reads, txs[i].span, std::cmp::Reverse(i)))
        .expect("a component has a member")
}

/// The exons `(start, end, transcript)` of the junction-bearing transcripts of one contig and strand by start, the starts,
/// and the running maximum end (non-decreasing): [`attach_single_exon_units`]'s overlap index.
type ExonIndex = (Vec<(i64, i64, usize)>, Vec<i64>, Vec<i64>);

/// Step 3 of the native regroup: every single-exon UNIT (`is_unit`) joins the junction-bearing transcript of its
/// contig and strand with which it shares most exonic bases, ties to the lower index; `(unit, target)` in unit order.
/// The index covers every transcript with >= 2 exons of the families input (units included), so a unit can join a
/// neighbouring gene's transcript. A unit with no such overlap stays alone. Targets are never single-exon, so the
/// attachments are independent of each other and never link two components.
fn attach_single_exon_units(txs: &[Tx], is_unit: &[bool]) -> Vec<(usize, usize)> {
    let mut by: HashMap<(&str, &str), Vec<(i64, i64, usize)>> = HashMap::new();
    for (i, t) in txs.iter().enumerate() {
        if t.exons.len() >= 2 {
            by.entry((t.chrom.as_str(), t.strand.as_str()))
                .or_default()
                .extend(t.exons.iter().map(|&(a, b)| (a, b, i)));
        }
    }
    // per key: the exons by start, their starts, and the running maximum end (non-decreasing)
    let index: HashMap<(&str, &str), ExonIndex> = by
        .into_iter()
        .map(|(k, mut v)| {
            v.sort_unstable();
            let starts = v.iter().map(|x| x.0).collect();
            let mut top = i64::MIN;
            let prefix_max = v
                .iter()
                .map(|x| {
                    top = top.max(x.1);
                    top
                })
                .collect();
            (k, (v, starts, prefix_max))
        })
        .collect();
    let mut moved = Vec::new();
    for (i, t) in txs.iter().enumerate() {
        if t.exons.len() != 1 || !is_unit[i] {
            continue;
        }
        let Some((v, starts, prefix_max)) = index.get(&(t.chrom.as_str(), t.strand.as_str()))
        else {
            continue;
        };
        let (a, b) = t.exons[0];
        let (lo, hi) = (
            prefix_max.partition_point(|&m| m < a),
            starts.partition_point(|&s| s <= b),
        );
        let mut shared: BTreeMap<usize, i64> = BTreeMap::new();
        for &(ea, eb, j) in v.get(lo..hi).unwrap_or(&[]) {
            if eb >= a {
                *shared.entry(j).or_insert(0) += b.min(eb) - a.max(ea) + 1;
            }
        }
        if let Some((&j, _)) = shared
            .iter()
            .max_by_key(|(&j, &n)| (n, std::cmp::Reverse(j)))
        {
            moved.push((i, j));
        }
    }
    moved
}

/// Step 4 of the native regroup, per component: among the components whose representative has the same input
/// gene_id `g`, the one with the best representative keeps `g` and the others are `<g>.nat<k>`, k = 2.. in the order of
/// their representatives' indices, skipping a name the input already has. A component whose representative has no
/// gene_id has none.
fn name_groups(txs: &[Tx], groups: &[Vec<usize>]) -> Vec<Option<String>> {
    let key = |i: usize| (txs[i].reads, txs[i].span, std::cmp::Reverse(i));
    let reps: Vec<usize> = groups.iter().map(|m| best_rep(txs, m)).collect();
    let mut taken: HashSet<String> = txs.iter().filter_map(|t| t.gene.clone()).collect();
    let mut order: Vec<&str> = Vec::new();
    let mut by_gene: HashMap<&str, Vec<usize>> = HashMap::new();
    for (g, &r) in reps.iter().enumerate() {
        if let Some(o) = txs[r].gene.as_deref() {
            by_gene
                .entry(o)
                .or_insert_with(|| {
                    order.push(o);
                    Vec::new()
                })
                .push(g);
        }
    }
    let mut name: Vec<Option<String>> = vec![None; groups.len()];
    for o in order {
        let gl = &by_gene[o];
        let keep = *gl
            .iter()
            .max_by_key(|&&g| key(reps[g]))
            .expect("a gene has a component");
        let mut by_rep = gl.clone();
        by_rep.sort_unstable_by_key(|&g| reps[g]);
        let mut k = 2;
        for g in by_rep {
            if g == keep {
                name[g] = Some(o.to_string());
            } else {
                while taken.contains(&format!("{o}.nat{k}")) {
                    k += 1;
                }
                let nm = format!("{o}.nat{k}");
                taken.insert(nm.clone());
                name[g] = Some(nm);
                k += 1;
            }
        }
    }
    name
}

/// The native regroup of a whole families input: [`native_components`], the single-exon units attached, the components
/// named ([`name_groups`], the representative taken AFTER the attachment). Per transcript its name, and the attachments.
fn native_names(txs: &[Tx], is_unit: &[bool]) -> (Vec<Option<String>>, Vec<(usize, usize)>) {
    let mut groups = native_components(txs);
    let mut gid = vec![0usize; txs.len()];
    for (g, m) in groups.iter().enumerate() {
        for &i in m {
            gid[i] = g;
        }
    }
    let moved = attach_single_exon_units(txs, is_unit);
    for &(i, j) in &moved {
        debug_assert_eq!(
            groups[gid[i]],
            vec![i],
            "a single-exon transcript is a component of its own before attachment"
        );
        groups[gid[i]].clear();
        groups[gid[j]].push(i);
    }
    groups.retain(|m| !m.is_empty());
    for m in groups.iter_mut() {
        m.sort_unstable();
    }
    let names = name_groups(txs, &groups);
    let mut out: Vec<Option<String>> = vec![None; txs.len()];
    for (m, nm) in groups.iter().zip(&names) {
        for &i in m {
            out[i] = nm.clone();
        }
    }
    (out, moved)
}

/// The families input of `f1units`: every transcript in line order, each cut transcript replaced, at its position, by
/// its units.
struct FamInput {
    txs: Vec<Tx>,
    /// per families-input transcript: (unit number from 1, units of its parent); `(0, 0)` for an original transcript
    unit: Vec<(usize, usize)>,
    /// per original transcript: its first families-input transcript (the next original's is its end)
    first: Vec<usize>,
}

fn family_input(txs: &[Tx], plan: &CutPlan) -> FamInput {
    let mut out = FamInput {
        txs: Vec::new(),
        unit: Vec::new(),
        first: Vec::new(),
    };
    for (i, t) in txs.iter().enumerate() {
        out.first.push(out.txs.len());
        let Some(cuts) = plan.cuts.get(&i) else {
            out.txs.push(t.clone());
            out.unit.push((0, 0));
            continue;
        };
        let introns: Vec<(i64, i64)> = cuts.iter().map(|c| (c.s, c.e)).collect();
        let runs = split_units(t, &introns);
        let n = runs.len();
        for (k, exons) in runs.into_iter().enumerate() {
            out.txs.push(Tx {
                tid: format!("{}.U{}", t.tid, k + 1),
                gene: t.gene.clone(),
                chrom: t.chrom.clone(),
                strand: t.strand.clone(),
                reads: t.reads,
                span: exons[exons.len() - 1].1 - exons[0].0 + 1,
                exons,
            });
            out.unit.push((k + 1, n));
        }
    }
    out
}

/// The families-input `gene_id` of every families-input transcript (`None` = it has none and is never regrouped) and
/// the single-exon attachments `(unit, target)`. SCOPED: the native regroup names the transcripts of the input
/// gene_ids that hold a unit; every other gene_id keeps RG3's names ([`rg3_pieces`]: what the default arm writes), and
/// a unit attached to a transcript of such a gene_id takes that transcript's name.
fn regroup_units(fam: &FamInput) -> (Vec<Option<String>>, Vec<(usize, usize)>) {
    let txs = &fam.txs;
    let touched: HashSet<&str> = (0..txs.len())
        .filter(|&i| fam.unit[i].0 > 0)
        .filter_map(|i| txs[i].gene.as_deref())
        .collect();
    let mut names: Vec<Option<String>> = vec![None; txs.len()];
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
    for g in genes.iter().filter(|g| !touched.contains(*g)) {
        let (comps, name) = rg3_pieces(txs, g, &members[g]);
        for (c, nm) in comps.iter().zip(&name) {
            for &t in c {
                names[t] = Some(nm.clone());
            }
        }
    }
    let is_unit: Vec<bool> = fam.unit.iter().map(|u| u.0 > 0).collect();
    let (native, moved) = native_names(txs, &is_unit);
    for (i, t) in txs.iter().enumerate() {
        if t.gene.as_deref().is_some_and(|g| touched.contains(g)) {
            names[i] = native[i].clone();
        }
    }
    for &(i, j) in &moved {
        names[i] = names[j].clone();
    }
    (names, moved)
}

/// The products of the units pass besides the rewritten lines.
pub struct UnitsOutcome {
    /// the families input: the GTF lines, each cut transcript replaced at its position by its units' lines
    pub families_lines: Vec<String>,
    /// `<out>.bridge_units.tsv`: one row per unit
    pub table_tsv: String,
    /// `f1` or `list:<file name>`
    pub detector: String,
}

const UNITS_HEADER: &str =
    "unit_tid\tparent_tid\tinput_gene\tnew_gene\tunit\tof\tchrom\tstrand\tn_exon\texonic_bp\tstart\tend\treads\tlabel\tcuts\tevidence";

/// The counts of the units pass (folded into [`Stats`]).
#[derive(Default)]
struct UnitCounts {
    transcripts: usize,
    cuts: usize,
    units: usize,
    single_exon: usize,
    attached: usize,
    gene_ids_touched: usize,
    families_gene_ids: usize,
}

/// The attributes with the first `gene_id "<old>"` replaced by `new` (the scripts' rewrite): `None` when the line has no
/// gene_id, the transcript has no new name, or the two are equal (the line is unchanged).
fn with_gene(attrs: &str, old: Option<&str>, new: Option<&str>) -> Option<String> {
    match (old, new) {
        (Some(o), Some(n)) if o != n => {
            Some(attrs.replacen(&format!("gene_id \"{o}\""), &format!("gene_id \"{n}\""), 1))
        }
        _ => None,
    }
}

/// The units pass on the LAST-but-one state of the lines (before `finish` rewrites them for `<out>.gtf`): the
/// families input, the units table and the counts.
fn build_units(lines: &[String], txs: &[Tx], plan: &CutPlan) -> Result<(UnitsOutcome, UnitCounts)> {
    let fam = family_input(txs, plan);
    let mut seen: HashSet<&str> = HashSet::new();
    for t in &fam.txs {
        anyhow::ensure!(seen.insert(t.tid.as_str()), "--bridge-regroup f1units: the transcript id {} occurs twice (a unit id `<T>.U<i>` collides)", t.tid);
    }
    let (names, moved) = regroup_units(&fam);
    // the pre-split locus of every input gene_id: its first transcript's contig, the span of ALL its transcripts
    let mut locus: HashMap<&str, (&str, i64, i64)> = HashMap::new();
    for t in txs {
        if let Some(g) = t.gene.as_deref() {
            let (lo, hi) = (t.exons[0].0, t.exons[t.exons.len() - 1].1);
            let l = locus.entry(g).or_insert((t.chrom.as_str(), lo, hi));
            l.1 = l.1.min(lo);
            l.2 = l.2.max(hi);
        }
    }
    let id_of: HashMap<&str, usize> = txs
        .iter()
        .enumerate()
        .map(|(i, t)| (t.tid.as_str(), i))
        .collect();
    let mut out: Vec<String> = Vec::with_capacity(lines.len() + 2 * plan.cuts.len());
    let mut table = String::from(UNITS_HEADER);
    let mut counts = UnitCounts {
        transcripts: plan.cuts.len(),
        cuts: plan.cuts.values().map(Vec::len).sum(),
        attached: moved.len(),
        ..Default::default()
    };
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        let ti = if line.is_empty() || line.starts_with('#') || f.len() < 9 {
            None
        } else {
            attr(f[8], "transcript_id")
                .and_then(|t| id_of.get(t))
                .copied()
        };
        let Some(ti) = ti else {
            out.push(line.clone());
            continue;
        };
        let Some(cuts) = plan.cuts.get(&ti) else {
            match with_gene(f[8], attr(f[8], "gene_id"), names[fam.first[ti]].as_deref()) {
                Some(attrs) => {
                    let mut g = f.clone();
                    g[8] = &attrs;
                    out.push(g.join("\t"));
                }
                None => out.push(line.clone()),
            }
            continue;
        };
        if f[2] != "transcript" {
            continue; // the exon lines of a cut transcript: the units' own follow its transcript line
        }
        let t = &txs[ti];
        let junctions: Vec<String> = cuts
            .iter()
            .map(|c| format!("{}:{}-{}:{}", t.chrom, c.s, c.e, t.strand))
            .collect();
        let evidence: Option<String> = cuts
            .iter()
            .map(|c| c.evidence.as_deref())
            .collect::<Option<Vec<&str>>>()
            .map(|v| v.join(","));
        let cut_col = cuts
            .iter()
            .map(|c| format!("{}-{}", c.s, c.e))
            .collect::<Vec<_>>()
            .join(";");
        let (start, n_units) = (fam.first[ti], fam.unit[fam.first[ti]].1);
        for (fi, (u, name)) in fam
            .txs
            .iter()
            .zip(&names)
            .enumerate()
            .skip(start)
            .take(n_units)
        {
            let new_gene = name.as_deref();
            let (lo, hi) = (u.exons[0].0, u.exons[u.exons.len() - 1].1);
            let mut attrs = f[8].replacen(
                &format!("transcript_id \"{}\"", t.tid),
                &format!("transcript_id \"{}\"", u.tid),
                1,
            );
            attrs = with_gene(&attrs, t.gene.as_deref(), new_gene).unwrap_or(attrs);
            let mut tags = format!(
                " fusion_of \"{}\"; fusion_unit \"{}/{}\"; fusion_junction \"{}\";",
                t.tid,
                fam.unit[fi].0,
                fam.unit[fi].1,
                junctions.join(",")
            );
            if let Some(g) = t.gene.as_deref() {
                let (c, a, b) = locus[g];
                tags.push_str(&format!(
                    " fusion_locus \"{c}:{a}-{b}\"; fusion_gene \"{g}\";"
                ));
            }
            tags.push_str(&format!(" fusion_detector \"{}\";", plan.detector));
            if let Some(ev) = &evidence {
                tags.push_str(&format!(" fusion_evidence \"{ev}\";"));
            }
            let (lo_s, hi_s) = (lo.to_string(), hi.to_string());
            let full = format!("{}{tags}", attrs.trim_end());
            let mut g = f.clone();
            g[3] = &lo_s;
            g[4] = &hi_s;
            g[8] = &full;
            out.push(g.join("\t"));
            for (k, &(a, b)) in u.exons.iter().enumerate() {
                let gene = new_gene
                    .map(|n| format!("gene_id \"{n}\"; "))
                    .unwrap_or_default();
                out.push(format!(
                    "{}\t{}\texon\t{a}\t{b}\t.\t{}\t.\t{gene}transcript_id \"{}\"; exon_number \"{}\";",
                    t.chrom,
                    f[1],
                    t.strand,
                    u.tid,
                    k + 1
                ));
            }
            counts.units += 1;
            counts.single_exon += usize::from(u.exons.len() == 1);
            table.push_str(&format!(
                "\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{lo}\t{hi}\t{}\t{}\t{cut_col}\t{}",
                u.tid,
                t.tid,
                t.gene.as_deref().unwrap_or("."),
                new_gene.unwrap_or("."),
                fam.unit[fi].0,
                fam.unit[fi].1,
                t.chrom,
                t.strand,
                u.exons.len(),
                u.exons.iter().map(|&(a, b)| b - a + 1).sum::<i64>(),
                t.reads,
                plan.detector,
                evidence.as_deref().unwrap_or(".")
            ));
        }
    }
    table.push('\n');
    counts.gene_ids_touched = (0..fam.txs.len())
        .filter(|&i| fam.unit[i].0 > 0)
        .filter_map(|i| fam.txs[i].gene.as_deref())
        .collect::<HashSet<_>>()
        .len();
    counts.families_gene_ids = names.iter().flatten().collect::<HashSet<_>>().len();
    Ok((
        UnitsOutcome {
            families_lines: out,
            table_tsv: table,
            detector: plan.detector.clone(),
        },
        counts,
    ))
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
    /// the bridge junctions acted on (F1: all of F1's; F1v2: those passing MINORITY; f1units with a list: the listed)
    pub bridge_junctions: usize,
    pub bridge_transcripts: usize,
    pub gene_ids_with_bridge: usize,
    /// input gene_ids with a bridge or split into >= 2 pieces
    pub gene_ids_split: usize,
    pub gene_ids_after: usize,
    /// gene_ids of the families input (bridges excluded; under f1units: the units included)
    pub families_gene_ids: usize,
    /// lines whose gene_id changed
    pub lines_changed: usize,
    /// lines of bridge transcripts (absent from the families input; under f1units: replaced by the units' lines)
    pub family_lines_dropped: usize,
    /// f1units: transcripts cut into units
    pub unit_transcripts: usize,
    /// f1units: introns cut (over all transcripts)
    pub unit_cuts: usize,
    /// f1units: units written
    pub units: usize,
    /// f1units: units with one exon
    pub unit_single_exon: usize,
    /// f1units: single-exon units attached to a junction-bearing transcript's gene_id
    pub unit_attached: usize,
    /// f1units: input gene_ids holding a unit (the native regroup applies to them only)
    pub unit_gene_ids_touched: usize,
    /// f1units with a list: listed transcripts the GTF does not have
    pub unit_list_unmatched: usize,
}

/// The pass's products besides the rewritten lines.
pub struct Outcome {
    pub stats: Stats,
    /// transcript ids of the bridge transcripts
    pub bridges: HashSet<String>,
    /// every STRUCTURAL junction with its evidence (`f1_bridge.py`'s `junctions.tsv`); empty when no evidence was read
    pub junctions_tsv: String,
    /// was F1's read evidence consulted (every arm but `f1units` with a list)? Only then is `junctions_tsv` a table
    pub evidence_used: bool,
    /// every F1 bridge junction with its reads and MINORITY decision (`f1v2.py`'s `bridges.tsv`); F1v2 only
    pub bridges_tsv: Option<String>,
    /// `f1units`: the families input with the units, and the units table
    pub units: Option<UnitsOutcome>,
}

impl Outcome {
    /// Is `line` part of the families input? Every line except those of a bridge transcript. (Under `f1units` the
    /// families input is [`UnitsOutcome::families_lines`]: the bridges are replaced by their units, not dropped.)
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

/// The detector's own counts: F1's junction table gives them all, a units list only the junctions it names.
#[derive(Default)]
struct DetectorCounts {
    structural_junctions: usize,
    up_proof: usize,
    down_proof: usize,
    f1_bridge_junctions: usize,
    /// the bridge junctions acted on
    bridge_junctions: usize,
}

/// Run `--bridge-regroup` on the final lines of the emitted GTF: rewrite the `gene_id`s in place, append the relation
/// attributes to the bridge transcripts' `transcript` lines, and return the tables. `evidence` is consumed contig by
/// contig; `genome_of` gives a contig's sequence (uppercase); `end_clusters` is `--polish-tes`'s 3' cluster rule on
/// sorted oriented ends. `Mode::F1Units` acts on every F1 bridge junction and also returns the units ([`Outcome::units`]).
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
            .map(|&(m, n, p, u)| {
                format!(
                    "{m}:{n}:{}",
                    if p {
                        "P"
                    } else if u {
                        "u"
                    } else {
                        "-"
                    }
                )
            })
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
            if clusters.is_empty() {
                ".".to_string()
            } else {
                clusters.join(";")
            },
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
        let mut out = String::from(if rows.is_empty() {
            "gene"
        } else {
            BRIDGES_HEADER
        });
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
    let plan = (mode == Mode::F1Units).then(|| f1_cut_plan(&txs, &junctions, &kept));
    let counts = DetectorCounts {
        structural_junctions: junctions.len(),
        up_proof: junctions.iter().filter(|j| j.up_proof).count(),
        down_proof: junctions.iter().filter(|j| j.down_proof).count(),
        f1_bridge_junctions: junctions.iter().filter(|j| j.bridge).count(),
        bridge_junctions: kept.iter().filter(|&&k| k).count(),
    };
    finish(
        lines,
        &txs,
        &bridges,
        &bj,
        plan.as_ref(),
        &counts,
        (junctions_tsv, true),
        bridges_tsv,
    )
}

/// `--bridge-regroup f1units --bridge-units-list FILE`: the same units execution as [`run`]'s `Mode::F1Units`, with the
/// transcripts and introns to cut read from `list` instead of F1's evidence (which is not consulted: no junction table).
/// `<out>.gtf` is F1's regrouping with the listed transcripts as the bridges.
pub fn run_list(lines: &mut [String], list: &UnitsList) -> Result<Outcome> {
    let txs = parse(lines, "--bridge-regroup")?;
    let (plan, unmatched) = list_cut_plan(&txs, list)?;
    let bridges: HashSet<usize> = plan.cuts.keys().copied().collect();
    let mut bj: HashSet<(&str, i64, i64, &str)> = HashSet::new();
    for (&t, cuts) in &plan.cuts {
        for c in cuts {
            bj.insert((txs[t].chrom.as_str(), c.s, c.e, txs[t].strand.as_str()));
        }
    }
    let counts = DetectorCounts {
        bridge_junctions: bj.len(),
        ..Default::default()
    };
    let mut out = finish(
        lines,
        &txs,
        &bridges,
        &bj,
        Some(&plan),
        &counts,
        (String::new(), false),
        None,
    )?;
    out.stats.unit_list_unmatched = unmatched;
    Ok(out)
}

/// The rest of a pass, whatever named the bridges: the units (when `plan` is given; from the lines as they are, before
/// they are rewritten), the regrouped names, the rewrite of `lines`, the statistics. `tables` = (the junction table,
/// whether F1's evidence was read).
#[allow(clippy::too_many_arguments)]
fn finish(
    lines: &mut [String],
    txs: &[Tx],
    bridges: &HashSet<usize>,
    bj: &HashSet<(&str, i64, i64, &str)>,
    plan: Option<&CutPlan>,
    counts: &DetectorCounts,
    tables: (String, bool),
    bridges_tsv: Option<String>,
) -> Result<Outcome> {
    let units = match plan {
        Some(p) => Some(build_units(lines, txs, p)?),
        None => None,
    };
    let (new, rel) = regroup(txs, bridges, bj);

    // the scripts' `rewrite`: gene_id "<old>" -> gene_id "<new>" (first occurrence) on every line of a relabelled
    // transcript, then the relation appended to a bridge's `transcript` line
    let id_of: HashMap<&str, usize> = txs
        .iter()
        .enumerate()
        .map(|(i, t)| (t.tid.as_str(), i))
        .collect();
    let (mut lines_changed, mut family_lines_dropped) = (0usize, 0usize);
    for line in lines.iter_mut() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            continue;
        }
        let Some(&i) = attr(f[8], "transcript_id").and_then(|t| id_of.get(t)) else {
            continue;
        };
        let mut attrs = with_gene(f[8], attr(f[8], "gene_id"), new[i].as_deref());
        lines_changed += usize::from(attrs.is_some());
        if let Some((fusion_of, junction)) = rel.get(&i) {
            family_lines_dropped += 1;
            if f[2] == "transcript" {
                let base = attrs.take().unwrap_or_else(|| f[8].to_string());
                attrs = Some(format!(
                    "{} fusion_of \"{fusion_of}\"; fusion_junction \"{junction}\";",
                    base.trim_end()
                ));
            }
        }
        if let Some(a) = attrs {
            let mut out = f.clone();
            out[8] = &a;
            *line = out.join("\t");
        }
    }

    let bgenes: HashSet<&str> = bridges
        .iter()
        .filter_map(|&t| txs[t].gene.as_deref())
        .collect();
    let mut split: HashSet<&str> = txs
        .iter()
        .enumerate()
        .filter(|(i, t)| {
            !bridges.contains(i) && t.gene.is_some() && new[*i].as_deref() != t.gene.as_deref()
        })
        .filter_map(|(_, t)| t.gene.as_deref())
        .collect();
    split.extend(bgenes.iter().copied());
    let mut stats = Stats {
        transcripts: txs.len(),
        gene_ids: txs
            .iter()
            .filter_map(|t| t.gene.as_deref())
            .collect::<HashSet<_>>()
            .len(),
        structural_junctions: counts.structural_junctions,
        up_proof: counts.up_proof,
        down_proof: counts.down_proof,
        f1_bridge_junctions: counts.f1_bridge_junctions,
        bridge_junctions: counts.bridge_junctions,
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
        ..Default::default()
    };
    let units = units.map(|(outcome, c)| {
        stats.families_gene_ids = c.families_gene_ids;
        stats.unit_transcripts = c.transcripts;
        stats.unit_cuts = c.cuts;
        stats.units = c.units;
        stats.unit_single_exon = c.single_exon;
        stats.unit_attached = c.attached;
        stats.unit_gene_ids_touched = c.gene_ids_touched;
        outcome
    });
    Ok(Outcome {
        stats,
        bridges: bridges.iter().map(|&t| txs[t].tid.clone()).collect(),
        junctions_tsv: tables.0,
        evidence_used: tables.1,
        bridges_tsv,
        units,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A transcript's GTF lines (`transcript` + `exon`s; 1-based closed exons).
    fn tx(
        chrom: &str,
        gene: &str,
        tid: &str,
        reads: i64,
        exons: &[(i64, i64)],
        strand: &str,
    ) -> Vec<String> {
        let a =
            format!("gene_id \"{gene}\"; transcript_id \"{tid}\"; reads \"{reads}\"; TPM \"1.0\";");
        let (lo, hi) = (exons[0].0, exons[exons.len() - 1].1);
        let mut v = vec![format!(
            "{chrom}\trustle\ttranscript\t{lo}\t{hi}\t.\t{strand}\t.\t{a}"
        )];
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
            vec![EndCluster {
                mode: *ends.last().unwrap(),
                n: ends.len(),
                proven: pas,
                unprimed: true,
            }]
        }
    }

    fn genome(_: &str) -> Result<Arc<GenomeIndex>> {
        Ok(Arc::new(GenomeIndex::from_seqs(&[("c1", &[b'C'; 3000])])))
    }

    /// The `+` locus: X (reads `rx`) ends inside J = [351, 1300]; the bridge B (reads `rb`) joins X's first exon to
    /// Y's last two exons through J; Y (reads 8) starts inside J at its own promoter.
    fn plus_locus(rx: i64, rb: i64) -> Vec<String> {
        let mut v = tx("c1", "G", "X", rx, &[(101, 200), (301, 400)], "+");
        v.extend(tx(
            "c1",
            "G",
            "B",
            rb,
            &[(101, 200), (301, 350), (1301, 1400), (1501, 1600)],
            "+",
        ));
        v.extend(tx(
            "c1",
            "G",
            "Y",
            8,
            &[(1101, 1200), (1301, 1400), (1501, 1600)],
            "+",
        ));
        v
    }

    /// X's reads (3' ends 398-400 inside J, deduplicated to three) and Y's (three at one start: a real cluster, first
    /// exon inside J), mirrored onto the `-` strand when `minus` (x -> 2001 - x).
    fn locus_evidence(minus: bool, with_y: bool) -> BridgeEvidence {
        let m = |ex: &[(u64, u64)]| -> Vec<(u64, u64)> {
            if minus {
                ex.iter()
                    .rev()
                    .map(|&(a, b)| (2001 - b, 2001 - a))
                    .collect()
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
                push(
                    &mut ev,
                    &m(&[(1101, 1200), (1301, 1400), (1501, 1600)]),
                    minus,
                );
            }
        }
        ev
    }

    fn run_on(lines: &mut [String], mode: Mode, mut ev: BridgeEvidence, pas: bool) -> Outcome {
        run(lines, mode, &mut ev, &genome, &fake_clusters(pas)).expect("run")
    }

    fn gene_of(lines: &[String], tid: &str) -> String {
        let l = lines
            .iter()
            .find(|l| l.contains("\ttranscript\t") && attr(l, "transcript_id") == Some(tid))
            .unwrap();
        attr(l, "gene_id").unwrap().to_string()
    }

    #[test]
    fn evidence_takes_every_n_and_ends_one_past_the_last_reference_base() {
        let mut ev = BridgeEvidence::default();
        let (m, d, n, i, s) = (
            Kind::Match,
            Kind::Deletion,
            Kind::Skip,
            Kind::Insertion,
            Kind::SoftClip,
        );
        // 10M 5D 100N 3I 20M 2S from 1-based 1000: intron [1015, 1114], last base 1134, first donor 1014
        let ops = [
            Op::new(m, 10),
            Op::new(d, 5),
            Op::new(n, 100),
            Op::new(i, 3),
            Op::new(m, 20),
            Op::new(s, 2),
        ];
        ev.push("c1", false, 999, &ops, || true);
        assert_eq!(ev.contigs["c1"].rows[0], vec![(1000, 1134, 1014)]);
        assert_eq!(ev.contigs["c1"].ends[0], vec![(1134, 1000)]);
        // a leading N and N-I-N are introns of their own (the script's walk, not the assembler's exon blocks)
        let ops = [
            Op::new(n, 50),
            Op::new(m, 10),
            Op::new(n, 20),
            Op::new(i, 2),
            Op::new(n, 30),
            Op::new(m, 10),
        ];
        ev.push("c2", false, 0, &ops, || true);
        assert_eq!(ev.contigs["c2"].rows[0], vec![(1, 120, 0)]);
        // `=` / `X` advance the reference; an unspliced record is never kept and its `ts` never read
        ev.push(
            "c3",
            false,
            0,
            &[
                Op::new(Kind::SequenceMatch, 5),
                Op::new(Kind::SequenceMismatch, 1),
            ],
            || panic!("ts read"),
        );
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
        assert_eq!(
            (c.ends[0].len(), c.ends[1].len()),
            (1, 1),
            "one U record per strand after deduplication"
        );
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
        let (mut all, mut region, mut other) = (
            BridgeEvidence::default(),
            locus_evidence(false, true),
            BridgeEvidence::default(),
        );
        region.seal();
        assert!(region.contigs["c1"].seen.is_empty());
        assert_eq!(
            (region.contigs["c1"].ends[0].len(), region.len()),
            (4, 7),
            "the evidence itself is kept"
        );
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
            noodles_sam::alignment::RecordBuf::builder()
                .set_data(data)
                .build()
        };
        assert!(ts_is_plus(&rec(None)), "absent = +");
        assert!(ts_is_plus(&rec(Some(BufValue::Character(b'+')))));
        assert!(!ts_is_plus(&rec(Some(BufValue::Character(b'-')))));
        assert!(
            ts_is_plus(&rec(Some(BufValue::String("+".into())))),
            "ts:Z:+ compares by value"
        );
        assert!(
            !ts_is_plus(&rec(Some(BufValue::Int32(43)))),
            "another type is not '+'"
        );
    }

    #[test]
    fn attr_is_the_first_key_quote_match() {
        let a = "gene_id \"G\"; matched_reads \"0\"; reads \"7\";";
        assert_eq!(
            attr(a, "reads"),
            Some("0"),
            "as the scripts' attr: the first `reads \"` wins"
        );
        assert_eq!(
            attr("reads \"7\"; matched_reads \"0\";", "reads"),
            Some("7")
        );
        assert_eq!(
            attr(
                "transcript_idx \"a\"; transcript_id \"b\";",
                "transcript_id"
            ),
            Some("b")
        );
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
        assert_eq!(
            (s(0, "+"), s(1, "+"), s(2, "+"), s(3, "+")),
            (Side::Up, Side::Down, Side::Straddle, Side::Inside)
        );
        assert_eq!(
            (s(0, "-"), s(1, "-"), s(2, "-"), s(3, "-")),
            (Side::Down, Side::Up, Side::Straddle, Side::Inside)
        );
    }

    #[test]
    fn components_join_same_strand_exon_overlap_only() {
        let mut lines = tx("c1", "G", "a", 1, &[(100, 200), (500, 600)], "+");
        lines.extend(tx("c1", "G", "b", 1, &[(300, 400)], "+")); // inside a's intron: not overlap
        lines.extend(tx("c1", "G", "c", 1, &[(600, 700)], "+")); // one shared base with a
        lines.extend(tx("c1", "G", "d", 1, &[(650, 800)], "-")); // other strand
        lines.extend(tx("c1", "G", "e", 1, &[(790, 900)], "+")); // chains through nothing on its strand
        let txs = parse(&lines, "--bridge-regroup").unwrap();
        assert_eq!(
            components(&txs, &[0, 1, 2, 3, 4]),
            vec![vec![0, 2], vec![1], vec![3], vec![4]]
        );
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
        assert_eq!(
            names,
            vec!["G.rg2".to_string(), "G.rg3".to_string(), "G".to_string()]
        );
        assert_eq!(
            rg3_pieces(&txs, "G", &[0, 3]),
            (vec![vec![0, 3]], vec!["G".to_string()])
        );
        assert_eq!(rg3_pieces(&txs, "G", &[]), (Vec::new(), Vec::new()));
    }

    /// The bridge: X keeps the name (10 reads > 8), Y becomes `.rg2`, B `.fus1` with its relation; the junction table
    /// holds the one STRUCTURAL junction; nothing but gene_id and the two appended attributes changes.
    #[test]
    fn f1_splits_a_proven_bridge_and_keeps_it_as_a_relation() {
        let src = plus_locus(10, 2);
        let mut lines = src.clone();
        let out = run_on(&mut lines, Mode::F1, locus_evidence(false, true), true);
        assert_eq!(
            (
                gene_of(&lines, "X"),
                gene_of(&lines, "Y"),
                gene_of(&lines, "B")
            ),
            ("G".into(), "G.rg2".into(), "G.fus1".into())
        );
        assert_eq!(
            lines[3],
            "c1\trustle\ttranscript\t101\t1600\t.\t+\t.\tgene_id \"G.fus1\"; transcript_id \"B\"; reads \"2\"; TPM \"1.0\"; \
             fusion_of \"G,G.rg2\"; fusion_junction \"c1:351-1300:+\";"
        );
        assert_eq!(lines.len(), src.len());
        for (a, b) in src.iter().zip(&lines) {
            assert_eq!(
                a.split('\t').take(8).collect::<Vec<_>>(),
                b.split('\t').take(8).collect::<Vec<_>>(),
                "coordinates never move"
            );
        }
        assert_eq!(
            out.junctions_tsv,
            format!("{JUNCTIONS_HEADER}\nG\tc1\t351\t1300\t+\t1\t2\t1\t1\t0\t3\t400:3:P\tTrue\t3\tTrue\tTrue\tB\n")
        );
        assert!(out.bridges_tsv.is_none());
        let fam: Vec<&String> = lines.iter().filter(|l| out.in_families(l)).collect();
        assert_eq!(
            fam.len(),
            src.len() - 5,
            "the bridge's transcript and 4 exon lines leave the families input"
        );
        assert!(fam.iter().all(|l| attr(l, "transcript_id") != Some("B")));
        let st = &out.stats;
        assert_eq!(
            (
                st.structural_junctions,
                st.f1_bridge_junctions,
                st.bridge_junctions,
                st.bridge_transcripts
            ),
            (1, 1, 1, 1)
        );
        assert_eq!(
            (
                st.gene_ids,
                st.gene_ids_split,
                st.gene_ids_after,
                st.families_gene_ids
            ),
            (1, 1, 3, 2)
        );
        assert_eq!(
            (st.lines_changed, st.family_lines_dropped),
            (9, 5),
            "Y: 4 lines relabelled, B: 5"
        );
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
        let m = |ex: &[(i64, i64)]| -> Vec<(i64, i64)> {
            ex.iter()
                .rev()
                .map(|&(a, b)| (2001 - b, 2001 - a))
                .collect()
        };
        let mut lines = tx("c1", "G", "X", 10, &m(&[(101, 200), (301, 400)]), "-");
        lines.extend(tx(
            "c1",
            "G",
            "B",
            2,
            &m(&[(101, 200), (301, 350), (1301, 1400), (1501, 1600)]),
            "-",
        ));
        lines.extend(tx(
            "c1",
            "G",
            "Y",
            8,
            &m(&[(1101, 1200), (1301, 1400), (1501, 1600)]),
            "-",
        ));
        let out = run_on(&mut lines, Mode::F1, locus_evidence(true, true), true);
        assert!(
            lines[3].ends_with("fusion_of \"G,G.rg2\"; fusion_junction \"c1:701-1650:-\";"),
            "{}",
            lines[3]
        );
        assert_eq!(gene_of(&lines, "Y"), "G.rg2");
        assert_eq!(
            out.junctions_tsv.lines().nth(1).unwrap(),
            "G\tc1\t701\t1650\t-\t1\t2\t1\t1\t0\t3\t1601:3:P\tTrue\t3\tTrue\tTrue\tB"
        );
    }

    /// F1v2: the link must carry fewer reads than EACH side (share < 1/2); a tie with a side abstains.
    #[test]
    fn f1v2_keeps_minority_bridges_and_a_tie_abstains() {
        let row = |rb: i64| {
            let mut lines = plus_locus(10, rb);
            let out = run_on(&mut lines, Mode::F1v2, locus_evidence(false, true), true);
            (
                out.bridges_tsv.unwrap().lines().nth(1).unwrap().to_string(),
                gene_of(&lines, "B"),
                out.stats,
            )
        };
        let (r, g, st) = row(2);
        assert_eq!(
            r,
            "G\tc1\t351\t1300\t+\t1\t2\t2\t1\t10\t10\t1\t8\t8\t0.2\tTrue\tB"
        );
        assert_eq!(
            (g.as_str(), st.f1_bridge_junctions, st.bridge_junctions),
            ("G.fus1", 1, 1)
        );
        let (r, g, st) = row(7);
        assert!(r.ends_with("\t0.4667\tTrue\tB"), "{r}");
        assert_eq!((g.as_str(), st.bridge_junctions), ("G.fus1", 1));
        // 8 = Y's 8: fewer than X's 10 but not fewer than Y's, share exactly 1/2 -> not a bridge; nothing changes
        let (r, g, st) = row(8);
        assert!(r.ends_with("\t0.5\tFalse\tB"), "{r}");
        assert_eq!(
            (
                g.as_str(),
                st.f1_bridge_junctions,
                st.bridge_junctions,
                st.bridge_transcripts
            ),
            ("G", 1, 0, 0)
        );
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
        let cases = [
            (1.0, "1.0"),
            (0.0, "0.0"),
            (0.5, "0.5"),
            (1.0 / 3.0, "0.3333"),
            (2.0 / 3.0, "0.6667"),
        ];
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
        lines.extend(tx(
            "c1",
            "G",
            "B",
            2,
            &[(101, 200), (301, 350), (1301, 1400), (1501, 1600)],
            "+",
        ));
        lines.push(
            "c1\trustle\ttranscript\t5\t9\t.\t+\t.\ttranscript_id \"N\"; reads \"3\";".to_string(),
        );
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
        let orphan = vec![
            "c1\trustle\texon\t1\t10\t.\t+\t.\tgene_id \"G\"; transcript_id \"X\";".to_string(),
        ];
        assert!(p(&orphan).is_err());
        let bare = vec![
            "c1\trustle\ttranscript\t1\t10\t.\t+\t.\tgene_id \"G\"; transcript_id \"X\";"
                .to_string(),
        ];
        assert!(p(&bare).is_err(), "a transcript without exon lines");
        for (exons, what) in [
            (vec![(100, 200), (100, 200)], "a duplicated exon line"),
            (vec![(100, 200), (150, 300)], "overlapping exons"),
            (vec![(100, 200), (200, 300)], "exons sharing a base"),
            (
                vec![(100, 200), (300, 250)],
                "an exon whose end precedes its start",
            ),
        ] {
            let e = p(&tx("c1", "G", "X", 1, &exons, "+"))
                .err()
                .unwrap_or_else(|| panic!("{what} must be an error"));
            assert!(
                e.to_string().contains("--bridge-regroup: transcript X has"),
                "{what}: {e}"
            );
        }
        assert!(
            p(&tx("c1", "G", "X", 1, &[(100, 200), (201, 300)], "+")).is_ok(),
            "adjacent exons are ordered"
        );
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
        let err = run(&mut lines, Mode::F1, &mut ev, &genome, &fake_clusters(true))
            .err()
            .expect("an error");
        assert!(
            err.to_string()
                .contains("transcript D has overlapping or duplicated exons 100-200 and 100-200"),
            "{err}"
        );
        assert_eq!(lines, src, "nothing was rewritten");
    }

    // ================================================================================================ units (f1units)

    fn mk(
        tid: &str,
        gene: &str,
        chrom: &str,
        strand: &str,
        reads: i64,
        exons: &[(i64, i64)],
    ) -> Tx {
        Tx {
            tid: tid.into(),
            gene: Some(gene.into()),
            chrom: chrom.into(),
            strand: strand.into(),
            reads,
            span: exons[exons.len() - 1].1 - exons[0].0 + 1,
            exons: exons.to_vec(),
        }
    }

    /// Each transcript's component representative, the shape `collapse_loci_groups` returns: the production helpers
    /// ([`native_components`], [`best_rep`]) the families input is named with.
    fn reps_of(txs: &[Tx]) -> Vec<usize> {
        let mut rep = vec![0usize; txs.len()];
        for members in native_components(txs) {
            let r = best_rep(txs, &members);
            for &m in &members {
                rep[m] = r;
            }
        }
        rep
    }

    /// xorshift64: the property test's fixed-seed generator
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// `native_components` + `best_rep` ARE `family_detect::collapse_loci_groups`, the rule the assembler names `gene_id`s
    /// with: on random transcripts over a small junction grid (so junctions are shared), on two contigs and both strands (a
    /// junction shared across strands joins them), with raw-tid collisions (equal ids, different components), read and span
    /// ties. `name_groups`, which picks each component's representative itself, names a component (one gene per transcript
    /// here) after the representative `collapse_loci_groups` chose.
    #[test]
    fn native_components_equal_collapse_loci_groups() {
        use crate::family::family_detect::{collapse_loci_groups, DenovoTranscript};
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        let (mut shared, mut cross_strand, mut collisions, mut ties) =
            (0usize, 0usize, 0usize, 0usize);
        for round in 0..500 {
            let n = 4 + rng.below(40) as usize;
            let mut txs: Vec<Tx> = Vec::new();
            for i in 0..n {
                let chrom = if rng.below(4) == 0 { "c2" } else { "c1" };
                let strand = if rng.below(2) == 0 { "+" } else { "-" };
                let nex = 1 + rng.below(4) as usize;
                let mut pos: Vec<i64> = Vec::new();
                while pos.len() < 2 * nex {
                    let p = 100 * (1 + rng.below(12)) as i64;
                    if !pos.contains(&p) {
                        pos.push(p);
                    }
                }
                pos.sort_unstable();
                let exons: Vec<(i64, i64)> = pos.chunks(2).map(|c| (c[0] + 1, c[1])).collect();
                let tid = format!("DN_{chrom}_{}_{nex}", exons[0].0 - 1);
                txs.push(mk(
                    &tid,
                    &format!("g{i}"),
                    chrom,
                    strand,
                    1 + rng.below(4) as i64,
                    &exons,
                ));
            }
            let den: Vec<DenovoTranscript> = txs
                .iter()
                .map(|t| DenovoTranscript {
                    tid: t.tid.clone(),
                    chrom: t.chrom.clone(),
                    start: (t.exons[0].0 - 1) as u64,
                    end: t.exons[t.exons.len() - 1].1 as u64,
                    n_reads: t.reads as u32,
                    strand: t.strand.chars().next().unwrap(),
                    introns: t
                        .exons
                        .windows(2)
                        .map(|w| (w[0].1 as u64, (w[1].0 - 1) as u64))
                        .collect(),
                    ..Default::default()
                })
                .collect();
            let want = collapse_loci_groups(&den);
            let got = reps_of(&txs);
            assert_eq!(got, want, "round {round}");
            let comps = native_components(&txs);
            for (m, name) in comps.iter().zip(name_groups(&txs, &comps)) {
                assert_eq!(
                    name.as_deref(),
                    Some(format!("g{}", want[m[0]]).as_str()),
                    "round {round}: the component of {m:?}"
                );
            }
            let reps: std::collections::BTreeSet<usize> = got.iter().copied().collect();
            shared += usize::from(reps.len() < txs.len());
            cross_strand += usize::from(
                (0..n).any(|i| (0..n).any(|j| got[i] == got[j] && txs[i].strand != txs[j].strand)),
            );
            collisions += usize::from(
                reps.iter()
                    .any(|&a| reps.iter().any(|&b| a < b && txs[a].tid == txs[b].tid)),
            );
            ties += usize::from((0..n).any(|i| {
                (0..n).any(|j| {
                    i < j
                        && got[i] == got[j]
                        && txs[i].reads == txs[j].reads
                        && txs[i].span == txs[j].span
                })
            }));
        }
        assert!(
            shared > 400 && cross_strand > 100 && collisions > 10 && ties > 10,
            "{shared} {cross_strand} {collisions} {ties}"
        );
    }

    #[test]
    fn native_components_by_hand() {
        // a and b share the intron 201-299: one component; its representative is the most-read transcript (b)
        let t = [
            mk("a", "g1", "c1", "+", 5, &[(100, 200), (300, 400)]),
            mk(
                "b",
                "g1",
                "c1",
                "+",
                9,
                &[(100, 200), (300, 450), (600, 700)],
            ),
            mk("c", "g2", "c1", "+", 9, &[(1000, 1100), (1200, 1300)]),
        ];
        assert_eq!(native_components(&t), vec![vec![0, 1], vec![2]]);
        assert_eq!(reps_of(&t), vec![1, 1, 2]);
        // reads tie -> longest span -> earliest index
        let t = [
            mk("a", "g", "c1", "+", 5, &[(100, 200), (300, 400)]),
            mk("b", "g", "c1", "+", 5, &[(100, 200), (300, 400)]),
        ];
        assert_eq!(reps_of(&t), vec![0, 0]);
        let t = [
            mk("a", "g", "c1", "+", 5, &[(100, 200), (300, 400)]),
            mk(
                "b",
                "g",
                "c1",
                "+",
                5,
                &[(100, 200), (300, 500), (700, 800)],
            ),
        ];
        assert_eq!(reps_of(&t), vec![1, 1], "the longer span wins the read tie");
        // strand-blind, contig-keyed
        let mut t = [
            mk("a", "g1", "c1", "+", 5, &[(100, 200), (300, 400)]),
            mk("b", "g2", "c1", "-", 5, &[(100, 200), (300, 400)]),
        ];
        assert_eq!(native_components(&t), vec![vec![0, 1]]);
        t[1].chrom = "c2".into();
        assert_eq!(native_components(&t), vec![vec![0], vec![1]]);
        // a single-exon transcript has no junction: a component of its own even inside another's exon
        let t = [
            mk("a", "g", "c1", "+", 5, &[(100, 200), (300, 400)]),
            mk("s", "g", "c1", "+", 5, &[(120, 180)]),
        ];
        assert_eq!(native_components(&t), vec![vec![0], vec![1]]);
    }

    #[test]
    fn split_units_cut_at_the_introns_in_transcription_order() {
        let ex = [(100, 200), (300, 400), (500, 600), (700, 800)];
        let plus = mk("t", "g", "c1", "+", 4, &ex);
        let minus = mk("t", "g", "c1", "-", 4, &ex);
        // m = 1
        let runs = split_units(&plus, &[(401, 499)]);
        assert_eq!(
            runs,
            vec![vec![(100, 200), (300, 400)], vec![(500, 600), (700, 800)]]
        );
        let runs = split_units(&minus, &[(401, 499)]);
        assert_eq!(
            runs,
            vec![vec![(500, 600), (700, 800)], vec![(100, 200), (300, 400)]],
            "`-`: U1 is the genomic right"
        );
        // m = 2: a run with no cut inside stays whole, the cuts are given in any order
        let runs = split_units(&plus, &[(601, 699), (201, 299)]);
        assert_eq!(
            runs,
            vec![
                vec![(100, 200)],
                vec![(300, 400), (500, 600)],
                vec![(700, 800)]
            ]
        );
        let runs = split_units(&minus, &[(601, 699), (201, 299)]);
        assert_eq!(
            runs,
            vec![
                vec![(700, 800)],
                vec![(300, 400), (500, 600)],
                vec![(100, 200)]
            ]
        );
        // m = 3: every exon is a unit
        let all = [(201, 299), (401, 499), (601, 699)];
        assert_eq!(split_units(&plus, &all).concat(), ex.to_vec());
        assert_eq!(split_units(&plus, &all).len(), 4);
        assert_eq!(
            split_units(&minus, &all).concat(),
            ex.iter().rev().copied().collect::<Vec<_>>()
        );
        // every exon of T is in exactly one unit, in every case
        for cuts in [&all[..1], &all[1..], &all[..2], &all[..]] {
            for t in [&plus, &minus] {
                let mut got: Vec<(i64, i64)> = split_units(t, cuts).concat();
                got.sort_unstable();
                assert_eq!(got, ex.to_vec());
                assert_eq!(split_units(t, cuts).len(), cuts.len() + 1);
            }
        }
        // a cut that is not an intron of the transcript cuts nothing
        assert_eq!(split_units(&plus, &[(1, 2)]).len(), 1);
    }

    #[test]
    fn the_units_of_a_cut_transcript_take_its_place_its_reads_and_its_gene() {
        let txs = vec![
            mk("A", "g", "c1", "+", 30, &[(100, 200), (300, 400)]),
            mk(
                "F",
                "g",
                "c1",
                "-",
                7,
                &[(100, 200), (300, 400), (2000, 2100)],
            ),
            mk("B", "g", "c1", "+", 30, &[(2000, 2100), (2200, 2300)]),
        ];
        let plan = CutPlan {
            detector: "x".into(),
            cuts: BTreeMap::from([(
                1,
                vec![Cut {
                    s: 401,
                    e: 1999,
                    evidence: None,
                }],
            )]),
        };
        let fam = family_input(&txs, &plan);
        let ids: Vec<&str> = fam.txs.iter().map(|t| t.tid.as_str()).collect();
        assert_eq!(ids, ["A", "F.U1", "F.U2", "B"]);
        assert_eq!(fam.unit, vec![(0, 0), (1, 2), (2, 2), (0, 0)]);
        assert_eq!(fam.first, vec![0, 1, 3]);
        // `-`: U1 is the genomic right (the single exon), U2 the two left exons; reads and gene are T's, the span the unit's
        assert_eq!(fam.txs[1].exons, vec![(2000, 2100)]);
        assert_eq!(fam.txs[2].exons, vec![(100, 200), (300, 400)]);
        assert!(fam.txs[1..3]
            .iter()
            .all(|t| t.reads == 7 && t.gene.as_deref() == Some("g") && t.strand == "-"));
        assert_eq!((fam.txs[1].span, fam.txs[2].span), (101, 301));
    }

    #[test]
    fn a_single_exon_unit_joins_the_junction_bearing_transcript_it_overlaps_most() {
        // F's single-exon unit (index 3) overlaps B's two exons by 51 + 51 bases and nothing else
        let txs = vec![
            mk("A", "gL", "c1", "+", 30, &[(100, 200), (300, 400)]),
            mk(
                "B",
                "gR",
                "c1",
                "+",
                30,
                &[(2000, 2100), (2200, 2300), (2400, 2500)],
            ),
            mk("F.U1", "gL", "c1", "+", 10, &[(100, 200), (300, 400)]),
            mk("F.U2", "gL", "c1", "+", 10, &[(2050, 2250)]),
        ];
        let unit = [false, false, true, true];
        assert_eq!(attach_single_exon_units(&txs, &unit), vec![(3, 1)]);
        // an ORIGINAL single-exon transcript is never attached (it is not a unit)
        let t2 = vec![
            mk("B", "gR", "c1", "+", 30, &[(2000, 2100), (2200, 2300)]),
            mk("S", "gS", "c1", "+", 3, &[(2050, 2090)]),
        ];
        assert!(attach_single_exon_units(&t2, &[false, false]).is_empty());
        assert_eq!(
            attach_single_exon_units(&t2, &[false, true]),
            vec![(1, 0)],
            "the same transcript as a unit attaches"
        );
        // the other strand / contig / no overlap: alone
        let t3 = vec![
            mk("B", "gR", "c1", "+", 30, &[(2000, 2100), (2200, 2300)]),
            mk("S", "gS", "c1", "-", 3, &[(2050, 2090)]),
        ];
        assert!(attach_single_exon_units(&t3, &[false, true]).is_empty());
        let t4 = vec![
            mk("B", "gR", "c1", "+", 30, &[(2000, 2100), (2200, 2300)]),
            mk("S", "gS", "c1", "+", 3, &[(2101, 2199)]),
        ];
        assert!(
            attach_single_exon_units(&t4, &[false, true]).is_empty(),
            "an intron is not an exon base"
        );
        // equal overlap with two transcripts: the lower index; more overlap beats a lower index
        let t5 = vec![
            mk("B1", "g1", "c1", "+", 5, &[(2000, 2100), (2200, 2300)]),
            mk("B2", "g2", "c1", "+", 50, &[(2000, 2100), (2500, 2600)]),
            mk("S", "gS", "c1", "+", 3, &[(2000, 2100)]),
        ];
        assert_eq!(
            attach_single_exon_units(&t5, &[false, false, true]),
            vec![(2, 0)]
        );
        let t5b = vec![
            mk("B1", "g1", "c1", "+", 5, &[(2000, 2100), (2200, 2300)]),
            mk("B2", "g2", "c1", "+", 50, &[(2085, 2260), (2500, 2600)]),
            mk("S", "gS", "c1", "+", 3, &[(2090, 2250)]),
        ];
        assert_eq!(
            attach_single_exon_units(&t5b, &[false, false, true]),
            vec![(2, 1)],
            "161 shared bases against 62"
        );
        // it can join a transcript of another gene_id and another unit of the same families input
        let t6 = vec![
            mk("U", "gA", "c1", "+", 5, &[(10, 20), (30, 40)]),
            mk("S", "gB", "c1", "+", 3, &[(35, 38)]),
        ];
        assert_eq!(attach_single_exon_units(&t6, &[true, true]), vec![(1, 0)]);
        // the prefix-max index finds a long exon that starts far left of the unit
        let t7 = vec![
            mk("L", "g", "c1", "+", 5, &[(100, 9000), (9500, 9600)]),
            mk("M", "g", "c1", "+", 5, &[(200, 300), (400, 500)]),
            mk("S", "g", "c1", "+", 3, &[(8000, 8100)]),
        ];
        assert_eq!(
            attach_single_exon_units(&t7, &[false, false, true]),
            vec![(2, 0)]
        );
    }

    /// Names: the best component of an input gene keeps it, the others are `<g>.nat<k>` in representative order, a taken
    /// name is skipped, and the representative is taken AFTER an attached unit joins (a unit carries its parent's reads).
    #[test]
    fn components_are_named_as_the_assembler_names_them() {
        let names = |txs: &[Tx]| {
            native_names(txs, &vec![false; txs.len()])
                .0
                .into_iter()
                .map(|n| n.unwrap())
                .collect::<Vec<_>>()
        };
        let t = [
            mk("A", "g", "c1", "+", 30, &[(100, 200), (300, 400)]),
            mk("B", "g", "c1", "+", 20, &[(2000, 2100), (2200, 2300)]),
        ];
        assert_eq!(names(&t), ["g", "g.nat2"]);
        // a name the input already has is skipped
        let t = [
            mk("A", "g", "c1", "+", 30, &[(100, 200), (300, 400)]),
            mk("B", "g", "c1", "+", 20, &[(2000, 2100), (2200, 2300)]),
            mk("C", "g.nat2", "c1", "+", 5, &[(5000, 5100), (5200, 5300)]),
        ];
        assert_eq!(names(&t), ["g", "g.nat3", "g.nat2"]);
        // different input genes do not compete: each keeps its own
        let t = [
            mk("A", "gL", "c1", "+", 3, &[(100, 200), (300, 400)]),
            mk("B", "gR", "c1", "+", 30, &[(2000, 2100), (2200, 2300)]),
        ];
        assert_eq!(names(&t), ["gL", "gR"]);
        // an attached unit (reads 50) becomes the representative of the component it joins and so decides who keeps `g`
        let txs = [
            mk("A", "g", "c1", "+", 30, &[(100, 200), (300, 400)]),
            mk("B", "g", "c1", "+", 20, &[(2000, 2100), (2200, 2300)]),
            mk("F.U2", "g", "c1", "+", 50, &[(2050, 2090)]),
        ];
        let (n, moved) = native_names(&txs, &[false, false, true]);
        assert_eq!(moved, vec![(2, 1)]);
        let n: Vec<String> = n.into_iter().map(|x| x.unwrap()).collect();
        assert_eq!(
            n,
            ["g.nat2", "g", "g"],
            "B's component (rep F.U2, 50 reads) now beats A's (30)"
        );
    }

    /// An attachment joins a unit to ONE component and never merges two: a single-exon unit over the exons of two components
    /// goes to the one it shares more with (here a tie: the lower index) and the two components stay apart.
    #[test]
    fn an_attachment_never_links_two_components() {
        let txs = [
            mk("A", "g", "c1", "+", 30, &[(100, 200), (300, 400)]),
            mk("B", "g", "c1", "+", 20, &[(1000, 1100), (1200, 1300)]),
            // 51 shared bases with A's second exon (350-400) and 51 with B's first (1000-1050)
            mk("S", "g", "c1", "+", 50, &[(350, 1050)]),
        ];
        let (n, moved) = native_names(&txs, &[false, false, true]);
        assert_eq!(moved, vec![(2, 0)]);
        let n: Vec<String> = n.into_iter().map(|x| x.unwrap()).collect();
        assert_eq!(
            n,
            ["g", "g.nat2", "g"],
            "S joins A's component (its reads make it the representative); B stays apart"
        );
        assert_eq!(
            native_components(&txs).len(),
            3,
            "before the attachment there are three components, after it two: no merge"
        );
    }

    #[test]
    fn a_gene_less_transcript_is_never_regrouped() {
        let mut t = mk("A", "g", "c1", "+", 30, &[(100, 200), (300, 400)]);
        t.gene = None;
        assert_eq!(native_names(&[t], &[false]).0, vec![None]);
    }

    /// X keeps `G` (10 reads), B is the F1 bridge (3 reads), Y (8 reads): the bridge's units join X's and Y's components.
    fn arm_lines(mode: Mode, rb: i64, pas: bool, minus: bool) -> (Vec<String>, Outcome) {
        let mut lines = if minus {
            let m = |ex: &[(i64, i64)]| -> Vec<(i64, i64)> {
                ex.iter()
                    .rev()
                    .map(|&(a, b)| (2001 - b, 2001 - a))
                    .collect()
            };
            let mut v = tx("c1", "G", "X", 10, &m(&[(101, 200), (301, 400)]), "-");
            v.extend(tx(
                "c1",
                "G",
                "B",
                rb,
                &m(&[(101, 200), (301, 350), (1301, 1400), (1501, 1600)]),
                "-",
            ));
            v.extend(tx(
                "c1",
                "G",
                "Y",
                8,
                &m(&[(1101, 1200), (1301, 1400), (1501, 1600)]),
                "-",
            ));
            v
        } else {
            plus_locus(10, rb)
        };
        let out = run_on(&mut lines, mode, locus_evidence(minus, true), pas);
        (lines, out)
    }

    fn line_of<'a>(lines: &'a [String], tid: &str) -> &'a String {
        lines
            .iter()
            .find(|l| l.contains("\ttranscript\t") && attr(l, "transcript_id") == Some(tid))
            .unwrap()
    }

    #[test]
    fn f1units_replaces_the_bridge_by_its_units_in_the_families_input_only() {
        let (f1_gtf, f1) = arm_lines(Mode::F1, 2, true, false);
        let (lines, out) = units_arm();
        // <out>.gtf is what f1 writes, every table too
        assert_eq!(lines, f1_gtf);
        assert_eq!(out.junctions_tsv, f1.junctions_tsv);
        assert!(out.evidence_used && out.bridges_tsv.is_none());
        assert_eq!(out.stats.bridge_transcripts, 1);
        // the families input: X, B.U1, B.U2, Y at B's line position; native names over {X, B.U1}, {B.U2, Y}
        let u = out.units.as_ref().expect("units");
        let fam = &u.families_lines;
        assert_eq!(u.detector, "f1");
        let genes: Vec<(String, String)> = fam
            .iter()
            .filter(|l| l.contains("\ttranscript\t"))
            .map(|l| {
                (
                    attr(l, "transcript_id").unwrap().to_string(),
                    attr(l, "gene_id").unwrap().to_string(),
                )
            })
            .collect();
        let want = |v: &[(&str, &str)]| {
            v.iter()
                .map(|(a, b)| (a.to_string(), b.to_string()))
                .collect::<Vec<_>>()
        };
        assert_eq!(
            genes,
            want(&[
                ("X", "G"),
                ("B.U1", "G"),
                ("B.U2", "G.nat2"),
                ("Y", "G.nat2")
            ])
        );
        assert_eq!(
            line_of(fam, "B.U1"),
            "c1\trustle\ttranscript\t101\t350\t.\t+\t.\tgene_id \"G\"; transcript_id \"B.U1\"; reads \"2\"; TPM \"1.0\"; \
             fusion_of \"B\"; fusion_unit \"1/2\"; fusion_junction \"c1:351-1300:+\"; fusion_locus \"c1:101-1600\"; \
             fusion_gene \"G\"; fusion_detector \"f1\"; fusion_evidence \"2;10;8;0.2\";"
        );
        let i = fam
            .iter()
            .position(|l| l.contains("transcript_id \"B.U2\""))
            .unwrap();
        assert_eq!(fam[i].split('\t').nth(3), Some("1301"));
        assert_eq!(
            fam[i + 1],
            "c1\trustle\texon\t1301\t1400\t.\t+\t.\tgene_id \"G.nat2\"; transcript_id \"B.U2\"; exon_number \"1\";"
        );
        assert_eq!(fam.len(), lines.len() - 5 + 3 + 3, "B's transcript line and 4 exon lines leave, two units of a transcript and 2 exon lines each come");
        assert_eq!(
            out.units_table_rows(),
            vec![
                "B.U1\tB\tG\tG\t1\t2\tc1\t+\t2\t150\t101\t350\t2\tf1\t351-1300\t2;10;8;0.2"
                    .to_string(),
                "B.U2\tB\tG\tG.nat2\t2\t2\tc1\t+\t2\t200\t1301\t1600\t2\tf1\t351-1300\t2;10;8;0.2"
                    .to_string(),
            ]
        );
        let st = &out.stats;
        assert_eq!(
            (
                st.unit_transcripts,
                st.unit_cuts,
                st.units,
                st.unit_single_exon,
                st.unit_attached,
                st.unit_gene_ids_touched
            ),
            (1, 1, 2, 0, 0, 1)
        );
        assert_eq!(
            st.families_gene_ids, 2,
            "G and G.nat2 (F1's own count of the families input is 2 as well)"
        );
    }

    fn units_arm() -> (Vec<String>, Outcome) {
        arm_lines(Mode::F1Units, 2, true, false)
    }

    impl Outcome {
        /// the units table rows without the header (tests)
        fn units_table_rows(&self) -> Vec<String> {
            self.units
                .as_ref()
                .unwrap()
                .table_tsv
                .lines()
                .skip(1)
                .map(str::to_string)
                .collect()
        }
    }

    /// `-`: the same locus mirrored; U1 is the genomic right (the part of the bridge that is Y's), U2 the left.
    #[test]
    fn f1units_on_the_minus_strand_numbers_the_units_in_transcription_order() {
        let (lines, out) = arm_lines(Mode::F1Units, 2, true, true);
        let (f1_gtf, _) = f1_minus();
        assert_eq!(lines, f1_gtf, "<out>.gtf is f1's");
        let fam = out.units.as_ref().unwrap().families_lines.clone();
        // B on `-` (mirror of 101-200, 301-350, 1301-1400, 1501-1600 by 2001 - x): exons 401-500, 601-700, 1651-1700, 1801-1900;
        // the bridge intron is the mirror of 351-1300 = 701-1650
        let u1 = line_of(&fam, "B.U1");
        assert_eq!(u1.split('\t').nth(3), Some("1651"));
        assert!(
            u1.contains("fusion_unit \"1/2\"") && u1.contains("fusion_junction \"c1:701-1650:-\""),
            "{u1}"
        );
        assert_eq!(line_of(&fam, "B.U2").split('\t').nth(4), Some("700"));
        assert_eq!(out.stats.units, 2);
    }

    fn f1_minus() -> (Vec<String>, Outcome) {
        arm_lines(Mode::F1, 2, true, true)
    }

    /// No bridge, no units: the families input is exactly f1's (the scoped form does not touch what no cut touches).
    #[test]
    fn without_a_bridge_the_families_input_is_the_default_arms() {
        for (pas, with_y) in [(false, true), (true, false)] {
            let mk_lines = || plus_locus(10, 2);
            let (mut a, mut b) = (mk_lines(), mk_lines());
            let f1 = run_on(&mut a, Mode::F1, locus_evidence(false, with_y), pas);
            let un = run_on(&mut b, Mode::F1Units, locus_evidence(false, with_y), pas);
            assert_eq!(a, b, "<out>.gtf, pas {pas} own starts {with_y}");
            let u = un.units.unwrap();
            let fam_f1: Vec<String> = a.iter().filter(|l| f1.in_families(l)).cloned().collect();
            assert_eq!(u.families_lines, fam_f1);
            assert_eq!(u.table_tsv, format!("{UNITS_HEADER}\n"));
            assert_eq!(un.stats.units, 0);
        }
        // RG3's names on the untouched gene: two exon-disjoint pieces keep `G` and `G.rg2` whatever the units pass does
        let mut lines = tx("c1", "G", "p", 5, &[(100, 200), (300, 400)], "+");
        lines.extend(tx("c1", "G", "q", 9, &[(5000, 5100), (5200, 5300)], "+"));
        let mut again = lines.clone();
        let f1 = run_on(&mut lines, Mode::F1, BridgeEvidence::default(), true);
        let un = run_on(&mut again, Mode::F1Units, BridgeEvidence::default(), true);
        assert_eq!(gene_of(&lines, "p"), "G.rg2");
        assert_eq!(
            un.units.unwrap().families_lines,
            lines
                .iter()
                .filter(|l| f1.in_families(l))
                .cloned()
                .collect::<Vec<_>>()
        );
    }

    /// The scoped form: a gene_id with no cut keeps RG3's names even when the native rule would name it otherwise
    /// (strand-blind junction sharing, single-exon transcripts), and a gene_id with a cut is named natively.
    #[test]
    fn only_the_gene_ids_that_hold_a_cut_are_regrouped_natively() {
        let mut lines = tx("c1", "G", "X", 10, &[(101, 200), (301, 400)], "+");
        lines.extend(tx(
            "c1",
            "G",
            "B",
            2,
            &[(101, 200), (301, 350), (1301, 1400), (1501, 1600)],
            "+",
        ));
        lines.extend(tx(
            "c1",
            "G",
            "Y",
            8,
            &[(1101, 1200), (1301, 1400), (1501, 1600)],
            "+",
        ));
        // H: untouched; its two strands share a junction (native: one gene; RG3: one piece per strand)
        lines.extend(tx("c1", "H", "h1", 10, &[(5101, 5200), (5301, 5400)], "+"));
        lines.extend(tx("c1", "H", "h2", 11, &[(5101, 5200), (5301, 5400)], "-"));
        let mut ev = locus_evidence(false, true);
        let mut l = lines.clone();
        let out = run(
            &mut l,
            Mode::F1Units,
            &mut ev,
            &genome,
            &fake_clusters(true),
        )
        .unwrap();
        let fam = out.units.unwrap().families_lines;
        let g = |tid: &str| attr(line_of(&fam, tid), "gene_id").unwrap().to_string();
        assert_eq!(
            (g("h2"), g("h1")),
            ("H".to_string(), "H.rg2".to_string()),
            "RG3: one piece per strand, the better keeps H"
        );
        assert_eq!(
            (g("B.U1"), g("B.U2")),
            ("G".to_string(), "G.nat2".to_string())
        );
    }

    // ---- the units list

    #[test]
    fn a_units_list_is_tid_and_junctions_with_a_header() {
        let l = UnitsList::parse(
            "dir/bridges.tsv",
            "# c\ntid\tjunctions\nB\t351-1300\nB\t351-1300,2001-2100\n",
        )
        .unwrap();
        assert_eq!(l.rows.len(), 1);
        assert_eq!(l.rows[0].0, "B");
        let cuts: Vec<(i64, i64)> = l.rows[0].1.iter().map(|c| (c.s, c.e)).collect();
        assert_eq!(
            cuts,
            [(351, 1300), (2001, 2100)],
            "a transcript listed twice gets the union, sorted"
        );
        assert_eq!(
            l.detector(),
            "list:bridges.tsv",
            "the file name, not its directory"
        );
        // both token forms, both separators, other columns ignored
        let l = UnitsList::parse(
            "l",
            "transcript_id\tx\tjunction\nT\t9\tc1:5-9:+;7-8,c2:1-4:-\n",
        )
        .unwrap();
        assert_eq!(
            l.rows[0].1,
            vec![
                ListedCut {
                    s: 1,
                    e: 4,
                    contig: Some("c2".into()),
                    strand: Some("-".into())
                },
                ListedCut {
                    s: 5,
                    e: 9,
                    contig: Some("c1".into()),
                    strand: Some("+".into())
                },
                ListedCut {
                    s: 7,
                    e: 8,
                    contig: None,
                    strand: None
                },
            ]
        );
        for (text, what) in [
            ("", "an empty file"),
            ("tid\tjunctions\n", "no row"),
            ("name\tjunctions\nB\t1-2\n", "no tid column"),
            ("tid\tjunctions\nB\t\n", "a row without a junction"),
            ("tid\tjunctions\nB\t12\n", "a token that is not S-E"),
            ("tid\tjunctions\nB\t9-3\n", "E before S"),
            ("tid\tjunctions\nB\t0-3\n", "S = 0"),
            ("tid\tjunctions\nB\tc1:1-3:x\n", "a bad strand"),
            ("tid\tjunctions\n\t1-3\n", "an empty id"),
            ("tid\tjunctions\nB\n", "a missing column"),
        ] {
            assert!(UnitsList::parse("l", text).is_err(), "{what}");
        }
    }

    #[test]
    fn a_list_names_the_cuts_and_f1s_evidence_is_not_read() {
        let mut lines = plus_locus(10, 2);
        let list = UnitsList::parse(
            "bridges.tsv",
            "tid\tjunctions\nB\tc1:351-1300:+\nNOPE\t1-2\n",
        )
        .unwrap();
        let out = run_list(&mut lines, &list).unwrap();
        // the list names B, so <out>.gtf is f1's regrouping with B as the bridge: same as the evidence arm's
        let (f1_gtf, _) = arm_lines(Mode::F1, 2, true, false);
        assert_eq!(lines, f1_gtf);
        assert!(
            !out.evidence_used && out.junctions_tsv.is_empty(),
            "no junction table without evidence"
        );
        let u = out.units.as_ref().unwrap();
        assert_eq!(u.detector, "list:bridges.tsv");
        let line = line_of(&u.families_lines, "B.U1");
        assert!(
            line.ends_with("fusion_gene \"G\"; fusion_detector \"list:bridges.tsv\";"),
            "no evidence attribute: {line}"
        );
        assert_eq!(
            out.units_table_rows()[0],
            "B.U1\tB\tG\tG\t1\t2\tc1\t+\t2\t150\t101\t350\t2\tlist:bridges.tsv\t351-1300\t."
        );
        assert_eq!(
            (
                out.stats.unit_list_unmatched,
                out.stats.bridge_junctions,
                out.stats.bridge_transcripts
            ),
            (1, 1, 1)
        );
        assert_eq!(out.stats.structural_junctions, 0);
        // the families input is the evidence arm's, except for the detector and the evidence attributes
        let (_, ev_out) = units_arm();
        let strip = |l: &String| {
            let mut s = l.clone();
            for k in ["fusion_detector", "fusion_evidence"] {
                if let Some(i) = s.find(&format!(" {k} \"")) {
                    let j = s[i..].find("\";").unwrap() + i + 2;
                    s.replace_range(i..j, "");
                }
            }
            s
        };
        let a: Vec<String> = u.families_lines.iter().map(strip).collect();
        let b: Vec<String> = ev_out
            .units
            .unwrap()
            .families_lines
            .iter()
            .map(strip)
            .collect();
        assert_eq!(a, b);
    }

    #[test]
    fn a_list_that_does_not_fit_the_gtf_is_an_error_and_writes_nothing() {
        let run_l = |text: &str| {
            let mut lines = plus_locus(10, 2);
            let src = lines.clone();
            let err = run_list(&mut lines, &UnitsList::parse("l.tsv", text).unwrap())
                .err()
                .map(|e| e.to_string());
            if err.is_some() {
                assert_eq!(lines, src, "nothing is rewritten on an error");
            }
            err
        };
        assert!(run_l("tid\tjunctions\nB\t352-1300\n")
            .unwrap()
            .contains("352-1300 is not an intron of that transcript"));
        assert!(run_l("tid\tjunctions\nB\tc9:351-1300:+\n")
            .unwrap()
            .contains("another contig or strand"));
        assert!(run_l("tid\tjunctions\nB\tc1:351-1300:-\n")
            .unwrap()
            .contains("another contig or strand"));
        assert!(run_l("tid\tjunctions\nNOPE\t1-2\n")
            .unwrap()
            .contains("none of its 1 transcripts is in this GTF"));
        assert!(run_l("tid\tjunctions\nB\t351-1300\n").is_none());
    }

    // ---- the prototype

    fn fixture(name: &str) -> String {
        let p = format!(
            "{}/tests/fixtures/bridge_units/{name}",
            env!("CARGO_MANIFEST_DIR")
        );
        std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("read {p}: {e}"))
    }

    /// Remove the attributes the dev prototype (`units2.py`) does not write, to compare with its GTF byte for byte.
    fn strip_extras(line: &str) -> String {
        let mut s = line.to_string();
        for k in [
            "fusion_locus",
            "fusion_gene",
            "fusion_detector",
            "fusion_evidence",
        ] {
            if let Some(i) = s.find(&format!(" {k} \"")) {
                let j = s[i..].find("\";").expect("a closed attribute") + i + 2;
                s.replace_range(i..j, "");
            }
        }
        s
    }

    /// Run `run_list` on a fixture (`plain`, `cuts`) and compare with the prototype's files (`expected` GTF and units table): the
    /// families GTF line for line once the four attributes the prototype does not write are removed, the units table column
    /// for column up to `cuts` (the prototype's `oracle_label` is benchmark-only; ours ends with `evidence`).
    fn assert_equals_prototype(
        plain: &str,
        cuts: &str,
        expected_gtf: &str,
        expected_tsv: &str,
    ) -> (Vec<String>, Outcome) {
        let mut lines: Vec<String> = fixture(plain).lines().map(str::to_string).collect();
        let list = UnitsList::read(&format!(
            "{}/tests/fixtures/bridge_units/{cuts}",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let out = run_list(&mut lines, &list).unwrap();
        let u = out.units.as_ref().unwrap();
        assert_eq!(u.detector, format!("list:{cuts}"));
        let got: Vec<String> = u.families_lines.iter().map(|l| strip_extras(l)).collect();
        let want: Vec<String> = fixture(expected_gtf).lines().map(str::to_string).collect();
        assert_eq!(got.len(), want.len());
        for (i, (a, b)) in got.iter().zip(&want).enumerate() {
            assert_eq!(a, b, "{expected_gtf} line {}", i + 1);
        }
        let cols = |t: &str| -> Vec<Vec<String>> {
            t.lines()
                .map(|l| l.split('\t').take(15).map(str::to_string).collect())
                .collect()
        };
        let (mine, theirs) = (cols(&u.table_tsv), cols(&fixture(expected_tsv)));
        assert_eq!(mine.len(), theirs.len());
        for (a, b) in mine.iter().zip(&theirs) {
            // the `label` column holds the detector: `list:<file>` here, the `--label` given to the prototype there
            assert_eq!([&a[..13], &a[14..]].concat(), [&b[..13], &b[14..]].concat());
        }
        (lines, out)
    }

    /// The Rust units execution equals the prototype's (`units2.py --regroup scoped --attach units`, 9fea69a1) on the
    /// synthetic fixture of `tests/fixtures/bridge_units` (m = 1, 2, 3; both strands; single-exon attachments across genes
    /// and on a tie; `.nat<k>`, `.rg<k>` and RG3 names; contigs): gene_id for gene_id, line for line, and its units table.
    #[test]
    fn units_equal_the_python_prototype_on_the_fixture() {
        let (lines, out) = assert_equals_prototype(
            "plain.gtf",
            "cuts.tsv",
            "expected.units.gtf",
            "expected.units.tsv",
        );
        let st = &out.stats;
        assert_eq!(
            (
                st.unit_transcripts,
                st.unit_cuts,
                st.units,
                st.unit_single_exon,
                st.unit_attached
            ),
            (5, 8, 13, 2, 2)
        );
        // the pre-split locus of each cut transcript is carried for the relation records
        let f3 = line_of(&out.units.as_ref().unwrap().families_lines, "F3.U2");
        assert!(
            f3.contains("fusion_locus \"c1:20000-24100\"; fusion_gene \"g3\";"),
            "{f3}"
        );
        // the assembled GTF keeps every transcript whole, as f1 writes it
        assert_eq!(
            lines.iter().filter(|l| l.contains("fusion_of \"")).count(),
            5
        );
        assert!(lines.iter().all(|l| !l.contains("fusion_unit")));
    }

    /// The same on real gorilla loci: 34 transcripts of 6 gene_ids of the simulation's f = .5 assembly (the four genes holding
    /// an F1 bridge, on both strands, and two others) with the four bridges F1 found (`R.cuts.tsv` of the units study).
    #[test]
    fn units_equal_the_python_prototype_on_real_gorilla_loci() {
        let (_, out) = assert_equals_prototype(
            "s_f0.5.plain.gtf",
            "s_f0.5.cuts.tsv",
            "s_f0.5.expected.units.gtf",
            "s_f0.5.expected.units.tsv",
        );
        let st = &out.stats;
        assert_eq!(
            (
                st.unit_transcripts,
                st.unit_cuts,
                st.units,
                st.unit_gene_ids_touched
            ),
            (4, 4, 8, 4)
        );
    }

    /// The scoped form's invariant: a single-exon unit attached to a transcript of an UNTOUCHED gene lands in that
    /// transcript's gene_id, also when the target sits in the gene's non-keeper RG3 piece (`H.rg2`), where the unit's own
    /// native component is named otherwise (`G.nat2`). Without the override (`names[i] = names[j]`) the unit would sit in
    /// another locus than the transcript it was attached to.
    #[test]
    fn an_attached_unit_lands_in_its_targets_gene_id() {
        let mut lines = tx("c1", "G", "A", 20, &[(100, 200), (300, 400)], "+");
        lines.extend(tx(
            "c1",
            "G",
            "F",
            7,
            &[(100, 200), (300, 400), (5050, 5090)],
            "+",
        ));
        lines.extend(tx("c1", "H", "H1", 10, &[(1000, 1100), (1200, 1300)], "+"));
        lines.extend(tx("c1", "H", "H2", 5, &[(5000, 5100), (5200, 5300)], "+"));
        let list = UnitsList::parse("l.tsv", "tid\tjunctions\nF\t401-5049\n").unwrap();
        let out = run_list(&mut lines, &list).unwrap();
        let fam = &out.units.as_ref().unwrap().families_lines;
        let gene = |t: &str| attr(line_of(fam, t), "gene_id").unwrap().to_string();
        assert_eq!(
            (gene("H1"), gene("H2")),
            ("H".to_string(), "H.rg2".to_string()),
            "RG3's pieces of the untouched gene H"
        );
        assert_eq!(
            (gene("A"), gene("F.U1")),
            ("G".to_string(), "G".to_string()),
            "the other unit is natively in A's component"
        );
        assert_eq!(out.stats.unit_attached, 1);
        assert_eq!(
            gene("F.U2"),
            "H.rg2",
            "the single-exon unit is in its target's gene_id, not in `G.nat2`"
        );
        // and so it is a member of H2's locus in the units table, whose `new_gene` column says the same
        assert!(
            out.units_table_rows()[1].starts_with("F.U2\tF\tG\tH.rg2\t2\t2\t"),
            "{}",
            out.units_table_rows()[1]
        );
    }

    /// One bridge transcript with TWO bridge introns (X - J1 - Y - J2 - Z): three units in transcription order, every unit line
    /// carries both cuts (genomic order) and both evidence strings comma-joined in cut order, the two single-exon units
    /// (the middle one and the last) are attached to Y's and Z's components, and the table's `cuts` / `evidence` columns say
    /// the same (the relations' 18th column is this string).
    #[test]
    fn two_cuts_of_one_transcript_join_their_junctions_and_evidence() {
        let mut lines = tx("c1", "G", "X", 10, &[(101, 200), (301, 400)], "+");
        lines.extend(tx(
            "c1",
            "G",
            "T",
            2,
            &[(101, 200), (301, 350), (1301, 1350), (2301, 2400)],
            "+",
        ));
        lines.extend(tx("c1", "G", "Y", 10, &[(1101, 1200), (1301, 1400)], "+"));
        lines.extend(tx("c1", "G", "Z", 10, &[(2101, 2200), (2301, 2400)], "+"));
        let mut ev = BridgeEvidence::default();
        for end in [398, 399, 400, 400] {
            push(&mut ev, &[(101, 200), (301, end)], false); // X: 3' ends inside J1 = 351-1300
        }
        for end in [1398, 1399, 1400] {
            push(&mut ev, &[(1101, 1200), (1301, end)], false); // Y: starts inside J1, 3' ends inside J2 = 1351-2300
        }
        for _ in 0..3 {
            push(&mut ev, &[(2101, 2200), (2301, 2400)], false); // Z: starts inside J2
        }
        let out = run_on(&mut lines, Mode::F1Units, ev, true);
        assert_eq!(
            (out.stats.f1_bridge_junctions, out.stats.bridge_transcripts),
            (2, 1),
            "T uses both bridge introns"
        );
        let fam = &out.units.as_ref().unwrap().families_lines;
        let evidence = "2;10;20;0.1667,2;20;10;0.1667"; // J1: link 2 reads, UP = X 10, DOWN = Y + Z 20; J2: UP = X + Y 20, DOWN = Z 10
        for (k, (tid, span)) in [
            ("T.U1", (101, 350)),
            ("T.U2", (1301, 1350)),
            ("T.U3", (2301, 2400)),
        ]
        .iter()
        .enumerate()
        {
            let l = line_of(fam, tid);
            assert!(l.contains(&format!("fusion_unit \"{}/3\"", k + 1)), "{l}");
            assert!(
                l.contains("fusion_junction \"c1:351-1300:+,c1:1351-2300:+\";"),
                "{l}"
            );
            assert!(
                l.ends_with(&format!("fusion_evidence \"{evidence}\";")),
                "{l}"
            );
            assert_eq!(
                l.split('\t').nth(3).and_then(|x| x.parse::<i64>().ok()),
                Some(span.0),
                "{l}"
            );
        }
        let genes: Vec<String> = ["X", "T.U1", "T.U2", "Y", "T.U3", "Z"]
            .iter()
            .map(|t| attr(line_of(fam, t), "gene_id").unwrap().to_string())
            .collect();
        assert_eq!(
            genes,
            ["G", "G", "G.nat2", "G.nat2", "G.nat3", "G.nat3"],
            "U2 joins Y's component, U3 Z's, U1 X's"
        );
        let st = &out.stats;
        assert_eq!(
            (
                st.unit_transcripts,
                st.unit_cuts,
                st.units,
                st.unit_single_exon,
                st.unit_attached
            ),
            (1, 2, 3, 2, 2)
        );
        let rows = out.units_table_rows();
        assert!(
            rows[1].ends_with(&format!("\t2\tf1\t351-1300;1351-2300\t{evidence}")),
            "{}",
            rows[1]
        );
    }

    /// An intron listed twice (as `S-E`, as `CONTIG:S-E:STRAND`, on two rows of one transcript) is cut once: one cut, two
    /// units, one junction in `fusion_junction`.
    #[test]
    fn a_cut_listed_twice_is_cut_once() {
        let mut lines = plus_locus(10, 2);
        let list = UnitsList::parse(
            "l.tsv",
            "tid\tjunctions\nB\t351-1300,c1:351-1300:+;351-1300\nB\t351-1300\n",
        )
        .unwrap();
        let out = run_list(&mut lines, &list).unwrap();
        assert_eq!(
            (
                out.stats.unit_cuts,
                out.stats.units,
                out.stats.bridge_junctions
            ),
            (1, 2, 1)
        );
        let u = line_of(&out.units.as_ref().unwrap().families_lines, "B.U1");
        assert!(u.contains("fusion_junction \"c1:351-1300:+\"; "), "{u}");
        assert!(
            out.units_table_rows()[0].contains("\t351-1300\t"),
            "{}",
            out.units_table_rows()[0]
        );
    }

    /// A list with a header and no transcript (or no header at all) is an ERROR, not "cut nothing": a list written for a
    /// sample with no annotated fusion would otherwise run as the default arm under another name. The message names the way out.
    #[test]
    fn a_header_only_list_is_an_error_that_says_what_to_do() {
        for text in [
            "tid\tjunctions\n",
            "# only a comment\ntid\tjunctions\n\n",
            "",
        ] {
            let e = UnitsList::parse("l.tsv", text).unwrap_err().to_string();
            assert!(
                e.contains("empty list: run without --bridge-units-list"),
                "{text:?}: {e}"
            );
        }
    }
}
