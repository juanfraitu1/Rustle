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
use crate::vg_family::family_container::{gtf_attr as attr, key_str, parse_key, py_int, Key};

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
