//! One family, start to end (o3port.py `run_family`, `isocon_side`, `isocon_keep`): piles -> clusters -> haplotypes -> merged
//! haplotypes -> DADA2 explanation -> allele slots -> certified copies, putative copies, and the assignment of tied reads.

use super::cluster::{absorb, cluster_pile, conflicting_set, haplotype, Cluster};
use super::pile::{blocks, cmp_base, observations, Piles, Prim, RecSummary};
use super::putative::{foreign_clipping, foreign_coverage, foreign_splicing, junction_usage, putative_copies, Usage};
use super::pyjson::Json;
use super::slots::{b_matching, consensus, divergence, explanation, forced_edges, pile_consensus, private_alleles, private_detail, Explain, Matching, Merged, Options, Parent};
use super::tied_reads::{disfavoured, homology, one_site, semi_tied, SiteCache};
use super::{Copy, Hap, Params, ReadId, SlotOrder, MIN_READS, MIN_SUPPORT};
use crate::genome::IndexedFasta;
use anyhow::{Context, Result};
use indexmap::IndexMap;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

/// What `run_family` returns (the Python's tuple), plus the family's dump object when asked for.
pub struct FamilyOut {
    /// (home copy, reads, contig) of each certified haplotype.
    pub certified: Vec<(String, BTreeSet<ReadId>, Vec<u8>)>,
    /// {copy: (pile reads, clusters, conflicting)} of the piles phased.
    pub per_copy: IndexMap<String, (usize, usize, usize)>,
    pub n_h: usize,
    /// Counter of the haplotypes' fates, in first-seen order.
    pub why: Vec<(&'static str, usize)>,
    pub piled: BTreeSet<ReadId>,
    /// (slot copy, reads, contig, multi fraction, the other holders' multi fractions) of each putative copy.
    pub looks: Vec<(String, BTreeSet<ReadId>, Vec<u8>, f64, Vec<f64>)>,
    /// (read, assigned, piles holding it) of every tied read.
    pub o2: Vec<(ReadId, String, usize)>,
    /// The dump object (family, copies, clusters, merged), when asked for.
    pub dump: Option<Vec<(String, Json)>>,
}

fn upper_base(genome: &IndexedFasta, chrom: &str, pos: i64) -> u8 {
    genome.fetch(chrom, pos, pos + 1).and_then(|v| v.first().copied()).map_or(b'N', |b| b.to_ascii_uppercase())
}

fn upper(v: Option<Vec<u8>>) -> Vec<u8> {
    let mut v = v.unwrap_or_default();
    v.make_ascii_uppercase();
    v
}

fn bump(why: &mut Vec<(&'static str, usize)>, k: &'static str) {
    match why.iter_mut().find(|(w, _)| *w == k) {
        Some(e) => e.1 += 1,
        None => why.push((k, 1)),
    }
}

fn b2s(b: u8) -> Json {
    Json::Str((b as char).to_string())
}

/// The haplotype's fate: (fate, by).
enum By {
    None,
    Copy(String),
    Join(String, String),
}

/// `certificate_matching` (`--certify-all-records`, register 1414): the slot matching a certificate must survive. A haplotype
/// whose reads' primaries mostly lie on no copy of the family is explained by that reference locus and claims no slot; the others
/// may take a slot on every family copy most of their reads have ANY record on, not only the copies inside the tie window.
#[allow(clippy::too_many_arguments)]
fn certificate_matching(
    demand: &[usize],
    merged: &[Merged],
    options: &Options,
    piles: &Piles,
    copies: &[Copy],
    records: &HashMap<ReadId, Vec<RecSummary>>,
    capacity: &BTreeMap<String, usize>,
    order: &[usize],
) -> Matching {
    let on = |ch: &str, s0: i64, e0: i64, c: &Copy| ch == c.chrom && s0 < c.end && e0 > c.start;
    let mut sorted_copies: Vec<&Copy> = copies.iter().collect();
    sorted_copies.sort_by(|a, b| a.tid.cmp(&b.tid));
    let none = Vec::new();
    let recs = |r: &ReadId| records.get(r).unwrap_or(&none);
    let mut copts = Options::new();
    for &i in demand {
        let reads = &merged[i].reads;
        let out = reads
            .iter()
            .flat_map(|r| recs(r).iter())
            .filter(|(ch, s0, e0, _, sec, _)| !*sec && !copies.iter().any(|c| on(ch, *s0, *e0, c)))
            .count();
        if out * 2 > reads.len() {
            continue;
        }
        let mut o: BTreeSet<String> = options[&i].iter().cloned().collect();
        for c in &sorted_copies {
            let k = reads.iter().filter(|r| recs(r).iter().any(|(ch, s0, e0, ..)| on(ch, *s0, *e0, c))).count();
            if k > 0 && piles.contains_key(&c.tid) && k * 2 >= reads.len() {
                o.insert(c.tid.clone());
            }
        }
        copts.insert(i, o.into_iter().collect());
    }
    let ord: Vec<usize> = order.iter().copied().filter(|i| copts.contains_key(i)).collect();
    b_matching(&copts, capacity, Some(&ord))
}

/// `run_family` on the family's piles (`pile_records`, built by the caller). `names` maps read ids to names for the dump.
#[allow(clippy::too_many_arguments)]
pub fn run_family(
    fam: &str,
    copies: &[Copy],
    piles: &Piles,
    prim: &HashMap<ReadId, Prim>,
    genome: &IndexedFasta,
    p: &Params,
    tmpdir: &Path,
    names: Option<&[String]>,
    records: &HashMap<ReadId, Vec<RecSummary>>,
) -> Result<FamilyOut> {
    let span: HashMap<&str, (i64, i64)> = copies.iter().map(|c| (c.tid.as_str(), (c.start, c.end))).collect();
    let chrom: HashMap<&str, &str> = copies.iter().map(|c| (c.tid.as_str(), c.chrom.as_str())).collect();
    let name = |r: ReadId| Json::Str(names.map_or(String::new(), |n| n[r as usize].clone()));
    let read_list = |rs: &BTreeSet<ReadId>| Json::List(rs.iter().map(|&r| name(r)).collect());
    let mut haps: Vec<(String, Cluster, f64, Hap)> = Vec::new();
    let mut per_copy: IndexMap<String, (usize, usize, usize)> = IndexMap::new();
    let mut allhaps: Vec<(String, usize, Hap)> = Vec::new();
    let mut cons: HashMap<String, IndexMap<i64, u8>> = HashMap::new();
    let mut d_copies: Vec<(String, Json)> = Vec::new();
    let mut d_clusters: Vec<Json> = Vec::new();
    let dumping = names.is_some();
    if dumping {
        for (c, pl) in piles {
            d_copies.push((c.clone(), Json::Obj(vec![("pile".into(), read_list(&pl.keys().copied().collect()))])));
        }
    }
    for (cid, pile) in piles {
        if pile.len() < 2 * MIN_READS {
            continue;
        }
        let (obs, _) = observations(pile.values(), prim, span[cid.as_str()]);
        let reads: BTreeSet<ReadId> = pile.keys().copied().collect();
        let (mut clusters, top) = cluster_pile(&reads, &obs, p.clusterer, p.correct);
        let mut hs: Vec<Hap> = clusters.iter().map(|c| haplotype(c, &obs, &top)).collect();
        let keep: Vec<usize> = if p.floor_first {
            let big: Vec<usize> = (0..clusters.len()).filter(|&i| clusters[i].len() >= MIN_SUPPORT).collect();
            if big.is_empty() {
                Vec::new()
            } else {
                let bh: Vec<Hap> = big.iter().map(|&i| hs[i].clone()).collect();
                conflicting_set(&bh).into_iter().map(|k| big[k]).collect()
            }
        } else if hs.is_empty() {
            Vec::new()
        } else {
            conflicting_set(&hs)
        };
        let joined = if p.absorb { absorb(&hs, &keep) } else { BTreeMap::new() };
        for (&j, &k) in &joined {
            let add = clusters[j].clone();
            clusters[k].extend(add);
        }
        for &k in joined.values().collect::<BTreeSet<_>>() {
            hs[k] = haplotype(&clusters[k], &obs, &top);
        }
        per_copy.insert(cid.clone(), (pile.len(), clusters.len(), keep.len()));
        let ch = chrom[cid.as_str()];
        let rf: BTreeMap<i64, u8> = top.keys().map(|&col| (col, upper_base(genome, ch, col))).collect();
        cons.insert(cid.clone(), pile_consensus(&obs));
        for (i, c) in clusters.iter().enumerate() {
            allhaps.push((cid.clone(), c.len(), hs[i].clone()));
        }
        for &i in &keep {
            haps.push((cid.clone(), clusters[i].clone(), divergence(&hs[i], &rf), hs[i].clone()));
        }
        if dumping {
            if let Some((_, Json::Obj(kv))) = d_copies.iter_mut().find(|(c, _)| c == cid) {
                kv.push(("columns".into(), Json::Int(top.len() as i64)));
            }
            for (i, c) in clusters.iter().enumerate() {
                d_clusters.push(Json::Obj(vec![
                    ("copy".into(), Json::str(cid.clone())),
                    ("reads".into(), read_list(c)),
                    ("hap".into(), Json::Obj(hs[i].iter().map(|(k, &v)| (k.to_string(), b2s(v))).collect())),
                    ("ref".into(), Json::Obj(hs[i].keys().map(|k| (k.to_string(), b2s(rf[k]))).collect())),
                    ("kept".into(), Json::Bool(keep.contains(&i))),
                    ("joined".into(), Json::Bool(joined.contains_key(&i))),
                ]));
            }
        }
    }
    // one haplotype found in several piles (tied reads) is one haplotype: merge by read overlap (>= half of the smaller)
    haps.sort_by(|a, b| b.1.len().cmp(&a.1.len()).then_with(|| a.0.cmp(&b.0)));
    let mut merged: Vec<Merged> = Vec::new();
    for (cid, reads_h, dv, hp) in haps {
        let hit = merged.iter_mut().find(|m| 2 * m.reads.intersection(&reads_h).count() >= m.reads.len().min(reads_h.len()));
        match hit {
            Some(m) => {
                if p.absorb {
                    m.reads.extend(reads_h.iter().copied());
                }
                m.piles.insert(cid.clone());
                m.div = m.div.min(dv);
                m.hap.entry(cid).or_insert(hp);
            }
            None => merged.push(Merged { reads: reads_h, piles: BTreeSet::from([cid.clone()]), div: dv, hap: IndexMap::from([(cid, hp)]) }),
        }
    }
    let in_pile = |m: &Merged, c: &str| m.reads.iter().filter(|r| piles[c].contains_key(*r)).count();
    let mut options: Options = Options::new();
    for (i, m) in merged.iter().enumerate() {
        // sorted, as the Python: b_matching takes the first option with a free slot, so the order picks the slot copy
        let mut opt: Vec<String> = piles.keys().filter(|c| 2 * in_pile(m, c) >= m.reads.len()).cloned().collect();
        opt.sort();
        options.insert(i, opt);
    }
    // DADA2's rule before the slots
    let mut maps_of: HashMap<String, HashMap<String, HashMap<i64, (i64, bool)>>> = HashMap::new();
    let mut why: Vec<(&'static str, usize)> = Vec::new();
    let mut fate: BTreeMap<usize, (&'static str, By)> = BTreeMap::new();
    let mut demand: Vec<usize> = Vec::new();
    for (i, m) in merged.iter().enumerate() {
        if m.reads.len() < MIN_SUPPORT {
            bump(&mut why, "below the floor");
            fate.insert(i, ("below the floor", By::None));
            continue;
        }
        let home = m.hap.keys().max_by(|a, b| (in_pile(m, a), a.as_str()).cmp(&(in_pile(m, b), b.as_str()))).unwrap().clone();
        let hp = &m.hap[&home];
        if !maps_of.contains_key(&home) {
            let hc = copies.iter().find(|c| c.tid == home).unwrap();
            let others: Vec<&Copy> = copies.iter().filter(|c| c.tid != home).collect();
            maps_of.insert(home.clone(), homology(hc, &others, genome, tmpdir)?);
        }
        let mp = &maps_of[&home];
        let hch = chrom[home.as_str()];
        let mut parents: Vec<Parent> = vec![(format!("ref:{home}"), hp.keys().map(|&c| (c, Some(upper_base(genome, hch, c)))).collect())];
        for c in copies {
            if c.tid == home || !mp.contains_key(&c.tid) {
                continue;
            }
            let mut pb = BTreeMap::new();
            for &col in hp.keys() {
                if let Some(&(z, plus)) = mp[&c.tid].get(&col) {
                    let b = upper_base(genome, &c.chrom, z);
                    pb.insert(col, Some(if plus { b } else { cmp_base(b) }));
                }
            }
            parents.push((format!("ref:{}", c.tid), pb));
        }
        for (k, (pid, n, g)) in allhaps.iter().enumerate() {
            if *n <= m.reads.len() {
                continue;
            }
            let mut pb = BTreeMap::new();
            for &col in hp.keys() {
                if *pid == home {
                    pb.insert(col, g.get(&col).copied().or_else(|| cons[pid].get(&col).copied()));
                } else if let Some(&(z, plus)) = mp.get(pid).and_then(|x| x.get(&col)) {
                    let b = g.get(&z).copied().or_else(|| cons[pid].get(&z).copied()).unwrap_or_else(|| upper_base(genome, chrom[pid.as_str()], z));
                    pb.insert(col, Some(if plus { b } else { cmp_base(b) }));
                }
            }
            parents.push((format!("hap{k}"), pb));
        }
        match explanation(hp, &parents) {
            None => {
                demand.push(i);
                bump(&mut why, "unexplained");
                fate.insert(i, ("unexplained", By::None));
            }
            Some(Explain::Single(n)) if n.starts_with("ref:") => {
                let z = n[4..].to_string();
                let o = options.get_mut(&i).unwrap();
                if !o.contains(&z) && piles.contains_key(&z) {
                    o.push(z.clone());
                    o.sort();
                }
                demand.push(i);
                bump(&mut why, "equals a reference copy");
                fate.insert(i, ("equals a reference copy", By::Copy(z)));
            }
            Some(Explain::Single(n)) => {
                bump(&mut why, "duplicate");
                fate.insert(i, ("duplicate", By::Copy(n)));
            }
            Some(Explain::Join(a, b)) => {
                bump(&mut why, "chimera");
                fate.insert(i, ("chimera", By::Join(a, b)));
            }
        }
    }
    if p.prefer {
        let mut pc = SiteCache::default();
        for &i in &demand {
            let drop = disfavoured(&merged[i].reads, &options[&i], piles, prim, genome, &mut pc);
            if !drop.is_empty() && drop.len() < options[&i].len() {
                options.get_mut(&i).unwrap().retain(|c| !drop.contains(c));
            }
        }
    }
    let capacity: BTreeMap<String, usize> = piles.keys().map(|c| (c.clone(), 2)).collect();
    let mut order = demand.clone();
    let n_reads = |i: usize| -(merged[i].reads.len() as i64);
    match p.slot_order {
        SlotOrder::Support => order.sort_by(|&a, &b| (n_reads(a), merged[a].div, a).partial_cmp(&(n_reads(b), merged[b].div, b)).unwrap()),
        SlotOrder::Divergence => order.sort_by(|&a, &b| (merged[a].div, n_reads(a), a).partial_cmp(&(merged[b].div, n_reads(b), b)).unwrap()),
    }
    let dopts: Options = demand.iter().map(|i| (*i, options[i].clone())).collect();
    let matched: Matching = b_matching(&dopts, &capacity, Some(&order));
    let forced = forced_edges(&dopts, &capacity, &matched, &order);
    let mut n_piles: IndexMap<ReadId, usize> = IndexMap::new();
    for pl in piles.values() {
        for &r in pl.keys() {
            *n_piles.entry(r).or_default() += 1;
        }
    }
    let usage: HashMap<&str, Usage> = per_copy.keys().map(|c| (c.as_str(), junction_usage(&piles[c]))).collect();
    let blk: HashMap<&str, HashMap<ReadId, Vec<(i64, i64)>>> =
        per_copy.keys().map(|c| (c.as_str(), piles[c].iter().map(|(&n, rec)| (n, blocks(rec))).collect())).collect();
    let foreign = |i: usize, x: &str, co: &[usize]| -> bool {
        let Some(u) = usage.get(x) else { return false };
        let ch = chrom[x];
        let site = |j: (i64, i64)| (upper(genome.fetch(ch, j.0, j.0 + 2)), upper(genome.fetch(ch, j.1 - 2, j.1)));
        foreign_splicing(&merged[i].reads, u, site)
            || (p.clipping && foreign_clipping(&merged[i].reads, &piles[x]))
            || (p.coverage && {
                let co_reads: BTreeSet<ReadId> = co.iter().flat_map(|j| merged[*j].reads.iter().copied()).collect();
                foreign_coverage(&merged[i].reads, &co_reads, &blk[x])
            })
    };
    let reads_of: BTreeMap<usize, BTreeSet<ReadId>> = demand.iter().map(|&i| (i, merged[i].reads.clone())).collect();
    let n_piles_map: HashMap<ReadId, usize> = n_piles.iter().map(|(&r, &k)| (r, k)).collect();
    let mut looks = Vec::new();
    let mut look_of: Vec<usize> = Vec::new();
    let fref: &dyn Fn(usize, &str, &[usize]) -> bool = &foreign;
    for (i, x, mi, mc) in putative_copies(&demand, &matched, &reads_of, &n_piles_map, if p.structure { Some(fref) } else { None }) {
        if merged[i].reads.len() < MIN_SUPPORT {
            continue;
        }
        let (_, steps) = observations(merged[i].reads.iter().filter_map(|r| piles[&x].get(r)), prim, span[x.as_str()]);
        looks.push((x, merged[i].reads.clone(), consensus(steps.values()), mi, mc));
        look_of.push(i);
    }
    let cmatch = p.certify_all.then(|| certificate_matching(&demand, &merged, &options, piles, copies, records, &capacity, &order));
    let mut certified = Vec::new();
    let mut cert_of: Vec<usize> = Vec::new();
    for &i in &demand {
        let m = &merged[i];
        let survives = cmatch.as_ref().is_none_or(|cm| matches!(cm.get(&i), Some(None)));
        if matched[&i].is_none() && private_alleles(i, &merged, &matched, &options) >= MIN_READS && survives {
            let home = m.piles.iter().max_by(|a, b| (in_pile(m, a), a.as_str()).cmp(&(in_pile(m, b), b.as_str()))).unwrap().clone();
            let (_, steps) = observations(m.reads.iter().filter_map(|r| piles[&home].get(r)), prim, span[home.as_str()]);
            certified.push((home, m.reads.clone(), consensus(steps.values())));
            cert_of.push(i);
        }
    }
    let mut d_merged = Vec::new();
    if dumping {
        for (i, m) in merged.iter().enumerate() {
            let (f, by) = &fate[&i];
            let by = match by {
                By::None => Json::Null,
                By::Copy(c) => Json::str(c.clone()),
                By::Join(a, b) => Json::List(vec![Json::str(a.clone()), Json::str(b.clone())]),
            };
            let in_demand = demand.contains(&i);
            let slot = if in_demand { matched[&i].clone().map_or(Json::Null, Json::Str) } else { Json::Null };
            let mut row = vec![
                ("reads".to_string(), read_list(&m.reads)),
                ("piles".into(), Json::List(m.piles.iter().map(|c| Json::str(c.clone())).collect())),
                ("div".into(), Json::Float(m.div)),
                ("options".into(), Json::List(options[&i].iter().map(|c| Json::str(c.clone())).collect())),
                ("fate".into(), Json::str(*f)),
                ("by".into(), by),
                ("slot".into(), slot),
                ("certified".into(), Json::Bool(cert_of.contains(&i))),
                ("putative".into(), Json::Bool(look_of.contains(&i))),
            ];
            if in_demand && matched[&i].is_none() {
                let (cid, pv) = private_detail(i, &merged, &matched, &options);
                let cid = cid.unwrap_or_default();
                let ch = chrom.get(cid.as_str()).copied().unwrap_or("");
                row.push((
                    "private".into(),
                    Json::List(pv.iter().map(|&(col, b)| Json::List(vec![Json::str(cid.clone()), Json::Int(col), b2s(b), b2s(upper_base(genome, ch, col))])).collect()),
                ));
            }
            d_merged.push(Json::Obj(row));
        }
    }
    let piled: BTreeSet<ReadId> = per_copy.keys().flat_map(|c| piles[c].keys().copied()).collect();
    // each AS-tied read (in >= 2 piles) goes to its haplotype's forced copy, or to the certified copy it belongs to; else abstain
    let mut home_of: HashMap<ReadId, usize> = HashMap::new();
    for &i in &demand {
        for &r in &merged[i].reads {
            home_of.entry(r).or_insert(i);
        }
    }
    let mut cert_reads: HashMap<ReadId, usize> = HashMap::new();
    if p.reconstruct {
        for (k, (_, rs, _)) in certified.iter().enumerate() {
            for &r in rs {
                cert_reads.insert(r, k);
            }
        }
    }
    let mut unsure: BTreeSet<String> = looks.iter().map(|l| l.0.clone()).collect();
    for &i in &look_of {
        unsure.extend(options[&i].iter().cloned());
    }
    let (mut o2, mut semi_cache, mut hap_cache) = (Vec::new(), SiteCache::default(), HashMap::new());
    for (&r, &k) in &n_piles {
        if k < 2 {
            continue;
        }
        let on: BTreeSet<&String> = piles.iter().filter(|(_, pl)| pl.contains_key(&r)).map(|(c, _)| c).collect();
        let a = if let Some(c) = cert_reads.get(&r) {
            format!("new:{c}")
        } else if p.guard && on.iter().any(|c| unsure.contains(*c)) {
            "ABSTAIN".into()
        } else if let Some(x) = home_of.get(&r).and_then(|h| forced.get(h).filter(|_| on.iter().all(|c| options[h].contains(c)))) {
            x.clone()
        } else if p.semi {
            let mut a = semi_tied(r, piles, prim, genome, &mut semi_cache, p);
            if a == "ABSTAIN" && p.one_site {
                if let Some(&h) = home_of.get(&r) {
                    a = one_site(r, &merged[h].reads, piles, prim, genome, &mut semi_cache, &mut hap_cache, h);
                }
            }
            a
        } else {
            "ABSTAIN".into()
        };
        o2.push((r, a, k));
    }
    let dump = dumping.then(|| {
        vec![
            ("family".to_string(), Json::str(fam)),
            ("copies".into(), Json::Obj(d_copies)),
            ("clusters".into(), Json::List(d_clusters)),
            ("merged".into(), Json::List(d_merged)),
        ]
    });
    Ok(FamilyOut { certified, per_copy, n_h: merged.len(), why, piled, looks, o2, dump })
}

/// One flagged candidate of the shipped stage: (candidate, member read names, contig, members that joined by attribution).
pub type IsoconCandidate = (String, BTreeSet<String>, Vec<u8>, usize);

/// `isocon_side`: {family: [flagged candidate of the stage]} from the stage's products at `prefix` and the net file.
pub fn isocon_side(prefix: &str, net_path: &str) -> Result<BTreeMap<String, Vec<IsoconCandidate>>> {
    let tsv = |path: String| -> Result<Vec<HashMap<String, String>>> {
        let text = std::fs::read_to_string(&path).with_context(|| format!("reading {path}"))?;
        let mut lines = text.lines();
        let head: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        Ok(lines.filter(|l| !l.is_empty()).map(|l| head.iter().zip(l.split('\t')).map(|(h, v)| (h.to_string(), v.to_string())).collect()).collect())
    };
    let via: HashMap<(String, String), String> =
        tsv(net_path.to_string())?.into_iter().map(|r| ((r["family"].clone(), r["read"].clone()), r["via"].clone())).collect();
    let cand: BTreeMap<String, HashMap<String, String>> = tsv(format!("{prefix}.candidates.tsv"))?.into_iter().map(|r| (r["candidate"].clone(), r)).collect();
    let owner: HashMap<String, String> = tsv(format!("{prefix}.clusters.tsv"))?
        .into_iter()
        .filter(|r| cand.contains_key(&r["candidate"]))
        .map(|r| (r["cluster"].clone(), r["candidate"].clone()))
        .collect();
    let mut members: HashMap<String, BTreeSet<(String, String)>> = HashMap::new();
    for r in tsv(format!("{prefix}.reads.tsv"))? {
        if let Some(c) = owner.get(&r["cluster"]) {
            members.entry(c.clone()).or_default().insert((r["family"].clone(), r["read"].clone()));
        }
    }
    let path = format!("{prefix}.contigs.fa");
    let text = std::fs::read_to_string(&path).with_context(|| format!("reading {path}"))?;
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut cur: Option<String> = None;
    for ln in text.lines() {
        let ln = ln.trim();
        if let Some(h) = ln.strip_prefix('>') {
            let n = h.split_whitespace().next().unwrap_or("").to_string();
            seqs.insert(n.clone(), Vec::new());
            cur = Some(n);
        } else if let Some(n) = &cur {
            seqs.get_mut(n).unwrap().extend_from_slice(ln.as_bytes());
        }
    }
    let mut out: BTreeMap<String, Vec<IsoconCandidate>> = BTreeMap::new();
    for (c, r) in &cand {
        let Some(m) = members.get(c) else { continue };
        if r["flagged"] == "1" && !m.is_empty() {
            if let Some(seq) = seqs.get(c) {
                let n_b = m.iter().filter(|x| via.get(*x).map(String::as_str) == Some("B")).count();
                out.entry(r["family"].clone()).or_default().push((c.clone(), m.iter().map(|x| x.1.clone()).collect(), seq.clone(), n_b));
            }
        }
    }
    Ok(out)
}

/// `isocon_keep`: v5 (no-pile) >= half of the members in no phased pile; v4 (via-b) >= half joined the net by attribution.
pub fn isocon_keep(members: &BTreeSet<String>, n_via_b: usize, piled: &BTreeSet<String>, via_b_rule: bool) -> bool {
    if via_b_rule {
        return 2 * n_via_b >= members.len();
    }
    2 * members.difference(piled).count() >= members.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set(xs: &[&str]) -> BTreeSet<String> {
        xs.iter().map(|x| x.to_string()).collect()
    }

    #[test]
    fn isocon_keep_rules() {
        let members = set(&["a", "b", "c", "d"]);
        assert!(isocon_keep(&members, 0, &set(&["a", "b"]), false));
        assert!(!isocon_keep(&members, 0, &set(&["a", "b", "c"]), false));
        assert!(isocon_keep(&members, 2, &set(&["a", "b", "c", "d"]), true));
        assert!(!isocon_keep(&members, 1, &set(&[]), true));
    }

    #[test]
    fn empty_family_row() {
        // a family whose reads land on no copy: nothing phased, nothing certified, no tied read
        let dir = std::env::temp_dir().join(format!("copy_phasing_empty_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("g.fa"), ">c\nACGT\n").unwrap();
        std::fs::write(dir.join("g.fa.fai"), "c\t4\t3\t4\t5\n").unwrap();
        let g = IndexedFasta::open(dir.join("g.fa").to_str().unwrap()).unwrap();
        let copies = vec![Copy { tid: "F:0".into(), chrom: "c".into(), start: 0, end: 4 }];
        let piles: Piles = [("F:0".to_string(), Default::default())].into_iter().collect();
        let out = run_family("F", &copies, &piles, &HashMap::new(), &g, &Params::default(), &dir, Some(&[]), &HashMap::new()).unwrap();
        assert!(out.certified.is_empty() && out.looks.is_empty() && out.o2.is_empty() && out.per_copy.is_empty());
        assert_eq!(out.n_h, 0);
        assert!(out.why.is_empty());
        assert_eq!(
            Json::Obj(out.dump.unwrap()).dumps(),
            r#"{"family": "F", "copies": {"F:0": {"pile": []}}, "clusters": [], "merged": []}"#
        );
        let _ = std::fs::remove_dir_all(dir);
    }
}
