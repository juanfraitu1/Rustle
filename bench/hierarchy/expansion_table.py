#!/usr/bin/env python3
"""expansion_table — the T2 E-table (Phase E, expansions; DNA side) of docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md, section 7.

An *expansion* in gorilla is a set of at least two annotated records (gene or pseudogene, exon union >= 200 bp) of one POOL component that are
co-orthologs of ONE locus that is single-copy in each of the outgroups H = CHM13v2.0, C = mPanTro3, O = mPonPyg2. This module holds everything that is
not a Liftoff run (that is liftoff_cross.py):

  pool        the gorilla self-lift table (figures/_liftoff.py `liftoff_loci.tsv`): copy rows (extra_copy, dropped_M1, dropped_M2) join annotated records
              through the note column ('overlaps gene-X' / 'overlaps copy of X') or to new `copy@contig:start0-end<strand>` nodes; union-find
              (`build_pool`); hubs (more than 51 members) are not called; the called set is the non-hub components with >= 2 such records
              (`called_components`)
  rank plan   run j holds the j-th record (sorted by contig as a Python string, start0, end, id) of every component of the plan, so no run holds two
              records of one component (`rank_plan`, `plan_run_shards`, `split_shard`)
  E-table     placement rule, outgroup-locus rule, triples, co-orthology as an equivalence relation, expansions, siblings, control pool, the outputs of
              step 5 (`classify_record`, `classify_components`, `call_expansions`, `write_outputs`). A record with a reverse copy row keeps the locus of
              its primary placement and is flagged `multiplicity`; a class (>= 2 records) that contains such a record is dropped as a whole and counted
              (the registered rule of step 3; Amendment 2 readings 5 and 6 restated)
  Gate 7      the LRPAP1 positive control and the planted negative control (`lrpap_control`, `split_control`, `check_positive_control`, `plant_negative`,
              `check_negative`). The two controls are lifted and checked FIRST, in calls of their own, before any verdict-set record is lifted: the
              LRPAP1 component and the fragment component leave the verdict plan and have a rank plan of their own (liftoff_cross.py)

CLI (light):  python3 -B bench/hierarchy/expansion_table.py pool|plan [--hub-mode nodes|records] [--out DIR]
"""
from __future__ import annotations

import argparse
import bisect
import collections
import csv
import json
import math
import os
import random
import re
import sys

LIFT_ROOT = "/mnt/linuxdisk/tmp/rustle_figures/liftoff"
COPY_CLASSES = ("extra_copy", "dropped_M1", "dropped_M2")
PRIMARY_CLASSES = ("in_place", "moved", "partial", "unmapped")
PLACED_CLASSES = ("in_place", "moved")
SHORT_BP = 200                 # exon union below this is never a member of a called expansion (figures/_liftoff.py SHORT_BP)
HUB_MAX = 51                   # more than this many members = a hub, not called (inherited: where Liftoff itself truncates whole-body copies)
MATCH_FRAC = 0.5               # a placement maps to an outgroup locus iff it shares >= 50% of the SHORTER exon union
RANK_CAP = 2500                # records per Liftoff shard (figures/_liftoff.py LIFTOFF_BLOCK_RECORDS)
MIN_COVERAGE = 0.5             # Liftoff's -a
MIN_SEQ_ID = 0.5               # Liftoff's -s
SEED_T2 = 20261007             # the registration date, the T2 seed
OVL = re.compile(r"overlaps (?:copy of )?(\S+)")
OUTGROUPS = (("H", "human"), ("C", "chimpanzee"), ("O", "orangutan"))


# ----------------------------------------------------------------------------------------------- intervals
def parse_iv(s: str) -> list[tuple[int, int]]:
    """'5-10,20-30' -> [(5, 10), (20, 30)] (0-based half-open, as the self-lift tables write them)"""
    if s in ("", "-", "NA"):
        return []
    out = []
    for part in s.split(","):
        a, b = part.split("-")
        out.append((int(a), int(b)))
    return out


def merge_iv(iv) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def iv_len(iv) -> int:
    return sum(e - s for s, e in iv)


def iv_inter(a, b) -> int:
    """bases shared by two merged, sorted interval lists"""
    i = j = n = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            n += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return n


def iv_str(iv) -> str:
    return ",".join(f"{s}-{e}" for s, e in iv)


def record_key(r: dict) -> tuple:
    return (r["contig"], int(r["start0"]), int(r["end"]), r["id"])


def record_label(r: dict) -> str:
    """record_id of every output table: contig:start0-end (0-based half-open)"""
    return f"{r['contig']}:{int(r['start0'])}-{int(r['end'])}"


# ----------------------------------------------------------------------------------------------- readers
def read_selflift(path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_records(path) -> dict[str, dict]:
    """records.tsv of a self-lift work dir: one row per gene / pseudogene record with its exon union"""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["start0"], r["end"], r["exonic_len"], r["idx"] = int(r["start0"]), int(r["end"]), int(r["exonic_len"]), int(r["idx"])
            r["exons"] = parse_iv(r["exons"])
            out[r["id"]] = r
    return out


# ----------------------------------------------------------------------------------------------- the pool
class Component:
    """a union-find component of the pool: annotated records (nodes) and unannotated copy sites (`copy@...` nodes)"""

    def __init__(self, cid: str):
        self.id = cid
        self.record_ids: list[str] = []
        self.sites: list[str] = []
        self.copy_rows: list[dict] = []
        self.recs: dict[str, dict] = {}

    def size(self, mode: str = "nodes") -> int:
        """members for the hub rule: 'nodes' = annotated records of any length + copy sites (the text of section 7), 'records' = annotated records only"""
        if mode == "nodes":
            return len(self.record_ids) + len(self.sites)
        if mode == "records":
            return len(self.record_ids)
        raise ValueError(mode)

    def hub(self, mode: str = "nodes", hub_max: int = HUB_MAX) -> bool:
        return self.size(mode) > hub_max

    def such(self) -> list[dict]:
        """placed annotated records (in_place or moved) with an exon union >= 200 bp, sorted by (contig, start0, end, id)"""
        rs = [r for r in self.recs.values() if r.get("cls") in PLACED_CLASSES and r["exonic_len"] >= SHORT_BP]
        return sorted(rs, key=record_key)


class Pool:
    def __init__(self):
        self.components: dict[str, Component] = collections.OrderedDict()
        self.node_component: dict[str, str] = {}
        self.records: dict[str, dict] = {}
        self.n_copy_rows = 0


def _record_from_row(r: dict) -> dict:
    return dict(id=r["source_id"], name=r.get("name", r["source_id"]), rtype=r.get("rtype", "gene"), biotype=r.get("biotype", "-"),
                contig=r["src_contig"], start0=int(r["src_start0"]), end=int(r["src_end"]), strand=r["src_strand"],
                exons=parse_iv(r["exons"]), exonic_len=int(r["exonic_len"]), idx=0)


def build_pool(rows: list[dict], records: dict[str, dict]) -> Pool:
    """union-find over the copy rows of a self-lift table (the registered rule); components of >= 2 nodes are kept"""
    parent: dict[str, str] = {}

    def node(x):
        parent.setdefault(x, x)
        return x

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    cls_of: dict[str, str] = {}
    first_row: dict[str, dict] = {}
    for r in rows:
        node(r["source_id"])
        first_row.setdefault(r["source_id"], r)
        if r["cls"] in PRIMARY_CLASSES:
            cls_of[r["source_id"]] = r["cls"]
    partner_of: list[tuple[dict, str]] = []
    n_copy = 0
    for r in rows:
        if r["cls"] not in COPY_CLASSES:
            continue
        n_copy += 1
        m = OVL.search(r["note"])
        if m:
            p = node(m.group(1))
        else:
            p = node(f"copy@{r['contig']}:{r['start0']}-{r['end']}{r['strand']}")
        union(r["source_id"], p)
        partner_of.append((r, p))
    pool = Pool()
    pool.n_copy_rows = n_copy
    members: dict[str, set] = collections.defaultdict(set)
    for x in list(parent):
        members[find(x)].add(x)
    comps = []
    for root, mem in members.items():
        if len(mem) < 2:
            continue
        ann = sorted(x for x in mem if not x.startswith("copy@"))
        sites = sorted(x for x in mem if x.startswith("copy@"))
        recs = {}
        for x in ann:
            rec = dict(records[x]) if x in records else _record_from_row(first_row[x])
            rec["cls"] = cls_of.get(x, "unmapped")
            recs[x] = rec
        cid = "PC_" + min(ann)
        c = Component(cid)
        c.record_ids, c.sites, c.recs = ann, sites, recs
        comps.append((root, c))
    comps.sort(key=lambda t: t[1].id)
    by_root = {}
    for root, c in comps:
        pool.components[c.id] = c
        by_root[root] = c
        for x in list(c.record_ids) + list(c.sites):
            pool.node_component[x] = c.id
    for r, p in partner_of:
        c = by_root.get(find(r["source_id"]))
        if c is not None:
            c.copy_rows.append(dict(source_id=r["source_id"], partner=p, cls=r["cls"], sequence_id=r["sequence_id"], note=r["note"],
                                    contig=r["contig"], start0=int(r["start0"]), end=int(r["end"]), strand=r["strand"], coverage=r["coverage"]))
    for c in pool.components.values():
        pool.records.update(c.recs)
    return pool


def called_components(pool: Pool, hub_mode: str = "nodes", hub_max: int = HUB_MAX) -> list[Component]:
    """non-hub components with at least two such records (placed, annotated, exon union >= 200 bp)"""
    return [c for c in pool.components.values() if not c.hub(hub_mode, hub_max) and len(c.such()) >= 2]


# ----------------------------------------------------------------------------------------------- the rank-sharded plan
def record_free_cuts(rs: list[dict]) -> list[int]:
    """positions k (1..n-1) of a list sorted by (contig, start0, end, id) before which no record of the same contig reaches: a record-free cut.
    A change of contig is always one (the generalisation of figures/_liftoff.py::_apply_splits to a run that holds records of many contigs)."""
    cuts, max_end, contig = [], -1, None
    for k, r in enumerate(rs):
        if r["contig"] != contig:
            if k:
                cuts.append(k)
            contig, max_end = r["contig"], -1
        elif r["start0"] >= max_end:
            cuts.append(k)
        max_end = max(max_end, r["end"])
    return cuts


class NoSplit(RuntimeError):
    pass


def _sorted_records(recs: list[dict]) -> list[dict]:
    return sorted(recs, key=lambda r: (r["contig"], r["start0"], r["end"], r["id"]))


def plan_run_shards(label: str, recs: list[dict], cap: int = RANK_CAP) -> list[tuple[str, list[str]]]:
    """one shard per run (label), or blocks `label.b1`... of <= cap records cut where no record spans (a block may exceed cap only when no
    record-free cut exists), balanced like figures/_liftoff.py::_plan_blocks"""
    rs = _sorted_records(recs)
    if len(rs) <= cap:
        return [(label, [r["id"] for r in rs])]
    target = math.ceil(len(rs) / math.ceil(len(rs) / cap))
    cuts = set(record_free_cuts(rs))
    blocks, cur = [], []
    for k, r in enumerate(rs):
        if cur and len(cur) >= target and k in cuts:
            blocks.append(cur)
            cur = []
        cur.append(r)
    if cur:
        blocks.append(cur)
    return [(f"{label}.b{i}", [r["id"] for r in b]) for i, b in enumerate(blocks, 1)]


def split_shard(label: str, ids: list[str], by_id: dict[str, dict]) -> list[tuple[str, list[str]]]:
    """a shard that did not finish inside the hard limit is split in two at the record-free cut nearest its middle (labels L.s1 / L.s2)"""
    rs = _sorted_records([by_id[i] for i in ids])
    cuts = record_free_cuts(rs)
    if not cuts:
        raise NoSplit(f"shard {label} ({len(rs)} records) has no record-free cut: stop rule, run it on the cluster")
    k = min(cuts, key=lambda c: (abs(c - len(rs) / 2), c))
    return [(f"{label}.s1", [r["id"] for r in rs[:k]]), (f"{label}.s2", [r["id"] for r in rs[k:]])]


def rank_plan(called: list[Component], cap: int = RANK_CAP, label_prefix: str = "r", pad: int | None = None) -> list[dict]:
    """run j (j = 1..J) holds the j-th such record (sorted by contig, start0) of every called component with at least j such records.
    Labels are r01.. for the called set and ctl1, ctl2 for the planted negative control (label_prefix='ctl')."""
    comps = sorted(called, key=lambda c: c.id)
    such = {c.id: c.such() for c in comps}
    J = max((len(v) for v in such.values()), default=0)
    width = pad if pad is not None else (max(2, len(str(J))) if label_prefix == "r" else 1)
    plan = []
    for j in range(1, J + 1):
        records = []
        for c in comps:
            if len(such[c.id]) >= j:
                r = such[c.id][j - 1]
                records.append(dict(id=r["id"], component=c.id, contig=r["contig"], start0=r["start0"], end=r["end"], strand=r["strand"],
                                    idx=r.get("idx", 0), rtype=r.get("rtype", "gene"), biotype=r.get("biotype", "-"), exonic_len=r["exonic_len"],
                                    rank=j))
        label = f"{label_prefix}{j:0{width}d}"
        plan.append(dict(run=j, label=label, records=records, shards=plan_run_shards(label, records, cap)))
    return plan


def plan_summary(plan: list[dict]) -> dict:
    per = [len(r["records"]) for r in plan]
    return dict(J=len(plan), records=sum(per), records_per_run=per, runs_with_at_most_10=sum(1 for n in per if n <= 10),
                shards=sum(len(r["shards"]) for r in plan))


# ----------------------------------------------------------------------------------------------- the E-table: placement and outgroup locus
def placement_ok(r: dict) -> bool:
    """a PLACEMENT is the primary row (extra_copy_number 0) of a record with coverage >= 0.5, sequence_ID >= 0.5 and neither the partial_mapping nor the
    low_identity flag set (whole-genome Liftoff writes rows from coverage 0.05 and marks the weaker ones with those flags); any other row is no placement"""
    try:
        cov, sid = float(r["coverage"]), float(r["sequence_id"])
    except (ValueError, KeyError):
        return False
    return int(r["tag"]) == 0 and cov >= MIN_COVERAGE and sid >= MIN_SEQ_ID and str(r["partial"]) == "0" and str(r["low_identity"]) == "0"


class OutgroupIndex:
    """the in_place records of one outgroup's own self-lift (liftoff/<species>/liftoff_loci.tsv) by contig, with their exon unions.

    `match(contig, exon_iv)`: the candidates are the in_place records on the contig that share >= 50% of the SHORTER of the two exon unions with the query
    (strand ignored). One candidate that has no copy row is the locus (status ok); a candidate that has a copy row is an outgroup_copy_row; two or more
    single-copy candidates are ambiguous; none is no_outgroup_locus. "Has a copy row" = a row of class extra_copy / dropped_M1 / dropped_M2 whose source
    is that record (the reading that gives the 19,841 / 21,970 / 21,704 single-copy protein-coding anchors of tier A); strict_single_copy also counts a
    record that is the landing site named by the note of another record's copy row.
    """

    def __init__(self, rows: list[dict], name: str = "", strict_single_copy: bool = False):
        self.name = name
        copy_src = {r["source_id"] for r in rows if r["cls"] in COPY_CLASSES}
        if strict_single_copy:
            for r in rows:
                if r["cls"] in COPY_CLASSES:
                    m = OVL.search(r.get("note", ""))
                    if m:
                        copy_src.add(m.group(1))
        self.multi = copy_src
        by = collections.defaultdict(list)
        for r in rows:
            if r["cls"] == "in_place":
                ex = merge_iv(parse_iv(r["exons"]))
                by[r["contig"]].append((int(r["start0"]), int(r["end"]), ex, r["source_id"]))
        self.by_contig, self.starts, self.maxlen = {}, {}, {}
        for contig, lst in by.items():
            lst.sort(key=lambda t: (t[0], t[1], t[3]))
            self.by_contig[contig] = lst
            self.starts[contig] = [t[0] for t in lst]
            self.maxlen[contig] = max(t[1] - t[0] for t in lst)

    def match(self, contig: str, exon_iv, frac: float = MATCH_FRAC) -> dict:
        q = merge_iv(exon_iv)
        lq = iv_len(q)
        cand = []
        if lq and contig in self.by_contig and q:
            lst, starts = self.by_contig[contig], self.starts[contig]
            lo = bisect.bisect_left(starts, q[0][0] - self.maxlen[contig])
            for i in range(lo, len(lst)):
                s0, e0, ex, rid = lst[i]
                if s0 >= q[-1][1]:
                    break
                if e0 <= q[0][0]:
                    continue
                inter = iv_inter(q, ex)
                short = min(lq, iv_len(ex))
                if short and inter >= frac * short - 1e-9:
                    cand.append(rid)
        cand = sorted(set(cand))
        if not cand:
            return dict(status="no_outgroup_locus", locus=None, candidates=[])
        multi = [c for c in cand if c in self.multi]
        if multi:
            return dict(status="outgroup_copy_row", locus=multi[0], candidates=cand)
        if len(cand) == 1:
            return dict(status="ok", locus=cand[0], candidates=cand)
        return dict(status="ambiguous", locus=None, candidates=cand)


REASON_ORDER = ("outgroup_copy_row", "no_placement", "no_outgroup_locus")


def has_multiplicity(info: dict) -> bool:
    """does the record have more than one placement in this outgroup? True when classify_record flagged it (a reverse copy row of any coverage, or more than one
    qualifying placement); for calls assembled by hand also a copy_rows count or an n_place above one. A record that has no placement has no multiplicity."""
    if info.get("status") == "no_placement":
        return False
    return bool(info.get("multiplicity")) or info.get("copy_rows", 0) > 0 or info.get("n_place", 0) > 1


def classify_record(rid: str, lifted: dict, unmapped: dict, indexes: dict) -> dict:
    """the call of one gorilla record: per outgroup the placement, the outgroup locus and the counts, and the triple (h, c, o) when all three are defined.
    lifted[o][rid] = the rows (primary and copy rows) of the record in outgroup o; unmapped[o] = ids listed unmapped.

    The registered rule (prereg section 7 step 3, Amendment 2 reading 5 restated): the number of placements of a record is recorded and a class with a member that
    has more than one is dropped and counted (call_expansions). So a record with a reverse copy row (any row with extra_copy_number > 0, whatever its coverage) KEEPS
    the locus of its primary placement as a candidate and is flagged `multiplicity` in that outgroup; it is not turned into a record without a triple. A record whose
    primary row is not a placement (partial, weak, flagged, absent) stays without one. The `outgroup_copy_row` status is a different thing: the OUTGROUP's own
    record at the locus has a copy row in the outgroup's self-lift (reading 6)."""
    per = {}
    for o, _ in OUTGROUPS:
        rows = lifted[o].get(rid, [])
        n_place = sum(1 for r in rows if _cov(r) >= MIN_COVERAGE)
        copies = sum(1 for r in rows if int(r["tag"]) > 0)
        primary = next((r for r in rows if int(r["tag"]) == 0 and r["feature_id"] == rid), None)
        info = dict(status=None, detail="", locus=None, n_place=n_place, copy_rows=copies, candidates=[], multiplicity=False)
        if primary is None:
            info["status"] = "no_placement"
            info["detail"] = "unmapped" if rid in unmapped.get(o, ()) else ("no_row" if not rows else "no_primary_row")
        elif not placement_ok(primary):
            # "partial_only": the record's only row is flagged partial_mapping (the count the pre-registration prints); anything else that fails the rule
            # (low_identity, coverage or identity below one half, or other rows besides) is "not_a_placement"
            info["status"] = "no_placement"
            info["detail"] = "partial_only" if (str(primary["partial"]) != "0" and len(rows) == 1) else "not_a_placement"
        else:
            m = indexes[o].match(primary["contig"], parse_iv(primary["exons"]))
            info["status"], info["locus"], info["candidates"] = m["status"], m["locus"], m["candidates"]
            info["detail"] = {"ok": "", "outgroup_copy_row": "outgroup_locus_has_copy_row", "ambiguous": "ambiguous", "no_outgroup_locus": ""}[m["status"]]
            if copies or n_place > 1:
                info["multiplicity"] = True
                if m["status"] == "ok":
                    info["detail"] = "reverse_copy_row"
        per[o] = info
    statuses = [per[o]["status"] for o, _ in OUTGROUPS]
    flagged = [o for o, _ in OUTGROUPS if per[o]["multiplicity"]]
    if all(s == "ok" for s in statuses):
        return dict(triple=tuple(per[o]["locus"] for o, _ in OUTGROUPS), reason=None, per=per, multiplicity=flagged)
    reason = next((r for r in REASON_ORDER if r in statuses), "no_outgroup_locus")      # an ambiguous locus is reported as no_outgroup_locus (detail: ambiguous)
    return dict(triple=None, reason=reason, per=per, multiplicity=flagged)


def _cov(r: dict) -> float:
    try:
        return float(r["coverage"])
    except (ValueError, KeyError):
        return -1.0


def classify_components(components: list[Component], lifted: dict, unmapped: dict, indexes: dict) -> dict:
    """the call of every such record of the components: {record id: classify_record(...)}"""
    return {r["id"]: classify_record(r["id"], lifted, unmapped, indexes) for c in components for r in c.such()}


# ----------------------------------------------------------------------------------------------- expansions, siblings, the control pool
ALL_OUTGROUPS = tuple(o for o, _ in OUTGROUPS)


def class_key(call: dict, use: tuple) -> tuple | None:
    """the locus tuple that decides co-orthology: the triple (all three outgroups), or the loci of the outgroups in `use` when each is ok (the human-only
    sensitivity of Amendment 3 uses ('H',))"""
    if tuple(use) == ALL_OUTGROUPS:
        return call["triple"]
    loci = tuple(call["per"][o].get("locus") if call["per"][o].get("status") == "ok" else None for o in use)
    return loci if all(loci) else None


def call_expansions(called: list[Component], calls: dict, use: tuple = ALL_OUTGROUPS, prefix: str = "EX") -> dict:
    """co-orthology is an equivalence relation: two such records of a component are co-orthologous iff they have the same defined triple; an expansion is
    the set of >= 2 such records of one component with the same triple (no chaining); a class (>= 2 members) with a member that has more than one placement in an
    outgroup of `use` (has_multiplicity: a reverse copy row, or more than one placement of coverage >= 0.5) is DROPPED AS A WHOLE and counted, never just that member;
    expansions of different components with an identical triple are siblings and are not merged.
    `use` = ('H',) reduces the triple to the human locus (the same lifts; ids start with `prefix`)."""
    expansions, dropped, comp_info = [], [], {}
    reasons: collections.Counter = collections.Counter()
    by_triple_exp: dict = collections.defaultdict(list)
    n = 0
    for c in sorted(called, key=lambda c: c.id):
        such = c.such()
        groups: dict = collections.defaultdict(list)
        for r in such:
            call = calls[r["id"]]
            key = class_key(call, use)
            if key is None:
                reasons[call["reason"] or "no_placement"] += 1
            else:
                groups[key].append(r)
        contributes = False
        for triple, members in sorted(groups.items(), key=lambda kv: record_key(kv[1][0])):
            if len(members) < 2:
                continue
            if any(has_multiplicity(calls[m["id"]]["per"][o]) for m in members for o in use):
                dropped.append(dict(component=c.id, triple=triple, members=[m["id"] for m in members], reason="multiplicity"))
                continue
            n += 1
            ex = dict(id=f"{prefix}{n:04d}", component=c.id, triple=triple, members=members)
            expansions.append(ex)
            by_triple_exp[triple].append(ex["id"])
            contributes = True
        comp_info[c.id] = dict(n_such=len(such), different_loci=len(groups), contributes=contributes)
    siblings = [dict(triple=t, expansions=ids) for t, ids in sorted(by_triple_exp.items()) if len(ids) >= 2]
    return dict(expansions=expansions, dropped_classes=dropped, siblings=siblings, components=comp_info, reasons=dict(reasons),
                records_in_expansions=sum(len(e["members"]) for e in expansions))


CONTROL_REASONS = ("outgroup_copy_row", "different_triple", "no_placement")      # the reasons control_pool.tsv lists (section 7, step 5)
NOT_ELIGIBLE = ("no_placement",)                                                  # "a record with reason no_placement is not eligible as a member of a control set"


def control_pool_rows(called: list[Component], calls: dict, res: dict, joining: dict) -> list[dict]:
    """the control pool: non-hub components that contribute no expansion and either hold a record that places on an outgroup locus with a copy row (in
    the outgroup's own self-lift: the OUTGROUP's record has a copy row), or a record flagged multiplicity (the reverse lift gave it a second place: its class is
    dropped in call_expansions), or hold records with different triples. One row per record with one of the three registered reasons: outgroup_copy_row,
    different_triple (a record with a triple in a component with two or more triples), no_placement. A FLAGGED record is listed as outgroup_copy_row whatever its
    call says (before D2 revised a reverse copy row gave the record that reason, ahead of no_placement, REASON_ORDER; listing it keeps the pool what it was).
    `eligible` is 0 for reason no_placement (nothing is known about its locus) and 1 for the others. A record placed on no annotated outgroup locus
    (no_outgroup_locus, ambiguous loci included) and a record whose triple is the only one of its component have none of the three reasons and are not listed
    (they stay in record_calls.tsv)."""
    flagged = lambda r: any(has_multiplicity(calls[r["id"]]["per"][o]) for o in ALL_OUTGROUPS)  # noqa: E731
    rows = []
    for c in sorted(called, key=lambda c: c.id):
        info = res["components"][c.id]
        if info["contributes"]:
            continue
        such = c.such()
        if not (any(calls[r["id"]]["reason"] == "outgroup_copy_row" or flagged(r) for r in such) or info["different_loci"] >= 2):
            continue
        for r in such:
            call = calls[r["id"]]
            if flagged(r):
                reason = "outgroup_copy_row"
            elif call["triple"] is None:
                reason = call["reason"] or "no_placement"
            else:
                reason = "different_triple" if info["different_loci"] >= 2 else None
            if reason not in CONTROL_REASONS:
                continue
            jr = ";".join(f"{p}:{cls}:{sid}" for p, cls, sid in joining.get(r["id"], []))
            rows.append(dict(pool_component_id=c.id, record_id=record_label(r), reason=reason, joining_copy_rows=jr,
                             eligible=0 if reason in NOT_ELIGIBLE else 1))
    return rows


def joining_copy_rows(pool: Pool) -> dict:
    """for every annotated record, the copy rows that join it to another annotated record: (partner, class, sequence_id)"""
    out: dict = collections.defaultdict(list)
    for c in pool.components.values():
        for cr in c.copy_rows:
            if cr["partner"].startswith("copy@"):
                continue
            out[cr["source_id"]].append((cr["partner"], cr["cls"], cr["sequence_id"]))
            out[cr["partner"]].append((cr["source_id"], cr["cls"], cr["sequence_id"]))
    return {k: sorted(v) for k, v in out.items()}


# ----------------------------------------------------------------------------------------------- reading the Liftoff runs
_COPY_ID = re.compile(r"(.+)_(\d+)")


def source_of(fid: str, known) -> str:
    """the gorilla record a Liftoff row belongs to: the feature id itself (the primary row, tag 0, keeps the record id) or the id without the `_<copy tag>`
    that Liftoff appends to the further copies"""
    return fid if fid in known else fid.rsplit("_", 1)[0]


def check_no_id_collisions(ids):
    """Liftoff names the copies of a record `<id>_<n>`, so a record whose id is `<another record's id>_<digits>` would be taken for a copy of it when the rows
    are read back (source_of). Raises ValueError for such ids; none exists in the 41,193 gorilla record ids."""
    known = set(ids)
    clashes = sorted(i for i in known if (m := _COPY_ID.fullmatch(i)) and m.group(1) in known)
    if clashes:
        raise ValueError(f"{len(clashes)} record id(s) look like the copy of another record of the plan ({', '.join(clashes[:5])}): rename them or read the rows another way")


def read_lifted_runs(outgroup_dir, known_ids, labels=None):
    """the finished runs of one outgroup (`<outgroup_dir>/runs/<label>/loci.tsv`, done.key, unmapped.txt) -> ({gorilla record: its rows}, {unmapped ids}).
    With `labels`, exactly those runs are read and each must be finished; without, every finished run is."""
    base = os.path.join(outgroup_dir, "runs")
    present = sorted(os.listdir(base)) if os.path.isdir(base) else []
    done = [lab for lab in present if os.path.exists(os.path.join(base, lab, "done.key"))]
    if labels is not None:
        missing = [lab for lab in labels if lab not in done]
        if missing:
            raise RuntimeError(f"{len(missing)} run(s) of {outgroup_dir} are not finished: {missing[:6]}")
        done = sorted(set(labels))
    lifted: dict = collections.defaultdict(list)
    unmapped: set = set()
    known = set(known_ids)
    for lab in done:
        d = os.path.join(base, lab)
        loci = os.path.join(d, "loci.tsv")
        try:
            with open(loci) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    src = source_of(r["feature_id"], known)
                    if src in known:
                        lifted[src].append(r)
        except (OSError, UnicodeDecodeError, csv.Error, KeyError) as e:
            raise RuntimeError(f"{loci} cannot be read ({type(e).__name__}: {e}): this finished run is damaged") from e
        up = os.path.join(d, "unmapped.txt")
        if not os.path.exists(up):
            raise RuntimeError(f"{d} has no unmapped.txt: Liftoff writes one for every run (empty when nothing is unmapped), so this is not a finished run")
        try:
            with open(up) as fh:
                unmapped.update(x.strip() for x in fh if x.strip())
        except (OSError, UnicodeDecodeError) as e:
            raise RuntimeError(f"{up} cannot be read ({type(e).__name__}: {e}): this finished run is damaged") from e
    return dict(lifted), unmapped


ANOMALY_DETAILS = ("no_row", "no_primary_row")


def anomaly_counts(calls: dict) -> dict:
    """{outgroup: {detail: n}} for the records whose rows are inconsistent: no row at all and not listed unmapped (no_row), or copy rows without a primary row
    (no_primary_row). A run is accepted only if every record is lifted or listed unmapped, and in the four real self-lifts no source with copy rows lacks its
    primary row, so any such record signals a parse or run anomaly, not biology."""
    out = {o: collections.Counter() for o, _ in OUTGROUPS}
    for c in calls.values():
        for o, _ in OUTGROUPS:
            d = c["per"][o].get("detail")
            if d in ANOMALY_DETAILS:
                out[o][d] += 1
    return {o: dict(v) for o, v in out.items() if v}


# ----------------------------------------------------------------------------------------------- outputs
def _w(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def _sha(path) -> str:
    import hashlib
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


ETABLE_HEADER = ["expansion_id", "pool_component_id", "record_id", "biotype", "exon_union_bp", "h", "c", "o", "kind"]


def etable_rows(expansions: list[dict]) -> list[list]:
    """E-table rows; a class defined by fewer than three outgroups leaves the other columns blank"""
    return [[e["id"], e["component"], record_label(m), m.get("biotype", "-"), m["exonic_len"], *(list(e["triple"]) + [""] * (3 - len(e["triple"]))), "annotated"]
            for e in expansions for m in sorted(e["members"], key=record_key)]


def component_status(c: Component, hub_mode: str, verdict_ids: set, control_ids: set) -> str:
    if c.hub(hub_mode):
        return "hub"
    if c.id in control_ids:
        return "control"
    if c.id in verdict_ids:
        return "verdict"
    return "too_few_such_records" if len(c.such()) < 2 else "not_called"


def write_outputs(res: dict, pool: Pool, called: list[Component], calls: dict, out_dir, rows=None, extra: dict | None = None, hub_mode: str = "nodes",
                  control=()) -> dict:
    """the files of step 5: etable.tsv (the deciding E-table), member_exons.tsv, control_pool.tsv, extension.tsv, and siblings.tsv, record_calls.tsv,
    dropped_classes.tsv, components.tsv (every pool component with its status: hub, control, verdict, too_few_such_records), etable_human_only.tsv (the
    sensitivity of Amendment 3: the triple reduced to the human locus, the same lifts), summary.json (with the sha256 of every file and `extra`).
    `called` is the verdict set; member_exons.tsv, extension.tsv and components.tsv list the whole pool."""
    os.makedirs(out_dir, exist_ok=True)
    P = lambda n: os.path.join(out_dir, n)  # noqa: E731
    _w(P("etable.tsv"), ETABLE_HEADER, etable_rows(res["expansions"]))
    res_h = call_expansions(called, calls, use=("H",), prefix="EH")
    _w(P("etable_human_only.tsv"), ETABLE_HEADER, etable_rows(res_h["expansions"]))
    verdict_ids, control_ids = {c.id for c in called}, {c.id for c in control}
    comp_rows = []
    for c in pool.components.values():
        n_exp = sum(1 for e in res["expansions"] if e["component"] == c.id)
        comp_rows.append([c.id, c.size("nodes"), len(c.record_ids), len(c.sites), len(c.such()), component_status(c, hub_mode, verdict_ids, control_ids), n_exp])
    _w(P("components.tsv"), ["pool_component_id", "members_nodes", "annotated_records", "copy_sites", "such_records", "status", "expansions"], comp_rows)
    me = []
    for c in pool.components.values():
        for rid in sorted(c.record_ids, key=lambda x: record_key(c.recs[x])):
            r = c.recs[rid]
            me.append([record_label(r), c.id, rid, iv_str(merge_iv(r["exons"]))])
    _w(P("member_exons.tsv"), ["record_id", "pool_component_id", "gene_id", "exons"], me)
    joining = joining_copy_rows(pool)
    cp = control_pool_rows(called, calls, res, joining)
    _w(P("control_pool.tsv"), ["pool_component_id", "record_id", "reason", "joining_copy_rows", "eligible"],
       [[r["pool_component_id"], r["record_id"], r["reason"], r["joining_copy_rows"], r["eligible"]] for r in cp])
    ext = []
    for c in pool.components.values():
        for cr in c.copy_rows:
            if cr["partner"].startswith("copy@"):
                ext.append([cr["partner"], c.id, cr["contig"], cr["start0"], cr["end"], cr["strand"], cr["source_id"], cr["sequence_id"], cr["cls"]])
    _w(P("extension.tsv"), ["site_id", "pool_component_id", "contig", "start0", "end", "strand", "source_id", "sequence_id", "cls"], sorted(ext))
    _w(P("siblings.tsv"), ["h", "c", "o", "expansions"], [[*s["triple"], ",".join(s["expansions"])] for s in res["siblings"]])
    rc = []
    n_multi = {o: 0 for o, _ in OUTGROUPS}
    for c in sorted(called, key=lambda c: c.id):
        for r in c.such():
            k = calls[r["id"]]
            per = k["per"]
            flagged = [o for o, _ in OUTGROUPS if has_multiplicity(per[o])]
            for o in flagged:
                n_multi[o] += 1
            rc.append([record_label(r), r["id"], c.id, k["reason"] or "", ",".join(k["triple"]) if k["triple"] else ""] +
                      [x for o, _ in OUTGROUPS for x in (per[o].get("status") or "", per[o].get("detail") or "", per[o]["n_place"], per[o].get("copy_rows", 0),
                                                         (per[o].get("locus") or "") if per[o].get("status") == "ok" else "")] + [",".join(flagged)])
    hdr = ["record_id", "gene_id", "pool_component_id", "reason", "triple"] + [f"{o}_{f}" for o, _ in OUTGROUPS
                                                                                for f in ("status", "detail", "n_place", "copy_rows", "locus")] + ["multiplicity"]
    _w(P("record_calls.tsv"), hdr, rc)
    _w(P("dropped_classes.tsv"), ["pool_component_id", "h", "c", "o", "members", "reason"],
       [[d["component"], *d["triple"], ",".join(d["members"]), d["reason"]] for d in res["dropped_classes"]])
    detail = {o: collections.Counter() for o, _ in OUTGROUPS}
    for c in called:
        for r in c.such():
            for o, _ in OUTGROUPS:
                if calls[r["id"]]["per"][o].get("status") == "no_placement":
                    detail[o][calls[r["id"]]["per"][o].get("detail") or "-"] += 1
    dl = [info["different_loci"] for info in res["components"].values()]
    summary = dict(called_components=len(called), expansions=len(res["expansions"]), records_in_expansions=res["records_in_expansions"],
                   sibling_groups=len(res["siblings"]), dropped_classes=len(res["dropped_classes"]),
                   records_in_dropped_classes=sum(len(d["members"]) for d in res["dropped_classes"]), records_with_reverse_copy_rows=n_multi,
                   records_without_triple_by_reason=res["reasons"],
                   no_placement_by_detail={o: dict(sorted(c.items())) for o, c in detail.items()},
                   human_only_expansions=len(res_h["expansions"]), human_only_records_in_expansions=res_h["records_in_expansions"],
                   human_only_dropped_classes=len(res_h["dropped_classes"]),
                   hubs=sum(1 for r in comp_rows if r[5] == "hub"),
                   hubs_with_two_such_records=sum(1 for c in pool.components.values() if c.hub(hub_mode) and len(c.such()) >= 2),
                   components_with_a_triple=sum(1 for n in dl if n >= 1), components_with_different_loci=sum(1 for n in dl if n >= 2),
                   classes_in_components_with_different_loci=sum(n for n in dl if n >= 2),
                   control_pool_rows=len(cp), control_pool_rows_eligible=sum(r["eligible"] for r in cp))
    summary.update(extra or {})
    names = ["etable.tsv", "member_exons.tsv", "control_pool.tsv", "extension.tsv", "siblings.tsv", "record_calls.tsv", "dropped_classes.tsv", "components.tsv",
             "etable_human_only.tsv"]
    summary["sha256"] = {n: _sha(P(n)) for n in names}
    with open(P("summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True)
    return {n: P(n) for n in names + ["summary.json"]}


# ----------------------------------------------------------------------------------------------- Gate 7: the controls of the E-table
LRPAP1_CONTROL = dict(reference="gene-LRPAP1", not_required=["gene-LOC129523503", "gene-LOC129526389"],
                      fragments=["gene-LOC134756368", "gene-LOC115932954", "gene-LOC115932756"], min_share=0.5,
                      expected_triple=("gene-LRPAP1", "gene-LRPAP1", "gene-LRPAP1"))      # the gene-LRPAP1 anchor of H, C, O (in_place, single-copy in the three self-lifts)


def lrpap_control(pool: Pool, rows: list[dict], spec: dict | None = None) -> dict:
    """The positive control: the gene-LRPAP1 pool component. A record qualifies iff at least 50% of its exon bases lie inside an LRPAP1 copy interval: the
    landing interval of one of the rows of the reference record in the gorilla self-lift (the in_place row and the seven dropped_M1 landing sites). The
    two long gene models that contain a copy as a minor part are NOT REQUIRED (either may or may not reach coverage 0.5 in the reverse lift; placed on the
    outgroup LRPAP1 locus it joins the class without a failure). The three 5' fragments form another component and can appear only as a sibling class.
    The spec's `expected_triple` (required) is carried into the control: check_positive_control passes only a class that sits on it."""
    spec = spec or LRPAP1_CONTROL
    expected_triple = tuple(spec["expected_triple"])               # a specification without one is an error: the control would silently judge less
    ref = spec["reference"]
    if ref not in pool.node_component:
        raise ValueError(f"{ref} is not in a pool component")
    comp = pool.components[pool.node_component[ref]]
    ivs: dict = collections.defaultdict(list)
    for r in rows:
        if r["source_id"] == ref:
            ivs[r["contig"]].append((int(r["start0"]), int(r["end"])))
    ivs = {k: merge_iv(v) for k, v in ivs.items()}
    such = comp.such()
    shares = {}
    for r in such:
        ex = merge_iv(r["exons"])
        tot = iv_len(ex)
        shares[r["id"]] = iv_inter(ex, ivs.get(r["contig"], [])) / tot if tot else 0.0
    qualifying = [r["id"] for r in such if shares[r["id"]] >= spec["min_share"]]
    non = [r["id"] for r in such if r["id"] not in qualifying]
    named = [x for x in spec["not_required"] if x in non]
    frag_comps = {pool.node_component.get(f) for f in spec["fragments"]}
    if len(frag_comps) != 1 or None in frag_comps:
        raise ValueError(f"the fragments {spec['fragments']} are not in one pool component: {sorted(map(str, frag_comps))}")
    return dict(component=comp.id, reference=ref, shares=shares, qualifying=qualifying, not_required=named,
                other_non_qualifying=[x for x in non if x not in spec["not_required"]],
                spec_problems=[x for x in spec["not_required"] if x not in non],
                fragments=dict(component=next(iter(frag_comps)), records=list(spec["fragments"])), expected_triple=expected_triple)


def control_view(ctrl: dict, called: list[Component], calls: dict) -> tuple[dict, list[dict]]:
    """(the calls the positive control is judged on, the flagged records of the control components). The control is judged on the multiplicity flags of the
    QUALIFYING records only: a record of the control components that is not one of them (a not-required long gene model, a lncRNA flank, a fragment) and carries a
    reverse copy row is LEFT OUT of the control's classes (its triple is set to None, reason "multiplicity_left_out") so that it can neither fail the control nor take
    the class down with it. It stays flagged in `calls`, which are not changed, and the registered rule of call_expansions (a class with a flagged member is dropped
    as a whole) is not touched. `flagged` lists, sorted by record id, every flagged record of the control components: record, outgroups (H, C, O order), qualifying,
    left_out. When nothing is left out the view IS `calls`."""
    qualifying = set(ctrl["qualifying"])
    view, flagged = calls, []
    for c in sorted(called, key=lambda c: c.id):
        for r in c.such():
            call = calls[r["id"]]
            outs = [o for o in ALL_OUTGROUPS if has_multiplicity(call["per"][o])]
            if not outs:
                continue
            left_out = r["id"] not in qualifying
            flagged.append(dict(record=r["id"], outgroups=outs, qualifying=not left_out, left_out=left_out))
            if left_out:
                if view is calls:
                    view = dict(calls)
                view[r["id"]] = dict(call, triple=None, reason="multiplicity_left_out")
    return view, sorted(flagged, key=lambda d: d["record"])


def check_positive_control(ctrl: dict, res: dict, calls: dict) -> dict:
    """PASS iff ONE expansion contains every qualifying record (the bar is those records, 6 for LRPAP1) AND the triple of that expansion is the EXPECTED triple
    of the control (the gene-LRPAP1 anchors of the three outgroups; stricter than the registered text, which asks for the one class only). Missing records are
    listed with the reason; `triple_ok` is True or False when a class holds every qualifying record and None otherwise (nothing to compare); the not-required
    records are reported (in the class or not, and why not) and never fail the control; fragments_sibling says whether an expansion of the fragment component
    has the same triple (informative only)."""
    Q = list(ctrl["qualifying"])
    expected = tuple(ctrl["expected_triple"])
    members = {e["id"]: {m["id"] for m in e["members"]} for e in res["expansions"]}
    dropped = {m for d in res.get("dropped_classes", []) for m in d["members"]}
    whole = [e for e in res["expansions"] if set(Q) <= members[e["id"]]]
    target = None
    if whole:
        target = whole[0]
    else:
        scored = [(len(set(Q) & members[e["id"]]), -i, e) for i, e in enumerate(res["expansions"])]
        scored = [t for t in scored if t[0] > 0]
        target = max(scored, key=lambda t: (t[0], t[1]))[2] if scored else None
    missing = []
    if not whole:
        tm = members[target["id"]] if target else set()
        for q in Q:
            if q in tm:
                continue
            c = calls[q]
            if c["triple"] is None:
                reason = c["reason"] or "no_placement"
            elif any(q in members[e["id"]] for e in res["expansions"]):
                reason = "in another class"
            elif q in dropped:
                reason = "class dropped for multiplicity"
            else:
                reason = "no partner with the same triple"
            missing.append(dict(record=q, reason=reason))
    tmem = members[target["id"]] if target else set()
    exp = []
    for x in ctrl["not_required"]:
        c = calls.get(x)
        exp.append(dict(record=x, in_class=x in tmem, reason="" if x in tmem else ((c["reason"] or "no_partner") if c else "not called")))
    sib = None
    if target is not None:
        sib = any(e["triple"] == target["triple"] and e["component"] == ctrl["fragments"]["component"] for e in res["expansions"])
    triple_ok = (tuple(whole[0]["triple"]) == expected) if whole else None
    return {"pass": bool(whole) and triple_ok is True, "expansion": target["id"] if target else None, "class_members": sorted(tmem), "missing": missing,
            "not_required": exp, "fragments_sibling": sib, "expected_triple": list(expected), "triple_ok": triple_ok}


def split_control(called: list[Component], ctrl: dict) -> tuple[list[Component], list[Component]]:
    """(verdict set, control components): the LRPAP1 component and the fragment component are lifted in the control calls, before the verdict set, and are not
    part of the verdict set ('LRPAP1 is spent, is excluded from the verdict set and is reported separately'). The positive-control component must be called;
    the fragment component is taken along when it is."""
    ids = [ctrl["component"], ctrl["fragments"]["component"]]
    by_id = {c.id: c for c in called}
    if ctrl["component"] not in by_id:
        raise ValueError(f"the positive-control component {ctrl['component']} is not in the called set (a hub, or fewer than two such records)")
    control = [by_id[i] for i in dict.fromkeys(ids) if i in by_id]
    keep = {c.id for c in control}
    return [c for c in called if c.id not in keep], control


def copy_sources(rows: list[dict]) -> set:
    """the gorilla records that are the source of a copy row (class extra_copy, dropped_M1 or dropped_M2)"""
    return {r["source_id"] for r in rows if r["cls"] in COPY_CLASSES}


def eligible_for_planting(rows: list[dict], recs: dict, pool: Pool) -> list[str]:
    """in_place records with an exon union >= 200 bp that have no copy row as source and are in no pool component (neither as source nor as the target of
    another record's copy row)"""
    src = copy_sources(rows)
    ok = [r["source_id"] for r in rows if r["cls"] == "in_place" and int(r["exonic_len"]) >= SHORT_BP and r["source_id"] not in src
          and r["source_id"] not in pool.node_component]
    return sorted(set(ok))


def planting_counts(rows: list[dict], pool: Pool) -> dict:
    """how the pool of candidates for the planted negative control shrinks: in_place records with an exon union >= 200 bp; those without a copy row as
    source; of those, the targets of another record's copy row (the note names them); those left once the pool components are removed"""
    src = copy_sources(rows)
    targets = set()
    for r in rows:
        if r["cls"] in COPY_CLASSES:
            m = OVL.search(r["note"])
            if m:
                targets.add(m.group(1))
    ge = {r["source_id"] for r in rows if r["cls"] == "in_place" and int(r["exonic_len"]) >= SHORT_BP}
    no_src = ge - src
    return dict(in_place_ge200=len(ge), no_copy_row_as_source=len(no_src), target_of_a_copy_row=len(no_src & targets),
                eligible=len({x for x in no_src if x not in pool.node_component}))


def plant_negative(rows: list[dict], recs: dict, pool: Pool, seed: int = SEED_T2, n_records: int = 100) -> list[Component]:
    """The planted negative control: n_records eligible records drawn with random.Random(seed).sample over the sorted ids and paired at random on DIFFERENT
    contigs into n_records / 2 fake two-record components (NEG01..); they are lifted by the same recipe (runs ctl1, ctl2); at most 1 may be called."""
    el = eligible_for_planting(rows, recs, pool)
    if len(el) < n_records:
        raise ValueError(f"only {len(el)} eligible records for {n_records} planted ones")
    rng = random.Random(seed)
    chosen = rng.sample(el, n_records)
    left = list(chosen)
    out = []
    first = {r["source_id"]: r for r in rows if r["cls"] in PRIMARY_CLASSES}
    while left:
        a = left.pop(0)
        j = next((k for k, b in enumerate(left) if (recs[b] if b in recs else _record_from_row(first[b]))["contig"] !=
                  (recs[a] if a in recs else _record_from_row(first[a]))["contig"]), None)
        if j is None:
            raise ValueError("the last planted records share a contig: choose another seed in an Amendment")
        b = left.pop(j)
        comp = Component(f"NEG{len(out) + 1:02d}")
        for x in sorted((a, b)):
            r = dict(recs[x]) if x in recs else _record_from_row(first[x])
            r["cls"] = "in_place"
            comp.recs[x] = r
            comp.record_ids.append(x)
        out.append(comp)
    return out


def negative_counts(planted: list[Component], calls: dict, res: dict) -> dict:
    """what the planted negative control actually saw, for the reader (no threshold is attached to these numbers): how many of the planted records got a triple in
    all three outgroups, how many planted pairs have a triple for both members, and how many planted classes were dropped for multiplicity (a dropped pair is not an
    expansion and so not counted as called)"""
    such = [(c, c.such()) for c in planted]
    return dict(records_with_triple=sum(1 for _, rs in such for r in rs if calls[r["id"]]["triple"] is not None),
                pairs_with_both_triples=sum(1 for _, rs in such if rs and all(calls[r["id"]]["triple"] is not None for r in rs)),
                dropped_pairs=len(res["dropped_classes"]))


def check_negative(res: dict, n_components: int, bar: int = 1) -> dict:
    called = len(res["expansions"])
    return {"called": called, "of": n_components, "bar": bar, "pass": called <= bar}


# ----------------------------------------------------------------------------------------------- the input-only report and the CLI
# The registered numbers (section 7, step 1 and 2, after the errata of 2026-10-07): the hub rule counts the annotated records of any length AND the copy@
# sites of the union-find component ("nodes"). The older counting, annotated records only, is the "records" reading; two components with 74 and 52 copy@
# sites separate them.
EXPECTED = dict(components_size_ge2=609, annotated_records_in_components=2976, annotated_ge200bp=2249, components_ge2_such=307,
                hubs=7, called=300, records=1403, by_rtype={"gene": 907, "pseudogene": 496}, J=44, records_per_run_first5=[300, 300, 165, 116, 81],
                runs_with_at_most_10=26, supported={"gorilla_KB3781": 80, "gorilla_OR6737": 139, "either": 157},
                supported_hubs={"gorilla_KB3781": 4, "gorilla_OR6737": 5})
EXPECTED_RECORDS_READING = dict(EXPECTED, hubs=5, called=302, records=1446, by_rtype={"gene": 944, "pseudogene": 502},
                                records_per_run_first5=[302, 302, 167, 118, 83])


def read_support(path) -> dict:
    """reads_capped of the placed records of a support table (figures/_liftoff.py: >= 2 primary reads with an aligned block on the exon union)"""
    d = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cls"] in PLACED_CLASSES:
                d[r["source_id"]] = int(r["reads_capped"])
    return d


def pool_report(lift_root=None, support_libs=("gorilla_KB3781", "gorilla_OR6737")) -> dict:
    """input-only: the pool, the hubs under both readings of the hub rule, the called set, the rank plan and the supported components.
    `lift_root` None = LIFT_ROOT as it is when the function is called (not as it was when the module was imported)"""
    lift_root = LIFT_ROOT if lift_root is None else lift_root
    rows = read_selflift(os.path.join(lift_root, "gorilla", "liftoff_loci.tsv"))
    recs = read_records(os.path.join(lift_root, "gorilla", "records.tsv"))
    pool = build_pool(rows, recs)
    comps = list(pool.components.values())
    comp2 = [c for c in comps if len(c.such()) >= 2]
    rep = dict(copy_rows=pool.n_copy_rows, components_size_ge2=len(comps), annotated_records_in_components=sum(len(c.record_ids) for c in comps),
               annotated_ge200bp=sum(1 for c in comps for r in c.recs.values() if r["exonic_len"] >= SHORT_BP), components_ge2_such=len(comp2),
               unannotated_sites=sum(len(c.sites) for c in comps), planting=planting_counts(rows, pool))
    for mode in ("nodes", "records"):
        hubs = [c for c in comp2 if c.hub(mode)]
        called = [c for c in comp2 if not c.hub(mode)]
        plan = rank_plan(called)
        by_rtype: collections.Counter = collections.Counter(r["rtype"] for c in called for r in c.such())
        rep[mode] = dict(hubs=len(hubs), hub_sizes=sorted(c.size(mode) for c in hubs), called=len(called), records=sum(by_rtype.values()),
                         by_rtype=dict(by_rtype), J=len(plan), plan=plan_summary(plan))
    sup = {}
    for lib in support_libs:
        sup[lib] = read_support(os.path.join(lift_root, "support", lib + ".tsv"))
    rep["supported"] = {}
    if sup:
        def n_supported(c, libs):
            return sum(1 for r in c.such() if any(sup[l].get(r["id"], 0) >= 2 for l in libs))
        for name, libs in list(((l, [l]) for l in sup)) + [("either", list(sup))]:
            sc = [c for c in comp2 if n_supported(c, libs) >= 2]
            rep["supported"][name] = dict(components=len(sc), hubs_nodes=sum(1 for c in sc if c.hub("nodes")), hubs_records=sum(1 for c in sc if c.hub("records")))
    return rep


def compare_expected(rep: dict, mode: str) -> list[str]:
    """differences between a report and the numbers quoted in the pre-registration (section 7, step 1 and 2): under the hub reading 'nodes' (the registered
    rule) against EXPECTED, under 'records' (the older counting, for reference) against EXPECTED_RECORDS_READING"""
    exp = EXPECTED if mode == "nodes" else EXPECTED_RECORDS_READING
    diffs = []
    for k in ("components_size_ge2", "annotated_records_in_components", "annotated_ge200bp", "components_ge2_such"):
        if rep[k] != exp[k]:
            diffs.append(f"{k}: expected {exp[k]}, got {rep[k]}")
    m = rep[mode]
    for k in ("hubs", "called", "records", "by_rtype", "J"):
        if m[k] != exp[k]:
            diffs.append(f"{k}: expected {exp[k]}, got {m[k]}")
    if m["plan"]["records_per_run"][:5] != exp["records_per_run_first5"]:
        diffs.append(f"records_per_run_first5: expected {exp['records_per_run_first5']}, got {m['plan']['records_per_run'][:5]}")
    if m["plan"]["runs_with_at_most_10"] != exp["runs_with_at_most_10"]:
        diffs.append(f"runs_with_at_most_10: expected {exp['runs_with_at_most_10']}, got {m['plan']['runs_with_at_most_10']}")
    for lib, n in exp["supported"].items():
        got = rep.get("supported", {}).get(lib)
        if got is not None and got["components"] != n:
            diffs.append(f"supported {lib}: expected {n} components, got {got['components']}")
        if got is not None and lib in exp["supported_hubs"] and got["hubs_" + mode] != exp["supported_hubs"][lib]:
            diffs.append(f"supported {lib}: expected {exp['supported_hubs'][lib]} hubs, got {got['hubs_' + mode]}")
    return diffs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["pool"])
    ap.add_argument("--lift-root", default=LIFT_ROOT)
    ap.add_argument("--out", default=None, help="write the report as JSON here")
    a = ap.parse_args(argv)
    rep = pool_report(a.lift_root)
    rep["compare"] = {mode: compare_expected(rep, mode) for mode in ("nodes", "records")}
    text = json.dumps(rep, indent=1, sort_keys=True)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
