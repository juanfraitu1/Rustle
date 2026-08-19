#!/usr/bin/env python3
"""Evaluate a transcription-orientation guard for RNA O1's E_r edge relation.

Rustle's RNA O1 exon-sum and RNA-locus genomic-span sequences are normalized into
transcript (5'->3') orientation.  Consequently a qualifying PAF record on strand
'-' is reverse-complement homology, not colinear transcript homology.  This script
scores the rule

    E_r^+(a,b) iff at least one existing E_r-qualifying record has PAF strand '+'

on the frozen O1 TP/FP arms.  The historical seven-family panel is rebuilt first,
because its Python sequence builder omitted Rustle's minus-locus reverse complement.
The correction changes orientation only, never sequence content or E_r membership.

This does NOT apply to ``gw_family_catalog --from-genome``: those read-free DNA
sequences deliberately remain in genomic-plus orientation, where '-' is a valid
inverted duplication.  This is an offline hypothesis test, not a production run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import pickle
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


FAMILIES = ("NPIP", "TBC1D3", "RABL2", "APOBEC3", "MAGEA", "GSTM", "HERC2")
REF_CIGAR = {"M", "=", "X", "D", "N"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name: str | None = None
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {name: "".join(parts) for name, parts in seqs.items()}


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]


def longest_forward_orf(seq: str) -> tuple[int, int]:
    """Coordinates of the longest forward-frame ATG-to-stop/end ORF."""
    best = (0, 0)
    for frame in range(3):
        start = None
        for pos in range(frame, len(seq) - 2, 3):
            codon = seq[pos : pos + 3].upper()
            if start is None:
                if codon == "ATG":
                    start = pos
            elif codon in STOP_CODONS:
                if pos + 3 - start > best[1] - best[0]:
                    best = (start, pos + 3)
                start = None
        if start is not None and len(seq) - start > best[1] - best[0]:
            best = (start, len(seq))
    return best


def overlap_fraction(start: int, end: int, interval: tuple[int, int]) -> float:
    return max(0, min(end, interval[1]) - max(start, interval[0])) / max(end - start, 1)


def reference_end(pos0: int, cigar: str) -> int:
    pos = pos0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in REF_CIGAR:
            pos += int(length)
    return pos


def cached_sam(cache: Path, bam: Path, chrom: str, start: int, end: int) -> str:
    region = f"{chrom}:{start + 1}-{end}"
    key = hashlib.sha1(f"{bam}|{region}".encode()).hexdigest()[:20]
    cached = cache / f"{key}.pkl"
    if cached.exists():
        with cached.open("rb") as fh:
            return pickle.load(fh)
    proc = subprocess.run(
        ["samtools", "view", "-F", "2308", str(bam), region],
        check=True,
        capture_output=True,
        text=True,
    )
    return proc.stdout


def dominant_transcript_strand(sam: str, start: int, end: int) -> tuple[str, int, int]:
    forward = reverse = 0
    for line in sam.splitlines():
        fields = line.split("\t", 6)
        if len(fields) < 6 or "N" not in fields[5]:
            continue
        pos0 = int(fields[3]) - 1
        if pos0 < start or reference_end(pos0, fields[5]) > end:
            continue
        if int(fields[1]) & 16:
            reverse += 1
        else:
            forward += 1
    if forward == reverse:
        raise RuntimeError(f"no dominant transcript strand ({forward}+/{reverse}-)")
    return ("-" if reverse > forward else "+"), forward, reverse


def rebuild_panel(
    panel_root: Path, bam: Path, minimap2: str, work: Path
) -> dict[str, Path]:
    cache = panel_root / "work" / ".cache"
    work.mkdir(parents=True, exist_ok=True)
    pafs: dict[str, Path] = {}
    for family in FAMILIES:
        rows: dict[str, tuple[str, int, int]] = {}
        with (panel_root / "nodes" / f"{family}.tsv").open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rows[row["node"]] = (row["chrom"], int(row["start0"]), int(row["end"]))
        seqs = read_fasta(panel_root / "work" / f"{family}.rnaseq.fa")
        corrected: dict[str, str] = {}
        strand_counts = Counter()
        for name, seq in seqs.items():
            chrom, start, end = rows[name]
            strand, forward, reverse = dominant_transcript_strand(
                cached_sam(cache, bam, chrom, start, end), start, end
            )
            strand_counts[strand] += 1
            corrected[name] = reverse_complement(seq) if strand == "-" else seq
        fasta = work / f"{family}.transcript_oriented.fa"
        with fasta.open("w") as fh:
            for name, seq in corrected.items():
                fh.write(f">{name}\n{seq}\n")
        paf = work / f"{family}.transcript_oriented.paf"
        with paf.open("w") as out:
            subprocess.run(
                [minimap2, "-c", "-X", "--no-long-join", "-t", "4", "-k", "11", "-w", "5", str(fasta), str(fasta)],
                check=True,
                stdout=out,
                stderr=subprocess.DEVNULL,
            )
        pafs[family] = paf
        print(
            f"panel_rebuild\t{family}\tsequences={len(corrected)}\t"
            f"plus={strand_counts['+']}\tminus_rc={strand_counts['-']}"
        )
    return pafs


def load_scored(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r["scored_core"] == "1"]


def qualifying_strands(path: Path, wanted: set[tuple[str, str]]) -> dict[tuple[str, str], set[str]]:
    found: dict[tuple[str, str], set[str]] = defaultdict(set)
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[0] == f[5]:
                continue
            pair = tuple(sorted((f[0], f[5])))
            if pair not in wanted:
                continue
            qlen, qstart, qend = int(f[1]), int(f[2]), int(f[3])
            tlen, tstart, tend = int(f[6]), int(f[7]), int(f[8])
            de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
            identity = 1.0 - de if de is not None else int(f[9]) / max(int(f[10]), 1)
            span = qend - qstart if qlen <= tlen else tend - tstart
            coverage = span / max(min(qlen, tlen), 1)
            if identity >= 0.60 and coverage >= 0.50:
                found[pair].add(f[4])
    return found


def qualifying_edges(path: Path) -> tuple[set[tuple[str, str]], set[tuple[str, str]]]:
    """Return (all qualifying pairs, pairs with a qualifying '+' record)."""
    found: dict[tuple[str, str], set[str]] = defaultdict(set)
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[0] == f[5]:
                continue
            pair = tuple(sorted((f[0], f[5])))
            qlen, qstart, qend = int(f[1]), int(f[2]), int(f[3])
            tlen, tstart, tend = int(f[6]), int(f[7]), int(f[8])
            de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
            identity = 1.0 - de if de is not None else int(f[9]) / max(int(f[10]), 1)
            span = qend - qstart if qlen <= tlen else tend - tstart
            coverage = span / max(min(qlen, tlen), 1)
            if identity >= 0.60 and coverage >= 0.50:
                found[pair].add(f[4])
    return set(found), {pair for pair, pair_strands in found.items() if "+" in pair_strands}


def qualifying_coding_edges(path: Path, sequences: dict[str, str]) -> set[tuple[str, str]]:
    """C1: R0 plus transcript orientation and an ORF-intersection repeat guard.

    A record survives when it is forward and either identity is at least 0.90 or at
    least 20% of its aligned block on one side intersects that sequence's longest
    forward ORF. Pair semantics remain ANY-surviving-record.
    """
    orfs = {name: longest_forward_orf(seq) for name, seq in sequences.items()}
    kept = set()
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[0] == f[5] or f[0] not in orfs or f[5] not in orfs:
                continue
            qlen, qstart, qend = int(f[1]), int(f[2]), int(f[3])
            tlen, tstart, tend = int(f[6]), int(f[7]), int(f[8])
            de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
            identity = 1.0 - de if de is not None else int(f[9]) / max(int(f[10]), 1)
            span = qend - qstart if qlen <= tlen else tend - tstart
            coverage = span / max(min(qlen, tlen), 1)
            orf_overlap = max(
                overlap_fraction(qstart, qend, orfs[f[0]]),
                overlap_fraction(tstart, tend, orfs[f[5]]),
            )
            if (
                identity >= 0.60
                and coverage >= 0.50
                and f[4] == "+"
                and (identity >= 0.90 or orf_overlap >= 0.20)
            ):
                kept.add(tuple(sorted((f[0], f[5]))))
    return kept


def components(nodes: set[str], edges: set[tuple[str, str]]) -> list[set[str]]:
    adjacent: dict[str, set[str]] = {node: set() for node in nodes}
    for a, b in edges:
        if a in adjacent and b in adjacent:
            adjacent[a].add(b)
            adjacent[b].add(a)
    out = []
    unseen = set(nodes)
    while unseen:
        seed = min(unseen)
        stack = [seed]
        component = set()
        while stack:
            node = stack.pop()
            if node not in unseen:
                continue
            unseen.remove(node)
            component.add(node)
            stack.extend(adjacent[node] & unseen)
        out.append(component)
    return sorted(out, key=lambda c: (-len(c), min(c)))


def report_panel_impact(panel_pafs: dict[str, Path]) -> None:
    print(
        "\npanel\tnodes\tedges_r0\tedges_plus\tedges_coding\t"
        "components_r0\tcomponents_plus\tcomponents_coding"
    )
    for family, path in panel_pafs.items():
        sequences = read_fasta(path.with_suffix(".fa"))
        nodes = set(sequences)
        all_edges, plus_edges = qualifying_edges(path)
        coding_edges = qualifying_coding_edges(path, sequences)
        before = ",".join(str(len(c)) for c in components(nodes, all_edges))
        after_components = components(nodes, plus_edges)
        after = ",".join(str(len(c)) for c in after_components)
        coding = ",".join(str(len(c)) for c in components(nodes, coding_edges))
        print(
            f"{family}\t{len(nodes)}\t{len(all_edges)}\t{len(plus_edges)}\t{len(coding_edges)}\t"
            f"{before}\t{after}\t{coding}"
        )
        if before != after:
            for component in after_components[1:]:
                print(f"PANEL_SPLIT\t{family}\t{','.join(sorted(component))}")


def load_catalog_truth(arms: Path) -> dict[tuple[str, tuple[str, str]], tuple[str, str, str]]:
    truth = {}
    for filename, label in (
        ("fp_set.tsv", "FP"),
        ("tp_set.tsv", "TP"),
        ("grey_not_scored.tsv", "GREY"),
        ("excluded_deep_paralog.tsv", "EXCLUDED"),
    ):
        with (arms / filename).open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                pair = tuple(sorted((row["seq_a"], row["seq_b"])))
                truth[(row["source_key"], pair)] = (
                    label,
                    row.get("scored_core", "0"),
                    row["klass"],
                )
    return truth


def report_catalog_impact(root: Path, arms: Path) -> None:
    truth = load_catalog_truth(arms)
    print(
        "\ncatalog\tfamilies\tnodes\tedges_r0\tedges_plus\tedges_rejected\t"
        "families_changed\tfamilies_disconnected\tnodes_isolated"
    )
    for species in ("GGO", "HSA"):
        nodes_by_family: dict[str, set[str]] = defaultdict(set)
        with (root / f"{species}_copy_annot.tsv").open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                nodes_by_family[row["family_id"]].add(f'{row["family_id"]}~{row["copy_idx"]}')
        all_edges, plus_edges = qualifying_edges(root / f"{species}_copies.er.paf")
        within_all = {
            edge for edge in all_edges
            if edge[0].rsplit("~", 1)[0] == edge[1].rsplit("~", 1)[0]
        }
        within_plus = {
            edge for edge in plus_edges
            if edge[0].rsplit("~", 1)[0] == edge[1].rsplit("~", 1)[0]
        }
        all_by_family: dict[str, set[tuple[str, str]]] = defaultdict(set)
        plus_by_family: dict[str, set[tuple[str, str]]] = defaultdict(set)
        for edge in within_all:
            all_by_family[edge[0].rsplit("~", 1)[0]].add(edge)
        for edge in within_plus:
            plus_by_family[edge[0].rsplit("~", 1)[0]].add(edge)
        changed = disconnected = isolated = 0
        rejected_truth = Counter()
        affected_truth: dict[str, Counter] = defaultdict(Counter)
        details = []
        source = f"o1_errorcensus_{species}_copies.er"
        for edge in within_all - within_plus:
            label, scored, klass = truth.get((source, edge), ("UNLABELLED", "", ""))
            key = f"{label}_core" if scored == "1" else label
            rejected_truth[key] += 1
            affected_truth[edge[0].rsplit("~", 1)[0]][key] += 1
        for family, nodes in nodes_by_family.items():
            old = all_by_family[family]
            new = plus_by_family[family]
            if old != new:
                changed += 1
            old_components = components(nodes, old)
            new_components = components(nodes, new)
            if len(new_components) > len(old_components):
                disconnected += 1
                newly_isolated = sum(
                    len(c) == 1 and next(iter(c)) in {n for edge in old for n in edge}
                    for c in new_components
                )
                isolated += newly_isolated
                details.append(
                    (family, len(nodes), len(old), len(new),
                     ",".join(str(len(c)) for c in old_components),
                     ",".join(str(len(c)) for c in new_components))
                )
        print(
            f"{species}\t{len(nodes_by_family)}\t{sum(map(len, nodes_by_family.values()))}\t"
            f"{len(within_all)}\t{len(within_plus)}\t{len(within_all - within_plus)}\t"
            f"{changed}\t{disconnected}\t{isolated}"
        )
        print(
            f"CATALOG_REJECT_TRUTH\t{species}\t" +
            "\t".join(f"{key}={value}" for key, value in sorted(rejected_truth.items()))
        )
        for family, n, old_n, new_n, old_cs, new_cs in sorted(details):
            labels = ",".join(
                f"{key}:{value}" for key, value in sorted(affected_truth[family].items())
            )
            print(
                f"CATALOG_SPLIT\t{species}\t{family}\tn={n}\tedges={old_n}->{new_n}\t"
                f"components={old_cs}->{new_cs}\trejected_truth={labels}"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arms", type=Path, default=Path("/home/juanfra/winloci_scratch/o1_fix/arms"))
    parser.add_argument("--panel-root", type=Path, default=Path("/home/juanfra/winloci_scratch/o1_closure"))
    parser.add_argument("--bam", type=Path, default=Path("/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam"))
    parser.add_argument("--minimap2", default="/home/juanfra/miniforge3/bin/minimap2")
    parser.add_argument("--work", type=Path, default=Path("/tmp/rustle_o1_orientation"))
    parser.add_argument(
        "--catalog-root",
        type=Path,
        default=Path("/home/juanfra/winloci_scratch/o1_errorcensus"),
    )
    args = parser.parse_args()

    panel_pafs = rebuild_panel(args.panel_root, args.bam, args.minimap2, args.work)
    report_panel_impact(panel_pafs)
    report_catalog_impact(args.catalog_root, args.arms)
    groups = {"FP": load_scored(args.arms / "fp_set.tsv"), "TP": load_scored(args.arms / "tp_set.tsv")}
    wanted: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for rows in groups.values():
        for row in rows:
            wanted[row["source_key"]].add(tuple(sorted((row["seq_a"], row["seq_b"]))))

    strands: dict[tuple[str, tuple[str, str]], set[str]] = {}
    for source, pairs in wanted.items():
        path = args.arms / "records" / f"{source}.paf"
        got = qualifying_strands(path, pairs)
        strands.update(((source, pair), value) for pair, value in got.items())
    for family, path in panel_pafs.items():
        source = f"o1_fix_paf_panel7_{family}.er"
        got = qualifying_strands(path, wanted[source])
        strands.update(((source, pair), value) for pair, value in got.items())

    print("\narm\tspecies\tstratum\tn\tretain_plus\treject_minus_only")
    for arm, rows in groups.items():
        by_group: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            by_group[row["species"], row["stratum"]].append(row)
        for (species, stratum), subset in sorted(by_group.items()):
            keep = 0
            for row in subset:
                pair = tuple(sorted((row["seq_a"], row["seq_b"])))
                keep += "+" in strands.get((row["source_key"], pair), set())
            print(f"{arm}\t{species}\t{stratum}\t{len(subset)}\t{keep}\t{len(subset) - keep}")

    fp = groups["FP"]
    rejected = []
    for row in fp:
        pair = tuple(sorted((row["seq_a"], row["seq_b"])))
        if "+" not in strands.get((row["source_key"], pair), set()):
            rejected.append(row)
    mechanisms = {row["mechanism_id"] for row in rejected}
    catalog_rows = [row for row in rejected if row["shipped_catalog_same_family"] == "1"]
    print(
        f"\nFP_TOTAL\trejected={len(rejected)}/{len(fp)}\t"
        f"mechanisms={len(mechanisms)}\tcatalog_rows={len(catalog_rows)}"
    )
    for row in rejected:
        print(
            "FP_REJECT\t{species}\t{stratum}\t{gene_pair_key}\t{mechanism_id}".format(**row)
        )

    tp = groups["TP"]
    lost = []
    for row in tp:
        pair = tuple(sorted((row["seq_a"], row["seq_b"])))
        if "+" not in strands.get((row["source_key"], pair), set()):
            lost.append(row)
    print(f"\nTP_TOTAL\tlost={len(lost)}/{len(tp)}")
    for row in lost:
        print(
            "TP_LOST\t{species}\t{stratum}\t{gene_pair_key}\t{seq_a}\t{seq_b}\t"
            "{truth_evidence}".format(**row)
        )


if __name__ == "__main__":
    main()
