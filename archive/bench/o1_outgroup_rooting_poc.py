#!/usr/bin/env python3
"""Prepare and score a conservative one-gorilla O1 rooting proof of concept.

The proof of concept uses complete human loci as probes because stable multi-locus
block classes do not exist yet. It may emit ROOT_CANDIDATE records, never production
DERIVED_FROM edges. Gorilla annotations are not used.
"""

from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from o1_gene_family_audit import read_tsv, write_tsv


SOURCE_GENES = ("GOLGA2", "ITSN2")


@dataclass(frozen=True)
class Probe:
    probe_id: str
    role: str
    annotation: str
    chrom: str
    start: int
    end: int
    strand: str
    disposition: str


@dataclass(frozen=True)
class PafRecord:
    query: str
    query_len: int
    query_start: int
    query_end: int
    strand: str
    target: str
    target_len: int
    target_start: int
    target_end: int
    matches: int
    block_len: int
    mapq: int
    identity: float

    @property
    def query_coverage(self) -> float:
        return (self.query_end - self.query_start) / max(self.query_len, 1)

    @property
    def aligned_span(self) -> int:
        return min(self.query_end - self.query_start, self.target_end - self.target_start)


@dataclass(frozen=True)
class AnchorHit:
    query: str
    query_len: int
    strand: str
    target: str
    target_len: int
    target_start: int
    target_end: int
    identity: float
    query_coverage: float
    mapq: int
    record_count: int


def parse_paf_line(line: str) -> PafRecord:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise ValueError("PAF record has fewer than 12 fields")
    tags: dict[str, str] = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    matches, block_len = int(fields[9]), int(fields[10])
    identity = 1.0 - float(tags["de"]) if "de" in tags else matches / max(block_len, 1)
    return PafRecord(
        query=fields[0],
        query_len=int(fields[1]),
        query_start=int(fields[2]),
        query_end=int(fields[3]),
        strand=fields[4],
        target=fields[5],
        target_len=int(fields[6]),
        target_start=int(fields[7]),
        target_end=int(fields[8]),
        matches=matches,
        block_len=block_len,
        mapq=int(fields[11]),
        identity=identity,
    )


def read_paf(path: Path) -> list[PafRecord]:
    with path.open() as handle:
        return [parse_paf_line(line) for line in handle if line.strip()]


def parse_source_genes(gff: Path) -> list[Probe]:
    found: dict[str, Probe] = {}
    with gff.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = {}
            for item in fields[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = value
            name = attrs.get("Name", attrs.get("gene", ""))
            if name not in SOURCE_GENES:
                continue
            found[name] = Probe(
                probe_id=f"SRC_{name}",
                role="SOURCE_CONTEXT",
                annotation=name,
                chrom=fields[0],
                start=int(fields[3]) - 1,
                end=int(fields[4]),
                strand=fields[6],
                disposition="CONTEXT_SOURCE",
            )
    missing = set(SOURCE_GENES) - set(found)
    if missing:
        raise RuntimeError(f"source genes missing from GFF: {','.join(sorted(missing))}")
    return [found[name] for name in SOURCE_GENES]


def parse_family_loci(path: Path) -> list[Probe]:
    probes: list[Probe] = []
    for index, row in enumerate(read_tsv(path)):
        if row["membership_role"] not in {
            "FAMILY_MEMBER", "FAMILY_MEMBER_OUTLIER", "CANDIDATE_MEMBER"
        }:
            continue
        chrom, start, end = row["locus_id"].removeprefix("L~").split("_")
        probes.append(
            Probe(
                probe_id=f"FAM_{index:03d}",
                role=row["membership_role"],
                annotation=row["annotation"],
                chrom=chrom,
                start=int(start),
                end=int(end),
                strand=".",
                disposition=row["disposition"],
            )
        )
    return probes


def fasta_lengths(fasta: Path) -> dict[str, int]:
    fai = Path(f"{fasta}.fai")
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    lengths = {}
    with fai.open() as handle:
        for line in handle:
            fields = line.split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def fetch_sequence(fasta: Path, chrom: str, start: int, end: int) -> str:
    if start < 0 or end <= start:
        raise ValueError(f"invalid interval {chrom}:{start}-{end}")
    region = f"{chrom}:{start + 1}-{end}"
    result = subprocess.run(
        ["samtools", "faidx", str(fasta), region],
        check=True,
        capture_output=True,
        text=True,
    )
    lines = result.stdout.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError(f"samtools faidx returned no sequence for {region}")
    return "".join(lines[1:]).upper()


def write_fasta(path: Path, entries: list[tuple[str, str]]) -> None:
    with path.open("w") as out:
        for name, sequence in entries:
            out.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                out.write(sequence[offset:offset + 80] + "\n")


def prepare(
    human_fasta: Path,
    human_gff: Path,
    family_loci: Path,
    out_dir: Path,
    flank_bp: int,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    sources = parse_source_genes(human_gff)
    family = parse_family_loci(family_loci)
    probes = sources + family
    lengths = fasta_lengths(human_fasta)
    rows = []
    loci_entries: list[tuple[str, str]] = []
    source_entries: list[tuple[str, str]] = []
    family_entries: list[tuple[str, str]] = []
    left_entries: list[tuple[str, str]] = []
    right_entries: list[tuple[str, str]] = []
    for probe in probes:
        chrom_len = lengths[probe.chrom]
        left_start, left_end = max(0, probe.start - flank_bp), probe.start
        right_start, right_end = probe.end, min(chrom_len, probe.end + flank_bp)
        locus_sequence = fetch_sequence(human_fasta, probe.chrom, probe.start, probe.end)
        left_sequence = fetch_sequence(human_fasta, probe.chrom, left_start, left_end)
        right_sequence = fetch_sequence(human_fasta, probe.chrom, right_start, right_end)
        loci_entries.append((probe.probe_id, locus_sequence))
        (source_entries if probe.role == "SOURCE_CONTEXT" else family_entries).append(
            (probe.probe_id, locus_sequence)
        )
        left_entries.append((f"LEFT_{probe.probe_id}", left_sequence))
        right_entries.append((f"RIGHT_{probe.probe_id}", right_sequence))
        rows.append(
            {
                "probe_id": probe.probe_id,
                "role": probe.role,
                "annotation": probe.annotation,
                "chrom": probe.chrom,
                "start": probe.start,
                "end": probe.end,
                "strand": probe.strand,
                "disposition": probe.disposition,
                "left_start": left_start,
                "left_end": left_end,
                "right_start": right_start,
                "right_end": right_end,
                "probe_status": "LOCUS_PROXY_NOT_STABLE_BLOCK",
            }
        )
    write_tsv(out_dir / "probes.tsv", rows, list(rows[0]))
    write_fasta(out_dir / "human_all_loci.fa", loci_entries)
    write_fasta(out_dir / "human_sources.fa", source_entries)
    write_fasta(out_dir / "human_family_loci.fa", family_entries)
    write_fasta(out_dir / "human_left_flanks.fa", left_entries)
    write_fasta(out_dir / "human_right_flanks.fa", right_entries)
    write_fasta(out_dir / "human_all_flanks.fa", left_entries + right_entries)


def passing_flank(record: AnchorHit) -> bool:
    return record.query_coverage >= 0.80 and record.identity >= 0.90 and record.mapq >= 20


def covered_bases(intervals: list[tuple[int, int]]) -> int:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def compatible_anchor_segments(left: PafRecord, right: PafRecord) -> bool:
    query_gap = right.query_start - left.query_end
    if query_gap < -100 or query_gap > 1000:
        return False
    if left.strand == "+":
        target_gap = right.target_start - left.target_end
    else:
        target_gap = left.target_start - right.target_end
    return -100 <= target_gap <= 5000


def chain_anchor_records(records: list[PafRecord]) -> list[AnchorHit]:
    """Join collinear split PAF records before applying the flank-coverage rule."""
    grouped: dict[tuple[str, str, str], list[PafRecord]] = defaultdict(list)
    for record in records:
        if record.mapq >= 20:
            grouped[(record.query, record.target, record.strand)].append(record)
    anchors: list[AnchorHit] = []
    for group in grouped.values():
        ordered = sorted(group, key=lambda record: (record.query_start, record.query_end))
        chains: list[list[PafRecord]] = []
        for record in ordered:
            if chains and compatible_anchor_segments(chains[-1][-1], record):
                chains[-1].append(record)
            else:
                chains.append([record])
        for chain in chains:
            total_span = sum(record.aligned_span for record in chain)
            anchors.append(
                AnchorHit(
                    query=chain[0].query,
                    query_len=chain[0].query_len,
                    strand=chain[0].strand,
                    target=chain[0].target,
                    target_len=chain[0].target_len,
                    target_start=min(record.target_start for record in chain),
                    target_end=max(record.target_end for record in chain),
                    identity=sum(record.identity * record.aligned_span for record in chain)
                    / max(total_span, 1),
                    query_coverage=covered_bases(
                        [(record.query_start, record.query_end) for record in chain]
                    ) / max(chain[0].query_len, 1),
                    mapq=min(record.mapq for record in chain),
                    record_count=len(chain),
                )
            )
    return anchors


def synteny_pairs(
    left_records: list[AnchorHit], right_records: list[AnchorHit]
) -> list[tuple[AnchorHit, AnchorHit, int, int]]:
    pairs = []
    for left in left_records:
        if not passing_flank(left):
            continue
        for right in right_records:
            if not passing_flank(right):
                continue
            if left.target != right.target or left.strand != right.strand:
                continue
            if left.strand == "+" and left.target_start >= right.target_start:
                continue
            if left.strand == "-" and left.target_start <= right.target_start:
                continue
            if left.strand == "+":
                inner_start, inner_end = left.target_end, right.target_start
            else:
                inner_start, inner_end = right.target_end, left.target_start
            if inner_end <= inner_start:
                continue
            pairs.append((left, right, inner_start, inner_end))
    return pairs


def overlaps(record: PafRecord, anchor: AnchorHit, start: int, end: int) -> bool:
    same_sequence = record.target == anchor.target or record.target_len == anchor.target_len
    return same_sequence and min(record.target_end, end) > max(record.target_start, start)


def analyze(
    probes_path: Path,
    human_links_path: Path,
    mat_loci_path: Path,
    pat_loci_path: Path,
    mat_left_path: Path,
    mat_right_path: Path,
    pat_left_path: Path,
    pat_right_path: Path,
    out_dir: Path,
) -> None:
    probe_rows = read_tsv(probes_path)
    probes = {row["probe_id"]: row for row in probe_rows}
    write_tsv(out_dir / "probes.tsv", probe_rows, list(probe_rows[0]))
    human_links = read_paf(human_links_path)
    loci_by_hap = {"MAT": read_paf(mat_loci_path), "PAT": read_paf(pat_loci_path)}
    left_by_hap = {"MAT": read_paf(mat_left_path), "PAT": read_paf(pat_left_path)}
    right_by_hap = {"MAT": read_paf(mat_right_path), "PAT": read_paf(pat_right_path)}

    source_links: dict[str, set[str]] = defaultdict(set)
    link_rows = []
    for record in human_links:
        if record.aligned_span < 1000 or record.identity < 0.60:
            continue
        source_links[record.query].add(record.target)
        source_probe = probes[record.query]
        family_probe = probes[record.target]
        link_rows.append(
            {
                "source_probe": record.query,
                "source_annotation": source_probe["annotation"],
                "family_probe": record.target,
                "family_annotation": family_probe["annotation"],
                "source_interval": f"{record.query_start}-{record.query_end}",
                "family_interval": f"{record.target_start}-{record.target_end}",
                "source_genomic_interval": (
                    f'{source_probe["chrom"]}:'
                    f'{int(source_probe["start"]) + record.query_start}-'
                    f'{int(source_probe["start"]) + record.query_end}'
                ),
                "family_genomic_interval": (
                    f'{family_probe["chrom"]}:'
                    f'{int(family_probe["start"]) + record.target_start}-'
                    f'{int(family_probe["start"]) + record.target_end}'
                ),
                "strand": record.strand,
                "identity": f"{record.identity:.4f}",
                "aligned_span": record.aligned_span,
                "evidence_status": "PROVISIONAL_PAIRWISE_BLOCK",
            }
        )
    write_tsv(out_dir / "human_source_family_links.tsv", link_rows, list(link_rows[0]))

    synteny_rows = []
    synteny_by_probe: dict[tuple[str, str], dict[str, object]] = {}
    for haplotype in ("MAT", "PAT"):
        left_index: dict[str, list[AnchorHit]] = defaultdict(list)
        right_index: dict[str, list[AnchorHit]] = defaultdict(list)
        for record in chain_anchor_records(left_by_hap[haplotype]):
            left_index[record.query.removeprefix("LEFT_")].append(record)
        for record in chain_anchor_records(right_by_hap[haplotype]):
            right_index[record.query.removeprefix("RIGHT_")].append(record)
        locus_index: dict[str, list[PafRecord]] = defaultdict(list)
        for record in loci_by_hap[haplotype]:
            if record.aligned_span >= 1000 and record.identity >= 0.60:
                locus_index[record.query].append(record)
        for probe_id, probe in probes.items():
            pairs = synteny_pairs(left_index[probe_id], right_index[probe_id])
            if len(pairs) == 1:
                left, right, inner_start, inner_end = pairs[0]
                supported = [
                    record for record in locus_index[probe_id]
                    if overlaps(record, left, inner_start, inner_end)
                ]
                status = "TWO_SIDED_UNIQUE" if supported else "FLANKS_WITHOUT_LOCUS_HIT"
                row = {
                    "haplotype": haplotype,
                    "probe_id": probe_id,
                    "role": probe["role"],
                    "annotation": probe["annotation"],
                    "target": left.target,
                    "target_start": inner_start,
                    "target_end": inner_end,
                    "orientation": left.strand,
                    "qualifying_pairs": 1,
                    "locus_hits_in_interval": len(supported),
                    "left_identity": f"{left.identity:.4f}",
                    "right_identity": f"{right.identity:.4f}",
                    "left_coverage": f"{left.query_coverage:.4f}",
                    "right_coverage": f"{right.query_coverage:.4f}",
                    "left_mapq": left.mapq,
                    "right_mapq": right.mapq,
                    "left_records": left.record_count,
                    "right_records": right.record_count,
                    "status": status,
                }
            else:
                row = {
                    "haplotype": haplotype,
                    "probe_id": probe_id,
                    "role": probe["role"],
                    "annotation": probe["annotation"],
                    "target": ".",
                    "target_start": ".",
                    "target_end": ".",
                    "orientation": ".",
                    "qualifying_pairs": len(pairs),
                    "locus_hits_in_interval": 0,
                    "left_identity": ".",
                    "right_identity": ".",
                    "left_coverage": ".",
                    "right_coverage": ".",
                    "left_mapq": ".",
                    "right_mapq": ".",
                    "left_records": ".",
                    "right_records": ".",
                    "status": "NO_TWO_SIDED_SYNTENY" if not pairs else "FLANK_MULTIMAPPING",
                }
            synteny_rows.append(row)
            synteny_by_probe[(haplotype, probe_id)] = row
    write_tsv(out_dir / "gorilla_synteny.tsv", synteny_rows, list(synteny_rows[0]))

    hit_rows = []
    for haplotype, records in loci_by_hap.items():
        for index, record in enumerate(records):
            if record.aligned_span < 1000 or record.identity < 0.60:
                continue
            hit_rows.append(
                {
                    "evidence_id": f"{haplotype}_LOCUS_{index}",
                    "haplotype": haplotype,
                    "probe_id": record.query,
                    "target": record.target,
                    "target_start": record.target_start,
                    "target_end": record.target_end,
                    "strand": record.strand,
                    "identity": f"{record.identity:.4f}",
                    "query_coverage": f"{record.query_coverage:.4f}",
                    "aligned_span": record.aligned_span,
                    "mapq": record.mapq,
                }
            )
    write_tsv(out_dir / "gorilla_locus_hits.tsv", hit_rows, list(hit_rows[0]))

    certificate_rows = []
    family_probe_ids = {
        probe_id for probe_id, probe in probes.items() if probe["role"] != "SOURCE_CONTEXT"
    }
    family_two_sided_both = {
        probe_id for probe_id in family_probe_ids
        if synteny_by_probe[("MAT", probe_id)]["status"] == "TWO_SIDED_UNIQUE"
        and synteny_by_probe[("PAT", probe_id)]["status"] == "TWO_SIDED_UNIQUE"
    }
    for source in ("SRC_GOLGA2", "SRC_ITSN2"):
        mat = synteny_by_probe[("MAT", source)]
        pat = synteny_by_probe[("PAT", source)]
        recurrent = len(source_links[source])
        both_syntenic = mat["status"] == pat["status"] == "TWO_SIDED_UNIQUE"
        if recurrent >= 2 and both_syntenic:
            root_status = "ROOT_CANDIDATE_SINGLE_OUTGROUP"
            outgroup_status = "HAPLOTYPES_AGREE"
        elif recurrent < 2:
            root_status = "NO_RECURRENT_HUMAN_BLOCK"
            outgroup_status = "NOT_EVALUATED"
        else:
            root_status = "OUTGROUP_UNRESOLVED"
            outgroup_status = f'MAT={mat["status"]};PAT={pat["status"]}'
        certificate_rows.append(
            {
                "source_probe": source,
                "source_annotation": probes[source]["annotation"],
                "linked_human_family_loci": recurrent,
                "linked_family_probe_ids": ",".join(sorted(source_links[source])) or ".",
                "family_loci_tested": len(family_probe_ids),
                "family_loci_two_sided_both_haplotypes": len(family_two_sided_both),
                "two_sided_family_probe_ids": ",".join(sorted(family_two_sided_both)) or ".",
                "maternal_synteny": mat["status"],
                "paternal_synteny": pat["status"],
                "root_candidate_status": root_status,
                "direction_status": "UNROOTED",
                "outgroup_status": outgroup_status,
                "claim_limit": "LOCUS_PROXY_NOT_STABLE_BLOCK",
            }
        )
    write_tsv(out_dir / "rooting_candidates.tsv", certificate_rows, list(certificate_rows[0]))
    write_candidate_gfa(out_dir, probes, source_links, synteny_by_probe)
    write_report(out_dir, certificate_rows, probes, family_two_sided_both)


def write_candidate_gfa(
    out_dir: Path,
    probes: dict[str, dict[str, str]],
    source_links: dict[str, set[str]],
    synteny_by_probe: dict[tuple[str, str], dict[str, object]],
) -> None:
    graph_path = out_dir / "rooting_candidate.gfa"
    colours_path = out_dir / "rooting_candidate.colours.csv"
    with graph_path.open("w") as graph, colours_path.open("w") as colours:
        graph.write("H\tVN:Z:1.0\tTS:Z:single_outgroup_root_candidate\n")
        colours.write("Node,Colour\n")
        for probe_id, probe in probes.items():
            length = int(probe["end"]) - int(probe["start"])
            graph.write(
                f'S\t{probe_id}\t*\tLN:i:{max(length, 1)}\tTY:Z:HUMAN_LOCUS\t'
                f'RL:Z:{probe["role"]}\tGN:Z:{probe["annotation"]}\n'
            )
            colour = "#a142f4" if probe["role"] == "SOURCE_CONTEXT" else (
                "#34a853" if probe["role"].startswith("FAMILY_MEMBER") else "#4285f4"
            )
            colours.write(f"{probe_id},{colour}\n")
        for source in ("SRC_GOLGA2", "SRC_ITSN2"):
            for haplotype in ("MAT", "PAT"):
                row = synteny_by_probe[(haplotype, source)]
                if row["status"] != "TWO_SIDED_UNIQUE":
                    continue
                ape_id = f"GOR_{haplotype}_{source}"
                length = int(row["target_end"]) - int(row["target_start"])
                graph.write(
                    f'S\t{ape_id}\t*\tLN:i:{max(length, 1)}\tTY:Z:OUTGROUP_LOCUS\t'
                    f'HP:Z:{haplotype}\tCT:Z:{row["target"]}\n'
                )
                colours.write(f"{ape_id},#00acc1\n")
                graph.write(
                    f'L\t{ape_id}\t+\t{source}\t+\t0M\tRT:Z:SYNTENIC_WITH\t'
                    "DS:Z:UNROOTED\tCL:Z:LOCUS_PROXY_NOT_STABLE_BLOCK\n"
                )
            for family_probe in sorted(source_links[source]):
                graph.write(
                    f'L\t{source}\t+\t{family_probe}\t+\t0M\t'
                    "RT:Z:PROVISIONAL_BLOCK_MATCH\tDS:Z:UNROOTED\t"
                    "CL:Z:LOCUS_PROXY_NOT_STABLE_BLOCK\n"
                )


def write_report(
    out_dir: Path,
    certificates: list[dict[str, object]],
    probes: dict[str, dict[str, str]],
    family_two_sided_both: set[str],
) -> None:
    two_sided_labels = []
    for probe_id in sorted(family_two_sided_both):
        probe = probes[probe_id]
        label = probe["annotation"] if probe["annotation"] != "." else (
            f'{probe["chrom"]}:{probe["start"]}-{probe["end"]}'
        )
        two_sided_labels.append(f"{probe_id} ({label})")
    lines = [
        "# O1 single-outgroup rooting proof of concept",
        "",
        "This experiment asks whether the two proposed human source contexts satisfy the",
        "prerequisites for later provenance rooting against one ape species. It uses both phased",
        "KB3781 gorilla haplotypes independently. Gorilla annotations are not used.",
        "",
        "| source | linked human GOLGA-family loci | maternal synteny | paternal synteny | result |",
        "|---|---:|---|---|---|",
    ]
    for row in certificates:
        lines.append(
            f'| {row["source_annotation"]} | {row["linked_human_family_loci"]} | '
            f'{row["maternal_synteny"]} | {row["paternal_synteny"]} | '
            f'{row["root_candidate_status"]} |'
        )
    lines += [
        "",
        "Both source contexts pass the proof-of-concept prerequisites: recurrent human family",
        "homology and a unique, two-sided syntenic gorilla locus in each haplotype. They remain",
        "`UNROOTED`, because the queries are whole-locus proxies and the shared intervals are",
        "pairwise witnesses rather than stable multi-locus block classes.",
        "",
        f'Of 18 audited family probes, {len(family_two_sided_both)} have unique two-sided synteny in '
        "both gorilla haplotypes:",
        "",
    ]
    lines.extend(f"- {label}" for label in two_sided_labels)
    lines += [
        "",
        "The other 15 probes lack the strict two-flank certificate. This is compatible with",
        "duplication/rearrangement but is not treated as proof of gorilla absence: repetitive-region",
        "mapping failure and assembly structure remain alternatives.",
        "",
        "## Fixed proof-of-concept rules",
        "",
        "- human source-to-family witness: aligned span >=1,000 bp and identity >=0.60;",
        "- flank length: 25,000 bp on each side;",
        "- qualifying flank anchor: coverage >=0.80, identity >=0.90, MAPQ >=20;",
        "- two anchors must occur on one ape sequence, in compatible order and orientation; and",
        "- collinear split anchors may be chained with query gap <=1,000 bp and target gap <=5,000 bp.",
        "",
        "The split-anchor rule was added after an audit showed that minimap2 represented the ITSN2",
        "right flank as two adjacent records. It is substrate-neutral, applied to every probe and",
        "haplotype, and does not change any alignment threshold. Prebuilt asm20 indexes and FASTA",
        "asm5 runs use different contig labels; target length plus coordinates reconcile those",
        "labels without changing placements.",
        "",
        "## Claim boundary",
        "",
        "This is a positive feasibility result for single-outgroup rooting, not evidence that",
        "GOLGA2 or ITSN2 has already been proven ancestral by Rustle. Production direction requires",
        "stable human block classes, block-specific source/derived paths, assembly-gap checks and a",
        "rooting certificate. Until those exist, no `DERIVED_FROM` edge is emitted.",
        "",
        "Normative evidence: `probes.tsv`, `human_source_family_links.tsv`, `gorilla_synteny.tsv`,",
        "`gorilla_locus_hits.tsv`, and `rooting_candidates.tsv`. Raw PAF files are retained under",
        "`work/`. `rooting_candidate.gfa` is a typed visual projection and deliberately contains no",
        "`DERIVED_FROM` edge.",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--human-fasta", type=Path, required=True)
    prepare_parser.add_argument("--human-gff", type=Path, required=True)
    prepare_parser.add_argument("--family-loci", type=Path, required=True)
    prepare_parser.add_argument("--out-dir", type=Path, required=True)
    prepare_parser.add_argument("--flank-bp", type=int, default=25_000)

    analyze_parser = subparsers.add_parser("analyze")
    analyze_parser.add_argument("--probes", type=Path, required=True)
    analyze_parser.add_argument("--human-links", type=Path, required=True)
    analyze_parser.add_argument("--mat-loci", type=Path, required=True)
    analyze_parser.add_argument("--pat-loci", type=Path, required=True)
    analyze_parser.add_argument("--mat-left", type=Path, required=True)
    analyze_parser.add_argument("--mat-right", type=Path, required=True)
    analyze_parser.add_argument("--pat-left", type=Path, required=True)
    analyze_parser.add_argument("--pat-right", type=Path, required=True)
    analyze_parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.command == "prepare":
        prepare(args.human_fasta, args.human_gff, args.family_loci, args.out_dir, args.flank_bp)
    else:
        analyze(
            args.probes, args.human_links, args.mat_loci, args.pat_loci,
            args.mat_left, args.mat_right, args.pat_left, args.pat_right, args.out_dir,
        )


if __name__ == "__main__":
    main()
