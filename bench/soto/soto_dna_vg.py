#!/usr/bin/env python3
"""DNA variation-graph ceiling vs Soto SDs (visual demo). Per flagship Soto family, build a base-level
variation graph (abpoa MSA -> GFA) from the member GENOMIC sequences and colour each copy green
(RNA-recovered) / red (DNA-only). The graph REPRESENTS all copies (the DNA ceiling); it does not
"detect" families. See docs/superpowers/specs/2026-07-18-soto-dna-vg-ceiling-design.md."""
import os
import sys
from collections import defaultdict

# ---- paths ----
BED = "bench/soto/80_fams.chr.bed"
MEMFA = "/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_members.fa"
DETECT = "bench/soto/soto_member_detection.tsv"
CAUSES = "bench/soto/soto_floor_decomposition.tsv"
OUT = "/home/juanfra/winloci_scratch/soto_vg"

WINDOW_BP = 20000              # memory-safe abpoa window on this 19 GB box (30000 still OOM'd; full SD spans OOM)
MM2 = "/home/juanfra/miniforge3/bin/minimap2"


def apply_window(seqs, member_spans, window_bp):
    """Pure: extract each member's homologous window. member_spans[i] = (qs, qe) coords in member i that
    align to the anchor's window, or None (no hit -> fallback first window_bp). Clamp each to window_bp."""
    out = []
    for i, s in enumerate(seqs):
        sp = member_spans[i]
        if sp is None:
            out.append(s[:window_bp])
        else:
            qs, qe = sp
            out.append(s[qs:min(qe, qs + window_bp)])
    return out


def window_members(seqs, window_bp=WINDOW_BP):
    """If any member exceeds window_bp, window ALL members to the region homologous to the SHORTEST member's
    first window_bp (via minimap2), so the graph stays a shared-backbone view that fits abpoa's memory.
    Returns (windowed_seqs, windowed: bool, window_span: int|None) — window_span is the ACTUAL window size
    (min(window_bp, len(shortest member))), or None when not windowed. No-op (seqs unchanged, False, None)
    when all fit."""
    if all(len(s) <= window_bp for s in seqs):
        return list(seqs), False, None
    import subprocess, tempfile
    anchor_idx = min(range(len(seqs)), key=lambda i: len(seqs[i]))
    anchor_win = seqs[anchor_idx][:window_bp]
    with tempfile.TemporaryDirectory() as td:
        aq, tf = f"{td}/anchor.fa", f"{td}/members.fa"
        open(aq, "w").write(f">anchor\n{anchor_win}\n")
        open(tf, "w").write("".join(f">m{i}\n{s}\n" for i, s in enumerate(seqs)))
        # ref=anchor, query=members  => PAF query is m{i}; member coords are qstart/qend (cols 2,3)
        paf = subprocess.run([MM2, "-x", "asm20", "-c", aq, tf],
                             capture_output=True, text=True, timeout=600).stdout
    best = {}
    for line in paf.splitlines():
        c = line.split("\t")
        if len(c) < 11:
            continue
        i = int(c[0][1:])                       # qname "m{i}"
        qs, qe, block = int(c[2]), int(c[3]), int(c[10])
        if i not in best or block > best[i][2]:
            best[i] = (qs, qe, block)
    spans = [(best[i][0], best[i][1]) if i in best else None for i in range(len(seqs))]
    return apply_window(seqs, spans, window_bp), True, len(anchor_win)


# ---- colours (match copy_graph.rs) ----
GREEN = "#1e8e3e"   # RNA-recovered
RED = "#d93025"     # DNA-only (K=0 / silent / coverage)
GREY = "#9aa0a6"    # shared / conserved


def parse_family_members(bed_lines, family_id):
    """(gene, chrom, start, end) for every member of family_id, in file order. BED col4 = GENE|ID_k."""
    out = []
    for ln in bed_lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        name = f[3]
        if "|" not in name:
            continue
        gene, fam = name.rsplit("|", 1)
        if fam == family_id:
            out.append((gene, f[0], int(f[1]), int(f[2])))
    return out


def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(buf)
            cur = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def member_seq(fa, chrom, start, end):
    """Genomic sequence for a BED member. soto_members.fa headers are 1-based (start+1)."""
    return fa[f"{chrom}:{start + 1}-{end}"]


def abpoa_msa(seqs):
    """Aligned gapped rows for the members (str inputs; bytes silently yield all-N). Returns the first
    len(seqs) rows of abpoa's column-MSA (excludes any appended consensus row)."""
    import pyabpoa
    aligner = pyabpoa.msa_aligner()
    res = aligner.msa([s.upper() for s in seqs], out_cons=False, out_msa=True)
    rows = list(res.msa_seq)[:len(seqs)]
    return [r.upper() for r in rows]


def msa_to_gfa(rows, names):
    """Column-MSA -> base-level variation-graph GFA. rows: equal-length uppercase gapped strings
    ('-' = gap), one per name. Returns (gfa_text, paths). Deterministic:
      - a maximal run of columns where all rows share the SAME non-gap base -> one shared node (all paths);
      - a maximal run of variant/gap columns -> one node per distinct gap-stripped allele (sorted); a
        member whose slice is all-gap skips the region."""
    assert rows and len({len(r) for r in rows}) == 1, "rows must be non-empty and equal length"
    assert len(rows) == len(names)
    m, L = len(rows), len(rows[0])

    def invariant(j):
        c0 = rows[0][j]
        return c0 != "-" and all(rows[i][j] == c0 for i in range(m))

    # segment columns into maximal same-class runs
    segments, j = [], 0
    while j < L:
        kind = "inv" if invariant(j) else "var"
        k = j + 1
        while k < L and ("inv" if invariant(k) else "var") == kind:
            k += 1
        segments.append((kind, j, k))
        j = k

    nodes, paths, nid = [], {n: [] for n in names}, 0
    for kind, a, b in segments:
        if kind == "inv":
            nid += 1
            node = str(nid)
            nodes.append((node, rows[0][a:b]))
            for n in names:
                paths[n].append(node)
        else:
            allele = {n: rows[i][a:b].replace("-", "") for i, n in enumerate(names)}
            node_of = {}
            for s in sorted(set(v for v in allele.values() if v)):
                nid += 1
                node_of[s] = str(nid)
                nodes.append((str(nid), s))
            for n in names:
                if allele[n]:
                    paths[n].append(node_of[allele[n]])

    links = set()
    for n in names:
        p = paths[n]
        for x, y in zip(p, p[1:]):
            links.add((x, y))

    out = ["H\tVN:Z:1.0"]
    for node, seq in nodes:
        out.append(f"S\t{node}\t{seq}")
    for x, y in sorted(links, key=lambda t: (int(t[0]), int(t[1]))):
        out.append(f"L\t{x}\t+\t{y}\t+\t0M")
    for n in names:
        out.append(f"P\t{n}\t{'+,'.join(paths[n])}+\t*")
    return "\n".join(out) + "\n", paths


def load_detection(tsv_lines):
    """(chrom,start,end) -> (detected: bool, recovered_by: str)."""
    out = {}
    for i, ln in enumerate(tsv_lines):
        if i == 0:
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        out[(f[2], int(f[3]), int(f[4]))] = (f[5] == "Y", f[6] if len(f) > 6 else "")
    return out


def load_causes(tsv_lines):
    """(chrom,start,end) -> miss cause, from soto_floor_decomposition.tsv (cols: 3 chrom,4 start,5 end,12 cause)."""
    out = {}
    for i, ln in enumerate(tsv_lines):
        if i == 0:
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        out[(f[2], int(f[3]), int(f[4]))] = f[11]
    return out


def node_colour(members_through_node, detected_by_name):
    """green if every member through the node is RNA-recovered, red if none is, grey otherwise."""
    flags = {bool(detected_by_name.get(n, False)) for n in members_through_node}
    if flags == {True}:
        return GREEN
    if flags == {False}:
        return RED
    return GREY


def colours_csv(paths, detected_by_name):
    node_members = defaultdict(set)
    for n, p in paths.items():
        for node in p:
            node_members[node].add(n)
    lines = ["Node,Colour"]
    for node in sorted(node_members, key=int):
        lines.append(f"{node},{node_colour(node_members[node], detected_by_name)}")
    return "\n".join(lines) + "\n"


def legend_tsv(members, detected_by_name, recovered_by_name, cause_by_name):
    lines = ["gene\tlocus\tdetected\trecovered_by\tcause\tcolour"]
    for gene, chrom, start, end in members:
        det = bool(detected_by_name.get(gene, False))
        lines.append("\t".join([
            gene, f"{chrom}:{start}-{end}", "Y" if det else "N",
            recovered_by_name.get(gene, ""),
            "" if det else cause_by_name.get(gene, ""),      # red members: their decomposition cause
            GREEN if det else RED]))
    return "\n".join(lines) + "\n"


FLAGSHIPS = [("ID_462", "SRGAP2"), ("ID_8", "PMS2P"), ("ID_131", "AMYLASE (AMY1/AMY2)")]

CAPTION = (
    "DNA variation graph = the identifiability CEILING: all {n} Soto copies present as paths "
    "(corroborated by Soto's INDEPENDENT DNA-read-depth catalog). green = RNA recovered; red = NOT "
    "recovered (per-member cause in the legend) — K=0 exon-identity copies are the identifiability floor "
    "(RNA cannot separate exon-identical copies; DNA can). The VG REPRESENTS what is given; it does not "
    "'detect' families. RNA recovers 76.2% of this ceiling genome-wide.{window}"
)


def build_family(family_id, members, fa, detection, causes):
    """Extract member seqs (skip+log those absent from fa), abpoa MSA, GFA, colours, legend, presence check."""
    present, missing, seqs, names = [], [], [], []
    for gene, chrom, start, end in members:
        try:
            seqs.append(member_seq(fa, chrom, start, end))
            names.append(gene)
            present.append((gene, chrom, start, end))
        except KeyError:
            missing.append(gene)
    if not present:
        return dict(family_id=family_id, n_members=len(members), n_present=0, gfa="", colours="",
                    legend="", missing=missing, windowed=False, window_span=None)
    seqs, windowed, window_span = window_members(seqs)
    rows = abpoa_msa(seqs) if len(seqs) > 1 else [seqs[0]]
    gfa, paths = msa_to_gfa(rows, names)
    det_by = {g: detection.get((c, s, e), (False, ""))[0] for g, c, s, e in present}
    rec_by = {g: detection.get((c, s, e), (False, ""))[1] for g, c, s, e in present}
    cause_by = {g: causes.get((c, s, e), "") for g, c, s, e in present}
    colours = colours_csv(paths, det_by)
    legend = legend_tsv(present, det_by, rec_by, cause_by)
    n_present = sum(1 for l in gfa.splitlines() if l.startswith("P\t"))   # checked, not assumed
    return dict(family_id=family_id, n_members=len(members), n_present=n_present, gfa=gfa,
                colours=colours, legend=legend, missing=missing, windowed=windowed, window_span=window_span)


def main():
    os.makedirs(OUT, exist_ok=True)
    bed = open(BED).read().splitlines()
    fa = read_fasta(MEMFA)
    detection = load_detection(open(DETECT).read().splitlines())
    causes = load_causes(open(CAUSES).read().splitlines())
    index = ["# Soto DNA variation-graph ceiling — flagship families\n"]
    for family_id, label in FLAGSHIPS:
        members = parse_family_members(bed, family_id)
        try:
            r = build_family(family_id, members, fa, detection, causes)
        except Exception as ex:
            print(f"{family_id} ({label}): RENDER FAILED — {type(ex).__name__}: {ex}")
            index.append(f"## {family_id} — {label}: RENDER FAILED ({type(ex).__name__}) — "
                         f"full genomic span exceeds this box's abpoa memory even windowed")
            continue
        base = f"{OUT}/{family_id}"
        open(f"{base}.gfa", "w").write(r["gfa"])
        open(f"{base}.colours.csv", "w").write(r["colours"])
        open(f"{base}.legend.tsv", "w").write(r["legend"])
        if r.get("windowed"):
            kb = (r["window_span"] or 0) // 1000
            wnote = (f" NOTE: this graph is WINDOWED to the ~{kb} kb region homologous to the anchor "
                     f"member — the full genomic span exceeds this box's abpoa memory.")
            wtag = f" [windowed ~{kb} kb]"
        else:
            wnote, wtag = "", ""
        miss = f"  (MISSING from graph: {r['missing']})" if r["missing"] else ""
        print(f"{family_id} ({label}): {r['n_present']}/{r['n_members']} copies as paths{wtag}{miss}")
        index.append(f"## {family_id} — {label}: {r['n_present']}/{r['n_members']} copies present as paths{wtag}")
        index.append(f"`{family_id}.gfa` + `{family_id}.colours.csv` (Bandage). " + CAPTION.format(n=r["n_present"], window=wnote))
        if r["missing"]:
            index.append(f"> honesty: {len(r['missing'])} member(s) absent from the graph: {r['missing']}")
    open(f"{OUT}/index.md", "w").write("\n\n".join(index) + "\n")
    print(f"wrote {OUT}/ (gfa + colours.csv + legend.tsv per family + index.md)")


if __name__ == "__main__":
    main()
