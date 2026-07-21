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


def main_flagships():
    """Legacy 3-flagship demo -> scratch OUT (kept for the original ceiling-demo entrypoint)."""
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


# ============================================================================
# ALL-83-FAMILIES orchestration (2026-07-20/21). Each family is rendered in a
# SEPARATE child process (`--build-one`) so a hard crash (abpoa OOM-kill /
# segfault — NOT always a catchable Python exception) only kills that one
# child; the parent orchestrator sees a non-zero/negative return code, logs
# it, and moves on to the next family. Never touches the flagship code above.
# ============================================================================
DNA_GRAPHS_DIR = "bench/soto/dna_graphs"

INDEX_COLS = ["family_id", "gene", "n_members", "n_green_rna_recovered",
              "n_red_dna_only", "n_grey", "gfa_path", "status", "skip_reason"]

SIGNAL_REASONS = {
    6: "SIGABRT (abpoa/abort)", 9: "SIGKILL (likely OOM-killed by the kernel)",
    11: "SIGSEGV (segfault, likely abpoa OOM/crash)", 15: "SIGTERM",
}


def all_family_ids(bed_lines):
    """Unique family IDs, in first-BED-appearance order then numerically re-sorted by the trailing
    integer (ID_8 before ID_35 before ID_131) for a stable, readable INDEX."""
    seen, order = set(), []
    for ln in bed_lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4 or "|" not in f[3]:
            continue
        _, fam = f[3].rsplit("|", 1)
        if fam not in seen:
            seen.add(fam)
            order.append(fam)

    def fam_num(fid):
        try:
            return int(fid.split("_", 1)[1])
        except (IndexError, ValueError):
            return 10 ** 9
    return sorted(order, key=fam_num)


def family_label(family_id, members):
    """Best-effort readable gene-family stem from BED gene names (e.g. SRGAP2C/SRGAP2D/SRGAP2B ->
    'SRGAP2'; GOLGA8CP/GOLGA8EP/... -> 'GOLGA8'). Falls back to family_id when the members' gene names
    share no usable stem (e.g. a family of lone accessions / unrelated symbols) -- 'else the ID'."""
    import re
    from collections import Counter
    genes = [g for g, _c, _s, _e in members]
    if not genes:
        return family_id
    cleaned = [re.sub(r"\.\d+$", "", g) for g in genes]      # strip accession version e.g. AL669831.1
    prefix = os.path.commonprefix(cleaned)
    if len(prefix) >= 3:
        return prefix

    def alpha_prefix(g):
        m = re.match(r"^[A-Za-z]+", g)
        return m.group(0) if m else ""
    stems = [s for s in (alpha_prefix(g) for g in cleaned) if len(s) >= 3]
    if stems:
        stem, count = Counter(stems).most_common(1)[0]
        if count >= max(2, len(genes) // 2):
            return stem
    return family_id


def sanitize_part(s):
    out = "".join(c if c.isalnum() else "_" for c in s).strip("_")
    return out or "x"


def filename_stem(family_id, label):
    return family_id if label == family_id else f"{family_id}_{sanitize_part(label)}"


def _blank_row(family_id, label, n_members, reason):
    return dict(family_id=family_id, gene=label, n_members=n_members, n_green_rna_recovered=0,
                n_red_dna_only=0, n_grey=0, gfa_path="", status="skipped", skip_reason=reason)


def build_and_write_one(family_id, label, bed_lines, fa, detection, causes, outdir):
    """Build + write one family's DNA VG to outdir. NEVER raises: every catchable failure (no BED
    members, no member sequence in soto_members.fa, a single-member family that can't form a
    comparative graph, minimap2/abpoa errors) is caught and returned as a 'skipped' row with a reason
    instead of propagating. (Uncatchable failures -- an OOM SIGKILL of the whole process -- are the
    caller's job, at the subprocess-isolation layer in _parent_run_all.)"""
    try:
        members = parse_family_members(bed_lines, family_id)
    except Exception as ex:
        return _blank_row(family_id, label, 0, f"parse_family_members failed: {type(ex).__name__}: {ex}")
    if not members:
        return _blank_row(family_id, label, 0, "no members in BED for this family_id")
    present = []
    for gene, chrom, start, end in members:
        try:
            member_seq(fa, chrom, start, end)
            present.append((gene, chrom, start, end))
        except KeyError:
            pass
    if not present:
        return _blank_row(family_id, label, len(members), "no member sequence found in soto_members.fa")
    if len(present) < 2:
        return _blank_row(family_id, label, len(members),
                           f"single-member family (n_present={len(present)}) -- no comparative graph to draw")
    try:
        r = build_family(family_id, members, fa, detection, causes)
    except Exception as ex:
        return _blank_row(family_id, label, len(members), f"{type(ex).__name__}: {ex}")
    if r["n_present"] < 2:
        return _blank_row(family_id, label, len(members),
                           f"fewer than 2 members reached the graph (n_present={r['n_present']}, "
                           f"missing={r['missing']})")
    try:
        stem = filename_stem(family_id, label)
        os.makedirs(outdir, exist_ok=True)
        base = f"{outdir}/{stem}"
        open(f"{base}.gfa", "w").write(r["gfa"])
        open(f"{base}.colours.csv", "w").write(r["colours"])
        open(f"{base}.legend.tsv", "w").write(r["legend"])
    except Exception as ex:
        return _blank_row(family_id, label, len(members), f"write failed: {type(ex).__name__}: {ex}")
    det_by = {g: detection.get((c, s, e), (False, ""))[0] for g, c, s, e in present}
    n_green = sum(1 for v in det_by.values() if v)
    n_red = sum(1 for v in det_by.values() if not v)
    n_grey_nodes = sum(1 for ln in r["colours"].splitlines()[1:] if ln.endswith("," + GREY))
    return dict(family_id=family_id, gene=label, n_members=len(members),
                n_green_rna_recovered=n_green, n_red_dna_only=n_red, n_grey=n_grey_nodes,
                gfa_path=f"{base}.gfa", status="rendered",
                skip_reason=(f"missing from graph: {r['missing']}" if r["missing"] else ""),
                windowed=r.get("windowed", False), window_span=r.get("window_span"))


def _child_build_one(family_id, label):
    """`--build-one` entrypoint: runs in its OWN process so a hard abpoa crash only takes this
    process down. Prints exactly one `RESULT_JSON:{...}` line on success or on any caught failure."""
    import json
    bed = open(BED).read().splitlines()
    fa = read_fasta(MEMFA)
    detection = load_detection(open(DETECT).read().splitlines())
    causes = load_causes(open(CAUSES).read().splitlines())
    row = build_and_write_one(family_id, label, bed, fa, detection, causes, DNA_GRAPHS_DIR)
    print("RESULT_JSON:" + json.dumps(row))


def write_index(rows, outdir):
    lines = ["\t".join(INDEX_COLS)]
    for r in rows:
        lines.append("\t".join(str(r.get(k, "")) for k in INDEX_COLS))
    open(f"{outdir}/INDEX.tsv", "w").write("\n".join(lines) + "\n")


def write_readme(rows, outdir):
    rendered = [r for r in rows if r["status"] == "rendered"]
    skipped = [r for r in rows if r["status"] != "rendered"]
    tot_green = sum(r["n_green_rna_recovered"] for r in rendered)
    tot_red = sum(r["n_red_dna_only"] for r in rendered)
    tot_grey = sum(r["n_grey"] for r in rendered)
    tot_members_rendered = sum(r["n_members"] for r in rendered)
    lines = [
        "# Soto DNA variation graphs -- all families",
        "",
        f"{len(rendered)}/{len(rows)} Soto families rendered as a DNA variation graph "
        f"(base-level abpoa MSA of the members' genomic sequences -> GFA); "
        f"{len(skipped)} skipped (logged below and in `INDEX.tsv`, never silently dropped).",
        "",
        "## What the colours mean",
        "",
        "- **green** (`#1e8e3e`) -- a copy (member) that RNA (IsoSeq A119b) actually RECOVERED "
        "(`soto_member_detection.tsv`, `detected=Y`).",
        "- **red** (`#d93025`) -- a copy that RNA did NOT recover (DNA-only / K=0-floor); its per-member "
        "cause (from `soto_floor_decomposition.tsv`) is in each family's `<stem>.legend.tsv`.",
        "- **grey** (`#9aa0a6`) -- a shared/conserved backbone **node** in the graph that BOTH a green "
        "and a red copy pass through (mixed). Grey is a **node**-level colour only -- every member "
        "(copy) itself is always definitively green or red per RNA detection; a node is grey where "
        "green and red copies' sequences coincide.",
        "",
        "The DNA variation graph is the identifiability CEILING: every Soto DNA-catalog copy is a path "
        "in the graph regardless of whether RNA recovered it. The graph REPRESENTS what is given (all "
        "copies, from Soto's independent DNA-read-depth catalog); it does not 'detect' families.",
        "",
        "## How to open a GFA",
        "",
        "Open any `bench/soto/dna_graphs/<family_id>_<gene>.gfa` in "
        "[Bandage](https://rrwick.github.io/Bandage/) (`File > Load graph`), then "
        "`File > Load CSV` the matching `.colours.csv` and colour nodes by the `Colour` column "
        "(Bandage: Graph drawing -> Node colours -> `Custom colours based on...` -> the loaded CSV "
        "column). The `.legend.tsv` alongside each GFA lists, per member: detected Y/N, "
        "recovered_by, cause (for red members), colour.",
        "",
        "Look up any family fast in `INDEX.tsv` (`family_id`, `gene`, `n_members`, "
        "`n_green_rna_recovered`, `n_red_dna_only`, `n_grey`, `gfa_path`, `status`, `skip_reason`).",
        "",
        "## Summary counts",
        "",
        f"- **{len(rendered)}/{len(rows)}** families rendered, **{len(skipped)}/{len(rows)}** skipped.",
        f"- Across rendered families: **{tot_members_rendered}** members total -- "
        f"**{tot_green} green** (RNA-recovered) / **{tot_red} red** (DNA-only) / "
        f"**{tot_grey} grey graph-nodes** (shared/mixed backbone, node-level not member-level).",
        "- These counts are over RENDERED families only -- members of a SKIPPED family (e.g. the 7 "
        "single-member families that can't form a comparative graph) are not pictured anywhere and are "
        "excluded from the green/red totals above; see `INDEX.tsv` `n_members` for every family "
        "including skipped ones, and `bench/soto/soto_member_detection.tsv` for the ground-truth "
        "green/red split across ALL 362 Soto members genome-wide (the 76.2% RNA-recovery figure).",
        "",
        "## Skipped families",
        "",
    ]
    if skipped:
        lines.append("| family_id | gene | n_members | reason |")
        lines.append("|---|---|---|---|")
        for r in skipped:
            lines.append(f"| {r['family_id']} | {r['gene']} | {r['n_members']} | {r['skip_reason']} |")
    else:
        lines.append("(none)")
    lines += [
        "",
        "## Regenerating",
        "",
        "GFAs are regenerable (deterministic given the same inputs): "
        "`/home/juanfra/miniforge3/bin/python3 bench/soto/soto_dna_vg.py` (needs `pyabpoa`; run from the "
        "repo root). Each family renders in its own child process so one abpoa OOM/crash cannot take "
        "down the run.",
        "",
    ]
    open(f"{outdir}/README.md", "w").write("\n".join(lines) + "\n")


def _parent_run_all(outdir=DNA_GRAPHS_DIR, timeout_s=900):
    """Orchestrator: iterate ALL Soto families, render each in its own child process (`--build-one`),
    never let one family's crash/OOM/timeout kill the run. Foreground, serial -- deliberately NOT
    parallel (memory-bound abpoa, WSL2 box)."""
    import subprocess
    os.makedirs(outdir, exist_ok=True)
    bed = open(BED).read().splitlines()
    fam_ids = all_family_ids(bed)
    rows = []
    for i, fid in enumerate(fam_ids, 1):
        members = parse_family_members(bed, fid)
        label = family_label(fid, members)
        print(f"[{i}/{len(fam_ids)}] {fid} ({label}, {len(members)} members) ... ", end="", flush=True)
        try:
            proc = subprocess.run(
                [sys.executable, os.path.abspath(__file__), "--build-one", fid, label],
                capture_output=True, text=True, timeout=timeout_s)
        except subprocess.TimeoutExpired:
            row = _blank_row(fid, label, len(members),
                              f"child process timed out after {timeout_s}s (likely abpoa hang/OOM-thrash)")
            rows.append(row)
            print(f"SKIPPED ({row['skip_reason']})")
            continue
        row = None
        if proc.returncode == 0:
            for line in proc.stdout.splitlines():
                if line.startswith("RESULT_JSON:"):
                    import json
                    row = json.loads(line[len("RESULT_JSON:"):])
                    break
        if row is None:
            rc = proc.returncode
            reason = SIGNAL_REASONS.get(-rc, f"child process exited {rc}") if rc < 0 else f"child process exited {rc}"
            stderr_tail = [ln for ln in (proc.stderr or "").strip().splitlines()[-3:]]
            if stderr_tail:
                reason += " -- " + " | ".join(stderr_tail)
            row = _blank_row(fid, label, len(members), reason)
        rows.append(row)
        print(f"{row['status']}" + (f" ({row['skip_reason']})" if row["skip_reason"] else ""))
    write_index(rows, outdir)
    write_readme(rows, outdir)
    n_rendered = sum(1 for r in rows if r["status"] == "rendered")
    print(f"\n{n_rendered}/{len(rows)} families rendered -> {outdir}/ (INDEX.tsv + README.md)")
    return rows


def main():
    if len(sys.argv) >= 4 and sys.argv[1] == "--build-one":
        _child_build_one(sys.argv[2], sys.argv[3])
    elif len(sys.argv) >= 2 and sys.argv[1] == "--flagships-only":
        main_flagships()
    else:
        _parent_run_all()


if __name__ == "__main__":
    main()
