#!/usr/bin/env python3
"""Mark files that are LIKELY dead or LIKELY superseded, for a later human-approved cleanup.

READ-ONLY: this moves, edits and deletes nothing. It writes two files:
  docs/cleanup_candidates.tsv   one row per file (class, confidence, evidence)
  docs/CLEANUP_CANDIDATES.md    method + summary, regenerated on every run

A file is ANCHORED (kept) if a path-like token naming it appears in an anchor source:
the ledger, the negative-results register, PREREGs, any top-level docs/*.md or docs/experiments
write-up, the auto-memory directory, Rust sources/tests, or Cargo.toml -- or, transitively, in a
script that is itself anchored (and a script that names anchored data is anchored too, so the
generator of a cited number is kept). Everything else is classified by explicit rules
(see RULES below); the rules are the whole method, nothing is judged by eye.

Usage:  python3 tools/audit_cleanup_candidates.py [--repo .] [--memory-dir DIR] [--today YYYY-MM-DD]
"""
import argparse
import collections
import datetime as dt
import fnmatch
import hashlib
import os
import re
import subprocess
import sys

EXCLUDE_PREFIXES = (".git/", "target/", ".worktrees/", "tools/stringtie/", ".remember/",
                    ".superpowers/", ".claude/")
CACHE_DIR_NAMES = ("__pycache__", ".pytest_cache")
TEXT_EXT = {"md", "txt", "py", "sh", "bash", "js", "mjs", "R", "r", "rs", "toml", "awk", "Rmd",
            "yaml", "yml", "cfg", "ini", "mk", "Makefile", "nf", "smk", "pl"}
SCRIPT_EXT = {"py", "sh", "bash", "js", "mjs", "R", "r", "awk", "pl", "nf", "smk", "Rmd"}
DATA_EXT = {"tsv", "csv", "json", "log", "txt", "gtf", "gff", "gff3", "fa", "fasta", "fai", "bed",
            "paf", "sam", "bam", "bai", "gfa", "png", "pdf", "svg", "html", "npz", "rds",
            "gffcmp", "tmap", "refmap", "stats", "loci", "tracking", "vcf", "npy", "pkl", "out"}
SIDECAR_EXT = {"fai", "bai", "csi", "tbi", "gzi", "crai"}
LEGACY_KEYWORDS = ("stringtie", "bundle", "gffcompare", "transfrag", "keeptrf", "rlink", "bnode",
                   "parity", "needy", "snapshot", "trace_bootstrap", "max_flow", "flow_iter",
                   "splice_graph", "subbundle")
SCOPE_DROPPED = re.compile(r"(^|[/_])(asj|allele_specific_junction)", re.I)
PORT_RE = re.compile(r"\b(port of|ported|faithful (rust )?port|native rust port|mirrors|migrat\w*"
                     r"|replaces|superseded by|supersedes)\b", re.I)
STALE_RE = re.compile(r"\b(stale|retracted|obsolete|deprecated|superseded|dead[- ]end|dead code|deleted|removed)\b", re.I)
WRITE_RE = re.compile(r"(open\([^)]*['\"][wa]b?['\"]|to_csv|to_json|json\.dump|savefig|ggsave|saveRDS|write\w*\(|"
                      r"\bout(put)?\w*\s*=|OUT\w*\s*=|>\s*\S+\s*$|\s-o\s|--out)", re.I)
PARTIAL_PORT_RE = re.compile(r"^[`'\"]?(::|'s\b|\s+(loaders?|helpers?|functions?|parser|core|leg|step)\b)", re.I)
IMPORT_RE = re.compile(r"^\s*(?:from\s+([A-Za-z_][\w.]*)\s+import\b|import\s+([A-Za-z_][\w.]*(?:\s+as\s+\w+)?"
                       r"(?:\s*,\s*[A-Za-z_][\w.]*(?:\s+as\s+\w+)?)*))")
TOKEN_RE = re.compile(r"[A-Za-z0-9_.+\-/*~${}]+")
EXT_OK = re.compile(r"\.[A-Za-z][A-Za-z0-9]{0,6}$")
HEADING_RE = re.compile(r"^#{1,4} (§[0-9]+[a-z0-9]*)?")
LEDGER_SECTION_RE = re.compile(r"§[0-9]+[a-z0-9]*")

RULES = [
    ("PROTECTED", "-", "build/config files, Rust sources and tests (Rust reachability is "
     "docs/MODULE_STATUS.md's job, enforced by module_status_tests), and the anchor docs "
     "themselves. Never a candidate."),
    ("TEMP", "high", "Python/pytest caches; untracked or git-ignored files at the repo root; "
     "untracked *.log / *err* / *out* / *.patch.txt / checkpoint files that nothing cites."),
    ("REFUTED-MODULE", "medium", "Rust module whose `//! **STATUS:**` header is REFUTED and that no other "
     "file names (a REFUTED module that is still imported, e.g. collapse_gate.rs, stays PROTECTED)."),
    ("SUPERSEDED-PORTED", "medium", "Python script that a Rust source line declares it ports "
     "('Port of', 'Faithful Rust port of', 'Mirrors', 'migration'); 'low' when only a function or "
     "part is ported (`x.py::f`, `x.py loaders`). The Python may still serve as a parity oracle "
     "or golden-fixture generator -- check tests before deleting."),
    ("SUPERSEDED-VERSION", "medium", "Older member of a version series in the same directory "
     "(_v1.._vN, foo/foo2/foo3, dated _YYYY-MM-DD copies, foo vs foo_fix/_final/_new) that "
     "nothing anchors. 'high' when the newest member IS anchored."),
    ("SUPERSEDED-CITED", "low", "Older member of a version series that IS anchored: provenance "
     "for a recorded result -- archive, don't delete."),
    ("LEGACY-ASSEMBLER", "medium", "Not anchored, and its path or first 200 lines name "
     "StringTie-era assembler machinery (bundle/transfrag/parity/gffcompare/...) -- the "
     "assembler layer was retired (docs/RETIREMENT_AND_MIGRATION.md)."),
    ("AMBIGUOUS-CITE", "low", "Not anchored; an anchor source names it only by a bare basename "
     "shared by several files, a directory too large (>40 files) or a glob matching several files -- may or may not "
     "mean this copy."),
    ("ORPHAN", "medium", "Not anchored; named only by files that are themselves not anchored "
     "(e.g. a figure named only by its un-cited plotting script)."),
    ("UNCITED-STALE", "medium", "Not named by anything (wide globs / big directories / scripts' "
     "ambiguous basenames ignored), last touched more than 30 days ago."),
    ("PROBABLE-PROVENANCE", "-", "Not cited by name, but it sits in an experiment directory (below "
     "bench/, docs/, ...) holding anchored files, or under a folder whose anchored README/write-up covers it -- usually an output of that experiment written under "
     "a computed name and cited as a folder. Verification judged 5/6 such files provenance: not a candidate."),
    ("UNCITED-RECENT", "low", "Not cited, touched within 30 days: may be work in progress."),
    ("KEEP-CITED", "-", "Anchored directly or transitively. Not a candidate."),
]
CANDIDATE_CLASSES = [r[0] for r in RULES if r[0] not in ("PROTECTED", "KEEP-CITED", "PROBABLE-PROVENANCE")]


def git(repo, *args):
    return subprocess.run(["git", "-C", repo, *args], capture_output=True, check=True).stdout


def list_files(repo):
    status = {}
    for p in git(repo, "ls-files", "-z").decode().split("\0"):
        if p:
            status[p] = "tracked"
    for p in git(repo, "ls-files", "-z", "--others", "--exclude-standard").decode().split("\0"):
        if p:
            status.setdefault(p, "untracked")
    for p in git(repo, "ls-files", "-z", "--others", "--ignored", "--exclude-standard").decode().split("\0"):
        if p:
            status.setdefault(p, "ignored")
    out, caches = {}, collections.Counter()
    for p, s in status.items():
        if p.startswith(EXCLUDE_PREFIXES) or p == "tools/stringtie":
            continue
        parts = p.split("/")
        cache = next((i for i, c in enumerate(parts) if c in CACHE_DIR_NAMES), None)
        if cache is not None:
            caches["/".join(parts[:cache + 1]) + "/"] += 1
            continue
        if p.endswith((".pyc", ".pyo")):
            continue
        if os.path.isfile(os.path.join(repo, p)):
            out[p] = s
    return out, caches


def git_last_dates(repo):
    last = {}
    date = None
    for line in git(repo, "log", "--no-renames", "--format=>%cs", "--name-only").decode(errors="replace").splitlines():
        if line.startswith(">"):
            date = line[1:]
        elif line and line not in last:
            last[line] = date
    return last


def ext_of(path):
    base = path.rsplit("/", 1)[-1]
    return base.rsplit(".", 1)[-1] if "." in base else base


def citer_tier(path, memory=False):
    """anchor = counts toward keeping; weak = evidence only."""
    if memory:
        return "memory"
    base = path.rsplit("/", 1)[-1]
    if path == "docs/o1_ledger.md":
        return "ledger"
    if path == "docs/NEGATIVE_RESULTS_REGISTER.md":
        return "register"
    if base.startswith("PREREG_"):
        return "prereg"
    if (path.startswith("src/") and path.endswith(".rs")) or path.startswith("tests/") or path == "Cargo.toml":
        return "code"
    if path == "README.md" or (path.startswith("docs/") and path.count("/") == 1 and path.endswith(".md")):
        return "doc"
    if path.startswith("docs/experiments/") and path.endswith(".md"):
        return "doc"
    return "weak"


ANCHOR_TIERS = {"ledger", "register", "prereg", "code", "doc", "memory"}


class Resolver:
    """Maps a path-like token found in text to the repo files it names."""

    def __init__(self, files, series_members=()):
        self.files = set(files)
        stems = collections.defaultdict(set)
        for f in series_members:
            stem = re.sub(r"[_-]?\d{4}-\d{2}-\d{2}$", "", f.rsplit("/", 1)[-1].split(".", 1)[0])
            if len(stem) >= 6 and re.search(r"\d$", stem):
                stems[stem].add(f)
        self.stems = {k: next(iter(v)) for k, v in stems.items() if len(v) == 1}
        self.suffix = collections.defaultdict(set)
        self.dotprefix = collections.defaultdict(set)
        self.dirs = collections.defaultdict(set)
        for f in files:
            parts = f.split("/")
            for i in range(len(parts)):
                self.suffix["/".join(parts[i:])].add(f)
            base = parts[-1]
            for m in re.finditer(r"\.", base):
                if m.start() > 0:
                    self.dotprefix["/".join(parts[:-1] + [base[:m.start()]])].add(f)
            for i in range(2, len(parts)):
                d = "/".join(parts[:i])
                for j in range(i - 1):  # every suffix of the dir with >= 2 components
                    self.dirs["/".join(parts[j:i])].add(f)
        self.memo = {}

    @staticmethod
    def normalize(tok):
        tok = re.sub(r"^\*\*|\*\*$", "", tok.split("::", 1)[0])
        if "/Rustle/" in tok:
            tok = tok.rsplit("/Rustle/", 1)[1]
        tok = re.sub(r"^\$\{?\w+\}?/", "", tok)
        tok = re.sub(r"^(\.\.?/)+", "", tok)
        tok = re.sub(r":\d+(-\d+)?$", "", tok)
        return tok.strip(".,;:)(-")

    def resolve_import(self, citer, module):
        """A Python import -> the repo .py file(s) it loads (citer's dir and its ancestors, as with sys.path hacks)."""
        rel = module.replace(".", "/")
        d = citer.rpartition("/")[0]
        while True:
            for cand in (f"{d}/{rel}.py" if d else f"{rel}.py", f"{d}/{rel}/__init__.py" if d else f"{rel}/__init__.py"):
                if cand in self.files:
                    return [cand]
            if not d:
                break
            d = d.rpartition("/")[0]
        hits = self.suffix.get(rel + ".py", ())
        return list(hits) if len(hits) == 1 else []

    def expand(self, tok):
        m = re.search(r"\{([^{}]*,[^{}]*)\}", tok)
        if not m:
            return [tok]
        return [tok[:m.start()] + alt + tok[m.end():] for alt in m.group(1).split(",")][:20]

    def resolve(self, raw):
        if raw not in self.memo:
            self.memo[raw] = list(self._resolve(raw))
        return self.memo[raw]

    def _resolve(self, raw):
        """Yield (file, kind), kind in exact|suffix|basename-ambiguous|dir|dir-large|prefix|glob."""
        tok = self.normalize(raw)
        if len(tok) < 4 or len(tok) > 300 or tok.startswith(("http", "/")):
            return
        if tok in self.stems:
            yield self.stems[tok], "stem"
            return
        if "/" not in tok and not EXT_OK.search(tok):
            return
        for t in self.expand(tok):
            if "*" in t:
                lit = re.sub(r"[*?]", "", t.rsplit("/", 1)[-1])
                if ("/" not in t and len(lit) < 6) or t.startswith("*"):
                    continue
                n = t.count("/") + 1
                hits = [f for f in self.files
                        if f.count("/") + 1 >= n and fnmatch.fnmatchcase("/".join(f.split("/")[-n:]), t)]
                kind = "glob" if len(hits) == 1 else "glob-wide"
                for f in hits:
                    yield f, kind
                continue
            t = t.rstrip("/")
            hits = self.suffix.get(t)
            if hits:
                if "/" in t or len(hits) == 1:
                    kind = "exact" if t in self.files else "suffix"
                    for f in hits:
                        yield f, kind
                else:
                    for f in hits:
                        yield f, "basename-ambiguous"
                continue
            if "/" in t:
                if t in self.dotprefix:
                    for f in self.dotprefix[t]:
                        yield f, "prefix"
                elif t in self.dirs and len(self.dirs[t]) <= 600:
                    kind = "dir" if len(self.dirs[t]) <= 40 else "dir-large"
                    for f in self.dirs[t]:
                        yield f, kind


def read_text(path):
    try:
        with open(path, "rb") as fh:
            data = fh.read(8_000_000)
    except OSError:
        return ""
    if b"\0" in data[:4096]:
        return ""
    return data.decode("utf-8", errors="replace")


def version_series(files):
    """Return {older_file: (newest_file, series_kind)}."""
    groups = collections.defaultdict(list)
    by_dir = collections.defaultdict(set)
    for f in files:
        d, _, base = f.rpartition("/")
        by_dir[d].add(base)
    for f in files:
        d, _, base = f.rpartition("/")
        stem, dot, ext = base.partition(".")
        m = re.match(r"^(.*?)[_-]?(\d{4}-\d{2}-\d{2})(.*)$", stem)
        if m:
            groups[(d, "date", m.group(1), m.group(3), ext)].append((m.group(2), f))
        m = re.match(r"^(.*?)[_-]v(\d+)(.*)$", stem)
        if m:
            groups[(d, "v", m.group(1), m.group(3), ext)].append((int(m.group(2)), f))
            plain = m.group(1) + m.group(3) + dot + ext
            if plain in by_dir[d]:
                groups[(d, "v", m.group(1), m.group(3), ext)].append((0, (d + "/" if d else "") + plain))
        m = re.match(r"^(.*[A-Za-z_])(\d)$", stem)
        if m and len(m.group(1)) >= 3 and (m.group(1) + dot + ext) in by_dir[d] and not re.search(r"(chr|fam|mcl|top|next|STRG_|_[a-z])$", m.group(1)):
            groups[(d, "n", m.group(1), "", ext)].append((int(m.group(2)), f))
            groups[(d, "n", m.group(1), "", ext)].append((1, (d + "/" if d else "") + m.group(1) + dot + ext))
        m = re.match(r"^(.*)_(old|bak|backup|orig)$", stem)
        if m and (m.group(1) + dot + ext) in by_dir[d]:
            groups[(d, "old", m.group(1), "", ext)] += [(0, f), (1, (d + "/" if d else "") + m.group(1) + dot + ext)]
        m = re.match(r"^(.*)_(fix|fixed|final|new|corrected)$", stem)
        if m and (m.group(1) + dot + ext) in by_dir[d]:
            groups[(d, "fix", m.group(1), "", ext)] += [(1, f), (0, (d + "/" if d else "") + m.group(1) + dot + ext)]
    older = {}
    for key, members in groups.items():
        members = sorted(set(members))
        if len({f for _, f in members}) < 2:
            continue
        if key[1] == "v" and max(n for n, _ in members) < 2:
            continue  # a lone `_v1` beside a base is a label (e.g. 'from_v1_sedef'), not a series
        newest = members[-1][1]
        for _, f in members[:-1]:
            if f != newest and f not in older:
                older[f] = (newest, key[1])
    return older


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--memory-dir", default=os.path.expanduser(
        "~/.claude/projects/-mnt-c-Users-jfris-Desktop/memory"))
    ap.add_argument("--today", default=dt.date.today().isoformat())
    ap.add_argument("--out-tsv", default="docs/cleanup_candidates.tsv")
    ap.add_argument("--out-md", default="docs/CLEANUP_CANDIDATES.md")
    a = ap.parse_args()
    repo = a.repo
    today = dt.date.fromisoformat(a.today)
    stale_cut = (today - dt.timedelta(days=30)).isoformat()
    head = git(repo, "rev-parse", "--short", "HEAD").decode().strip()
    branch = git(repo, "rev-parse", "--abbrev-ref", "HEAD").decode().strip()

    status, caches = list_files(repo)
    own_outputs = {a.out_tsv, a.out_md}
    for p in own_outputs:
        status.pop(p, None)
    files = sorted(status)
    last = git_last_dates(repo)
    info = {}
    for f in files:
        full = os.path.join(repo, f)
        st = os.stat(full)
        touched = last.get(f) if status[f] == "tracked" else None
        touched = touched or dt.date.fromtimestamp(st.st_mtime).isoformat()
        info[f] = {"size": st.st_size, "touched": touched}

    series = version_series(files)
    digests = collections.defaultdict(list)
    for f in files:
        if 0 < info[f]["size"] <= 20_000_000:
            with open(os.path.join(repo, f), "rb") as fh:
                digests[hashlib.md5(fh.read()).hexdigest()].append(f)
    identical = {f: [g for g in grp if g != f] for grp in digests.values() if len(grp) > 1 for f in grp}
    res = Resolver(files, {f for f in series} | {n for n, _ in series.values()})
    # ---- citation scan -------------------------------------------------------------------
    cites = collections.defaultdict(lambda: collections.defaultdict(set))  # target -> citer -> kinds
    tiers = {}
    port_of = collections.defaultdict(set)
    stale_mention = collections.defaultdict(set)  # target -> doc/memory citers whose naming line says stale/retracted/...
    writes = collections.defaultdict(set)  # data file -> scripts whose naming line looks like a write
    sections = collections.defaultdict(set)
    dead_refs = {}
    head_text = {}
    citers = [(f, citer_tier(f)) for f in files if ext_of(f) in TEXT_EXT or f.rsplit("/", 1)[-1] in TEXT_EXT]
    if os.path.isdir(a.memory_dir):
        for name in sorted(os.listdir(a.memory_dir)):
            if name.endswith(".md"):
                citers.append((os.path.join(a.memory_dir, name), "memory"))
    for citer, tier in citers:
        is_mem = tier == "memory"
        text = read_text(citer if is_mem else os.path.join(repo, citer))
        if not text:
            continue
        key = "memory:" + os.path.basename(citer) if is_mem else citer
        tiers[key] = tier
        if not is_mem:
            head_text[citer] = "\n".join(text.splitlines()[:200]).lower()
        section = None
        n_refs = n_dead = 0
        for line in text.splitlines():
            if citer.endswith("o1_ledger.md") and line.startswith("#"):
                m = LEDGER_SECTION_RE.search(line)
                section = m.group(0) if m else line.strip("# ")[:25]
            if citer.endswith(".py"):
                m = IMPORT_RE.match(line)
                if m:
                    names = [m.group(1)] if m.group(1) else [x.split()[0] for x in m.group(2).split(",")]
                    for name in names:
                        for tgt in res.resolve_import(citer, name):
                            if tgt != citer:
                                cites[tgt][key].add("import")
            if "/" not in line and "." not in line:
                continue
            porting = bool(PORT_RE.search(line)) and citer.endswith(".rs")
            for raw in TOKEN_RE.findall(line):
                hit = False
                for tgt, kind in res.resolve(raw):
                    if tgt == citer:
                        continue
                    hit = True
                    cites[tgt][key].add(kind)
                    if tier in ("doc", "memory") and STALE_RE.search(line):
                        stale_mention[tgt].add(key)
                    if kind in ("exact", "suffix", "prefix") and ext_of(citer) in SCRIPT_EXT and WRITE_RE.search(line):
                        writes[tgt].add(citer)
                    if section and tier == "ledger":
                        sections[tgt].add(section)
                    if porting and tgt.endswith(".py") and kind in ("exact", "suffix"):
                        after = line[line.find(raw) + len(raw):]
                        partial = "::" in raw or PARTIAL_PORT_RE.match(after)
                        port_of[tgt].add((citer, "partial" if partial else "whole"))
                if not is_mem:
                    t = Resolver.normalize(raw)
                    if re.match(r"^(bench|src|docs|scripts|tools|tests|analysis|figures)/", t) and EXT_OK.search(t) and "*" not in t and "{" not in t:
                        n_refs += 1
                        n_dead += 0 if hit else 1
        if not is_mem:
            dead_refs[citer] = (n_dead, n_refs)

    # ---- anchoring fixpoint ----------------------------------------------------------------
    refuted = set()
    for f in files:
        if f.endswith(".rs") and re.search(r"^//! \*\*STATUS:\*\* REFUTED", head_text.get(f, ""), re.M | re.I):
            refuted.add(f)
    for f in refuted:
        tiers[f] = "weak"
        mod = f.rsplit("/", 1)[-1][:-3]
        for g in files:
            if g != f and g.endswith(".rs") and re.search(rf"\b{mod}::", read_text(os.path.join(repo, g))):
                cites[f][g].add("rust-use")

    def strong(citer, kinds):
        return bool(kinds - {"basename-ambiguous", "dir-large", "glob-wide"})

    anchored = {}
    for tgt, by in cites.items():
        why = sorted(c for c, k in by.items() if tiers.get(c) in ANCHOR_TIERS and strong(c, k))
        if why:
            anchored[tgt] = why
    changed = True
    while changed:
        changed = False
        for tgt, by in cites.items():
            if tgt in anchored:
                continue
            via = sorted(c for c, k in by.items()
                         if c in anchored and (ext_of(c) in SCRIPT_EXT or c.endswith(".md")) and strong(c, k))
            if via:
                anchored[tgt] = ["via " + v for v in via]
                changed = True
        # index sidecars (.fai/.bai/...) inherit their parent's anchor
        for f in files:
            parent, _, sfx = f.rpartition(".")
            if f not in anchored and sfx in SIDECAR_EXT and parent in anchored:
                anchored[f] = ["sidecar of " + parent]
                changed = True
        # the generator of anchored data is anchored too (a mere reader is not)
        for tgt, gens in writes.items():
            if tgt not in anchored or ext_of(tgt) not in DATA_EXT:
                continue
            for c in gens:
                if c not in anchored:
                    anchored[c] = ["writes anchored data " + tgt]
                    changed = True

    anchored_dirs = {f.rpartition("/")[0] for f in anchored if f.count("/") >= 2}  # not bench/ or docs/ roots
    # an experiment folder whose own write-up (e.g. bench/x/README.md) is anchored covers its subfolders too
    documented_dirs = {f.rpartition("/")[0] for f in anchored if f.count("/") >= 2 and f.endswith(".md")}

    def in_documented(f):
        d = f.rpartition("/")[0]
        while d.count("/") >= 1:
            if d in documented_dirs:
                return True
            d = d.rpartition("/")[0]
        return False

    # ---- classify -------------------------------------------------------------------------
    rows = []
    klass = {}
    for f in files:
        ext = ext_of(f)
        base = f.rsplit("/", 1)[-1]
        by = cites.get(f, {})
        strong_by = {c for c, k in by.items() if strong(c, k)}
        weak_anchor_by = {c for c, k in by.items() if not strong(c, k) and tiers.get(c) in ANCHOR_TIERS}
        anchor = anchored.get(f, [])
        weak = sorted(c for c in by if c not in anchor and tiers.get(c) not in ANCHOR_TIERS)
        ev = []
        superseded_by = ""
        touched = info[f]["touched"]
        if f in refuted and not by:
            c, conf = "REFUTED-MODULE", "medium"
        elif f in refuted:
            c, conf = "PROTECTED", "-"
            ev.append("status-REFUTED-but-still-referenced-by:" + sorted(by)[0])
        elif (f.startswith(("src/", "tests/")) and ext == "rs") or f in (
                "Cargo.toml", "Cargo.lock", ".gitignore", ".gitmodules", "LICENSE", "README.md") \
                or f.startswith(".cargo/") or f == "tools/audit_cleanup_candidates.py" \
                or tiers.get(f) in ("ledger", "register", "prereg") or f == "docs/MEMORY_DIGEST.md":
            c, conf = "PROTECTED", "-"
        elif f in port_of:
            whole = any(k == "whole" for _, k in port_of[f])
            c, conf = "SUPERSEDED-PORTED", "medium" if whole else "low"
            superseded_by = ",".join(sorted({r for r, _ in port_of[f]}))
            if not whole:
                ev.append("partial-port")
            users = sorted(c for c, k in by.items() if c in anchored and ext_of(c) in SCRIPT_EXT and strong(c, k))
            if users:  # e.g. a gen_*_fixture.py still imports it to build the Rust golden fixture
                conf = "low"
                ev.append("still-used-by:" + users[0])
        elif f in series:
            newest, kind = series[f]
            superseded_by = newest
            ev.append(f"series:{kind}")
            if anchor:
                c, conf = "SUPERSEDED-CITED", "low"
            else:
                c, conf = "SUPERSEDED-VERSION", "high" if newest in anchored else "medium"
        elif anchor:
            c, conf = "KEEP-CITED", "-"
        elif status[f] != "tracked" and not strong_by and not weak_anchor_by and (
                "/" not in f or ext in ("log", "npz") or re.search(r"(^|[_.-])(err|out|tmp|scratch)([_.-]|$)", base, re.I)
                or base.endswith(".patch.txt")):
            c, conf = "TEMP", "high"
        elif any(k in f.lower() or (ext in SCRIPT_EXT and k in head_text.get(f, "")) for k in LEGACY_KEYWORDS) \
                and sum(k in (f.lower() + head_text.get(f, "")) for k in LEGACY_KEYWORDS) >= 2:
            c, conf = "LEGACY-ASSEMBLER", "medium" if touched < stale_cut else "low"
        elif weak_anchor_by:
            c, conf = "AMBIGUOUS-CITE", "low"
        elif strong_by:
            c, conf = "ORPHAN", "medium" if touched < stale_cut else "low"
            ev.append("only-cited-by-candidates")
        elif (f.rpartition("/")[0] in anchored_dirs or in_documented(f)) and ext in DATA_EXT:
            c, conf = "PROBABLE-PROVENANCE", "-"
        elif touched < stale_cut:
            c, conf = "UNCITED-STALE", "medium"
        else:
            c, conf = "UNCITED-RECENT", "low"
        if f in stale_mention and anchor:
            ev.append("named-beside-stale-language-by:" + sorted(stale_mention[f])[0])
        if f in identical:
            ev.append("byte-identical-to:" + identical[f][0] + (f"(+{len(identical[f]) - 1})" if len(identical[f]) > 1 else ""))
        if SCOPE_DROPPED.search(f):
            ev.append("asj-objective-dropped")
        if f in dead_refs and dead_refs[f][1] >= 5 and dead_refs[f][0] / dead_refs[f][1] > 0.5:
            ev.append(f"dead-path-refs:{dead_refs[f][0]}/{dead_refs[f][1]}")
        kinds = sorted({k for ks in by.values() for k in ks})
        if kinds and kinds != ["exact"]:
            ev.append("match:" + "+".join(kinds))
        klass[f] = c
        rows.append([f, status[f], info[f]["size"], touched, c, conf, superseded_by,
                     len(anchor), ";".join(anchor[:4]), ",".join(sorted(sections.get(f, []))[-6:]),
                     ";".join(weak[:4]), ";".join(ev)])
    for cache, n in sorted(caches.items()):
        rows.append([cache, "ignored", 0, "-", "TEMP", "high", "", 0, "", "", "", f"cache-dir:{n}-files"])

    # ---- write ----------------------------------------------------------------------------
    header = ["path", "git_status", "bytes", "last_touched", "class", "confidence", "superseded_by",
              "n_anchor_citers", "anchor_citers", "ledger_sections", "weak_citers", "evidence"]
    with open(os.path.join(repo, a.out_tsv), "w") as fh:
        fh.write(f"# generated {a.today} by tools/audit_cleanup_candidates.py at {branch}@{head}; "
                 f"read-only audit, nothing was moved\n")
        fh.write("\t".join(header) + "\n")
        order = {c: i for i, (c, _, _) in enumerate(RULES)}
        conf_order = {"high": 0, "medium": 1, "low": 2, "-": 3}
        for r in sorted(rows, key=lambda r: (order[r[4]], conf_order[r[5]], r[0])):
            fh.write("\t".join(str(x).replace("\t", " ") for x in r) + "\n")
    write_md(os.path.join(repo, a.out_md), rows, a, branch, head, stale_cut)
    n_cand = sum(r[4] in CANDIDATE_CLASSES for r in rows)
    print(f"{len(rows)} rows, {n_cand} candidates -> {a.out_tsv}, {a.out_md}", file=sys.stderr)


def human(n):
    for unit in ("B", "K", "M", "G"):
        if n < 1024:
            return f"{n:.0f}{unit}"
        n /= 1024
    return f"{n:.1f}T"


def write_md(path, rows, a, branch, head, stale_cut):
    counts = collections.Counter((r[4], r[5]) for r in rows)
    by_class = collections.defaultdict(list)
    for r in rows:
        by_class[r[4]].append(r)
    L = []
    L.append("# Cleanup candidates — likely dead / likely superseded files\n")
    L.append(f"Generated {a.today} at `{branch}@{head}` by `tools/audit_cleanup_candidates.py` "
             "(re-run it; this file is overwritten). **Read-only: nothing was moved, edited or deleted.** "
             f"Full per-file table: `{a.out_tsv}` (filter on `class` and `confidence`).\n")
    L.append("⚠ A mark is a *candidate*, not a verdict. Before deleting anything: (1) grep the path once more, "
             "(2) check the `ledger_sections` / `anchor_citers` columns, (3) prefer `git mv` into an "
             "archive directory over `rm` for anything tracked, (4) remember `bench/` data can be slow "
             "to regenerate (AGENTS.md §2).\n")
    if os.path.exists(os.path.join(os.path.dirname(path), "CLEANUP_CANDIDATES_VERIFICATION.md")):
        L.append("Measured precision of these marks (blind, pre-registered spot-check): "
                 "`docs/CLEANUP_CANDIDATES_VERIFICATION.md`.\n")
    L.append("## How a file is kept\n")
    L.append("A file is **anchored** if a path-like token naming it (full path, unique path suffix, "
             "unique basename, output prefix, a glob matching exactly 1 file, an enclosing directory of ≤40 files, "
             "or a Python `import` resolved against the importing script's directory and its parents) appears in "
             "an anchor source: `docs/o1_ledger.md`, `docs/NEGATIVE_RESULTS_REGISTER.md`, "
             "`docs/PREREG_*.md`, top-level `docs/*.md`, `docs/experiments/*.md`, `README.md`, the "
             "auto-memory directory, Rust `src/`/`tests/`, or `Cargo.toml`. Anchoring then propagates: "
             "a file named by an anchored script or anchored markdown write-up is anchored, index sidecars "
             "(`.fai`, `.bai`) follow their parent, and a script that names anchored data "
             "(its generator — the naming line looks like a write) is anchored; a script that merely reads anchored data is not. An ambiguous bare basename (e.g. `reads.fa`, "
             "found in many directories) never anchors on its own. Citations from `docs/archive/`, "
             "`docs/superpowers/`, `AGENTS.md` and un-anchored scripts are *weak*: recorded, not counted.\n")
    L.append(f"## Rules (applied in this order; stale cutoff = {stale_cut})\n")
    L.append("| class | confidence | rule | files | bytes |")
    L.append("|---|---|---|---:|---:|")
    for c, conf, rule in RULES:
        rs = by_class.get(c, [])
        L.append(f"| **{c}** | {conf} | {rule} | {len(rs)} | {human(sum(int(r[2]) for r in rs))} |")
    L.append("")
    cand = [r for r in rows if r[4] in CANDIDATE_CLASSES]
    L.append(f"**{len(cand)} candidates** of {len(rows)} files "
             f"({sum(1 for r in cand if r[5] == 'high')} high, {sum(1 for r in cand if r[5] == 'medium')} medium, "
             f"{sum(1 for r in cand if r[5] == 'low')} low confidence).\n")
    L.append("## Candidates by directory\n")
    L.append("| directory | " + " | ".join(c for c in CANDIDATE_CLASSES) + " | kept |")
    L.append("|---|" + "---:|" * (len(CANDIDATE_CLASSES) + 1))
    dirs = collections.defaultdict(collections.Counter)
    for r in rows:
        parts = r[0].split("/")
        d = "/".join(parts[:2]) + "/" if len(parts) > 2 else (parts[0] + "/" if len(parts) > 1 else "(root)")
        dirs[d][r[4]] += 1
    for d in sorted(dirs, key=lambda d: -sum(dirs[d][c] for c in CANDIDATE_CLASSES)):
        cnt = dirs[d]
        if not any(cnt[c] for c in CANDIDATE_CLASSES):
            continue
        L.append(f"| `{d}` | " + " | ".join(str(cnt[c] or "") for c in CANDIDATE_CLASSES)
                 + f" | {cnt['KEEP-CITED'] + cnt['PROTECTED'] + cnt['PROBABLE-PROVENANCE']} |")
    L.append("")
    for c in ("TEMP", "REFUTED-MODULE", "SUPERSEDED-PORTED", "SUPERSEDED-VERSION"):
        rs = sorted(by_class.get(c, []), key=lambda r: r[0])
        if not rs:
            continue
        L.append(f"## {c} ({len(rs)})\n")
        for r in rs[:200]:
            extra = f" → superseded by `{r[6]}`" if r[6] else ""
            ev = f" ({r[11]})" if r[11] else ""
            L.append(f"- `{r[0]}` [{r[1]}, {r[5]}]{extra}{ev}")
        if len(rs) > 200:
            L.append(f"- … {len(rs) - 200} more in the TSV")
        L.append("")
    flagged = [r for r in rows if "asj-objective-dropped" in r[11]]
    if flagged:
        L.append(f"## Overlay: files tied to the dropped ASJ objective ({len(flagged)})\n")
        L.append("Not classified as dead by this flag alone — ASJ was dropped as an objective (memory, "
                 "2026-08-07), but the binaries still build. Decide as a scope question, not a cleanup one.\n")
        for r in sorted(flagged):
            L.append(f"- `{r[0]}` — {r[4]}")
        L.append("")
    L.append("## Known blind spots\n")
    L.append("- Paths built at run time (`f\"bench/{name}.tsv\"`, shell loops) are invisible to a token scan: "
             "such outputs land in UNCITED-IN-ANCHORED-DIR or UNCITED-*, never KEEP.")
    L.append("- A citation proves a file was *named*, not that the naming text is still true; an anchor in a "
             "retracted ledger section still anchors (see `ledger_sections`).")
    L.append("- Data outside the repo (`/mnt/linuxdisk`, BAMs) is not audited; `.worktrees/`, `tools/stringtie` "
             "(submodule), `.remember/`, `.superpowers/` are excluded.")
    L.append("- Rust module reachability is not re-derived here — see `docs/MODULE_STATUS.md`.")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
