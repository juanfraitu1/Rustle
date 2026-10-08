#!/usr/bin/env python3
"""Generate docs/INDEX.md: one line per document under docs/ and docs/archive/*/.

Usage (repo root):  python3 tools/docs_index.py > docs/INDEX.md
Columns: file, the document's own first heading, last commit date, and who cites it
(rows of docs/NEGATIVE_RESULTS_REGISTER.md; memory topic files if the memory dir exists).
"""
import datetime, pathlib, re, subprocess, sys, collections
ROOT = pathlib.Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
MEM = pathlib.Path.home() / ".claude/projects/-mnt-c-Users-jfris-Desktop/memory"

def title(p):
    for ln in p.read_text(encoding="utf-8", errors="ignore").splitlines():
        if ln.startswith("#"):
            return ln.lstrip("# ").strip().replace("|", "/")[:120]
    return ""

def last_commits():
    out = subprocess.run(["git", "log", "--format=%ad|%H", "--date=short", "--name-only", "--", "docs/"],
                         cwd=ROOT, capture_output=True, text=True).stdout
    last, cur = {}, None
    for ln in out.splitlines():
        if re.match(r"\d{4}-\d{2}-\d{2}\|", ln):
            cur = ln.split("|")[0]
        elif ln.startswith("docs/") and cur:
            last.setdefault(pathlib.Path(ln).name, cur)
    return last

def main():
    files = sorted(DOCS.glob("*.md")) + sorted(DOCS.glob("archive/*/*.md"))
    files = [f for f in files if f.name != "INDEX.md"]
    reg = (DOCS / "NEGATIVE_RESULTS_REGISTER.md").read_text(encoding="utf-8", errors="ignore")
    mem = "\n".join(p.read_text(encoding="utf-8", errors="ignore") for p in MEM.glob("*.md")) if MEM.exists() else ""
    last = last_commits()
    def row(f):
        rel = f.relative_to(DOCS).as_posix()
        cites = []
        if f.name != "NEGATIVE_RESULTS_REGISTER.md" and (c := reg.count(f.name)):
            cites.append(f"register ×{c}")
        if (c := mem.count(f.name)):
            cites.append(f"memory ×{c}")
        return f"| [`{f.name}`]({rel}) | {title(f)} | {last.get(f.name, 'uncommitted')} | {', '.join(cites) or '—'} |"
    hdr = "| file | title | last commit | cited by |\n|---|---|---|---|"
    print(f"# docs/ index — every document, one line each (generated {datetime.date.today()})\n")
    print("Regenerate with `python3 tools/docs_index.py > docs/INDEX.md` after adding or archiving a document. "
          "Archived files keep their names; only the directory changed, so grepping a filename still finds it. "
          "Verdicts of archived studies are in the register rows that cite them; pre-registrations are kept verbatim "
          "(the commit that introduced each one is the proof it preceded the result).\n")
    top = [f for f in files if f.parent == DOCS]
    living = [f for f in top if not re.search(r"20\d\d-\d\d-\d\d", f.name)]
    dated = sorted((f for f in top if re.search(r"20\d\d-\d\d-\d\d", f.name)),
                   key=lambda f: re.search(r"20\d\d-\d\d-\d\d", f.name).group(0), reverse=True)
    print(f"## Living documents ({len(living)}, undated, maintained in place)\n\n{hdr}")
    print("\n".join(row(f) for f in living))
    print(f"\n## Open studies and current handoff ({len(dated)})\n\n{hdr}")
    print("\n".join(row(f) for f in dated))
    groups = collections.defaultdict(list)
    for f in files:
        if f.parent != DOCS:
            groups[f.parent.relative_to(DOCS).as_posix()].append(f)
    for g in sorted(groups, reverse=True):
        fs = sorted(groups[g], key=lambda f: (re.search(r"20\d\d-\d\d-\d\d", f.name).group(0) if re.search(r"20\d\d-\d\d-\d\d", f.name) else last.get(f.name, ""), f.name), reverse=True)
        print(f"\n## `{g}/` — {len(fs)} closed documents\n\n{hdr}")
        print("\n".join(row(f) for f in fs))

if __name__ == "__main__":
    main()
