"""Verify each registry heuristic against its source file:line, then emit a TSV.
Anti-drift: if a constant's literal value is no longer on its recorded line, fail loudly."""
import argparse, re, sys, tomllib, pathlib

COLS = ["stage", "tier", "kind", "name", "value", "file", "line", "rationale"]

def load_registry(toml_path):
    with open(toml_path, "rb") as f:
        return tomllib.load(f).get("heuristic", [])

def verify_entry(entry, repo_root):
    p = pathlib.Path(repo_root) / entry["file"]
    if not p.exists():
        return f'{entry["name"]}: file not found: {entry["file"]}'
    lines = p.read_text(errors="replace").splitlines()
    ln = entry["line"]
    if ln < 1 or ln > len(lines):
        return f'{entry["name"]}: line {ln} out of range in {entry["file"]}'
    val = str(entry["value"])
    line = lines[ln - 1]
    if re.match(r'^[0-9.]+$', val):
        found = re.search(r'(?<![0-9.])' + re.escape(val) + r'(?![0-9.])', line) is not None
    else:
        found = val in line
    if not found:
        return (f'{entry["name"]}: value {entry["value"]!r} not on '
                f'{entry["file"]}:{ln} — got: {lines[ln-1].strip()!r}')
    return None

def emit_tsv(entries, out_path):
    with open(out_path, "w") as f:
        f.write("\t".join(COLS) + "\n")
        for e in sorted(entries, key=lambda e: (e["stage"], e["tier"], e["name"])):
            f.write("\t".join(str(e.get(c, "")) for c in COLS) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", required=True)
    ap.add_argument("--toml", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    entries = load_registry(a.toml)
    errors = [msg for e in entries if (msg := verify_entry(e, a.repo_root))]
    if errors:
        print("DRIFT DETECTED — registry does not match source:", file=sys.stderr)
        for m in errors:
            print("  " + m, file=sys.stderr)
        sys.exit(1)
    emit_tsv(entries, a.out)
    print(f"verified {len(entries)} heuristics -> {a.out}")

if __name__ == "__main__":
    main()
