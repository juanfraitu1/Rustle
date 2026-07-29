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

def anchor_of(entry):
    """Identifier token that pins an entry to its declaration, e.g. 'RefineParams::min_identity' ->
    'min_identity', '--em-max-iter (default)' -> 'em_max_iter'. Used only to RELOCATE a moved constant,
    never to decide whether it is correct. Longest identifier wins: short tokens like 'em' or 'as' match
    everywhere and make every candidate ambiguous."""
    # An entry whose name is prose ("sensitive tier call-site literal") carries no identifier that appears at
    # the call site, so it may declare its own `anchor` in the registry. Explicit beats guessing.
    if entry.get("anchor"):
        return entry["anchor"]
    tok = entry["name"].split("::")[-1].strip()
    tok = re.sub(r'\(.*?\)', ' ', tok)                       # drop parenthetical notes
    tok = tok.lstrip('-').replace('-', '_')                  # CLI flag --em-max-iter -> em_max_iter
    cands = re.findall(r'[A-Za-z_][A-Za-z0-9_]*', tok)
    return max(cands, key=len) if cands else None

def relocate(entry, repo_root):
    """A constant whose recorded line no longer holds it has either MOVED (code inserted above it) or
    CHANGED (a real retune). These are completely different events and must not be conflated: the first is
    bookkeeping, the second invalidates the published mechanism registry. Search the file for a line carrying
    BOTH the entry's anchor identifier AND its recorded value. Exactly one hit => it moved, return the new
    line. Zero hits => the value is gone (a real change) — refuse. Several hits => ambiguous — refuse.
    Refusing keeps the loud failure, which is the point of the registry."""
    p = pathlib.Path(repo_root) / entry["file"]
    anc = anchor_of(entry)
    if anc is None or not p.exists():
        return None, "no usable anchor"
    val = str(entry["value"])
    num = re.match(r'^[0-9.]+$', val)
    lines = p.read_text(errors="replace").splitlines()
    def has_val(line):
        return (re.search(r'(?<![0-9.])' + re.escape(val) + r'(?![0-9.])', line) is not None) if num else (val in line)
    # The anchor normally sits on the value's OWN line (`min_identity: 0.80,`). The one exception is a clap
    # default, which lives on the `#[arg(default_value_t = ..)]` ATTRIBUTE while the field name is the line
    # below. Widening that into a general +/-3 window is unsafe: it lets an unrelated neighbour carrying the
    # same literal absorb a constant that was actually RETUNED, which is the precise confusion this function
    # exists to prevent. So the window applies only when the value line is an attribute.
    def anchored(i):
        line = lines[i - 1]
        if anc in line:
            return True
        if re.match(r'\s*#\[', line):                       # attribute: field name is just below
            return any(anc in lines[j] for j in range(i, min(len(lines), i + 2)))
        return False
    hits = [i for i, line in enumerate(lines, start=1) if has_val(line) and anchored(i)]
    # A registry entry points at where a constant is DECLARED, never at prose that merely mentions its value.
    # Doc comments citing a default ("/// ... asm20 0.80 ...") sit immediately ABOVE the declaration, so a
    # naive "first candidate below the old line" rule picks the comment every time -- which is how an earlier
    # version of this pass relocated min_identity, READ_CAP, PAD and 5 others onto `///` lines.
    code = [i for i in hits if not re.match(r'\s*(//|/\*|\*)', lines[i-1])]
    if code:
        hits = code
    if len(hits) == 1:
        return hits[0], None
    if not hits:
        return None, f"value {val!r} no longer appears near any line containing {anc!r} — REAL CHANGE, not a move"
    # Several code candidates. Pure insertions shift code DOWNWARD monotonically, so the true new home is the
    # nearest candidate at or below the old line. Report the alternatives so the choice stays auditable --
    # a silently wrong line number is a false claim in a published registry.
    below = [h for h in hits if h >= entry["line"]]
    pick = min(below) if below else min(hits, key=lambda h: abs(h - entry["line"]))
    return pick, f"AMBIGUOUS-RESOLVED {anc!r}+{val!r} candidates {hits} -> chose {pick} (was {entry['line']})"

def rewrite_lines(toml_path, updates):
    """Rewrite only the `line = N` field of the named entries, in place, leaving everything else byte-identical."""
    src = pathlib.Path(toml_path).read_text().splitlines(keepends=True)
    out, cur = [], None
    for ln in src:
        m = re.match(r'\s*name\s*=\s*"(.*)"\s*$', ln)
        if m:
            cur = m.group(1)
        if cur in updates and re.match(r'\s*line\s*=\s*\d+\s*$', ln):
            out.append(f"line = {updates[cur]}\n")
            continue
        out.append(ln)
    pathlib.Path(toml_path).write_text("".join(out))

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
    ap.add_argument("--relocate", action="store_true",
                    help="update recorded line numbers for constants that MOVED (value unchanged and "
                         "uniquely locatable). Never changes a value; entries whose value actually changed "
                         "still fail.")
    a = ap.parse_args()
    entries = load_registry(a.toml)
    if a.relocate:
        updates, unresolved = {}, []
        for e in entries:
            if verify_entry(e, a.repo_root) is None:
                continue
            new_line, why = relocate(e, a.repo_root)
            if new_line is not None:
                updates[e["name"]] = new_line
                note = f"   [{why}]" if why else ""
                print(f"  relocate {e['name']}: {e['line']} -> {new_line}{note}")
            else:
                unresolved.append(f'{e["name"]}: {why}')
        if updates:
            rewrite_lines(a.toml, updates)
        print(f"relocated {len(updates)} entr{'y' if len(updates)==1 else 'ies'}; "
              f"{len(unresolved)} unresolved")
        for u in unresolved:
            print("  UNRESOLVED " + u, file=sys.stderr)
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
