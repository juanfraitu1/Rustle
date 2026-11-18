#!/usr/bin/env python3
"""Family transcripts: the transcripts of an assembly joined across loci by tied molecules (driver stage `famwise`).

    famwise.py GTF FW_PREFIX [--eq-loci EQ.loci.tsv] [--window 10] [--min-link 1]

Reads FW_PREFIX.reads.tsv and FW_PREFIX.placed.tsv (written by `eq_classes.py --famwise FW_PREFIX`: the records of
every tied molecule with both scores, and the untied molecules per member). A node is a (locus, intron chain) the GTF
emits; a record with no intron or with a chain the locus does not emit sits at the locus's remainder (`U`, `other`).
Nodes that one molecule carries at tied records are one family transcript; each molecule is counted once. Per node:
lo = molecules placed only there, hi = lo + every tied molecule that reaches it; a status from the two bands (exact
score ties, and scores within --window of the best; bands on minimap2's ms). Nothing is assigned, added or dropped.
Writes FW_PREFIX.tx.tsv, .classes.tsv, .families.tsv and .gtf (the input GTF with fw_* attributes)."""
import argparse, collections, csv, itertools, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gtf_names import chain, opn

UNSPLICED, OTHER, OUTSIDE, UNSEEN = "U", "other", "outside", "unseen"
FW_WINDOW = 20  # the reads table keeps records this close to the molecule's best, so narrower windows need no BAM pass

def record_chain(pos, cigartuples):
    """intron chain of a record, in gtf_names.chain's convention; pos = 0-based reference start. Only N opens an intron."""
    out, cur = [], pos
    for op, n in cigartuples:
        if op == 3:
            out.append((cur, cur + n + 1))
        if op in (0, 2, 3, 7, 8):
            cur += n
    return tuple(out)

def load_nodes(gtf):
    """-> nodes {locus: {chain: node_id}}, tx {node_id: [transcript ids]}; multi-exon transcripts, locus = gene_id."""
    ex = collections.defaultdict(list)
    for l in opn(gtf):
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            ex[(re.search(r'gene_id "([^"]+)"', f[8]).group(1), re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
    by = collections.defaultdict(list)
    for (g, t), e in ex.items():
        if len(e) > 1:
            by[(g, chain(e))].append(t)
    nodes, tx = collections.defaultdict(dict), {}
    for (g, c), ts in by.items():
        nodes[g][c] = min(ts); tx[min(ts)] = sorted(ts)
    return dict(nodes), tx

def member_of(locus, chain_, nodes):
    return "." if locus == OUTSIDE else UNSPLICED if not chain_ else nodes.get(locus, {}).get(chain_, OTHER)

def bands(records, window=10):
    """records {read: [(locus, member, gap, primary)]} -> {"exact": {read: frozenset of (locus, member)}, "wide": ...}.
    The primary record is in both bands whatever its score; a read with no record in a band is absent from it."""
    out = {"exact": {}, "wide": {}}
    for q, rs in records.items():
        for band, ok in (("exact", lambda g: g == 0), ("wide", lambda g: g <= window)):
            s = frozenset((l, m) for l, m, g, p in rs if p or ok(g))
            if s:
                out[band][q] = s
    return out

def is_node(m):
    return m[1] not in (UNSPLICED, OTHER, ".")

def layer(sets, base, repeat=frozenset(), min_link=1):
    """One band. sets {read: frozenset of members}, base Counter {member: untied molecules} -> dict:
    placed, classes, ft {node: key of its family transcript}, inside / open {key: molecules}, lo / hi {member: n},
    class_ft {class: (keys it touches, inside?)}, collision (keys of family transcripts holding two nodes of one locus)."""
    placed, classes = collections.Counter(base), collections.Counter()
    for s in sets.values():
        if len(s) == 1:
            placed[next(iter(s))] += 1
        else:
            classes[s] += 1
    linkable = lambda m: is_node(m) and m[0] not in (OUTSIDE, UNSEEN) and m[0] not in repeat
    parent, pair = {}, collections.Counter()
    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for m in list(placed) + [m for c in classes for m in c]:
        if is_node(m):
            find(m)
    for c, n in classes.items():
        ms = sorted(m for m in c if linkable(m))
        for i, a in enumerate(ms):
            for b in ms[i + 1:]:
                pair[(a, b)] += n
    for (a, b), n in pair.items():
        if n >= min_link:
            parent[find(a)] = find(b)
    comp = collections.defaultdict(list)
    for m in list(parent):
        comp[find(m)].append(m)
    ft = {m: min(ms) for ms in comp.values() for m in ms}
    inside, opn_, hi, class_ft = collections.Counter(), collections.Counter(), collections.Counter(placed), {}
    for c, n in classes.items():
        keys = {ft[m] for m in c if m in ft}
        class_ft[c] = (sorted(keys), len(keys) == 1 and all(linkable(m) for m in c))
        if class_ft[c][1]:
            inside[next(iter(keys))] += n
        else:
            for k in keys:
                opn_[k] += n
        for m in c:
            hi[m] += n
    collision = {min(ms) for ms in comp.values() if len({m[0] for m in ms}) < len(ms)}
    return {"placed": placed, "classes": classes, "ft": ft, "inside": inside, "open": opn_, "lo": placed, "hi": hi, "class_ft": class_ft, "collision": collision}

def locus_counts(sets, base):
    """-> placed, max per locus, as eq_classes counts them: molecules with every band record on the locus, and with any.
    `unseen` members are dropped first (eq_classes never saw them); `outside` is a locus that gets no entry."""
    pl, mx = collections.Counter(), collections.Counter()
    for (l, _), n in base.items():
        if l != OUTSIDE:
            pl[l] += n; mx[l] += n
    for s in sets.values():
        loci = {l for l, _ in s if l != UNSEEN}
        if len(loci) == 1:
            pl.update(loci - {OUTSIDE})
        mx.update(loci - {OUTSIDE})
    return pl, mx

def statuses(exact, wide, all_nodes):
    """per node, first match: certified (a molecule fits nowhere else within the window), best_only (a molecule scores
    best only here), tied_only (best-scoring records, all tied elsewhere), near_only (only records below their
    molecule's best), no_band_record."""
    return {n: "certified" if wide["lo"][n] else "best_only" if exact["lo"][n] else "tied_only" if exact["hi"][n]
            else "near_only" if wide["hi"][n] else "no_band_record" for n in all_nodes}

def min_copies(required, class_loci, cap=200000):
    """fewest loci that contain `required` and meet every set of class_loci -> (n, exact?). A lower bound on the copies
    expressed, never a call. Past `cap` tries: required + a greedy packing of disjoint unmet sets (still a lower bound)."""
    unmet = sorted({s for s in class_loci if not s & required}, key=lambda s: (len(s), sorted(s)))
    unmet = [s for s in unmet if not any(t < s for t in unmet)]
    if not unmet:
        return len(required), True
    cand, tries = sorted(set().union(*unmet)), 0
    for k in range(1, len(cand) + 1):
        for combo in itertools.combinations(cand, k):
            tries += 1
            if tries > cap:
                picked = []
                for s in unmet:
                    if not any(s & t for t in picked):
                        picked.append(s)
                return len(required) + len(picked), False
            if all(s & set(combo) for s in unmet):
                return len(required) + k, True

STATUSES = ("certified", "best_only", "tied_only", "near_only", "no_band_record")

def rows(path):
    with open(path) as fh:
        next(fh)
        for l in fh:
            yield l.rstrip("\n").split("\t")

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("gtf"); p.add_argument("prefix"); p.add_argument("--eq-loci"); p.add_argument("--window", type=float, default=10); p.add_argument("--min-link", type=int, default=1)
    a = p.parse_args()
    nodes, tx = load_nodes(a.gtf)
    lines, loc = list(opn(a.gtf)), {}
    for l in lines:
        f = l.split("\t")
        if len(f) > 8 and f[2] == "exon":
            c = loc.setdefault(re.search(r'transcript_id "([^"]+)"', f[8]).group(1), [f[0], int(f[3]), int(f[4]), f[6]])
            c[1], c[2] = min(c[1], int(f[3])), max(c[2], int(f[4]))
    all_nodes = sorted(((g, n) for g, d in nodes.items() for n in d.values()), key=lambda m: (loc[m[1]][0], loc[m[1]][1], m))
    records, base = collections.defaultdict(list), collections.Counter()
    for f in rows(a.prefix + ".reads.tsv"):
        records[f[0]].append((f[1], f[2], float(f[3]), f[6] == "1"))
    for f in rows(a.prefix + ".placed.tsv"):
        base[(f[0], f[1])] += int(f[2])
    known = set(all_nodes)
    bad = next((m for m in itertools.chain(base, ((l, m) for rs in records.values() for l, m, _, _ in rs)) if is_node(m) and m not in known), None)
    if bad:
        sys.exit(f"famwise: {a.prefix}.reads.tsv / .placed.tsv name transcript {bad[1]} of locus {bad[0]}, which {a.gtf} does not hold: were they made with another GTF?")
    eq = {r["locus"]: r for r in csv.DictReader(open(a.eq_loci), delimiter="\t")} if a.eq_loci else {}
    sets = bands(records, a.window)
    L = {b: layer(sets[b], base, frozenset(l for l, r in eq.items() if "repeat" in (r["label"], r["label_w"])), a.min_link) for b in sets}
    for b, sfx in (("exact", ""), ("wide", "_w")) if a.eq_loci else ():  # P1: every molecule counted once, exactly as eqclasses counts it per locus
        pl, mx = locus_counts(sets[b], base)
        for l in sorted(set(eq) | set(mx)):
            want = (int(eq[l]["placed" + sfx]), int(eq[l]["max" + sfx])) if l in eq else (0, 0)
            if want != (pl[l], mx[l]):
                sys.exit(f"famwise: P1 fails at locus {l} ({b}): eq placed/max {want[0]}/{want[1]}, famwise {pl[l]}/{mx[l]} (is {a.eq_loci} of the same run and --window?)")
    st, attr = statuses(L["exact"], L["wide"], all_nodes), {}
    with open(a.prefix + ".tx.tsv", "w") as ftx, open(a.prefix + ".classes.tsv", "w") as fcl, open(a.prefix + ".families.tsv", "w") as ffa:
        ftx.write("band\tnode\tlocus\tchrom\tstart\tend\tstrand\ttranscripts\tft\tlo\thi\tstatus\n")
        fcl.write("band\treads\tsize\tmembers\tft\tkind\n")
        ffa.write("band\tft\tnodes\tloci\tcopies\tmolecules\tplaced\tinside\topen\t" + "\t".join(STATUSES) + "\tmin_copies\tmin_copies_exact\tcollision\trep\n")
        for b in ("exact", "wide"):
            l, comp, inside_loci = L[b], collections.defaultdict(list), collections.defaultdict(list)
            for n in all_nodes:
                comp[l["ft"].get(n, n)].append(n)  # a transcript no band record reaches is a family transcript of its own
            mol = {k: sum(l["lo"][n] for n in ns) + l["inside"][k] for k, ns in comp.items()}
            order = sorted(comp, key=lambda k: (-mol[k], k))
            name = {k: f"FT{i + 1}" for i, k in enumerate(order)}
            for c, n in sorted(l["classes"].items(), key=lambda kv: (-kv[1], sorted(kv[0]))):
                keys, ins = l["class_ft"][c]
                if ins:
                    inside_loci[keys[0]].append(frozenset(m[0] for m in c))
                fcl.write(f"{b}\t{n}\t{len(c)}\t{','.join('|'.join(m) for m in sorted(c))}\t{','.join(name[k] for k in keys) or '.'}\t{'inside' if ins else 'open'}\n")
            for k in order:
                ns, cnt = comp[k], collections.Counter(st[n] for n in comp[k])
                loci = sorted({n[0] for n in ns})
                mc, exact = min_copies({n[0] for n in ns if st[n] in STATUSES[:2]}, inside_loci[k])
                rep = min(ns, key=lambda n: (-l["lo"][n], loc[n[1]][1], n[1]))
                ffa.write("\t".join(map(str, [b, name[k], len(ns), len(loci), ",".join(loci), mol[k], mol[k] - l["inside"][k], l["inside"][k], l["open"][k],
                                               *(cnt[x] for x in STATUSES), mc, "yes" if exact else "no", "yes" if k in l["collision"] else "no", rep[1]])) + "\n")
                for n in ns:
                    ftx.write("\t".join(map(str, [b, n[1], n[0], *loc[n[1]], ",".join(tx[n[1]]), name[k], l["lo"][n], l["hi"][n], st[n]])) + "\n")
                    if b == "wide":
                        for t in tx[n[1]]:
                            attr[t] = f' fw_tx "{name[k]}"; fw_status "{st[n]}"; fw_lo "{l["lo"][n]}"; fw_hi "{l["hi"][n]}"; fw_copies "{len(loci)}"; fw_rep "{"yes" if n == rep else "no"}";'
            multi = [k for k in order if len(comp[k]) > 1]
            ins, opn_ = (sum(n for c, n in l["classes"].items() if l["class_ft"][c][1] == x) for x in (True, False))
            print(f"{b}{'' if b == 'exact' else f' (within {a.window:g} of best)'}: {len(multi)} family transcripts join {sum(len(comp[k]) for k in multi)} transcripts of "
                  f"{len({n[0] for k in multi for n in comp[k]})} loci; {ins} tied molecules inside one family transcript, {opn_} open; {len(l['collision'])} with two transcripts of one locus")
    print("transcripts by status: " + ", ".join(f"{x} {sum(v == x for v in st.values())}" for x in STATUSES))
    with open(a.prefix + ".gtf", "w") as fh:
        for l in lines:
            f = l.split("\t")
            t = re.search(r'transcript_id "([^"]+)"', f[8]) if len(f) > 8 and f[2] == "transcript" else None
            fh.write(l.rstrip("\n") + attr[t.group(1)] + "\n" if t and t.group(1) in attr else l)

if __name__ == "__main__":
    main()
