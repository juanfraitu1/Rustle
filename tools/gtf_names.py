#!/usr/bin/env python3
"""Name the transcripts of an assembled GTF after an annotation (driver stage `names`, option --names ANNOT).

    gtf_names.py ANNOT.gff3|gtf[.gz] IN.gtf OUT.gtf

The annotation is a source of names only: no transcript is added, dropped or changed, ids stay as they are. Appended to
the transcript and exon lines of every transcript that overlaps an annotated gene on its strand:
  ref_gene_name / ref_gene_id   the gene; for an exact match the matched transcript's gene, else the gene sharing the
                                most exon bases with the model ('.' strand models are compared with both strands)
  ref_transcript_id             only when the intron chain equals an annotated transcript's (multi-exon models)
  ref_match                     exact | overlap
  ref_genes                     every overlapped gene name, when more than one (a model across two genes)
Prints one summary line (the driver logs it)."""
import bisect, collections, gzip, re, sys

def opn(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def attrs(s):
    """-> (is_gtf, {key: value})"""
    if re.match(r'\s*\w+ "', s):
        return True, dict(re.findall(r'(\S+) "([^"]*)"', s))
    return False, dict(x.split("=", 1) for x in s.strip().split(";") if "=" in x)

def chain(ex):
    ex = sorted(ex)
    return tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))

def load_annot(path):
    """-> chains {(chrom, strand, chain): (gene_key, transcript_id)}, genes {gene_key: (id, name)},
    gex {(chrom, strand): [(start, end, gene_key)]} sorted exon intervals."""
    gname, parent, tex, tloc = {}, {}, collections.defaultdict(list), {}
    for l in opn(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or l[0] == "#":
            continue
        is_gtf, a = attrs(f[8])
        if is_gtf:
            g = a.get("gene_id", "")
            gname.setdefault(g, a.get("gene_name") or a.get("gene") or g)
            if a.get("gene_name"):
                gname[g] = a["gene_name"]
            if f[2] == "exon":
                t = a.get("transcript_id", g)
                parent[t] = g; tex[t].append((int(f[3]), int(f[4]))); tloc[t] = (f[0], f[6])
        elif f[2] == "exon":
            for t in a.get("Parent", "").split(","):
                tex[t].append((int(f[3]), int(f[4]))); tloc[t] = (f[0], f[6])
        elif "ID" in a:
            if "Parent" in a:
                parent[a["ID"]] = a["Parent"]
            else:
                gname[a["ID"]] = a.get("gene_name") or a.get("Name") or a.get("gene") or a["ID"]
    chains, gex = {}, collections.defaultdict(list)
    for t, ex in tex.items():
        g = parent.get(t, t)  # RefSeq pseudogene exons hang off the gene itself
        gname.setdefault(g, g)
        c, st = tloc[t]
        if len(ex) > 1:
            chains.setdefault((c, st, chain(ex)), (g, t))
        gex[(c, st)].extend((s, e, g) for s, e in ex)
    for k in gex:
        gex[k] = sorted(set(gex[k]))
    return chains, {g: (re.sub(r"^gene[-:]", "", g), n) for g, n in gname.items()}, gex

def overlap_genes(ex, ivs, starts, longest):
    """exon bases shared with each gene: {gene_key: bp}. ivs sorted by start; longest = the longest annotated exon."""
    out = collections.Counter()
    for s, e in ex:
        i = bisect.bisect_left(starts, s - longest)
        while i < len(ivs) and ivs[i][0] <= e:
            o = min(e, ivs[i][1]) - max(s, ivs[i][0]) + 1
            if o > 0:
                out[ivs[i][2]] += o
            i += 1
    return out

def main(annot, src, dst):
    chains, genes, gex = load_annot(annot)
    idx = {k: ([s for s, _, _ in v], max(e - s for s, e, _ in v)) for k, v in gex.items()}
    lines, tx, loc = opn(src).readlines(), collections.defaultdict(list), {}
    tid = lambda f: re.search(r'transcript_id "([^"]+)"', f[8])
    for l in lines:
        f = l.split("\t")
        if len(f) > 8 and f[2] == "exon" and tid(f):
            t = tid(f).group(1); tx[t].append((int(f[3]), int(f[4]))); loc[t] = (f[0], f[6])
    tag, n = {}, collections.Counter()
    for t, ex in tx.items():
        c, st = loc[t]
        strands = (st,) if st in "+-" else ("+", "-")
        hit = next((chains[(c, s, chain(ex))] for s in strands if len(ex) > 1 and (c, s, chain(ex)) in chains), None)
        ov = collections.Counter()
        for s in strands:
            if (c, s) in gex:
                ov.update(overlap_genes(ex, gex[(c, s)], *idx[(c, s)]))
        if hit:
            g, rt = hit
        elif ov:
            g, rt = min(ov, key=lambda k: (-ov[k], genes[k][1])), None
        else:
            n["unnamed"] += 1
            continue
        n["exact" if rt else "overlap"] += 1
        a = f' ref_gene_name "{genes[g][1]}"; ref_gene_id "{genes[g][0]}";'
        if rt:
            a += f' ref_transcript_id "{re.sub(r"^(rna|transcript)[-:]", "", rt)}";'
        a += f' ref_match "{"exact" if rt else "overlap"}";'
        names = sorted({genes[k][1] for k in ov} | {genes[g][1]})
        if len(names) > 1:
            a += f' ref_genes "{",".join(names)}";'
        tag[t] = a
    with open(dst, "w") as out:
        for l in lines:
            f = l.split("\t")
            m = tid(f) if len(f) > 8 else None
            out.write(l.rstrip("\n") + tag[m.group(1)] + "\n" if m and m.group(1) in tag else l)
    print(f"{len(tx)} transcripts: {n['exact']} exact intron chain, {n['overlap']} named by exon overlap, {n['unnamed']} outside the annotation")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
