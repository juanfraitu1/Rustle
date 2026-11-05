#!/usr/bin/env python3
"""Opt-in family HMM layer (docs/superpowers/specs/2026-11-05-family-hmm-design.md): one DNA and one protein
profile per multi-copy family, Pfam labels, and a pre-registered concordance table against the MCL families.
Not a member-discovery or assignment emitter (unmapped-rescue Amendment 5 registered that arm)."""
import collections
import os
import re
import subprocess
import sys


def read_fasta(path):
    """[(header without '>', seq)]"""
    out, hdr, seq = [], None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if hdr is not None:
                    out.append((hdr, "".join(seq)))
                hdr, seq = ln[1:].strip(), []
            else:
                seq.append(ln.strip())
    if hdr is not None:
        out.append((hdr, "".join(seq)))
    return out


def copies_by_family_records(records):
    """{family: [(header, seq)]} keyed by the header's first '|' field"""
    fams = collections.defaultdict(list)
    for h, s in records:
        fams[h.split("|")[0]].append((h, s))
    return dict(fams)


def copies_by_family(copies_fa):
    return copies_by_family_records(read_fasta(copies_fa))


def longest_orf(records):
    """esl-translate records of ONE copy -> (header, frame, 0-based nt offset, aa) of the longest, or None"""
    if not records:
        return None
    hdr, aa = max(records, key=lambda r: len(r[1]))
    m = re.search(r"coords=(\d+)\.\.(\d+) .*frame=(\d+)", hdr)
    return hdr, int(m.group(3)), int(m.group(1)) - 1, aa


def esl_translate(copy_fa, min_aa):
    """raw FASTA text of the top-strand ORFs of at least min_aa residues"""
    return subprocess.run(["esl-translate", "-l", str(min_aa), "--watson", copy_fa],
                          capture_output=True, text=True, check=True).stdout


def _run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, check=True, **kw)


def write_fasta(path, records):
    with open(path, "w") as fh:
        for h, s in records:
            fh.write(">%s\n%s\n" % (h, s))


def build_profile(fasta, hmm_out, name, alphabet, threads=1):
    """mafft --auto then hmmbuild --<alphabet> -n name"""
    aln = hmm_out + ".aln.fa"
    with open(aln, "w") as fh:
        subprocess.run(["mafft", "--auto", "--quiet", "--thread", str(threads), fasta], stdout=fh, check=True)
    _run(["hmmbuild", "--" + alphabet, "-n", name, hmm_out, aln])
    os.remove(aln)


def select_families(fams, only, min_copies):
    """family ids with >= min_copies copies, restricted to `only` (ids absent from the table go to stderr)"""
    if only is not None:
        for f in sorted(only - set(fams)):
            print("skip %s: not in copies" % f, file=sys.stderr)
    return sorted(f for f, c in fams.items() if len(c) >= min_copies and (only is None or f in only))


def _concat_press(out, sub, name):
    paths = sorted(os.path.join(out, sub, f) for f in os.listdir(os.path.join(out, sub)) if f.endswith(".hmm"))
    db = os.path.join(out, name)
    with open(db, "w") as fh:
        for p in paths:
            fh.write(open(p).read())
    if paths:
        _run(["hmmpress", "-f", db])
    return db


def build(prefix, out, only=None, min_copies=2, min_orf_aa=50, threads=1):
    """DNA + protein profile per family -> out/{dna,prot}/FAM.hmm, out/orf.tsv, out/build.tsv, pressed out/families.{dna,prot}.hmm"""
    import tempfile
    fams = copies_by_family(prefix + ".fam.copies.fa")
    for sub in ("dna", "prot"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    rows, orf_rows = [], []
    with tempfile.TemporaryDirectory() as tmp:
        for fam in select_families(fams, only, min_copies):
            copies = fams[fam]
            nt_fa = os.path.join(tmp, fam + ".nt.fa")
            write_fasta(nt_fa, copies)
            build_profile(nt_fa, os.path.join(out, "dna", fam + ".hmm"), fam, "dna", threads)
            prots = []
            for h, s in copies:
                one = os.path.join(tmp, "copy.fa")
                write_fasta(one, [(h.split("|")[2], s)])
                recs = _parse_fasta_text(esl_translate(one, min_orf_aa))
                o = longest_orf(recs)
                if o is None:
                    continue
                _hdr, frame, off, aa = o
                orf_rows.append((fam, h.split("|")[2], frame, off, len(aa)))
                prots.append((h, aa))
            status = "ok" if len(prots) >= 2 else "none"
            if status == "ok":
                aa_fa = os.path.join(tmp, fam + ".aa.fa")
                write_fasta(aa_fa, prots)
                build_profile(aa_fa, os.path.join(out, "prot", fam + ".hmm"), fam, "amino", threads)
            rows.append({"family": fam, "n_copies": len(copies), "prot_status": status})
    with open(os.path.join(out, "orf.tsv"), "w") as fh:
        fh.write("family\ttid\tframe\tnt_offset\taa_len\n")
        for r in orf_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(out, "build.tsv"), "w") as fh:
        fh.write("family\tn_copies\tprot_status\n")
        for r in rows:
            fh.write("%s\t%d\t%s\n" % (r["family"], r["n_copies"], r["prot_status"]))
    _concat_press(out, "dna", "families.dna.hmm")
    _concat_press(out, "prot", "families.prot.hmm")
    return rows


def _parse_fasta_text(text):
    out, hdr, seq = [], None, []
    for ln in text.splitlines():
        if ln.startswith(">"):
            if hdr is not None:
                out.append((hdr, "".join(seq)))
            hdr, seq = ln[1:].strip(), []
        elif hdr is not None:
            seq.append(ln.strip())
    if hdr is not None:
        out.append((hdr, "".join(seq)))
    return out


WS = "/anvme/workspace/v138da13-rustle"
PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"


def parse_domtblout(lines, max_evalue=1e-5):
    """hmmscan --domtblout -> {query: [(ali_from, pfam_name, i_evalue)] sorted by ali_from}"""
    out = collections.defaultdict(list)
    for ln in lines:
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split()
        iev = float(f[12])
        if iev <= max_evalue:
            out[f[3]].append((int(f[17]), f[0], iev))
    return {q: sorted(v) for q, v in out.items()}


def pfam_path(arg=None):
    """Pfam-A.hmm path or exit 2 with the download command"""
    p = arg or os.path.join(WS, "pfam", "Pfam-A.hmm")
    if os.path.exists(p + ".h3m"):
        return p
    d = os.path.dirname(p)
    print("Pfam-A not found at %s. Download once:\n  wget -P %s %s && gunzip %s/Pfam-A.hmm.gz && hmmpress %s"
          % (p, d, PFAM_URL, d, p), file=sys.stderr)
    sys.exit(2)


def scan(out, pfam, max_evalue=1e-5, threads=1):
    """Pfam domain labels of each family's protein consensus -> out/pfam_labels.tsv"""
    cons = os.path.join(out, "prot.consensus.fa")
    with open(cons, "w") as fh:
        subprocess.run(["hmmemit", "-c", os.path.join(out, "families.prot.hmm")], stdout=fh, check=True)
    dom = os.path.join(out, "pfam.domtblout")
    _run(["hmmscan", "--cpu", str(threads), "--noali", "--domtblout", dom, "-o", os.devnull, pfam, cons])
    labels = parse_domtblout(open(dom), max_evalue)
    with open(os.path.join(out, "pfam_labels.tsv"), "w") as fh:
        fh.write("family\tdomains\tbest_evalue\tprot_status\n")
        for ln in open(os.path.join(out, "build.tsv")).read().splitlines()[1:]:
            fam, _n, status = ln.split("\t")
            hits = labels.get(fam + "-consensus", [])
            fh.write("%s\t%s\t%s\t%s\n" % (fam, ";".join(h[1] for h in hits),
                                           min((h[2] for h in hits), default="NA"), status))
