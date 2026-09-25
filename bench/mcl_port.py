#!/usr/bin/env python3
"""Python MCL comparator — now a thin shim over the bit-faithful Rust bin `mcl_port` (§6z3, r1047).

Kept because the bench scorers call it (since wave 7 through `lib.mcl`) as `mcl(edges, ...)`; their numbers were designed
against THIS algorithm (register 917: deliberately not bit-identical to the shipped
`annotation_families::mcl`, 245 vs 168 clusters on chr2). The Rust bin reproduces the former
numpy/scipy implementation exactly — same clusters, same cluster order, same member order — verified on
20 graphs (10 real-weight + 10 uniform-weight, the tie-heavy case): 20/20 identical. numpy and scipy are
no longer imported anywhere on this path.

Semantics (unchanged): self-loops of weight 1, column normalisation, expansion M*M, inflation, absolute
prune, convergence at 1e-7, at most `max_iter` iterations; clusters = columns joined to their heaviest row.
"""
import os
import subprocess
import tempfile

_BIN = os.environ.get('RUSTLE_MCL_PORT_BIN',
                      '/mnt/linuxdisk/home/juanfraitu/rustle_target/release/mcl_port')


def mcl(edges, inflation=2.8, prune=1e-9, max_iter=100):
    """edges: {(a, b): weight} over STRING node ids (undirected). Returns clusters (lists of node ids), size >= 1,
    in the order the numpy implementation returned them (ascending smallest member; members ascending)."""
    if not edges:
        return []
    fd, path = tempfile.mkstemp(suffix='.mcl.tsv', dir=os.environ.get('TMPDIR'))
    try:
        with os.fdopen(fd, 'w') as fh:
            for (a, b), w in edges.items():
                if not (isinstance(a, str) and isinstance(b, str)):
                    raise TypeError(f'mcl_port.mcl: node ids must be str, got {type(a).__name__}/{type(b).__name__}')
                if '\t' in a or '\n' in a or '\t' in b or '\n' in b:
                    raise ValueError('mcl_port.mcl: node ids may not contain tab or newline')
                fh.write(f'{a}\t{b}\t{float(w)!r}\n')   # repr round-trips the f64 exactly
        out = subprocess.run([_BIN, '--graph', path, '--inflation', repr(float(inflation)),
                              '--prune', repr(float(prune)), '--max-iter', str(int(max_iter))],
                             capture_output=True, text=True, check=True).stdout
    finally:
        os.unlink(path)
    return [line.split('\t') for line in out.split('\n') if line]
