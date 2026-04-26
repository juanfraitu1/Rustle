import csv
from pathlib import Path

OUT = Path(__file__).parent
NODE_FIELDS = ["source","chrom","bdstart","bdend","strand","node_idx","start","end","cov","hardstart","hardend","nchildren","nparents"]
EDGE_FIELDS = ["source","chrom","bdstart","bdend","strand","from_idx","to_idx","from_start","from_end","to_start","to_end","edge_kind","from_cov","to_cov","bottleneck_cov"]


def make_chain_graph(copy: str, offset: int):
    """4-node linear-with-skip chain, all junctions, fixed shape."""
    nodes = []
    for i, (s, e, cov) in enumerate([(100,200,10),(400,500,11),(700,800,9),(1000,1100,12)]):
        nodes.append({
            "source":"rustle","chrom":copy,"bdstart":offset,"bdend":offset+1100,"strand":"+",
            "node_idx":i,"start":offset+s,"end":offset+e,"cov":cov,
            "hardstart":1 if i==0 else 0,"hardend":1 if i==3 else 0,
            "nchildren":2 if i==1 else (1 if i<3 else 0),
            "nparents":2 if i==3 else (1 if i>0 else 0),
        })
    edges = []
    for (a,b) in [(0,1),(1,2),(1,3),(2,3)]:
        edges.append({
            "source":"rustle","chrom":copy,"bdstart":offset,"bdend":offset+1100,"strand":"+",
            "from_idx":a,"to_idx":b,
            "from_start":nodes[a]["start"],"from_end":nodes[a]["end"],
            "to_start":nodes[b]["start"],"to_end":nodes[b]["end"],
            "edge_kind":"junction",
            "from_cov":nodes[a]["cov"],"to_cov":nodes[b]["cov"],
            "bottleneck_cov":min(nodes[a]["cov"], nodes[b]["cov"]),
        })
    return nodes, edges


def write(name, fields, rows):
    with (OUT/name).open("w") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)


# Iso triple
all_n, all_e = [], []
for copy in ("copy1","copy2","copy3"):
    n,e = make_chain_graph(copy, 0)
    all_n += n; all_e += e
write("iso_triple_nodes.tsv", NODE_FIELDS, all_n)
write("iso_triple_edges.tsv", EDGE_FIELDS, all_e)

# Non-iso pair: copyA = 4-node chain, copyB = 5-node linear (no skip)
n_a, e_a = make_chain_graph("copyA", 0)

# 5-node linear copyB
n_b = []
for i,(s,e,cov) in enumerate([(100,200,10),(300,400,11),(500,600,9),(700,800,10),(900,1000,12)]):
    n_b.append({
        "source":"rustle","chrom":"copyB","bdstart":0,"bdend":1000,"strand":"+",
        "node_idx":i,"start":s,"end":e,"cov":cov,
        "hardstart":1 if i==0 else 0,"hardend":1 if i==4 else 0,
        "nchildren":1 if i<4 else 0,
        "nparents":1 if i>0 else 0,
    })
e_b = []
for a,b in [(0,1),(1,2),(2,3),(3,4)]:
    e_b.append({
        "source":"rustle","chrom":"copyB","bdstart":0,"bdend":1000,"strand":"+",
        "from_idx":a,"to_idx":b,
        "from_start":n_b[a]["start"],"from_end":n_b[a]["end"],
        "to_start":n_b[b]["start"],"to_end":n_b[b]["end"],
        "edge_kind":"junction",
        "from_cov":n_b[a]["cov"],"to_cov":n_b[b]["cov"],
        "bottleneck_cov":min(n_b[a]["cov"], n_b[b]["cov"]),
    })
write("non_iso_nodes.tsv", NODE_FIELDS, n_a + n_b)
write("non_iso_edges.tsv", EDGE_FIELDS, e_a + e_b)
