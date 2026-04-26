import csv
from pathlib import Path

OUT = Path(__file__).parent

NODES = [
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","node_idx":0,"start":1000,"end":1200,"cov":12.5,"hardstart":1,"hardend":0,"nchildren":1,"nparents":0},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","node_idx":1,"start":2000,"end":2300,"cov":15.0,"hardstart":0,"hardend":0,"nchildren":2,"nparents":1},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","node_idx":2,"start":3000,"end":3150,"cov":11.0,"hardstart":0,"hardend":0,"nchildren":1,"nparents":1},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","node_idx":3,"start":4000,"end":4500,"cov":13.0,"hardstart":0,"hardend":1,"nchildren":0,"nparents":2},
]
EDGES = [
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","from_idx":0,"to_idx":1,"from_start":1000,"from_end":1200,"to_start":2000,"to_end":2300,"edge_kind":"junction","from_cov":12.5,"to_cov":15.0,"bottleneck_cov":12.5},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","from_idx":1,"to_idx":2,"from_start":2000,"from_end":2300,"to_start":3000,"to_end":3150,"edge_kind":"junction","from_cov":15.0,"to_cov":11.0,"bottleneck_cov":11.0},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","from_idx":1,"to_idx":3,"from_start":2000,"from_end":2300,"to_start":4000,"to_end":4500,"edge_kind":"junction","from_cov":15.0,"to_cov":13.0,"bottleneck_cov":13.0},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","from_idx":2,"to_idx":3,"from_start":3000,"from_end":3150,"to_start":4000,"to_end":4500,"edge_kind":"junction","from_cov":11.0,"to_cov":13.0,"bottleneck_cov":11.0},
]

def write(name, fields, rows):
    p = OUT / name
    with p.open("w") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

write("small_nodes.tsv", list(NODES[0].keys()), NODES)
write("small_edges.tsv", list(EDGES[0].keys()), EDGES)
