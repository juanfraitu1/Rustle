import csv
from pathlib import Path

OUT = Path(__file__).parent

ROWS = [
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","bundle_run":1,"iter_idx":1,"path_len":4,"path_nodes":"S,0,1,3,T","bottleneck":12.5,"total_flow_after":12.5},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","bundle_run":1,"iter_idx":2,"path_len":5,"path_nodes":"S,0,1,2,3,T","bottleneck":2.0,"total_flow_after":14.5},
]

with (OUT/"flow_3iter.tsv").open("w") as f:
    w = csv.DictWriter(f, fieldnames=list(ROWS[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(ROWS)
