#!/usr/bin/env python3
"""Phase 0b: build the RUSTLE_LOWINTRON_ORACLE file from ST pred_intron_low events.
File format (per line, whitespace-separated): "<strand> <comma-FREE mask> <chain>".
The loader (transcript_filter.rs:1549) parses the mask PER-CHARACTER, so commas must be stripped.
Usage: python3 bench/build_lowintron_oracle.py /tmp/ri_st.jsonl /tmp/lowintron.oracle"""
import json, sys

def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "/tmp/ri_st.jsonl"
    out = sys.argv[2] if len(sys.argv) > 2 else "/tmp/lowintron.oracle"
    n = 0
    with open(out, "w") as fh:
        for line in open(src):
            if '"step":"pred_intron_low"' not in line: continue
            try: r = json.loads(line)
            except ValueError: continue
            p = r.get("payload", {})
            chain = p.get("chain", "")
            mask = p.get("intron_low", "")
            if not chain or not mask: continue
            mask_nc = mask.replace(",", "")
            if len(mask_nc) != chain.count(",") + 1:
                continue
            fh.write(f"{r['strand']} {mask_nc} {chain}\n")
            n += 1
    print(f"wrote {n} oracle entries to {out}")

if __name__ == "__main__":
    main()
