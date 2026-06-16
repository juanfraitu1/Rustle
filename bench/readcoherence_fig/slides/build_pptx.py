#!/usr/bin/env python3
"""
Assemble the rendered slide PNGs + speaker notes into a 16:9 PPTX.

Each slideN.png is placed FULL-BLEED on a 16:9 slide; speaker notes are
attached per slide.

Run with: /home/juanfra/miniforge3/bin/python build_pptx.py
"""
import os
from pptx import Presentation
from pptx.util import Inches

HERE = os.path.dirname(os.path.abspath(__file__))

# 16:9 at 13.33 x 7.5 inches (matches the PNG figsize)
SLIDE_W = Inches(13.333)
SLIDE_H = Inches(7.5)

# Speaker notes per slide (1-indexed -> list index 0..6)
NOTES = [
    # 1 - title
    "",
    # 2 - fig2 mechanism
    "Flow assembly builds one splice graph from all reads and picks the "
    "FEWEST paths that explain coverage (parsimony) — so real but "
    "low-abundance / overlapping isoforms get merged or dropped. "
    "Read-coherence instead emits each molecule's own spliced path.",
    # 3 - fig1 igv locus
    "Real chr19 locus (NC_073243.2). The four StringTie variants splice a "
    "long intron straight over the internal exons; the read-coherence isoform "
    "retains them, and 15 spanning long-read alignments are the evidence. "
    "Eyeball-verifiable in IGV.",
    # 4 - rigor
    "Clustering needs a distance + threshold and can merge distinct "
    "molecules; this is a discrete equivalence on junction chains — "
    "parameter-free. Junctions are exact; only the 5'/3' ends are soft, "
    "handled by the outer-envelope + fragment-collapse so truncation never "
    "mints spurious short isoforms.",
    # 5 - fig4 gate collapse
    "Two precision mechanisms, both annotation-free. The realness gate keeps "
    "a read-chain isoform only if every junction is canonical, it has no "
    "RT-switch repeat signature, and read-depth >= min_cov. Degradation-aware "
    "collapse folds 5'/3'/internal fragments into their full-length parent "
    "(ISM -> FSM).",
    # 6 - fig3 architecture
    "Read-chain is held out of every per-bundle and global filter that could "
    "displace a flow transcript, then unioned back only as NOVEL intron chains "
    "and passed through the gate. So the output is a provable superset of the "
    "flow baseline: genome-wide gd\\rcg = 0 and st\\rcg = 0 (never loses a "
    "flow find or a StringTie guide).",
    # 7 - fig5 results
    "Genome-wide (25 chromosomes). At the shipped default min_cov=3 (the "
    "precision/recall knee), read-coherence recovers +1,735 exact (FSM) "
    "annotated isoforms StringTie misses — about 3x the raw per-molecule "
    "approach. Of the gated read-chain extras: 100% canonical, 63% real, "
    "12% noise.",
]


def main():
    prs = Presentation()
    prs.slide_width = SLIDE_W
    prs.slide_height = SLIDE_H

    blank_layout = prs.slide_layouts[6]  # fully blank

    for i in range(1, 8):
        png = os.path.join(HERE, f"slide{i}.png")
        if not os.path.exists(png):
            raise FileNotFoundError(png)

        slide = prs.slides.add_slide(blank_layout)

        # full-bleed image
        slide.shapes.add_picture(png, 0, 0, width=SLIDE_W, height=SLIDE_H)

        # speaker notes
        note = NOTES[i - 1]
        if note:
            slide.notes_slide.notes_text_frame.text = note

    out = os.path.join(HERE, "read_coherence_deck.pptx")
    prs.save(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
