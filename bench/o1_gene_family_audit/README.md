# O1 typed gene-family audit

This audit treats Soto families as segmental-duplication discovery strata, not as gene-family truth.
RNA forward homology supplies primary membership; DNA supplies corroboration and can rescue an
independently named RNA-null member, but DNA-only SD homology cannot recruit an anonymous locus.

## Known multi-copy families

| family | named members | supported | primary graph (nodes/edges/components) | RNA edges | DNA edges | RNA-only / both / DNA-only | status |
|---|---:|---:|---:|---:|---:|---:|---|
| NPIP | 19 | 19 | 26/280/26 | 235 | 379 | 3 / 232 / 147 | PASS |
| TBC1D3 | 9 | 9 | 11/55/11 | 55 | 57 | 0 / 55 / 2 | PASS |
| RABL2 | 2 | 2 | 2/1/2 | 1 | 1 | 0 / 1 / 0 | PASS |
| APOBEC3 | 7 | 6 | 7/6/6+1 | 0 | 6 | 0 / 0 / 6 | PARTIAL_KNOWN_OUTLIER |
| MAGEA | 12 | 12 | 12/62/12 | 60 | 58 | 16 / 44 / 14 | PASS |
| GSTM | 4 | 3 | 4/3/3+1 | 1 | 3 | 0 / 1 / 2 | PARTIAL_KNOWN_OUTLIER |
| HERC2 | 1 | 1 | 10/16/10 | 16 | 23 | 6 / 10 / 13 | ANNOTATION_LIMITED |

The panel removes 3 adjudicated contaminants and withholds 1 independently
conflicting gene annotations from primary membership. GSTM3 and APOBEC3H remain named
members but are reported as graph false negatives rather than being relabeled or deleted.
HERC2 has a connected ten-node candidate graph but only one independently named member in
this annotation, so it is annotation-limited rather than counted as a validated pass.

## Soto-missing loci

Independent annotation supports 4 unique Soto-missing gene copies.
A further 3 loci are novel gene-copy candidates, not confirmed copies, because they
lack independent same-family annotation.

| locus | SD family | independent annotation | disposition |
|---|---|---|---|
| chr1:120449518-120478614 | ID_400 | NBPF8 | SUPPORTED_SOTO_MISSING_GENE_COPY |
| chr1:148601552-148620972 | ID_400 | NBPF19 | SUPPORTED_SOTO_MISSING_GENE_COPY |
| chr15:20812032-20814214 | ID_116 | GOLGA6L1 | SUPPORTED_SOTO_MISSING_GENE_COPY |
| chr9:77863432-77881287 | ID_280 | LOC128966611 | SUPPORTED_SOTO_MISSING_GENE_COPY |
| chr1:16478340-16479683 | ID_393 | . | NOVEL_GENE_COPY_CANDIDATE |
| chr15:21025273-21028497 | ID_71 | LOC102723564 | SUPPORTED_NONCODING_FAMILY_LOCUS |
| chr15:28174361-28187051 | ID_179 | . | REJECT_GENE_FAMILY_ASSIGNMENT |
| chr5:70456936-70462706 | ID_163 | . | NOVEL_GENE_COPY_CANDIDATE |
| chr7:76085482-76092920 | ID_245 | . | REJECT_NOT_MULTICOPY |
| chr9:77926044-77949112 | ID_280 | . | NOVEL_GENE_COPY_CANDIDATE |

The full TSV retains SD-only and partial-homology rows so they remain reviewable without
inflating the gene-family copy count. `<family>.gene_family.gfa` contains only primary members
and membership edges; `<family>.audit.gfa` retains every SD-discovery node for inspection. GFA
link tag `EV` records RNA/DNA evidence and audit tag `MB:i:1` marks admitted edges. The companion
Bandage colour CSV uses green for known members, blue for RNA-recruited candidates, orange for
known false negatives, grey/
purple for review nodes, and red for adjudicated contaminants.
