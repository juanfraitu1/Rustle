# Avoiding Assembly Artifacts in Multi-Copy Gene Families

## The Problem: Copies That Are Too Close Together

When gene copies share >95% sequence identity, standard assemblers face
two failure modes:

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '13px'}}}%%
flowchart TD
    subgraph PROBLEM["Two Failure Modes"]
        direction TB
        subgraph MERGE["Failure 1: Erroneous Merging"]
            direction LR
            MA["Copy A\n(10 reads)"] 
            MB["Copy B\n(8 reads)"]
            MC["Merged Chimera\n(18 reads)\nWRONG: mixes exons\nfrom both copies"]
            MA -->|"shared exons\noverlap"| MC
            MB -->|"shared exons\noverlap"| MC
        end
        subgraph MISS["Failure 2: Copy Dropout"]
            direction LR
            DA["Copy A\n(10 reads)\nassembled"]
            DB["Copy B\n(8 reads)\n!!! MISSED !!!\n(reads assigned to A)"]
            DA -.-|"multi-mappers\nstolen by A"| DB
        end
    end

    style MERGE fill:#ffebee,stroke:#c62828,stroke-width:2px
    style MISS fill:#ffebee,stroke:#c62828,stroke-width:2px
    style MC fill:#ffcdd2,stroke:#c62828
    style DB fill:#ffcdd2,stroke:#c62828
```

---

## How VG Mode Prevents Erroneous Merging

### Step 1: Family Discovery (Separate, Don't Merge)

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '12px'}}}%%
flowchart LR
    subgraph STANDARD["Standard Assembler"]
        direction TB
        S1["All reads at chr7:100-250kb\ngo into ONE big bundle"]
        S2["Single splice graph\nmixes both copies"]
        S3["Chimeric transcripts\n(wrong junctions combined)"]
        S1 --> S2 --> S3
    end

    subgraph VG_MODE["Rustle VG Mode"]
        direction TB
        V1["Reads clustered into\nSEPARATE bundles:\n  Bundle A (chr7:100-150kb)\n  Bundle B (chr7:200-250kb)"]
        V2["Bundles LINKED as family\nvia shared multi-mapped reads\n(not merged!)"]
        V3["Each copy assembled\nindependently with\ncorrected read weights"]
        V1 --> V2 --> V3
    end

    style STANDARD fill:#ffebee,stroke:#c62828,stroke-width:2px
    style VG_MODE fill:#e8f5e9,stroke:#2E7D32,stroke-width:2px
    style S3 fill:#ffcdd2,stroke:#c62828
```

**Key principle:** VG mode **links** related bundles but keeps them as
**separate assembly units**. Multi-mapped reads are redistributed with
corrected weights — they are never physically merged into one bundle.

### Step 2: Junction-Based Copy Discrimination

Even when copies share most exons, they often differ at one or two
splice junctions. VG mode uses these differences:

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '12px'}}}%%
flowchart LR
    subgraph COPIES["Two copies differing by one exon"]
        direction TB
        subgraph CA["Copy A"]
            direction LR
            CA1["E1"] --- CA2["E2"] --- CA3["<b>E3a</b>"] --- CA4["E4"]
        end
        subgraph CB["Copy B"]
            direction LR
            CB1["E1"] --- CB2["E2"] --- CB3["<b>E3b</b>\n(3bp shorter)"] --- CB4["E4"]
        end
    end

    subgraph READS_J["Read Junction Evidence"]
        direction TB
        RU["Multi-mapped reads\n(E1-E2: shared junction)\nCannot distinguish copies"]
        RD1["Unique to A:\nE2---E3a junction\n(donor at pos 5000)"]
        RD2["Unique to B:\nE2---E3b junction\n(donor at pos 5003)"]
    end

    subgraph ASSIGN["VG Assignment"]
        direction TB
        A1["Reads with E2-E3a junction\n--> 100% to Copy A"]
        A2["Reads with E2-E3b junction\n--> 100% to Copy B"]
        A3["Reads with E1-E2 only\n--> split by EM/MFLP/Flow\nproportional to unique evidence"]
    end

    COPIES --> READS_J --> ASSIGN

    style CA3 fill:#c8e6c9,stroke:#2E7D32,stroke-width:3px
    style CB3 fill:#bbdefb,stroke:#1565C0,stroke-width:3px
    style RD1 fill:#c8e6c9,stroke:#2E7D32
    style RD2 fill:#bbdefb,stroke:#1565C0
    style RU fill:#fff9c4,stroke:#F57F17
```

**Even a single nucleotide difference** at a splice junction creates a
unique junction coordinate that unambiguously assigns a read to one copy.
The compatibility score in EM/MFLP uses these junction differences to
weight assignments.

---

## Recovering Missing Gene Family Members

### How Unmapped Reads Find Their Family

Reads from a paralog absent from the reference fail to align. But they
share sequence (k-mers) with known family members. VG mode exploits this:

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '12px'}}}%%
flowchart TD
    subgraph BAM["BAM File (post-alignment)"]
        direction LR
        MAPPED["Mapped reads\n(primary + supplementary)"]
        UNMAPPED["Unmapped reads\n(flag 0x4)"]
    end

    subgraph STEP1["Step 1: Build Family K-mer Index"]
        direction TB
        KI["For each gene family:\n  extract k-mers from\n  all assembled transcripts\n  across all copies"]
    end

    subgraph STEP2["Step 2: Scan Unmapped Reads"]
        direction TB
        SCAN["For each unmapped read:\n  compute k-mer overlap\n  with each family's index"]
        MATCH["If overlap >= threshold\n(e.g., Jaccard >= 0.3):\n  assign read to family"]
    end

    subgraph STEP3["Step 3: Cluster & Assemble"]
        direction TB
        CLUSTER["Cluster rescued reads\nby approximate mapping\nposition (supplementary\nalignment coordinates)"]
        NOVEL["If cluster has >= N reads\nwith novel junctions:\n  create synthetic bundle\n  --> assemble as new copy"]
    end

    BAM --> STEP1 --> STEP2 --> STEP3

    style UNMAPPED fill:#fff9c4,stroke:#F57F17,stroke-width:2px
    style NOVEL fill:#fff9c4,stroke:#F57F17,stroke-width:3px
    style STEP2 fill:#fff3e0,stroke:#FF9800
```

**The reads never passed through the mapper for this locus** — they were
unmapped or mapped elsewhere as supplementary alignments. VG mode
re-discovers them through sequence similarity (k-mer matching) to the
family's assembled transcripts, bypassing the aligner's reference bias.

---

## Summary: VG Mode Artifact Prevention

```
    Artifact              | Cause                    | VG Prevention
    ----------------------|--------------------------|---------------------------
    Chimeric transcripts  | Merging reads from       | Keep copies as separate
                          | different copies into    | bundles; link via family
                          | one splice graph         | grouping, not merging
                          |                          |
    Copy dropout          | Multi-mappers stolen     | EM/MFLP/Flow redistribute
                          | by dominant copy         | weights using junction
                          | (winner-takes-all)       | compatibility evidence
                          |                          |
    Missing paralogs      | Novel copy absent from   | K-mer scan of unmapped
                          | reference; reads can't   | reads against family VG;
                          | align                    | create novel bundles
                          |                          |
    SNP-identical copies  | Copies differ only by    | --vg-snp mode: parse MD
                          | point mutations, not     | tag for per-base variants;
                          | splice junctions         | build diagnostic SNP
                          |                          | profiles per copy
                          |                          |
    Haplotype confusion   | Two haplotypes of same   | --vg-phase: use HP/PS
                          | copy create false        | tags to split reads by
                          | isoform diversity        | haplotype before assembly
```
