# Variation Graph: Novel Copy Discovery and Read Rescue

## Gene Family Variation Graph with Unmapped Read Rescue

In a linear reference, reads from paralogs absent in the assembly have
nowhere to map and are lost. A **variation graph** encodes all known copies
of a gene family as alternate paths through shared and divergent exons,
allowing unmapped reads to be rescued by aligning to the graph.

### Conceptual Overview

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '14px'}, 'flowchart': {'nodeSpacing': 20, 'rankSpacing': 40}}}%%
flowchart LR

    subgraph LINEAR["Linear Reference (one copy only)"]
        direction LR
        L1["Exon 1"] --- li1((" ")) --- L2["Exon 2"] --- li2((" ")) --- L3["Exon 3"]
        li1 -.-|"intron"| li1
        li2 -.-|"intron"| li2
    end

    subgraph UNMAPPED["Unmapped Reads"]
        direction TB
        U1["Read X ????????\n(no alignment)"]
        U2["Read Y ????????\n(no alignment)"]
        U3["Read Z ????????\n(no alignment)"]
    end

    LINEAR -.->|"aligner finds\nno match"| UNMAPPED

    style LINEAR fill:#ffebee,stroke:#c62828,stroke-width:2px
    style UNMAPPED fill:#ffebee,stroke:#c62828,stroke-width:2px
    style U1 fill:#ffcdd2,stroke:#c62828
    style U2 fill:#ffcdd2,stroke:#c62828
    style U3 fill:#ffcdd2,stroke:#c62828
    style li1 fill:none,stroke:none
    style li2 fill:none,stroke:none
```

**Problem:** The reference genome contains only one copy of the gene.
Reads originating from a missing paralog cannot map and are discarded.

---

### The Variation Graph Solution

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '13px'}, 'flowchart': {'nodeSpacing': 15, 'rankSpacing': 50, 'curve': 'basis'}}}%%
flowchart LR

    START((" "))

    subgraph VG["Gene Family Variation Graph"]
        direction LR

        E1_shared["<b>Exon 1</b>\n(shared)"]

        E2_A["<b>Exon 2a</b>\nCopy A\n...ACGTGC..."]
        E2_B["<b>Exon 2b</b>\nCopy B\n...ACATGC...\n(1 SNP)"]
        E2_C["<b>Exon 2c</b>\nNovel Copy\n...ACTTGC...\n(2 SNPs)"]

        E3_shared["<b>Exon 3</b>\n(shared)"]

        E4_A["<b>Exon 4a</b>\nCopy A only"]
        E4_B["<b>Exon 4b</b>\nCopy B only"]

        E5_shared["<b>Exon 5</b>\n(shared)"]

        E1_shared -->|"intron"| E2_A
        E1_shared -->|"intron"| E2_B
        E1_shared -->|"intron"| E2_C

        E2_A -->|"intron"| E3_shared
        E2_B -->|"intron"| E3_shared
        E2_C -->|"intron"| E3_shared

        E3_shared -->|"intron"| E4_A
        E3_shared -->|"intron"| E4_B

        E4_A -->|"intron"| E5_shared
        E4_B -->|"intron"| E5_shared
    end

    STOP((" "))

    START --> E1_shared
    E5_shared --> STOP

    style VG fill:#e3f2fd,stroke:#1565C0,stroke-width:3px
    style E1_shared fill:#bbdefb,stroke:#1565C0
    style E3_shared fill:#bbdefb,stroke:#1565C0
    style E5_shared fill:#bbdefb,stroke:#1565C0
    style E2_A fill:#c8e6c9,stroke:#2E7D32
    style E2_B fill:#c8e6c9,stroke:#2E7D32
    style E2_C fill:#fff9c4,stroke:#F57F17,stroke-width:3px
    style E4_A fill:#c8e6c9,stroke:#2E7D32
    style E4_B fill:#c8e6c9,stroke:#2E7D32
    style START fill:#000,stroke:#000
    style STOP fill:#000,stroke:#000
```

**Legend:**
- **Blue nodes** = exons shared across all copies (identical sequence)
- **Green nodes** = exons specific to known copies A and B (in reference)
- **Yellow node (bold border)** = exon from a **novel copy** not in the reference

---

### Read Rescue: Mapping Unmapped Reads to the VG

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '13px'}, 'flowchart': {'nodeSpacing': 15, 'rankSpacing': 50, 'curve': 'basis'}}}%%
flowchart LR

    subgraph VG["Gene Family Variation Graph"]
        direction LR

        E1["<b>Exon 1</b>\n(shared)"]

        E2_A["<b>Exon 2a</b>\nCopy A"]
        E2_B["<b>Exon 2b</b>\nCopy B"]
        E2_C["<b>Exon 2c</b>\nNovel Copy"]

        E3["<b>Exon 3</b>\n(shared)"]

        E4_A["<b>Exon 4a</b>\nCopy A"]
        E4_B["<b>Exon 4b</b>\nCopy B"]

        E5["<b>Exon 5</b>\n(shared)"]

        E1 --> E2_A --> E3
        E1 --> E2_B --> E3
        E1 --> E2_C --> E3
        E3 --> E4_A --> E5
        E3 --> E4_B --> E5
    end

    subgraph READS["Unmapped Reads (rescued)"]
        direction TB
        R1["<b>Read X</b>\nExon1 - Exon2c - Exon3\n(matches novel copy path)"]
        R2["<b>Read Y</b>\nExon2c - Exon3 - Exon4b\n(novel Exon2 + Copy B Exon4)"]
        R3["<b>Read Z</b>\nExon1 - Exon2c - Exon3 - Exon4a\n(novel Exon2 + Copy A Exon4)"]
    end

    R1 -.->|"k-mer\nmatch"| E2_C
    R2 -.->|"k-mer\nmatch"| E2_C
    R3 -.->|"k-mer\nmatch"| E2_C

    style VG fill:#e3f2fd,stroke:#1565C0,stroke-width:3px
    style READS fill:#e8f5e9,stroke:#2E7D32,stroke-width:3px
    style E1 fill:#bbdefb,stroke:#1565C0
    style E3 fill:#bbdefb,stroke:#1565C0
    style E5 fill:#bbdefb,stroke:#1565C0
    style E2_A fill:#c8e6c9,stroke:#2E7D32
    style E2_B fill:#c8e6c9,stroke:#2E7D32
    style E2_C fill:#fff9c4,stroke:#F57F17,stroke-width:3px
    style E4_A fill:#c8e6c9,stroke:#2E7D32
    style E4_B fill:#c8e6c9,stroke:#2E7D32
    style R1 fill:#a5d6a7,stroke:#2E7D32
    style R2 fill:#a5d6a7,stroke:#2E7D32
    style R3 fill:#a5d6a7,stroke:#2E7D32
```

### How the Rescue Works (step by step)

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '14px'}}}%%
flowchart TD
    subgraph STEP1["Step 1: Standard Assembly"]
        S1A["Assemble each locus\nindependently\n(StringTie core)"]
        S1B["Copy A: 2 transcripts\nCopy B: 1 transcript"]
    end

    subgraph STEP2["Step 2: Family Discovery"]
        S2A["Detect that Copy A and\nCopy B share multi-mapped\nreads and k-mer similarity"]
        S2B["Group into\nGene Family F1"]
    end

    subgraph STEP3["Step 3: Build Family VG"]
        S3A["Merge exon graphs\nfrom all copies into\na single variation graph"]
        S3B["Shared exons become\nshared nodes;\ndivergent exons become\nalternate paths"]
    end

    subgraph STEP4["Step 4: Rescue Unmapped Reads"]
        S4A["Scan unmapped reads\nfor k-mer matches\nagainst family VG"]
        S4B["Reads with sufficient\nk-mer overlap are\nassigned to the family"]
        S4C["If enough reads cluster\ntogether with novel sequence:\ncreate NEW BUNDLE\n(novel paralog)"]
    end

    subgraph STEP5["Step 5: Reassemble with Rescued Reads"]
        S5A["Rerun assembly\nwith updated weights\nand novel bundles"]
        S5B["Final output includes\ntranscripts from all copies\nINCLUDING novel paralogs"]
    end

    S1A --> S1B --> S2A --> S2B --> S3A --> S3B --> S4A --> S4B --> S4C --> S5A --> S5B

    style STEP1 fill:#e8f4fd,stroke:#2196F3
    style STEP2 fill:#e8f4fd,stroke:#2196F3
    style STEP3 fill:#e3f2fd,stroke:#1565C0,stroke-width:2px
    style STEP4 fill:#fff3e0,stroke:#FF9800,stroke-width:3px
    style STEP5 fill:#e8f5e9,stroke:#4CAF50,stroke-width:2px
    style S4C fill:#fff9c4,stroke:#F57F17,stroke-width:3px
```

### Why This Matters: Reference Bias in Gene Families

```
    LINEAR REFERENCE                    VARIATION GRAPH
    (1 copy only)                       (all copies)

    =====[E1]==[E2a]==[E3]==[E4a]====   E1 is shared by all copies
                                         |
    Reads from novel copy:               E2 branches into 3 variants:
    "I have exon 2c, never seen          a (Copy A), b (Copy B), c (Novel)
     in this reference"                  |
         |                               E3 is shared
         v                               |
    UNMAPPED / DISCARDED                 E4 branches: a (Copy A), b (Copy B)
    (reference bias!)                    |
                                         Reads from novel copy:
                                         "E2c matches a VG path!"
                                              |
                                              v
                                         RESCUED --> new bundle
                                         --> novel paralog assembled
```

> **Key insight:** Unmapped reads are not necessarily low-quality or
> erroneous. In gene families, they often originate from paralogs that
> are simply absent from the reference assembly. By building a variation
> graph from all known copies and scanning unmapped reads against it,
> Rustle can **rescue** these reads and even **discover novel gene copies**
> that the reference genome is missing, directly reducing reference bias.
