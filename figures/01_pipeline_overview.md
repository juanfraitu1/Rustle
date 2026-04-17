# Rustle Pipeline Overview

## How StringTie works (and by extension, Rustle)

Rustle is a Rust reimplementation of the StringTie transcript assembly algorithm,
extended with a **Variation Graph (VG) mode** for multi-copy gene family resolution.

### Core StringTie/Rustle Pipeline

```mermaid
flowchart TD
    subgraph INPUT["1. Input"]
        BAM["BAM/SAM\n(aligned reads)"]
        REF["Reference GTF\n(optional guide)"]
        FASTA["Genome FASTA\n(splice consensus)"]
    end

    subgraph BUNDLE["2. Bundle Detection"]
        CLUSTER["Cluster overlapping\nreads into loci"]
        JUNC["Extract & filter\nsplice junctions"]
        CPAS["Detect poly-A\ncut sites (long reads)"]
    end

    subgraph GRAPH["3. Splice Graph Construction"]
        NODES["Segment exonic regions\ninto graph nodes"]
        EDGES["Add collinear +\njunction edges"]
        SRC_SINK["Add source & sink\nnodes"]
        LONGTRIM["Longtrim: split nodes\nat read-boundary peaks\n(long-read mode)"]
    end

    subgraph MAPFLOW["4-5. Read Mapping & Flow"]
        MAP["Map reads onto\ngraph paths (transfrags)"]
        ABSORB["Merge compatible\ntransfrags"]
        MAXFLOW["Edmonds-Karp\nmax-flow decomposition"]
    end

    subgraph EXTRACT["6. Path Extraction"]
        SEED["Seed transcripts from:\n- guide matches\n- junction paths\n- long-read paths"]
        EXTEND["Extend paths\nsource <-> sink"]
        ABUND["Compute per-transcript\nabundance (cov, FPKM, TPM)"]
    end

    subgraph FILTER["7. Transcript Filtering"]
        ISOFRAC["Isofrac: remove\nlow-fraction isoforms"]
        PAIRWISE["Pairwise: remove\ncontained transcripts"]
        DEDUP["Dedup: remove\nintron-chain subsets"]
        RUNOFF["Filter polymerase\nrun-off / run-on"]
    end

    subgraph OUTPUT["8. Output"]
        GTF["Transcript GTF"]
        ABUND_TAB["Gene abundance\ntable"]
        BALL["Ballgown tables\n(optional)"]
    end

    INPUT --> BUNDLE --> GRAPH --> MAPFLOW --> EXTRACT --> FILTER --> OUTPUT

    style INPUT fill:#e8f4fd,stroke:#2196F3
    style BUNDLE fill:#e8f4fd,stroke:#2196F3
    style GRAPH fill:#e8f4fd,stroke:#2196F3
    style MAPFLOW fill:#e8f4fd,stroke:#2196F3
    style EXTRACT fill:#e8f4fd,stroke:#2196F3
    style FILTER fill:#e8f4fd,stroke:#2196F3
    style OUTPUT fill:#e8f4fd,stroke:#2196F3
```

### Rustle's VG Extension: Variation Graph Mode

The standard StringTie pipeline (above) treats each locus independently.
Rustle adds a **Variation Graph (VG) layer** that links related loci
(paralogs, tandem duplicates) and resolves multi-mapping reads across them.

```mermaid
flowchart TD
    subgraph STRINGTIE_CORE["StringTie Core Pipeline (per-locus)"]
        direction LR
        A["Bundles"] --> B["Splice\nGraphs"] --> C["Flow\nDecomposition"] --> D["Transcripts"]
    end

    subgraph VG_LAYER["Rustle VG Extension  --vg"]
        direction TB

        DISCOVER["VG Discovery\n- Index multi-mapped reads\n  across bundles\n- K-mer Jaccard similarity\n  for tandem duplicates\n- Union-Find to group\n  bundles into families"]

        subgraph SOLVERS["Multi-mapping Resolution"]
            direction LR
            EM["EM\nIterative\nreweighting"]
            MFLP["MFLP\nLinear\nprogramming"]
            FLOW_S["Flow\nTwo-pass\nredistribution"]
        end

        REWEIGHT["Update read weights\nper gene-family copy"]
    end

    subgraph APPLICATIONS["Applications Enabled by VG"]
        direction TB

        APP1["<b>1. New Copies</b>\n--vg-discover-novel\nMap unmapped reads to graphs\nto discover paralogs absent\nfrom the reference genome\n(reduces reference bias)"]

        APP2["<b>2. Multi-mapping Resolution</b>\nResolve ambiguous read assignments\nacross gene family copies using:\n  MFLP (LP optimization)\n  EM (expectation-maximization)\n  Flow (coverage redistribution)"]

        APP3["<b>3. SNP-based Assignment</b>\n--vg-snp\nUse sequence variants (MD tag)\nto distinguish gene copies\nby their unique SNP profiles"]

        APP4["<b>4. Phased Assembly</b>\n--vg-phase\nUse haplotype phase tags (HP)\nto split bundles by haplotype\nand co-assemble related copies\nwith haplotype awareness"]
    end

    STRINGTIE_CORE --> VG_LAYER
    DISCOVER --> SOLVERS --> REWEIGHT
    VG_LAYER --> STRINGTIE_CORE
    VG_LAYER -.- APPLICATIONS

    style STRINGTIE_CORE fill:#e8f4fd,stroke:#2196F3
    style VG_LAYER fill:#fff3e0,stroke:#FF9800,stroke-width:3px
    style APPLICATIONS fill:#e8f5e9,stroke:#4CAF50
    style APP1 fill:#e8f5e9,stroke:#4CAF50,text-align:left
    style APP2 fill:#e8f5e9,stroke:#4CAF50,text-align:left
    style APP3 fill:#e8f5e9,stroke:#4CAF50,text-align:left
    style APP4 fill:#e8f5e9,stroke:#4CAF50,text-align:left
    style EM fill:#fff3e0,stroke:#FF9800
    style MFLP fill:#fff3e0,stroke:#FF9800
    style FLOW_S fill:#fff3e0,stroke:#FF9800
```

### How the VG layer interacts with the core pipeline

```mermaid
flowchart LR
    subgraph PASS1["Pass 1: Independent Assembly"]
        B1["Bundle A\n(copy 1)"]
        B2["Bundle B\n(copy 2)"]
        B3["Bundle C\n(copy 3)"]
    end

    subgraph VG["VG: Link & Resolve"]
        LINK["Discover family\n(shared reads +\nk-mer similarity)"]
        SOLVE["Resolve multi-maps\n(EM / MFLP / Flow)"]
    end

    subgraph PASS2["Pass 2: Reweighted Assembly"]
        R1["Bundle A\n(adjusted weights)"]
        R2["Bundle B\n(adjusted weights)"]
        R3["Bundle C\n(adjusted weights)"]
    end

    subgraph NOVEL["Novel Copy Discovery"]
        UNMAPPED["Unmapped reads\nwith k-mer matches"]
        NEWBUNDLE["New bundle\n(novel paralog)"]
    end

    B1 --> LINK
    B2 --> LINK
    B3 --> LINK
    LINK --> SOLVE
    SOLVE --> R1
    SOLVE --> R2
    SOLVE --> R3

    UNMAPPED --> NEWBUNDLE
    NEWBUNDLE --> VG

    style PASS1 fill:#e8f4fd,stroke:#2196F3
    style VG fill:#fff3e0,stroke:#FF9800,stroke-width:3px
    style PASS2 fill:#c8e6c9,stroke:#4CAF50
    style NOVEL fill:#fce4ec,stroke:#E91E63
```

### Key Insight

> **StringTie** assembles each locus in isolation.
> **Rustle** adds a graph-aware layer that recognizes when multiple loci
> belong to the same gene family and uses variation-graph principles
> (multi-mapping resolution, SNP discrimination, haplotype phasing,
> novel copy discovery) to produce more accurate, less reference-biased
> transcript assemblies for duplicated gene families.
