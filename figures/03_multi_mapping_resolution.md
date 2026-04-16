# Multi-Mapping Resolution: EM, MFLP, and Flow

## The Multi-Mapping Problem in Gene Families

When a read maps equally well to multiple gene copies, the aligner reports it
as "multi-mapped" (NH tag > 1). Standard assemblers either discard these reads
or split their weight equally (1/NH), both of which distort transcript abundance
and can merge distinct copies or miss low-expressed copies entirely.

### Why Multi-Mappers Are Real (Not Artifacts)

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '13px'}}}%%
flowchart TB
    subgraph GENOME["Two Gene Copies (99% identical)"]
        direction LR
        subgraph COPY_A["Copy A (chr7:100kb)"]
            A1["Exon 1"] --- A2["Exon 2\n...ACGTGCAT..."] --- A3["Exon 3\n(copy-specific)"]
        end
        subgraph COPY_B["Copy B (chr7:250kb)"]
            B1["Exon 1"] --- B2["Exon 2\n...ACGTGCAT..."] --- B3["Exon 3'\n(different)"]
        end
    end

    subgraph READS["Sequenced Reads"]
        R1["Read 1: Exon1-Exon2\n(maps to BOTH copies)\nNH=2"]
        R2["Read 2: Exon2-Exon3\n(maps to Copy A ONLY)\nNH=1"]
        R3["Read 3: Exon2-Exon3'\n(maps to Copy B ONLY)\nNH=1"]
        R4["Read 4: Exon1-Exon2\n(maps to BOTH copies)\nNH=2"]
    end

    R1 -.->|"0.5"| COPY_A
    R1 -.->|"0.5"| COPY_B
    R2 -->|"1.0"| COPY_A
    R3 -->|"1.0"| COPY_B
    R4 -.->|"0.5"| COPY_A
    R4 -.->|"0.5"| COPY_B

    style COPY_A fill:#c8e6c9,stroke:#2E7D32
    style COPY_B fill:#bbdefb,stroke:#1565C0
    style R1 fill:#fff9c4,stroke:#F57F17
    style R4 fill:#fff9c4,stroke:#F57F17
    style R2 fill:#c8e6c9,stroke:#2E7D32
    style R3 fill:#bbdefb,stroke:#1565C0
```

**The core question:** Read 1 and Read 4 map to both copies. Do they truly come
from both? Or from one copy only? The answer depends on the _expression context_.

---

## Three Resolution Strategies

### Strategy 1: EM (Expectation-Maximization)

Iteratively refines read weights based on junction compatibility scores.
Each round updates the probability that a read belongs to each copy,
converging to a stable assignment.

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '12px'}}}%%
flowchart LR
    subgraph INIT["Initialization"]
        I["All multi-mappers\nget weight = 1/NH\n(uniform)"]
    end

    subgraph ITER["Iteration (repeat until stable)"]
        direction TB
        E["<b>E-step:</b>\nFor each multi-mapped read,\ncompute compatibility score\nwith each copy's junctions"]
        M["<b>M-step:</b>\nRe-estimate copy abundances\nfrom weighted read counts"]
        U["<b>Update:</b>\nSet read weight at each copy\n= abundance(copy) * compat(read,copy)\n/ sum over all copies"]
    end

    subgraph RESULT["After Convergence"]
        direction TB
        RA["Read 1 -> Copy A: 0.7\nRead 1 -> Copy B: 0.3\n(Copy A has more\njunction evidence)"]
        RB["Copy A abundance: 15.4\nCopy B abundance: 8.2"]
    end

    INIT --> ITER --> RESULT
    E --> M --> U
    U -->|"not converged"| E

    style INIT fill:#e8f4fd,stroke:#2196F3
    style ITER fill:#fff3e0,stroke:#FF9800,stroke-width:2px
    style RESULT fill:#e8f5e9,stroke:#4CAF50
```

**Key property:** EM allows **fractional assignment** — a multi-mapper can
contribute weight to multiple copies simultaneously. This is correct when
a read truly could come from either copy (e.g., exon 1-2 reads in a family
where both copies are expressed). EM naturally handles the "read belongs in
several places" scenario by distributing weight proportional to evidence.

---

### Strategy 2: MFLP (Minimum Flow Linear Program)

Single-shot optimization that assigns each multi-mapper to copies by
maximizing total compatibility, subject to the constraint that weights
sum to 1 per read.

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '12px'}}}%%
flowchart TD
    subgraph FORMULATION["LP Formulation"]
        direction TB
        VAR["<b>Variables:</b>\nw[read][copy] in [0, 1]\nweight of read r at copy k"]
        CON["<b>Constraints:</b>\nFor each read r:\n  sum_k w[r][k] = 1\n(weights must sum to 1)"]
        OBJ["<b>Objective:</b>\nMaximize sum_{r,k}\n  w[r][k] * compat(r, k)\n(total compatibility)"]
    end

    subgraph SOLVE["Solver (MicroLP)"]
        LP["Linear program:\nSingle-shot optimization\n(no iteration needed)"]
    end

    subgraph OUTPUT["Result"]
        direction TB
        O1["w[Read1][CopyA] = 1.0\nw[Read1][CopyB] = 0.0\n(MFLP chose best match)"]
        O2["w[Read4][CopyA] = 0.0\nw[Read4][CopyB] = 1.0\n(assigned to other copy)"]
    end

    FORMULATION --> SOLVE --> OUTPUT

    style FORMULATION fill:#e8f4fd,stroke:#2196F3
    style SOLVE fill:#fff3e0,stroke:#FF9800,stroke-width:2px
    style OUTPUT fill:#e8f5e9,stroke:#4CAF50
```

**Key property:** MFLP tends toward **hard assignment** (0 or 1 weights) because
the LP optimum sits at vertices of the feasible polytope. This is appropriate when
you believe each read truly originates from a single copy and want the globally
optimal assignment. However, when a read genuinely comes from a region shared
between expressed copies, MFLP will still assign it to one — potentially
underestimating the weaker copy.

---

### Strategy 3: Flow-Based Redistribution

Two-pass approach: assemble first with uniform weights, then use the
assembled transcripts' coverage to redistribute multi-mapper weights,
then reassemble.

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '12px'}}}%%
flowchart TD
    subgraph PASS1["Pass 1: Initial Assembly"]
        direction TB
        P1A["Assemble all copies\nwith uniform 1/NH weights"]
        P1B["Copy A: 2 transcripts\n(cov = 10.2, 3.1)\nCopy B: 1 transcript\n(cov = 5.7)"]
    end

    subgraph REDISTRIBUTE["Flow Redistribution"]
        direction TB
        FR1["For each multi-mapped read:\n  check which assembled transcripts\n  it supports at each copy"]
        FR2["Redistribute weight proportional\nto transcript coverage:\n  w(Copy A) = 10.2 / (10.2+5.7) = 0.64\n  w(Copy B) = 5.7 / (10.2+5.7) = 0.36"]
    end

    subgraph PASS2["Pass 2: Reweighted Assembly"]
        direction TB
        P2A["Reassemble with\nupdated weights"]
        P2B["Copy A: 2 transcripts\n(cov = 12.1, 3.8)\nCopy B: 1 transcript\n(cov = 4.9)"]
    end

    PASS1 --> REDISTRIBUTE --> PASS2

    style PASS1 fill:#e8f4fd,stroke:#2196F3
    style REDISTRIBUTE fill:#fff3e0,stroke:#FF9800,stroke-width:2px
    style PASS2 fill:#e8f5e9,stroke:#4CAF50
```

**Key property:** Flow uses the **assembler's own output** as evidence for
redistribution. Like EM, it allows fractional weights, but the evidence comes
from assembled transcripts rather than raw junction compatibility. This makes
it robust to cases where junction structure alone is ambiguous but coverage
patterns distinguish the copies.

---

## When Does a Multi-Mapper Really Belong in Several Places?

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '13px'}}}%%
flowchart TD
    subgraph QUESTION["A multi-mapped read maps to Copy A and Copy B"]
        Q["Does this read truly come\nfrom one copy or both?"]
    end

    subgraph SCENARIOS["Three Biological Scenarios"]
        direction TB

        S1["<b>Scenario 1: One copy expressed</b>\nOnly Copy A is transcribed.\nThe read ONLY comes from A.\nCorrect: assign 100% to A."]

        S2["<b>Scenario 2: Both copies expressed</b>\nBoth copies make identical exon 1-2.\nThe read could be from either.\nCorrect: split proportional to\ncopy abundance (e.g., 60/40)."]

        S3["<b>Scenario 3: Novel/divergent copy</b>\nThe read has SNPs matching a\nthird copy not in the reference.\nCorrect: assign to novel copy\n(rescued via VG)."]
    end

    subgraph METHODS["How Each Method Handles This"]
        direction TB
        ME["<b>EM:</b> naturally splits proportional\nto abundance (handles Scenario 2).\nConverges to ~100% for Scenario 1."]
        MM["<b>MFLP:</b> picks the single best copy\n(handles Scenario 1 well).\nMay undercount in Scenario 2."]
        MF["<b>Flow:</b> splits proportional to\nassembled coverage (handles Scenario 2).\nTwo-pass catches Scenario 1."]
        MS["<b>SNP mode (--vg-snp):</b>\nUses per-base variants to distinguish.\nHandles Scenario 3 directly."]
    end

    Q --> SCENARIOS
    SCENARIOS --> METHODS

    style QUESTION fill:#fff3e0,stroke:#FF9800,stroke-width:2px
    style S1 fill:#c8e6c9,stroke:#2E7D32
    style S2 fill:#bbdefb,stroke:#1565C0
    style S3 fill:#fff9c4,stroke:#F57F17
    style ME fill:#fff3e0,stroke:#FF9800
    style MM fill:#fff3e0,stroke:#FF9800
    style MF fill:#fff3e0,stroke:#FF9800
    style MS fill:#fff3e0,stroke:#FF9800
```

### Controlling False Multi-Mapping: How Many Are Real?

The question "how many multi-mappers are real?" is controlled by two mechanisms:

1. **Compatibility scoring:** A read is only assigned weight at a copy where
   its splice junctions are compatible. If Read 1 has junctions A-B-C and
   Copy B only has junctions A-B-D, the compatibility score is low and the
   read gets near-zero weight at Copy B regardless of method.

2. **Abundance feedback:** In EM and Flow, the abundance of each copy acts
   as a prior. If Copy B has very few uniquely-mapped reads, multi-mappers
   get less weight there — the method "learns" that Copy B is weakly expressed
   and avoids over-assigning reads to it.

3. **The fractional truth:** A multi-mapper that maps to two expressed copies
   _does_ belong in both places. EM and Flow correctly split it. This is not
   an error — it reflects the biological reality that the sequencer cannot
   distinguish which copy produced the molecule when the copies are identical
   in that region. The fractional weight (e.g., 0.6/0.4) is the honest answer.

```
    Method    | Fractional? | Handles "belongs in both" | Best for
    ----------|-------------|---------------------------|------------------
    EM        | Yes         | Yes (proportional split)  | General use
    MFLP      | Tends to 0/1| No (picks one)           | Clear-cut cases
    Flow      | Yes         | Yes (coverage-based)      | Complex families
    SNP       | N/A         | Distinguishes copies      | Divergent copies
```
