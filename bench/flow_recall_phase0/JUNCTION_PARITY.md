# Junction-acceptance parity rustle vs StringTie

**Date:** 2026-06-09. Both tools, junction_accept parity logs, NC_073227.2 whole chrom
(any-accepted semantics: a junction is "good" if accepted in any log entry).

## Result: rustle's good-junctions are a STRICT SUBSET of ST's
- rustle accepted **10,008** · ST accepted **24,952** · both **10,007** · rustle-only **1** · ST-only **14,945**.
- ST's 14,945 extra junctions: **86% have <=1 supporting read** (median support = 1). ST keeps
  essentially EVERY junction it sees, including single-read ones.
- Every junction rustle rejects (5,535), ST accepts. Of ST's extra: ~5,535 rustle saw+rejected
  (its filter) + ~9,410 rustle never saw (stricter read->junction collection).

## Interpretation — deliberate filter-STAGE difference, not a bug
- ST: permissive junctions (keep all) -> strict DOWNSTREAM path/transcript filtering.
- rustle: strict junctions early -> different downstream filtering.
Making junctions "exactly the same as ST" = accepting ~15k single-read junctions/chrom WITHOUT
ST's exact downstream filter = precision catastrophe. Matches the st_shadow_mode finding (literal
ST-faithful junction acceptance REGRESSES). rustle's stricter junction set is CORRECT.

## Key: junction acceptance is NOT where the recoverable signal is
The near-misses that matter (altsplice, etc.) are NOT junction-acceptance differences — at the
traced altsplice case BOTH tools accepted the same donors (40-read & 28-read); the divergence was
PATH SELECTION (parse_trflong) on junctions both already accept. So junction parity is a dead end;
the real divergence is one stage later.
