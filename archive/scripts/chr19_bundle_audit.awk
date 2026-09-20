BEGIN {
    OFS = "\t"
}

FNR == 1 && FILENAME ~ /GGO_19_familytrace_full\.bundles\.tsv$/ {
    next
}

FNR == 1 && FILENAME ~ /GGO_19_familytrace_full\.junctions\.tsv$/ {
    next
}

FILENAME ~ /GGO_19_familytrace_full\.bundles\.tsv$/ {
    id = $1
    bchrom[id] = $2
    bstrand[id] = $3
    bstart[id] = $4 + 0
    bend[id] = $5 + 0
    breads[id] = $6 + 0
    brawfam[id] = $17 + 0
    bcorrfam[id] = $18 + 0
    bfamchg[id] = $19 + 0
    bundle_ids[++nb] = id
    next
}

FILENAME ~ /GGO_19_familytrace_full\.junctions\.tsv$/ {
    id = $1
    cur_bad[id] += ($13 + 0)
    par_bad[id] += ($14 + 0)
    par_only_bad[id] += (($14 + 0) == 1 && ($13 + 0) == 0) ? 1 : 0
    if ($20 != "none") {
        cur_redir[id]++
    }
    if ($21 != "none") {
        par_redir[id]++
    }
    next
}

FILENAME ~ /GGO_19\.log$/ {
    if (match($0, /^\[[^]]+\]>bundle ([^:]+):([0-9]+)-([0-9]+) \[([0-9]+) alignments \(([0-9]+) distinct\), ([0-9]+) junctions, ([0-9]+) guides\] begins processing\.\.\./, m)) {
        log_n++
        lchrom[log_n] = m[1]
        lstart[log_n] = m[2] + 0
        lend[log_n] = m[3] + 0
        lalign[log_n] = m[4] + 0
        ldistinct[log_n] = m[5] + 0
        ljunc[log_n] = m[6] + 0
        lguides[log_n] = m[7] + 0
        next
    }

    if (match($0, /^\[[^]]+\]\^bundle ([^:]+):([0-9]+)-([0-9]+) done \(([0-9]+) processed potential transcripts\)\./, m)) {
        ldone_n++
        ltx[ldone_n] = m[4] + 0
        next
    }

    next
}

END {
    print "log_bundle_id", "chrom", "start", "end", "stringtie_alignments", "stringtie_distinct", "stringtie_junctions", "stringtie_guides", "stringtie_transcripts", "our_trace_bundles", "our_trace_reads", "our_current_bad_junctions", "our_parity_bad_junctions", "our_parity_only_bad_junctions", "our_current_redirects", "our_parity_redirects", "our_raw_families", "our_corrected_families", "our_family_key_changes", "our_bundle_ids", "our_bundle_strands"

    for (i = 1; i <= log_n; i++) {
        trace_bundles = 0
        trace_reads = 0
        current_bad = 0
        parity_bad = 0
        parity_only = 0
        cur_r = 0
        par_r = 0
        rawf = 0
        corrf = 0
        famchg = 0
        ids = ""
        strands = ""

        for (k = 1; k <= nb; k++) {
            id = bundle_ids[k]
            if (bchrom[id] != lchrom[i]) {
                continue
            }
            if (bstart[id] <= lend[i] && bend[id] >= lstart[i]) {
                trace_bundles++
                trace_reads += breads[id]
                current_bad += cur_bad[id]
                parity_bad += par_bad[id]
                parity_only += par_only_bad[id]
                cur_r += cur_redir[id]
                par_r += par_redir[id]
                rawf += brawfam[id]
                corrf += bcorrfam[id]
                famchg += bfamchg[id]
                ids = ids (ids ? "," : "") id
                strands = strands (strands ? "," : "") bstrand[id]
            }
        }

        print i, lchrom[i], lstart[i], lend[i], lalign[i], ldistinct[i], ljunc[i], lguides[i], ltx[i], trace_bundles, trace_reads, current_bad, parity_bad, parity_only, cur_r, par_r, rawf, corrf, famchg, ids, strands
    }
}
