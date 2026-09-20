BEGIN {
    FS = OFS = "\t"
}

FNR == NR {
    if (FNR == 1) {
        next
    }
    n++
    bid[n] = $1
    bstart[n] = $4
    bend[n] = $5
    raw[n] = $9
    corr[n] = $10
    red[n] = $13
    frac[n] = $14
    next
}

FNR == 1 {
    print $0, "max_raw_transfrags", "max_corrected_transfrags", "max_transfrag_reduction", "max_transfrag_reduction_frac", "max_transfrag_bundle_id"
    next
}

{
    bestRed = -1
    bestFrac = -1
    bestRaw = ""
    bestCorr = ""
    bestId = ""

    targetStart = $3
    targetEnd = $4

    for (i = 1; i <= n; i++) {
        if (bstart[i] <= targetEnd && bend[i] >= targetStart) {
            if (red[i] > bestRed || (red[i] == bestRed && frac[i] > bestFrac)) {
                bestRed = red[i]
                bestFrac = frac[i]
                bestRaw = raw[i]
                bestCorr = corr[i]
                bestId = bid[i]
            }
        }
    }

    print $0, bestRaw, bestCorr, bestRed, bestFrac, bestId
}
