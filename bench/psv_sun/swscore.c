/* swscore — best LOCAL alignment score (Smith-Waterman, Gotoh affine gaps) for each pair of sequences.
 * Diagnostic for docs/PSV_SUN_PILEUP_2026-10-12.md (are the unaligned stretches between copies homologous?).
 * BLASTN scoring: match +2, mismatch -3, gap open 5, extend 2 (a gap of length L costs 5 + 2L). Score only, O(n) memory.
 *
 *   gcc -O2 -o swscore swscore.c && ./swscore < pairs.tsv > scores.tsv   (pairs.tsv: id <TAB> seqA <TAB> seqB)
 * Output: id <TAB> score
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

enum { MATCH = 2, MISMATCH = -3, GAP_OPEN = 5, GAP_EXT = 2 };

static int sw(const char *a, int n, const char *b, int m, int *H, int *E) {
    /* H[j]: best score ending at (i, j); E[j]: best ending at (i, j) in a gap in a (vertical); F: horizontal */
    int best = 0;
    for (int j = 0; j <= m; j++) { H[j] = 0; E[j] = -1000000; }
    for (int i = 1; i <= n; i++) {
        int diag = 0, F = -1000000;
        H[0] = 0;
        for (int j = 1; j <= m; j++) {
            int e = E[j] - GAP_EXT, eo = H[j] - GAP_OPEN - GAP_EXT;
            E[j] = e > eo ? e : eo;
            int f = F - GAP_EXT, fo = H[j - 1] - GAP_OPEN - GAP_EXT;
            F = f > fo ? f : fo;
            int s = diag + (toupper((unsigned char)a[i - 1]) == toupper((unsigned char)b[j - 1]) ? MATCH : MISMATCH);
            int h = s;
            if (E[j] > h) h = E[j];
            if (F > h) h = F;
            if (h < 0) h = 0;
            diag = H[j];
            H[j] = h;
            if (h > best) best = h;
        }
    }
    return best;
}

int main(void) {
    size_t cap = 1 << 16;
    char *line = malloc(cap);
    int hcap = 0, *H = NULL, *E = NULL;
    ssize_t len;
    while ((len = getline(&line, &cap, stdin)) > 0) {
        line[strcspn(line, "\r\n")] = 0;
        char *id = strtok(line, "\t"), *a = strtok(NULL, "\t"), *b = strtok(NULL, "\t");
        if (!id || !a || !b) continue;
        int n = (int)strlen(a), m = (int)strlen(b);
        if (m + 1 > hcap) {
            hcap = m + 1;
            H = realloc(H, hcap * sizeof *H);
            E = realloc(E, hcap * sizeof *E);
        }
        printf("%s\t%d\n", id, sw(a, n, b, m, H, E));
    }
    return 0;
}
