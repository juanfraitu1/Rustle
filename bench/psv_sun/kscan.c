/* kscan — count genome occurrences of a list of candidate k-mers (k <= 31), both strands (canonical).
 * docs/PREREG_psv_sun_pileup_2026-10-12.md M1. Memory ~ 2 x 16 B per candidate; the genome is streamed.
 *
 *   gcc -O2 -o kscan kscan.c && ./kscan K kmers.txt genome.fa > counts.tsv   (kmers.txt: one ACGT k-mer per line)
 * Output: kmer <TAB> count (canonical form as given, count saturates at 65535). Non-ACGT bases break k-mers.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int K;
static uint64_t MASK, *keys;
static uint16_t *cnt;
static size_t cap;
static const uint64_t EMPTY = ~0ULL;

static int code(int c) {
    switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
    }
}
static uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; return x; }
static size_t slot(uint64_t k) {
    size_t i = mix(k) & (cap - 1);
    while (keys[i] != EMPTY && keys[i] != k) i = (i + 1) & (cap - 1);
    return i;
}
static uint64_t enc(const char *s, uint64_t *rc) {
    uint64_t f = 0, r = 0;
    for (int i = 0; i < K; i++) {
        int c = code(s[i]);
        if (c < 0) return EMPTY;
        f = (f << 2) | (uint64_t)c;
        r |= (uint64_t)(3 - c) << (2 * i);
    }
    *rc = r;
    return f;
}

int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: kscan K kmers.txt genome.fa\n"); return 2; }
    K = atoi(argv[1]);
    if (K < 1 || K > 31) { fprintf(stderr, "K must be 1..31\n"); return 2; }
    MASK = (K == 32) ? ~0ULL : ((1ULL << (2 * K)) - 1);
    FILE *kf = fopen(argv[2], "r");
    if (!kf) { perror(argv[2]); return 1; }
    char line[256];
    size_t n = 0;
    while (fgets(line, sizeof line, kf)) n++;
    rewind(kf);
    cap = 1024;
    while (cap < 2 * n + 16) cap <<= 1;
    keys = malloc(cap * sizeof *keys);
    cnt = calloc(cap, sizeof *cnt);
    memset(keys, 0xff, cap * sizeof *keys);
    char **names = calloc(n, sizeof *names);
    uint64_t *canon = calloc(n, sizeof *canon);
    size_t m = 0;
    while (fgets(line, sizeof line, kf)) {
        line[strcspn(line, "\r\n")] = 0;
        if ((int)strlen(line) != K) continue;
        uint64_t rc, f = enc(line, &rc);
        if (f == EMPTY) continue;
        uint64_t c = f < rc ? f : rc;
        keys[slot(c)] = c;
        names[m] = strdup(line);
        canon[m++] = c;
    }
    fclose(kf);
    FILE *g = fopen(argv[3], "r");
    if (!g) { perror(argv[3]); return 1; }
    static char buf[1 << 20];
    uint64_t f = 0, r = 0;
    int len = 0, header = 0;
    size_t got;
    while ((got = fread(buf, 1, sizeof buf, g)) > 0) {
        for (size_t i = 0; i < got; i++) {
            char ch = buf[i];
            if (ch == '>') { header = 1; len = 0; continue; }
            if (ch == '\n') { header = 0; continue; }
            if (header || ch == '\r') continue;
            int c = code(ch);
            if (c < 0) { len = 0; continue; }
            f = ((f << 2) | (uint64_t)c) & MASK;
            r = (r >> 2) | ((uint64_t)(3 - c) << (2 * (K - 1)));
            if (++len >= K) {
                uint64_t cc = f < r ? f : r;
                size_t s = slot(cc);
                if (keys[s] == cc && cnt[s] < 65535) cnt[s]++;
            }
        }
    }
    fclose(g);
    for (size_t i = 0; i < m; i++) printf("%s\t%u\n", names[i], cnt[slot(canon[i])]);
    return 0;
}
