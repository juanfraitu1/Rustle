// kmerhit: count occurrences of a small set of canonical 21-mers in a FASTA stream.
// usage: kmerhit qset.u64 out.i32 < reads.fa      qset.u64 = sorted unique canonical 21-mers (uint64 little endian, 2 bits per base A0 C1 G2 T3)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define K 21
static const uint64_t KM = (1ULL << (2 * K)) - 1;
static uint64_t *tab; static int32_t *val; static size_t mask;
static inline size_t h(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x & mask; }
int main(int argc, char **argv) {
  if (argc < 3) return 2;
  FILE *f = fopen(argv[1], "rb"); fseek(f, 0, SEEK_END); long nb = ftell(f); fseek(f, 0, SEEK_SET); long nq = nb / 8;
  uint64_t *q = malloc(nb); if (fread(q, 8, nq, f) != (size_t)nq) return 3; fclose(f);
  size_t cap = 1; while (cap < (size_t)nq * 4) cap <<= 1; mask = cap - 1;
  tab = malloc(cap * 8); val = malloc(cap * 4); memset(tab, 0xff, cap * 8);
  for (long i = 0; i < nq; i++) { size_t p = h(q[i]); while (tab[p] != ~0ULL) p = (p + 1) & mask; tab[p] = q[i]; val[p] = (int32_t)i; }
  int32_t *cnt = calloc(nq, 4);
  static int8_t code[256]; memset(code, -1, 256); code['A'] = code['a'] = 0; code['C'] = code['c'] = 1; code['G'] = code['g'] = 2; code['T'] = code['t'] = 3;
  char *line = NULL; size_t cl = 0; ssize_t n; long recs = 0, bases = 0;
  while ((n = getline(&line, &cl, stdin)) > 0) {
    if (line[0] == '>') { recs++; continue; }
    uint64_t fw = 0, rv = 0; int len = 0;
    for (ssize_t i = 0; i < n; i++) {
      int c = code[(unsigned char)line[i]]; if (c < 0) { len = 0; fw = rv = 0; continue; }
      fw = ((fw << 2) | c) & KM; rv = (rv >> 2) | ((uint64_t)(3 - c) << (2 * (K - 1))); bases++;
      if (++len >= K) { uint64_t x = fw < rv ? fw : rv; size_t p = h(x); while (tab[p] != ~0ULL) { if (tab[p] == x) { cnt[val[p]]++; break; } p = (p + 1) & mask; } }
    }
  }
  f = fopen(argv[2], "wb"); fwrite(cnt, 4, nq, f); fclose(f);
  fprintf(stderr, "records=%ld bases=%ld\n", recs, bases); return 0;
}
