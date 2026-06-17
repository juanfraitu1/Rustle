export const meta = {
  name: 'hidden-collapse-verify',
  description: 'Adversarially verify the hidden-collapse scan: pressure-test the method against FP modes + verify every TIER A/B/C hit from BAM evidence, then synthesize an honest headroom.',
  phases: [
    { title: 'Methodology', detail: 'judge panel: each agent refutes one false-positive mode' },
    { title: 'Verify', detail: 'per-hit BAM evidence verification (batched, read from TSV)' },
    { title: 'Synthesize', detail: 'honest count + completeness critic' },
  ],
}

// args = { n_hits: <int>, dist: "<verdict totals + n_coseg histogram text>" }
const N = (args && args.n_hits) || 0
const dist = (args && args.dist) || '(no distribution provided)'
const SCAN = '/mnt/c/Users/jfris/Desktop/Rustle/bench/hidden_collapse_scan.py'
const EVID = '/mnt/c/Users/jfris/Desktop/Rustle/bench/hidden_collapse_evidence.py'
const HITS = '/mnt/c/Users/jfris/Desktop/Rustle/bench/hidden_collapse_hits.tsv'
// hits TSV columns: tier chrom start end n_isoforms n_fsm n_psv n_coseg n_copies frac_sec frac_mq0 n_aln minor_frac

const COMMON = `
Context: task (c) tests whether read-coherence recall isoforms at ANNOTATED single-copy loci secretly
sit on HIDDEN collapsed/cross-mapped paralog copies with copy-discriminating PSV signal (which the
annotation-based geometric headroom probe could not see). A detector calls PSVs de novo from the GGO
HiFi BAM pileup and tests CO-SEGREGATION (the validated copy_split identifiability signal). The honest
headroom = loci that are REAL hidden collapsed paralogs, NOT false-positive modes.

Detector source: ${SCAN}
Per-locus evidence tool (run it yourself via Bash):  python3 ${EVID} --region <chrom:start-end>
It prints: candidate PSV columns; frac_mq0 (this locus's PRIMARY reads multimap => a copy exists
elsewhere); frac_sec (reads from elsewhere spill here); base-change spectrum (A>G+T>C fraction >0.8 =>
RNA-editing; mixed 12 substitution types with Ti:Tv~2 => paralog divergence); the co-segregating copy
haplotypes + read counts + minor-copy fraction; primary MAPQ histogram; reference low-complexity frac.

Tier meaning: A = >=3 allele groups (cannot be one diploid het locus); B = exactly 2 groups but >=8
linked PSV columns with low multimapping (far above diploid-het exonic expectation ~3/transcript);
C = multimapping-driven (could be a segdup/repeat/pseudogene, softest). Genome-wide distribution:
${dist}
`

// ---------------- Phase 1: methodology judge panel ----------------
phase('Methodology')
const FP_MODES = [
  { key: 'het', desc: 'DIPLOID HETEROZYGOSITY: a highly polymorphic gene whose two haplotypes carry many linked het SNPs co-segregating into 2 read groups, mimicking a 2-copy collapse (the main threat to TIER B).' },
  { key: 'rna_editing', desc: 'RNA EDITING: clustered A-to-I (read A>G on +, T>C on -) sites that co-segregate by molecule, mimicking PSVs.' },
  { key: 'lowcomplexity', desc: 'LOW-COMPLEXITY / STR / homopolymer alignment artifacts: aligners pile correlated mismatches at repeats, faking a linked block.' },
  { key: 'segdup', desc: 'SEGMENTAL DUPLICATION / pseudogene cross-mapping that is real but NOT an isoform-copy-resolvable transcript pair (the main threat to TIER C).' },
  { key: 'paralog_resolved', desc: 'PARALOG ALREADY RESOLVED by minimap2: copies map to distinct loci; what remains here is residual cross-map noise, not a genuine collapse at THIS locus.' },
]
const METH_SCHEMA = {
  type: 'object',
  required: ['fp_mode', 'can_inflate_AB', 'severity', 'discriminating_check', 'verdict'],
  properties: {
    fp_mode: { type: 'string' },
    can_inflate_AB: { type: 'boolean' },
    severity: { type: 'string', enum: ['none', 'low', 'medium', 'high'] },
    discriminating_check: { type: 'string' },
    recommended_change: { type: 'string' },
    verdict: { type: 'string' },
  },
}
const methodology = await parallel(FP_MODES.map((fm) => () =>
  agent(
    `${COMMON}\nYou are an adversarial reviewer. Argue HARD that the FP mode below inflates the COLLAPSED_LIKE count, then assess honestly whether the detector's guards (de-novo PSV freq floor f>=0.2 & >=3 reads; >=K=3 co-segregation; >=3-group test for TIER A; >=8 linked-column threshold for TIER B; MAPQ-0/secondary signature) control it. Read the detector source; run the evidence tool on 2-3 example hits if useful.\n\nFP MODE: ${fm.desc}`,
    { label: `method:${fm.key}`, phase: 'Methodology', schema: METH_SCHEMA }
  ).then((r) => ({ ...r, fp_mode: fm.key }))
)).then((a) => a.filter(Boolean))
log(`methodology: ${methodology.length} FP modes assessed; ${methodology.filter((m) => m.can_inflate_AB).length} can inflate TIER A/B`)

// ---------------- Phase 2: per-hit verification (batched, read TSV by index) ----------------
phase('Verify')
const VERIFY_SCHEMA = {
  type: 'object',
  required: ['results'],
  properties: {
    results: {
      type: 'array',
      items: {
        type: 'object',
        required: ['region', 'tier', 'verdict', 'n_isoforms', 'n_fsm', 'reasoning'],
        properties: {
          region: { type: 'string' },
          tier: { type: 'string' },
          n_isoforms: { type: 'integer' },
          n_fsm: { type: 'integer' },
          verdict: { type: 'string', enum: ['REAL_HIDDEN_COLLAPSE', 'HET', 'RNA_EDITING', 'REPEAT_ARTIFACT', 'TOO_THIN', 'UNCERTAIN'] },
          reasoning: { type: 'string' },
        },
      },
    },
  },
}
const BATCH = Math.max(8, Math.ceil(N / 24))
const ranges = []
for (let lo = 1; lo <= N; lo += BATCH) ranges.push([lo, Math.min(lo + BATCH - 1, N)])
log(`verifying ${N} hits in ${ranges.length} batches of ~${BATCH}`)

const verified = await parallel(ranges.map(([lo, hi], bi) => () =>
  agent(
    `${COMMON}\nVERIFY a batch of hits. First get your rows (1-based data rows ${lo}..${hi}) from the hits TSV:\n` +
    `  awk -F'\\t' 'NR>1 && (NR-1)>=${lo} && (NR-1)<=${hi}' ${HITS}\n` +
    `Columns: tier chrom start end n_isoforms n_fsm n_psv n_coseg n_copies frac_sec frac_mq0 n_aln minor_frac.\n` +
    `For EACH row, run \`python3 ${EVID} --region <chrom:start-end>\`, read the evidence, and classify verdict:\n` +
    `- REAL_HIDDEN_COLLAPSE: >=2 clean balanced co-segregating copies, MIXED substitution spectrum (NOT A>G/T>C-dominated), low reference low-complexity, AND (>=3 groups OR many linked columns OR genuine multimapping: low primary MAPQ / high frac_mq0/frac_sec). A real paralog copy the annotation misses.\n` +
    `- HET: 2 groups, modest linked SNPs, HIGH primary MAPQ, no multimapping (diploid heterozygosity). A highly polymorphic gene with 8-15 clean biallelic linked SNPs and high MAPQ is HET, not collapse.\n` +
    `- RNA_EDITING: A>G+T>C fraction > 0.8.\n` +
    `- REPEAT_ARTIFACT: high reference low-complexity/homopolymer at the block.\n` +
    `- TOO_THIN: too few reads/columns to judge.\n` +
    `- UNCERTAIN: genuinely ambiguous.\n` +
    `Be skeptical: default to the FP explanation unless REAL_HIDDEN_COLLAPSE evidence is clear. Carry tier, n_isoforms, n_fsm from the TSV row into each result. Set region as "chrom:start-end".`,
    { label: `verify:b${bi}[${lo}-${hi}]`, phase: 'Verify', schema: VERIFY_SCHEMA }
  )
)).then((a) => a.filter(Boolean))

const results = verified.flatMap((v) => (v && v.results) || [])
const isReal = (r) => r.verdict === 'REAL_HIDDEN_COLLAPSE'
const sum = (arr, f) => arr.reduce((s, x) => s + (Number(f(x)) || 0), 0)
const ab = results.filter((r) => r.tier === 'A_ge3groups' || r.tier === 'B_dense2copy')
const realAB = ab.filter(isReal)
const realAll = results.filter(isReal)
const breakdown = {}
for (const r of results) breakdown[r.verdict] = (breakdown[r.verdict] || 0) + 1
log(`verification: ${realAll.length}/${results.length} REAL_HIDDEN_COLLAPSE; ` +
    `confound-controlled (A+B) real = ${realAB.length} loci / ${sum(realAB, (x) => x.n_isoforms)} isoforms / ${sum(realAB, (x) => x.n_fsm)} FSM; ` +
    `breakdown ${JSON.stringify(breakdown)}`)

// ---------------- Phase 3: synthesis + completeness critic ----------------
phase('Synthesize')
const SYNTH_SCHEMA = {
  type: 'object',
  required: ['honest_headroom_loci', 'honest_headroom_isoforms', 'go_no_go', 'summary', 'completeness_gaps'],
  properties: {
    honest_headroom_loci: { type: 'integer' },
    honest_headroom_isoforms: { type: 'integer' },
    honest_headroom_fsm: { type: 'integer' },
    go_no_go: { type: 'string', enum: ['GO', 'NO_GO', 'WEAK_GO'] },
    summary: { type: 'string' },
    completeness_gaps: { type: 'array', items: { type: 'string' } },
    caveats: { type: 'array', items: { type: 'string' } },
  },
}
const synthesis = await agent(
  `${COMMON}\nSynthesize an honest go/no-go on the undercount caveat.\n\n` +
  `METHODOLOGY (FP modes): ${JSON.stringify(methodology)}\n\n` +
  `VERIFICATION: ${results.length} hits judged. Verdict breakdown ${JSON.stringify(breakdown)}.\n` +
  `Confound-controlled (TIER A+B) confirmed REAL_HIDDEN_COLLAPSE = ${realAB.length} loci / ` +
  `${sum(realAB, (x) => x.n_isoforms)} recall isoforms / ${sum(realAB, (x) => x.n_fsm)} FSM.\n` +
  `All-tier confirmed REAL = ${realAll.length} loci / ${sum(realAll, (x) => x.n_isoforms)} isoforms / ${sum(realAll, (x) => x.n_fsm)} FSM.\n\n` +
  `Recall: the geometric (annotation-based) headroom probe found 0 recall isoforms at confirmed-real ` +
  `collapsed paralogs. This task asks whether the BAM reveals HIDDEN ones. Decide whether threading PSVs ` +
  `through the molecule graph would yield REAL copy-resolution on recall loci (GO/WEAK_GO/NO_GO), give the ` +
  `honest confound-controlled headroom numbers, and list COMPLETENESS gaps the detector structurally cannot ` +
  `see (e.g. identical-sequence copies => no PSVs => invisible; copies whose reads never reach this locus).`,
  { label: 'synthesize', phase: 'Synthesize', schema: SYNTH_SCHEMA }
)

return {
  methodology,
  verification_breakdown: breakdown,
  confirmed_real_AB: { loci: realAB.length, isoforms: sum(realAB, (x) => x.n_isoforms), fsm: sum(realAB, (x) => x.n_fsm) },
  confirmed_real_all: { loci: realAll.length, isoforms: sum(realAll, (x) => x.n_isoforms), fsm: sum(realAll, (x) => x.n_fsm) },
  confirmed_regions: realAB.map((r) => ({ region: r.region, tier: r.tier, reasoning: r.reasoning })),
  synthesis,
}
