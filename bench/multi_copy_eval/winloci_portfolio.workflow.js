export const meta = {
  name: 'vg-winloci-portfolio',
  description: 'Per-locus: does rustle+VG beat StringTie via multimapping reads (vs RefSeq truth)? Evaluate each paralog candidate, then adversarially verify each apparent win is genuinely read-backed (not a DAZ3-style phantom).',
  phases: [
    { title: 'Eval', detail: 'StringTie vs rustle-baseline vs rustle-VG(win stack), scored against RefSeq per locus' },
    { title: 'Verify', detail: 'strict-NM read-backing of each apparent win — the ZFY-real vs DAZ3-phantom discriminant' },
  ],
}

const A = (typeof args === 'string') ? JSON.parse(args) : (args || {})
const cands = A.candidates || []
log(`parsed ${cands.length} candidates (args was ${typeof args})`)
const VGFLAGS = A.vgflags || 'RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1'
const GI = A.gi || '/mnt/c/Users/jfris/Desktop/Rustle/bench/paralog_secondary_scan/scan_out/gene_introns.tsv'
const FLANK = A.flank || 15000
const SH = '/mnt/c/Users/jfris/Desktop/Rustle/bench/multi_copy_eval'

const EVAL_SCHEMA = {
  type: 'object',
  properties: {
    gene: { type: 'string' },
    region: { type: 'string' },
    strand: { type: 'string' },
    n_primary_reads: { type: 'integer' },
    n_secondary: { type: 'integer' },
    st_tx: { type: 'integer' }, base_tx: { type: 'integer' }, vg_tx: { type: 'integer' },
    st_eq: { type: 'array', items: { type: 'string' } },
    base_eq: { type: 'array', items: { type: 'string' } },
    vg_eq: { type: 'array', items: { type: 'string' } },
    vg_gain_vs_st: { type: 'array', items: { type: 'string' } },
    vg_gain_vs_base: { type: 'array', items: { type: 'string' } },
    vg_loss_vs_st: { type: 'array', items: { type: 'string' } },
    classification: { type: 'string', enum: ['win_vs_st', 'win_vs_base_only', 'tie', 'regress', 'error'] },
  },
  required: ['gene', 'classification', 'vg_gain_vs_st', 'vg_eq', 'st_eq', 'base_eq'],
}

const VERDICT_SCHEMA = {
  type: 'object',
  properties: {
    gene: { type: 'string' },
    n_chain_reads: { type: 'integer' },
    n_strict_owner: { type: 'integer' },
    n_pure_owner: { type: 'integer' },
    n_tied: { type: 'integer' },
    n_sibling_favored: { type: 'integer' },
    owner_frac: { type: 'number' },
    verdict: { type: 'string', enum: ['real', 'phantom', 'ambiguous', 'weak', 'error'] },
    reason: { type: 'string' },
  },
  required: ['gene', 'verdict'],
}

function evalPrompt(c) {
  const sib = (c.dom_sibling && String(c.dom_sibling).indexOf(':') === -1) ? c.dom_sibling : ''
  return `Run this EXACT command (it extracts the locus + sibling from the genome BAM, runs StringTie + rustle-baseline + rustle-VG, scores each vs RefSeq, and prints ONE JSON object on stdout):

GI='${GI}' FLANK=${FLANK} VGFLAGS='${VGFLAGS}' bash ${SH}/winloci_eval.sh '${c.gene_id}' '${c.chrom}' ${c.start} ${c.end} '${c.strand}' '${sib}'

Return the JSON object it prints, verbatim, via StructuredOutput. Do NOT invent or alter field values. If the script errors or prints no JSON, return {gene:'${c.gene_id}', classification:'error', vg_gain_vs_st:[], vg_eq:[], st_eq:[], base_eq:[]}.`
}

function verifyPrompt(ev, c) {
  return `At locus ${c.gene_id}, rustle-VG matched RefSeq transcripts that StringTie missed: ${JSON.stringify(ev.vg_gain_vs_st)} (a candidate win via multi-mapping reads). Adversarially verify this is genuinely read-backed at THIS copy and not a DAZ3-style phantom borrowed from a sibling paralog. Run this EXACT command (it reuses the locus.bam already extracted by the eval step and computes strict edit-distance decisiveness — how many reads' sequence STRICTLY prefers this copy over siblings):

GI='${GI}' bash ${SH}/winloci_verify.sh '${c.gene_id}' '${c.chrom}' ${c.start} ${c.end} '${c.strand}'

Return its JSON verdict object via StructuredOutput. verdict meanings: 'real' = enough reads strictly anchor here (genuine recovery StringTie can't make); 'phantom' = no read prefers this copy, sequence favors a sibling (fabrication, DAZ3 pattern — NOT a real win); 'ambiguous'/'weak' = inconclusive. Return verbatim.`
}

log(`Evaluating ${cands.length} paralog candidates with VG config: ${VGFLAGS}`)

const results = await pipeline(
  cands,
  (c) => agent(evalPrompt(c), { label: `eval:${c.gene_id}`, phase: 'Eval', schema: EVAL_SCHEMA }),
  (ev, c) => {
    if (!ev || ev.classification !== 'win_vs_st' || !(ev.vg_gain_vs_st || []).length) return ev
    return agent(verifyPrompt(ev, c), { label: `verify:${c.gene_id}`, phase: 'Verify', schema: VERDICT_SCHEMA })
      .then(v => ({ ...ev, verdict: v }))
  }
)

const ok = results.filter(Boolean)
const winsVsSt = ok.filter(r => r.classification === 'win_vs_st')
const winsVsBase = ok.filter(r => r.classification === 'win_vs_base_only')
const ties = ok.filter(r => r.classification === 'tie')
const regress = ok.filter(r => r.classification === 'regress')
const errors = ok.filter(r => r.classification === 'error')
const confirmed = winsVsSt.filter(r => r.verdict && r.verdict.verdict === 'real')
const phantom = winsVsSt.filter(r => r.verdict && r.verdict.verdict === 'phantom')
const inconclusive = winsVsSt.filter(r => r.verdict && (r.verdict.verdict === 'ambiguous' || r.verdict.verdict === 'weak'))

log(`DONE: ${winsVsSt.length} candidate wins-vs-ST → ${confirmed.length} CONFIRMED read-backed, ${phantom.length} phantom, ${inconclusive.length} inconclusive | ${winsVsBase.length} win-vs-baseline-only | ${ties.length} tie | ${regress.length} regress | ${errors.length} err`)

return {
  total: cands.length,
  evaluated: ok.length,
  summary: {
    win_vs_st: winsVsSt.length,
    confirmed_real: confirmed.length,
    phantom_wins: phantom.length,
    inconclusive_wins: inconclusive.length,
    win_vs_base_only: winsVsBase.length,
    tie: ties.length,
    regress: regress.length,
    error: errors.length,
  },
  confirmed: confirmed.map(r => ({ gene: r.gene, region: r.region, gains: r.vg_gain_vs_st, verdict: r.verdict })),
  phantom: phantom.map(r => ({ gene: r.gene, gains: r.vg_gain_vs_st, verdict: r.verdict })),
  inconclusive: inconclusive.map(r => ({ gene: r.gene, gains: r.vg_gain_vs_st, verdict: r.verdict })),
  win_vs_base_only: winsVsBase.map(r => ({ gene: r.gene, region: r.region, gains: r.vg_gain_vs_base })),
  regressions: regress.map(r => ({ gene: r.gene, lost: r.vg_loss_vs_st })),
}
