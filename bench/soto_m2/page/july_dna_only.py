"""DNA-only versions of the two July tabs (2026-09-30, user request: compare DNA to DNA, as Soto did; the A119b Iso-Seq results are
removed from the meeting page). Detection opens and stays in DNA mode (--from-genome: genome self-alignment, no reads); the members
table loses its RNA tile, filter and column. Each edit must match exactly once."""
import re


def _rep(s, a, b):
    assert s.count(a) == 1, (s.count(a), a[:80])
    return s.replace(a, b)


def det(markup, js):
    markup = _rep(markup, '<button class="modebtn" data-m="rna" aria-pressed="true">RNA mode</button>'
                          '<button class="modebtn" data-m="dna" aria-pressed="false">DNA mode</button>',
                  '<button class="modebtn" data-m="dna" aria-pressed="true">DNA mode</button>')
    markup = _rep(markup, "Rustle × Soto 2025 · human A119b · T2T-CHM13 v2.0", "Rustle × Soto 2025 · DNA mode, genome only · T2T-CHM13 v2.0")
    lede = re.search(r'<p class="lede">.*?</p>', markup, re.S).group(0)
    markup = _rep(markup, lede, '<p class="lede">Rustle\'s family engine run on the genome alone (DNA mode: self-alignment + '
                  'γ-quasi-cliques, no reads, no annotation) against the Soto segmental-duplication truth set — <b>83 families, 362 '
                  'members</b>. Every member becomes a locus grouped by homology, and the method also surfaces <b>extra copies</b> '
                  '(unlisted paralogs) that lower raw precision. Expand a family to review members and candidates; uncheck segdup '
                  'pieces or candidates you judge real — and use the size-outlier lever on the extra copies. Both metrics update '
                  'live.</p>')
    markup = _rep(markup, "The 5 members marked <b>· top-up</b> are missed on the real BAM but found when the member's own reads are "
                          "resampled to ideal coverage (a simulated top-up); they are <b>counted as found</b> (coverage-limited, not a "
                          "fundamental identifiability limit).", "")
    markup = _rep(markup, "<span><i>· top-up</i> = found at ideal coverage (counted as found)</span>", "")
    js = _rep(js, "let mode='rna';", "let mode='dna';")
    return markup, js


def mem(markup, js):
    markup = _rep(markup, " <b>RNA</b> = the Iso-Seq pipeline recovered it.", "")
    markup = _rep(markup, '<div class="stat rn"><div class="k">RNA pipeline</div><div class="v num">313 / 362</div><div class="s num">'
                          '86.5% — a strict subset of the DNA method</div></div>', "")
    markup = _rep(markup, '<button class="chip" data-f="rnamiss" aria-pressed="false">RNA missed <span class="num">49</span></button>', "")
    markup = _rep(markup, '<th data-s="rna">RNA</th>', "")
    markup, n = re.subn(r"RNA: verified per-member attribution\.\s*RNA is a strict subset \(both 313, RNA-only 0\)\.\s*", "", markup)
    assert n == 1, n
    js = _rep(js, "cell(ex)+cell(pipe)+cell(rna)", "cell(ex)+cell(pipe)")
    return markup, js
