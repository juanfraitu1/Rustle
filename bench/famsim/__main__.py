#!/usr/bin/env python3
"""python3 bench/famsim <command> — see bench/FAMSIM.md.

  spec [RUNG]                         print a scenario template (or one rung of the ladder on the template spec)
  make SPEC.json --out DIR            build genome + truth + reads
  verify DIR                          prove the condition from the products -> verify.tsv
  align DIR [--threads 2]             map the reads with the shipped minimap2 command
  run DIR [--stages a,b] [--bin DIR]  run the pipeline stages (assemble,denovo,guided,assign,flag)
  score DIR                           score every stage -> score.tsv, summary.md
  all SPEC.json --out DIR             make + verify + align + run + score
  ladder --template-spec T.json --out DIR [--only r1,r2] [--background ...]   the whole built-in ladder -> ladder.tsv
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from famsim import chromosome, evaluate, pipeline, reads, scenarios, verify  # noqa: E402
from famsim.model import GeneModel, seeded  # noqa: E402
from famsim.template import resolve_decoys, resolve_template  # noqa: E402


def say(msg):
    print(f"[famsim] {msg}", file=sys.stderr, flush=True)


def cmd_make(spec, out):
    planted, man = chromosome.build(spec, out, say)
    reads.simulate(planted, spec.get("reads", {}), out, int(spec.get("seed", 1)), say)
    with open(os.path.join(out, "spec.json"), "w") as fh:
        json.dump(spec, fh, indent=1)
    return planted


def cmd_all(spec, out, a):
    cmd_make(spec, out)
    R = verify.run(out, threads=1, allow_background_homology=a.allow_background_homology, log=say)
    if R.failed and not a.keep_going:
        say(f"{len(R.failed)} verify FAIL(s): stopping before alignment (use --keep-going to continue)")
        return R
    pipeline.align(out, a.threads, say)
    pipeline.run(out, a.stages.split(","), a.bin, a.threads, a.copy_table, say)
    evaluate.run(out, say)
    return R


def cmd_ladder(a):
    base = json.load(open(a.template_spec))
    os.makedirs(a.out, exist_ok=True)
    seed = int(base.get("seed", 1))
    tj, dj = os.path.join(a.out, "template.json"), os.path.join(a.out, "decoys.json")
    if not os.path.exists(tj) or a.force:
        t = resolve_template(base["template"], seeded(seed, "template"))
        t.save(tj)
        ds = resolve_decoys(base.get("decoys"), t, seeded(seed, "decoys"))
        with open(dj, "w") as fh:
            json.dump([d.to_json() for d in ds], fh)
    t = GeneModel.load(tj)
    say(f"ladder template: {t.summary()} [{t.provenance.get('gene') or t.provenance.get('source')}], decoys {len(json.load(open(dj)))}")
    only = set(a.only.split(",")) if a.only else None
    specs = scenarios.ladder_specs(base, tj, dj, len(t.exons), only)
    skipped = [n for n, _, _ in scenarios.rungs(99) if n not in {s["name"] for s in specs} and (not only or n in only)]
    if skipped:
        say(f"rungs not applicable to a {len(t.exons)}-exon template (or not selected): {skipped}")
    rows, hdr = [], None
    lad = os.path.join(a.out, "ladder.tsv")
    if only and os.path.exists(lad):           # a partial rerun replaces its rungs inside the existing table
        with open(lad) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            rows = [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]
    for s in specs:
        d = os.path.join(a.out, s["name"])
        say(f"=== {s['name']}: {s['condition']}")
        try:
            R = cmd_all(s, d, a)
            ver = "PASS" if not R.failed else f"FAIL:{len(R.failed)}"
            S = evaluate.Score(s["name"])
            if os.path.exists(os.path.join(d, "score.tsv")):
                for line in open(os.path.join(d, "score.tsv")):
                    f = line.rstrip("\n").split("\t")
                    if f[0] != "scenario":
                        S.rows.append(tuple(f))
            h = evaluate.headline(S)
        except Exception as e:  # noqa: BLE001  one broken rung must not hide the others
            say(f"    {s['name']} FAILED: {e}")
            ver, h = f"ERROR", {}
        row = {"rung": s["name"], "condition": s["condition"], "verify": ver, **h}
        hdr = hdr or list(row)
        rows = [r for r in rows if r.get("rung") != s["name"]] + [row]
        order = [n for n, _, _ in scenarios.rungs(99)]
        rows.sort(key=lambda r: order.index(r["rung"]) if r.get("rung") in order else 99)
        with open(os.path.join(a.out, "ladder.tsv"), "w") as fh:
            fh.write("\t".join(hdr) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(k, "")) for k in hdr) + "\n")
    say(f"ladder: {len(rows)} rungs -> {os.path.join(a.out, 'ladder.tsv')}")
    for r in rows:
        print("\t".join(str(r.get(k, "")) for k in hdr))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="famsim", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("spec"); p.add_argument("rung", nargs="?")
    for name in ("make", "all"):
        p = sub.add_parser(name); p.add_argument("spec"); p.add_argument("--out", required=True)
    for name in ("verify", "align", "run", "score"):
        p = sub.add_parser(name); p.add_argument("dir")
    p = sub.add_parser("ladder"); p.add_argument("--template-spec", required=True, help="a scenario JSON whose template/background/decoys/reads blocks are shared by every rung")
    p.add_argument("--out", required=True); p.add_argument("--only", default=""); p.add_argument("--force", action="store_true")
    for name in ("run", "all", "ladder"):
        p = sub.choices[name]
        p.add_argument("--stages", default="assemble,denovo,guided,assign,flag")
        p.add_argument("--bin", default=pipeline.DEFAULT_BIN); p.add_argument("--copy-table", default="families", choices=("families", "catalog"))
    for name in ("align", "run", "all", "ladder"):
        sub.choices[name].add_argument("--threads", type=int, default=2)
    for name in ("verify", "all", "ladder"):
        sub.choices[name].add_argument("--allow-background-homology", action="store_true")
    for name in ("all", "ladder"):
        sub.choices[name].add_argument("--keep-going", action="store_true", help="align and run even when verify reports a FAIL")
    a = ap.parse_args(argv)
    if a.cmd == "spec":
        if a.rung:
            base = dict(scenarios.TEMPLATE_SPEC)
            r = next((x for x in scenarios.rungs(6) if x[0] == a.rung), None)
            if r is None:
                sys.exit(f"unknown rung {a.rung}; known: {[x[0] for x in scenarios.rungs(6)]}")
            base["name"], base["condition"], base["copies"] = r[0], r[1], r[2]
            print(json.dumps(base, indent=1))
        else:
            print(json.dumps(scenarios.TEMPLATE_SPEC, indent=1))
    elif a.cmd == "make":
        cmd_make(json.load(open(a.spec)), a.out)
    elif a.cmd == "verify":
        R = verify.run(a.dir, threads=1, allow_background_homology=a.allow_background_homology, log=say)
        sys.exit(1 if R.failed else 0)
    elif a.cmd == "align":
        pipeline.align(a.dir, a.threads, say)
    elif a.cmd == "run":
        pipeline.run(a.dir, a.stages.split(","), a.bin, a.threads, a.copy_table, say)
    elif a.cmd == "score":
        evaluate.run(a.dir, say)
    elif a.cmd == "all":
        R = cmd_all(json.load(open(a.spec)), a.out, a)
        sys.exit(1 if R.failed and not a.keep_going else 0)
    elif a.cmd == "ladder":
        cmd_ladder(a)


if __name__ == "__main__":
    main()
