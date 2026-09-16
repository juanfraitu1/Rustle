#!/usr/bin/env python3
"""NPIP/TBC1D3 layer order (integration) — apply the verification fixes to the light/heavy builds.

Writes corrected tables NEXT TO the originals (suffix .corrected.tsv; originals untouched):
  light/members.corrected.tsv            member rule applied the same way to both families (fix: members asymmetry)
  light/P.groups.corrected.tsv           in_P genes + every other U gene with a §6ko protein as 'P|other:<gene>' (audit
                                         2026-09-16: PKD1, DHX40); non-coding members are NOT P singletons (fix: P encoding)
  light/P.members_status.corrected.tsv   every corrected member, group_id 'P|NA' when not in P
  light/D.groups.corrected.tsv           layer_dna.py rerun on the corrected members, duplicate-coordinate node fix,
                                         catalog-absent members -> 'D|NA' (fix: per-chromosome catalogs, cheapest route);
                                         + every other U gene that is a catalog node, with its catalog cluster (audit
                                         2026-09-16: TBC1D26 -> MCL24, status clustered_nonmember_group)
  light/D.edges.corrected.tsv            same edge rule as layer_dna.py (pre-MCL --dump-graph edges touching member groups)
  light/C.groups.corrected.tsv           TBC1D3 positional clusters moved out of clade_L1 (fix: positional != clade);
                                         members absent from the reference trees -> 'NA' (not assessed)
  light/universe.corrected.tsv           U = members + P + D + C (S1 dropped: deferred layer), per-layer universe flags
  light/truth_hgnc.corrected.tsv, truth_soto.corrected.tsv, truth_guided.corrected.tsv (E0/E1 filled for every row),
  light/truth_literature_subfamilies.corrected.tsv
  heavy/EXPR.counts.corrected.tsv        + rows for universe genes EXPR lacked (integrate_slim/expr_recount.tsv)

Corrected member rule (both families):
  description  RefSeq description matches /nuclear pore complex[- ]interacting protein/i  (NPIP)
               or /^TBC1 domain family member 3( |[A-Z]|$)/                                (TBC1D3; excludes TBC1D30-32)
  readthrough  a record whose description says 'readthrough' and whose name/description names NPIP or TBC1D3 is a
               member only if its span overlaps NO same-strand description-rule record (one record per copy: the
               readthrough is kept only when its family part has no gene record of its own).
"""
import collections
import csv
import re
import sys

sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light/scripts")
import soto_map  # noqa: E402

sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench")
from protein_families import excluded  # noqa: E402  (§6ko r2 filter)

LIGHT = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light"
HEAVY = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/heavy"
INT = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/integrate_slim"
H = "/mnt/linuxdisk/home/juanfraitu/o1_falsemerge"
HGNC = "/mnt/linuxdisk/home/juanfraitu/winloci_data/hgnc/hgnc_complete_set.txt"
LIT = "/mnt/c/Users/jfris/Desktop/Rustle/docs/lit_subclusters_npip_tbc1d3_truth.tsv"
CATALOGS = {
    "c15_17_22": dict(kind="regions", path=f"{H}/human2/genes.regions", D=f"{LIGHT}/work/D/c15_17_22.e1",
                      E0=f"{H}/human2/guided", E1=f"{H}/lit/aj_dev/refseq_e1"),
    "c16_19_20": dict(kind="nodes", path=f"{H}/lit/aj_ho/refseq/nodes.tsv", D=f"{LIGHT}/work/D/c16_19_20.e1",
                      E0=f"{H}/lit/aj_ho/refseq/e0", E1=f"{H}/lit/aj_ho/refseq/e1"),
}
DESC_NPIP = re.compile(r"nuclear pore complex[- ]interacting protein", re.I)
DESC_TBC = re.compile(r"^TBC1 domain family member 3( |[A-Z]|$)")
FAMNAME = re.compile(r"NPIP|TBC1D3")
SYM = re.compile(r"^(NPIP.*|TBC1D3(?:$|[A-Z]|P\d|-).*)$")  # the original symbol rule (members.py), for the audit line


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def write(path, cols, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def key_of(s):
    c, r = s.rsplit(":", 1)
    a, b = r.split("-")
    return (c, int(a), int(b))


def catalog_keys(by_coord, kind, path):
    """record key -> RefSeq gene rows. Fix vs layer_dna.py/truths_universe.py: duplicate-coordinate node rows accumulate
    (setdefault/extend) instead of the last row overwriting the others."""
    key2genes = {}
    if kind == "nodes":
        names = {r["idx"]: r["name"] for r in tsv(path + ".names.tsv")}
        for r in tsv(path):
            k = (r["chrom"], int(r["start"]) + 1, int(r["end"]))
            lst = key2genes.setdefault(k, [])
            for g in by_coord.get(k, []):
                if g["name"] == names[r["idx"]] and g not in lst:
                    lst.append(g)
    else:
        for line in open(path):
            k = key_of(line.strip())
            key2genes[k] = list(by_coord.get(k, []))
    return key2genes


def membership(key2genes, prefix):
    rep = {key_of(r["annotation"]): key_of(r["representative"]) for r in tsv(prefix + ".loci.tsv")}
    cl = {(r["chrom"], int(r["start"]), int(r["end"])): r["cluster_id"] for r in tsv(prefix + ".clusters.tsv")}
    out = {}
    for k, gs in key2genes.items():
        r = rep.get(k, k)
        for g in gs:
            out[g["gene_id"]] = (cl.get(r, ""), f"{k[0]}:{k[1]}-{k[2]}", f"{r[0]}:{r[1]}-{r[2]}" if r != k else "")
    return out


def main():
    genes = {r["gene_id"]: r for r in tsv(f"{LIGHT}/work/refseq/genes.tsv")}
    by_coord = collections.defaultdict(list)
    for g in genes.values():
        by_coord[(g["chrom"], int(g["start0"]) + 1, int(g["end"]))].append(g)

    # ------------------------------------------------------------------ A. members
    old_members = {r["gene_id"]: r for r in tsv(f"{LIGHT}/members.tsv")}
    desc = {}
    for gid, g in genes.items():
        if DESC_NPIP.search(g["description"]):
            desc[gid] = "NPIP"
        elif DESC_TBC.search(g["description"]):
            desc[gid] = "TBC1D3"
    rt_cand = {gid: ("NPIP" if "NPIP" in (g["name"] + " " + g["description"]) else "TBC1D3")
               for gid, g in genes.items()
               if "readthrough" in g["description"] and FAMNAME.search(g["name"] + " " + g["description"])}
    members, rt_log = {}, []
    for gid, fam in desc.items():
        members[gid] = (fam, "description" + ("; symbol" if SYM.match(genes[gid]["name"]) else ""))
    for gid, fam in sorted(rt_cand.items()):
        g = genes[gid]
        s, e = int(g["start0"]), int(g["end"])
        ov = sorted(genes[d]["name"] for d in desc if genes[d]["chrom"] == g["chrom"] and genes[d]["strand"] == g["strand"]
                    and int(genes[d]["start0"]) < e and s < int(genes[d]["end"]))
        rt_log.append(f"{g['name']} ({g['description']}): overlaps same-strand description-rule records {ov or '-'} -> "
                      f"{'NOT member' if ov else 'member'}")
        if not ov:
            members[gid] = (fam, "readthrough (family part has no own gene record)")
    sym_only = [genes[g]["name"] for g in genes if SYM.match(genes[g]["name"]) and g not in members and g not in rt_cand]
    lit = {r["name"]: r for r in tsv(LIT)}
    db = soto_map.load()
    exons = {r["gene_id"]: soto_map.parse_blocks(r["exons"]) for r in tsv(f"{LIGHT}/work/refseq/exons.tsv")}
    mrows = []
    for gid, (fam, basis) in members.items():
        g = genes[gid]
        m = soto_map.map_gene(db, g["name"], g["chrom"], g["strand"], exons.get(gid) or [(int(g["start0"]), int(g["end"]))])
        mrows.append({"gene_id": gid, "name": g["name"], "biotype": g["biotype"], "chrom": g["chrom"], "start": g["start0"],
                      "end": g["end"], "strand": g["strand"], "family": fam, "member_basis": basis,
                      "in_lit_truth_31": "yes" if g["name"] in lit else "no", "description": g["description"],
                      "in_original_members": "yes" if gid in old_members else "no", **m})
    mrows.sort(key=lambda x: (x["family"], x["chrom"], int(x["start"])))
    # audit (2026-09-16, minor): the readthrough restriction is not a general 'one record per copy' rule. Flag every member
    # whose span lies inside another same-strand member record (non-readthrough containment), and report the member count
    # under the alternative 'fold the enclosed record into its encloser' option. No layer is changed by this column.
    for x in mrows:
        s, e = int(x["start"]), int(x["end"])
        enc = sorted(y["name"] for y in mrows if y is not x and y["chrom"] == x["chrom"] and y["strand"] == x["strand"]
                     and int(y["start"]) <= s and e <= int(y["end"]) and (int(y["end"]) - int(y["start"])) > (e - s))
        x["span_inside_member_record"] = ",".join(enc)
    n_enc = collections.Counter(x["family"] for x in mrows if x["span_inside_member_record"])
    print(f"[members] records whose span lies inside another same-strand member record: "
          f"{[(x['name'], x['span_inside_member_record']) for x in mrows if x['span_inside_member_record']]}; member count if "
          f"they are folded into the encloser: NPIP {sum(x['family'] == 'NPIP' for x in mrows) - n_enc['NPIP']}, TBC1D3 "
          f"{sum(x['family'] == 'TBC1D3' for x in mrows) - n_enc['TBC1D3']}")
    mcols = ["gene_id", "name", "biotype", "chrom", "start", "end", "strand", "family", "member_basis", "in_lit_truth_31",
             "description", "in_original_members", "span_inside_member_record"] + soto_map.FIELDS
    write(f"{LIGHT}/members.corrected.tsv", mcols, mrows)
    M = {r["gene_id"]: r for r in mrows}
    fam_of = {g: r["family"] for g, r in M.items()}
    print(f"[members] corrected {len(M)}: NPIP {sum(v == 'NPIP' for v in fam_of.values())}, TBC1D3 "
          f"{sum(v == 'TBC1D3' for v in fam_of.values())}; added {sorted(genes[g]['name'] for g in set(M) - set(old_members))}; "
          f"removed {sorted(genes[g]['name'] for g in set(old_members) - set(M))}; symbol-rule names outside the rule {sym_only or '-'}")
    for x in rt_log:
        print("   readthrough:", x)

    # ------------------------------------------------------------------ B. P
    pidx = tsv(f"{LIGHT}/work/P/proteins.index.tsv")
    P_universe = {r["gene_id"] for r in pidx}
    Prow = [r for r in tsv(f"{LIGHT}/P.groups.tsv") if r["p_status"] == "in_P"]
    for r in Prow:
        r["is_member"] = "yes" if r["gene_id"] in M else "no"
        r["member_family"] = fam_of.get(r["gene_id"], "")
    pcols = list(tsv(f"{LIGHT}/P.groups.tsv")[0].keys())
    write(f"{LIGHT}/P.groups.corrected.tsv", pcols, Prow)
    Pg = {r["gene_id"]: r["group_id"] for r in Prow}
    old_status = {r["gene_id"]: r for r in tsv(f"{LIGHT}/P.members_status.tsv")}
    cds = {r["gene_id"] for r in tsv(f"{LIGHT}/work/refseq/cds.tsv")}
    srows = []
    for gid, r in sorted(M.items(), key=lambda kv: (kv[1]["family"], kv[1]["name"])):
        if gid in Pg:
            st, grp = "in_P", Pg[gid]
        elif gid in old_status:
            st, grp = old_status[gid]["p_status"], "P|NA"
        else:
            st = ";".join(x for x in ("excluded_r2_biotype" if excluded(r["biotype"], 2) else "",
                                      "no_CDS_in_RefSeq" if gid not in cds else "") if x) or "CHECK"
            grp = "P|NA"
        srows.append({"gene_id": gid, "name": r["name"], "biotype": r["biotype"], "family": r["family"], "p_status": st,
                      "group_id": grp, "in_P_universe": "yes" if gid in P_universe else "no"})
    write(f"{LIGHT}/P.members_status.corrected.tsv",
          ["gene_id", "name", "biotype", "family", "p_status", "group_id", "in_P_universe"], srows)
    assert all((r["group_id"] != "P|NA") == (r["in_P_universe"] == "yes") for r in srows), "member in P universe but not in_P"
    print(f"[P] in_P rows {len(Prow)} (members {sum(r['is_member'] == 'yes' for r in Prow)}); members not in P "
          f"{sum(r['group_id'] == 'P|NA' for r in srows)}: {dict(collections.Counter(r['p_status'] for r in srows if r['group_id'] == 'P|NA'))}")

    # ------------------------------------------------------------------ C. D (layer_dna.py logic, corrected members)
    D_all, E0_all, E1_all, D_catalog_of = {}, {}, {}, {}
    drows, erows = [], []
    for tag, c in CATALOGS.items():
        k2g = catalog_keys(by_coord, c["kind"], c["path"])
        mem = membership(k2g, c["D"])
        e0 = membership(k2g, c["E0"])
        e1 = membership(k2g, c["E1"])
        assert all(mem[g][0] == e1[g][0] for g in mem), f"D != E1 on {tag}"
        for g in mem:
            D_catalog_of[g] = tag
            D_all[g] = mem[g]
            E0_all[g] = e0[g]
            E1_all[g] = e1[g]
        mclus = {mem[m][0] for m in M if m in mem and mem[m][0]}
        in_groups = sorted((g for g in mem if mem[g][0] in mclus),
                           key=lambda x: (mem[x][0], genes[x]["chrom"], int(genes[x]["start0"]), x))
        for gid in in_groups:
            g = genes[gid]
            cl, k, f = mem[gid]
            drows.append({"gene_id": gid, "group_id": f"D|{tag}|{cl}", "name": g["name"], "biotype": g["biotype"],
                          "chrom": g["chrom"], "start": g["start0"], "end": g["end"], "strand": g["strand"],
                          "is_member": "yes" if gid in M else "no", "member_family": fam_of.get(gid, ""),
                          "catalog": tag, "cluster_id": cl, "record_key": k, "folded_into": f, "status": "clustered"})
        for m in M:
            if m in mem and not mem[m][0]:
                g = genes[m]
                drows.append({"gene_id": m, "group_id": f"D|{tag}|singleton:{m}", "name": g["name"], "biotype": g["biotype"],
                              "chrom": g["chrom"], "start": g["start0"], "end": g["end"], "strand": g["strand"],
                              "is_member": "yes", "member_family": fam_of[m], "catalog": tag, "cluster_id": "",
                              "record_key": mem[m][1], "folded_into": mem[m][2], "status": "in_catalog_unclustered"})
        keyset = {key_of(mem[g][1]) for g in in_groups}
        grp = {g: f"D|{tag}|{mem[g][0]}" for g in mem if mem[g][0]}
        for line in open(c["D"] + ".graph.tsv"):
            f = line.rstrip("\n").split("\t")
            if len(f) != 3:
                continue
            a, b = key_of(f[0]), key_of(f[1])
            if a not in keyset and b not in keyset:
                continue
            for ga in k2g.get(a, []) or [{"gene_id": "?", "name": "?"}]:
                for gb in k2g.get(b, []) or [{"gene_id": "?", "name": "?"}]:
                    erows.append({"catalog": tag, "u_gene_id": ga["gene_id"], "u_name": ga["name"],
                                  "v_gene_id": gb["gene_id"], "v_name": gb["name"], "weight": f[2],
                                  "u_group": grp.get(ga["gene_id"], ""), "v_group": grp.get(gb["gene_id"], ""),
                                  "u_key": f[0], "v_key": f[1]})
        print(f"[D] {tag}: keys {len(k2g)} (genes mapped {len(mem)}); member clusters {sorted(mclus)}; genes in them "
              f"{len(in_groups)} (members {sum(g in M for g in in_groups)})")
    for m in M:
        if m not in D_all:
            g = genes[m]
            drows.append({"gene_id": m, "group_id": "D|NA", "name": g["name"], "biotype": g["biotype"], "chrom": g["chrom"],
                          "start": g["start0"], "end": g["end"], "strand": g["strand"], "is_member": "yes",
                          "member_family": fam_of[m], "catalog": "", "cluster_id": "", "record_key": "", "folded_into": "",
                          "status": "not_in_any_catalog (outside D's universe)"})
    dcols = ["gene_id", "group_id", "name", "biotype", "chrom", "start", "end", "strand", "is_member", "member_family",
             "catalog", "cluster_id", "record_key", "folded_into", "status"]
    write(f"{LIGHT}/D.groups.corrected.tsv", dcols, drows)
    write(f"{LIGHT}/D.edges.corrected.tsv",
          ["catalog", "u_gene_id", "u_name", "v_gene_id", "v_name", "weight", "u_group", "v_group", "u_key", "v_key"], erows)
    old_d = {(r["gene_id"], r["group_id"]) for r in tsv(f"{LIGHT}/D.groups.tsv") if r["status"] == "clustered"}
    new_d = {(r["gene_id"], r["group_id"]) for r in drows if r["status"] == "clustered"}
    old_e = {tuple(r[c] for c in ("u_gene_id", "v_gene_id", "weight")) for r in tsv(f"{LIGHT}/D.edges.tsv")}
    new_e = {tuple(r[c] for c in ("u_gene_id", "v_gene_id", "weight")) for r in erows}
    print(f"[D] clustered (gene, group) rows original {len(old_d)} corrected {len(new_d)} identical {old_d == new_d}; "
          f"edge rows original {len(old_e)} corrected {len(new_e)} identical {old_e == new_e}; not in any catalog: "
          f"{sorted(genes[r['gene_id']]['name'] for r in drows if r['group_id'] == 'D|NA')}")
    Dg = {r["gene_id"]: r["group_id"] for r in drows if r["group_id"] != "D|NA"}

    # ------------------------------------------------------------------ D. C
    oc = {r["gene_id"]: r for r in tsv(f"{LIGHT}/C.groups.tsv")}
    crows = []
    for gid, r in sorted(M.items(), key=lambda kv: (kv[1]["family"], kv[1]["chrom"], int(kv[1]["start"]))):
        o = oc.get(gid)
        if o is None or o["in_reference_trees"] != "yes":
            crows.append({"gene_id": gid, "name": r["name"], "family": r["family"], "in_reference_trees": "no",
                          "C_fine": "C|NA", "C_mid": "Cmid|NA", "C_L1": "CL1|NA", "positional_cluster": "",
                          "evidence_fine": "", "evidence_mid": "", "evidence_L1": "",
                          "note": "not in the 31-record reference trees: not assessed (outside C's universe)"})
            continue
        l1, ev1, pos, note = o["clade_L1"], o["evidence_L1"], "", ""
        if r["family"] == "TBC1D3" and not l1.startswith("CL1|singleton"):
            pos = f"{l1[4:]} (intron tree {ev1.split('SH-aLRT/UFBoot ')[-1]}; POSITIONAL, not a clade: o1_ledger §6jp/§6jr/§6js, L14917-14923)"
            l1, ev1, note = f"CL1|singleton:{gid}", "", "clade_L1 cleared: TBC1D3 cluster1/2 are positional"
        crows.append({"gene_id": gid, "name": r["name"], "family": r["family"], "in_reference_trees": "yes",
                      "C_fine": o["clade_id"], "C_mid": o["clade_mid"], "C_L1": l1, "positional_cluster": pos,
                      "evidence_fine": o["evidence_fine"], "evidence_mid": o["evidence_mid"], "evidence_L1": ev1,
                      "note": note})
    ccols = ["gene_id", "name", "family", "in_reference_trees", "C_fine", "C_mid", "C_L1", "positional_cluster",
             "evidence_fine", "evidence_mid", "evidence_L1", "note"]
    write(f"{LIGHT}/C.groups.corrected.tsv", ccols, crows)
    Cr = {r["gene_id"]: r for r in crows if r["in_reference_trees"] == "yes"}
    print(f"[C] rows {len(crows)}; in trees {len(Cr)}; TBC1D3 clade_L1 cleared "
          f"{sum(1 for r in crows if r['positional_cluster'])}")

    # ------------------------------------------------------------------ E. universe
    memP = {Pg[g] for g in Pg if g in M}
    memD = {Dg[g] for g in Dg if g in M and "singleton" not in Dg[g]}
    U = set(M) | {g for g in Pg if Pg[g] in memP} | {g for g in Dg if Dg[g] in memD} | set(Cr)
    side = {}
    for g in U:
        fams = set()
        if g in M:
            fams.add(fam_of[g])
        for Lg in (Pg, Dg):
            if g in Lg:
                fams |= {fam_of[m] for m in M if Lg.get(m) == Lg[g]}
        assert len(fams) == 1, (g, fams)
        side[g] = fams.pop()
    expr = {r["gene_id"]: r for r in tsv(f"{HEAVY}/EXPR.counts.tsv")}
    rec = {r["gene_id"]: r for r in tsv(f"{INT}/expr_recount.tsv")}
    hg = tsv(HGNC)
    by_id = {r["hgnc_id"]: r for r in hg}
    by_sym = {r["symbol"]: r for r in hg}
    dbx = dict(line.rstrip("\n").split("\t") for line in open(f"{LIGHT}/work/refseq/gene_dbxref.tsv"))
    old_soto = {r["gene_id"]: r for r in tsv(f"{LIGHT}/truth_soto.tsv")}
    urows, hrows, srows2, grows = [], [], [], []
    order = sorted(U, key=lambda g: (genes[g]["chrom"], int(genes[g]["start0"]), g))
    for g in order:
        r = genes[g]
        name = r["name"]
        ex = rec.get(name)
        assert ex is not None, name
        ids = [x[5:] for x in dbx.get(g, "").split(",") if x.startswith("HGNC:")]
        h, how = None, "none"
        if ids and ids[0] in by_id:
            h, how = by_id[ids[0]], "dbxref"
        elif name in by_sym:
            h, how = by_sym[name], "symbol"
        hrows.append({"gene_id": g, "name": name, "is_member": "yes" if g in M else "no", "match": how,
                      "hgnc_id": h["hgnc_id"] if h else "", "hgnc_symbol": h["symbol"] if h else "",
                      "locus_group": h["locus_group"] if h else "", "gene_group": h["gene_group"] if h else "",
                      "gene_group_id": h["gene_group_id"] if h else ""})
        if g in old_soto:
            s = dict(old_soto[g])
            s["is_member"] = "yes" if g in M else "no"
        else:
            m = soto_map.map_gene(db, name, r["chrom"], r["strand"], exons.get(g) or [(int(r["start0"]), int(r["end"]))])
            flag = ("unmatched" if not m["soto_gene_id"] else "weak_match" if m["soto_match_quality"] == "weak"
                    else "ambiguous_multi_family" if m["soto_ambiguous"] == "yes" else "ok")
            s = {"gene_id": g, "name": name, "is_member": "yes" if g in M else "no", "flag": flag, **m}
        srows2.append(s)
        c = Cr.get(g)
        layers = [L for L, ok in (("P", Pg.get(g) in memP), ("D", Dg.get(g) in memD), ("C", g in Cr)) if ok]
        urows.append({"gene_id": g, "name": name, "biotype": r["biotype"], "chrom": r["chrom"], "start": r["start0"],
                      "end": r["end"], "strand": r["strand"], "is_member": "yes" if g in M else "no",
                      "member_family": fam_of.get(g, ""), "family_side": side[g],
                      "in_P_universe": "yes" if g in P_universe else "no",
                      "P_group": Pg.get(g, f"P|other:{g}" if g in P_universe else "P|NA"),
                      "in_D_universe": "yes" if g in D_all else "no",
                      "D_group": (f"D|{D_catalog_of[g]}|{D_all[g][0]}" if D_all[g][0] else f"D|{D_catalog_of[g]}|singleton:{g}")
                      if g in D_all else "D|NA",
                      "D_folded_into": D_all[g][2] if g in D_all else "",
                      "in_C_universe": "yes" if c else "no",
                      "C_L1": c["C_L1"] if c else "CL1|NA", "C_mid": c["C_mid"] if c else "Cmid|NA",
                      "C_fine": c["C_fine"] if c else "C|NA", "positional_cluster": c["positional_cluster"] if c else "",
                      "layers_placing_gene_with_member": ",".join(layers),
                      "E0_group": (f"E0|{D_catalog_of[g]}|{E0_all[g][0]}" if E0_all[g][0] else f"E0|{D_catalog_of[g]}|singleton:{g}")
                      if g in D_all else "",
                      "hgnc_gene_group_id": h["gene_group_id"] if h else "", "hgnc_gene_group": h["gene_group"] if h else "",
                      "soto_families": s["soto_families"], "soto_flag": s["flag"],
                      "n_reads_any": ex["n_reads_any"], "n_reads_unique": ex["n_reads_unique"],
                      "description": r["description"]})
    ucols = list(urows[0].keys())
    write(f"{LIGHT}/universe.corrected.tsv", ucols, urows)
    # audit (2026-09-16, important): the group tables must carry every U gene of the layer's universe, not only the genes
    # of member groups, so that containment computed from P.groups / D.groups equals containment computed from the
    # universe table. P universe rule (spec, 'P defined on coding genes'): every U gene with a §6ko protein
    # (work/P/proteins.index.tsv), including coding genes P does not place with a member ('P|other:<gene>').
    # D universe: every U gene that is a node (record) of one of the two catalogs, with its catalog cluster.
    plen = {r["gene_id"]: r["length_aa"] for r in pidx}
    have_p = {r["gene_id"] for r in Prow}
    extra_p = []
    for g in order:
        if g in P_universe and g not in have_p:
            r = genes[g]
            extra_p.append({"gene_id": g, "group_id": f"P|other:{g}", "name": r["name"], "biotype": r["biotype"],
                            "chrom": r["chrom"], "start": r["start0"], "end": r["end"], "strand": r["strand"],
                            "protein_len_aa": plen[g], "is_member": "yes" if g in M else "no",
                            "member_family": fam_of.get(g, ""), "p_status": "in_P_universe_not_with_member"})
    write(f"{LIGHT}/P.groups.corrected.tsv", pcols, Prow + extra_p)
    have_d = {r["gene_id"] for r in drows}
    extra_d = []
    for g in order:
        if g in D_all and g not in have_d:
            r = genes[g]
            cl, k, f = D_all[g]
            extra_d.append({"gene_id": g, "group_id": f"D|{D_catalog_of[g]}|{cl}" if cl else f"D|{D_catalog_of[g]}|singleton:{g}",
                            "name": r["name"], "biotype": r["biotype"], "chrom": r["chrom"], "start": r["start0"],
                            "end": r["end"], "strand": r["strand"], "is_member": "yes" if g in M else "no",
                            "member_family": fam_of.get(g, ""), "catalog": D_catalog_of[g], "cluster_id": cl,
                            "record_key": k, "folded_into": f,
                            "status": "clustered_nonmember_group" if cl else "in_catalog_unclustered"})
    write(f"{LIGHT}/D.groups.corrected.tsv", dcols, drows + extra_d)
    print(f"[P/D group tables] rows added for U genes of the layer universe outside member groups: P "
          f"{[(genes[x['gene_id']]['name'], x['group_id']) for x in extra_p]}; D "
          f"{[(genes[x['gene_id']]['name'], x['group_id'], x['folded_into']) for x in extra_d]}")
    pu_tab = {r["gene_id"]: r["group_id"] for r in Prow + extra_p}
    du_tab = {r["gene_id"]: r["group_id"] for r in drows + extra_d if r["group_id"] != "D|NA"}
    assert pu_tab == {u["gene_id"]: u["P_group"] for u in urows if u["in_P_universe"] == "yes"}, "P table != universe"
    assert du_tab == {u["gene_id"]: u["D_group"] for u in urows if u["in_D_universe"] == "yes"}, "D table != universe"
    write(f"{LIGHT}/truth_hgnc.corrected.tsv", list(hrows[0].keys()), hrows)
    write(f"{LIGHT}/truth_soto.corrected.tsv", ["gene_id", "name", "is_member", "flag"] + soto_map.FIELDS, srows2)
    # guided truth: every U gene in a catalog + every gene sharing an E0/E1 cluster with a member, cluster ALWAYS filled
    mcl0 = {(D_catalog_of[m], E0_all[m][0]) for m in M if m in E0_all and E0_all[m][0]}
    mcl1 = {(D_catalog_of[m], E1_all[m][0]) for m in M if m in E1_all and E1_all[m][0]}
    gl = [g for g in D_all if g in U or (D_catalog_of[g], E0_all[g][0]) in mcl0 or (D_catalog_of[g], E1_all[g][0]) in mcl1]
    for g in sorted(gl, key=lambda x: (genes[x]["chrom"], int(genes[x]["start0"]), x)):
        r = genes[g]
        grows.append({"gene_id": g, "name": r["name"], "biotype": r["biotype"], "chrom": r["chrom"], "start": r["start0"],
                      "end": r["end"], "is_member": "yes" if g in M else "no", "in_universe": "yes" if g in U else "no",
                      "catalog": D_catalog_of[g],
                      "E0_cluster": f"{D_catalog_of[g]}|{E0_all[g][0]}" if E0_all[g][0] else "unclustered",
                      "E0_folded_into": E0_all[g][2],
                      "E1_cluster": f"{D_catalog_of[g]}|{E1_all[g][0]}" if E1_all[g][0] else "unclustered",
                      "E1_folded_into": E1_all[g][2], "record_key": E0_all[g][1]})
    write(f"{LIGHT}/truth_guided.corrected.tsv", list(grows[0].keys()), grows)
    # literature
    lrows = []
    for gid, r in sorted(M.items(), key=lambda kv: (kv[1]["family"], kv[1]["chrom"], int(kv[1]["start"]))):
        t = lit.get(r["name"])
        lrows.append({"gene_id": gid, "name": r["name"], "family": r["family"], "in_literature_truth": "yes" if t else "no",
                      "level1": t["level1"] if t else "", "level2": t["level2"] if t else "",
                      "level1_kind": ("" if not t else "phylogenetic subfamily (Dishuck 2025)" if r["family"] == "NPIP"
                                      else "POSITIONAL genomic cluster (Guitart 2024; o1_ledger L14917-14923)"),
                      "npipb_named_subfamily": ("yes" if r["name"] in {"NPIPB3", "NPIPB4", "NPIPB5", "NPIPB11", "NPIPB12",
                                                                        "NPIPB13"} else "no") if r["family"] == "NPIP" and t else ""})
    write(f"{LIGHT}/truth_literature_subfamilies.corrected.tsv", list(lrows[0].keys()), lrows)
    # EXPR corrected
    xrows = []
    names_U = {genes[g]["name"]: g for g in U}
    for n, r in sorted(rec.items()):
        o = expr.get(n)
        gid = "gene-" + n
        row = dict(o) if o else {"gene_id": n, "name": n, "chrom": genes[gid]["chrom"], "n_reads_any": r["n_reads_any"],
                                 "n_reads_unique": r["n_reads_unique"],
                                 "expressed_ge3_any": "yes" if int(r["n_reads_any"]) >= 3 else "no",
                                 "expressed_ge3_unique": "yes" if int(r["n_reads_unique"]) >= 3 else "no",
                                 "seed_family": "", "biotype": genes[gid]["biotype"], "in_S2": "no", "S2_group_id": ""}
        row["member_family_corrected"] = fam_of.get(gid, "")
        row["in_universe_corrected"] = "yes" if n in names_U else "no"
        row["row_source"] = "EXPR.counts.tsv" if o else "integrate_slim/expr_recount.tsv (same rule)"
        xrows.append(row)
    write(f"{HEAVY}/EXPR.counts.corrected.tsv", list(xrows[0].keys()), xrows)
    print(f"[U] {len(U)} genes (members {len(M)}); per layer universe: P {sum(u['in_P_universe'] == 'yes' for u in urows)}, "
          f"D {sum(u['in_D_universe'] == 'yes' for u in urows)}, C {sum(u['in_C_universe'] == 'yes' for u in urows)}; "
          f"placed with a member by P {sum('P' in u['layers_placing_gene_with_member'].split(',') for u in urows)}, "
          f"D {sum('D' in u['layers_placing_gene_with_member'].split(',') for u in urows)}; family side "
          f"{dict(collections.Counter(side.values()))}; chromosomes {dict(collections.Counter(genes[g]['chrom'] for g in U))}")
    print(f"[truths] HGNC match {dict(collections.Counter(h['match'] for h in hrows))}; Soto flags "
          f"{dict(collections.Counter(s['flag'] for s in srows2))}; guided rows {len(grows)}; EXPR corrected rows {len(xrows)}")


if __name__ == "__main__":
    main()
