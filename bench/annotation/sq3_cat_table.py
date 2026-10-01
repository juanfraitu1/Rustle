# Same columns as bench/SQANTI3_POLISH.md, RefSeq (sq3/) vs CAT (sq3_cat/) per chromosome and pooled
import csv, collections
B = "/mnt/linuxdisk/home/juanfraitu/bakeoff"
ART = {"genic", "antisense", "intergenic", "genic_intron", "fusion"}
def stats(c, a, sub):
    s = collections.Counter()
    for r in csv.DictReader(open(f"{B}/human_chr{c}/{sub}/filt_{a}/{a}_RulesFilter_result_classification.txt"), delimiter="\t"):
        k = r["structural_category"]; p = r["filter_result"] == "Isoform"
        s["n"] += 1; s["FSM"] += k == "full-splice_match"; s["ISM"] += k == "incomplete-splice_match"
        s["NIC"] += k == "novel_in_catalog"; s["NNC"] += k == "novel_not_in_catalog"; s["art"] += k in ART
        s["PASS"] += p; s["FSMp"] += p and k == "full-splice_match"
    return s
chroms = (20, 11, 7, 14, 5, 9); arms = ("SHIP", "stringtie", "flair")
pool = {(a, sub): collections.Counter() for a in arms for sub in ("sq3", "sq3_cat")}
rows = []
for c in chroms:
    for a in arms:
        for sub in ("sq3", "sq3_cat"):
            s = stats(c, a, sub); pool[(a, sub)].update(s); rows.append((f"chr{c}", a, sub, s))
for c, a, sub, s in rows + [("pooled", a, sub, pool[(a, sub)]) for a in arms for sub in ("sq3", "sq3_cat")]:
    n = s["n"]; pc = lambda k: f'{s[k]} ({100*s[k]/n:.1f}%)'
    print("\t".join([c, a, "RefSeq" if sub == "sq3" else "CAT", str(n), pc("FSM"), pc("ISM"), pc("NIC"), pc("NNC"), pc("art"), pc("PASS"), str(s["FSMp"])]))
