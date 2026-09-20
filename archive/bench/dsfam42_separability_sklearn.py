#!/usr/bin/env python3
"""Stage 2 (sklearn): can ANYTHING separate the PSV-tied reads by copy? Two tests:
 (1) SUPERVISED on the ASSIGNED reads (copy known via PSV): can non-PSV BAM features predict the copy?
     Run TWICE — with all features (upper bound), and EXCLUDING AS/NM/de/mapq (the metrics that are
     TIED for the unassignable reads). If the second is at chance, nothing orthogonal to AS/PSV
     separates copies → the tied reads are genuinely unassignable.
 (2) UNSUPERVISED on the TIED reads: do they cluster into copy-like groups (silhouette), and does any
     clustering recover the minimap2 placement (adjusted Rand index)?
Run with the sqanti3 env python (has sklearn)."""
import csv, numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.preprocessing import StandardScaler

rows = list(csv.DictReader(open("/tmp/dsfam42_features.tsv"), delimiter="\t"))
FEATS = ["exonic","n_exons","gspan","softclip_frac","sc5","sc3","is_rev","mapq","AS","NM","de",
         "gc","qual_mean","qual_min","seq_entropy","offset5","offset3"]
ORTHO = [f for f in FEATS if f not in ("AS","NM","de","mapq")]   # tied reads have AS/PSV tied
def mat(rs, feats): return np.array([[float(r[f]) for f in feats] for r in rs])

asn = [r for r in rows if r["status"]=="assigned"]
tied = [r for r in rows if r["status"]=="tied"]
ya = np.array([int(r["assigned_copy"]) for r in asn])
print(f"assigned={len(asn)} (across {len(set(ya))} copies)  tied={len(tied)}")
from collections import Counter
print("  assigned-copy balance:", dict(Counter(ya.tolist())))

def supervised(feats, tag):
    X = StandardScaler().fit_transform(mat(asn, feats)); y = ya
    cv = StratifiedKFold(5, shuffle=True, random_state=0)
    rf = RandomForestClassifier(n_estimators=300, random_state=0, class_weight="balanced")
    acc = cross_val_score(rf, X, y, cv=cv, scoring="accuracy").mean()
    dum = cross_val_score(DummyClassifier(strategy="most_frequent"), X, y, cv=cv, scoring="accuracy").mean()
    f1 = cross_val_score(rf, X, y, cv=cv, scoring="f1_macro").mean()
    rf.fit(X, y)
    imp = sorted(zip(feats, rf.feature_importances_), key=lambda x:-x[1])[:6]
    print(f"\n[{tag}] features={len(feats)}")
    print(f"  RandomForest CV accuracy = {acc:.3f}  |  chance (most-frequent) = {dum:.3f}  |  1/k = {1/len(set(y)):.3f}  |  macro-F1 = {f1:.3f}")
    print(f"  top features: " + ", ".join(f"{f}={i:.2f}" for f,i in imp))
    return acc, dum

print("\n=== TEST 1: SUPERVISED — can non-PSV features predict the copy (on ASSIGNED reads)? ===")
a_all,d_all = supervised(FEATS, "ALL features (incl AS/NM/de — upper bound, NOT usable for tied)")
a_o,d_o = supervised(ORTHO, "ORTHOGONAL only (excl AS/NM/de/mapq — what COULD help tied reads)")

print("\n=== TEST 1b: SUPERVISED (non-confounded) — do SEQUENCE/QUALITY features alone predict which")
print("    copy LOCUS minimap2 placed a read at (9-class, ALL reads)? If chance, the copies are")
print("    sequence-indistinguishable -> the floor. (Position/score features excluded.) ===")
SEQ = ["gc","qual_mean","qual_min","seq_entropy","n_exons","exonic","softclip_frac","sc5","sc3","is_rev"]
allr = [r for r in rows if r["status"] != "ambiguous"]
Xs = StandardScaler().fit_transform(mat(allr, SEQ))
ys = np.array([int(r["best_copy"]) for r in allr])
cv = StratifiedKFold(5, shuffle=True, random_state=0)
rf = RandomForestClassifier(n_estimators=300, random_state=0, class_weight="balanced")
acc_s = cross_val_score(rf, Xs, ys, cv=cv, scoring="accuracy").mean()
f1_s = cross_val_score(rf, Xs, ys, cv=cv, scoring="f1_macro").mean()
dum_s = cross_val_score(DummyClassifier(strategy="most_frequent"), Xs, ys, cv=cv, scoring="accuracy").mean()
print(f"  predict minimap2 copy-locus (k={len(set(ys))}) from SEQUENCE/QUALITY only: acc={acc_s:.3f} "
      f"f1_macro={f1_s:.3f} | chance={dum_s:.3f} | 1/k={1/len(set(ys)):.3f}")
print(f"  -> {'AT CHANCE: no sequence feature separates the copies (floor)' if f1_s < 0.2 else 'above chance: residual signal — a lead'}")

print("\n=== TEST 2: UNSUPERVISED — do the TIED reads cluster into copy-like groups? ===")
Xt = StandardScaler().fit_transform(mat(tied, ORTHO))
k = len(set(ya))
km = KMeans(n_clusters=k, n_init=10, random_state=0).fit(Xt)
sil = silhouette_score(Xt, km.labels_)
bc = np.array([int(r["best_copy"]) for r in tied])
ari = adjusted_rand_score(bc, km.labels_)
print(f"  KMeans k={k} on tied reads: silhouette = {sil:.3f}  (>0.5 = real clusters; ~0 = no structure)")
print(f"  adjusted Rand index vs minimap2 placement = {ari:.3f}  (1=perfect, 0=random)")

print("\n=== VERDICT ===")
if a_o - d_o < 0.05:
    print(f"  ORTHOGONAL-feature accuracy ({a_o:.3f}) ~ chance ({d_o:.3f}): NO BAM feature other than AS/PSV")
    print(f"  separates the copies. The tied reads (AS+PSV tied) carry NO recoverable copy signal — confirms")
    print(f"  the information-theoretic floor. EM/ML cannot help: there is nothing to learn.")
else:
    print(f"  ORTHOGONAL accuracy ({a_o:.3f}) > chance ({d_o:.3f}): a residual non-AS signal exists — a LEAD.")
print(f"  (AS-inclusive upper bound {a_all:.3f}: any gain over orthogonal is AS-driven = useless for tied.)")
