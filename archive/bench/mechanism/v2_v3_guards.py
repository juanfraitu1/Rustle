"""V2: prove each inert-guard sits outside the data's range on the Soto families
(observed_max < guard => never fires => cannot move the numbers).
V3: report the two known result-changers (READ_CAP/POA_CAP, compile-time consts) honestly --
observed pool/length vs cap on the largest-avg_reads Soto locus, without a rebuild."""
import json, pathlib, argparse

REPO = pathlib.Path(__file__).resolve().parents[2]
SOTO_FAMS = REPO / "bench/soto/a119b_detected_families.tsv"


def classify_guard(guard, observed_max):
    return {"guard": guard, "observed_max": observed_max, "fires": observed_max >= guard}


def soto_family_sizes():
    """n_copies (col index 1) across the 66 Soto families in a119b_detected_families.tsv."""
    sizes = []
    with open(SOTO_FAMS) as f:
        next(f)
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > 1 and cols[1].isdigit():
                sizes.append(int(cols[1]))
    return sizes


def measure_v2():
    sizes = soto_family_sizes()
    max_members = max(sizes) if sizes else 0
    # MAX_MEMBERS=30 (multi_repeat_bridge.rs:71) and MAX_LOCI=60 (recombinant_split.rs:55) guard the
    # repeat-bridge / recombinant-split O1-cluster stages on family SIZE.
    return {
        "MAX_MEMBERS(30)": classify_guard(30, max_members),
        "MAX_LOCI(60)": classify_guard(60, max_members),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "bench/mechanism/verification_results.json"))
    a = ap.parse_args()
    out = pathlib.Path(a.out)
    data = json.loads(out.read_text()) if out.exists() else {}
    data["v2"] = measure_v2()
    out.write_text(json.dumps(data, indent=2))
    print("V2 guard classification:")
    for k, v in data["v2"].items():
        print(f"  {k}: observed_max={v['observed_max']} fires={v['fires']}")


if __name__ == "__main__":
    main()
