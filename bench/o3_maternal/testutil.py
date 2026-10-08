"""Synthetic BAMs for the o3_maternal tests."""
import os

import pysam


def write_bam(path, refs, segs):
    """refs: {name: length}; segs: [dict(name, flag=0, ref, start, cigar='100M', mapq=60, AS=100, de=0.01, seq='A'*100)].
    Written coordinate-sorted and indexed."""
    tmp = path + ".unsorted"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": n, "LN": l} for n, l in refs.items()]}
    ids = {n: i for i, n in enumerate(refs)}
    with pysam.AlignmentFile(tmp, "wb", header=header) as f:
        for s in segs:
            a = pysam.AlignedSegment(f.header)
            a.query_name = s["name"]
            a.flag = s.get("flag", 0)
            unm = bool(a.flag & 4)
            a.reference_id = -1 if unm else ids[s["ref"]]
            a.reference_start = -1 if unm else s["start"]
            a.mapping_quality = s.get("mapq", 60)
            if not unm:
                a.cigarstring = s.get("cigar", "100M")
            a.query_sequence = s.get("seq", "A" * 100)
            if not unm:
                a.set_tag("AS", s.get("AS", 100))
                a.set_tag("de", s.get("de", 0.01))
            f.write(a)
    pysam.sort("-o", path, tmp)
    pysam.index(path)
    os.remove(tmp)
