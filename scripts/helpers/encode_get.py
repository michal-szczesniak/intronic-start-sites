#!/usr/bin/env python3
"""
Download processed ENCODE files (GRCh38) for the second cell line.

Usage:
  encode_get.py signal  NAME ACCESSION DIR [--to-counts]      plus/minus strand bigWigs (unique reads),
                                                              replicates merged -> NAME_plus.bw, NAME_minus.bw
  encode_get.py peaks   NAME ACCESSION DIR                    PINTS peaks -> NAME_bi.bed.gz, NAME_uni.bed.gz
  encode_get.py histone NAME ACCESSION DIR                    "replicated peaks" -> NAME.bed.gz
  encode_get.py bam     NAME ACCESSION DIR                    "alignments" (all replicates) -> NAME_repN.bam
When a file exists in several versions, takes the one flagged as default, otherwise the newest.
--to-counts: divides the signal by the smallest non-zero value in the file (normalised signal -> counts).
"""
import os
import subprocess
import sys

import numpy as np
import requests

B = "https://www.encodeproject.org"
H = {"Accept": "application/json"}
HERE = os.path.dirname(os.path.abspath(__file__))


def files(acc, fmt):
    p = [("type", "File"), ("dataset", f"/experiments/{acc}/"), ("assembly", "GRCh38"), ("file_format", fmt),
         ("status", "released"), ("format", "json"), ("limit", "all"), ("field", "accession"), ("field", "output_type"),
         ("field", "biological_replicates"), ("field", "preferred_default"), ("field", "href"), ("field", "date_created")]
    fs = requests.get(B + "/search/", params=p, headers=H, timeout=60).json().get("@graph", [])
    fs.sort(key=lambda f: f.get("date_created", ""), reverse=True)
    fs.sort(key=lambda f: 0 if f.get("preferred_default") else 1)
    return fs


def get(f, path):
    if not os.path.exists(path):
        subprocess.run(["wget", "-nv", "-O", path, B + f["href"]], check=True)
    print(f"  {f['accession']} ({f['output_type']}, rep {f.get('biological_replicates')}) -> {path}", flush=True)


def to_counts(path):
    import pyBigWig
    bw = pyBigWig.open(path)
    mn = np.inf
    for c in list(bw.chroms())[:5]:
        iv = bw.intervals(c)
        if iv:
            v = np.abs(np.array([x[2] for x in iv]))
            v = v[v > 0]
            if len(v):
                mn = min(mn, v.min())
    out = path.replace(".bigWig", ".counts.bw")
    o = pyBigWig.open(out, "w")
    o.addHeader(list(bw.chroms().items()))
    for c in bw.chroms():
        iv = bw.intervals(c)
        if iv:
            a = np.array(iv, dtype=float)
            o.addEntries([c] * len(a), a[:, 0].astype(int).tolist(), ends=a[:, 1].astype(int).tolist(),
                         values=(np.abs(a[:, 2]) / mn).tolist())
    o.close()
    print(f"  {path}: smallest non-zero value {mn:.4g} -> counts in {out}", flush=True)
    return out


def main():
    kind, name, acc, out = sys.argv[1:5]
    counts = "--to-counts" in sys.argv
    os.makedirs(out, exist_ok=True)
    print(f"== {kind} {name} {acc} ==", flush=True)
    if kind == "signal":
        best = {}
        for f in files(acc, "bigWig"):
            if "signal of unique reads" not in f["output_type"]:
                continue
            st = "plus" if f["output_type"].startswith("plus") else "minus"
            best.setdefault((tuple(f.get("biological_replicates", [])), st), f)
        for st in ("plus", "minus"):
            paths = []
            for (reps, s), f in sorted(best.items()):
                if s != st:
                    continue
                p = os.path.join(out, f["accession"] + ".bigWig")
                get(f, p)
                paths.append(to_counts(p) if counts else p)
            subprocess.run([sys.executable, os.path.join(HERE, "merge_bw.py"), os.path.join(out, f"{name}_{st}.bw")] + paths,
                           check=True)
    elif kind == "peaks":
        for ot, suf in (("bidirectional peaks", "bi"), ("unidirectional peaks", "uni")):
            fs = [f for f in files(acc, "bed") if f["output_type"] == ot]
            if fs:
                get(fs[0], os.path.join(out, f"{name}_{suf}.bed.gz"))
            else:
                print(f"  no files: {ot}")
    elif kind == "histone":
        fs = [f for f in files(acc, "bed") if f["output_type"] == "replicated peaks"]
        get(fs[0], os.path.join(out, f"{name}.bed.gz"))
    elif kind == "bam":
        fs = [f for f in files(acc, "bam") if f["output_type"] == "alignments"]
        seen = set()
        for f in fs:
            r = tuple(f.get("biological_replicates", []))
            if r in seen:
                continue
            seen.add(r)
            get(f, os.path.join(out, f"{name}_rep{'_'.join(map(str, r))}.bam"))


if __name__ == "__main__":
    main()
