#!/usr/bin/env python3
"""Downloads the RNA-seq bigWigs (unique reads, GRCh38) of an ENCODE experiment and its control and merges replicates.
For each replicate and strand takes the file flagged as default, otherwise the newest.
Usage: encode_rnaseq_bw.py NAME KO_ACCESSION CONTROL_ACCESSION DIR"""
import os
import subprocess
import sys

import requests

B = "https://www.encodeproject.org"
H = {"Accept": "application/json"}
name, ko, ctrl, out = sys.argv[1:5]
os.makedirs(out, exist_ok=True)
here = os.path.dirname(os.path.abspath(__file__))


def pick(acc):
    p = [("type", "File"), ("dataset", f"/experiments/{acc}/"), ("assembly", "GRCh38"), ("file_format", "bigWig"),
         ("status", "released"), ("format", "json"), ("limit", "all"), ("field", "accession"), ("field", "output_type"),
         ("field", "biological_replicates"), ("field", "preferred_default"), ("field", "href"), ("field", "date_created")]
    fs = [f for f in requests.get(B + "/search/", params=p, headers=H, timeout=60).json().get("@graph", [])
          if "signal of unique reads" in f.get("output_type", "")]
    fs.sort(key=lambda f: f.get("date_created", ""), reverse=True)
    fs.sort(key=lambda f: 0 if f.get("preferred_default") else 1)
    best = {}
    for f in fs:
        strand = "plus" if f["output_type"].startswith("plus") else "minus"
        k = (tuple(f.get("biological_replicates", [])), strand)
        best.setdefault(k, f)
    return best


for label, acc in (("KO", ko), ("CTRL", ctrl)):
    best = pick(acc)
    for strand in ("plus", "minus"):
        paths = []
        for (reps, st), f in sorted(best.items()):
            if st != strand:
                continue
            path = os.path.join(out, f["accession"] + ".bigWig")
            if not os.path.exists(path):
                subprocess.run(["wget", "-nv", "-O", path, B + f["href"]], check=True)
            paths.append(path)
            print(f"{name} {label} {acc} {strand} rep{list(reps)} -> {f['accession']}", flush=True)
        subprocess.run([sys.executable, os.path.join(here, "merge_bw.py"),
                        os.path.join(out, f"{name}_{label}_{strand}.bw")] + paths, check=True)
