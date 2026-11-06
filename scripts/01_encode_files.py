#!/usr/bin/env python3
"""
iSTART: list (and optionally download) the files of ENCODE experiments through the REST API.

Examples:
  # 1. see which files are available (always do this first):
  python scripts/01_encode_files.py ENCSR261KBX ENCSR114HGS ENCSR220XSM
  # 2. download selected files (filters are combined with AND; --output-type is a case-insensitive name fragment):
  python scripts/01_encode_files.py ENCSR261KBX --download --assembly GRCh38 --format bigWig \
         --output-type "strand signal" --outdir data/procap

List columns: experiment, file accession, format, output type, assembly, replicates, size, URL.
"""
import argparse
import os
import sys

import requests

BASE = "https://www.encodeproject.org"
FIELDS = ["accession", "file_format", "file_type", "output_type", "assembly", "href",
          "biological_replicates", "technical_replicates", "file_size", "status"]


def list_files(exp):
    params = [("type", "File"), ("dataset", f"/experiments/{exp}/"), ("status", "released"),
              ("format", "json"), ("limit", "all")] + [("field", f) for f in FIELDS]
    r = requests.get(f"{BASE}/search/", params=params, headers={"Accept": "application/json"}, timeout=60)
    if r.status_code == 404:
        return []
    r.raise_for_status()
    return r.json().get("@graph", [])


def keep(f, a):
    if a.assembly and f.get("assembly") != a.assembly:
        return False
    if a.format and f.get("file_format") != a.format:
        return False
    if a.output_type and a.output_type.lower() not in f.get("output_type", "").lower():
        return False
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("experiments", nargs="+")
    ap.add_argument("--assembly")
    ap.add_argument("--format")
    ap.add_argument("--output-type")
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--outdir", default="data/encode")
    a = ap.parse_args()

    print("\t".join(["experiment", "accession", "format", "output_type", "assembly",
                     "bio_rep", "tech_rep", "size_MB", "url"]))
    todo = []
    for exp in a.experiments:
        files = list_files(exp)
        if not files:
            print(f"[WARNING] no files for {exp}", file=sys.stderr)
        for f in sorted(files, key=lambda x: (x.get("file_format", ""), x.get("output_type", ""))):
            if not keep(f, a):
                continue
            url = BASE + f.get("href", "")
            print("\t".join(str(x) for x in [
                exp, f.get("accession"), f.get("file_format"), f.get("output_type"), f.get("assembly", "NA"),
                ",".join(map(str, f.get("biological_replicates", []))),
                ",".join(map(str, f.get("technical_replicates", []))),
                round((f.get("file_size") or 0) / 1e6, 1), url]))
            todo.append((exp, f, url))

    if a.download:
        for exp, f, url in todo:
            d = os.path.join(a.outdir, exp)
            os.makedirs(d, exist_ok=True)
            path = os.path.join(d, os.path.basename(url))
            if os.path.exists(path):
                continue
            print(f"[DOWNLOAD] {url}", file=sys.stderr)
            with requests.get(url, stream=True, timeout=600) as r:
                r.raise_for_status()
                with open(path + ".part", "wb") as fo:
                    for chunk in r.iter_content(1 << 20):
                        fo.write(chunk)
            os.rename(path + ".part", path)


if __name__ == "__main__":
    main()
