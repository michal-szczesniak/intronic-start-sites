#!/usr/bin/env python3
"""
List the samples of a GEO series with their SRA runs (NCBI runinfo); samples are selected by a title pattern.

Usage: geo_runs.py GSE148844 output.tsv [--include REGEX] [--exclude REGEX]
By default selects wild-type samples (title contains wt/wild) and excludes KO, rescue, siRNA and ELOF1/CSB.
Output: one row per run: gsm, title, selected (tak/nie = yes/no), time_min, rep, run, layout, source, read count.
The FASTQ files of the selected runs are then downloaded from ENA (jobs/hct116_drb.sh).
"""
import argparse
import re
import time

import requests


def ncbi_runs(srx):
    """Runs of an SRA experiment from NCBI (runinfo): [(run, layout, 'SRA', read_count)]."""
    import csv
    import io
    text = ""
    for _ in range(3):
        try:
            r = requests.get("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo", params={"acc": srx}, timeout=120)
            if r.status_code == 200 and r.text.startswith("Run,"):
                text = r.text
                break
            ids = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                               params={"db": "sra", "term": srx, "retmode": "json"}, timeout=120).json()["esearchresult"]["idlist"]
            r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                             params={"db": "sra", "id": ",".join(ids), "rettype": "runinfo", "retmode": "text"}, timeout=120)
            if r.status_code == 200 and r.text.lstrip().startswith("Run,"):
                text = r.text.lstrip()
                break
        except Exception:
            pass
        time.sleep(5)
    out = []
    for row in csv.DictReader(io.StringIO(text)):
        if row.get("Run", "").startswith(("SRR", "ERR", "DRR")) and row.get("Experiment", srx) == srx:
            out.append((row["Run"], row.get("LibraryLayout", "?"), "SRA", row.get("spots", "?")))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gse")
    ap.add_argument("out")
    ap.add_argument("--include", default=r"(^|[^a-z])(wt|wild)")
    ap.add_argument("--exclude", default=r"(ko\b|_ko|-/-|knock|rescue|re-?exp|sirna|elof1|csb|ercc|uvssa)")
    a = ap.parse_args()
    inc, exc = re.compile(a.include, re.I), re.compile(a.exclude, re.I)

    txt = requests.get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
                       params={"acc": a.gse, "targ": "gsm", "form": "text", "view": "brief"}, timeout=120).text
    samples, cur = [], None
    for line in txt.splitlines():
        if line.startswith("^SAMPLE"):
            cur = {"gsm": line.split("=", 1)[1].strip(), "title": "", "srx": ""}
            samples.append(cur)
        elif cur is not None and line.startswith("!Sample_title"):
            cur["title"] = line.split("=", 1)[1].strip()
        elif cur is not None and line.startswith("!Sample_relation") and "SRA" in line:
            m = re.search(r"(SRX\d+)", line)
            if m:
                cur["srx"] = m.group(1)
    print(f"{a.gse}: samples {len(samples)}")

    rows = []
    any_inc = any(inc.search(s["title"]) and not exc.search(s["title"]) for s in samples)
    if not any_inc:
        print("  (no title matches the --include pattern: selecting all samples that are not excluded)")
    for s in samples:
        t = s["title"]
        sel = (bool(inc.search(t)) or not any_inc) and not exc.search(t)
        m = re.search(r"(\d+)\s*(min|m)(?![a-z])", t, re.I)
        tmin = int(m.group(1)) if m else (0 if re.search(r"(no.?release|0h|unreleased)", t, re.I) else -1)
        mr = re.search(r"(rep|replicate|_r)\s*_?(\d)", t, re.I)
        rep = mr.group(2) if mr else "?"
        runs = ncbi_runs(s["srx"]) if (s["srx"] and sel) else []
        # column "selected" holds "tak" (yes) / "nie" (no); jobs/hct116_drb.sh tests $3=="tak"
        for run, layout, url, nreads in runs or [("NA", "?", "", "0")]:
            rows.append([s["gsm"], t, "tak" if sel else "nie", str(tmin), rep, run, layout, url, nreads])
        print(f"  {'*' if sel else ' '} {s['gsm']}  {t:<60} time={tmin:>3}  rep={rep}  runs={len(runs)}")
    with open(a.out, "w") as fo:
        fo.write("gsm\ttitle\tselected\ttime_min\trep\trun\tlayout\turl_r1\treads\n")
        for r in rows:
            fo.write("\t".join(r) + "\n")
    n = sum(1 for r in rows if r[2] == "tak")
    print(f"selected runs: {n}; table: {a.out}")


if __name__ == "__main__":
    main()
