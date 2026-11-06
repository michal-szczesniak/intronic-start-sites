#!/usr/bin/env python3
"""
iSTART, description of the start sites (core version): chromatin and (optionally) conservation.

For each group of sites in the module A table (*.drb.tsv), computes the percentage of sites with a peak of each
histone mark within ±--window nt and assigns the sites to types (column typ):
  promoter-like     = H3K4me3 present
  enhancer-like     = no H3K4me3, H3K27ac or H3K4me1 present
  no marks          = none of the three
Reference: host start sites used in K562 (promoters). In addition, the DRB effect (gap_index)
by chromatin type for capped iTSS in active hosts.

Input:
  --drb       *.drb.tsv from 20_drb_gap.py
  --peaks     NAME=file.bed.gz (repeatable), e.g. H3K4me3=... H3K27ac=... H3K4me1=...
  --phylop    (optional) phyloP bigWig (file or URL); mean over ±100 nt
"""
import argparse
import bisect
import gzip
import sys
from collections import defaultdict

import numpy as np
import pandas as pd


def load_peaks(path):
    idx = defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track")) or not line.strip():
                continue
            f = line.split("\t")
            idx[f[0]].append((int(f[1]), int(f[2])))
    out = {}
    for c, v in idx.items():
        v.sort()
        out[c] = ([s for s, _ in v], [e for _, e in v], max(e - s for s, e in v))
    return out


def has_peak(pk, chrom, pos, w):
    if chrom not in pk:
        return False
    st, en, mx = pk[chrom]
    j0 = bisect.bisect_left(st, pos - w - mx)
    j1 = bisect.bisect_right(st, pos + w)
    return any(en[j] > pos - w for j in range(j0, j1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--drb", required=True)
    ap.add_argument("--peaks", nargs="+", required=True)
    ap.add_argument("--window", type=int, default=1000)
    ap.add_argument("--phylop")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    d = pd.read_csv(a.drb, sep="\t", low_memory=False)
    for col in ("host_active", "tss_supported", "testable", "eligible"):
        d[col] = d[col].astype(str) == "True"
    peaks = {}
    for spec in a.peaks:
        name, path = spec.split("=", 1)
        peaks[name] = load_peaks(path)
        print(f"[PEAKS] {name}: {sum(len(v[0]) for v in peaks[name].values())}", file=sys.stderr)
    names = list(peaks)

    # reference: host start sites used in K562
    hp = d[d.host_active].dropna(subset=["host_tss_used"]).drop_duplicates("host_gene_id")
    hp = pd.DataFrame(dict(site_id="host:" + hp.host_gene_id, site_type="host_promoter", chrom=hp.chrom,
                           strand=hp.strand, pos=hp.host_tss_used.astype(int), host_active=True,
                           tss_supported=True, testable=False, eligible=False))
    allsites = pd.concat([d, hp], ignore_index=True)
    for n in names:
        allsites[n] = [has_peak(peaks[n], r.chrom, int(r.pos), a.window) for r in allsites.itertuples()]
    k4me3 = allsites[names[0]] if "H3K4me3" not in names else allsites["H3K4me3"]
    enh = [n for n in names if n != "H3K4me3"]
    # chromatin type labels written to the table: promoter-like / enhancer-like / no marks
    allsites["typ"] = np.where(k4me3, "promotor-podobny",
                               np.where(allsites[enh].any(axis=1), "enhancer-podobny", "bez znaczników"))

    if a.phylop:
        import pyBigWig
        bw = pyBigWig.open(a.phylop)
        def ph(r):
            try:
                v = bw.values(r.chrom, int(r.pos) - 100, int(r.pos) + 101, numpy=True)
                return float(np.nanmean(v))
            except Exception:
                return np.nan
        sel = allsites.site_type.isin(["candidate", "internal_promoter", "host_promoter"]) | \
            ((allsites.site_type == "random_control") & (allsites.index % 3 == 0))
        allsites["phylop"] = np.nan
        allsites.loc[sel, "phylop"] = [ph(r) for r in allsites[sel].itertuples()]

    cand = allsites.site_type == "candidate"
    act = allsites.host_active
    sup = allsites.tss_supported
    groups = [
        ("host start sites (promoters, reference)", allsites.site_type == "host_promoter"),
        ("capped internal promoters", (allsites.site_type == "internal_promoter") & act & sup),
        ("capped iTSS, active host", cand & act & sup),
        ("capped iTSS, inactive host", cand & ~act & sup),
        ("uncapped iTSS", cand & ~sup),
        ("random controls", (allsites.site_type == "random_control") & act),
    ]
    log = [f"Chromatin marks within ±{a.window} nt (% of sites with a peak) and chromatin type:"]
    hdr = f"  {'group':<44}{'n':>6}" + "".join(f"{n:>9}" for n in names) + \
          f"{'promoter':>10}{'enhancer':>10}{'none':>8}" + (f"{'phyloP':>8}" if a.phylop else "")
    log.append(hdr)
    for gname, m in groups:
        g = allsites[m]
        if not len(g):
            continue
        row = f"  {gname:<44}{len(g):>6}" + "".join(f"{g[n].mean() * 100:>8.0f}%" for n in names)
        for t in ("promotor-podobny", "enhancer-podobny", "bez znaczników"):
            row += f"{(g.typ == t).mean() * 100:>9.0f}%"
        if a.phylop:
            row += f"{g.phylop.mean():>8.2f}"
        log.append(row)

    # DRB effect by chromatin type (capped iTSS, active host, testable)
    ctrl = allsites[(allsites.site_type == "random_control") & act & allsites.testable].gap_index
    if len(ctrl):
        thr = np.percentile(ctrl, 95)
        log.append("")
        log.append(f"DRB effect by chromatin type (capped iTSS, active host; control threshold {thr:.2f}):")
        t = allsites[cand & act & sup & allsites.testable]
        for typ, g in t.groupby("typ"):
            log.append(f"  {typ:<20} n={len(g):>4}   above threshold {(g.gap_index > thr).mean() * 100:5.1f}%   "
                       f"median gap_index {g.gap_index.median():5.2f}")
        log.append(f"  {'random controls':<20} n={len(ctrl):>4}   above threshold   5.0%   median gap_index {ctrl.median():5.2f}")

    allsites.to_csv(f"{a.out_prefix}.describe.tsv", sep="\t", index=False)
    txt = "\n".join(log)
    open(f"{a.out_prefix}.describe.txt", "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
