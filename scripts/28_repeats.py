#!/usr/bin/env python3
"""
iSTART: do intronic start sites arise from repeat elements (transposons)?

For each site: the RepeatMasker element (UCSC rmsk) covering the start (±--pad nt), its class,
family, name and orientation relative to the direction of transcription from the start
(orientacja: zgodna = same, przeciwna = opposite).
Comparison with random controls in the same genes; the most enriched elements;
starts in transposons vs the rest: conservation, stable RNA, DRB class (if given).
Usage: 28_repeats.py --rmsk rmsk.txt.gz --set NAME=describe.tsv[,phylop.tab[,expr.tsv[,classes.tsv]]] ...
"""
import argparse
import bisect
import gzip
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def load_rmsk(path):
    idx = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            idx[f[5]].append((int(f[6]), int(f[7]), f[9], f[10], f[11], f[12]))
    out = {}
    for c, v in idx.items():
        v.sort()
        out[c] = (v, [x[0] for x in v], max(x[1] - x[0] for x in v))
    return out


def hit(rm, chrom, pos, pad):
    if chrom not in rm:
        return None
    v, st, mx = rm[chrom]
    best = None
    for j in range(bisect.bisect_left(st, pos - pad - mx), bisect.bisect_right(st, pos + pad)):
        s, e, strand, name, cls, fam = v[j]
        if e > pos - pad and s <= pos + pad:
            if best is None or (e - s) > (best[1] - best[0]):
                best = v[j]
    return best


def simple(cls, fam):
    if cls in ("SINE",):
        return "SINE/" + ("Alu" if fam == "Alu" else "inne")
    if cls == "LINE":
        return "LINE/" + ("L1" if fam == "L1" else "inne")
    if cls == "LTR":
        return "LTR/" + ("ERV" if fam.startswith("ERV") else fam)
    if cls == "DNA":
        return "DNA"
    if cls in ("Simple_repeat", "Low_complexity"):
        return "proste/niska złożoność"
    return "inne"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rmsk", required=True)
    ap.add_argument("--set", nargs="+", required=True)
    ap.add_argument("--pad", type=int, default=0)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    rm = load_rmsk(a.rmsk)
    print(f"RepeatMasker elements: {sum(len(v[0]) for v in rm.values())}")

    for spec in a.set:
        name, files = spec.split("=", 1)
        files = files.split(",")
        d = pd.read_csv(files[0], sep="\t", low_memory=False).reset_index()
        for c in ("host_active", "tss_supported"):
            d[c] = d[c].astype(str) == "True"
        if len(files) > 1 and files[1]:
            t = pd.read_csv(files[1], sep="\t", header=None, names=["n", "size", "covered", "sum", "mean0", "mean"])
            t["index"] = t.n.str[1:].astype(int)
            d = d.merge(t[["index", "mean"]].rename(columns={"mean": "phylop"}), on="index", how="left")
        if len(files) > 2 and files[2]:
            e = pd.read_csv(files[2], sep="\t", low_memory=False)
            ea = e.host_active.astype(str) == "True"
            thr = np.nanpercentile(e[(e.site_type == "random_control") & ea].rna_step, 95)
            d = d.merge(e[["site_id", "rna_step"]], on="site_id", how="left")
            d["stabilne_RNA"] = d.rna_step > thr
        if len(files) > 3 and files[3]:
            k = pd.read_csv(files[3], sep="\t", low_memory=False)[["site_id", "klasa"]]
            d = d.merge(k, on="site_id", how="left")
        rows = []
        for r in d.itertuples():
            h = hit(rm, r.chrom, int(r.pos), a.pad)
            if h is None:
                rows.append(("brak", "brak", "", ""))
            else:
                s, e_, strand, nm, cls, fam = h
                orient = "zgodna" if strand == r.strand else "przeciwna"
                rows.append((simple(cls, fam), f"{cls}/{fam}", nm, orient))
        d["grupa_rep"], d["rodzina"], d["element"], d["orientacja"] = zip(*rows)
        d.to_csv(f"{a.out_prefix}.{name}.tsv", sep="\t", index=False)

        it = d[(d.site_type == "candidate") & d.host_active & d.tss_supported]
        ct = d[(d.site_type == "random_control") & d.host_active]
        hp = d[d.site_type == "host_promoter"]
        print(f"\n===== {name}: intronic start sites (capped, active host) n={len(it)}, "
              f"random controls n={len(ct)}, host promoters n={len(hp)} =====")
        print(f"  {'element class':<26}{'starts':>9}{'controls':>10}{'promoters':>11}{'fold':>10}{'p (Fisher)':>12}")
        for g in sorted(set(d.grupa_rep)):
            a1, b1 = int((it.grupa_rep == g).sum()), len(it)
            a2, b2 = int((ct.grupa_rep == g).sum()), len(ct)
            p = fisher_exact([[a1, b1 - a1], [a2, b2 - a2]])[1]
            fold = (a1 / b1) / max(a2 / b2, 1e-9)
            print(f"  {g:<26}{a1 / b1 * 100:>8.1f}%{a2 / b2 * 100:>9.1f}%{(hp.grupa_rep == g).mean() * 100:>10.1f}%{fold:>10.2f}{p:>12.1e}")
        print("  element orientation relative to transcription from the start (intronic start sites):")
        for g in ("LTR/ERV", "LINE/L1", "SINE/Alu"):
            x = it[it.grupa_rep == g]
            if len(x):
                print(f"    {g:<10} n={len(x):>4}  same orientation {(x.orientacja == 'zgodna').mean() * 100:.0f}%")
        # most enriched elements (min. 8 starts)
        ci, cc = it.element.value_counts(), ct.element.value_counts()
        tab = []
        for el, n1 in ci.items():
            if el == "" or n1 < 8:
                continue
            n2 = int(cc.get(el, 0))
            p = fisher_exact([[n1, len(it) - n1], [n2, len(ct) - n2]], alternative="greater")[1]
            tab.append((el, it[it.element == el].rodzina.iloc[0], n1, n2, (n1 / len(it)) / max((n2 + 0.5) / len(ct), 1e-9), p))
        tab.sort(key=lambda x: x[5])
        print("  most enriched elements (starts / controls, fold, p):")
        for el, fam, n1, n2, f, p in tab[:15]:
            print(f"    {el:<16}{fam:<22}{n1:>5}{n2:>7}{f:>8.1f}{p:>10.1e}")
        # starts in transposons vs the rest
        it = it.copy()
        it["TE"] = ~it.grupa_rep.isin(["brak", "proste/niska złożoność"])
        print(f"  starts in transposons: {it.TE.mean() * 100:.0f}% (controls: {(~ct.grupa_rep.isin(['brak', 'proste/niska złożoność'])).mean() * 100:.0f}%)")
        for col, lab in (("phylop", "phyloP > 1"), ("stabilne_RNA", "stable RNA")):
            if col in it:
                v = (it[col] > 1) if col == "phylop" else it[col].astype(bool)
                print(f"    {lab}: in transposons {v[it.TE].mean() * 100:.0f}%, outside them {v[~it.TE].mean() * 100:.0f}%")
        if "klasa" in it:
            x = it[it.klasa.isin(["produktywny", "krótki", "bez odpowiedzi"])]
            print("    DRB class (in transposons / outside them): " + ", ".join(
                f"{k} {(x[x.TE].klasa == k).mean() * 100:.0f}% / {(x[~x.TE].klasa == k).mean() * 100:.0f}%"
                for k in ("produktywny", "krótki", "bez odpowiedzi")))


if __name__ == "__main__":
    main()
