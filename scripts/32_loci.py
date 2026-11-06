#!/usr/bin/env python3
"""
iSTART: example loci for Figure 6, drawn with pyGenomeTracks.

Picks, from the final tables:
  nested   – named known nested genes (default NSMCE3, CHMP1B) from results/known_genes_wide.tsv
  short    – intronic start sites with a short DRB burst in HCT116 (class "short")
  product  – intronic start sites with productive transcription in HCT116 (class "productive")
  internal – internal promoters of host genes, productive in HCT116
Intronic sites and internal promoters must have a confirmed capped start in both cell lines and an active host
gene in both; they are ranked by the smaller of the two PRO-cap counts. The top --n of each kind are drawn.
Window: 10 kb before to 20 kb after the site in the direction of transcription (nested genes: gene +-10 kb).
Tracks: genes (canonical transcripts), PRO-cap both strands (K562, HCT116), RNA-seq on the site's strand
(K562, HCT116), nascent RNA 10 and 40 min after DRB (HCT116), H3K27ac and H3K4me3 peaks (both lines), phyloP.
"""
import argparse
import gzip
import importlib.util
import os
import subprocess

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
s = importlib.util.spec_from_file_location("kg", os.path.join(HERE, "30_known_genes.py"))
kg = importlib.util.module_from_spec(s)
s.loader.exec_module(kg)

TRACKS = dict(
    cap_K=("data/cap_plus.bw", "data/cap_minus.bw"),
    cap_H=("data/hct116/CTCFline_plus.bw", "data/hct116/CTCFline_minus.bw"),
    rna_K=("data/rnaseq/K562_total_plus.bw", "data/rnaseq/K562_total_minus.bw"),
    rna_H=("data/hct116/RNA_plus.bw", "data/hct116/RNA_minus.bw"),
    tt10=("data/hct116_drb/TT_10min.read_plus.bw", "data/hct116_drb/TT_10min.read_minus.bw"),
    tt40=("data/hct116_drb/TT_40min.read_plus.bw", "data/hct116_drb/TT_40min.read_minus.bw"),
    peaks=dict(K562=dict(H3K27ac="data/encode/histone/H3K27ac.bed.gz", H3K4me3="data/encode/histone/H3K4me3.bed.gz"),
               HCT116=dict(H3K27ac="data/hct116/H3K27ac.bed.gz", H3K4me3="data/hct116/H3K4me3.bed.gz")),
)


def near(df, chrom, strand, pos, w=200):
    x = df[(df.chrom == chrom) & (df.strand == strand)]
    x = x[(x.position - pos).abs() <= w]
    return x.iloc[(x.position - pos).abs().argsort()[:1]] if len(x) else x


def pick(a):
    S1 = pd.read_csv("results/final.supp_table1_intronic_starts.tsv", sep="\t", low_memory=False)
    for c in ("host_active", "capped_start_confirmed", "shared_with_other_cell_line"):
        S1[c] = S1[c].astype(str) == "True"
    K, H = S1[S1.cell_line == "K562"], S1[S1.cell_line == "HCT116"]
    picks = []
    for cls, kind in (("short", "short"), ("productive", "product")):
        x = H[(H.drb_class_10min == cls) & H.capped_start_confirmed & H.host_active & H.shared_with_other_cell_line]
        rows = []
        for r in x.itertuples():
            k = near(K, r.chrom, r.strand, r.position)
            if len(k) and bool(k.capped_start_confirmed.iloc[0]) and bool(k.host_active.iloc[0]):
                rows.append((min(r.procap_reads, k.procap_reads.iloc[0]), r))
        for sc, r in sorted(rows, key=lambda t: -t[0])[:a.n]:
            picks.append(dict(kind=kind, name=f"{kind}_{r.host_gene}_{r.chrom}_{r.position}", chrom=r.chrom, strand=r.strand,
                              pos=int(r.position), label=f"intronic start in {r.host_gene}", cap_min=sc, host=r.host_gene))
    # internal promoters: same site ids in both lines
    ch = pd.read_csv("results/hct116/classes_10.tsv", sep="\t", low_memory=False)
    dk = pd.read_csv("results/v2/itss.describe.tsv", sep="\t", low_memory=False)
    for d in (ch, dk):
        for c in ("host_active", "tss_supported"):
            d[c] = d[c].astype(str) == "True"
    ch = ch[(ch.site_type == "internal_promoter") & ch.host_active & ch.tss_supported & (ch.klasa == "produktywny")]
    dk = dk[(dk.site_type == "internal_promoter") & dk.host_active & dk.tss_supported][["site_id", "cap_sum"]]
    m = ch.merge(dk, on="site_id", suffixes=("_H", "_K"))
    m["sc"] = m[["cap_sum_H", "cap_sum_K"]].min(axis=1)
    hn = "host_gene_name" if "host_gene_name" in m else ("host_name" if "host_name" in m else None)
    for r in m.sort_values("sc", ascending=False).head(a.n).itertuples():
        host = getattr(r, hn) if hn else "host"
        picks.append(dict(kind="internal", name=f"internal_{host}_{r.chrom}_{int(r.pos)}", chrom=r.chrom, strand=r.strand,
                          pos=int(r.pos), label=f"internal promoter of {host}", cap_min=r.sc, host=host))
    kgw = pd.read_csv("results/known_genes_wide.tsv", sep="\t", low_memory=False)
    for g in a.nested:
        x = kgw[kgw.name == g]
        if not len(x):
            print(f"nested gene {g} not found")
            continue
        r = x.iloc[0]
        tss = [int(t) for t in str(r.tss).split(";")]
        picks.append(dict(kind="nested", name=f"nested_{g}", chrom=r.chrom, strand=r.strand, pos=tss[0],
                          label=f"{g} nested in {r.host}", start=int(r.start) - 10000, end=int(r.end) + 10000,
                          cap_min=min(r.cap_K562, r.cap_HCT116), host=r.host))
    for p in picks:
        if "start" not in p:
            p["start"], p["end"] = ((p["pos"] - 10000, p["pos"] + 20000) if p["strand"] == "+" else (p["pos"] - 20000, p["pos"] + 10000))
    return pd.DataFrame(picks)


def genes_bed(genes, tx, exons, chrom, s, e, path, host="", host_path=None):
    n_host = 0
    with open(path, "w") as f, open(host_path or os.devnull, "w") as fh:
        for tid, t in tx.items():
            g = genes.get(t["gene_id"])
            if g is None or g["chrom"] != chrom or t["end"] < s or t["start"] > e:
                continue
            if "Ensembl_canonical" not in t["tags"] or g["type"] not in ("protein_coding", "lncRNA"):
                continue
            ex = sorted(exons[tid])
            sizes = ",".join(str(b - a) for a, b in ex)
            starts = ",".join(str(a - t["start"]) for a, b in ex)
            line = f"{chrom}\t{t['start']}\t{t['end']}\t{g['name']}\t0\t{t['strand']}\t{t['start']}\t{t['end']}\t31,120,180\t{len(ex)}\t{sizes}\t{starts}\n"
            if g["name"] == host:
                fh.write(line)
                n_host += 1
            else:
                f.write(line)
        # pyGenomeTracks 3.9 applies all_labels_inside to every feature of a track except the first one; this unnamed
        # 1-nt white feature at the left edge of the window comes first, so that a gene label that would run past the
        # right edge is drawn to the left of its gene.
        s0 = max(0, s)
        f.write(f"{chrom}\t{s0}\t{s0 + 1}\t\t0\t+\t{s0}\t{s0 + 1}\t255,255,255\t1\t1\t0\n")
    return n_host


def smooth_bg(src, chrom, s, e, path, w=301, step=50, sign=1):
    """PRO-cap reads summed in a sliding window (display only), written as bedGraph; returns the maximum."""
    import pyBigWig
    bw = pyBigWig.open(src)
    L = bw.chroms(chrom) or 0
    s2, e2 = max(0, s), min(L, e)
    v = np.abs(np.nan_to_num(bw.values(chrom, s2, e2, numpy=True))) if e2 > s2 else np.zeros(0)
    c = np.convolve(v, np.ones(w), "same") if len(v) else v
    with open(path, "w") as f:
        for i in range(0, len(c), step):
            f.write(f"{chrom}\t{s2 + i}\t{min(s2 + i + step, e2)}\t{sign * c[i:i + step].max():.1f}\n")
    return float(c.max()) if len(c) else 0.0


def peaks_bed(src, chrom, s, e, path):
    with gzip.open(src, "rt") as fi, open(path, "w") as fo:
        for line in fi:
            f = line.split("\t")
            if f[0] == chrom and int(f[2]) > s and int(f[1]) < e:
                fo.write("\t".join(f[:3]) + "\n")


def ini(p, d, phylop, capmax, has_host=False):
    st = p["strand"]
    k = 0 if st == "+" else 1           # file of the site's strand (RNA strand)
    ko = 1 - k                          # nascent RNA reads lie on the opposite strand
    sname = "+" if st == "+" else "−"
    L = ["[x-axis]\nwhere = top\nfontsize = 9\n",
         (f"[host]\nfile = {d}/host.bed\ntitle = host gene: {p['host']} ({sname} strand)\nheight = 0.6\nstyle = flybase\n"
          f"gene_rows = 1\nlabels = false\ncolor = #2b6cb0\nborder_color = none\nfile_type = bed\n" if has_host else ""),
         f"[genes]\nfile = {d}/genes.bed\ntitle = other genes (canonical)\nheight = 1.8\nfontsize = 8\nstyle = flybase\ngene_rows = 3\nlabels = true\nall_labels_inside = true\ncolor = bed_rgb\nborder_color = none\nfile_type = bed\n",
         f"[site]\nfile = {d}/site.bed\ntitle = start site\nheight = 0.5\ncolor = #e53e3e\nlabels = false\nfile_type = bed\n",
         "[spacer]\nheight = 0.3\n"]
    for line, key in (("K562", "cap_K"), ("HCT116", "cap_H")):
        M = max(capmax[key], 1.0)
        L.append(f"[procap {line} +]\nfile = {d}/{key}_plus.bedgraph\ntitle = PRO-cap {line} (+ up, − down; reads/300 nt)\n"
                 f"height = 2.0\ncolor = #2b6cb0\nmin_value = {-M:.0f}\nmax_value = {M:.0f}\n")
        L.append(f"[procap {line} -]\nfile = {d}/{key}_minus.bedgraph\ncolor = #c05621\noverlay_previous = share-y\n")
        L.append(f"[spacer]\nheight = 0.3\n")
    for line, key in (("K562", "rna_K"), ("HCT116", "rna_H")):
        L.append(f"[rna {line}]\nfile = {os.path.abspath(TRACKS[key][k])}\ntitle = RNA-seq {line} ({sname})\nheight = 1.2\n"
                 f"color = #4a5568\nnumber_of_bins = 700\nmin_value = 0\n")
    for key, lab in (("tt10", "10 min"), ("tt40", "40 min")):
        L.append(f"[{key}]\nfile = {os.path.abspath(TRACKS[key][ko])}\ntitle = nascent RNA HCT116, {lab} after DRB ({sname})\n"
                 f"height = 1.2\ncolor = #2f855a\nnumber_of_bins = 700\nmin_value = 0\n")
    for line in ("K562", "HCT116"):
        for mark in ("H3K27ac", "H3K4me3"):
            L.append(f"[{mark} {line}]\nfile = {d}/{mark}_{line}.bed\ntitle = {mark} {line}\nheight = 0.4\n"
                     f"color = #C9D1DA\nborder_color = none\ndisplay = collapsed\nlabels = false\nfile_type = bed\n")
    if phylop:
        L.append(f"[phylop]\nfile = {phylop}\ntitle = phyloP (100 vertebrates)\nheight = 1.0\ncolor = #805ad5\n"
                 f"number_of_bins = 700\nmin_value = -2\n")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", default="ref/gencode.v47.annotation.gtf.gz")
    ap.add_argument("--phylop", help="local copy of hg38.phyloP100way.bw (optional)")
    ap.add_argument("--nested", nargs="+", default=["NSMCE3", "CHMP1B"])
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--out", default="results/figures/loci")
    ap.add_argument("--pgt", default="pyGenomeTracks")
    ap.add_argument("--only", nargs="*", help="draw only these loci (names from picked_loci.tsv)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    P = pick(a)
    P.to_csv(f"{a.out}/picked_loci.tsv", sep="\t", index=False)
    print(P[["kind", "label", "chrom", "strand", "pos", "cap_min"]].to_string(index=False))
    genes, tx, exons = kg.parse_gtf(a.gtf)
    for p in P.to_dict("records"):
        if a.only and p["name"] not in a.only:
            continue
        d = os.path.abspath(f"{a.out}/{p['name']}")
        os.makedirs(d, exist_ok=True)
        n_host = genes_bed(genes, tx, exons, p["chrom"], p["start"], p["end"], f"{d}/genes.bed",
                           host=p.get("host", ""), host_path=f"{d}/host.bed")
        for line, marks in TRACKS["peaks"].items():
            for mark, src in marks.items():
                peaks_bed(src, p["chrom"], p["start"], p["end"], f"{d}/{mark}_{line}.bed")
        capmax = {}
        for key in ("cap_K", "cap_H"):
            capmax[key] = max(smooth_bg(TRACKS[key][i], p["chrom"], p["start"], p["end"], f"{d}/{key}_{sname}.bedgraph",
                                        sign=(1 if i == 0 else -1)) for i, sname in ((0, "plus"), (1, "minus")))
        with open(f"{d}/site.bed", "w") as f:
            f.write(f"{p['chrom']}\t{max(0, p['pos'] - 150)}\t{p['pos'] + 150}\tsite\t0\t{p['strand']}\n")
        with open(f"{d}/tracks.ini", "w") as f:
            f.write(ini(p, d, a.phylop, capmax, has_host=n_host > 0))
        region = f"{p['chrom']}:{max(0, p['start'])}-{p['end']}"
        for ext in ("pdf", "png"):
            r = subprocess.run([a.pgt, "--tracks", f"{d}/tracks.ini", "--region", region, "--outFileName",
                                f"{a.out}/{p['name']}.{ext}", "--width", "18", "--dpi", "300"], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            if r.returncode != 0:
                raise SystemExit(f"pyGenomeTracks failed for {p['name']}:\n{r.stderr[-2000:]}")
        print(f"drawn: {a.out}/{p['name']}.pdf")


if __name__ == "__main__":
    main()
