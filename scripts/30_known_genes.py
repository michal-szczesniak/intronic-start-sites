#!/usr/bin/env python3
"""
iSTART: does the start-site test find known genes that lie inside introns? (positive control)

Nested genes: GENCODE protein-coding and lncRNA genes whose whole span lies inside one intron of the
canonical transcript (tag Ensembl_canonical) of another protein-coding gene on the SAME strand; hosts with
any readthrough transcript are excluded. Well-established: protein-coding with MANE_Select, or lncRNA with
gene level 1-2 and a transcript of support level (TSL) 1-2. Broadly expressed (optional, GTEx median TPM
table): TPM >= 1 in >= 80% of tissues.
Stand-alone genes (reference): protein-coding and lncRNA genes not overlapping any other protein-coding gene.

Per cell line (--line, repeatable): "NAME|cap_plus,cap_minus|uncap_plus,uncap_minus|rna_plus,rna_minus|sites.tsv"
(uncap and sites may be empty). For each gene:
  expression = mean RNA-seq coverage over its exons minus mean coverage of its intronic surroundings
               (2-5 kb on both sides, same strand), so host pre-mRNA does not count as expression
  detected   = the cap test passes near any annotated start (offsets -300..+300 nt, step 50):
               >= 5 capped reads in +-100 nt, Poisson p < 1e-3 against flanks 500-2000 nt,
               capped/uncapped ratio >= 2 when uncapped data are given
Summary: detection by expression decile; observed vs expected nested genes detected (expected from stand-alone
genes of the same type and decile); well-established and broadly expressed nested genes; overlap with the de novo
intronic start sites; examples detected in all cell lines.
"""
import argparse
import bisect
import gzip
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import poisson

ATTR = re.compile(r'(\w+) (?:"([^"]*)"|([^;\s"]+))')


def parse_gtf(path):
    genes, tx, exons = {}, {}, defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("gene", "transcript", "exon"):
                continue
            a, tags = {}, []
            for k, v1, v2 in ATTR.findall(f[8]):
                v = v1 or v2
                if k == "tag":
                    tags.append(v)
                else:
                    a[k] = v
            gid = a["gene_id"].split(".")[0]
            if f[2] == "gene":
                genes[gid] = dict(gene_id=gid, name=a.get("gene_name", gid), type=a.get("gene_type", ""), chrom=f[0],
                                  start=int(f[3]) - 1, end=int(f[4]), strand=f[6], level=int(a.get("level", 3)))
            elif f[2] == "transcript":
                tx[a["transcript_id"]] = dict(gene_id=gid, start=int(f[3]) - 1, end=int(f[4]), strand=f[6], tags=set(tags),
                                              tsl=a.get("transcript_support_level", "NA"))
            else:
                exons[a["transcript_id"]].append((int(f[3]) - 1, int(f[4])))
    return genes, tx, exons


class Cov:
    def __init__(self, spec):
        p, m = spec.split(",")
        self.bw = {"+": pyBigWig.open(p), "-": pyBigWig.open(m)}

    def vals(self, chrom, s, e, strand):
        bw = self.bw[strand]
        L = bw.chroms(chrom)
        out = np.zeros(max(e - s, 0))
        if L is None or e <= 0:
            return out
        s2, e2 = max(0, s), min(L, e)
        if e2 > s2:
            v = bw.values(chrom, s2, e2, numpy=True)
            out[s2 - s:e2 - s] = np.abs(np.nan_to_num(v))
        return out

    def mean(self, chrom, s, e, strand):
        v = self.vals(chrom, s, e, strand)
        return float(v.mean()) if len(v) else 0.0


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--gtex", help="GTEx gene median TPM table (.gct.gz), optional")
    ap.add_argument("--line", action="append", required=True)
    ap.add_argument("--orient", choices=["same", "both"], default="same",
                    help="nested genes on the host's strand only, or on either strand")
    ap.add_argument("--introns", choices=["canonical", "any"], default="canonical",
                    help="introns of the host's canonical transcript only, or of any of its transcripts")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    genes, tx, exons = parse_gtf(a.gtf)
    by_gene = defaultdict(list)
    for tid, t in tx.items():
        by_gene[t["gene_id"]].append(tid)

    # --- hosts: protein-coding, no readthrough transcript, canonical introns ---
    host_introns = defaultdict(list)  # (chrom, strand) -> (intron_s, intron_e, host_id)
    host_exons_all = defaultdict(list)  # host_id -> all exons of all transcripts
    for gid, g in genes.items():
        if g["type"] != "protein_coding":
            continue
        tids = by_gene[gid]
        if any("readthrough_transcript" in tx[t]["tags"] for t in tids):
            continue
        for t in tids:
            host_exons_all[gid] += exons[t]
        can = [t for t in tids if "Ensembl_canonical" in tx[t]["tags"]]
        use = can[:1] if a.introns == "canonical" else tids
        seen_i = set()
        for t in use:
            ex = sorted(exons[t])
            for (s1, e1), (s2, e2) in zip(ex[:-1], ex[1:]):
                if s2 > e1 and (e1, s2) not in seen_i:
                    seen_i.add((e1, s2))
                    host_introns[g["chrom"]].append((e1, s2, gid, g["strand"]))
    def index(d):
        out = {}
        for k, v in d.items():
            v.sort()
            out[k] = (v, [x[0] for x in v], max(x[1] - x[0] for x in v))
        return out

    HI = index(host_introns)
    pc_spans = defaultdict(list)
    for gid, g in genes.items():
        if g["type"] == "protein_coding":
            pc_spans[g["chrom"]].append((g["start"], g["end"], gid))
    PS = index(pc_spans)

    def overlapping(idx, key, s, e):
        if key not in idx:
            return []
        v, st, mx = idx[key]
        lo = bisect.bisect_left(st, s - mx)
        hi = bisect.bisect_left(st, e)
        return [x for x in v[lo:hi] if x[1] > s]

    def tss_list(gid):
        out = set()
        for t in by_gene[gid]:
            out.add(tx[t]["start"] if tx[t]["strand"] == "+" else tx[t]["end"] - 1)
        return sorted(out)

    rows = []
    for gid, g in genes.items():
        if g["type"] not in ("protein_coding", "lncRNA"):
            continue
        host, orient = None, ""
        for s, e, hid, hstr in overlapping(HI, g["chrom"], g["start"], g["end"]):
            if hid == gid or not (s <= g["start"] and g["end"] <= e):
                continue
            if hstr != g["strand"] and a.orient == "same":
                continue
            if host is None or (orient == "anti" and hstr == g["strand"]):
                host, orient = hid, ("same" if hstr == g["strand"] else "anti")
        tss = tss_list(gid)
        if host is not None:
            hx = host_exons_all[host]
            if any(min(abs(p - hs), abs(p - he)) < 500 or hs <= p < he for p in tss for hs, he in hx):
                continue  # start too close to a host exon to be judged
            kind = "nested"
        else:
            overl = any(hid != gid for s, e, hid in overlapping(PS, g["chrom"], g["start"], g["end"]))
            if overl:
                continue
            kind = "standalone"
        tids = by_gene[gid]
        if g["type"] == "protein_coding":
            well = any("MANE_Select" in tx[t]["tags"] for t in tids)
        else:
            well = g["level"] <= 2 and any(tx[t]["tsl"] in ("1", "2") for t in tids)
        ex = merge([e for t in tids for e in exons[t]])
        rows.append(dict(gene_id=gid, name=g["name"], type=g["type"], kind=kind, orient=orient if host else "", chrom=g["chrom"], strand=g["strand"],
                         start=g["start"], end=g["end"], host=genes[host]["name"] if host else "", well_established=well,
                         tss=";".join(map(str, tss)), exons=ex))
    G = pd.DataFrame(rows)
    print(f"genes: nested {int((G.kind == 'nested').sum())} (protein-coding {int(((G.kind == 'nested') & (G.type == 'protein_coding')).sum())}, "
          f"lncRNA {int(((G.kind == 'nested') & (G.type == 'lncRNA')).sum())}); stand-alone {int((G.kind == 'standalone').sum())}")
    for o in ("same", "anti"):
        x = G[(G.kind == "nested") & (G.orient == o)]
        if len(x):
            print(f"  nested on the {'host strand' if o == 'same' else 'opposite strand'}: protein-coding "
                  f"{int((x.type == 'protein_coding').sum())}, lncRNA {int((x.type == 'lncRNA').sum())}")

    if a.gtex:
        t = pd.read_csv(a.gtex, sep="\t", skiprows=2)
        t["gid"] = t.Name.str.split(".").str[0]
        tis = t.columns[2:-1]
        t["broad"] = (t[tis] >= 1).mean(axis=1) >= 0.8
        G = G.merge(t[["gid", "broad"]].rename(columns={"gid": "gene_id"}), on="gene_id", how="left")
        G["broad"] = G.broad.fillna(False).astype(bool)
    else:
        G["broad"] = np.nan

    lines = []
    for spec in a.line:
        name, capf, uncf, rnaf, sitesf = (spec.split("|") + ["", "", "", ""])[:5]
        cap, rna = Cov(capf), Cov(rnaf)
        unc = Cov(uncf) if uncf else None
        cap_tot = sum(abs(cap.bw[s].header()["sumData"]) for s in "+-")
        unc_tot = sum(abs(unc.bw[s].header()["sumData"]) for s in "+-") if unc else None
        itss = None
        if sitesf:
            d = pd.read_csv(sitesf, sep="\t", low_memory=False)
            d = d[(d.site_type == "candidate")]
            itss = defaultdict(list)
            for r in d.itertuples():
                itss[(r.chrom, r.strand)].append(int(r.pos))
            itss = {k: np.sort(v) for k, v in itss.items()}
        expr, det, capmax, in_itss, exm, flm = [], [], [], [], [], []
        for r in G.itertuples():
            ex = r.exons
            L = sum(e - s for s, e in ex)
            ex_mean = sum(rna.mean(r.chrom, s, e, r.strand) * (e - s) for s, e in ex) / max(L, 1)
            fl = (rna.mean(r.chrom, r.start - 5000, r.start - 2000, r.strand) + rna.mean(r.chrom, r.end + 2000, r.end + 5000, r.strand)) / 2
            expr.append(max(ex_mean - fl, 0.0))
            exm.append(ex_mean)
            flm.append(fl)
            best, ok = 0.0, False
            for p in map(int, r.tss.split(";")):
                v = cap.vals(r.chrom, p - 2400, p + 2401, r.strand)
                u = unc.vals(r.chrom, p - 400, p + 401, r.strand) if unc else None
                for off in range(-300, 301, 50):
                    c = 2400 + off
                    s = v[c - 100:c + 101].sum()
                    best = max(best, s)
                    if s < 5 or ok:
                        continue
                    flank = (v[c - 2000:c - 500].mean() + v[c + 500:c + 2000].mean()) / 2 if c >= 2000 else v[c + 500:c + 2000].mean()
                    if poisson.sf(s - 1, max(flank, 1e-3) * 201) >= 1e-3:
                        continue
                    if unc is not None:
                        uu = u[400 + off - 100:400 + off + 101].sum()
                        if ((s + 1) / cap_tot) / ((uu + 1) / unc_tot) < 2:
                            continue
                    ok = True
            det.append(ok)
            capmax.append(best)
            hit = False
            if itss is not None and r.kind == "nested":
                v = itss.get((r.chrom, r.strand))
                if v is not None:
                    for p in map(int, r.tss.split(";")):
                        j = np.searchsorted(v, p)
                        if any(0 <= jj < len(v) and abs(v[jj] - p) <= 500 for jj in (j - 1, j)):
                            hit = True
            in_itss.append(hit)
        G[f"expr_{name}"], G[f"det_{name}"], G[f"cap_{name}"], G[f"itss_{name}"] = expr, det, capmax, in_itss
        G[f"exon_{name}"], G[f"flank_{name}"] = exm, flm
        lines.append(name)

    for name in lines:
        e, dcol = f"expr_{name}", f"det_{name}"
        print(f"\n===== {name} =====")
        ref = G[(G.kind == "standalone") & (G[e] > 0)]
        edges = np.quantile(ref[e], [0.2, 0.4, 0.6, 0.8])
        G[f"bin_{name}"] = np.where(G[e] > 0, np.searchsorted(edges, G[e], side="right") + 1, 0)
        # clearly expressed: exon signal at least 3x its intronic surroundings, upper two fifths
        G[f"clear_{name}"] = (G[f"exon_{name}"] >= 3 * G[f"flank_{name}"]) & (G[f"bin_{name}"] >= 4)
        print("Detection by expression bin (0 = no signal above surroundings; 1-5 = fifths of expressed stand-alone genes):")
        print(f"  {'bin':<6}" + "".join(f"{t:>26}" for t in ("stand-alone pc", "nested pc", "stand-alone lncRNA", "nested lncRNA")))
        for dd in range(0, 6):
            cells = []
            for kind, ty in (("standalone", "protein_coding"), ("nested", "protein_coding"), ("standalone", "lncRNA"), ("nested", "lncRNA")):
                x = G[(G.kind == kind) & (G.type == ty) & (G[f"bin_{name}"] == dd)]
                cells.append(f"{x[dcol].mean() * 100:>16.0f}% (n={len(x):>5})" if len(x) else f"{'–':>26}")
            print(f"  {dd:<6}" + "".join(cells))
        print("Observed vs expected nested genes detected (expected from stand-alone genes of the same type and bin):")
        for lab_, extra in (("all", None), ("clearly expressed only", f"clear_{name}")):
            for ty in ("protein_coding", "lncRNA"):
                base = G[(G.kind == "standalone") & (G.type == ty)]
                if extra:
                    base = base[base[extra]]
                rate = base.groupby(f"bin_{name}")[dcol].mean()
                for o in ("same", "anti"):
                    x = G[(G.kind == "nested") & (G.type == ty) & (G.orient == o)]
                    if extra:
                        x = x[x[extra]]
                    if not len(x):
                        continue
                    exp_ = x[f"bin_{name}"].map(rate).fillna(0).sum()
                    obs = int(x[dcol].sum())
                    print(f"  {lab_:<24}{ty:<16}{o + ' strand':<13} n={len(x):>5}  observed {obs:>5}  expected {exp_:>7.1f}  "
                          f"ratio {obs / exp_ if exp_ else float('nan'):.2f}")
        sel = (G.kind == "nested") & G.well_established
        lab = "well-established nested genes"
        if G.broad.notna().any():
            sel &= G.broad.astype(bool)
            lab = "well-established, broadly expressed nested genes"
        x = G[sel & G[f"clear_{name}"]]
        print(f"  {lab}, clearly expressed in {name}: {len(x)}, detected {int(x[dcol].sum())}"
              f" ({x[dcol].mean() * 100 if len(x) else float('nan'):.0f}%)")
        dn = G[(G.kind == "nested") & G[dcol] & (G.orient == "same")]
        if len(dn):
            print(f"  detected nested genes on the host strand also found among the de novo intronic start sites (within 500 nt): "
                  f"{int(dn[f'itss_{name}'].sum())} of {len(dn)} ({dn[f'itss_{name}'].mean() * 100:.0f}%)")
        x = G[(G.kind == "nested") & (G.type == "protein_coding") & G[f"clear_{name}"]]
        print(f"  clearly expressed nested protein-coding genes in {name}: {len(x)}, detected {int(x[dcol].sum())}"
              f" ({x[dcol].mean() * 100 if len(x) else float('nan'):.0f}%); same strand {int((x.orient == 'same').sum())}, "
              f"opposite strand {int((x.orient == 'anti').sum())}")
        for r in x.sort_values(e, ascending=False).itertuples():
            print(f"    {r.name:<14} host {r.host:<12} {r.orient:<5} {'detected' if getattr(r, dcol) else 'NOT detected':<13} "
                  f"capped reads {getattr(r, 'cap_' + name):>6.0f}  exon/flank {getattr(r, 'exon_' + name):.3g}/"
                  f"{getattr(r, 'flank_' + name):.3g}  {'well-established' if r.well_established else ''}"
                  f"{' broad' if r.broad == True else ''}")

    sel = (G.kind == "nested") & G.well_established
    if G.broad.notna().any():
        sel &= G.broad.astype(bool)
    allok = np.all([G[f"det_{n}"] for n in lines], axis=0)
    ex = G[sel & allok].copy()
    ex["capsum"] = sum(ex[f"cap_{n}"] for n in lines)
    print(f"\nExamples: well-established{' broadly expressed' if G.broad.notna().any() else ''} nested genes detected in all "
          f"cell lines ({len(ex)}); top 25 by capped reads:")
    print(f"  {'gene':<16}{'type':<16}{'host':<14}{'strand':<7}" + "".join(f"{'cap ' + n:>12}{'bin':>8}" for n in lines))
    for r in ex.sort_values("capsum", ascending=False).head(25).itertuples():
        print(f"  {r.name:<16}{r.type:<16}{r.host:<14}{r.orient:<7}" + "".join(
            f"{getattr(r, 'cap_' + n):>12.0f}{getattr(r, 'bin_' + n):>8}" for n in lines))

    G.drop(columns=["exons"]).to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
