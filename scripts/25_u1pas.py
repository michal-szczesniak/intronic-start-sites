#!/usr/bin/env python3
"""
iSTART: early-termination signals after the start site (U1–PAS axis).

A polymerase that, after the start site, meets a polyadenylation signal (PAS: AATAAA/ATTAAA) before it meets a
splice donor site recognised by U1 snRNP (5'SS) readily terminates transcription early; U1 bound
to a 5'SS protects against such cleavage. Classic pattern: in the gene direction (sense) 5'SS before PAS, in the
PROMPT direction (antisense) the reverse.

For each site, in the window of --win nt after the start site (on the transcribed strand):
  pas_n, ss_n          number of PAS and 5'SS motifs (consensus [AC]AG|GT[AG]AGT)
  protected            the first 5'SS lies before the first PAS (or there is no PAS)
Each site also has the opposite direction (anty, as for a PROMPT): the PAS sense/anti ratio compares two strands
with the same local base composition, which removes the effect of GC content. Groups: host promoters
(reference), internal promoters, iTSS (by chromatin type and host activity), random controls.
Optionally relates the "protection" of iTSS to stable RNA (rna_step) and to the response to MTREX knockout (delta).
"""
import argparse
import re

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
PAS = re.compile(r"(?=(AATAAA|ATTAAA))")
SS5 = re.compile(r"(?=([AC]AGGT[AG]AGT))")


class Fasta:
    def __init__(self, path):
        self.f = open(path, "rb")
        self.idx = {}
        for line in open(path + ".fai"):
            name, length, offset, lb, lw = line.split("\t")[:5]
            self.idx[name] = (int(length), int(offset), int(lb), int(lw))

    def get(self, chrom, s, e):
        if chrom not in self.idx:
            return ""
        L, off, lb, lw = self.idx[chrom]
        s, e = max(0, s), min(L, e)
        if e <= s:
            return ""
        a = off + (s // lb) * lw + s % lb
        b = off + (e // lb) * lw + e % lb
        self.f.seek(a)
        return self.f.read(b - a).decode().replace("\n", "").upper()


def downstream(fa, chrom, pos, strand, w):
    if strand == "+":
        return fa.get(chrom, pos, pos + w)
    return fa.get(chrom, pos - w + 1, pos + 1).translate(COMP)[::-1]


def scan(seq):
    pas = [m.start() for m in PAS.finditer(seq)]
    ss = [m.start() + 3 for m in SS5.finditer(seq)]
    first_pas = pas[0] if pas else None
    first_ss = ss[0] if ss else None
    protected = (first_pas is None) or (first_ss is not None and first_ss < first_pas)
    return len(pas), len(ss), protected


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True, help="itss.describe.tsv (step 21)")
    ap.add_argument("--fasta", required=True, help="genome FASTA with a .fai index")
    ap.add_argument("--expr", help="expr.itss_expression.tsv (step 22), column rna_step")
    ap.add_argument("--exo", help="exo.MTREX.tsv (step 24), column delta")
    ap.add_argument("--win", type=int, default=1000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    fa = Fasta(a.fasta)
    s = pd.read_csv(a.sites, sep="\t", low_memory=False)
    for col in ("host_active", "tss_supported"):
        s[col] = s[col].astype(str) == "True"
    rows = []
    for r in s.itertuples():
        p = int(r.pos)
        rec = dict(site_id=r.site_id)
        ok = True
        for lab, st in (("sens", r.strand), ("anty", "-" if r.strand == "+" else "+")):
            q = p if lab == "sens" else (p - 300 if r.strand == "+" else p + 300)  # antisense direction (as for a PROMPT)
            seq = downstream(fa, r.chrom, q, st, a.win)
            if len(seq) < a.win * 0.9:
                ok = False
                break
            n_pas, n_ss, prot = scan(seq)
            rec.update({f"pas_{lab}": n_pas, f"ss_{lab}": n_ss, f"prot_{lab}": prot,
                        f"gc_{lab}": (seq.count("G") + seq.count("C")) / len(seq)})
        if ok:
            rows.append(rec)
    d = s.merge(pd.DataFrame(rows), on="site_id")
    d["protected"] = d.prot_sens
    if a.expr:
        e = pd.read_csv(a.expr, sep="\t", low_memory=False)[["site_id", "rna_step"]]
        d = d.merge(e, on="site_id", how="left")
    if a.exo:
        x = pd.read_csv(a.exo, sep="\t", low_memory=False)[["site_id", "delta"]]
        d = d.merge(x, on="site_id", how="left")
    d.to_csv(a.out, sep="\t", index=False)

    cand = d.site_type == "candidate"
    groups = [("host promoters (anti = PROMPT)", d.site_type == "host_promoter"),
              ("capped internal promoters", (d.site_type == "internal_promoter") & d.host_active & d.tss_supported),
              ("capped iTSS, active host", cand & d.host_active & d.tss_supported),
              ("  of which promoter-like", cand & d.host_active & d.tss_supported & (d.typ == "promotor-podobny")),
              ("  of which enhancer-like", cand & d.host_active & d.tss_supported & (d.typ == "enhancer-podobny")),
              ("capped iTSS, inactive host", cand & ~d.host_active & d.tss_supported),
              ("random controls (intron background)", (d.site_type == "random_control") & d.host_active)]
    k = 1000 / a.win
    print(f"Signals in the first {a.win} nt after the start site: direction of transcription (sense) and opposite (anti):")
    print(f"  {'group':<40}{'n':>6}{'PAS/kb sense':>13}{'PAS/kb anti':>13}{'sense/anti':>11}{'GC sense':>9}{'protected*':>12}")
    for name, msk in groups:
        g = d[msk]
        if not len(g):
            continue
        ps, pa = g.pas_sens.mean() * k, g.pas_anty.mean() * k
        print(f"  {name:<40}{len(g):>6}{ps:>13.2f}{pa:>13.2f}{(ps / pa if pa > 0 else float('nan')):>11.2f}{g.gc_sens.mean() * 100:>8.0f}%"
              f"{g.prot_sens.mean() * 100:>11.0f}%")
    print("  * sense direction: first 5'SS before the first PAS, or no PAS in the window")
    print("  sense/anti < 1: PAS depletion on the transcribed strand vs the opposite strand of the same composition")

    it = d[cand & d.host_active & d.tss_supported]
    if "rna_step" in it:
        ctrl = d[(d.site_type == "random_control") & d.host_active]
        thr = np.nanpercentile(ctrl.rna_step, 95)
        print(f"\nStable RNA after iTSS (active host) by protection (control threshold {thr:.2f}):")
        tab = []
        for prot, g in it.groupby("protected"):
            k_ = int((g.rna_step > thr).sum())
            tab.append([k_, len(g) - k_])
            print(f"  {'protected' if prot else 'unprotected':<14} n={len(g):>4}  RNA rise > threshold {k_ / len(g) * 100:5.1f}%")
        if len(tab) == 2:
            print(f"  Fisher's exact test (protected vs unprotected): p = {fisher_exact(tab)[1]:.3f}")
    if "delta" in it:
        ctrl = d[(d.site_type == "random_control") & d.host_active]
        thr = np.nanpercentile(ctrl.delta.dropna(), 95)
        print(f"\nResponse to MTREX knockout after iTSS (active host) by protection (control threshold {thr:.2f}):")
        for prot, g in it.groupby("protected"):
            gg = g.dropna(subset=["delta"])
            print(f"  {'protected' if prot else 'unprotected':<14} n={len(gg):>4}  Δ > threshold {(gg.delta > thr).mean() * 100:5.1f}%"
                  f"  median Δ {gg.delta.median():5.2f}")


if __name__ == "__main__":
    main()
