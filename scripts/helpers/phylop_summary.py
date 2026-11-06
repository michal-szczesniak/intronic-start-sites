import pandas as pd
d = pd.read_csv("results/modA/itss.describe.tsv", sep="\t", low_memory=False).reset_index()
t = pd.read_csv("results/modA/sites_phylop.tab", sep="\t", header=None, names=["n", "size", "covered", "sum", "mean0", "mean"])
t["index"] = t.n.str[1:].astype(int)
d = d.merge(t[["index", "covered", "mean"]], on="index")
for c in ("host_active", "tss_supported"):
    d[c] = d[c].astype(str) == "True"
cand = d.site_type == "candidate"
G = [("host starts", d.site_type == "host_promoter"),
     ("capped internal promoters", (d.site_type == "internal_promoter") & d.host_active & d.tss_supported),
     ("capped iTSS, active host", cand & d.host_active & d.tss_supported),
     ("  of which promoter-like", cand & d.host_active & d.tss_supported & (d.typ == "promotor-podobny")),
     ("  of which enhancer-like", cand & d.host_active & d.tss_supported & (d.typ == "enhancer-podobny")),
     ("capped iTSS, inactive host", cand & ~d.host_active & d.tss_supported),
     ("uncapped iTSS", cand & ~d.tss_supported),
     ("random controls", (d.site_type == "random_control") & d.host_active)]
print(f"{'group':<38}{'n':>6}{'mean phyloP':>16}{'median':>10}{'% > 1':>8}")
for name, m in G:
    g = d[m & (d.covered > 0)]["mean"]
    print(f"{name:<38}{len(g):>6}{g.mean():>16.2f}{g.median():>10.2f}{(g > 1).mean() * 100:>7.0f}%")
