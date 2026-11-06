import importlib.util, sys
from collections import defaultdict, Counter
def load(n, f):
    s = importlib.util.spec_from_file_location(n, f); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
c10 = load("c10", "scripts/10_candidates.py"); c12 = load("c12", "scripts/12_itss.py")
GTF, UNI, BI = sys.argv[1], sys.argv[2], sys.argv[3]
tx, exons = c10.parse_gtf(GTF)
tss, ex, pc, oth = defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list)
for t in tx.values():
    k = (t["chrom"], t["strand"])
    tss[k].append(c10.tss_of(t))
    ex[k] += exons[t["transcript_id"]]
    (pc if t["gene_type"] == "protein_coding" else oth)[k].append((t["start"], t["end"]))
tss = {k: sorted(v) for k, v in tss.items()}
EX = {k: c10.IntervalSet(v) for k, v in ex.items()}
PC = {k: c10.IntervalSet(v) for k, v in pc.items()}
OT = {k: c10.IntervalSet(v) for k, v in oth.items()}
E = c10.IntervalSet([])
starts = [(c, s, p, "pojedynczy") for c, s, p, _, _ in c12.read_uni(UNI)]
for c, p5, p6, _ in c12.read_bi(BI):
    starts += [(c, "+", p5, "z pary"), (c, "-", p6, "z pary")]
def cat(c, s, p):
    k, o = (c, s), (c, "-" if s == "+" else "+")
    tl = tss.get(k, [])
    j = c10.bisect.bisect_left(tl, p - 500)
    if j < len(tl) and tl[j] <= p + 500:
        return "1 at an annotated start (±500 nt)"
    if EX.get(k, E).dist(p) == 0:
        return "2 in an exon, same strand"
    if PC.get(k, E).overlaps(p, p + 1):
        return "3 protein-coding intron, same strand, ≥1 kb from exons" if EX.get(k, E).dist(p) >= 1000 \
            else "4 protein-coding intron, same strand, <1 kb from exons"
    if OT.get(k, E).overlaps(p, p + 1):
        return "5 intron of another gene (e.g. lncRNA), same strand"
    if PC.get(o, E).overlaps(p, p + 1) or OT.get(o, E).overlaps(p, p + 1):
        return "6 in a gene on the opposite strand (antisense)"
    return "7 between genes"
cnt = Counter((cat(c, s, p), kind) for c, s, p, kind in starts)
tot = len(starts)
print(f"All PRO-cap starts: {tot}")
print(f"{'category':<58}{'total':>8}{'%':>7}{'single':>12}{'paired':>8}")
for k in sorted({k for k, _ in cnt}):
    a, b = cnt[(k, 'pojedynczy')], cnt[(k, 'z pary')]
    print(f"{k:<58}{a + b:>8}{(a + b) / tot * 100:>6.1f}%{a:>12}{b:>8}")
if len(sys.argv) > 4:
    import pandas as pd
    s = pd.read_csv(sys.argv[4], sep="\t", low_memory=False)
    s = s[s.site_type == "candidate"]
    near = 0
    for r in s.itertuples():
        tl = tss.get((r.chrom, r.strand), [])
        j = c10.bisect.bisect_left(tl, r.pos - 500)
        near += j < len(tl) and tl[j] <= r.pos + 500
    print(f"\nIntronic starts (iTSS, from step 12): {len(s)}; at any annotated start "
          f"(±500 nt, any gene or transcript on the same strand): {near} ({near / len(s) * 100:.1f}%)")
