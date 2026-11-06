#!/usr/bin/env python3
"""Sums several bigWigs (e.g. replicates) into one; takes absolute values,
so a minus strand stored as negative numbers becomes positive.
Usage: python scripts/helpers/merge_bw.py output.bw input1.bw input2.bw ..."""
import sys

import numpy as np
import pyBigWig

out, ins = sys.argv[1], sys.argv[2:]
bws = [pyBigWig.open(f) for f in ins]
chroms = {}
for bw in bws:
    for c, L in bw.chroms().items():
        chroms[c] = max(chroms.get(c, 0), L)
order = sorted(chroms)
o = pyBigWig.open(out, "w")
o.addHeader([(c, chroms[c]) for c in order])
total = 0.0
for c in order:
    pos_l, val_l = [], []
    for bw in bws:
        if c not in bw.chroms():
            continue
        iv = bw.intervals(c)
        if not iv:
            continue
        a = np.array(iv, dtype=np.float64)
        s, e, v = a[:, 0].astype(np.int64), a[:, 1].astype(np.int64), np.abs(a[:, 2])
        ln = e - s
        offs = np.arange(ln.sum()) - np.repeat(np.cumsum(ln) - ln, ln)
        pos_l.append(np.repeat(s, ln) + offs)
        val_l.append(np.repeat(v, ln))
    if not pos_l:
        continue
    pos = np.concatenate(pos_l)
    val = np.concatenate(val_l)
    u, inv = np.unique(pos, return_inverse=True)
    summed = np.bincount(inv, weights=val)
    keep = summed > 0
    u, summed = u[keep], summed[keep]
    o.addEntries([c] * len(u), u.tolist(), ends=(u + 1).tolist(), values=summed.tolist())
    total += summed.sum()
o.close()
print(f"{out}: total signal = {total:.0f}")
