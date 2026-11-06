#!/bin/bash
#SBATCH --job-name=istart_ltr12c_gata
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.log
# Do LTR12C copies that contain a confirmed intronic start site (active host) carry a GATA motif
# more often than the other LTR12C copies? Analogous to Brocks et al. 2017, Fig. 5e.
# Controls: NF-Y and SP1 motifs (also in LTR12C), full-length copies (>= 1200 bp), the HCT116 line.
set -eo pipefail
trap 'echo "ERROR: command at line $LINENO exited with code $?"' ERR
# project directory = the directory sbatch was called from (ISTART_DIR only for testing)
cd "${ISTART_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[ -d results/v2 ] && [ -d envs/homer ] || { echo "This is not the iSTART project directory ($PWD). Run sbatch from the iSTART directory."; exit 1; }
echo "==== start: $(date)  directory: $PWD"
export O=results/motifs/ltr12c; mkdir -p $O logs
GENOME=ref/GRCh38.primary_assembly.genome.fa
PY=$PWD/envs/istart/bin/python
export PATH=$PWD/envs/homer/bin:$PATH
for f in $GENOME ref/repeats/rmsk.txt.gz results/v2/itss.describe.tsv results/hct116/itss.describe.tsv $PY; do
  [ -e "$f" ] || { echo "Missing file: $f"; exit 1; }
done

# 1) Motifs from the HOMER library: GATA2 (test), NF-Y and SP1 (controls)
echo "==== 1) looking for HOMER motif files"
F=$(find "$PWD/envs/homer" -name 'gata2.motif' -print -quit 2>/dev/null || true)
if [ -z "$F" ]; then
  echo "gata2.motif not found in envs/homer. GATA/NF-Y/SP1 motif files that are there:"
  find "$PWD/envs/homer" -iname '*.motif' 2>/dev/null | grep -iE 'gata|nfy|sp1' | head -20 || true
  exit 1
fi
export MDIR=$(dirname "$F"); echo "     motif library: $MDIR"
for m in gata2 nfy sp1; do [ -s $MDIR/$m.motif ] || { echo "$m.motif missing in $MDIR"; exit 1; }; done
cat $MDIR/gata2.motif $MDIR/nfy.motif $MDIR/sp1.motif > $O/motifs.motif

# 2) LTR12C copies (RepeatMasker), flagging those that contain a confirmed start site (K562, HCT116)
echo "==== 2) LTR12C copies and the start sites inside them"
$PY - <<'PY'
import os, numpy as np, pandas as pd
O = os.environ["O"]
cols = ["bin", "sw", "dv", "dl", "ins", "chrom", "start", "end", "left", "strand",
        "name", "cls", "fam", "rs", "re", "rl", "rid"]
r = pd.read_csv("ref/repeats/rmsk.txt.gz", sep="\t", header=None, names=cols,
                usecols=["chrom", "start", "end", "strand", "name"])
r = r[(r.name == "LTR12C") & r.chrom.str.fullmatch(r"chr[0-9XY]+")].reset_index(drop=True)
r["id"], r["len"] = [f"LTR12C_{i}" for i in range(len(r))], r.end - r.start
for cell, path in (("K562", "results/v2/itss.describe.tsv"), ("HCT116", "results/hct116/itss.describe.tsv")):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    for c in ("host_active", "tss_supported"):
        d[c] = d[c].astype(str) == "True"
    d = d[(d.site_type == "candidate") & d.tss_supported & d.host_active]
    n_in, n_sense = np.zeros(len(r), int), np.zeros(len(r), int)
    for ch, g in d.groupby("chrom"):
        e = r[r.chrom == ch]
        for pos, st in zip(g.pos.astype(int), g.strand):
            hit = e.index[(e.start <= pos) & (e.end > pos)]
            n_in[hit] += 1
            n_sense[hit] += (e.loc[hit, "strand"].values == st).astype(int)
    r[f"{cell}_itss"], r[f"{cell}_sense"] = n_in, n_sense
    print(f"{cell}: {n_in.sum()} of {len(d)} start sites lie in LTR12C; copies with a start: {(n_in > 0).sum()} "
          f"of {len(r)}; starts in the same orientation as the LTR: {n_sense.sum()}/{n_in.sum()}")
r.to_csv(f"{O}/ltr12c_copies.tsv", sep="\t", index=False)
r.assign(sc=0)[["chrom", "start", "end", "id", "sc", "strand"]].to_csv(
    f"{O}/ltr12c.bed", sep="\t", header=False, index=False)
PY

# 3) Motif occurrences in whole LTR12C copies
echo "==== 3) HOMER: searching for motifs (log: $O/homer.log)"
findMotifsGenome.pl $O/ltr12c.bed $GENOME $O/homer -size given -find $O/motifs.motif \
    -p 8 -preparsedDir results/motifs/preparsed > $O/hits.txt 2> $O/homer.log

echo "     motif hits: $(($(wc -l < $O/hits.txt) - 1))"

# 4) Fraction of copies with the motif: with a start vs without (Fisher test)
echo "==== 4) comparison"
$PY - <<'PY'
import os, pandas as pd
from scipy.stats import fisher_exact
O, MDIR, MOTS = os.environ["O"], os.environ["MDIR"], ["GATA2", "NFY", "SP1"]
lab = {line.rstrip("\n").split("\t")[1]: m.upper()
       for m in ("gata2", "nfy", "sp1") for line in open(f"{MDIR}/{m}.motif") if line.startswith(">")}
h = pd.read_csv(f"{O}/hits.txt", sep="\t")
h["motif"] = h["Motif Name"].map(lab)
if h.motif.isna().all():
    raise SystemExit(f"Motif names not recognised: {h['Motif Name'].unique()[:5]}")
cnt = h.groupby(["PositionID", "motif"]).size().unstack(fill_value=0).reindex(columns=MOTS, fill_value=0)
r = pd.read_csv(f"{O}/ltr12c_copies.tsv", sep="\t").set_index("id").join(cnt).fillna({m: 0 for m in MOTS})
rows = []
for cell in ("K562", "HCT116"):
    for sub, keep in (("wszystkie kopie", r.len > 0), ("kopie >= 1200 pz", r.len >= 1200)):
        x = r[keep]
        has = x[f"{cell}_itss"] > 0
        for mot in MOTS:
            pos = x[mot] > 0
            t = [[int((has & pos).sum()), int((has & ~pos).sum())], [int((~has & pos).sum()), int((~has & ~pos).sum())]]
            odds, p = fisher_exact(t) if sum(t[0]) and sum(t[1]) else (float("nan"), float("nan"))
            rows.append({"linia": cell, "zbiór": sub, "motyw": mot,
                         "n_ze_startem": sum(t[0]), "%_ze_startem": 100 * t[0][0] / max(sum(t[0]), 1),
                         "n_bez": sum(t[1]), "%_bez": 100 * t[1][0] / max(sum(t[1]), 1),
                         "mediana_dł_ze_startem": x.len[has].median(), "OR": odds, "p": p})
s = pd.DataFrame(rows)
s.to_csv(f"{O}/summary.tsv", sep="\t", index=False)
print(s.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
PY
echo "==== done ===="