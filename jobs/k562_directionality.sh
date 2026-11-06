#!/bin/bash
#SBATCH --job-name=k562_directionality
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/k562_directionality_%x_%j.log
# Step 26 for K562; defaults written out explicitly.
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate envs/istart
for f in ref/gencode.v47.annotation.gtf.gz data/encode/ENCSR261KBX/ENCFF819CPA.bed.gz \
         data/cap_plus.bw data/cap_minus.bw results/v2/host_tss.tsv \
         results/brudrb/K562_Bru10.read_plus.bw results/brudrb/K562_Bru10.read_minus.bw \
         data/encode/histone/H3K36me3.bed.gz; do
  [ -e "$f" ] || { echo "Missing input: $f"; exit 1; }
done
python scripts/26_directionality.py --gtf ref/gencode.v47.annotation.gtf.gz \
  --bi data/encode/ENCSR261KBX/ENCFF819CPA.bed.gz \
  --cap-plus data/cap_plus.bw --cap-minus data/cap_minus.bw \
  --host-tss results/v2/host_tss.tsv \
  --bru-plus results/brudrb/K562_Bru10.read_plus.bw --bru-minus results/brudrb/K562_Bru10.read_minus.bw \
  --bru-flip yes --min-exon-dist 1000 \
  --h3k36me3 data/encode/histone/H3K36me3.bed.gz \
  --out results/v2/directionality.rerun.tsv
python - <<'PY'
import numpy as np
import pandas as pd
old = pd.read_csv("results/v2/directionality.tsv", sep="\t")
new = pd.read_csv("results/v2/directionality.rerun.tsv", sep="\t")
print(f"rows: old {len(old)}, rerun {len(new)}")
m = old.merge(new, on=["chrom", "mid"], suffixes=("", "_new"))
print(f"pairs in both: {len(m)}")
print(f"same class: {100 * (m.klasa == m.klasa_new).mean():.2f}%")
for c in ("sens", "anty", "DI", "host_signal", "host_tx", "k36"):
    x = pd.to_numeric(m[c], errors="coerce")
    y = pd.to_numeric(m[c + "_new"], errors="coerce")
    print(f"{c}: identical in {100 * np.isclose(x, y, equal_nan=True).mean():.2f}% of pairs")
PY
echo "==== done ===="
