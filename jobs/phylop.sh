#!/bin/bash
#SBATCH --job-name=istart_phylop
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate envs/istart
echo "== ±100 nt windows around the sites =="
python - <<'PY'
import pandas as pd
d = pd.read_csv("results/modA/itss.describe.tsv", sep="\t", low_memory=False).reset_index()
pd.DataFrame({"c": d.chrom, "s": (d.pos.astype(int) - 100).clip(lower=0), "e": d.pos.astype(int) + 101,
              "n": "r" + d["index"].astype(str)}).to_csv("results/modA/sites_pm100.bed", sep="\t", header=False, index=False)
print(len(d), "windows")
PY
T=${TMPDIR:-/tmp}/istart_$SLURM_JOB_ID
mkdir -p $T
echo "== copying phyloP to local disk =="
cp ref/phylop/hg38.phyloP100way.bw $T/
echo "== bigWigAverageOverBed =="
time bigWigAverageOverBed $T/hg38.phyloP100way.bw results/modA/sites_pm100.bed results/modA/sites_phylop.tab
rm -rf $T
python scripts/helpers/phylop_summary.py
echo "== done =="
