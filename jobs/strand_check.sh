#!/bin/bash
#SBATCH --job-name=istart_strand
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
{ source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate envs/istart; } || export PATH=$PWD/envs/istart/bin:$PATH
mkdir -p results/strand_check
python scripts/46_strand_check.py --genes results/known_genes_wide.tsv --gtf ref/gencode.v47.annotation.gtf.gz \
  --line "K562|data/rnaseq/K562_total.junctions.tsv|data/rnaseq/K562_total_plus.bw,data/rnaseq/K562_total_minus.bw" \
  --line "HCT116|data/hct116/RNA.junctions.tsv|data/hct116/RNA_plus.bw,data/hct116/RNA_minus.bw" \
  --out-prefix results/strand_check/strand
echo "==== done ===="
