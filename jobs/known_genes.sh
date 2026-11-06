#!/bin/bash
#SBATCH --job-name=istart_known
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
{ source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate envs/istart; } || export PATH=$PWD/envs/istart/bin:$PATH
mkdir -p ref/gtex
GX=ref/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
[ -s $GX ] || wget -nv -O $GX https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz || rm -f $GX
GT=""; [ -s $GX ] && GT="--gtex $GX" || echo "GTEx table not available: running without breadth of expression"
python scripts/30_known_genes.py --gtf ref/gencode.v47.annotation.gtf.gz $GT \
  --line "K562|data/cap_plus.bw,data/cap_minus.bw|data/uncap_plus.bw,data/uncap_minus.bw|data/rnaseq/K562_total_plus.bw,data/rnaseq/K562_total_minus.bw|results/v2/itss.describe.tsv" \
  --line "HCT116|data/hct116/CTCFline_plus.bw,data/hct116/CTCFline_minus.bw||data/hct116/RNA_plus.bw,data/hct116/RNA_minus.bw|results/hct116/itss.describe.tsv" \
  --out results/known_genes.tsv
echo "==== done ===="
