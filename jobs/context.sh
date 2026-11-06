#!/bin/bash
#SBATCH --job-name=istart_ctx
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
{ source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate envs/istart; } || export PATH=$PWD/envs/istart/bin:$PATH
G=ref/gencode.v47.annotation.gtf.gz; FA=ref/GRCh38.primary_assembly.genome.fa
python scripts/29_context.py --name K562 --gtf $G --describe results/v2/itss.describe.tsv \
  --uni data/encode/ENCSR261KBX/ENCFF636REV.bed.gz --bi data/encode/ENCSR261KBX/ENCFF819CPA.bed.gz \
  --cap data/cap_plus.bw,data/cap_minus.bw --uncap data/uncap_plus.bw,data/uncap_minus.bw \
  --rna data/rnaseq/K562_total_plus.bw,data/rnaseq/K562_total_minus.bw \
  --ko data/exosome/MTREX_KO_plus.bw,data/exosome/MTREX_KO_minus.bw \
  --ko-ctrl data/exosome/MTREX_CTRL_plus.bw,data/exosome/MTREX_CTRL_minus.bw \
  --fasta $FA --out results/context_K562.tsv 2>&1 | grep -v "^\[GTF\]"
python scripts/29_context.py --name HCT116 --gtf $G --describe results/hct116/itss.describe.tsv \
  --uni data/hct116/CTCFline_uni.bed.gz --bi data/hct116/CTCFline_bi.bed.gz \
  --cap data/hct116/CTCFline_plus.bw,data/hct116/CTCFline_minus.bw \
  --rna data/hct116/RNA_plus.bw,data/hct116/RNA_minus.bw \
  --drb data/hct116_drb/TT_10min.read_plus.bw,data/hct116_drb/TT_10min.read_minus.bw \
  --drb-ref data/hct116_drb/TT_40min.read_plus.bw,data/hct116_drb/TT_40min.read_minus.bw \
  --fasta $FA --out results/context_HCT116.tsv 2>&1 | grep -v "^\[GTF\]"
echo "==== done ===="
