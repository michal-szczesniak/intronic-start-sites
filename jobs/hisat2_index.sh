#!/bin/bash
#SBATCH --job-name=istart_idx
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate envs/istart
cd ref
pigz -dc GRCh38.primary_assembly.genome.fa.gz > GRCh38.primary_assembly.genome.fa
samtools faidx GRCh38.primary_assembly.genome.fa
mkdir -p hisat2
hisat2-build -p 8 GRCh38.primary_assembly.genome.fa hisat2/grch38
ls -lh hisat2
echo "== done =="
