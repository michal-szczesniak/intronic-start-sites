#!/bin/bash
#SBATCH --job-name=istart_loci
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
{ source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate envs/istart; } || export PATH=$PWD/envs/istart/bin:$PATH
T=${TMPDIR:-/tmp}/istart_$SLURM_JOB_ID; mkdir -p $T
cp ref/phylop/hg38.phyloP100way.bw $T/
python scripts/32_loci.py --phylop $T/hg38.phyloP100way.bw --pgt $PWD/envs/pgt/bin/pyGenomeTracks
rm -rf $T
echo "==== done ===="
