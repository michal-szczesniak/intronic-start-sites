#!/bin/bash
#SBATCH --job-name=istart_exo
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
{ source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate envs/istart; } || export PATH=$PWD/envs/istart/bin:$PATH
X=data/exosome; mkdir -p $X results/v2/exosome
for spec in "MTREX ENCSR300UBK ENCSR271ZAB" "EXOSC10 ENCSR299PRA ENCSR854OFV" "EXOSC5 ENCSR406SEH ENCSR341TTW"; do
  set -- $spec
  echo "== $1 =="
  python scripts/helpers/encode_rnaseq_bw.py $1 $2 $3 $X
  python scripts/24_exosome.py --sites results/v2/itss.describe.tsv \
    --ko-plus $X/$1_KO_plus.bw --ko-minus $X/$1_KO_minus.bw --ctrl-plus $X/$1_CTRL_plus.bw --ctrl-minus $X/$1_CTRL_minus.bw \
    --name $1 --out-prefix results/v2/exosome/exo
done
echo "== done =="
