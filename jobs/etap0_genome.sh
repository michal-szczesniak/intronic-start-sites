#!/bin/bash
#SBATCH --job-name=istart_e0
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate envs/istart
O=results/etap0
mkdir -p $O
CAP="--cap-plus data/cap_plus.bw --cap-minus data/cap_minus.bw"
UNCAP="--uncap-plus data/uncap_plus.bw --uncap-minus data/uncap_minus.bw"
PAS=ref/atlas.clusters.2.0.GRCh38.96.bed.gz

echo "== step 10: candidates from the annotation =="
python scripts/10_candidates.py --gtf ref/gencode.v47.annotation.gtf.gz --out $O

echo "== step 11: evidence for the annotation =="
python scripts/11_evidence.py --sites $O/sites.tsv --host-tss $O/host_tss.tsv $CAP $UNCAP --pas $PAS --out-prefix $O/annot

echo "== step 12: intronic start sites from PRO-cap =="
python scripts/12_itss.py --gtf ref/gencode.v47.annotation.gtf.gz \
  --uni data/encode/ENCSR261KBX/ENCFF636REV.bed.gz --bi data/encode/ENCSR261KBX/ENCFF819CPA.bed.gz \
  $CAP --pas $PAS --annot-sites $O/sites.tsv --out $O/sites_itss.tsv

echo "== step 11: evidence for iTSS =="
python scripts/11_evidence.py --sites $O/sites_itss.tsv --host-tss $O/host_tss.tsv $CAP $UNCAP --pas $PAS --out-prefix $O/itss
echo "== done =="
