#!/bin/bash
#SBATCH --job-name=istart_robust
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate envs/istart
O=results/robust500; mkdir -p $O
G=ref/gencode.v47.annotation.gtf.gz; PAS=ref/atlas.clusters.2.0.GRCh38.96.bed.gz
CAP="--cap-plus data/cap_plus.bw --cap-minus data/cap_minus.bw"
UNCAP="--uncap-plus data/uncap_plus.bw --uncap-minus data/uncap_minus.bw"
B=results/brudrb; H=data/encode/histone
DRB="--drb-plus $B/K562_DRB10.read_plus.bw --drb-minus $B/K562_DRB10.read_minus.bw --bru-plus $B/K562_Bru10.read_plus.bw --bru-minus $B/K562_Bru10.read_minus.bw"
# grep/sed filters below accept both the original Polish wording (logs and summaries written by the
# Polish version of the scripts) and the current English wording.
echo "== 500 nt: candidates and iTSS =="
python scripts/10_candidates.py --gtf $G --out $O --min-exon-dist 500 2>&1 | grep -E "CAND|CTRL"
python scripts/12_itss.py --gtf $G --uni data/encode/ENCSR261KBX/ENCFF636REV.bed.gz --bi data/encode/ENCSR261KBX/ENCFF819CPA.bed.gz \
  $CAP --pas $PAS --annot-sites $O/sites.tsv --out $O/sites_itss.tsv --min-exon-dist 500 2>&1 | grep -E "iTSS|OUT"
echo "== 500 nt: start-site evidence =="
python scripts/11_evidence.py --sites $O/sites_itss.tsv --host-tss $O/host_tss.tsv $CAP $UNCAP --pas $PAS \
  --out-prefix $O/itss --min-exon-dist 500 | grep -E "promotory|kontrole|kandydaci|klasa|G1|promoters|controls|candidates|Candidate genes|class [BC-]:"
echo "== 500 nt: DRB =="
python scripts/20_drb_gap.py --evidence $O/itss.evidence.tsv $DRB --win 8000 --gap 5000 --out-prefix $O/itss_w8 \
  | sed -n -E '/Kalibracja empiryczna|Empirical calibration/,/szacunek|estimate/p;/skok|step/p'
echo "== 500 nt: chromatin =="
python scripts/21_describe.py --drb $O/itss_w8.drb.tsv \
  --peaks H3K4me3=$H/H3K4me3.bed.gz H3K27ac=$H/H3K27ac.bed.gz H3K4me1=$H/H3K4me1.bed.gz --out-prefix $O/itss 2>/dev/null
for W in "6000 4000" "10000 6000"; do
  set -- $W
  echo "== main set, windows $1 / $2 =="
  python scripts/20_drb_gap.py --evidence results/etap0/itss.evidence.tsv $DRB --win $1 --gap $2 --out-prefix results/modA/itss_w$1 \
    | sed -n -E '/Kalibracja empiryczna|Empirical calibration/,/szacunek|estimate/p;/skok|step/p'
done
echo "== done =="
