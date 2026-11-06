#!/bin/bash
#SBATCH --job-name=istart_supp
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
{ source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate envs/istart; } || export PATH=$PWD/envs/istart/bin:$PATH
ANNOT=""; for L in "K562 results/v2" "HCT116 results/hct116"; do set -- $L; [ -s $2/annot.evidence.tsv ] && ANNOT="$ANNOT $1=$2/annot.evidence.tsv"; done
RP=$(ls data/hct116/*RAD21*plus*.bw 2>/dev/null | head -1); RM=$(ls data/hct116/*RAD21*minus*.bw 2>/dev/null | head -1)
WAVE="K562:ref=results/brudrb/K562_Bru10.read_plus.bw,results/brudrb/K562_Bru10.read_minus.bw K562:DRB10=results/brudrb/K562_DRB10.read_plus.bw,results/brudrb/K562_DRB10.read_minus.bw"
WAVE="$WAVE HCT116:ref=data/hct116_drb/TT_40min.read_plus.bw,data/hct116_drb/TT_40min.read_minus.bw"
for t in 10 20 30; do f=data/hct116_drb/TT_${t}min.read_plus.bw; [ -s $f ] && WAVE="$WAVE HCT116:${t}min=$f,${f/plus/minus}"; done
echo "annot: $ANNOT"; echo "RAD21: $RP $RM"; echo "wave: $WAVE"
python scripts/34_supp_numbers.py --evidence K562=results/v2/itss.evidence.tsv HCT116=results/hct116/itss.evidence.tsv \
  ${ANNOT:+--annot $ANNOT} ${RP:+--rad21 $RP,$RM} --wave $WAVE \
  --context K562=results/context_K562.tsv HCT116=results/context_HCT116.tsv --out-prefix results/supp/supp
echo; echo "==== S4: K562 DRB test (summary of 20_drb_gap) ===="
# the sed filter below accepts both the original Polish wording (summaries written by the Polish
# version of 20_drb_gap.py) and the current English wording.
sed -n -E '/^(Kalibracja|Calibration|Empirical calibration)/,/^$/p;/^Profil/,/skok|step/p' results/v2/itss_w8.summary.txt
echo; echo "==== S5: exosome knockout checks ===="
ls results/v2/exosome/
for f in results/v2/exosome/*summary*; do [ -s "$f" ] && { echo "== $f"; cat "$f"; }; done
echo "==== done ===="
