#!/bin/bash
#SBATCH --job-name=istart_v2
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate envs/istart
V=results/v2; mkdir -p $V
G=ref/gencode.v47.annotation.gtf.gz; PAS=ref/atlas.clusters.2.0.GRCh38.96.bed.gz
CAP="--cap-plus data/cap_plus.bw --cap-minus data/cap_minus.bw"
UNCAP="--uncap-plus data/uncap_plus.bw --uncap-minus data/uncap_minus.bw"
B=results/brudrb; H=data/encode/histone; R=data/rnaseq
DRB="--drb-plus $B/K562_DRB10.read_plus.bw --drb-minus $B/K562_DRB10.read_minus.bw --bru-plus $B/K562_Bru10.read_plus.bw --bru-minus $B/K562_Bru10.read_minus.bw"
# grep/sed filters below accept both the original Polish wording (logs and summaries written by the
# Polish version of the scripts) and the current English wording.
echo "== 10: candidates =="
python scripts/10_candidates.py --gtf $G --out $V 2>&1 | grep -E "HOST|CAND|CTRL"
echo "== 11: evidence for the catalogue =="
python scripts/11_evidence.py --sites $V/sites.tsv --host-tss $V/host_tss.tsv $CAP $UNCAP --pas $PAS --out-prefix $V/annot \
  | grep -E "promotory|kontrole|kandydaci|klasa|promoters|controls|candidates|Candidate genes|class [BC-]:"
echo "== 12: intronic start sites =="
python scripts/12_itss.py --gtf $G --uni data/encode/ENCSR261KBX/ENCFF636REV.bed.gz --bi data/encode/ENCSR261KBX/ENCFF819CPA.bed.gz \
  $CAP --pas $PAS --annot-sites $V/sites.tsv --out $V/sites_itss.tsv 2>&1 | grep -E "iTSS|OUT"
echo "== 11: evidence for iTSS =="
python scripts/11_evidence.py --sites $V/sites_itss.tsv --host-tss $V/host_tss.tsv $CAP $UNCAP --pas $PAS --out-prefix $V/itss \
  | grep -E "promotory|kontrole|kandydaci|klasa|G1|promoters|controls|candidates|Candidate genes|class [BC-]:"
echo "== 20: DRB =="
python scripts/20_drb_gap.py --evidence $V/itss.evidence.tsv $DRB --win 8000 --gap 5000 --out-prefix $V/itss_w8 \
  | sed -n -E '/Kalibracja empiryczna|Empirical calibration/,/szacunek|estimate/p;/skok|step/p'
echo "== 21: chromatin =="
python scripts/21_describe.py --drb $V/itss_w8.drb.tsv \
  --peaks H3K4me3=$H/H3K4me3.bed.gz H3K27ac=$H/H3K27ac.bed.gz H3K4me1=$H/H3K4me1.bed.gz --out-prefix $V/itss 2>/dev/null
echo "== phyloP =="
python - <<'PY'
import pandas as pd
d = pd.read_csv("results/v2/itss.describe.tsv", sep="\t", low_memory=False).reset_index()
pd.DataFrame({"c": d.chrom, "s": (d.pos.astype(int) - 100).clip(lower=0), "e": d.pos.astype(int) + 101,
              "n": "r" + d["index"].astype(str)}).to_csv("results/v2/sites_pm100.bed", sep="\t", header=False, index=False)
PY
T=${TMPDIR:-/tmp}/istart_$SLURM_JOB_ID; mkdir -p $T
cp ref/phylop/hg38.phyloP100way.bw $T/
bigWigAverageOverBed $T/hg38.phyloP100way.bw $V/sites_pm100.bed $V/sites_phylop.tab
rm -rf $T
sed 's#results/modA/#results/v2/#g' scripts/helpers/phylop_summary.py | python -
echo "== 22: stable RNA =="
python scripts/22_expression.py --gtf $G --cand-tx $V/candidates_tx.tsv --cand-genes $V/annot.genes.tsv \
  --sites $V/itss.describe.tsv --junctions $R/K562_total.junctions.tsv \
  --plus $R/K562_total_plus.bw --minus $R/K562_total_minus.bw --out-prefix $V/expr 2>&1 | grep -v "^\[GTF\]"
echo "== 23: candidates =="
python scripts/23_candidates.py --describe $V/itss.describe.tsv --phylop $V/sites_phylop.tab \
  --expr $V/expr.itss_expression.tsv --out $V/candidates_ranked.tsv --top 25
echo "== done =="
