#!/bin/bash
#SBATCH --job-name=istart_cls
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
{ source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate envs/istart; } || export PATH=$PWD/envs/istart/bin:$PATH
export D=data/hct116_drb; O=results/hct116
echo "==== 1. bigWigs of individual runs ===="
mk() { r=$1; for ST in + -; do NM=$([ "$ST" = "+" ] && echo plus || echo minus)
  [ -s $D/$r.read_$NM.bw ] && continue
  bedtools genomecov -ibam $D/$r.bam -bg -split -strand $ST | LC_ALL=C sort -S 2G -k1,1 -k2,2n > $D/$r.$NM.bg
  bedGraphToBigWig $D/$r.$NM.bg ref/grch38.chrom.sizes $D/$r.read_$NM.bw; rm $D/$r.$NM.bg; done; echo "$r ready"; }
export -f mk
printf "%s\n" SRR11563139 SRR11563142 SRR11563147 SRR11563150 SRR11563140 SRR11563148 | xargs -P 6 -I{} bash -c 'mk {}'
B() { echo $D/$1.read_plus.bw,$D/$1.read_minus.bw; }
VAL="--u1pas $O/u1pas.tsv --expr $O/expr.itss_expression.tsv --describe $O/itss.describe.tsv --phylop $O/sites_phylop.tab"
echo "==== 2. classification: 10 min vs 40 min ===="
python scripts/27_classify.py --evidence $O/itss.evidence.tsv --test $(B SRR11563139) $(B SRR11563147) \
  --ref $(B SRR11563142) $(B SRR11563150) --front 16000 $VAL --out $O/classes_10.tsv
echo "==== 3. classification: 20 min vs 40 min ===="
python scripts/27_classify.py --evidence $O/itss.evidence.tsv --test $(B SRR11563140) $(B SRR11563148) \
  --ref $(B SRR11563142) $(B SRR11563150) --front 44000 $VAL --out $O/classes_20.tsv
echo "==== 4. agreement of 10- and 20-min classes (starts classified in both) ===="
python - <<'PY'
import pandas as pd
a = pd.read_csv("results/hct116/classes_10.tsv", sep="\t", low_memory=False)[["site_id", "site_type", "klasa"]]
b = pd.read_csv("results/hct116/classes_20.tsv", sep="\t", low_memory=False)[["site_id", "klasa"]]
m = a.merge(b, on="site_id", suffixes=("_10", "_20"))
m = m[(m.site_type == "candidate") & ~m.klasa_10.str.startswith("nieroz") & ~m.klasa_20.str.startswith("nieroz")]
print(pd.crosstab(m.klasa_10, m.klasa_20, margins=True))
PY
echo "==== done ===="
