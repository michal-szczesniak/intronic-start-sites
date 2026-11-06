#!/bin/bash
#SBATCH --job-name=istart_rna
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate envs/istart
R=data/rnaseq; U=https://www.encodeproject.org/files
mkdir -p $R
echo "== download =="
for A in ENCFF932MJL ENCFF170FKF; do wget -nv -c -O $R/$A.bam $U/$A/@@download/$A.bam; done
for A in ENCFF312ZLI ENCFF272TZC ENCFF530FJG ENCFF123ORL; do wget -nv -c -O $R/$A.bigWig $U/$A/@@download/$A.bigWig; done
echo "== splice junctions (uniquely mapped reads) =="
for A in ENCFF932MJL ENCFF170FKF; do
  samtools view -@ 3 -F 0x904 -q 255 $R/$A.bam \
  | awk 'BEGIN{OFS="\t"} $6 ~ /N/ { pos=$4; cig=$6;
      while (match(cig, /^[0-9]+[MIDNSHP=X]/)) { len=substr(cig,1,RLENGTH-1)+0; op=substr(cig,RLENGTH,1);
        if (op=="N") print $3, pos-1, pos-1+len;
        if (op=="M"||op=="D"||op=="N"||op=="="||op=="X") pos+=len;
        cig=substr(cig,RLENGTH+1) } }' \
  | LC_ALL=C sort -S 8G -k1,1 -k2,2n -k3,3n | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1}' > $R/$A.junctions.tsv
  wc -l $R/$A.junctions.tsv
done
cat $R/ENCFF932MJL.junctions.tsv $R/ENCFF170FKF.junctions.tsv \
  | awk 'BEGIN{OFS="\t"}{k=$1"\t"$2"\t"$3; c[k]+=$4} END{for(k in c) print k, c[k]}' \
  | LC_ALL=C sort -k1,1 -k2,2n > $R/K562_total.junctions.tsv
wc -l $R/K562_total.junctions.tsv
echo "== merging replicate bigWigs =="
python scripts/helpers/merge_bw.py $R/K562_total_plus.bw  $R/ENCFF312ZLI.bigWig $R/ENCFF272TZC.bigWig
python scripts/helpers/merge_bw.py $R/K562_total_minus.bw $R/ENCFF530FJG.bigWig $R/ENCFF123ORL.bigWig
echo "== done =="
