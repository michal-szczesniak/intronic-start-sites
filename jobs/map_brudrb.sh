#!/bin/bash
#SBATCH --job-name=istart_map
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate envs/istart
D=data/brudrb
B=results/brudrb
mkdir -p $B/qc
cut -f1,2 ref/GRCh38.primary_assembly.genome.fa.fai > ref/grch38.chrom.sizes
ADAPT=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

for R in SRR1182352 SRR1182353 SRR1182354 SRR1182355; do
  echo "== $R =="
  fastp -i $D/$R.fastq.gz -o $B/$R.trim.fastq.gz -w 8 -a $ADAPT -l 20 \
        -j $B/qc/$R.fastp.json -h $B/qc/$R.fastp.html 2> $B/qc/$R.fastp.log
  hisat2 -p 8 -x ref/hisat2/grch38 -U $B/$R.trim.fastq.gz --no-unal 2> $B/qc/$R.hisat2.log \
    | samtools view -@ 2 -b -q 10 - | samtools sort -@ 4 -m 1G -o $B/$R.bam -
  samtools index $B/$R.bam
  rm $B/$R.trim.fastq.gz
done

echo "== merging runs into samples =="
samtools merge -f -@ 8 $B/K562_Bru10.bam $B/SRR1182352.bam $B/SRR1182353.bam
samtools merge -f -@ 8 $B/K562_DRB10.bam $B/SRR1182354.bam $B/SRR1182355.bam

echo "== coverage bigWigs by read strand =="
for S in K562_Bru10 K562_DRB10; do
  samtools index $B/$S.bam
  for ST in + -; do
    N=$([ "$ST" = "+" ] && echo plus || echo minus)
    bedtools genomecov -ibam $B/$S.bam -bg -split -strand $ST | LC_ALL=C sort -k1,1 -k2,2n > $B/$S.$N.bg
    bedGraphToBigWig $B/$S.$N.bg ref/grch38.chrom.sizes $B/$S.read_$N.bw
    rm $B/$S.$N.bg
  done
done

echo "== summary =="
for R in SRR1182352 SRR1182353 SRR1182354 SRR1182355; do
  echo "$R: $(grep 'overall alignment rate' $B/qc/$R.hisat2.log); after MAPQ>=10 filter: $(samtools view -c $B/$R.bam)"
done
ls -lh $B/*.bw
echo "== done =="
