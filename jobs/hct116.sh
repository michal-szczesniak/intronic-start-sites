#!/bin/bash
#SBATCH --job-name=istart_hct
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
{ source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate envs/istart; } || export PATH=$PWD/envs/istart/bin:$PATH
D=data/hct116; O=results/hct116; V=results/v2; mkdir -p $D $O
G=ref/gencode.v47.annotation.gtf.gz; PAS=ref/atlas.clusters.2.0.GRCh38.96.bed.gz; GET=scripts/helpers/encode_get.py
echo "==== 1. download ===="
python $GET signal CTCFline ENCSR165QKR $D --to-counts
python $GET peaks  CTCFline ENCSR165QKR $D
python $GET signal RAD21line ENCSR635SXX $D --to-counts
python $GET peaks  RAD21line ENCSR635SXX $D
for m in "H3K4me3 ENCSR333OPW" "H3K27ac ENCSR661KMA" "H3K4me1 ENCSR161MXP" "H3K36me3 ENCSR091QXP"; do
  set -- $m; python $GET histone $1 $2 $D
done
python $GET signal RNA ENCSR698RPL $D
python $GET bam RNA ENCSR698RPL $D
echo "==== 2. splice junctions ===="
for b in $D/RNA_rep*.bam; do
  samtools view -@ 3 -F 0x904 -q 255 $b \
  | awk 'BEGIN{OFS="\t"} $6 ~ /N/ { pos=$4; cig=$6;
      while (match(cig, /^[0-9]+[MIDNSHP=X]/)) { len=substr(cig,1,RLENGTH-1)+0; op=substr(cig,RLENGTH,1);
        if (op=="N") print $3, pos-1, pos-1+len;
        if (op=="M"||op=="D"||op=="N"||op=="="||op=="X") pos+=len;
        cig=substr(cig,RLENGTH+1) } }' \
  | LC_ALL=C sort -S 8G -k1,1 -k2,2n -k3,3n | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1}' > ${b%.bam}.junctions.tsv
done
cat $D/RNA_rep*.junctions.tsv | awk 'BEGIN{OFS="\t"}{k=$1"\t"$2"\t"$3; c[k]+=$4} END{for(k in c) print k, c[k]}' \
  | LC_ALL=C sort -k1,1 -k2,2n > $D/RNA.junctions.tsv
wc -l $D/RNA.junctions.tsv
CAP="--cap-plus $D/CTCFline_plus.bw --cap-minus $D/CTCFline_minus.bw"
# grep/sed filters below accept both the original Polish wording (logs and summaries written by the
# Polish version of the scripts) and the current English wording.
echo "==== 3. Result 1: catalogue ===="
python scripts/11_evidence.py --sites $V/sites.tsv --host-tss $V/host_tss.tsv $CAP --pas $PAS --out-prefix $O/annot \
  | grep -E "Orientacja|promotory|kontrole|kandydaci|klasa|Orientation|promoters|controls|candidates|Candidate genes|class [BC-]:"
echo "==== 4. Result 2: intronic start sites (CTCF line) ===="
python scripts/12_itss.py --gtf $G --uni $D/CTCFline_uni.bed.gz --bi $D/CTCFline_bi.bed.gz $CAP --pas $PAS \
  --annot-sites $V/sites.tsv --out $O/sites_itss.tsv 2>&1 | grep -E "BI\]|PEAKS|iTSS|OUT"
python scripts/11_evidence.py --sites $O/sites_itss.tsv --host-tss $V/host_tss.tsv $CAP --pas $PAS --out-prefix $O/itss \
  | grep -E "promotory|kontrole|kandydaci|klasa|G1|promoters|controls|candidates|Candidate genes|class [BC-]:"
echo "==== 4b. reproducibility: RAD21 line ===="
CAP2="--cap-plus $D/RAD21line_plus.bw --cap-minus $D/RAD21line_minus.bw"
python scripts/12_itss.py --gtf $G --uni $D/RAD21line_uni.bed.gz --bi $D/RAD21line_bi.bed.gz $CAP2 --pas $PAS \
  --annot-sites $V/sites.tsv --out $O/sites_itss_rad21.tsv 2>&1 | grep -E "iTSS"
python - <<'PY'
import pandas as pd, numpy as np
a = pd.read_csv("results/hct116/sites_itss.tsv", sep="\t"); b = pd.read_csv("results/hct116/sites_itss_rad21.tsv", sep="\t")
a = a[a.site_type == "candidate"]; b = b[b.site_type == "candidate"]
bk = {}
for r in b.itertuples(): bk.setdefault((r.chrom, r.strand), []).append(r.pos)
for k in bk: bk[k] = np.sort(bk[k])
hit = 0
for r in a.itertuples():
    v = bk.get((r.chrom, r.strand)); 
    if v is not None and len(v) and np.min(np.abs(v - r.pos)) <= 100: hit += 1
print(f"iTSS CTCF line: {len(a)}, RAD21 line: {len(b)}; shared (within 100 nt): {hit} ({hit / len(a) * 100:.0f}% of CTCF)")
PY
echo "==== 5. Result 4: chromatin (without DRB) ===="
python - <<'PY'
import pandas as pd
e = pd.read_csv("results/hct116/itss.evidence.tsv", sep="\t", low_memory=False)
e["eligible"] = False; e["testable"] = False; e["gap_index"] = float("nan")
e.to_csv("results/hct116/itss.nodrb.tsv", sep="\t", index=False)
PY
python scripts/21_describe.py --drb $O/itss.nodrb.tsv \
  --peaks H3K4me3=$D/H3K4me3.bed.gz H3K27ac=$D/H3K27ac.bed.gz H3K4me1=$D/H3K4me1.bed.gz --out-prefix $O/itss 2>/dev/null
echo "==== 6. Results 1 and 4: stable RNA ===="
python scripts/22_expression.py --gtf $G --cand-tx $V/candidates_tx.tsv --cand-genes $O/annot.genes.tsv \
  --sites $O/itss.describe.tsv --junctions $D/RNA.junctions.tsv --plus $D/RNA_plus.bw --minus $D/RNA_minus.bw \
  --out-prefix $O/expr 2>&1 | grep -v "^\[GTF\]"
echo "==== 7. stable RNA at matched start strength ===="
python - <<'PY'
import pandas as pd, numpy as np
e = pd.read_csv("results/hct116/expr.itss_expression.tsv", sep="\t", low_memory=False)
for c in ("host_active", "tss_supported"): e[c] = e[c].astype(str) == "True"
ctrl = e[(e.site_type == "random_control") & e.host_active]; thr = np.nanpercentile(ctrl.rna_step, 95)
x = e[e.host_active & e.tss_supported & e.site_type.isin(["candidate", "internal_promoter"])].copy()
edges = np.percentile(x.cap_sum, [0, 33.3, 66.7, 100])
x["sila"] = pd.cut(x.cap_sum, edges, include_lowest=True, labels=["1 weak", "2 medium", "3 strong"])
print(f"control threshold {thr:.2f}; PRO-cap tertile boundaries {edges.round(0)}")
for (s, t), g in x.groupby(["sila", "site_type"], observed=True):
    print(f"{s:<10} {t:<18} n={len(g):>4}  stable RNA > threshold {(g.rna_step > thr).mean() * 100:5.1f}%")
PY
echo "==== 8. Result 4: phyloP ===="
python - <<'PY'
import pandas as pd
d = pd.read_csv("results/hct116/itss.describe.tsv", sep="\t", low_memory=False).reset_index()
pd.DataFrame({"c": d.chrom, "s": (d.pos.astype(int) - 100).clip(lower=0), "e": d.pos.astype(int) + 101,
              "n": "r" + d["index"].astype(str)}).to_csv("results/hct116/sites_pm100.bed", sep="\t", header=False, index=False)
PY
T=${TMPDIR:-/tmp}/istart_$SLURM_JOB_ID; mkdir -p $T
cp ref/phylop/hg38.phyloP100way.bw $T/
bigWigAverageOverBed $T/hg38.phyloP100way.bw $O/sites_pm100.bed $O/sites_phylop.tab
rm -rf $T
sed 's#results/modA/#results/hct116/#g' scripts/helpers/phylop_summary.py | python -
echo "==== 9. Result 5: U1-PAS code ===="
python scripts/25_u1pas.py --sites $O/itss.describe.tsv --fasta ref/GRCh38.primary_assembly.genome.fa \
  --expr $O/expr.itss_expression.tsv --out $O/u1pas.tsv
echo "==== 10. Result 6: sense-antisense asymmetry ===="
python scripts/26_directionality.py --gtf $G --bi $D/CTCFline_bi.bed.gz $CAP --host-tss $V/host_tss.tsv \
  --bru-plus $D/RNA_plus.bw --bru-minus $D/RNA_minus.bw --bru-flip no --h3k36me3 $D/H3K36me3.bed.gz \
  --out $O/directionality.tsv 2>&1 | grep -v "^\[GTF\]\|^\[HOST\]\|^\[PEAKS\]"
echo "==== done ===="
