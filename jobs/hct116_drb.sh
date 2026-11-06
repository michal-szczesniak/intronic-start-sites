#!/bin/bash
#SBATCH --job-name=istart_ttdrb
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=36:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
{ source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate envs/istart; } || export PATH=$PWD/envs/istart/bin:$PATH
D=data/hct116_drb; O=results/hct116; mkdir -p $D/fq $D/qc $D/sra $D/tmp $O
# grep/sed filters below accept both the original Polish wording (logs and summaries written by the
# Polish version of the scripts) and the current English wording.
echo "==== 1. samples GSE148844 ===="
python scripts/helpers/geo_runs.py GSE148844 $D/samples.tsv | grep -E "^\s+\*|wybranych|selected runs"
# samples.tsv column 3 (selected): "tak" = yes, "nie" = no (values written by geo_runs.py)
N=$(awk -F'\t' 'NR>1 && $3=="tak" && $6!="NA"' $D/samples.tsv | wc -l)
if [ "$N" -lt 2 ] || [ "$N" -gt 30 ]; then echo "STOP: $N runs selected - check the sample titles"; exit 1; fi
echo "==== 2. download (ENA) and mapping (R1) ===="
awk -F'\t' 'NR>1 && $3=="tak" && $6!="NA" {print $6"\t"$4}' $D/samples.tsv | while IFS=$'\t' read run tmin; do
  if [ -s $D/$run.bam ]; then echo "$run: already done"; continue; fi
  echo "== $run (time $tmin min) =="
  U=https://ftp.sra.ebi.ac.uk/vol1/fastq/${run:0:6}/0${run: -2}/$run
  wget -nv -c -O $D/fq/$run.fastq.gz $U/$run.fastq.gz < /dev/null || wget -nv -c -O $D/fq/$run.fastq.gz $U/${run}_1.fastq.gz < /dev/null
  R1=$D/fq/$run.fastq.gz
  fastp -i $R1 -o $D/fq/$run.trim.fq.gz -w 16 -l 25 -j $D/qc/$run.fastp.json -h $D/qc/$run.fastp.html \
    2> $D/qc/$run.fastp.log < /dev/null
  hisat2 -p 16 -x ref/hisat2/grch38 -U $D/fq/$run.trim.fq.gz --no-unal 2> $D/qc/$run.hisat2.log < /dev/null \
    | samtools view -@ 4 -b -q 10 - | samtools sort -@ 8 -m 1G -o $D/$run.bam -
  samtools index $D/$run.bam
  rm -rf $D/fq/${run}* $D/sra/$run
  echo "$run: $(grep 'overall alignment rate' $D/qc/$run.hisat2.log); MAPQ>=10: $(samtools view -c $D/$run.bam)"
done
echo "==== 3. time points and bigWigs ===="
for t in $(awk -F'\t' 'NR>1 && $3=="tak" && $6!="NA" {print $4}' $D/samples.tsv | sort -n | uniq); do
  runs=$(awk -F'\t' -v t=$t -v d=$D 'NR>1 && $3=="tak" && $6!="NA" && $4==t {print d"/"$6".bam"}' $D/samples.tsv)
  samtools merge -f -@ 8 $D/TT_${t}min.bam $runs
  samtools index $D/TT_${t}min.bam
  for ST in + -; do
    NM=$([ "$ST" = "+" ] && echo plus || echo minus)
    bedtools genomecov -ibam $D/TT_${t}min.bam -bg -split -strand $ST | LC_ALL=C sort -k1,1 -k2,2n > $D/TT_${t}min.$NM.bg
    bedGraphToBigWig $D/TT_${t}min.$NM.bg ref/grch38.chrom.sizes $D/TT_${t}min.read_$NM.bw
    rm $D/TT_${t}min.$NM.bg
  done
  echo "TT_${t}min: $(samtools view -c $D/TT_${t}min.bam) uniquely mapped reads"
done
TIMES=$(ls $D/TT_*min.read_plus.bw | sed 's/.*TT_\(-*[0-9]*\)min.*/\1/' | sort -n | tr '\n' ' ')
REF=$(echo $TIMES | awk '{print $NF}')
echo "time points: $TIMES; reference (latest): $REF min"
echo "==== 4. wave check: ratio to the reference vs distance from the gene start ===="
python - $D $REF $TIMES <<'PY'
import sys, numpy as np, pandas as pd, pyBigWig
D, ref, times = sys.argv[1], sys.argv[2], [t for t in sys.argv[3:] if t != sys.argv[2]]
hg = pd.read_csv("results/v2/host_genes.tsv", sep="\t")
cap = {"+": pyBigWig.open("data/hct116/CTCFline_plus.bw"), "-": pyBigWig.open("data/hct116/CTCFline_minus.bw")}
def vals(bw, c, s, e):
    try: v = bw.values(c, max(0, s), e, numpy=True); return np.nan_to_num(np.abs(v)) if len(v) == e - max(0, s) else np.zeros(e - max(0, s))
    except Exception: return np.zeros(max(0, e - max(0, s)))
L, step = 120000, 4000
long = hg[(hg.end - hg.start) >= L].copy()
long["sig"] = [vals(cap[r.strand], r.chrom, r.ref_tss - 100, r.ref_tss + 101).sum() for r in long.itertuples()]
top = long.sort_values("sig", ascending=False).head(500)
def bw(t, k): return pyBigWig.open(f"{D}/TT_{t}min.read_{k}.bw")
def prof(t, same):
    b = {k: bw(t, k) for k in ("plus", "minus")}; out = np.zeros(L // step)
    for r in top.itertuples():
        k = "plus" if r.strand == "+" else "minus"
        if not same: k = "minus" if k == "plus" else "plus"
        s, e = (r.ref_tss, r.ref_tss + L) if r.strand == "+" else (r.ref_tss - L + 1, r.ref_tss + 1)
        v = vals(b[k], r.chrom, s, e)
        if r.strand == "-": v = v[::-1]
        v = np.pad(v, (0, L - len(v)))[:L]; out += v.reshape(-1, step).mean(1)
    return out
same = prof(ref, True).sum() >= prof(ref, False).sum()
print(f"read orientation: {'RNA strand' if same else 'opposite strand (swapped)'}")
pr = prof(ref, same)
print("kb\t" + "\t".join(f"{t}/{ref}" for t in times))
P = {t: prof(t, same) / np.maximum(pr, 1e-9) for t in times}
for i in range(L // step):
    print(f"{i * 4}-{i * 4 + 4}\t" + "\t".join(f"{P[t][i]:.2f}" for t in times))
PY
echo "==== 5. gap test: 10 min (start sites >= 40 kb from the host start) and 20 min (>= 70 kb) ===="
for spec in "10 40000" "20 70000"; do
  set -- $spec
  [ -s $D/TT_$1min.read_plus.bw ] || { echo "no $1 min time point"; continue; }
  echo "== $1 min vs $REF min =="
  python scripts/20_drb_gap.py --evidence $O/itss.evidence.tsv \
    --drb-plus $D/TT_$1min.read_plus.bw --drb-minus $D/TT_$1min.read_minus.bw \
    --bru-plus $D/TT_${REF}min.read_plus.bw --bru-minus $D/TT_${REF}min.read_minus.bw \
    --min-dist $2 --win 8000 --gap 5000 --out-prefix $O/ttdrb_$1 \
    | sed -n -E '/Orientacja|Orientation/p;/Kalibracja empiryczna|Empirical calibration/,/szacunek|estimate/p;/skok|step/p'
done
echo "==== 6. comparison at matched start strength (10 min) ===="
python - <<'PY'
import pandas as pd, numpy as np
d = pd.read_csv("results/hct116/ttdrb_10.drb.tsv", sep="\t", low_memory=False)
for c in ("host_active", "tss_supported", "testable"): d[c] = d[c].astype(str) == "True"
ctrl = d[(d.site_type == "random_control") & d.host_active & d.testable]
thr = np.percentile(ctrl.gap_index, 95)
x = d[d.host_active & d.tss_supported & d.testable & d.site_type.isin(["candidate", "internal_promoter"])].copy()
edges = np.percentile(x.cap_sum, [0, 33.3, 66.7, 100])
x["sila"] = pd.cut(x.cap_sum, edges, include_lowest=True, labels=["1 weak", "2 medium", "3 strong"])
print(f"control threshold {thr:.2f}; PRO-cap tertile boundaries {edges.round(0)}")
for (s, t), g in x.groupby(["sila", "site_type"], observed=True):
    print(f"{s:<10} {t:<18} n={len(g):>4}  DRB > threshold {(g.gap_index > thr).mean() * 100:5.1f}%  median {g.gap_index.median():5.2f}")
PY
echo "==== done ===="
