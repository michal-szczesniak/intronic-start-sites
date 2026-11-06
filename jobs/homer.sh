#!/bin/bash
#SBATCH --job-name=istart_homer
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.log
set -eo pipefail
M=results/motifs; mkdir -p $M/preparsed
$PWD/envs/istart/bin/python - <<'PY'
import pandas as pd
for name, path in (("K562", "results/v2/itss.describe.tsv"), ("HCT116", "results/hct116/itss.describe.tsv")):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    for c in ("host_active", "tss_supported"):
        d[c] = d[c].astype(str) == "True"
    for lab, m in (("itss", (d.site_type == "candidate") & d.tss_supported & d.host_active),
                   ("bg", (d.site_type == "random_control") & d.host_active)):
        x = d[m]
        pd.DataFrame({"c": x.chrom, "s": (x.pos.astype(int) - 150).clip(lower=0), "e": x.pos.astype(int) + 150,
                      "id": x.site_id, "sc": 0, "st": x.strand}).to_csv(f"results/motifs/{name}_{lab}.bed", sep="\t", header=False, index=False)
        print(name, lab, len(x))
PY
export PATH=$PWD/envs/homer/bin:$PATH
for L in K562 HCT116; do
  echo "==== $L: intronic start sites vs random intronic positions ===="
  findMotifsGenome.pl $M/${L}_itss.bed ref/GRCh38.primary_assembly.genome.fa $M/$L -size given \
    -bg $M/${L}_bg.bed -nomotif -p 8 -preparsedDir $M/preparsed > $M/$L.log 2>&1
  head -26 $M/$L/knownResults.txt | cut -f1,3,7,9 | sed 's#/Homer##; s#/HOMER##'
done
echo "==== done ===="
