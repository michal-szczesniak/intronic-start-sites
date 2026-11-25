# Intronic start sites

Code for the analyses and figures of:

> Szcześniak MW, Wanowska E. Transcription starting inside human introns is widespread but mostly yields short-lived RNA. Manuscript submitted for publication (2026).

The scripts map capped transcription start sites inside introns of human protein-coding genes from public PRO-cap data (K562 and HCT116 cells) and test whether these sites produce stable RNA (total RNA-seq), sustain transcription after DRB release (BruDRB-seq, DRB/TTchem-seq), give RNA that accumulates without the nuclear exosome (ENCODE knockouts), carry the splice-donor and polyadenylation-signal arrangement of promoters, and differ from start sites between genes.

## Contents

- `scripts/`: all analysis and figure scripts, numbered roughly in the order in which they were run; `scripts/helpers/`: download and file helpers. `SCRIPTS.md` gives a one-line description of each.
- `jobs/`: the SLURM batch jobs that ran the steps.
- `env/`: conda environments, installed Python packages and Python versions.
- Processed tables are archived in Zenodo ([10.5281/zenodo.23984488](https://doi.org/10.5281/zenodo.23984488)).

## Data

All input data are public: PRO-cap, total RNA-seq, histone ChIP-seq peaks and exosome-knockout RNA-seq from the ENCODE portal; BruDRB-seq (GEO GSE55534) and DRB/TTchem-seq (GEO GSE148844); the GRCh38 primary assembly with the GENCODE release 47 annotation; PolyASite 2.0 clusters; the UCSC phyloP (100 vertebrates) and RepeatMasker tracks; and the GTEx v8 median expression table. Accession numbers are listed in the supplementary data of the article. `scripts/01_encode_files.py` and the scripts in `scripts/helpers/` download the ENCODE and GEO files; the FASTQ files of the BruDRB-seq and Bru-seq runs (SRR1182352-SRR1182355) are expected in `data/brudrb/` (available from ENA).

## Running

The steps were run from the project directory, which holds `scripts/`, `data/`, `ref/` and `results/`, on a SLURM cluster with the jobs in `jobs/` (environment `envs/istart`; `envs/pgt` for pyGenomeTracks and the figures). Each script writes its tables under `results/`.

The figures can be redrawn from the archived tables without the raw data: extract the Zenodo data archive into this directory and run

```
python scripts/50_figure_data.py --no-bigwig --skip-gtf
python scripts/51_figures.py
```

with Python 3 and pandas, numpy, scipy and matplotlib (for example the environment in `env/pgt.yml`).

`--no-bigwig` takes the nascent RNA profiles from the tables of steps 20 and 24 instead of the signal files, and `--skip-gtf` leaves out the genome-wide breakdown of start sites, which needs the GENCODE annotation and the PRO-cap peak files.

## Figure and table numbering

The scripts of steps 01–45, 50 and 51 name the figures and tables as in the first preprint version of the article; step 46 and `SCRIPTS.md` follow the current version. In the current version, Figs. 1–5 and panels a and b of Fig. 6 keep their numbers, Fig. 6c is Fig. 7 and Fig. 7 is Fig. 8; Supp. Figs. S1–S6 are Supplementary Figs. 1–6, S7 is Supplementary Fig. 8, S8 is Supplementary Fig. 9, S9 is Supplementary Fig. 10 and S10 is Supplementary Fig. 7. Supplementary Table 1 of the scripts (all intronic start sites) is Supplementary Data 1 of the article, and Supplementary Table 1 of the article is the output of step 46.

## License

Code: MIT License (see `LICENSE`). Data in the Zenodo archive: CC BY 4.0.

## Citation

Please cite the article above and the Zenodo archive ([10.5281/zenodo.23984488](https://doi.org/10.5281/zenodo.23984488)).
