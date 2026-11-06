# Scripts

Scripts are numbered roughly in the order in which they were run. All paths are relative to the project directory.

## Analysis steps (`scripts/`)

| Script | Description |
|---|---|
| `scripts/01_encode_files.py` | Lists (and optionally downloads) the files of ENCODE experiments through the ENCODE REST API. |
| `scripts/10_candidates.py` | Candidate intronic sites on the host strand and control sets (internal promoters, random intronic positions). |
| `scripts/11_evidence.py` | Start (PRO-cap) and end (PolyASite) evidence for the sites from step 10. |
| `scripts/12_itss.py` | Intronic start sites called directly from PRO-cap peaks, with the cap test. |
| `scripts/20_drb_gap.py` | Group test of nascent RNA after DRB release for intronic start sites, internal promoters and random positions. |
| `scripts/21_describe.py` | Chromatin marks (and optionally conservation) at start sites. |
| `scripts/22_expression.py` | Stable RNA after sites from total RNA-seq; splicing of annotated intronic RNAs. |
| `scripts/23_candidates.py` | List of the strongest intronic start sites. |
| `scripts/24_exosome.py` | RNA after start sites in the ENCODE exosome knockouts (MTREX, EXOSC10, EXOSC5). |
| `scripts/25_u1pas.py` | Splice-donor sites and polyadenylation signals downstream of start sites. |
| `scripts/26_directionality.py` | Strand asymmetry of divergent start pairs inside genes and between genes. |
| `scripts/27_classify.py` | Single-site classes after DRB release from both replicates (HCT116, DRB/TTchem-seq). |
| `scripts/28_repeats.py` | Transposable elements at start sites. |
| `scripts/29_context.py` | Intronic start sites compared with start sites between genes; RNA per initiation event. |
| `scripts/30_known_genes.py` | Validation: detection of the start sites of genes nested in introns. |
| `scripts/31_final_tables.py` | Final per-site tables (Supplementary Table 1). |
| `scripts/32_loci.py` | Locus figures drawn with pyGenomeTracks. |
| `scripts/33_hm13_ipa.py` | Intronic polyadenylation of HM13 in the intron that holds MCTS2. |
| `scripts/34_supp_numbers.py` | Numbers for the supplementary figures, computed from the result tables. |
| `scripts/35_exon_distance.py` | Distance from intronic start sites to the next host exon, and initiation and transcript fate by distance. |
| `scripts/36_stable_far.py` | Stable RNA measured only for sites at least 3 kb upstream of the next host exon. |
| `scripts/37_stable_far_rest.py` | The 3-kb rule applied to the remaining results based on stable RNA. |
| `scripts/38_fig4a_fig5a.py` | Exosome knockout and RNA per initiation event with the 3-kb rule. |
| `scripts/39_fig5a.py` | Stable RNA per initiation event with the 3-kb rule, from the context tables. |
| `scripts/40_psi_both_lines.py` | Inclusion of the next host exon and initiation at intronic start sites, in both cell lines. |
| `scripts/41_register_checks.py` | Values reported in the text that no earlier step printed. |
| `scripts/42_extra_checks.py` | Additional controls reported in the text. |
| `scripts/43_fig6_panels.py` | Locus panels of single start sites (pyGenomeTracks). |
| `scripts/44_or_ci.py` | Confidence intervals of the Mantel–Haenszel odds ratios for shared and cell-type-specific sites. |
| `scripts/45_depth_matched_deficit.py` | Depth-matched deficit of the arm on the host's strand in divergent pairs, compared across definitions. |
| `scripts/45b_depth_matched_nn.py` | Depth-matched deficit of the arm on the host's strand with nearest-neighbour matching, for active and silent genes. |
| `scripts/50_figure_data.py` | Tables of the plotted values, one per figure panel, from the result tables. |
| `scripts/51_figures.py` | Figures drawn with matplotlib from the tables of step 50. |

## Helpers (`scripts/helpers/`)

| Script | Description |
|---|---|
| `scripts/helpers/encode_get.py` | Downloads processed ENCODE files (GRCh38) for a cell line. |
| `scripts/helpers/encode_rnaseq_bw.py` | Downloads the RNA-seq bigWigs (unique reads, GRCh38) of an ENCODE experiment and its control and merges the replicates. |
| `scripts/helpers/geo_runs.py` | Lists the samples of a GEO series with their SRA runs, selected by a pattern in the sample title. |
| `scripts/helpers/merge_bw.py` | Sums several bigWigs (for example replicates) into one, using absolute values. |
| `scripts/helpers/phylop_summary.py` | Summary of phyloP conservation at start sites. |
| `scripts/helpers/starts_breakdown.py` | Genome-wide breakdown of PRO-cap start sites by position relative to annotated genes. |

## SLURM jobs (`jobs/`)

| Job | Description |
|---|---|
| `jobs/etap0_genome.sh` | First genome-wide run of steps 10–12 (K562). |
| `jobs/core_v2.sh` | Main K562 run: steps 10–12 and 20–23, with phyloP conservation. |
| `jobs/hct116.sh` | Main HCT116 run: downloads of the ENCODE files, steps 11, 12, 21, 22, 25 and 26, with phyloP conservation. |
| `jobs/hct116_drb.sh` | Download and processing of the HCT116 DRB/TTchem-seq time course (GSE148844). |
| `jobs/map_brudrb.sh` | Download and mapping of K562 BruDRB-seq and Bru-seq (GSE55534). |
| `jobs/hisat2_index.sh` | HISAT2 index of the genome. |
| `jobs/rnaseq.sh` | Download of ENCODE total RNA-seq alignments and signal files, merging of replicates. |
| `jobs/exosome.sh` | Exosome knockout RNA-seq (ENCODE) and step 24. |
| `jobs/classify.sh` | Single-site classes after DRB release (step 27). |
| `jobs/context.sh` | Start sites between genes and RNA per initiation event (step 29). |
| `jobs/k562_directionality.sh` | Strand asymmetry of divergent pairs in K562 (step 26). |
| `jobs/homer.sh` | Known transcription factor motifs at start sites (HOMER). |
| `jobs/ltr12c_gata.sh` | GATA, NF-Y and other motifs in LTR12C copies with and without an intronic start site. |
| `jobs/phylop.sh` | phyloP conservation at start sites. |
| `jobs/known_genes.sh` | Validation on genes nested in introns (step 30). |
| `jobs/known_genes_wide.sh` | Validation on genes nested in introns, both strands and all host transcripts (step 30). |
| `jobs/loci.sh` | Locus figures (step 32). |
| `jobs/loci_final.sh` | Final locus figures (step 32). |
| `jobs/robust.sh` | Earlier robustness run with a 500-nt distance rule (steps 10–12, 20, 21). |
| `jobs/robust_v2.sh` | Robustness run of the main K562 analysis with a 500-nt distance rule. |
| `jobs/supp_numbers.sh` | Numbers for the supplementary figures (step 34). |
