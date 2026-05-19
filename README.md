# Cattle_Analysis

# Gut-Protective metabolic phenotype for diarrhoea remission caused by an environmental probiotic thermophile

## Overview

This repository contains analysis scripts, datasets, and
visualisation code accompanying the manuscript submitted
to *Microbiome*.

## Repository Structure

| # | Content | Format | Description |
|---|---------|--------|-------------|
| 1 | Fig 1  | zip | Diarrhoea_daily analysis and Time course series for SCFA etc. |
| 2 | Fig 2a | zip | Bacterial α diversity analysis |
| 3 | Fig 2b | zip | Bacterial β diversity (NMDS + PERMANOVA) |
| 4 | Fig 3a (incl. Fig S6) | zip | Metabolite statistical analysis |
| 5 | Fig 3b | zip | Correlation analysis |
| 6 | Fig 4a | zip | Cliff's delta bubble chart (Metabolite feature selection) |
| 7 | Fig 4bc | zip | Individual-level fold change heatmap |
| 8 | Fig 5 (incl. Table S2, S3, S4) | md | Genomic and proteomic analysis  | Butyrate cluster visualisation |
| 9 | Fig 5a | py | Butyrate cluster visualisation |
| 10 | Fig S4 | zip | Microbial statistical analysis |
| 11 | Fig S6 (incl. Fig.4)| zip | Permutation test (Metabolite feature selection) |
| 12 | Fig S7 (incl. Table S1) | zip | PPICRUSt2 functional prediction of EC numbers |

Each zip file contains both the dataset and the analysis
script/command used to generate the corresponding figure.

## Requirements

- Python 3.10+
  - pandas, numpy, scipy, scikit-learn, matplotlib, openpyxl
- R 5.1
  - scale, vegan, pairwiseAdonis, corrplot, ggplot2
- NCBI BLAST+ (for genomic analysis)

## Usage

Each zip folder is self-contained with both the dataset
and the analysis script. To reproduce a figure:

1. Extract the zip file
2. Open the script file (.py or .R) in the extracted folder
3. Run the script 

## Genome Data

Caldifermentibacillus hisashii N11 genome: [AP028807.1](https://www.ncbi.nlm.nih.gov/nuccore/AP028807.1)

## License

Please check the license file on this website.

## Acknowledgments

- NCBI for genome sequence and annotation
- antiSMASH developers for BGC prediction tools
- Analysis scripts in this repository were developed with
  the assistance of Claude (Anthropic). All code was
  critically reviewed, tested, and validated by the authors.
