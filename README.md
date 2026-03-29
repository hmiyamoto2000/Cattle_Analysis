# Cattle_Analysis

# Gut-Protective metabolic phenotype for diarrhoea remission caused by an environmental probiotic thermophile

## Overview

This repository contains analysis scripts, datasets, and
visualisation code accompanying the manuscript submitted
to *Microbiome*.

## Repository Structure

| # | Content | Format | Description |
|---|---------|--------|-------------|
| 1 | Fig 1 (and Fig S1) | zip | Kaplan-Meier survival analysis |
| 2 | Fig 2a | zip | α diversity analysis |
| 3 | Fig 2b | zip | β diversity (NMDS + PERMANOVA) |
| 4 | Fig 3a (incl. Fig S6) | zip | Metabolite feature selection |
| 5 | Fig 3b | zip | Cliff's delta bubble chart |
| 6 | Fig 4a | zip | ROC curve analysis |
| 7 | Fig 4b (incl. Fig S5) | zip | Statistical analysis (BH correction) |
| 8 | Fig 5a | zip | PICRUSt2 functional analysis |
| 9 | Fig 5b | zip | PICRUSt2 visualisation |
| 10 | Fig 6 | md | Genomic analysis procedure (incl. Table S1, S2) |
| 11 | Fig 6a | py | Butyrate cluster visualisation |
| 12 | Fig S4 | zip | Correlation analysis |
| 13 | Fig S7 | zip | Time-series analysis |

Each zip file contains both the dataset and the analysis
script/command used to generate the corresponding figure.

## Requirements

- Python 3.10+
  - pandas, numpy, scipy, scikit-learn, matplotlib, openpyxl
- R 5.1
  - survival, survminer, vegan, pairwiseAdonis, corrplot, ggplot2
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
