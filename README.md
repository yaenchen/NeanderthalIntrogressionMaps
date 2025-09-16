---
title: "NeanderthalIntrogressionMaps.md"
author: "Yaen Chen"
date: "2025-08-11"
output: github_document
---
# Code and data for "Comparing genomic Neanderthal introgression maps reveals core agreement but substantial heterogeneity"
   
## Overview

This repository contains data and scripts used in the Neanderthal introgression maps project. Please see the associated manuscript: "Comparing genomic Neanderthal introgression maps reveals core agreement but substantial heterogeneity"

The data directory contains data generated in the analysis, including introgression maps and deserts.

The bin directory contains scripts and Jupyter Notebooks used in the analysis for data generation (data_preparation.ipynb), data analysis, and data visualization (figures_and_analysis.ipynb, phenotypes_analysis.Rmd). 

The introgression_tools directory contains a custom tools.py containing functions used throughout the analysis.

Do not hesitate to reach out with questions or comments: yaen.chen@ucsf.edu and tony@capralab.org.


## Data Structure

Project directory: 

```bash
project_directory # The working directory
├── data # see README for additional file descriptions within each folder
│   ├── hg19
│   ├── hg38
│   └── T2T
├── introgression_tools
│   └── tools.py # custom package containing functions written for analysis
└── bin 
    ├── bash_scripts # custom scripts used in analysis
    ├── data_preparation.ipynb # data cleaning notebook
    ├── figures_and_analysis.ipynb # notebook for analysis and figure generation in the manuscript
    └── phenotypes_analysis.Rmd # R markdown file for analysis and figure generation in the manuscripts (phenotype enrichments)
```
