# Identifying Transcriptional Differences in Lung Cancer
20.440 Project (Cora Wendlandt and Rachel McGinn)

![Annotated Volcano Plot](https://github.com/CoraWendlandt/SexTranscriptionLungCancer/blob/1bd106da799adb20fd20d36f7fd25f1fb730a2d2/figures/ALLannotated_volcano.png)

To redo the creation of the figure, run volcano.py and then annotated_volcano.py. If you do not care about recreating the data parsing and just want the plot, run annotated_volcano.py only.
## Overview
Lung cancer affects people from a variety of backgrounds, but research has suggested that incidence and outcomes differ on the basis of sex. Specifically, [Stabellini et al. (2022)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8983352/) found that males have a higher incidence of lung cancer and that female sex was associated with higher surgical rates, lower immunotherapy use rates, higher rates of endocrinologic complications after immunotherapy use, and higher rates of psychological disorders. Accordingly, here we explore transcriptional differences in male and female lung cancer samples hoping to answer the following:
### What (if any) are the expression differences in genes for male versus female lung cancer patients?
Such findings would offer hypotheses as to the origin of differences in incidence and outcome across sex.

## Data
Data was procured from the [Genetic Data Commons Data Portal](https://portal.gdc.cancer.gov/). Specifically, we look at 3 RNAseq datasets that explored bronchus and lung cancer: [CPTAC-3](https://portal.gdc.cancer.gov/projects/CPTAC-3), [TCGA-LUAD](https://portal.gdc.cancer.gov/projects/TCGA-LUAD), and [TCGA-LUSC](https://portal.gdc.cancer.gov/projects/TCGA-LUSC). This gives 843 male cases and 512 female cases which we used for downstream analysis.

## Folder Structure
### Code
Holds the scripts to reproduce all figures and data analysis
### Data
Holds the data used for all downstream analysis
### Figures
Holds the figures from our analysis

## Installation
### Packages and Software
All the scripts were run using Python 3.7.5 and required packages can be found in code/requirements.txt
### Running the code
Clone the repo, unzip the tar.gz files in data, run scripts following the instructions commented at their top
