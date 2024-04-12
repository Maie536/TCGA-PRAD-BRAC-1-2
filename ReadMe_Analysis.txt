# TCGA Prostate Cancer Data Analysis

This script performs an analysis of TCGA (The Cancer Genome Atlas) data, focusing on the relationship between various genetic expressions and clinical attributes in prostate cancer samples. It utilizes data from the TCGA study published in Cell, 2015, to explore mRNA expressions of BRCA1 and BRCA2 genes, preoperative PSA levels, clinical Gleason scores, patient age, and race.

## Getting Started

These instructions will guide you through setting up your project and cleaning up the data needed for analysis.

## Prerequisites

Before you run this script, you need to have R installed on your computer along with the following R packages:
- `ggplot2`
- `gridExtra`
- `readr`
- `dplyr`

You can install these packages using R commands like:

```r
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("readr")
install.packages("dplyr")

## Dataset

The script analyzes data from the following study:

Cancer Genome Atlas Research Network (2015). The Molecular Taxonomy of Primary Prostate Cancer. Cell, 163(4), 1011–1025. https://doi.org/10.1016/j.cell.2015.10.025
Ensure you have ran the TCGA_Data_CleanUp to obtain the following file: full_mut_exp_data_TCGA.csv and update the path as necessary for your environment.

## Features

This script includes functions and visualizations for:

Converting specific dataset columns to numeric for analysis.
Filtering the dataset for BRCA1 and BRCA2 gene expressions.
Creating scatter plots and box plots to explore the relationship between PSA levels, mRNA expression, clinical Gleason scores, patient age, and mutation status.
Comparing gene expression levels in mutated versus unmutated samples.
Analyzing the distribution of age and race within the dataset.
Usage

To run this script, simply load it in your R environment and execute. Make sure the dataset path matches your local setup.

## Visualization

The script generates several plots, including:

Scatter plots for BRCA1 and BRCA2 gene expressions against preoperative PSA levels.
Box plots comparing clinical Gleason scores against PSA levels.
Violin plots for patient age distribution across different races.
Combined graphs for quick comparative analysis.
Citations

For citations and further reading, refer to the TCGA 2015 study linked above and additional resources provided within the script.

## Citation

Please cite the following if you use this dataset for your research:

Cancer Genome Atlas Research Network (2015). The Molecular Taxonomy of Primary Prostate Cancer. Cell, 163(4), 1011–1025. https://doi.org/10.1016/j.cell.2015.10.025

## Data Source

All data associated with this study can be found at cBioPortal.
Link: https://www.cbioportal.org/study/summary?id=prad_tcga_pub

## Feedback

For questions, suggestions, or feedback, please open an issue in this repository.

## License

This script is provided for educational and research use. Any commercial use of the script or the data it processes should comply with the terms provided by the data source.