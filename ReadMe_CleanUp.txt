# TCGA Prostate Cancer Data Analysis

This project aims to analyze The Cancer Genome Atlas (TCGA) data for prostate cancer, focusing on the relationship between BRCA1 and BRCA2 gene mutations and expression levels, and clinical outcomes.

## Getting Started

These instructions will guide you through setting up your project and cleaning up the data needed for analysis.

## Prerequisites

Before running this script, ensure you have R installed on your machine. You will also need to install the following R packages:
- `tidyverse`
- `dplyr`
- `ggplot2`

You can install these packages using the following R commands:

```r
install.packages("tidyverse")
install.packages("dplyr")
install.packages("ggplot2")

## Data Files
Make sure the following data files are placed in the specified directory:

data_clinical_patient.txt
data_clinical_sample.txt
data_mrna_seq_v2_rsem.txt
data_mutations.txt

These two files were obtained directly from cBioPortal virtual tool and will be added to the compiled list of files with the script.

altered_samples.txt
unaltered_samples.txt

Running the clean up Script
Open the R script in your R IDE or text editor.
Make sure to set your working directory to the folder containing your data files and update the file path on the script.
Run the script. It will perform data preprocessing, merge clinical and genomic data, and output combined datasets for further analysis.

## Output

The script will generate CSV files with merged clinical and genetic data for BRCA1 and BRCA2, categorized into altered and unaltered gene expressions:

Full_mut_data_TCGA.csv
unalt_brca1_exp_data_TCGA.csv
unalt_brca2_exp_data_TCGA.csv
Full_mut_exp_data_TCGA.csv

These files will be saved to specified paths you updated

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
