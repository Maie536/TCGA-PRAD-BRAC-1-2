# TCGA-PRAD-BRAC-1-2-
The Cancer Genome Atlas (TCGA) through cBioPortal, Prostate Adenocarcinoma data repository

# - INTRO -
Investigating the relationship between mutations in BRCA1 and BRCA2 genes and the development of prostate cancer. 
We sought to answer the pivotal question: How do the BRCA1 and BRCA2 genes influence the development and clinical manifestations of prostate cancer? 
Given the constraints imposed by the limited number of data points available for these specific genes, our findings suggest the necessity for broader research to determine the significance and reliability of the observed outcomes.


Using The Cancer Genome Atlas (TCGA) through cBioPortal we focused on the Prostate Adenocarcinoma data repository.  The data repository consisted of clinical, mutation, and gene expression data of patients diagnosed with prostate cancer. The TCGA-PRAD specifcally sequences genomic data from tumor cells taken from paitcents already diagnosed with prostate cancer, and differes from the Cancer Cell Line Encyclopedia (CCLE) where they are not grown in as cell lines. The data used is from primary tissue samples only. 

# - CODING INFORMATION -
This project runs on Rstudios version R 3.3.0+ and is coded in R using bash as a file directory. We used the `dplyr` package to wrangle the selected dataset as well as other library packages like `GenomicDataCommons`, `ggplot2`, `gridExtra`, and `readr` assist in data acquisition, wrangling, and analyzing. Additionally, the `tidyverse` collection of packages, including `tidyr` for tidying and reshaping the data, ensured its optimal structure for subsequent analyses.

The analysis primarily focused on exploring the association between BRCA1 and BRCA2 mutations, clinical variables (e.g., PSA levels, Gleason scores, age, race), and gene expression levels. Statistical tests, such as linear regression analysis, were conducted to assess the strength and significance of associations between variables. Logarithmic transformations were applied to normalize skewed distributions, and appropriate statistical models were fitted to the data. This approach aimed to provide a clear and comprehensive understanding of the patterns and correlations within the data, enhancing the visualization of our data analysis results. 


 # - DATA ANALYSIS -

Initial data analysis revealed a distant pattern for BRCA1 and BRCA2 mutant samples across the prostate cancer clinical manifestations as seen in PSA levels and Gleason scores. As seen in figure 1, BRCA2 mutations were observed to be associated with higher PSA levels stratified by Gleason scores sum compared to BRCA1 mutations, aligning with our hypothesis that BRCA2 alterations might play a more pronounced role in the pathogenesis of prostate cancer. 

Further stratification of mutation status by expression, PSA levels, and patient age revealed distinct patterns. Figure 2, Panel 1, shows that on average those with BRCA1 mutations had a higher expression level across both mutation status. A deeper analysis is required.

Utilizing linear regression analysis via the `lm` function in R, we identified a statistically significant positive correlation between PSA levels and BRCA2 gene expression (R² = 0.8746, p-value = 0.006518), indicating a potent link to prostate cancer aggressiveness. In contrast, while BRCA1 mutations also showed a positive correlation with PSA levels (R² = 0.9487), the association did not achieve statistical significance (p-value = 0.1455).


# - NOTES -

# - Source -

These files are derived from the following virtual study on the cBioPortal for Cancer Genomics:

- **Study Name**: Prostate Adenocarcinoma (TCGA, Cell 2015)
- **cBioPortal Study ID**: `prad_tcga_pub`
- **cBioPortal URL**: [https://www.cbioportal.org/study/summary?id=prad_tcga_pub](https://www.cbioportal.org/study/summary?id=prad_tcga_pub)

# - Intended Use -

These files are intended for research purposes, to facilitate the analysis of genetic alterations in prostate cancer patients and their potential impact on prognosis, treatment outcomes, or any other clinical correlation. Researchers might use these files to segregate patient data based on BRCA1/BRCA2 alteration status before conducting further genetic, molecular, or clinical analyses.

# - Citation -

When using these files for academic or research purposes, please cite the original study as follows:

- Cancer Genome Atlas Research Network (2015). The Molecular Taxonomy of Primary Prostate Cancer. Cell, 163(4), 1011–1025. https://doi.org/10.1016/j.cell.2015.10.025

Additionally, acknowledge the use of cBioPortal for sourcing patient data.

# - Disclaimer - 

The meta file and the associated patient ID lists do not contain any personal health information (PHI) or any data that can directly identify an individual patient. All data has been de-identified in accordance with applicable laws and regulations. The user assumes all responsibility for the ethical and legal use of this data.

# - BUG FIXES -



