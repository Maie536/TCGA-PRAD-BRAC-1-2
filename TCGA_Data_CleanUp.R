# Make sure to install required packages prior to running script
library(tidyverse)
library(dplyr)
library(ggplot2) 

#______________________Step 1a: Read in clinical and patient data files_______________________________

# Read the txt file of data_clinical_patient and data_clinical_sample for each data set
clin_data_TCGA <- read_tsv("/Users/benitaariaslara/Desktop/lsc586/R/TCGA/data/data_clinical_patient.txt")
sample_data_TCGA <- read_tsv("/Users/benitaariaslara/Desktop/lsc586/R/TCGA/data/data_clinical_sample.txt")

#________________________Step 1b: Clean up data to merge clinical files____________________________

# Set column names for clin_data using row 4
colnames(clin_data_TCGA) <- as.character(unlist(clin_data_TCGA[4, ]))

# Set column names for sample_data using row 4
colnames(sample_data_TCGA) <- as.character(unlist(sample_data_TCGA[4, ]))

# Remove the first 3 rows as they are no longer needed
clin_data_TCGA <- clin_data_TCGA[-(1:4), ]
sample_data_TCGA <- sample_data_TCGA[-(1:4), ]

# Merge the two datasets based on Patient Identifier
merged_data_TCGA <- merge(clin_data_TCGA, sample_data_TCGA, by = "PATIENT_ID")
head(merged_data_TCGA)

# Select only the desired columns
selected_columns_TCGA <- c("PATIENT_ID", "SAMPLE_ID", "AGE","RACE","PREOPERATIVE_PSA","CLINICAL_GLEASON_SUM")
merged_data_TCGA <- merged_data_TCGA[, selected_columns_TCGA]

#_____________________________Step 1c: Read mutation and expression data____________________________________

# Read the txt file of z-scores and raw expression values
raw_data_TCGA <- read_tsv("/Users/benitaariaslara/Desktop/lsc586/R/TCGA/data/data_mrna_seq_v2_rsem.txt")
data_mutations_TCGA <- read_tsv("/Users/benitaariaslara/Desktop/lsc586/R/TCGA/data/data_mutations.txt")

# extract sample ID and Hugo Symbol from mutations table
columns_mutations <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Mutation_Status","Variant_Classification")
data_mutations_TCGA <- data_mutations_TCGA[, columns_mutations]

# rotate the expression tables
raw_data_TCGA_long <- raw_data_TCGA %>%
  pivot_longer(cols = -c(Hugo_Symbol, Entrez_Gene_Id), 
               names_to = "SAMPLE_ID", 
               values_to = "Expression") %>%
  select(-Entrez_Gene_Id)

# Filter for BRCA1 and BRCA2 for unaltered
brca1_data_ <- raw_data_TCGA_long %>%
  filter(Hugo_Symbol %in% c("BRCA1"))
brca2_data_ <- raw_data_TCGA_long %>%
  filter(Hugo_Symbol %in% c("BRCA2"))
brca1_data_ <- brca1_data_ %>%
  select(-Hugo_Symbol)
brca2_data_ <- brca2_data_ %>%
  select(-Hugo_Symbol)

# Filter for brca121 for mutations
brca1_mut_data <- data_mutations_TCGA %>%
  filter(Hugo_Symbol %in% c("BRCA1"))
brca2_mut_data <- data_mutations_TCGA %>%
  filter(Hugo_Symbol %in% c("BRCA2"))

brca1_data <- raw_data_TCGA_long %>%
  filter(Hugo_Symbol %in% c("BRCA1"))
brca2_data <- raw_data_TCGA_long %>%
  filter(Hugo_Symbol %in% c("BRCA2"))

brca1_mut_data <- data_mutations_TCGA %>%
  filter(Hugo_Symbol %in% c("BRCA1"))
brca2_mut_data <- data_mutations_TCGA %>%
  filter(Hugo_Symbol %in% c("BRCA2"))

# Read in file with Altered Patient IDs
altered_data_TCGA <- read_delim("/Users/benitaariaslara/Desktop/lsc586/R/TCGA/altered_samples_TCGA.txt", delim = ":", col_names = FALSE)

# Filter expression data using altered IDs
brca1_data <- filter(brca1_data, `SAMPLE_ID` %in% altered_data_TCGA$X2)
brca2_data <- filter(brca2_data, `SAMPLE_ID` %in% altered_data_TCGA$X2)
brca1_data <- brca1_data %>%
  select(-Hugo_Symbol)
brca2_data <- brca2_data %>%
  select(-Hugo_Symbol)

# Rename Tumor_Sample_Barcode to Sample_ID
brca1_mut_data <- brca1_mut_data %>% rename(SAMPLE_ID = Tumor_Sample_Barcode)
brca2_mut_data <- brca2_mut_data %>% rename(SAMPLE_ID = Tumor_Sample_Barcode)

# Merge expression and mutation data with clinical data
brca1_mut_data_TCGA <- merge(merged_data_TCGA, brca1_mut_data, by = "SAMPLE_ID")
brca2_mut_data_TCGA <- merge(merged_data_TCGA, brca2_mut_data, by = "SAMPLE_ID")

# add col to say altered for mutated files
brca1_mut_data_TCGA$Status <- "Altered"
brca2_mut_data_TCGA$Status <- "Altered"
brca1_mut_data_TCGA$Study <- "Prostate Adenocarcinoma (TCGA, Cell 2015)"
brca2_mut_data_TCGA$Study <- "Prostate Adenocarcinoma (TCGA, Cell 2015)"

# combine mutation data
merged_exp_brca1_data_TCGA <- merge(brca1_mut_data_TCGA, brca1_data, by = "SAMPLE_ID")
merged_exp_brca2_data_TCGA <- merge(brca2_mut_data_TCGA, brca2_data, by = "SAMPLE_ID")

Full_mut_data_TCGA <- bind_rows(
  merged_exp_brca1_data_TCGA,
  merged_exp_brca2_data_TCGA
)

# write files into folders
write_csv(Full_mut_data_TCGA, "/Users/benitaariaslara/Desktop/lsc586/R/TCGA/Full_mut_data_TCGA.csv")

# merge brca1/2_data_ with clinical for unaltered
merged_brca1_data_ <- merge(merged_data_TCGA, brca1_data_, by = "SAMPLE_ID")
merged_brca2_data_ <- merge(merged_data_TCGA, brca2_data_, by = "SAMPLE_ID")

#_________________Step 3a: Filter by altered and unaltered, obtained from cBioPortal______________

# Read in file with unaltered Patient IDs
unaltered_data_TCGA <- read_delim("/Users/benitaariaslara/Desktop/lsc586/R/TCGA/unaltered_samples_TCGA.txt", delim = ":", col_names = FALSE)

# Filter merged data using unaltered patient IDs
unalt_brca1_exp_data_TCGA <- filter(merged_brca1_data_, `SAMPLE_ID` %in% unaltered_data_TCGA$X2)
unalt_brca2_exp_data_TCGA <- filter(merged_brca2_data_, `SAMPLE_ID` %in% unaltered_data_TCGA$X2)

#_________________________Step 3b: Merge 3 altered mutations data sets together by gene__________________________________

# Add col for Study Name
unalt_brca1_exp_data_TCGA$Study <- "Prostate Adenocarcinoma (TCGA, Cell 2015)"
unalt_brca2_exp_data_TCGA$Study <- "Prostate Adenocarcinoma (TCGA, Cell 2015)"
unalt_brca1_exp_data_TCGA$Status <- "Unaltered"
unalt_brca2_exp_data_TCGA$Status <- "Unaltered"

# write files into folders
write_csv(unalt_brca1_exp_data_TCGA, "/Users/benitaariaslara/Desktop/lsc586/R/TCGA/unalt_brca1_exp_data_TCGA.csv")
write_csv(unalt_brca2_exp_data_TCGA, "/Users/benitaariaslara/Desktop/lsc586/R/TCGA/unalt_brca2_exp_data_TCGA.csv")

# combine mut and exp data

Full_mut_exp_data_TCGA <- bind_rows(
  Full_mut_data_TCGA,
  unalt_brca1_exp_data_TCGA,
  unalt_brca2_exp_data_TCGA
)

# write files into folders
write_csv(Full_mut_exp_data_TCGA, "/Users/benitaariaslara/Desktop/lsc586/R/TCGA/Full_mut_exp_data_TCGA.csv")


############### DATA SET CITATION ###################

# Cancer Genome Atlas Research Network (2015). The Molecular 
# Taxonomy of Primary Prostate Cancer. Cell, 163(4), 1011–1025.
# https://doi.org/10.1016/j.cell.2015.10.025

# cBioPortal Link to all data associated with study:
# https://www.cbioportal.org/study/summary?id=prad_tcga_pub&plots_horz_selection=%7B%22selectedGeneOption%22%3A672%2C%22dataType%22%3A%22MRNA_EXPRESSION%22%2C%22selectedDataSourceOption%22%3A%22rna_seq_v2_mrna%22%2C%22logScale%22%3A%22false%22%7D&plots_vert_selection=%7B%22dataType%22%3A%22clinical_attribute%22%2C%22selectedDataSourceOption%22%3A%22PREOPERATIVE_PSA%22%2C%22logScale%22%3A%22false%22%7D&plots_coloring_selection=%7B%22colorByMutationType%22%3A%22true%22%2C%22colorByCopyNumber%22%3A%22true%22%2C%22colorBySv%22%3A%22false%22%7D



