library(ggplot2)
library(gridExtra)
library(readr)
library(dplyr)

# Read the tcga data
Full_mut_exp_data_TCGA <- read_csv(".../Full_mut_exp_data_TCGA.csv")

# Convert to numeric
Full_mut_exp_data_TCGA <- Full_mut_exp_data_TCGA %>%
  mutate(PREOPERATIVE_PSA = as.numeric(PREOPERATIVE_PSA),
         Expression = as.numeric(Expression))

# Replace NA values with "Non-altered BRCA1/2" in the data frame
Full_mut_exp_data_TCGA$Hugo_Symbol[is.na(Full_mut_exp_data_TCGA$Hugo_Symbol)] <- "Non-altered BRCA1/2"

# Subset data for BRCA1 and BRCA2
brca_data <- Full_mut_exp_data_TCGA %>%
  filter(Hugo_Symbol %in% c("BRCA1", "BRCA2"))

brca_data <- na.omit(brca_data)

# Function to create scatter plot psa and expression
create_scatter_plot <- function(data, gene) {
  p <- ggplot(data, aes(x = PREOPERATIVE_PSA, y = Expression, color = Hugo_Symbol)) +
    geom_point() +
    stat_smooth(method = "lm", se = FALSE, aes(group = 1)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = paste("(TCGA, Cell 2015)"),
         x = "PSA (ng/mL)", y = paste(gene, ": mRNA expression (RNA Seq V2 RSEM)")) +
    scale_color_manual(name = "Gene Mutation",
                       values = c("BRCA1" = "blue", "BRCA2" = "orange")) +
    theme_minimal()
  
  # Calculate R-squared and p-value
  lm_model <- lm(Expression ~ PREOPERATIVE_PSA, data = data)
  rsq <- summary(lm_model)$r.squared
  p_val <- summary(lm_model)$coefficients[2, 4]
  
  # Add annotations as legend
  p + 
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
             label = paste("R-squared:", round(rsq, 4)), size = 4) +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 2.5,
             label = paste("p-value:", format.pval(p_val, digits = 4)), size = 4)
}

# Generating plots
psa1 <- create_scatter_plot(filter(brca_data, Hugo_Symbol == "BRCA1"), "BRCA1")
psa2 <- create_scatter_plot(filter(brca_data, Hugo_Symbol == "BRCA2"), "BRCA2")

# Combining both graphs into one plot
grid.arrange(psa1, psa2, ncol = 2)

# Convert CLINICAL_GLEASON_SUM to a factor
Full_mut_exp_data_TCGA$CLINICAL_GLEASON_SUM <- factor(Full_mut_exp_data_TCGA$CLINICAL_GLEASON_SUM)

# Create the boxplot for Gleason Score and PSA Levels
ggplot(Full_mut_exp_data_TCGA, aes(x = CLINICAL_GLEASON_SUM, y = PREOPERATIVE_PSA, color = Hugo_Symbol)) +
  geom_boxplot() +
  labs(title = "(TCGA, Cell 2015)",
       x = "Clinical Gleason Sum", y = "PSA (ng/mL)") +
  scale_color_manual(name = "Gene Mutation",
                     values = c("BRCA1" = "blue", "BRCA2" = "orange", "Non-altered BRCA1/2" = "black")) +  
  theme_minimal()

# Convert Race as factor and Age as numeric
Full_mut_exp_data_TCGA$RACE <- factor(Full_mut_exp_data_TCGA$RACE)
Full_mut_exp_data_TCGA$AGE <- as.numeric(Full_mut_exp_data_TCGA$AGE)

# Create the violin plot for Race and Age
ggplot(Full_mut_exp_data_TCGA, aes(x = RACE, y = AGE)) +
  geom_violin(alpha = 0.7, color = "Grey", fill = "#87CEEB") +
  geom_point(aes(color = Status, shape = Status, size = Status), position = position_jitter(width = 0.2)) + 
  labs(title = "(TCGA, Cell 2015)",
       x = "RACE", y = "AGE") +
  scale_shape_manual(values = c("Altered" = 15, "Unaltered" = 16)) +  # Set shapes for "Altered" and "Unaltered"
  scale_color_manual(values = c("Altered" = "red", "Unaltered" = "black")) +  # Set colors for "Altered" and "Unaltered"
  scale_size_manual(values = c("Altered" = 1.9, "Unaltered" = 1.4)) +  # Set sizes for "Altered" and "Unaltered"
  theme_minimal()

# Function to create scatter plot for age and expression
create_scatter_plot2 <- function(data, gene) {
  p <- ggplot(data, aes(x = AGE, y = Expression, color = Hugo_Symbol)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = paste("(TCGA, Cell 2015)"),
         x = "AGE (Years)", y = paste(gene, ": mRNA expression (RNA Seq V2 RSEM)")) +
    scale_color_manual(name = "Gene Mutation",
                       values = c("BRCA1" = "blue", "BRCA2" = "orange")) + 
    theme_minimal()
  
  # Calculate R-squared and p-value
  lm_model <- lm(Expression ~ AGE, data = data)
  rsq <- summary(lm_model)$r.squared
  p_val <- summary(lm_model)$coefficients[2, 4]
  
  # Add annotations as legend
  p + 
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
             label = paste("R-squared:", round(rsq, 4)), size = 4) +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 2.5,
             label = paste("p-value:", format.pval(p_val, digits = 4)), size = 4)
}

# Generating plots
age1 <- create_scatter_plot2(filter(brca_data, Hugo_Symbol == "BRCA1"), "BRCA1")
age2 <- create_scatter_plot2(filter(brca_data, Hugo_Symbol == "BRCA2"), "BRCA2")

# Combining both graphs into one plot
grid.arrange(age1, age2, ncol = 2)

# Function to create scatter plot for age and psa
create_scatter_plot3 <- function(data, gene) {
  p <- ggplot(data, aes(x = AGE, y = PREOPERATIVE_PSA, color = Hugo_Symbol)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = paste("(TCGA, Cell 2015)"),
         x = "AGE (Years)", y = paste(gene, ": PSA (ng/mL)")) +
    scale_color_manual(name = "Gene Mutation",
                       values = c("BRCA1" = "blue", "BRCA2" = "orange")) + 
    theme_minimal()
  
  # Calculate R-squared and p-value
  lm_model <- lm(AGE ~ PREOPERATIVE_PSA, data = data)
  rsq <- summary(lm_model)$r.squared
  p_val <- summary(lm_model)$coefficients[2, 4]
  
  # Add annotations as legend
  p + 
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
             label = paste("R-squared:", round(rsq, 4)), size = 4) +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 2.5,
             label = paste("p-value:", format.pval(p_val, digits = 4)), size = 4)
}

# Generating plots
age_psa1 <- create_scatter_plot3(filter(brca_data, Hugo_Symbol == "BRCA1"), "BRCA1")
age_psa2 <- create_scatter_plot3(filter(brca_data, Hugo_Symbol == "BRCA2"), "BRCA2")

# Combining both graphs into one plot
grid.arrange(age_psa1, age_psa2, ncol = 2)

# Create the boxplot for mutation status and expression
plot1 <- ggplot(brca_data, aes(x = Mutation_Status, y = log(Expression), color = Hugo_Symbol)) +
  geom_boxplot() +
  labs(title = "(TCGA, Cell 2015)",
       x = "Mutation Status", y = "mRNA expression (RNA Seq V2 RSEM)") +
  scale_color_manual(name = "Gene Mutation",
                     values = c("BRCA1" = "blue", "BRCA2" = "orange")) +
  theme_minimal()

# Create the boxplot for PSA and mutation status
plot2 <- ggplot(brca_data, aes(x = Mutation_Status, y = log(PREOPERATIVE_PSA), color = Hugo_Symbol)) +
  geom_boxplot() +
  labs(title = "(TCGA, Cell 2015)",
       x = "Mutation Status", y = "PSA (ng/mL)") +
  scale_color_manual(name = "Gene Mutation",
                     values = c("BRCA1" = "blue", "BRCA2" = "orange")) +
  theme_minimal()

# Create the boxplot for mutation status and AGE
plot3 <- ggplot(brca_data, aes(x = Mutation_Status, y = AGE, color = Hugo_Symbol)) +
  geom_boxplot() +
  labs(title = "(TCGA, Cell 2015)",
       x = "Mutation Status", y = "AGE (Years)") +
  scale_color_manual(name = "Gene Mutation",
                     values = c("BRCA1" = "blue", "BRCA2" = "orange")) +
  theme_minimal()

# Combining both graphs into one plot
grid.arrange(plot1, plot2, plot3, ncol = 3)

############### DATA SET CITATION ###################

# Cancer Genome Atlas Research Network (2015). The Molecular 
# Taxonomy of Primary Prostate Cancer. Cell, 163(4), 1011–1025.
# https://doi.org/10.1016/j.cell.2015.10.025

# cBioPortal Link to all data associated with study:
# https://www.cbioportal.org/study/summary?id=prad_tcga_pub&plots_horz_selection=%7B%22selectedGeneOption%22%3A672%2C%22dataType%22%3A%22MRNA_EXPRESSION%22%2C%22selectedDataSourceOption%22%3A%22rna_seq_v2_mrna%22%2C%22logScale%22%3A%22false%22%7D&plots_vert_selection=%7B%22dataType%22%3A%22clinical_attribute%22%2C%22selectedDataSourceOption%22%3A%22PREOPERATIVE_PSA%22%2C%22logScale%22%3A%22false%22%7D&plots_coloring_selection=%7B%22colorByMutationType%22%3A%22true%22%2C%22colorByCopyNumber%22%3A%22true%22%2C%22colorBySv%22%3A%22false%22%7D


