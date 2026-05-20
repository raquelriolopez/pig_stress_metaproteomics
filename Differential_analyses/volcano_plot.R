#################Volcano plot
library(tidyverse)
library(dplyr)
library(ggrepel)
library(ggplot2)
library("easyCODA")

###Inputs

metadata <- read.table("metadataMetaP.txt", header = T, sep = "\t") 

#MAGG
MAGG_abundance <- read.table("MAGG_abundance.tsv", header=T, sep="\t")  ##RAW Catalogue abundance

#MAGP
MAGP_abundance_raw <- read.table("summedLFQ_per_bin_raw_New_X.txt", header=T, sep="\t") ##RAW MAG abundance from proteins
MAGP_abundance_raw <- as.data.frame(MAGP_abundance_raw)

MAGP_abundance_log <- read.table("MAGP_abundance_log_New.txt", header=T, sep="\t") ##LOG MAG abundance from proteins
MAGP_abundance_log <- as.data.frame(MAGP_abundance_log)

#MAGP_abundance_log <- MAGP_abundance_raw
#MAGP_abundance_log[, -1] <- log2(MAGP_abundance_raw[, -1] + 1) ####to avoid log(0) we add +1
#write.table(MAGP_abundance_log, file = "MAGP_abundance_log_New.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

MAGP_abundance_clr <- read.table("MAGP_abundance_clr_New.txt", header=T, sep="\t") ##CLR MAG abundance from proteins
MAGP_abundance_clr <- as.data.frame(MAGP_abundance_clr)

#numeric_abundance <- MAGP_abundance_raw[, -1]
#MAGP_abundance_clr <- CLR(numeric_abundance + 0.0001, weight = TRUE)$LR
#MAGP_abundance_clr <- data.frame(Bin = MAGP_abundance_raw$Bin, MAGP_abundance_clr)
#write.table(MAGP_abundance_clr, file = "MAGP_abundance_clr_New.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#PROTEINS
protein_abundance_raw <- read.table("MetaP2_raw_filtfreq.tsv", header=T, sep="\t") ##RAW protein abundance
protein_abundance_raw <- as.data.frame(protein_abundance_raw)

#Just to be sure I've recalculated log, but it was correct, they were the same.
#protein_abundance_log <- log2(protein_abundance_raw + 1) ####to avoid log(0) we add +1
#write.table(protein_abundance_log, file = "MetaP2_log_filtfreq.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

protein_abundance_log <- read.table("MetaP2_log_filtfreq.tsv", header=T, sep="\t") ###LOG protein abundance
protein_abundance_log <- as.data.frame(protein_abundance_log)

#protein_abundance_clr <- CLR(protein_abundance_raw + 0.0001, weight = TRUE)$LR
#protein_abundance_clr <- as.data.frame(protein_abundance_clr)
#write.table(protein_abundance_clr, file = "MetaP2_clr_filtfreq.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

protein_abundance_clr <- read.table("MetaP2_clr_filtfreq.txt", header=T, sep="\t") ####CLR protein abundance
protein_abundance_clr <- as.data.frame(protein_abundance_clr)

#Visualisation

View(MAGG_abundance)
View(MAGP_abundance_raw)
View(MAGP_abundance_log)
View(MAGP_abundance_clr)

View(protein_abundance_raw)
View(protein_abundance_log)
View(protein_abundance_clr)

library(readxl)
MAG_taxonomy <- read_excel("taxonomy_MAG_C70C10_MetaP2.xlsx")
View(MAG_taxonomy)


###################Volcano plot for MAGG RAW DATA (relative abundance)

#Transpose the abundance table to have samples as rows and MAGs as columns
abundance_long <- MAGG_abundance %>%
  pivot_longer(cols = -Genomes, names_to = "SampleID", values_to = "Abundance")

#Combine with metadata
combined_data <- abundance_long %>%
  left_join(metadata, by = "SampleID")

View(combined_data)

#Perform statistical tests for each MAG
results <- combined_data %>%
  group_by(Genomes) %>%
  summarise(
    # Average per condition
    mean_control = mean(Abundance[Tto == "control"], na.rm = TRUE),
    mean_stress = mean(Abundance[Tto == "stress"], na.rm = TRUE),
    # fold change (en escala log2)
    log2Fold_change = log2(mean_stress / mean_control),
    # t-test
    p_value = t.test(Abundance[Tto == "stress"], Abundance[Tto == "control"])$p.value
    ) %>%
  ungroup() %>%
  mutate(
    # Calcular -log10(p-value) para el volcano plot
    SCORE_neg_log10_p = -log10(p_value),
    Significance = case_when(
      p_value < 0.05 & abs(log2Fold_change) > 1 ~ "Significant (p<0.05, |FC|>1)",
      p_value < 0.1 & abs(log2Fold_change) >3 ~ "Significant (p<0.1, |FC|>3)",
      TRUE ~ "Not significant"
    )
  )
View(results)
head(results)

# volcano plot
volcano_plot <- ggplot(results, aes(x = log2Fold_change, y = SCORE_neg_log10_p)) +
  geom_point(aes(color = Significance)) +
  scale_color_manual(values = c(
    "Significant (p<0.05, |FC|>1)" = "red",
    "Significant (p<0.1, |FC|>3)" = "orange",
    "Not significant" = "gray"
    )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.1), linetype = "dotted", color = "darkgreen") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-3, 3), linetype = "dotted", color = "darkgreen") +
  labs(x = "Log2 Fold Change (Stress/Control)",
       y = "-Log10(p-value)",
       title = "Volcano Plot MAG abundance: Control vs Stress",
       color = "Significance") +
  theme_minimal() +
  # Extreme significant points
  geom_text_repel(
    data = subset(results, Significance != "Not significant"),
    aes(label = Genomes),
    size = 3,
    box.padding = 0.5,
    max.overlaps = 20
  )

print(volcano_plot)

# p-value distribution
ggplot(results, aes(x = p_value)) + 
  geom_histogram(bins = 30) +
  labs(title = "Distribución de p-values")

# significant MAGs
cat("MAGs significativos:", sum(results$p_value < 0.05 & abs(results$log2Fold_change) > 1, na.rm = TRUE))
#6 significant MAGs


# p-value distribution and fold changes
par(mfrow = c(1, 2))
hist(results$p_value, main = "p-values ajustados")
hist(results$log2Fold_change, main = "Log2 Fold Changes")


### Generate a table of significant MAGs
significance_counts <- results %>%
  count(Significance) %>%
  mutate(Percentage = n / sum(n) * 100)
print(significance_counts)

significant_mags <- results %>%
  filter(Significance != "Not significant") %>%
  arrange(Significance, desc(abs(log2Fold_change)))
View(significant_mags)

write.csv(significant_mags, "significant_MAGs.csv", row.names = FALSE)

######Add taxonomy
library(readxl)
library(writexl)

taxonomy <- read_excel("taxonomia_MetaP2.xlsx")

significant_mags_with_taxonomy <- significant_mags %>%
  left_join(taxonomy %>% select(user_genome, classification), 
            by = c("Genomes" = "user_genome"))

View(significant_mags_with_taxonomy)
write.csv(significant_mags_with_taxonomy, "significant_mags_with_taxonomy.csv", row.names = FALSE)

write_xlsx(x = significant_mags_with_taxonomy, path = "significant_mags_with_taxonomy.xlsx", col_names = TRUE, format_headers = TRUE)

###########################Add taxonomy to the volcano plot

# add taxonomy
results_with_tax <- results %>%
  left_join(MAG_taxonomy %>% 
              select(user_genome, classification),
            by = c("Genomes" = "user_genome")) %>%
  mutate(
    # Extraer el último nivel taxonómico disponible
    last_taxon = sub(".*;(g__[^;]+).*", "\\1", classification),
    last_taxon = ifelse(grepl(";", last_taxon), 
                        sub(".*;(f__[^;]+).*", "\\1", classification),
                        last_taxon),
    last_taxon = ifelse(grepl(";", last_taxon), 
                        sub(".*;(o__[^;]+).*", "\\1", classification),
                        last_taxon),
    
    # Crear etiqueta de especie: s__ si existe, último taxón + (NS) si no
    Species_label = ifelse(grepl(";s__[^;]+", classification),
                           sub(".*;s__([^;]+)", "\\1", classification),
                           paste0(last_taxon, " (NS)")),
    
    # Versión corta para el gráfico
    Species_short = ifelse(nchar(Species_label) > 25,
                           paste0(substr(Species_label, 1, 22), "..."),
                           Species_label)
  )

# volcano plot again
volcano_plot <- ggplot(results_with_tax, aes(x = log2Fold_change, y = SCORE_neg_log10_p)) +
  geom_point(aes(color = Significance)) +
  scale_color_manual(values = c(
    "Significant (p<0.05, |FC|>1)" = "red",
    "Significant (p<0.1, |FC|>3)" = "orange",
    "Not significant" = "gray"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.1), linetype = "dotted", color = "darkgreen") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-3, 3), linetype = "dotted", color = "darkgreen") +
  labs(x = "Log2 Fold Change (Stress/Control)",
       y = "-Log10(p-value)",
       title = "Volcano Plot MAG abundance: Control vs Stress",
       color = "Significance") +
  theme_minimal() +
  # Extreme significant points with taxonomy labels
  geom_text_repel(
    data = subset(results_with_tax, Significance != "Not significant"),
    aes(label = Species_short),  # Usar la etiqueta taxonómica corta
    size = 3,
    box.padding = 0.5,
    max.overlaps = 20
  )

print(volcano_plot)












#####################################################Volcano plot for MAGs calculated from proteins (MAGP_raw, _log, _clr)
#INPUT

metadata <- read.table("metadataMetaP.txt", header = T, sep = "\t") 

#Transpose the abundance table to have samples as rows and MAGs as columns
abundance_long_magp <- MAGP_abundance_log %>%    #Change here MAGP_abundance_log or _clr (CLR hay que tratar los datos diferente, repensar)
  pivot_longer(cols = -Bin, names_to = "SampleID", values_to = "Abundance")  %>%
  mutate(Abundance = if_else(Abundance == 0, rnorm(n(), mean = 8, sd = 1), as.numeric(Abundance))) #MAGNUS correction in log

View(abundance_long_magp)

#Combine with metadata
combined_data_magp <- abundance_long_magp %>%
  left_join(metadata, by = "SampleID")

View(combined_data_magp)

#Perform statistical tests for each MAG
results <- combined_data_magp %>%
  group_by(Bin) %>%
  summarise(
    # Average per condition
    mean_control = mean(Abundance[Tto == "control"], na.rm = TRUE),
    mean_stress = mean(Abundance[Tto == "stress"], na.rm = TRUE),
    # fold change (en escala log2)
    log2Fold_change = mean_stress - mean_control, ###Change to log2(mean_stress / mean_control) in RAW
    # t-test
    p_value = t.test(Abundance[Tto == "stress"], Abundance[Tto == "control"])$p.value
  ) %>%
  ungroup() %>%
  mutate(
    # Calcular -log10(p-value) para el volcano plot
    SCORE_neg_log10_p = -log10(p_value),
    Significance = case_when(
      p_value < 0.05 & abs(log2Fold_change) > 1 ~ "Significant (p<0.05, |FC|>1)",
      p_value < 0.1 & abs(log2Fold_change) >3 ~ "Significant (p<0.1, |FC|>3)",
      TRUE ~ "Not significant"
    )
  )
View(results)
head(results)

# volcano plot
volcano_plot <- ggplot(results, aes(x = log2Fold_change, y = SCORE_neg_log10_p)) +
  geom_point(aes(color = Significance)) +
  scale_color_manual(values = c(
    "Significant (p<0.05, |FC|>1)" = "red",
    "Significant (p<0.1, |FC|>3)" = "orange",
    "Not significant" = "gray"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.1), linetype = "dotted", color = "darkgreen") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-3, 3), linetype = "dotted", color = "darkgreen") +
  labs(x = "Log2 Fold Change (Stress/Control)",
       y = "-Log10(p-value)",
       title = "Volcano Plot MAG abundance calculated from Proteins (log): Control vs Stress",
       color = "Significance") +
  theme_minimal() +
  # Extreme significant points
  geom_text_repel(
    data = subset(results, Significance != "Not significant"),
    aes(label = Bin),
    size = 3,
    box.padding = 0.5,
    max.overlaps = 20
  )

print(volcano_plot)


### Generar una tabla de los MAGs significativos
significance_counts <- results %>%
  count(Significance) %>%
  mutate(Percentage = n / sum(n) * 100)
print(significance_counts)

significant_mags <- results %>%
  filter(Significance != "Not significant") %>%
  arrange(Significance, desc(abs(log2Fold_change)))
View(significant_mags)

write.csv(significant_mags, "significant_MAGPs_clr_volcano.csv", row.names = FALSE)

######Unir la taxonomía para saber de quién se trata
library(readxl)
taxonomy <- read_excel("taxonomia_MetaP2.xlsx") %>%
  rename(Bin = user_genome) %>%
  select(Bin, classification)

View(taxonomy)

significant_mags_with_taxonomy <- significant_mags %>%
  left_join(taxonomy, by = "Bin")

View(significant_mags_with_taxonomy)
write.table(significant_mags_with_taxonomy, "significant_MAGPs_log_volcano_with_taxonomy.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

library(dplyr)

############# MAGP_log Volcano plot with taxonomy

# Transpose the abundance table to have samples as rows and MAGs as columns
abundance_long_magp <- MAGP_abundance_log %>%    
  pivot_longer(cols = -Bin, names_to = "SampleID", values_to = "Abundance")  %>%
  mutate(Abundance = if_else(Abundance == 0, rnorm(n(), mean = 8, sd = 1), as.numeric(Abundance))) #MAGNUS correction in log

# Combine with metadata
combined_data_magp <- abundance_long_magp %>%
  left_join(metadata, by = "SampleID")

# Perform statistical tests for each MAG
results <- combined_data_magp %>%
  group_by(Bin) %>%
  summarise(
    mean_control = mean(Abundance[Tto == "control"], na.rm = TRUE),
    mean_stress = mean(Abundance[Tto == "stress"], na.rm = TRUE),
    # fold change
    log2Fold_change = mean_stress - mean_control,
    # t-test
    p_value = t.test(Abundance[Tto == "stress"], Abundance[Tto == "control"])$p.value
  ) %>%
  ungroup() %>%
  mutate(
    SCORE_neg_log10_p = -log10(p_value),
    Significance = case_when(
      p_value < 0.05 & abs(log2Fold_change) > 1 ~ "Significant (p<0.05, |FC|>1)",
      p_value < 0.1 & abs(log2Fold_change) > 3 ~ "Significant (p<0.1, |FC|>3)",
      TRUE ~ "Not significant"
    )
  )

results_with_tax <- results %>%
  left_join(MAG_taxonomy %>% 
              select(user_genome, classification),
            by = c("Bin" = "user_genome")) %>%  # Cambiar a Bin
  mutate(
    last_taxon = sub(".*;(g__[^;]+).*", "\\1", classification),
    last_taxon = ifelse(grepl(";", last_taxon), 
                        sub(".*;(f__[^;]+).*", "\\1", classification),
                        last_taxon),
    last_taxon = ifelse(grepl(";", last_taxon), 
                        sub(".*;(o__[^;]+).*", "\\1", classification),
                        last_taxon),
    
    Species_label = ifelse(grepl(";s__[^;]+", classification),
                           sub(".*;s__([^;]+)", "\\1", classification),
                           paste0(last_taxon, " (NS)")),
    
    Species_short = ifelse(nchar(Species_label) > 25,
                           paste0(substr(Species_label, 1, 22), "..."),
                           Species_label)
  )

# volcano plot with taxonomy
volcano_plot <- ggplot(results_with_tax, aes(x = log2Fold_change, y = SCORE_neg_log10_p)) +
  geom_point(aes(color = Significance)) +
  scale_color_manual(values = c(
    "Significant (p<0.05, |FC|>1)" = "red",
    "Significant (p<0.1, |FC|>3)" = "orange",
    "Not significant" = "gray"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.1), linetype = "dotted", color = "darkgreen") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-3, 3), linetype = "dotted", color = "darkgreen") +
  labs(x = "Log2 Fold Change (Stress/Control)",
       y = "-Log10(p-value)",
       title = "Volcano Plot MAG abundance calculated from Proteins (log): Control vs Stress",
       color = "Significance") +
  theme_minimal() +
  # Extreme significant points with taxonomy labels
  geom_text_repel(
    data = subset(results_with_tax, Significance != "Not significant"),
    aes(label = Species_short),  # Usar la etiqueta taxonómica corta
    size = 3,
    box.padding = 0.5,
    max.overlaps = 20
  )

print(volcano_plot)

#########################################Volcano plot for RAW/LOG/CLR protein abundance (relative abundance)

metadata <- read.table("metadataMetaP.txt", header = T, sep = "\t") 

protein_abundance <- protein_abundance_log #change to log, to raw, to clr

View(protein_abundance)

#Transpose the abundance table to have samples as rows and proteins as columns
abundance_long_p <- protein_abundance %>%
  rownames_to_column(var = "Protein") %>%
  pivot_longer(cols = -Protein, names_to = "SampleID", values_to = "Abundance") 

combined_data_p <- abundance_long_p %>%
  left_join(metadata, by = "SampleID") %>%
  mutate(Abundance = if_else(Abundance == 0, rnorm(n(), mean = 6, sd = 1), as.numeric(Abundance))) #MAGNUS correction in log

View(combined_data_p)

#Perform statistical tests for each protein
results_p <- combined_data_p %>%
  group_by(Protein) %>%
  summarise(
    # Average per condition
    mean_control = mean(Abundance[Tto == "control"], na.rm = TRUE),
    mean_stress = mean(Abundance[Tto == "stress"], na.rm = TRUE),
    # fold change (en escala log2)
    log2Fold_change = mean_stress - mean_control, ###Change to log2(mean_stress / mean_control) in RAW
    # t-test
    p_value = t.test(Abundance[Tto == "stress"], Abundance[Tto == "control"])$p.value
  ) %>%
  ungroup() %>%
  mutate(
    # Calcular -log10(p-value) para el volcano plot
    SCORE_neg_log10_p = -log10(p_value)
  )
View(results)
head(results)

# volcano plot
volcano_plot <- ggplot(results_p, aes(x = log2Fold_change, y = SCORE_neg_log10_p)) +
  geom_point(aes(color = ifelse(p_value < 0.05 & abs(log2Fold_change) > 1, "Significant", "Not significant"))) +
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(x = "Log2 Fold Change (Stress/Control)",
       y = "-Log10(p-value)",
       title = "Volcano Plot Protein abundance log: Control vs Stress",
       color = "Significance") +
  theme_minimal()

print(volcano_plot)


############# Identification of significant proteins #############

significant_proteins <- results_p %>%
  filter(p_value < 0.05 & abs(log2Fold_change) > 1) %>%
  arrange(p_value, desc(abs(log2Fold_change)))

protein_counts <- results_p %>%
  summarise(
    Total_Proteins = n(),
    Significant = sum(p_value < 0.05 & abs(log2Fold_change) > 1),
    Not_Significant = Total_Proteins - Significant,
    Percentage_Significant = (Significant / Total_Proteins) * 100
  )

View(significant_proteins)

#Save results
write.csv(significant_proteins, "significant_proteins_log.csv", row.names = FALSE)
write.csv(protein_counts, "protein_significance_counts_log.csv", row.names = FALSE)


#Add the values to the volcano plot
volcano_plot <- volcano_plot +
  labs(subtitle = paste0("Significant proteins: ", protein_counts$Significant, 
                         " of ", protein_counts$Total_Proteins, 
                         " (", round(protein_counts$Percentage_Significant, 1), "%)"))

print(volcano_plot)




