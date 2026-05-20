library(tidyverse)
library(BiocManager)
library(ggrepel)
library(KEGGREST)
library(readxl)
library(writexl)

##########################################################EXTRACT LFQ VALUES (SKIP, this is done)!!!!!!!!!
#  911,040 genes masterfile with ko_id
masterfile <- read.delim("masterfile_genes_proteins_gutbrain.txt", header = TRUE, quote = "")
View(masterfile)
head(masterfile)

#metadata
metadata <- read.table("metadataMetaP.txt", header = T, sep = "\t") 
View(metadata)

#Get sample_columns

sample_columns <- grep("X[0-9]+", colnames(masterfile), value = TRUE)
View(sample_columns)

########################################## Get Function abundances (SKIP)!!!!!!!!!

summedLFQ_functions <- masterfile %>%
  mutate(across(all_of(sample_columns), ~2^.x)) %>% #log2 to raw (2^value)
  filter(Expressed == TRUE, !is.na(ko_id)) %>%
  select(ko_id, all_of(sample_columns)) %>%
  pivot_longer(-ko_id, names_to = "SampleID", values_to = "value") %>%
  group_by(ko_id, SampleID) %>%
  summarize(SummedLFQ = sum(value, na.rm=T))
      
write.table(summedLFQ_functions, "summedLFQ_functions.txt", sep = "\t", row.names = FALSE)

#Add metadata

summedLFQ_functions <- summedLFQ_functions %>%
  left_join(metadata, by = "SampleID")

head(summedLFQ_functions)

########################################## Filter ko_id by frequency 3 (SKIP)!!!!!!!!!!!

summedLFQ_functions_filtfreq <- summedLFQ_functions %>%
  group_by(ko_id) %>%
  filter(
    # 3 samples non-zero in control
    sum(Tto == "control" & SummedLFQ > 0) >= 3 |
    # 3 samples non-zero in stress
    sum(Tto == "stress" & SummedLFQ > 0) >= 3
  ) %>%
  ungroup()

write.table(summedLFQ_functions_filtfreq, "summedLFQ_functions_filtfreq.txt", sep = "\t", row.names = FALSE)

View(summedLFQ_functions)
View(summedLFQ_functions_filtfreq)

################################################################ Apply a PLS-DA (use the other script) and load the results

#import ko_id abundances per sample (log and clr) and discriminant ko_id (log and clr)

functions_filtfreq_clr <- read.table("functions_filtfreq_clr.txt", header = T, sep = "\t") 
functions_filtfreq_log <- read.table("functions_filtfreq_log.txt", header = T, sep = "\t") #1836 SUMMED KO_ID LOG ABUNDANCES

relevant_koid_log_names <- read.table("relevant_koid_log_names.txt", header = T, sep = "\t") #251 DISCRIMINANT SUMMED KO_ID LOG LIST NAMES
relevant_koid_clr_names <- read.table("relevant_koid_clr_names.txt", header = T, sep = "\t") 

View(functions_filtfreq_log)
View(relevant_koid_log_names)

#Filter abundances by discriminant ko_id

d_functions_log <- functions_filtfreq_log %>% #251 entries
  filter(ko_id %in% relevant_koid_log_names$ko_id)

d_functions_clr <- functions_filtfreq_clr %>%
  filter(ko_id %in% relevant_koid_clr_names$ko_id)

View(d_functions_log)

write.table(d_functions_log, "d_functions_log_abundance.txt", sep = "\t", row.names = FALSE)
write.table(d_functions_clr, "d_functions_clr_abundance.txt", sep = "\t", row.names = FALSE)


#Adapt the format for volcano plot and add metadata

summedLFQ_functions_d_log<- d_functions_log %>% #15,060 entries
  pivot_longer(-ko_id, names_to = "SampleID", values_to = "SummedLFQ") %>%
  group_by(ko_id, SampleID) %>%
  left_join(metadata, by = "SampleID")

summedLFQ_functions_d_clr<- d_functions_clr %>%
  pivot_longer(-ko_id, names_to = "SampleID", values_to = "SummedLFQ") %>%
  group_by(ko_id, SampleID)  %>%
  left_join(metadata, by = "SampleID")

View(summedLFQ_functions_d_log)
View(summedLFQ_functions_d_clr)

########################################################## VOLCANO PLOT

######################### with FDR


# Statistical analyses with FDR
results <- summedLFQ_functions_d_log %>%
  group_by(ko_id) %>%
  summarise(
    mean_control = mean(SummedLFQ[Tto == "control"], na.rm = TRUE),
    mean_stress = mean(SummedLFQ[Tto == "stress"], na.rm = TRUE),
    log2Fold_change = mean_stress - mean_control,
    p_value = t.test(SummedLFQ[Tto == "stress"], 
                     SummedLFQ[Tto == "control"])$p.value
  ) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(p_value, method = "fdr"),
    SCORE_neg_log10_p = -log10(p_value),
    Significance = case_when(
      FDR < 0.05 & abs(log2Fold_change) > 1 ~ "FDR < 0.05 & |FC| > 1",
      FDR < 0.1 & abs(log2Fold_change) > 3 ~ "FDR < 0.1 & |FC| > 3",
      p_value < 0.05 & abs(log2Fold_change) > 1 ~ "p < 0.05 & |FC| > 1",
      p_value < 0.1 & abs(log2Fold_change) > 3 ~ "p < 0.1 & |FC| > 3",
      TRUE ~ "Not significant")
  )

View(results) #251 entries

write.table(results, "d_functions_log_statistics_FDR.txt", sep = "\t", row.names = FALSE)
write.table(results, "d_functions_clr_statistics_FDR.txt", sep = "\t", row.names = FALSE)

summary(results$log2Fold_change)
summary(results$SCORE_neg_log10_p)

#VOLCANO PLOT
volcano_plot <- ggplot(results, aes(x = log2Fold_change, y = SCORE_neg_log10_p)) +
  geom_point(aes(color = Significance), alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c(
      "FDR < 0.05 & |FC| > 1" = "blue",
      "FDR < 0.1 & |FC| > 3" = "darkblue",
      "p < 0.05 & |FC| > 1" = "red",
      "p < 0.1 & |FC| > 3" = "orange",
      "Not significant" = "gray"
    ),
    drop = FALSE
  ) +
  # Líneas de referencia
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.1), linetype = "dotted", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-3, 3), linetype = "dotted", color = "black") +
  labs(
    x = "Log2 Fold Change (Stress/Control)",
    y = "-Log10(p-value)",
    title = "Discriminant KO_ID LOG: Control vs Stress",
    color = "Significance Level"
  ) +
  theme_minimal()
  

print(volcano_plot)

# significant ko_id identification

significant_ko <- results %>%
  filter(
    (FDR < 0.05 & abs(log2Fold_change) > 1) |      
      (FDR < 0.1 & abs(log2Fold_change) > 3) |       
      (p_value < 0.05 & abs(log2Fold_change) > 1) |  
      (p_value < 0.1 & abs(log2Fold_change) > 3)     
  ) %>%
  arrange(FDR, desc(abs(log2Fold_change)))         

write.table(significant_ko, file = "significant_ko_LOG_volcano.txt", sep = "\t", row.names = FALSE)

View(significant_ko)

###Strict significant
significant_ko_p005 <- results %>%
  filter(
    (p_value < 0.05 & abs(log2Fold_change) > 1))        

write.table(significant_ko, file = "significant_ko_LOG_volcano_p005.txt", sep = "\t", row.names = FALSE)

View(significant_ko_p005)

# Subtitle for volcano
ko_percentages <- results %>%
  summarise(
    Total = n(),
    p005_FC1_perc = round(sum(p_value < 0.05 & abs(log2Fold_change) > 1, na.rm = TRUE) / Total * 100, 1),
    p01_FC3_perc = round(sum(p_value < 0.1 & abs(log2Fold_change) > 3, na.rm = TRUE) / Total * 100, 1)
  )

plot_subtitle <- sprintf(
  "Discriminant KO_IDs that are significant: p<0.05 & |FC|>1 = %g%%, p<0.1 & |FC|>3 = %g%%",
  ko_percentages$p005_FC1_perc,
  ko_percentages$p01_FC3_perc
)

volcano_plot <- volcano_plot +
  labs(subtitle = plot_subtitle) +
  theme(
    plot.subtitle = element_text(size = 9, color = "black")
  )
    
print(volcano_plot)



###########################Get ko_id annotation

masterfile <- read.delim("masterfile_genes_proteins_gutbrain.txt", header = TRUE, quote = "")

significant_ko <- read.delim("significant_ko_LOG_volcano.txt", header = TRUE, quote = "") #LOG
significant_ko <- read.delim("significant_ko_CLR_volcano.txt", header = TRUE, quote = "") #CLR

significant_ko_annotated <- significant_ko %>%
  left_join(masterfile %>% 
              dplyr::select(ko_id, kegg_hit) %>%
              distinct(ko_id, .keep_all = TRUE),
            by = "ko_id") %>%
  rename(KEGG_annotation = kegg_hit)

write.table(significant_ko_annotated, file = "significant_ko_LOG_volcano_an.txt", sep = "\t", row.names = FALSE) #LOG
write.table(significant_ko_annotated, file = "significant_ko_CLR_volcano_an.txt", sep = "\t", row.names = FALSE) #CLR
View(significant_ko_annotated)




