library(readxl)
library(dplyr)
library(writexl)

#INPUTS
Functions_table <- read_excel("Functions_table.xlsx")

#KO_IDs that are enriched from the total discriminant proteins
ko_enrichment_control_discriminant <- read.delim("ko_enrichment_discriminant_control_significant.txt", header = TRUE, sep ="\t")
ko_enrichment_stress_discriminant <- read.delim("ko_enrichment_discriminant_stress_significant.txt", header = TRUE, sep ="\t")
#ko_enrichment_control_significant <- read.delim("ko_enrichment_significant_control_significant.txt", header = TRUE, sep ="\t")
ko_enrichment_stress_significant <- read.delim("ko_enrichment_significant_stress_significant.txt", header = TRUE, sep ="\t")

#Summed KO_IDs that are discriminant
koid_control_stress_PLSDA <- read.delim("koid_control_stress_NEW.txt", header = TRUE, sep ="\t")

#Summed KO_IDs that are not only discriminant but also differentially abundant
significant_ko_LOG_volcano<- read.delim("significant_ko_LOG_volcano.txt", header = TRUE, sep ="\t")
View(significant_ko_LOG_volcano)
#Subsets
PLSDA_koid_control <- subset(koid_control_stress_PLSDA, Association == "control")
PLSDA_koid_stress <- subset(koid_control_stress_PLSDA, Association == "stress")

View(ko_enrichment_control_discriminant)
View(ko_enrichment_stress_significant)

#Adding p-values from enrichment analysis

Functions_table <- Functions_table %>%
  left_join(
    ko_enrichment_control_discriminant %>% 
      select(ko_id, p_value_control = p_value),
    by = c("KO_ID" = "ko_id")
  ) %>%
  left_join(
    ko_enrichment_stress_discriminant %>% 
      select(ko_id, p_value_stress = p_value),
    by = c("KO_ID" = "ko_id")
  ) %>%
  mutate(
    Enrichment_pvalue = case_when(
      Condition == "Control" ~ p_value_control,
      Condition == "Stress" ~ p_value_stress,
      Condition == "Mixed*" ~ p_value_control,
      TRUE ~ NA_real_
    )
  ) %>%
  select(-p_value_control, -p_value_stress)

#Adding loading values from discriminant analysis

Functions_table <- Functions_table %>%
  left_join(
    PLSDA_koid_control %>% 
      select(ko_id, Loading_Value_control = Loading_Value),
    by = c("KO_ID" = "ko_id")
  ) %>%
  left_join(
    PLSDA_koid_stress %>% 
      select(ko_id, Loading_Value_stress = Loading_Value),
    by = c("KO_ID" = "ko_id")
  ) %>%
  mutate(
    Discriminant_Loading = case_when(
      Condition == "Control" ~ Loading_Value_control,
      Condition == "Stress" ~ Loading_Value_stress,
      Condition == "Mixed*" ~ Loading_Value_stress,
      TRUE ~ NA_real_
    )
  ) %>%
  select(-Loading_Value_control, -Loading_Value_stress)


# Adding p-values and log2Fold_change from differential analysis

Functions_table <- Functions_table %>%
  left_join(
    significant_ko_LOG_volcano %>%
      select(
        ko_id,
        Differential_abundance_pvalue = p_value,
        Differential_log2Fold_change = log2Fold_change
      ),
    by = c("KO_ID" = "ko_id")
  )


#################Add the taxonomy

discr_prot_tax <- read.delim("Discr_prot_MetaP_PLSDA_tax.txt", header = TRUE, sep ="\t")
View(discr_prot_tax)
View(masterfile)

discr_prot_tax <- discr_prot_tax %>%
  left_join(
    masterfile %>% 
      select(gene, ko_id),
    by = c("Protein" = "gene")
  )

###############################ADD MANNUALLY ARCHAEA CLASSIFICATIONS (2nd sheet taxonomy file)
proteins_discriminant_koids <- discr_prot_tax %>% 
  filter(ko_id %in% Functions_table$KO_ID)

write_xlsx(proteins_discriminant_koids, 
           "proteins_discriminant_koids.xlsx")

discr_prot_tax <- read_excel("proteins_discriminant_koids.xlsx")

View(discr_prot_tax)
##############################


Functions_table <- Functions_table %>%
  left_join(
    discr_prot_tax %>% select(ko_id, classification),
    by = c("KO_ID" = "ko_id"),
    relationship = "many-to-many"
  )


View(proteins_discriminant_koids)
View(Functions_table)

write.table(Functions_table, file = "Functions_table_expanded.txt", sep = "\t")


write_xlsx(Functions_table, 
           "Functions_table_definitiva.xlsx")

