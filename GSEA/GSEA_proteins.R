library(dplyr)
library(purrr)
library(tidyverse)
library(readxl)
library(writexl)
library(BiocManager)
BiocManager::install("KEGGREST")
library(KEGGREST)

packageVersion("KEGGREST")

######################INPUT

#  911,040 genes masterfile with ko_id
masterfile <- read.delim("masterfile_genes_proteins_gutbrain.txt", header = TRUE, quote = "")

# 42,696 proteins, rownames 1, 2, 3...  
legend <- read.delim("Leyenda_MetaP2.txt", header = TRUE, quote = "") %>% rename(gene = Protein)

# 26,868 proteins (filtered by frequency), rownames P2, P3,...
protein_abundance_log <- read.table("MetaP2_log_filtfreq.tsv", header=T, sep="\t") ###LOG protein abundance
protein_abundance_log <- as.data.frame(protein_abundance_log)

# DIFFERENTIALLY ABUNDANT 1,253 proteins first column named "Proteins" whose entries are "P21738"...
#significant_proteins <- read.csv("significant_proteins_log.csv", header = TRUE, sep = ",")

# 1001 DISCRIMINANT proteins
discriminant_proteins <- read.delim("discriminant_proteins_control_stress_plsda_NEW.txt", header = TRUE, quote = "")  %>% 
  rename(gene = Protein) %>%
  rename(Protein = code_protein) %>%
  mutate(Protein = gsub("^X", "P", Protein))

View(discriminant_proteins)
View(protein_abundance_log)
View(legend)
############################FORMAT

#Adapt the format and obtain the list of filtered proteins by frequency and the list of significant proteins with associated gene.

protein_abundance_log <- rownames_to_column(protein_abundance_log, var = "Protein")

legend <- rownames_to_column(legend, var = "Protein")
legend$Protein <- paste0("P", legend$Protein)

filt_prot <- legend %>%
  filter(Protein %in% protein_abundance_log$Protein)

#sign_prot <- legend %>%
  #filter(Protein %in% significant_proteins$Protein)

discr_prot <- legend %>%  #1001 discriminant proteins
  filter(Protein %in% discriminant_proteins$Protein)

View(filt_prot)
View(discr_prot)

###########################Get ko_id

masterfile_ko_id <- masterfile %>% 
  select(gene, ko_id)

filt_prot <- filt_prot %>% 
  left_join(masterfile_ko_id, by = "gene")

#sign_prot <- sign_prot %>% 
  #left_join(masterfile_ko_id, by = "gene")

discr_prot <- discr_prot %>% 
  left_join(masterfile_ko_id, by = "gene")

# Save tables
write.table(filt_prot, "filt_prot.txt", sep = "\t" , quote = FALSE)
#write.table(sign_prot, "sign_prot.txt", sep = "\t", quote = FALSE)
write.table(discr_prot, "discr_prot.txt", sep = "\t", quote = FALSE)

###########
#Filter proteins without ko_id

filt_prot_total_with_ko_id <- filt_prot %>%      #  From 26,868 proteins to 19,179 entries
  filter(!is.na(ko_id))

#sign_prot <- sign_prot %>%      # From 1,244 proteins to 852 entries
  #filter(!is.na(ko_id))

discr_prot_total_with_ko_id <- discr_prot %>%      # From 1,001 proteins to 678 entries
  filter(!is.na(ko_id))

View(filt_prot_total_with_ko_id)
View(discr_prot_total_with_ko_id)


############################################################
### KO ENRICHMENT — CONTROL vs STRESS ###
############################################################


### 1. Background counts (proteínas filtradas + con KO_ID)
background_counts <- filt_prot_total_with_ko_id %>% # # 1,783 unique KO_IDs
  count(ko_id, name = "background_count")

total_background <- nrow(filt_prot_total_with_ko_id)   # 19,179 proteins with KO_ID in background

View(background_counts)

############################################################
### CONTROL
############################################################

# Proteínas control con KO_ID
control_prots_with_ko <- discr_control %>%    #227 discriminant proteins for control condition with KO_ID
  left_join(dplyr::select(filt_prot, Protein, ko_id), by = "Protein") %>%
  filter(!is.na(ko_id))

View(control_prots_with_ko) 

# Tamaño correcto de la muestra
total_control <- nrow(control_prots_with_ko)  # 227 discriminant proteins for control condition with KO_ID

# Conteo por KO
control_counts <- control_prots_with_ko %>% # 150 unique KO_IDs from the previous 227 proteins
  count(ko_id, name = "control_count")

View(control_counts)

# Merge background vs control
ko_enrichment_control <- full_join(background_counts,
                                   control_counts,
                                   by = "ko_id") %>%
  mutate(
    background_count = tidyr::replace_na(background_count, 0),
    control_count    = tidyr::replace_na(control_count, 0)
  )

# Test hipergeométrico
ko_enrichment_control <- ko_enrichment_control %>%
  rowwise() %>%
  mutate(
    p_value = ifelse(control_count > 0,
                     phyper(
                       q = control_count - 1,
                       m = background_count,
                       n = total_background - background_count,
                       k = total_control,
                       lower.tail = FALSE
                     ),
                     NA_real_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(p_value))

# Adjust
ko_enrichment_control$fdr <- p.adjust(ko_enrichment_control$p_value, method = "BH")

# Significativos
discriminant_kos_control <- ko_enrichment_control %>% #30 enriched KO_IDs in control
  filter(p_value < 0.05) %>%
  arrange(p_value)

View(discriminant_kos_control)

############################################################
### STRESS
############################################################

# Stress proteins with ko_id
stress_prots_with_ko <- discr_stress %>% #451 discriminant proteins with ko_id
  left_join(dplyr::select(filt_prot, Protein, ko_id), by = "Protein") %>%
  filter(!is.na(ko_id))

View(stress_prots_with_ko)

# Tamaño correcto de la muestra
total_stress <- nrow(stress_prots_with_ko)  # 451

# Conteo por KO
stress_counts <- stress_prots_with_ko %>%
  count(ko_id, name = "stress_count")

View(stress_counts) #286

# Merge background vs stress
ko_enrichment_stress <- full_join(background_counts,
                                  stress_counts,
                                  by = "ko_id") %>%
  mutate(
    background_count = tidyr::replace_na(background_count, 0),
    stress_count     = tidyr::replace_na(stress_count, 0)
  )

# Test hipergeométrico
ko_enrichment_stress <- ko_enrichment_stress %>%
  rowwise() %>%
  mutate(
    p_value = ifelse(stress_count > 0,
                     phyper(
                       q = stress_count - 1,
                       m = background_count,
                       n = total_background - background_count,
                       k = total_stress,
                       lower.tail = FALSE
                     ),
                     NA_real_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(p_value))

# Adjust
ko_enrichment_stress$fdr <- p.adjust(ko_enrichment_stress$p_value, method = "BH")

# Significativos
discriminant_kos_stress <- ko_enrichment_stress %>% #40 enriched KO_IDs in stress
  filter(p_value < 0.05) %>%
  arrange(p_value)

View(discriminant_kos_stress)
############################################################
### Añadir anotaciones KEGG (igual que antes)
############################################################

masterfile_ko_id_hit <- masterfile %>% 
  dplyr::select(ko_id, kegg_hit) %>%
  distinct(ko_id, .keep_all = TRUE)

add_kegg_from_masterfile <- function(df) {
  df %>% 
    left_join(masterfile_ko_id_hit, by = "ko_id") %>% 
    rename(`KEGG hit` = kegg_hit)
}

discriminant_kos_control <- add_kegg_from_masterfile(discriminant_kos_control)
disc_ko_enrichment_control <- add_kegg_from_masterfile(ko_enrichment_control)
discriminant_kos_stress <- add_kegg_from_masterfile(discriminant_kos_stress) 
disc_ko_enrichment_stress <- add_kegg_from_masterfile(ko_enrichment_stress)

############################################################
### Exportar (igual que antes)
############################################################

write.table(disc_ko_enrichment_control, "disc_ko_enrichment_control_all.txt", sep = "\t", row.names = FALSE)
write.table(discriminant_kos_control, "ko_enrichment_control_discriminant.txt", sep = "\t", row.names = FALSE)
write.table(disc_ko_enrichment_stress, "disc_ko_enrichment_stress_all.txt", sep = "\t", row.names = FALSE)
write.table(discriminant_kos_stress, "ko_enrichment_stress_discriminant.txt", sep = "\t", row.names = FALSE)












###################################################################################KO_ID enrichment STRESS vs CONTROL

###############Subset significant_proteins into control or stress

# Control
sig_control <- significant_proteins %>%
  filter(log2Fold_change < 0)

# Stress
sig_stress <- significant_proteins %>%
  filter(log2Fold_change > 0)

nrow(sig_control) #467 proteins associated to control 
nrow(sig_stress)  #786 proteins associated to stress

###############Subset discriminant_proteins into control or stress

# Control
discr_control <- discriminant_proteins %>%
  filter(Association == "control")

# Stress
discr_stress <- discriminant_proteins %>%
  filter(Association == "stress")

View(discr_control)
View(discr_stress)
View(filt_prot)

nrow(discr_control) #345/1001 proteins that distinguish control
nrow(discr_stress) #656/1001 proteins that distinguish stress

########### ko_id enrichment for control

background_counts <- filt_prot %>%
  count(ko_id, name = "background_count")

#control_counts <- sig_control %>%
#  left_join(select(filt_prot, Protein, ko_id), by = "Protein") %>%
#  filter(!is.na(ko_id)) %>%  
#  count(ko_id, name = "control_count")

control_counts <- discr_control %>%
  left_join(dplyr::select(filt_prot, Protein, ko_id), by = "Protein") %>%
  filter(!is.na(ko_id)) %>%  
  count(ko_id, name = "control_count")


#Replace NA with 0
ko_enrichment_control <- full_join(
  background_counts,
  control_counts,
  by = "ko_id"
) %>%
  mutate(
    background_count = replace_na(background_count, 0),
    control_count = replace_na(control_count, 0)
  )

#Parameters
total_background <- nrow(filt_prot)       # Background proteins with ko_id: 19,179 
#total_control <- nrow(sig_control)        # Significant Control proteins with ko_id: 467
total_control <- nrow(control_prots_with_ko)  # Discriminant Control proteins with ko_id: 451

#Apply hypergeometric test to each ko_id

ko_enrichment_control <- ko_enrichment_control %>%
  rowwise() %>%
  mutate(
    p_value = ifelse(control_count > 0,  
                     phyper(
                       q = control_count - 1,          
                       m = background_count,           
                       n = total_background - background_count, 
                       k = total_control,              
                       lower.tail = FALSE
                     ),
                     NA_real_
    )
  ) %>%
  ungroup()  %>%
  filter(!is.na(p_value))

#Adjust P-values
ko_enrichment_control$fdr <- p.adjust(ko_enrichment_control$p_value, method = "BH")

#Significant ko_id control
#significant_kos_control <- ko_enrichment_control %>%
 # filter(p_value < 0.05) %>%
  #arrange(p_value)

#Discriminant ko_id control
discriminant_kos_control <- ko_enrichment_control %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)

View(discriminant_kos_control) #14 discriminant proteins for control enriched

########### ko_id enrichment for stress

# Count KO occurrences in background
background_counts <- filt_prot %>%
  count(ko_id, name = "background_count")

# Count KO occurrences in stress-enriched proteins
stress_counts <- discr_stress %>%  
  left_join(dplyr::select(filt_prot, Protein, ko_id), by = "Protein") %>%
  filter(!is.na(ko_id)) %>%      # Remove proteins without KO annotation
  count(ko_id, name = "stress_count")

View(stress_counts) #286

# Merge and replace NA with 0
ko_enrichment_stress <- full_join(
  background_counts,
  stress_counts,
  by = "ko_id"
) %>%
  mutate(
    background_count = replace_na(background_count, 0),
    stress_count = replace_na(stress_count, 0)
  )

# Fixed parameters
total_background <- nrow(filt_prot)       # Background proteins with KO: 19,179 
#total_stress <- nrow(sig_stress)          # Significant Stress-enriched proteins: 786 
total_stress <-  nrow(stress_prots_with_ko)          # Discriminant Stress-enriched proteins: 286

View(total_stress)
View(discr_stress)

# Hypergeometric test for each KO
ko_enrichment_stress <- ko_enrichment_stress %>%
  rowwise() %>%
  mutate(
    p_value = ifelse(stress_count > 0,    # Only test KOs present in stress
                     phyper(
                       q = stress_count - 1,            # Successes in sample (P(X ≥ q))
                       m = background_count,             # Successes in background
                       n = total_background - background_count, # Failures in background
                       k = total_stress,                 # Sample size (stress proteins)
                       lower.tail = FALSE
                     ),
                     NA_real_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(p_value))                 

# Multiple testing correction
ko_enrichment_stress$fdr <- p.adjust(ko_enrichment_stress$p_value, method = "BH")

# Significant stress-enriched KOs (p < 0.05)
#significant_kos_stress <- ko_enrichment_stress %>%
 # filter(p_value < 0.05) %>%
#  arrange(p_value)

# Discriminant stress-enriched KOs (p < 0.05)
discriminant_kos_stress <- ko_enrichment_stress %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)

View(discriminant_kos_stress)
########Get KEGG functional annotations for these ko_id
masterfile_ko_id_hit <- masterfile %>% 
  dplyr::select(ko_id, kegg_hit) %>%
  distinct(ko_id, .keep_all = TRUE)

add_kegg_from_masterfile <- function(df) {
  df %>% 
    left_join(masterfile_ko_id_hit, by = "ko_id") %>% 
    rename(`KEGG hit` = kegg_hit)  
}

#
#significant_kos_control <- add_kegg_from_masterfile(significant_kos_control)
#ko_enrichment_control <- add_kegg_from_masterfile(ko_enrichment_control)
#significant_kos_stress <- add_kegg_from_masterfile(significant_kos_stress) 
#ko_enrichment_stress <- add_kegg_from_masterfile(ko_enrichment_stress)

discriminant_kos_control <- add_kegg_from_masterfile(discriminant_kos_control)
disc_ko_enrichment_control <- add_kegg_from_masterfile(ko_enrichment_control)
discriminant_kos_stress <- add_kegg_from_masterfile(discriminant_kos_stress) 
disc_ko_enrichment_stress <- add_kegg_from_masterfile(ko_enrichment_stress)

# Save results

write.table(disc_ko_enrichment_control, "disc_ko_enrichment_control_all.txt", sep = "\t", row.names = FALSE)
write.table(discriminant_kos_control, "ko_enrichment_control_discriminant.txt", sep = "\t", row.names = FALSE)
write.table(disc_ko_enrichment_stress, "disc_ko_enrichment_stress_all.txt", sep = "\t", row.names = FALSE)
write.table(discriminant_kos_stress, "ko_enrichment_stress_discriminant.txt", sep = "\t", row.names = FALSE)

# View results

View(discriminant_kos_control)
View(disc_ko_enrichment_control)
View(discriminant_kos_stress)
View(disc_ko_enrichment_stress)

###Reload

ko_enrichment_stress <- read.table("ko_enrichment_stress_all.txt", header=T, sep="\t")
significant_kos_stress <- read.table("ko_enrichment_stress_significant.txt", header=T, sep="\t")
ko_enrichment_control <- read.table("ko_enrichment_control_all.txt", header=T, sep="\t")
significant_kos_control <- read.table("ko_enrichment_control_significant.txt", header=T, sep="\t")


write.table(discriminant_kos_control, "ko_enrichment_control_discriminant.txt", sep = "\t", row.names = FALSE)
write.table(discriminant_kos_stress, "ko_enrichment_stress_discriminant.txt", sep = "\t", row.names = FALSE)


#####Guardar excel para supplementary material

install.packages("openxlsx")

library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "Control")
addWorksheet(wb, "Stress")

writeData(wb, "Control", discriminant_kos_control)
writeData(wb, "Stress", discriminant_kos_stress)

saveWorkbook(wb, "ko_enrichment_discriminant.xlsx", overwrite = TRUE)






#######Other ideas


####################################################PATHWAY ENRICHMENT

## 1.  Get ko_id for each condition
#Get KOs for control-enriched proteins (log2FC < 0)
control_kos <- sig_control %>%
  left_join(select(filt_prot, Protein, ko_id), by = "Protein") %>%
  filter(!is.na(ko_id)) %>%
  pull(ko_id) %>%
  unique()

# Get KOs for stress-enriched proteins (log2FC > 0)
stress_kos <- sig_stress %>%
  left_join(select(filt_prot, Protein, ko_id), by = "Protein") %>%
  filter(!is.na(ko_id)) %>%
  pull(ko_id) %>%
  unique()

## 2. Get pathway information (same for both conditions)
background_kos <- unique(filt_prot$ko_id)

# Modified function to cache pathway data
get_kegg_pathways_cached <- function(ko_ids) {
  ko_pathways <- list()
  for (ko in ko_ids) {
    # Skip if already processed
    if (!is.null(ko_pathways[[ko]])) next
    
    pathways <- tryCatch({
      keggGet(paste0("ko:", ko))[[1]]$PATHWAY
    }, error = function(e) NULL)
    
    if (!is.null(pathways)) {
      ko_pathways[[ko]] <- names(pathways)
    }
    Sys.sleep(0.1) # Be nice to KEGG server
  }
  return(ko_pathways)
}

# Get pathways for background KOs
background_pathways <- get_kegg_pathways_cached(background_kos)
all_pathways <- unique(unlist(background_pathways))

## 3. Hypergeometric test function (modified for condition-specific analysis)

run_pathway_enrichment <- function(significant_kos, condition_name) {
  results <- lapply(all_pathways, function(pathway) {
    pathway_kos <- names(background_pathways)[sapply(background_pathways, function(x) pathway %in% x)]
    
    q <- sum(significant_kos %in% pathway_kos)
    if (q == 0) return(NULL) #exclude pathways without significant proteins
    m <- sum(background_kos %in% pathway_kos)
    n <- length(background_kos) - m
    k <- length(significant_kos)
    
    data.frame(
      pathway = pathway,
      p_value = phyper(q - 1, m, n, k, lower.tail = FALSE),
      count_in_set = q,
      count_in_background = m,
      condition = condition_name,
      stringsAsFactors = FALSE
    )
  }) %>% 
    bind_rows() %>%
    arrange(p_value) %>% 
    mutate(fdr = p.adjust(p_value, method = "BH"))  # Añadir FDR
}

## 4. Run analysis for both conditions
control_pathways <- run_pathway_enrichment(control_kos, "control")
stress_pathways <- run_pathway_enrichment(stress_kos, "stress")

write.table(control_pathways, "kegg_pathway_enrichment_control.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = "")       
write.table(stress_pathways, "kegg_pathway_enrichment_stress.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = "")

# 5. Get names

control_pathways$pathway_name <- sapply(control_pathways$pathway, function(p) {
  keggGet(p)[[1]]$NAME
})

stress_pathways$pathway_name <- sapply(stress_pathways$pathway, function(p) {
  keggGet(p)[[1]]$NAME
})

#Save objects

#All names for each pathway code
saveRDS(stress_pathways, file = "stress_pathways.RDS")
saveRDS(control_pathways, file = "control_pathways.RDS")

clean_dataframe <- function(df) { #Avoid listing pathway name entries
  df %>%
    mutate(across(where(is.list), ~ sapply(.x, function(x) paste(unlist(x), collapse = "|"))))
}

control_pathways_clean <- clean_dataframe(control_pathways)
stress_pathways_clean <- clean_dataframe(stress_pathways)

write.table(control_pathways_clean, "kegg_pathway_enrichment_control_simplified.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = "")       
write.table(stress_pathways_clean, "kegg_pathway_enrichment_stress_simplified.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = "")

# 7. Combined view
combined_pathways <- bind_rows(control_pathways, stress_pathways) %>%
  select(pathway, pathway_name, everything())

combined_pathways_clean <- clean_dataframe(combined_pathways)

write.table(combined_pathways_clean, "kegg_pathway_enrichment.txt", sep="\t", row.names=FALSE)

## View results
View(control_pathways)
View(stress_pathways)
View(combined_pathways_clean)

head(combined_pathways_clean)
head(control_pathways)





######################
# 1) Proteínas de stress originales (antes de join)
cat("nrow(discr_stress) — proteins labelled 'stress' in discriminant_proteins:\n")
print(nrow(discr_stress))

# 2) Proteínas de stress con ko_id tras el join
stress_prots_with_ko <- discr_stress %>%
  left_join(dplyr::select(filt_prot, Protein, ko_id), by = "Protein") %>%
  filter(!is.na(ko_id))

cat("nrow(stress_prots_with_ko) — proteins 'stress' that HAVE ko_id:\n")
print(nrow(stress_prots_with_ko))

cat("n_distinct(stress_prots_with_ko$ko_id) — distinct KO_IDs among those proteins:\n")
print(n_distinct(stress_prots_with_ko$ko_id))

# 3) Comparar con control
control_prots_with_ko <- discr_control %>%
  left_join(dplyr::select(filt_prot, Protein, ko_id), by = "Protein") %>%
  filter(!is.na(ko_id))

cat("control: nrow proteins with ko_id:", nrow(control_prots_with_ko), "\n")
cat("control: n_distinct ko_id:", n_distinct(control_prots_with_ko$ko_id), "\n")









############################################################### (SKIP IT) KO_ID enrichment (without splitting groups) -Not useful at all

# Count occurrences of each KO in the background and in the significant proteins
background_counts <- filt_prot %>%
  count(ko_id, name = "background_count")

significant_counts <- sign_prot %>%
  count(ko_id, name = "significant_count")

# Merge the data and replace NA with 0
ko_enrichment <- full_join(
  background_counts,
  significant_counts,
  by = "ko_id"
) %>%
  mutate(
    background_count = replace_na(background_count, 0),
    significant_count = replace_na(significant_count, 0)
  )

# Fixed parameters
total_background <- nrow(filt_prot)       # 19,179
total_significant <- nrow(sign_prot)      # 852

# Function to apply the hypergeometric test to each KO
ko_enrichment <- ko_enrichment %>%
  rowwise() %>%
  mutate(
    p_value = phyper(
      q = significant_count - 1,         # Number of successes in the sample (corrected for P(X ≥ q))
      m = background_count,              # Successes in the background
      n = total_background - background_count, # Failures in the background
      k = total_significant,             # Sample size
      lower.tail = FALSE                 # Enrichment test (right tail)
    )
  ) %>%
  ungroup()

# Adjust p-values by FDR
ko_enrichment$fdr <- p.adjust(ko_enrichment$p_value, method = "BH")

# Filter significant KOs (p < 0.05 or FDR < 0.05)
significant_kos <- ko_enrichment %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)

# Save results
write.table(ko_enrichment, "gsea_kos.txt", sep = "\t", row.names = FALSE)
write.table(significant_kos, "gsea_significant_kos.txt", sep = "\t", row.names = FALSE)

View(significant_kos)

View(ko_enrichment)
