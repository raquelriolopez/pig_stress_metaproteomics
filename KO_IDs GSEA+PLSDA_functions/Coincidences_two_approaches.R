library(dplyr)
library(writexl)
############KO_ID matches that are enriched and summedlfq discriminants

#Enrichment
ko_enrichment_control_discriminant <- read.delim("ko_enrichment_discriminant_control_significant.txt", header = TRUE, sep ="\t") #30 KO_IDs
ko_enrichment_stress_discriminant <- read.delim("ko_enrichment_discriminant_stress_significant.txt", header = TRUE, sep ="\t") #40 KO_IDs

#PLS-DA
koid_control_stress_PLSDA <- read.delim("koid_control_stress_NEW.txt", header = TRUE, sep ="\t") #251 KO_IDs PLS-DA

#ko_enrichment_control_significant <- read.delim("ko_enrichment_significant_control_significant.txt", header = TRUE, sep ="\t")
#ko_enrichment_stress_significant <- read.delim("ko_enrichment_significant_stress_significant.txt", header = TRUE, sep ="\t")

View(koid_control_stress_PLSDA)
View(ko_enrichment_control_significant)

#Subset Summed LFQ ko_id PLSDA control and stress

koid_control <- subset(koid_control_stress_PLSDA, Association == "control") #93 KO_IDs
koid_stress <- subset(koid_control_stress_PLSDA, Association == "stress") #158 KO_IDs

#Volcano KO_IDs

volcano_significant <- read.table("significant_ko_LOG_volcano_an.txt", header = T, sep = "\t")  #68 KO_IDs 

volcano_significant <- volcano_significant %>%
  mutate(
    condition = case_when(
      log2Fold_change < 0 ~ "control",
      log2Fold_change > 0 ~ "stress"
    )
  )

View(volcano_significant)

#MASTERFILE
masterfile <- read.delim("masterfile_genes_proteins_gutbrain.txt", header = TRUE, quote = "")

##########################################################

# Helper
extract_ko <- function(df, analysis_name, ko_col = "ko_id") {
  df %>%
    select(all_of(ko_col)) %>%
    distinct() %>%
    rename(ko_id = all_of(ko_col)) %>%
    mutate(analysis = analysis_name)
}

# Enrichment
enrich_control <- extract_ko(
  ko_enrichment_control_discriminant,
  "Enriched_from_discr_proteins_control"
)

enrich_stress <- extract_ko(
  ko_enrichment_stress_discriminant,
  "Enriched_from_discr_proteins_stress"
)

# PLS-DA
plsda_control <- extract_ko(
  koid_control,
  "Discriminant_summed_koid_control"
)

plsda_stress <- extract_ko(
  koid_stress,
  "Discriminant_summed_koid_stress"
)

# Differential (volcano)
diff_control <- extract_ko(
  volcano_significant %>% filter(condition == "control"),
  "Differentially_abundant_summed_koid_control"
)

diff_stress <- extract_ko(
  volcano_significant %>% filter(condition == "stress"),
  "Differentially_abundant_summed_koid_stress"
)


ko_long <- bind_rows(
  enrich_control,
  enrich_stress,
  plsda_control,
  plsda_stress,
  diff_control,
  diff_stress
)

ko_summary <- ko_long %>%
  group_by(ko_id) %>%
  summarise(
    counts = n_distinct(analysis),
    approach_result = paste(sort(unique(analysis)), collapse = " + "),
    .groups = "drop"
  )

masterfile_ko_id_hit <- masterfile %>%
  select(ko_id, kegg_hit) %>%
  distinct()

ko_final <- ko_summary %>%
  left_join(masterfile_ko_id_hit, by = "ko_id") %>%
  rename(`KEGG hit` = kegg_hit)

View(ko_final)

write_xlsx(
  ko_final,
  "ko_id_coincidences_twoapproaches_COMPLETE.xlsx"
)

##############WITH YES OR NO for each analysis

ko_presence <- ko_long %>%
  distinct(ko_id, analysis) %>%
  mutate(present = "YES") %>%
  tidyr::pivot_wider(
    names_from  = analysis,
    values_from = present,
    values_fill = "NO"
  )

analysis_cols <- setdiff(colnames(ko_presence), "ko_id")

ko_presence <- ko_presence %>%
  mutate(
    counts = rowSums(across(all_of(analysis_cols), ~ . == "YES"))
  )

ko_presence <- ko_presence %>%
  mutate(
    approach_result = apply(
      select(., all_of(analysis_cols)),
      1,
      function(x) {
        paste(names(x)[x == "YES"], collapse = " + ")
      }
    )
  )

ko_final <- ko_presence %>%
  left_join(masterfile_ko_id_hit, by = "ko_id") %>%
  rename(`KEGG hit` = kegg_hit)


write_xlsx(
  ko_final,
  "ko_id_coincidences_twoapproaches_WITH_FLAGS.xlsx"
)

View(ko_final)



###Conteos

enriched_ko <- ko_final %>%
  filter(
    Enriched_from_discr_proteins_control == "YES" |
      Enriched_from_discr_proteins_stress  == "YES"
  )

n_enriched <- nrow(enriched_ko)
n_enriched


enriched_and_discriminant <- enriched_ko %>%
  filter(
    Discriminant_summed_koid_control == "YES" |
      Discriminant_summed_koid_stress  == "YES"
  )

n_enriched_and_discriminant <- nrow(enriched_and_discriminant)
n_enriched_and_discriminant

View(enriched_and_discriminant)
enriched_discriminant_and_differential <- enriched_and_discriminant %>%
  filter(
    Differentially_abundant_summed_koid_control == "YES" |
      Differentially_abundant_summed_koid_stress  == "YES"
  )

n_enriched_discriminant_and_differential <- 
  nrow(enriched_discriminant_and_differential)

n_enriched_discriminant_and_differential

cat("Enriched KO_IDs (control or stress):", n_enriched, "\n")
cat("Of these, discriminatory:, n_enriched_and_discriminant", "\n")
cat("Of these, also differentially abundant:", 
    n_enriched_discriminant_and_differential, "\n")


enriched_discriminant_and_differential %>%
  select(
    ko_id,
    approach_result,
    `KEGG hit`
  ) %>%
  View()



