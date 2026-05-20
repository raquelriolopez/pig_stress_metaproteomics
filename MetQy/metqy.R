########METQY

# To calculate richness from proteomics data we need to count the number of quantifiable species per sample.
# You have the abundances already, you just need to transfer those to a presence/absence matrix, 1 if it was quantified, 0 if not. 
# Then from this we can calculate the richness measures. 
# Alpha diversity (Shannon index) and beta diversity (Bray-Curtis) can be calculated in R, either manually or using the package vegan. 
# I can share scripts for the manual approach but maybe the vegan package is simpler.

#Libraries
install("BiocManager")
a
library(BiocManager)
BiocManager::install("MetQy")
library(MetQy)
library(tidyverse)
install.packages("pheatmap")
library(pheatmap)
packageVersion("vegan")

#Load masterfile
masterfile <- readRDS(file="masterfile.rds")
View(masterfile)

# For MCF we need to convert the masterfile into a list of K-numbers per MAG like this:

bin_vs_KO <- masterfile %>%
  select(Bin, ko_id) %>%
  filter(!is.na(ko_id)) %>%
  separate_rows(ko_id, sep = ";") %>% #For some unknown reason, DRAM only return one K-number per gene...
  group_by(Bin) %>%
  summarise(K = paste0(ko_id, collapse = ";"))

View(bin_vs_KO)

# The center line is only if your data has multiple K-numbers per gene, but I don't think DRAM has that.

# Then, once we have this matrix, which essentially is just two columns, 
# one with a MAG name and one with a semicolon-separated list of K-numbers, we do the MetQy analysis:

all_metqy <- query_genomes_to_modules(as.data.frame(bin_vs_KO), GENOME_ID_COL = "Bin", GENES_COL = "K", META_OUT = T, ADD_OUT = T)

saveRDS(all_metqy, file="R_all_metqy.rds")

View(all_metqy)
# This takes about 15 minutes so take a coffee break there. 
# Then save the object so you don't have to do it again later. Next time you just rread in the object instead of recalculating it:

all_metqy <- readRDS(file="R_all_metqy.rds")

# Then for plotting I normally add some pretty names for the modules, and plot it with the pheatmap package:

mcf <- as.data.frame(all_metqy$MATRIX) %>%
  select(where(~ any(.x > 0.0))) #Filter the modules to at least >0 mcf in at least one microbe - This can be changed to some other value if you want.

# Create pretty names for all the M0000x ids

KO_short_names <- as.data.frame(cbind(all_metqy$METADATA[,1], all_metqy$METADATA[,3]))

colnames(KO_short_names) <- c("ID", "Short name")

KO_short_names$Name <- apply(KO_short_names[, 1:2], 1, paste, collapse = " - ")

# Add pretty names to the matrix columns

mcf_names <- as.data.frame(colnames(mcf)) %>%
  left_join(KO_short_names, by = c("colnames(mcf)" = "ID")) %>%
  select(Name)

colnames(mcf) <- mcf_names$Name

# Generate heatmap - transposed version

pheatmap(t(mcf), cellwidth = 10, cellheight = 8, fontsize = 6,  filename = "MetQy_mcf.pdf")

# Export txt file

write_delim(mcf %>% rownames_to_column(var="Bin"), file="mcf_metaG.txt", delim = "\t")


#Archaea
#945_SemiBin_201
#MEGAHIT-MetaBAT2-S939.26
#MEGAHIT-MetaBAT2-S946.39
#MEGAHIT-MetaBAT2-S968.197
#group6_bins.257

############SUBSET 3 Characterised MAGs

# Subset the mcf matrix to include only the 3 specific MAGs
target_mags <- c("935_SemiBin_192", "group6_bins.106", "MEGAHIT-MetaBAT2-S972.52")

# Filter the mcf matrix for the target MAGs
mcf_subset <- mcf[rownames(mcf) %in% target_mags, ]

# Generate the smaller heatmap - transposed version
pheatmap(t(mcf_subset), 
         cellwidth = 15, 
         cellheight = 8, 
         fontsize = 8,
         filename = "MetQy_mcf_3MAGs.pdf")

# If you want to keep the same order as your target list:
mcf_subset_ordered <- mcf[target_mags, ]
pheatmap(t(mcf_subset_ordered), 
         cellwidth = 15, 
         cellheight = 8, 
         fontsize = 8,
         filename = "MetQy_mcf_3MAGs_ordered.pdf")


##########SUBSET all the 95 MAGs that are discriminant by MAGG or MAGP

final_table <- read_tsv("all_RELEVANT_MAGs_unified_taxonomy_NEW.tsv")
View(final_table)

final_table <- read_tsv("all_RELEVANT_MAGs_unified_taxonomy_NEW.tsv") %>%
  mutate(
    taxon = case_when(
      str_detect(classification, "s__[^;]+") ~ str_extract(classification, "s__[^;]+"),
      str_detect(classification, "g__[^;]+") ~ str_extract(classification, "g__[^;]+"),
      str_detect(classification, "f__[^;]+") ~ str_extract(classification, "f__[^;]+"),
      TRUE ~ "Unclassified"
    ),
    short_label = paste0(taxon, " (", Genome, ")")
  )

mcf_filtered <- mcf[final_table$Genome, ]
rownames(mcf_filtered) <- final_table$short_label
mcf_filtered[is.na(mcf_filtered)] <- 0

pheatmap(t(mcf_filtered),
         cellwidth = 10,
         cellheight = 8,
         fontsize_row = 6,
         fontsize_col = 7,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         filename = "MetQy_mcf_classification_labels_LOG_NEW.pdf")


#############

##########SUBSET all the 95 MAGs that are discriminant by MAGG or MAGP

final_table <- read_tsv("all_RELEVANT_MAGs_unified_taxonomy_LOG_NEW.tsv") %>%
  mutate(
    taxon = case_when(
      str_detect(classification, "s__[^;]+") ~ str_extract(classification, "s__[^;]+"),
      str_detect(classification, "g__[^;]+") ~ str_extract(classification, "g__[^;]+"),
      str_detect(classification, "f__[^;]+") ~ str_extract(classification, "f__[^;]+"),
      TRUE ~ "Unclassified"
    ),
    short_label = paste0(taxon, " (", Genome, ")")
  )

# Sort by indicator type AND WITHIN EACH GROUP BY TAXONOMY (alphabetical)
final_table_ordered <- final_table %>%
  arrange(taxon, Consistent_Indicator)

mcf_ordered <- mcf[final_table_ordered$Genome, ]
rownames(mcf_ordered) <- final_table_ordered$short_label
mcf_ordered[is.na(mcf_ordered)] <- 0

# FILTER MODULES: keep only those that have at least one value > 0 in some MAG
mcf_filtered_modules <- mcf_ordered[, colSums(mcf_ordered > 0) > 0]

# You can also be stricter and filter modules with very little variability
# For example, keep only modules present in at least X MAGs:
min_mags_with_module <- 2  # adapt value, in this case (present in eat least 2 MAGs)
mcf_filtered_modules_strict <- mcf_ordered[, colSums(mcf_ordered > 0) >= min_mags_with_module]

# Create annotations for the COLUMNS (MAGs)
annotation_col <- data.frame(
  Indicator = final_table_ordered$Consistent_Indicator,
  row.names = final_table_ordered$short_label
)

# Define colors for each type of indicator
annotation_colors <- list(
  Indicator = c(
    "control" = "blue",
    "stress" = "red", 
    "mixed" = "purple"
  )
)

# Calculate gaps to separate the groups
indicator_groups <- rle(final_table_ordered$Consistent_Indicator)
gaps_col <- cumsum(indicator_groups$lengths)[-length(indicator_groups$lengths)]

# Create a vector of colors for the column labels
col_colors <- ifelse(final_table_ordered$Consistent_Indicator == "control", "blue",
                     ifelse(final_table_ordered$Consistent_Indicator == "stress", "red", "purple"))

# Heatmap with filtered modules (only those with at least one value > 0)
pheatmap(t(mcf_filtered_modules),
         cellwidth = 10,
         cellheight = 8,
         fontsize_row = 6,
         fontsize_col = 7,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         gaps_col = gaps_col,
         fontsize_col_color = 7,
         color_labels_col = col_colors,
         filename = "MetQy_mcf_classification_labels_LOG_filtered.pdf")

# Heatmap con filtro más estricto (módulos presentes en al menos 2 MAGs)
pheatmap(t(mcf_filtered_modules_strict),
         cellwidth = 10,
         cellheight = 8,
         fontsize_row = 6,
         fontsize_col = 7,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         gaps_col = gaps_col,
         fontsize_col_color = 7,
         color_labels_col = col_colors,
         filename = "MetQy_mcf_classification_labels_LOG_strict_filtered.png")

########## Add MAG origin (OPTIONAL, there is no need for this)

final_table <- read_tsv("all_RELEVANT_MAGs_unified_taxonomy_LOG.tsv") %>%
  mutate(
    taxon = case_when(
      str_detect(classification, "s__[^;]+") ~ str_extract(classification, "s__[^;]+"),
      str_detect(classification, "g__[^;]+") ~ str_extract(classification, "g__[^;]+"),
      str_detect(classification, "f__[^;]+") ~ str_extract(classification, "f__[^;]+"),
      TRUE ~ "Unclassified"
    ),
    short_label = paste0(taxon, " (", Genome, ")")
  )

# Sort FIRST by Indicator, then by taxonomy, and finally by origin 
# For gaps by Indicator to work correctly
final_table_ordered <- final_table %>%
  arrange(Consistent_Indicator, Origins, taxon)  # ¡Indicator PRIMERO!

mcf_ordered <- mcf[final_table_ordered$Genome, ]
rownames(mcf_ordered) <- final_table_ordered$short_label
mcf_ordered[is.na(mcf_ordered)] <- 0

# FILTER MODULES: keep only those present in at least 2 MAGs
min_mags_with_module <- 2
mcf_filtered_modules_strict <- mcf_ordered[, colSums(mcf_ordered > 0) >= min_mags_with_module]

# Create annotations for the COLUMNS (MAGs) - TWO INDICATORS
annotation_col <- data.frame(
  Indicator = final_table_ordered$Consistent_Indicator,
  Origin = final_table_ordered$Origins,
  row.names = final_table_ordered$short_label
)

# Define colors for both indicators 
annotation_colors <- list(
  Indicator = c(
    "control" = "blue",
    "stress" = "red", 
    "mixed" = "purple"
  ),
  Origin = c(
    "MAGG_PLSDA_log" = "lightgreen",
    "MAGG_volcano" = "darkgreen",
    "MAGG_PLSDA_log; MAGG_volcano" = "turquoise",
    "MAGP_PLSDA_log" = "orange",
    "MAGP_volcano" = "darkorange4",
    "MAGP_PLSDA_log; MAGP_volcano" = "yellow",
    "MAGG_PLSDA_log; MAGP_PLSDA_log" = "pink",
    "MAGG_PLSDA_log; MAGP_volcano" = "coral",
    "MAGG_PLSDA_log; MAGP_PLSDA_log; MAGP_volcano" = "magenta"
  )
)

# Calculate gaps to separate the groups by Indicator
indicator_groups <- rle(final_table_ordered$Consistent_Indicator)
gaps_col <- cumsum(indicator_groups$lengths)[-length(indicator_groups$lengths)]

# Create color vector for column labels (based on Consistent_Indicator)
col_colors <- ifelse(final_table_ordered$Consistent_Indicator == "control", "blue",
                     ifelse(final_table_ordered$Consistent_Indicator == "stress", "red", "purple"))

# Heatmap with gaps by Indicator
pheatmap(t(mcf_filtered_modules_strict),
         cellwidth = 10,
         cellheight = 8,
         fontsize_row = 6,
         fontsize_col = 7,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         annotation_names_col = TRUE,
         annotation_legend = TRUE,
         gaps_col = gaps_col,  # Gaps por indicator
         fontsize_col_color = 7,
         color_labels_col = col_colors,
         main = "KEGG Modules - Ordered by Indicator > Origin > Taxonomy",
         filename = "MetQy_mcf_ordered_indicator_origin_taxonomy.png")

############################# Analysis of relevant MAGs

#####################MAGs poorly described with little metabolic potential

final_table <- read_tsv("all_RELEVANT_MAGs_unified_taxonomy_LOG_NEW.tsv") %>%
  mutate(
    taxon = case_when(
      str_detect(classification, "s__[^;]+") ~ str_extract(classification, "s__[^;]+"),
      str_detect(classification, "g__[^;]+") ~ str_extract(classification, "g__[^;]+"),
      str_detect(classification, "f__[^;]+") ~ str_extract(classification, "f__[^;]+"),
      TRUE ~ "Unclassified"
    ),
    short_label = paste0(taxon, " (", Genome, ")")
  )
View(final_table)


# MAGs with low metabolic potential that cluster
low_potential_mags <- c(
  "MEGAHIT-MetaBAT2-S971.164", "MEGAHIT-MetaBAT2-S960.129", "MEGAHIT-MetaBAT2-S970.243",
  "MEGAHIT-MetaBAT2-S946.39", "MEGAHIT-MetaBAT2-S969.165", "MEGAHIT-MetaBAT2-S981.151",
  "MEGAHIT-MetaBAT2-S944.97", "group7_bins.369", "MEGAHIT-MetaBAT2-S962.71", 
  "MEGAHIT-MetaBAT2-S968.197", "962_SemiBin_32",
  "MEGAHIT-MetaBAT2-S988.19", "973_SemiBin_207", "MEGAHIT-MetaBAT2-S992.120",
  "MEGAHIT-MetaBAT2-S969.185", "group8_bins.297", "MEGAHIT-MetaBAT2-S946.206",
  "983_SemiBin_128", "MEGAHIT-MetaBAT2-S933.139"
)
# Important Treponema and MAGs - CORRECTED
important_mags <- c(
  "932_SemiBin_153", "935_SemiBin_192", "970_SemiBin_111", 
  "981_SemiBin_196", "986_SemiBin_166", "991_SemiBin_49",
  "group6_bins.106", "group9_bins.244", "MEGAHIT-MetaBAT2-S941.86",
  "MEGAHIT-MetaBAT2-S945.92", "MEGAHIT-MetaBAT2-S945.140",
  "MEGAHIT-MetaBAT2-S947.49", "MEGAHIT-MetaBAT2-S947.6",
  "MEGAHIT-MetaBAT2-S948.142", "MEGAHIT-MetaBAT2-S964.185",
  "MEGAHIT-MetaBAT2-S972.52", "MEGAHIT-MetaBAT2-S973.204",
  "MEGAHIT-MetaBAT2-S985.88", "MEGAHIT-MetaBAT2-S932.114",
  "983_SemiBin_321"
)

# Total

all_relevant_mags <- c(
  "970_SemiBin_111", "991_SemiBin_49", "981_SemiBin_196", "935_SemiBin_192",
  "group9_bins.244", "MEGAHIT-MetaBAT2-S948.142", "MEGAHIT-MetaBAT2-S947.49",
  "932_SemiBin_153", "group6_bins.106", "MEGAHIT-MetaBAT2-S973.204",
  "MEGAHIT-MetaBAT2-S972.52", "986_SemiBin_166", "MEGAHIT-MetaBAT2-S945.92",
  "MEGAHIT-MetaBAT2-S945.140", "MEGAHIT-MetaBAT2-S985.88", "MEGAHIT-MetaBAT2-S932.114",
  "MEGAHIT-MetaBAT2-S947.6", "MEGAHIT-MetaBAT2-S941.86", "MEGAHIT-MetaBAT2-S964.185",
  "MEGAHIT-MetaBAT2-S969.185", "983_SemiBin_128", "973_SemiBin_207", "MEGAHIT-MetaBAT2-S946.206",
  "MEGAHIT-MetaBAT2-S988.19", "MEGAHIT-MetaBAT2-S992.120", "MEGAHIT-MetaBAT2-S933.139", "983_SemiBin_321",
  "MEGAHIT-MetaBAT2-S971.164", "MEGAHIT-MetaBAT2-S960.129", "MEGAHIT-MetaBAT2-S970.243",
  "MEGAHIT-MetaBAT2-S946.39", "MEGAHIT-MetaBAT2-S969.165", "MEGAHIT-MetaBAT2-S981.151",
  "MEGAHIT-MetaBAT2-S944.97", "group7_bins.369", "MEGAHIT-MetaBAT2-S962.71", "MEGAHIT-MetaBAT2-S968.197", "962_SemiBin_32")

# Use the corrected list
module_ids<- c(
  "M00009", "M00010", "M00011", "M00015", "M00016", "M00017", 
  "M00018", "M00019", "M00021", "M00023", "M00026", "M00028", "M00029", 
  "M00037", "M00038", "M00045", "M00082", "M00083", "M00086", "M00093", 
  "M00118", "M00119", "M00124", "M00126", "M00131", "M00133", "M00134","M00135",
  "M00222", "M00237", "M00247", "M00250", "M00357", "M00358", "M00499",
  "M00554", "M00569", "M00572", "M00579", "M00609", "M00620", "M00627",
  "M00628", "M00631", "M00632", "M00647", "M00700", "M00706", "M00707", 
  "M00727", "M00728", "M00740", "M00742", "M00793", "M00840"
)

module_ids <- c(
  # Neurotransmitters and precursors
  "M00023",  # Tryptophan biosynthesis (precursor de serotonina)
  "M00037",  # Melatonin biosynthesis (triptófano -> serotonina -> melatonina)
  "M00038",  # Tryptophan metabolism (vía kynurenine)
  "M00135",  # GABA biosynthesis (neurotransmisor inhibitorio)
  
  # SCFA
  "M00082",  # Fatty acid biosynthesis, initiation
  "M00083",  # Fatty acid biosynthesis, elongation
  "M00086",  # beta-Oxidation, acyl-CoA synthesis
  "M00579",  # Acetate production (acetato es SCFA)
  
  # Neurotransmitter-related amino acid metabolism
  "M00015",  # Proline biosynthesis
  "M00016",  # Lysine biosynthesis
  "M00017",  # Methionine biosynthesis
  "M00018",  # Threonine biosynthesis
  "M00019",  # Valine/isoleucine biosynthesis
  "M00021",  # Cysteine biosynthesis
  "M00026",  # Histidine biosynthesis (precursor de histamina)
  "M00028",  # Ornithine biosynthesis (relacionado con urea/glutamato)
  "M00029",  # Urea cycle (detoxificación amonio)
  "M00045",  # Histidine degradation (a histamina/glutamato)
  
  # Cofactors and antioxidants
  "M00118",  # Glutathione biosynthesis (antioxidante clave)
  "M00119",  # Pantothenate biosynthesis (vitamina B5, cofactor)
  "M00124",  # Pyridoxal biosynthesis (vitamina B6, cofactor neurotransmisores)
  "M00126",  # Tetrahydrofolate biosynthesis (folato, metabolismo 1-carbon)
  "M00840",  # Tetrahydrofolate biosynthesis alternativa
  
  # Polyamines and related metabolism
  "M00131",  # Inositol phosphate metabolism
  "M00133",  # Polyamine biosynthesis (putrescine -> spermidine)
  "M00134",  # Polyamine biosynthesis alternativa
  
  # Transportation and stress
  "M00222",  # Phosphate transport system
  "M00237",  # Branched-chain amino acid transport system
  "M00247",  # Putative ABC transport system
  "M00250",  # Lipopolysaccharide transport system (LPS, inflamación)
  "M00499",  # HydH-HydG (metal tolerance, estrés)
  
  # Methanogenesis (affects SCFA/gas production)
  "M00357",  # Methanogenesis, acetate => methane
  "M00358"  # Coenzyme M biosynthesis
  
)

# Find the full names
all_module_names <- colnames(mcf)
specific_modules <- c()

for (module_id in module_ids) {
  pattern <- paste0("^", module_id, " - ")
  matches <- grep(pattern, all_module_names, value = TRUE)
  if (length(matches) > 0) {
    specific_modules <- c(specific_modules, matches[1])
  }
}

cat("Specific modules found:", length(specific_modules), "\n")

# Get full module names
all_module_names <- colnames(mcf)
specific_modules <- find_module_names(module_ids, all_module_names)

cat("Modules found:", length(specific_modules), "of", length(module_ids), "\n")
print(specific_modules)


# Heatmap 1: MAGs with low metabolic potential
low_potential_table <- final_table %>% filter(Genome %in% low_potential_mags)
View(low_potential_table)
if (nrow(low_potential_table) > 0 & length(specific_modules) > 0) {
  low_potential_mcf <- mcf[low_potential_table$Genome, specific_modules, drop = FALSE]
  low_potential_mcf[is.na(low_potential_mcf)] <- 0
  
  low_potential_table_ordered <- low_potential_table %>%
    arrange(Consistent_Indicator, Origins, taxon)
  rownames(low_potential_mcf) <- low_potential_table_ordered$short_label
  
  annotation_low <- data.frame(
    Indicator = low_potential_table_ordered$Consistent_Indicator,
    Origin = low_potential_table_ordered$Origins,
    row.names = low_potential_table_ordered$short_label
  )
  
  pheatmap(t(low_potential_mcf),
           cellwidth = 12,
           cellheight = 8,
           fontsize_row = 7,
           fontsize_col = 7,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           annotation_col = annotation_low,
           annotation_colors = annotation_colors,
           main = "Low metabolic potential MAGs - Specific modules",
           filename = "MetQy_low_potential_specific_GUTBRAIN.png")
} else {
  cat("There is not enough data for a low-potential heatmap\n")
}

# Heatmap 2: Important Treponema and MAGs
important_table <- final_table %>% filter(Genome %in% important_mags)
View(important_table)
if (nrow(important_table) > 0 & length(specific_modules) > 0) {
  important_mcf <- mcf[important_table$Genome, specific_modules, drop = FALSE]
  important_mcf[is.na(important_mcf)] <- 0
  
  important_table_ordered <- important_table %>%
    arrange(Consistent_Indicator, Origins, taxon)
  rownames(important_mcf) <- important_table_ordered$short_label
  
  annotation_important <- data.frame(
    Indicator = important_table_ordered$Consistent_Indicator,
    Origin = important_table_ordered$Origins,
    row.names = important_table_ordered$short_label
  )
  
  pheatmap(t(important_mcf),
           cellwidth = 12,
           cellheight = 8,
           fontsize_row = 7,
           fontsize_col = 7,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           annotation_col = annotation_important,
           annotation_colors = annotation_colors,
           main = "Treponema & Important MAGs - Specific modules",
           filename = "MetQy_treponema_specific_GUTBRAIN.pdf")
} else {
  cat("There is not enough data for the Treponema and important MAGs heatmap\n")
}


######with cluster

if (nrow(low_potential_table) > 0 & length(specific_modules) > 0) {
  low_potential_mcf <- mcf[low_potential_table$Genome, specific_modules, drop = FALSE]
  low_potential_mcf[is.na(low_potential_mcf)] <- 0
  
  low_potential_table_ordered <- low_potential_table %>%
    arrange(Consistent_Indicator, Origins, taxon)
  rownames(low_potential_mcf) <- low_potential_table_ordered$short_label
  
  annotation_low <- data.frame(
    Indicator = low_potential_table_ordered$Consistent_Indicator,
    Origin = low_potential_table_ordered$Origins,
    row.names = low_potential_table_ordered$short_label
  )
  
  pheatmap(t(low_potential_mcf),
           cellwidth = 12,
           cellheight = 8,
           fontsize_row = 7,
           fontsize_col = 7,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = annotation_low,
           annotation_colors = annotation_colors,
           main = "Low metabolic potential MAGs - Specific modules",
           filename = "MetQy_low_potential_specific_cluster_GUTBRAIN.pdf")
} else {
  cat("There is not enough data for a low-potential heatmap\n")
}


# Alternative version with clustering to see patterns
if (nrow(important_table) > 0 & length(specific_modules) > 0) {
  important_mcf <- mcf[important_table$Genome, specific_modules, drop = FALSE]
  important_mcf[is.na(important_mcf)] <- 0
  
  important_table_ordered <- important_table %>%
    arrange(Consistent_Indicator, Origins, taxon)
  rownames(important_mcf) <- important_table_ordered$short_label
  
  annotation_important <- data.frame(
    Indicator = important_table_ordered$Consistent_Indicator,
    Origin = important_table_ordered$Origins,
    row.names = important_table_ordered$short_label
  )
  
  pheatmap(t(important_mcf),
           cellwidth = 12,
           cellheight = 8,
           fontsize_row = 7,
           fontsize_col = 7,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = annotation_important,
           annotation_colors = annotation_colors,
           main = "Treponema & Important MAGs - Specific modules",
           filename = "MetQy_treponema_specific_cluster_GUTBRAIN.pdf")
} else {
  cat("There is not enough data for the Treponema and important MAGs heatmap\n")
}





########################Complete heatmaps for the subsets

# Heatmap 1: MAGs with low metabolic potential - ALL modules present
low_potential_table <- final_table %>% filter(Genome %in% low_potential_mags)

if (nrow(low_potential_table) > 0) {
  low_potential_mcf <- mcf[low_potential_table$Genome, ]
  low_potential_mcf[is.na(low_potential_mcf)] <- 0
  
  # Filter modules that appear in at least 1 MAG (values> 0)
  modules_present <- colSums(low_potential_mcf > 0) >= 1
  low_potential_filtered <- low_potential_mcf[, modules_present, drop = FALSE]
  
  cat("Low potential MAGs:", nrow(low_potential_filtered), "\n")
  cat("Modules present:", ncol(low_potential_filtered), "\n")
  
  low_potential_table_ordered <- low_potential_table %>%
    arrange(Consistent_Indicator, Origins, taxon)
  rownames(low_potential_filtered) <- low_potential_table_ordered$short_label
  
  annotation_low <- data.frame(
    Indicator = low_potential_table_ordered$Consistent_Indicator,
    Origin = low_potential_table_ordered$Origins,
    row.names = low_potential_table_ordered$short_label
  )
  
  pheatmap(t(low_potential_filtered),
           cellwidth = 12,
           cellheight = 8,
           fontsize_row = 6,
           fontsize_col = 7,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           annotation_col = annotation_low,
           annotation_colors = annotation_colors,
           main = "Low metabolic potential MAGs - All modules",
           filename = "MetQy_low_potential_all_modules.png")
} else {
  cat("No low-potential MAGs were found\n")
}

# Heatmap 2: Important Treponema and MAGs - ALL modules present
important_table <- final_table %>% filter(Genome %in% important_mags)

if (nrow(important_table) > 0) {
  important_mcf <- mcf[important_table$Genome, ]
  important_mcf[is.na(important_mcf)] <- 0
  
  # Filter modules that appear in at least 1 MAG (values ​​> 0)
  modules_present <- colSums(important_mcf > 0) >= 1
  important_filtered <- important_mcf[, modules_present, drop = FALSE]
  
  cat("Important MAGs:", nrow(important_filtered), "\n")
  cat("Modules present:", ncol(important_filtered), "\n")
  
  important_table_ordered <- important_table %>%
    arrange(Consistent_Indicator, Origins, taxon)
  rownames(important_filtered) <- important_table_ordered$short_label
  
  annotation_important <- data.frame(
    Indicator = important_table_ordered$Consistent_Indicator,
    Origin = important_table_ordered$Origins,
    row.names = important_table_ordered$short_label
  )
  
  pheatmap(t(important_filtered),
           cellwidth = 12,
           cellheight = 8,
           fontsize_row = 6,
           fontsize_col = 7,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           annotation_col = annotation_important,
           annotation_colors = annotation_colors,
           main = "Treponema & Important MAGs - All modules",
           filename = "MetQy_treponema_all_modules.pdf")
} else {
  cat("No relevant MAGs were found\n")
}

# Alternative version with clustering to see patterns
if (nrow(low_potential_table) > 0) {
  pheatmap(t(low_potential_filtered),
           cellwidth = 12,
           cellheight = 8,
           fontsize_row = 6,
           fontsize_col = 7,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = annotation_low,
           annotation_colors = annotation_colors,
           main = "Low metabolic potential MAGs - All modules (clustered)",
           filename = "MetQy_low_potential_all_modules_clustered.pdf")
}

if (nrow(important_table) > 0) {
  pheatmap(t(important_filtered),
           cellwidth = 12,
           cellheight = 8,
           fontsize_row = 6,
           fontsize_col = 7,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = annotation_important,
           annotation_colors = annotation_colors,
           main = "Treponema & Important MAGs - All modules (clustered)",
           filename = "MetQy_treponema_all_modules_clustered.pdf")
}

# Also export the data for analysis
if (nrow(low_potential_table) > 0) {
  write_tsv(as.data.frame(low_potential_filtered) %>% 
              rownames_to_column("MAG"),
            "low_potential_all_modules_data.tsv")
}

if (nrow(important_table) > 0) {
  write_tsv(as.data.frame(important_filtered) %>% 
              rownames_to_column("MAG"),
            "treponema_all_modules_data.tsv")
}


########################################RECOVERING THE HEATMAP FROM 99 MAGS, BUT THIS TIME CLUSTERING IT

########## Heatmap of the 99 MAGs with clustering

final_table <- read_tsv("all_RELEVANT_MAGs_unified_taxonomy_LOG.tsv") %>%
  mutate(
    taxon = case_when(
      str_detect(classification, "s__[^;]+") ~ str_extract(classification, "s__[^;]+"),
      str_detect(classification, "g__[^;]+") ~ str_extract(classification, "g__[^;]+"),
      str_detect(classification, "f__[^;]+") ~ str_extract(classification, "f__[^;]+"),
      TRUE ~ "Unclassified"
    ),
    short_label = paste0(taxon, " (", Genome, ")")
  )

# Sort FIRST by Indicator, then by taxonomy, and finally by origin
final_table_ordered <- final_table %>%
  arrange(Consistent_Indicator, Origins, taxon)

mcf_ordered <- mcf[final_table_ordered$Genome, ]
rownames(mcf_ordered) <- final_table_ordered$short_label
mcf_ordered[is.na(mcf_ordered)] <- 0

# FILTER MODULES: keep only those present in at least 2 MAGs
min_mags_with_module <- 2
mcf_filtered_modules_strict <- mcf_ordered[, colSums(mcf_ordered > 0) >= min_mags_with_module]

# Create annotations for the COLUMNS (MAGs) - TWO INDICATORS
annotation_col <- data.frame(
  Indicator = final_table_ordered$Consistent_Indicator,
  Origin = final_table_ordered$Origins,
  row.names = final_table_ordered$short_label
)

# Define colors for both indicators
annotation_colors <- list(
  Indicator = c(
    "control" = "blue",
    "stress" = "red", 
    "mixed" = "purple"
  ),
  Origin = c(
    "MAGG_PLSDA_log" = "lightgreen",
    "MAGG_volcano" = "darkgreen",
    "MAGG_PLSDA_log; MAGG_volcano" = "turquoise",
    "MAGP_PLSDA_log" = "orange",
    "MAGP_volcano" = "darkorange4",
    "MAGP_PLSDA_log; MAGP_volcano" = "yellow",
    "MAGG_PLSDA_log; MAGP_PLSDA_log" = "pink",
    "MAGG_PLSDA_log; MAGP_volcano" = "coral",
    "MAGG_PLSDA_log; MAGP_PLSDA_log; MAGP_volcano" = "magenta"
  )
)

# Heatmap with Clustering
pheatmap(t(mcf_filtered_modules_strict),
         cellwidth = 10,
         cellheight = 8,
         fontsize_row = 6,
         fontsize_col = 7,
         cluster_rows = TRUE,    # CLUSTERING de módulos
         cluster_cols = TRUE,    # CLUSTERING de MAGs
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         annotation_names_col = TRUE,
         annotation_legend = TRUE,
         show_colnames = TRUE,
         main = "KEGG Modules - 99 MAGs with clustering",
         filename = "MetQy_99MAGs_clustered.pdf")

# Version without clustering for comparison (maintaining the original order)
pheatmap(t(mcf_filtered_modules_strict),
         cellwidth = 8,
         cellheight = 6,
         fontsize_row = 5,
         fontsize_col = 6,
         cluster_rows = FALSE,   # Sin clustering
         cluster_cols = FALSE,   # Sin clustering
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         annotation_names_col = TRUE,
         annotation_legend = TRUE,
         show_colnames = TRUE,
         gaps_col = cumsum(rle(final_table_ordered$Consistent_Indicator)$lengths)[-length(rle(final_table_ordered$Consistent_Indicator)$lengths)],
         main = "KEGG Modules - 99 MAGs without clustering",
         filename = "MetQy_99MAGs_no_clustering.pdf")

# Export the clustering data
write_tsv(as.data.frame(mcf_filtered_modules_strict) %>% 
            rownames_to_column("MAG_label"),
          "99MAGs_module_data.tsv")

# Information about clustering
cat("Total MAGs:", nrow(mcf_filtered_modules_strict), "\n")
cat("Total modules (present in ≥2 MAGs):", ncol(mcf_filtered_modules_strict), "\n")
cat("MAGs by group:\n")
print(table(final_table_ordered$Consistent_Indicator))






########################99 MAGs but only modules of interest.

########## Heatmap of the 99 MAGs with clustering - GUT-BRAIN AXIS MODULES

final_table <- read_tsv("all_RELEVANT_MAGs_unified_taxonomy_LOG.tsv") %>%
  mutate(
    taxon = case_when(
      str_detect(classification, "s__[^;]+") ~ str_extract(classification, "s__[^;]+"),
      str_detect(classification, "g__[^;]+") ~ str_extract(classification, "g__[^;]+"),
      str_detect(classification, "f__[^;]+") ~ str_extract(classification, "f__[^;]+"),
      TRUE ~ "Unclassified"
    ),
    short_label = paste0(taxon, " (", Genome, ")")
  )

# Sort FIRST by Indicator, then by taxonomy, and finally by origin
final_table_ordered <- final_table %>%
  arrange(Consistent_Indicator, Origins, taxon)

mcf_ordered <- mcf[final_table_ordered$Genome, ]
rownames(mcf_ordered) <- final_table_ordered$short_label
mcf_ordered[is.na(mcf_ordered)] <- 0

# GUT-BRAIN AXIS
module_ids <- c(
  "M00023",  # Tryptophan biosynthesis (precursor de serotonina)
  "M00037",  # Melatonin biosynthesis (triptófano -> serotonina -> melatonina)
  "M00038",  # Tryptophan metabolism (vía kynurenine)
  "M00135",  # GABA biosynthesis (neurotransmisor inhibitorio)
  
  "M00082",  # Fatty acid biosynthesis, initiation
  "M00083",  # Fatty acid biosynthesis, elongation
  "M00086",  # beta-Oxidation, acyl-CoA synthesis
  "M00579",  # Acetate production (acetato es SCFA)
  
  "M00015",  # Proline biosynthesis
  "M00016",  # Lysine biosynthesis
  "M00017",  # Methionine biosynthesis
  "M00018",  # Threonine biosynthesis
  "M00019",  # Valine/isoleucine biosynthesis
  "M00021",  # Cysteine biosynthesis
  "M00026",  # Histidine biosynthesis (precursor de histamina)
  "M00028",  # Ornithine biosynthesis (relacionado con urea/glutamato)
  "M00029",  # Urea cycle (detoxificación amonio)
  "M00045",  # Histidine degradation (a histamina/glutamato)

  "M00118",  # Glutathione biosynthesis (antioxidante clave)
  "M00119",  # Pantothenate biosynthesis (vitamina B5, cofactor)
  "M00124",  # Pyridoxal biosynthesis (vitamina B6, cofactor neurotransmisores)
  "M00126",  # Tetrahydrofolate biosynthesis (folato, metabolismo 1-carbon)
  "M00840",  # Tetrahydrofolate biosynthesis alternativa
  
  "M00131",  # Inositol phosphate metabolism
  "M00133",  # Polyamine biosynthesis (putrescine -> spermidine)
  "M00134",  # Polyamine biosynthesis alternativa
  
  "M00222",  # Phosphate transport system
  "M00237",  # Branched-chain amino acid transport system
  "M00247",  # Putative ABC transport system
  "M00250",  # Lipopolysaccharide transport system (LPS, inflamación)
  "M00499",  # HydH-HydG (metal tolerance, estrés)
  
  "M00357",  # Methanogenesis, acetate => methane
  "M00358"   # Coenzyme M biosynthesis
)

# Find the full names of the modules
all_module_names <- colnames(mcf_ordered)
gutbrain_modules <- c()

for (module_id in module_ids) {
  pattern <- paste0("^", module_id, " - ")
  matches <- grep(pattern, all_module_names, value = TRUE)
  if (length(matches) > 0) {
    gutbrain_modules <- c(gutbrain_modules, matches[1])
  } else {
    cat("Módulo no encontrado:", module_id, "\n")
  }
}

cat("Gut-Brain Axis Modules found:", length(gutbrain_modules), "\n")

# Filter the array to show only these modules
mcf_gutbrain <- mcf_ordered[, gutbrain_modules, drop = FALSE]

# Filter modules that appear in at least 2 MAGs
mcf_gutbrain_filtered <- mcf_gutbrain[, colSums(mcf_gutbrain > 0) >= 2, drop = FALSE]

cat("Gut-Brain Axis modules present in ≥2 MAGs:", ncol(mcf_gutbrain_filtered), "\n")

# Create annotations for the COLUMNS (MAGs)
annotation_col <- data.frame(
  Indicator = final_table_ordered$Consistent_Indicator,
  Origin = final_table_ordered$Origins,
  row.names = final_table_ordered$short_label
)

# Define colors for both indicators
annotation_colors <- list(
  Indicator = c(
    "control" = "blue",
    "stress" = "red", 
    "mixed" = "purple"
  ),
  Origin = c(
    "MAGG_PLSDA_log" = "lightgreen",
    "MAGG_volcano" = "darkgreen",
    "MAGG_PLSDA_log; MAGG_volcano" = "turquoise",
    "MAGP_PLSDA_log" = "orange",
    "MAGP_volcano" = "darkorange4",
    "MAGP_PLSDA_log; MAGP_volcano" = "yellow",
    "MAGG_PLSDA_log; MAGP_PLSDA_log" = "pink",
    "MAGG_PLSDA_log; MAGP_volcano" = "coral",
    "MAGG_PLSDA_log; MAGP_PLSDA_log; MAGP_volcano" = "magenta"
  )
)

# Heatmap WITH CLUSTERING - GUT-BRAIN AXIS MODULES
pheatmap(t(mcf_gutbrain_filtered),
         cellwidth = 12,
         cellheight = 10, 
         fontsize_row = 7,
         fontsize_col = 7,
         cluster_rows = TRUE,    # CLUSTERING modules
         cluster_cols = TRUE,    # CLUSTERING MAGs
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         annotation_names_col = TRUE,
         annotation_legend = TRUE,
         show_colnames = TRUE,
         main = "Gut-Brain Axis Modules - 99 MAGs with clustering",
         filename = "MetQy_99MAGs_gutbrain_clustered.png")

# Non-clustered version for comparison
pheatmap(t(mcf_gutbrain_filtered),
         cellwidth = 12,
         cellheight = 10,
         fontsize_row = 7,
         fontsize_col = 7,
         cluster_rows = FALSE,  
         cluster_cols = FALSE,   
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         annotation_names_col = TRUE,
         annotation_legend = TRUE,
         show_colnames = TRUE,
         gaps_col = cumsum(rle(final_table_ordered$Consistent_Indicator)$lengths)[-length(rle(final_table_ordered$Consistent_Indicator)$lengths)],
         main = "Gut-Brain Axis Modules - 99 MAGs without clustering",
         filename = "MetQy_99MAGs_gutbrain_no_clustering.png")

# Export data
write_tsv(as.data.frame(mcf_gutbrain_filtered) %>% 
            rownames_to_column("MAG_label"),
          "99MAGs_gutbrain_module_data.tsv")

# Information about the Gut-Brain Axis modules
cat("Total MAGs:", nrow(mcf_gutbrain_filtered), "\n")
cat("Gut-Brain Axis Modules (present in ≥2 MAGs):", ncol(mcf_gutbrain_filtered), "\n")
cat("MAGs by group:\n")
print(table(final_table_ordered$Consistent_Indicator))

# Show included modules
cat("\nGut-Brain Axis modules included in the heatmap:\n")
print(colnames(mcf_gutbrain_filtered))

# Completeness statistics by group
gutbrain_completeness <- data.frame(
  MAG = rownames(mcf_gutbrain_filtered),
  Indicator = final_table_ordered$Consistent_Indicator,
  Total_Modules = ncol(mcf_gutbrain_filtered),
  Present_Modules = rowSums(mcf_gutbrain_filtered > 0),
  Completeness_Percent = rowSums(mcf_gutbrain_filtered > 0) / ncol(mcf_gutbrain_filtered) * 100
)

cat("\nAverage completion rate of Gut-Brain Axis modules per group:\n")
print(gutbrain_completeness %>% 
        group_by(Indicator) %>% 
        summarise(Mean_Completeness = mean(Completeness_Percent),
                  SD_Completeness = sd(Completeness_Percent)))




