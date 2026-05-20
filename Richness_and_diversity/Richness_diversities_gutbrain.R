######Richness
#To calculate richness from proteomics data we need to count the number of quantifiable species per sample.
#You have the abundances already, you just need to transfer those to a presence/absence matrix, 1 if it was quantified, 0 if not. 
#Then from this we can calculate the richness measures. 
#Alpha diversity (Shannon index) and beta diversity (Bray-Curtis) can be calculated in R, either manually or using the package vegan. 

#libraries
library(tidyverse)
library(readxl)
library(writexl)
library(vegan)
library(pheatmap)
library(ggplot2)
library(reshape2)

#MICROECO at the end

################################################################### VEGAN
#INPUTS
MAGG_abundance <- read.table("MAGG_abundance.tsv", header = TRUE, sep = "\t")
MAGP_abundance <- read.table("MAGP_abundance_raw_New.txt", header = TRUE, sep = "\t")
MAGP_abundance_log <- read.table("MAGP_abundance_log_New.txt", header = TRUE, sep = "\t") 
#MAGP_abundance_clr <- read.table("MAGP_abundance_clr_New.txt", header = TRUE, sep = "\t")

#Correction to avoid 0 in heatmap

MAGP_abundance_log_corrected <- MAGP_abundance_log %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = ~ if_else(.x == 0, rnorm(n(), mean = 8, sd = 1), as.numeric(.x))
  ))
View(MAGP_abundance_log_corrected)

metadata <-  read.table("metadataMetaP.txt", header = TRUE, sep = "\t") %>%
  rename(treatment = Tto) %>%    
  rename(pen = corral)

## 1. Calculate Richness ------------------------------------------------------
# Convert to presence/absence matrices (1/0)
MAGP_log_presence <- MAGP_abundance_log %>% mutate(across(-1, ~if_else(. > 0, 1, 0)))
MAGG_presence <- MAGG_abundance %>% mutate(across(-1, ~if_else(. > 0, 1, 0)))

# Calculate richness (number of present species per sample)
richness_proteomics_log <- MAGP_log_presence %>% 
  select(-1) %>% 
  colSums() %>% 
  enframe(name = "Sample", value = "Richness")

richness_genomics <- MAGG_presence %>% 
  select(-1) %>% 
  colSums() %>% 
  enframe(name = "Sample", value = "Richness")

## 2. Calculate Diversity Indices ---------------------------------------------
# Alpha diversity (Shannon index)
shannon_proteomics_log <- MAGP_abundance_log_corrected %>% 
  select(-1) %>% 
  t() %>% 
  diversity(index = "shannon") %>% 
  enframe(name = "Sample", value = "Shannon_Index")

shannon_genomics <- MAGG_abundance %>% 
  select(-1) %>% 
  t() %>% 
  diversity(index = "shannon") %>% 
  enframe(name = "Sample", value = "Shannon_Index")

# Beta diversity (Bray-Curtis dissimilarity)
bray_proteomics_raw <- MAGP_abundance %>%
  select(-1) %>% 
  t() %>% 
  vegdist(method = "bray") %>% 
  as.matrix() %>% 
  as.data.frame()

# Beta diversity (Bray-Curtis dissimilarity)
bray_proteomics_log <- MAGP_abundance_log_corrected %>% #USING CORRECTION
  select(-1) %>% 
  t() %>% 
  vegdist(method = "bray") %>% 
  as.matrix() %>% 
  as.data.frame()

#bray_proteomics_clr <- MAGP_abundance_clr %>% 
 # select(-1) %>% 
 # t() %>% 
 # vegdist(method = "bray") %>% 
  #as.matrix() %>% 
  #as.data.frame()

bray_genomics <- MAGG_abundance %>% 
  select(-1) %>% 
  t() %>% 
  vegdist(method = "bray") %>% 
  as.matrix() %>% 
  as.data.frame()

## 3. Save All Results --------------------------------------------------------
# Create a list with all results
results <- list(
  Proteomics_Richness_log = richness_proteomics_log,
  Genomics_Richness = richness_genomics,
  Proteomics_Shannon_log = shannon_proteomics_log,
  Genomics_Shannon = shannon_genomics,
  Proteomics_BrayCurtis_log = bray_proteomics_log,
  Genomics_BrayCurtis = bray_genomics
)

# Save to Excel file
write_xlsx(results, "Diversity_Analysis_Results_MAGGraw_MAGPlog_corrected.xlsx")

# Also save as RData for further analysis
save(results, file = "Diversity_Analysis_Results_MAGGraw_MAGPlog_corrected.RData")


## 4. Visualization ------------------------------------------------------------


################# SIMPLE HEATMAP

# For proteomics data
pheatmap(as.matrix(bray_proteomics_log),
         main = "MAGP LOG Sample Similarity (Bray-Curtis)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         clustering_distance_rows = "euclidean",  # Cluster using Euclidean on the dissimilarity matrix
         clustering_distance_cols = "euclidean",
         border_color = NA)

# For genomics data
pheatmap(as.matrix(bray_genomics),
         main = "MAGG Sample Similarity (Bray-Curtis)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         border_color = NA)

################### HEATMAP considering metadata

# Define a consistent color palette for both plots
my_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
annotation_colors <- list(
  treatment = c("control" = "green3", "stress" = "orange"),
  sex = c("F" = "pink", "M" = "purple"),
  pen = RColorBrewer::brewer.pal(max(length(unique(metadata$pen))), "Set3")
)

names(annotation_colors$pen) <- unique(metadata$pen)

# visualisation function
create_annotated_heatmap <- function(bray_matrix, title_prefix, metadata) {

  sample_names <- rownames(as.matrix(bray_matrix))

# Prepare metadata annotation
  annotation_df <- metadata %>%
    filter(SampleID %in% sample_names) %>%
    arrange(match(SampleID, sample_names)) %>%
    column_to_rownames("SampleID") %>%
    select(treatment, pen, sex)  
    
# Create heatmap
  pheatmap(
    as.matrix(bray_matrix),
    main = paste(title_prefix, "Sample Similarity"),
    color = my_palette,
    annotation_row = annotation_df,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    border_color = NA,
    fontsize = 9,
    cellwidth = 12,
    cellheight = 12
    )
  }
  
  # Create plots for both datasets
create_annotated_heatmap(bray_proteomics_log, "MAGP LOG (Proteomics)", metadata)
create_annotated_heatmap(bray_proteomics_raw, "MAGP (Proteomics)", metadata)
create_annotated_heatmap(bray_genomics, "MAGG (Genomics)", metadata)

create_annotated_heatmap(bray_proteomics_clr, "MAGP CLR (Proteomics)", metadata)

##OUTPUT

png("heatmap_proteomics_log.png", width = 1200, height = 1000)
create_annotated_heatmap(bray_proteomics_log, "MAGP LOG (Proteomics)", metadata)
dev.off()

png("heatmap_proteomics_raw.png", width = 1200, height = 1000)
create_annotated_heatmap(bray_proteomics_raw, "MAGP (Proteomics)", metadata)
dev.off()

png("heatmap_genomics.png", width = 1200, height = 1200)
create_annotated_heatmap(bray_genomics, "MAGG (Genomics)", metadata)
dev.off()

######################################################USING MICROECO
library(microeco)
library(magrittr)
library(tidyverse)
library(readxl)

#INPUT
MAGG_abundance <- read.table("MAGG_abundance.tsv", header = TRUE, sep = "\t")
MAGP_abundance <- read.table("MAGP_abundance_raw_New.txt", header = TRUE, sep = "\t")
MAGP_abundance_log <- read.table("MAGP_abundance_log_New.txt", header = TRUE, sep = "\t") 
#MAGP_abundance_clr <- read.table("MAGP_abundance_clr_New.txt", header = TRUE, sep = "\t")

metadata <-  read.table("metadataMetaP.txt", header = TRUE, sep = "\t") %>%
  rename(treatment = Tto) %>%    
  rename(pen = corral) %>%
  t()
taxonomy <- read_excel("C:/Users/RRIO/OneDrive - IRTA/Escritorio/DOCTORAT/SEGUNDO AÑO/Estancia/DATA for the stay/protein_legend, metadata, all_proteins abundance and MAG taxonomy/taxonomia_MetaP2.xlsx")

View(MAGP_abundance_log)
##### Prepare the data 
# Samples as rows, MAGs as columns

abundance_table <- MAGP_abundance %>% #####Change for MAGG_abundance or MAGP_abundance_log or MAGP_abundance_clr
  as.data.frame() %>% 
  remove_rownames() %>%
  column_to_rownames("Bin")

# Make sure sample names match between abundance and metadata
colnames(metadata) <- colnames(abundance_table) #recharge metadata and taxonomy before this if you've changed abundance input
metadata <- metadata[-1, ] %>%
  as.data.frame()

# Format taxonomy
taxonomy_table <- taxonomy %>% 
  select(user_genome, classification) %>% 
  separate(classification, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% 
  column_to_rownames("user_genome")

##### Create microtable object
dataset <- microtable$new(
  sample_table = metadata,
  otu_table = abundance_table,
  tax_table = taxonomy_table
)

##### Calculate alpha diversity indices (Shannon, Simpson, Chao1, etc.)

dataset$cal_alphadiv(PD = FALSE)  # Set PD = TRUE if you want phylogenetic diversity

head(dataset$alpha_diversity)

# Save results
dir.create("alpha_diversity", showWarnings = FALSE)
dataset$save_alphadiv(dirpath = "alpha_diversity")

##### Calculate beta diversity (Bray-Curtis, Jaccard, UniFrac)

dataset$cal_betadiv(unifrac = FALSE)  # Set unifrac = TRUE if you want UniFrac distances (needs phylogenetic tree)

# View the results (Bray-Curtis dissimilarity matrix)
dataset$beta_diversity$bray

# Save results
dir.create("beta_diversity", showWarnings = FALSE)
dataset$save_betadiv(dirpath = "beta_diversity")


#################Visualization 

print(colnames(abundance_table))
print(rownames(metadata))

# Reload metadata (without t())
metadata <- read.table("metadataMetaP.txt", header = TRUE, sep = "\t") %>%
  rename(treatment = Tto) %>%    
  rename(pen = corral)
rownames(metadata) <- metadata$SampleID

# Remake microtable object
dataset <- microtable$new(
  sample_table = metadata,
  otu_table = abundance_table,
  tax_table = taxonomy_table
)

# Alfa plots
t_alpha <- trans_alpha$new(dataset = dataset, group = "treatment")

# Statistical test (Kruskal-Wallis)
t_alpha$cal_diff(method = "KW")
head(t_alpha$res_diff)

# Plot Shannon diversity
t_alpha$plot_alpha(measure = "Shannon", order_x_mean = TRUE)

# Beta plots
##### Recalculate beta diversity (Bray-Curtis, Jaccard, UniFrac)
dataset$cal_betadiv(unifrac = FALSE)
dataset$beta_diversity$bray

t_beta <- trans_beta$new(dataset = dataset, group = "treatment", measure = "bray")
t_beta$cal_ordination(method = "PCoA")

# Plot PCoA
t_beta$plot_ordination(plot_color = "treatment", plot_shape = "sex")


