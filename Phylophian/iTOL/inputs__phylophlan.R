#######ITOL
library(dplyr)
library(readr)
library(tidyr)
library(tidyverse)
library(writexl)
library(readxl)
library(stringr)

#Import the data needed to create the iTol dataset
MAGG_abundance <- read.table("catalogue_MAGs_MetaP2_abundance.txt", header=T, sep="\t") ##RAW
MAGG_catalogue <- read.table("MAG_C70C10_MetaP2.txt", header=T, sep="\t")
metadata <- read.table("metadataMetaP.txt", header=T, sep="\t")
taxonomy <- read_excel("taxonomia_MetaP2.xlsx") %>%
  select(., user_genome, classification)

#Filter the abundance table by the MAGs catalog used for the metaP and for the tree
MAGG_catalogue$Genomes <- rownames(MAGG_catalogue)
MAGG_abundance <- MAGG_abundance[MAGG_abundance$Genomes %in% MAGG_catalogue$Genomes, ]
nrow(MAGG_abundance)

#Remove	IsollateEukC-trimmed-reads.medaka
MAGG_abundance <- MAGG_abundance[MAGG_abundance$Genomes != "IsollateEukC-trimmed-reads.medaka", ]
View(MAGG_abundance)

write.table(MAGG_abundance, "MAGG_abundance.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#I create a column with the ID of each MAG ####ONLY ONCE, USE THE MATCH SAVED IN MAGG_ID below
MAGG_abundance <- MAGG_abundance %>%
  mutate(
    MAG_ID = sprintf("MAG.%03d", 1:n()),
  ) %>%
  relocate(MAG_ID, .after = 1)


##### LOAD FROM HERE to avoid creating a new legend and losing the original MAG.000 code
##Sqrt transformation on MAGG_abundance USE SAME ID, DO NOT RENAME, USE MAGG_ID generated later
# MAGG_sqrt <- as.data.frame(MAGG_abundance)
# MAGG_sqrt[, -1] <- sqrt(MAGG_sqrt[, -1])
# View(MAGG_sqrt)
#write.table(MAGG_sqrt, "MAGG_sqrt.tsv", sep = "\t", row.names= FALSE, col.names = TRUE, quote = FALSE)

MAGP_abundance_raw <- read.table("MAGP_abundance_raw_New.txt", header=T, sep="\t")
View(MAGP_abundance_raw)

MAGP_sqrt <- as.data.frame(MAGP_abundance_raw)
MAGP_sqrt[, -1] <- sqrt(MAGP_sqrt[, -1])
View(MAGP_sqrt)
write.table(MAGP_sqrt, "MAGP_sqrt.tsv", sep = "\t", row.names= FALSE, col.names = TRUE, quote = FALSE)

MAGP_sqrt <- read.table("MAGP_sqrt.tsv", header=T, sep="\t")

MAGG_ID <- read.table("MAGG_ID.tsv", header=T, sep="\t") 


MAGG_abundance <- read.table("MAGG_abundance.txt", header=T, sep="\t")
MAGG_abundance$MAG_ID <- MAGG_ID$MAG_ID[match(MAGG_abundance$Genomes, MAGG_ID$Genomes)]
MAGG_sqrt <- read.table("MAGG_sqrt.txt", header=T, sep="\t")
MAGG_sqrt$MAG_ID <- MAGG_ID$MAG_ID[match(MAGG_sqrt$Genomes, MAGG_ID$Genomes)] 


#########Calculate the average abundance for control and for stress.

control_samples <- metadata %>%
  filter(Tto == "control") %>%
  pull(SampleID)

stress_samples <- metadata %>%
  filter(Tto == "stress") %>%
  pull(SampleID)

print(control_samples)
print(stress_samples)

MAGP_abundance_sqrt_means <- MAGP_sqrt%>%  ###You can change to MAGG_sqrt, MAGP_sqrt and vice versa
  mutate(
    average_control = rowMeans(select(., all_of(control_samples))),
    average_stress = rowMeans(select(., all_of(stress_samples)))
  ) %>%
  select(., Bin, average_control, average_stress) #Change Bin by Genomes if working with MAGG

View(MAGP_abundance_sqrt_means)
write.table(MAGP_abundance_sqrt_means, "MAGP_abundance_sqrt_means.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#I create one subset for control and another for stress so I can upload them as two different bars in the tree.

MAGG_control_means <- MAGG_abundance_means %>%
  select(., Genomes, average_control, MAG_ID
  )

MAGG_stress_means <- MAGG_abundance_means %>%
  select(., Genomes, average_stress, MAG_ID
  )
View(MAGG_stress_means)

#write.table(MAGG_control_means, "MAGG_control_sqrt_means.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(MAGG_stress_means, "MAGG_stress_sqrt_means.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#### ADJUST the input format for iTOL

MAGG_control_means_input <- MAGG_abundance_means %>%
  select(., Genomes, average_control) 

MAGG_stress_means_input <- MAGG_abundance_means %>%
  select(., Genomes, average_stress)

View(MAGG_stress_means_input)

write.table(MAGG_control_means_input, "MAGG_control_means_input.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(MAGG_stress_means_input, "MAGG_stress_means_input.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Create a table with only the MAG and the ID

MAGG_ID <- MAGG_control_means %>%
  select(., Genomes, MAG_ID)

write.table(MAGG_ID, "MAGG_ID.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE )

#To color the taxonomy and create the legend, we need a file with the names, IDs, and phyla.

tax_filum <- taxonomy %>%
  rename(Genomes = user_genome) %>%  
  mutate(
    Phylum = str_extract(classification, "p__([^;]+)", group = 1),
    Phylum = str_remove(Phylum, "p__")
  ) %>%
  select(Genomes, Phylum) 

View(tax_filum)

####Create the labels input
labels_tree <- MAGG_ID %>%
  left_join(tax_filum, by = "Genomes")

labels_tree <- labels_tree %>%
  mutate(
    Phylum = ifelse(is.na(Phylum), "Unknown", as.character(Phylum))
  )

write.table(labels_tree, "labels_tree.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE )

print(unique(labels_tree$Phylum))
View(labels_tree)
###Create the color_branches input

# 1. Define unique phyla (including NA)
phyla <- unique(labels_tree$Phylum)

# 2. Create a distinct color palette for each edge (using RColorBrewer + additional colors)
set.seed(123)
phyla_colours <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728","#CCCCCC", "#9467bd", "#8c564b", "#e377c2",
  "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
  "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#8c6d31"
)[1:length(phyla)] %>% 
  setNames(phyla)

# 3. Generate the file for iTOL

itol_colours <- labels_tree %>%
  mutate(
    color = ifelse(
      is.na(Phylum) | !Phylum %in% names(phyla_colours),
      "#CCCCCC",
      phyla_colours[as.character(Phylum)]
    ),
    label = ifelse(is.na(Phylum), "Unknown", as.character(Phylum))
  ) %>%
  select(Genomes, color, label)


View(itol_colours)

# 4. Create the LEGEND section
legend_entries <- data.frame(
  phylum = names(phyla_colours),
  color = unname(phyla_colours)) %>%
  mutate(shape = 1,label = ifelse(is.na(phylum), "Unknown", phylum))

# 5. Write the file in iTOL format
output <- c(
  "DATASET_COLORSTRIP",
  "SEPARATOR TAB",
  "DATASET_LABEL Phylum_Colors",
  "COLOR #ff0000",
  "COLOR_BRANCHES 1",
  paste("LEGEND_TITLE", "Phylum Colors"),
  paste("LEGEND_SHAPES", paste(rep(1, nrow(legend_entries)), collapse = " ")),
  paste("LEGEND_COLORS", paste(legend_entries$color, collapse = " ")),
  paste("LEGEND_LABELS", paste(legend_entries$label, collapse = " ")),
  "STRIP_WIDTH 25",
  "DATA",
  paste(itol_colours$Genomes, itol_colours$color, itol_colours$label, sep = "\t")
)

writeLines(output, "itol_phylum_colours.txt")

#Edit other parameters manually; in the generated script, replace spaces with \t





##########Differences between raw and sqrt

library(dplyr)

raw_sqrt_dif <- MAGG_abundance_sqrt_means %>%
  left_join(
    MAGG_abundance_means %>% select(-MAG_ID), 
    by = "Genomes") %>%
  mutate(
    sqrt_diff = sqrt_average_control - sqrt_average_stress,
    raw_diff = average_control - average_stress,
    .after = 4
  ) %>%
  relocate(raw_diff, .after = last_col())
View(raw_sqrt_dif)

write.table(raw_sqrt_dif, file = "raw_sqrt_dif.txt", sep = "\t", col.names = TRUE)
