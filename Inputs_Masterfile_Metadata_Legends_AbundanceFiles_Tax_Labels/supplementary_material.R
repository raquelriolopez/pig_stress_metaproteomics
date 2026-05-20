library(tidyverse)
library("dplyr")
library(readxl)
library(writexl)

######Create excel with info for 42,696 proteins

#Inputs
MetaP2_raw <- read.delim("MetaP2_raw.txt", header = TRUE, stringsAsFactors = FALSE)
MetaP2_raw[is.na(MetaP2_raw)] <- 0 #NA to 0

MetaP2_raw_filtfreq <- read.delim("MetaP2_raw_filtfreq.txt", header = TRUE, stringsAsFactors = FALSE)

Leyenda_MetaP2 <- read.delim("Leyenda_MetaP2.txt", header = TRUE, stringsAsFactors = FALSE)

metadataMetaP <- read.delim("metadataMetaP.txt", header = TRUE, stringsAsFactors = FALSE)
X484_mag_taxonomy <- read_excel("484_mag_taxonomy.xlsx")
masterfile_genes_proteins_gutbrain <- read.delim("masterfile_genes_proteins_gutbrain.txt", header = TRUE, stringsAsFactors = FALSE)

#Calculate average raw abaundance per protein and per condition

# Samples
muestras_en_metadatos <- intersect(names(MetaP2_raw), metadataMetaP$SampleID)
View(muestras_en_metadatos)

# Samples per condition
controles <- metadataMetaP$SampleID[metadataMetaP$Tto == "control"]
estres <- metadataMetaP$SampleID[metadataMetaP$Tto == "stress"]

# Samples that exist in MetaP2_raw (just to confirm)
controles <- intersect(controles, muestras_en_metadatos)
estres <- intersect(estres, muestras_en_metadatos)

# Calculate average
info_proteins <- MetaP2_raw #rename to avoid loosing steps when making mistakes

info_proteins$Average_control <- rowMeans(info_proteins[, controles, drop = FALSE])
info_proteins$Average_stress <- rowMeans(info_proteins[, estres, drop = FALSE])

View(info_proteins)

# Add legend
nrow(Leyenda_MetaP2) == nrow(info_proteins)  # debe ser TRUE
info_proteins$gene <- Leyenda_MetaP2$Protein

#Add annotation hits
columnas_master <- c("gene", "fasta", "ko_id", "kegg_hit", "pfam_hits", "cazy_ids")

columnas_presentes <- intersect(columnas_master, names(masterfile_genes_proteins_gutbrain))
if (length(columnas_presentes) < length(columnas_master)) {
  warning("Faltan algunas columnas en masterfile: ", 
          paste(setdiff(columnas_master, columnas_presentes), collapse = ", "))
}

info_proteins <- merge(info_proteins, 
                       masterfile_genes_proteins_gutbrain[, columnas_presentes, drop = FALSE], 
                       by = "gene", 
                       all.x = TRUE)

if ("Protein" %in% names(Leyenda_MetaP2) && "Description" %in% names(Leyenda_MetaP2)) {
  info_proteins <- merge(info_proteins, 
                         Leyenda_MetaP2[, c("Protein", "Description")], 
                         by.x = "ene", 
                         by.y = "Protein", 
                         all.x = TRUE)
} else {
  warning("Leyenda_MetaP2 no tiene las columnas 'Protein' y/o 'Description'")
}


names(info_proteins)[names(info_proteins) == "fasta"] <- "genome"

#Add taxonomy

info_proteins <- merge(info_proteins, 
                       X484_mag_taxonomy[, c("user_genome", "classification")], 
                       by.x = "genome", 
                       by.y = "user_genome", 
                       all.x = TRUE)

View(info_proteins)
colnames(info_proteins)

# Reorganise
order <- c("gene", "ko_id", "kegg_hit", "pfam_hits", "cazy_ids", 
           "Description", "genome", "classification", 
           "Average_control", "Average_stress")

# Reordenar el dataframe
info_proteins <- info_proteins[, c(order, setdiff(names(info_proteins), order))]

# Create codes for proteins
info_proteins$Protein <- paste0("P", 1:nrow(info_proteins))

info_proteins <- info_proteins[, c("Protein", setdiff(names(info_proteins), "Protein"))]

#Create column to specify if they passed the frequency filter

ids_filtfreq <- rownames(MetaP2_raw_filtfreq)
info_proteins$Kept <- ifelse(info_proteins$Protein %in% ids_filtfreq, "YES", "NO")
info_proteins <- info_proteins[, c("Protein", "Kept", setdiff(names(info_proteins), c("Protein", "Kept")))]


length(unique(info_proteins$genome)) # 483/484 from MAG catalogue, MEGAHIT-MetaBAT2-S959.102, CAG-533 sp004563275, had no proteins identified 

#save excel
write_xlsx(info_proteins, "info_proteins.xlsx")


############### Create excel with info for 

taxonomia_MetaP2_1<- read_excel("taxonomia_MetaP2.xlsx", sheet = 1)
taxonomia_MetaP2_4<- read_excel("taxonomia_MetaP2.xlsx", sheet = 4)
View(taxonomia_MetaP2_4)
View(taxonomia_MetaP2_1)
MAGG_raw <- read_tsv("MAGG_abundance_raw.tsv", col_names = TRUE)
MAGP_raw<- read.delim("MAGP_abundance_raw_NeW.txt", header = TRUE)

View(MAGG_raw)
view(MAGP_raw)
#Merge both dataframes and select columns

info_genomes_raw <- merge(taxonomia_MetaP2_1, taxonomia_MetaP2_4, by.x = "user_genome", by.y = "bins", all = TRUE)

info_genomes_raw <- info_genomes_raw[, c("user_genome", 
                                                  "Completeness", 
                                                  "Contamination", 
                                                  "Genome_Size", 
                                                  "classification", 
                                                  "closest_genome_ani", 
                                                  "closest_genome_reference", 
                                                  "closest_placement_ani", 
                                                  "closest_placement_reference", 
                                                  "classification_method")]

#Calculate average MAGG and MAGP for control and stress condition.

control_samples <- metadataMetaP$SampleID[metadataMetaP$Tto == "control"]  #The same as we did before
stress_samples <- metadataMetaP$SampleID[metadataMetaP$Tto == "stress"]

# Average MAGG_raw
MAGG_raw$Average_MAGG_control <- rowMeans(MAGG_raw[, control_samples, drop = FALSE])
MAGG_raw$Average_MAGG_stress <- rowMeans(MAGG_raw[, stress_samples, drop = FALSE])

# Average MAGP_raw
MAGP_raw$Average_MAGP_control <- rowMeans(MAGP_raw[, control_samples, drop = FALSE])
MAGP_raw$Average_MAGP_stress <- rowMeans(MAGP_raw[, stress_samples, drop = FALSE])

#### Merge averages and raw abundances

#to avoid missunderstanting we add MP or MG to sample columns from MAGP abundance file or MAGG abundance file
muestras_MAGP <- grep("^X", names(MAGP_raw), value = TRUE)
names(MAGP_raw)[names(MAGP_raw) %in% muestras_MAGP] <- paste0("MP", muestras_MAGP)

muestras_MAGG <- grep("^X", names(MAGG_raw), value = TRUE)
names(MAGG_raw)[names(MAGG_raw) %in% muestras_MAGG] <- paste0("MG", muestras_MAGG)

#Merge
info_temp <- merge(info_genomes_raw, MAGG_raw, 
                   by.x = "user_genome", by.y = "Genomes", 
                   all.x = TRUE)  # left join

info_genomes <- merge(info_temp, MAGP_raw, 
                    by.x = "user_genome", by.y = "Bin", 
                    all.x = TRUE)

# Reorganise
abundancias <- grep("^(MP|MG)", names(info_genomes), value = TRUE)
nombres_sin_abund <- setdiff(names(info_genomes), abundancias)
pos_size <- which(nombres_sin_abund == "Genome_Size")

nuevo_orden <- c(
  nombres_sin_abund[1:pos_size],          
  abundancias,                             
  nombres_sin_abund[(pos_size + 1):length(nombres_sin_abund)]  
)

info_genomes <- info_genomes[, nuevo_orden]

# Average cols
avg_cols <- c("Average_MAGG_control", "Average_MAGG_stress", "Average_MAGP_control", "Average_MAGP_stress")

info_genomes <- info_genomes[, c(
  names(info_genomes)[1:which(names(info_genomes) == "Genome_Size")],
  avg_cols,
  setdiff(names(info_genomes), c(names(info_genomes)[1:which(names(info_genomes) == "Genome_Size")], avg_cols))
)]


View(info_genomes)

#save
write_xlsx(info_genomes, "info_genomes_raw.xlsx")

#############Create table for genes_info

gene_cols <- c("gene", "gene_position", "start_position", "end_position", 
                       "strandedness", "rank", "ko_id", "kegg_hit", "pfam_hits", 
                       "cazy_ids", "cazy_hits")

genes_info <- masterfile_genes_proteins_gutbrain[, gene_cols]

View(genes_info)

#save
write_xlsx(genes_info, "info_genes.xlsx")


##########################CREATE one single excel for all supplementary tables

#Inputs manually
info_genomes_raw <- read_excel("info_genomes_raw.xlsx")
info_genes <- read_excel("info_genes.xlsx")
info_proteins_raw <- read_excel("info_proteins_raw.xlsx")
ST3_Volcano <- read_excel("ST3_Volcano_ MAG_and_MAG_protein_abundance_significant_mags.xlsx")
ST6_Enrichment <- read_excel("ST6_Enrichment_from_discriminant_proteins.xlsx")
ST9_Comparison <- read_excel("ST9_Comparison_two_approaches_functions.xlsx")

ST1_PLSDA_MAG <- read.delim("ST1_PLSDA_MAG_abundance_LOG_discriminant_mags.txt", 
                            header = TRUE, stringsAsFactors = FALSE, sep = "\t")
ST2_PLSDA_MAG_protein <- read.delim("ST2_PLSDA_MAG_protein_abundance_LOG_discriminant_mags.txt", 
                                    header = TRUE, stringsAsFactors = FALSE, sep = "\t")
ST4_PLSDA_proteins <- read.delim("ST4_PLSDA_LOG_discriminant_individual_proteins.txt", 
                                 header = TRUE, stringsAsFactors = FALSE, sep = "\t")
ST5_Volcano_proteins <- read.delim("ST5_Volcano_LOG_significant_individual_proteins.txt", 
                                   header = TRUE, stringsAsFactors = FALSE, sep = "\t")
ST7_PLSDA_functions <- read.delim("ST7_PLSDA_LFQ_KO_ID_LOG_discriminant_functions.txt", 
                                  header = TRUE, stringsAsFactors = FALSE, sep = "\t")
ST8_Volcano_functions <- read.delim("ST8_Volcano_LFQ_KO_ID_LOG_significant_functions.txt", 
                                    header = TRUE, stringsAsFactors = FALSE, sep = "\t")



#list
archivos_y_hojas <- list(
  c("info_genomes_raw.xlsx", "Table S1"),
  c("info_genes.xlsx", "Table S2"),
  c("info_proteins_raw.xlsx", "Table S3"),
  c("ST1_PLSDA_MAG_abundance_LOG_discriminant_mags.txt", "Table S4"),
  c("ST2_PLSDA_MAG_protein_abundance_LOG_discriminant_mags.txt", "Table S5"),
  c("ST3_Volcano_ MAG_and_MAG_protein_abundance_significant_mags.xlsx", "Table S6"),
  c("ST4_PLSDA_LOG_discriminant_individual_proteins.txt", "Table S7"),
  c("ST5_Volcano_LOG_significant_individual_proteins.txt", "Table S8"),
  c("ST6_Enrichment_from_discriminant_proteins.xlsx", "Table S9"),
  c("ST7_PLSDA_LFQ_KO_ID_LOG_discriminant_functions.txt", "Table S10"),
  c("ST8_Volcano_LFQ_KO_ID_LOG_significant_functions.txt", "Table S11"),
  c("ST9_Comparison_two_approaches_functions.xlsx", "Table S12")
)

# Lista para almacenar los dataframes
lista_dataframes <- list()

# Bucle para leer cada archivo
for (i in seq_along(archivos_y_hojas)) {
  archivo <- archivos_y_hojas[[i]][1]
  nombre_hoja <- archivos_y_hojas[[i]][2]
  
  cat("Procesando archivo:", archivo, "como", nombre_hoja, "\n")
  
  # Comprobar que el archivo existe
  if (!file.exists(archivo)) {
    stop("El archivo ", archivo, " no existe en el directorio de trabajo: ", getwd())
  }
  
  # Leer según la extensión
  if (grepl("\\.xlsx$", archivo, ignore.case = TRUE)) {
    # Archivo Excel
    df <- read_excel(archivo)
  } else if (grepl("\\.txt$", archivo, ignore.case = TRUE)) {
    # Archivo de texto (asumimos tabulado)
    df <- read.delim(archivo, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  } else {
    stop("Formato de archivo no soportado: ", archivo)
  }
  
  # Añadir a la lista con el nombre de hoja deseado
  lista_dataframes[[nombre_hoja]] <- df
}

# Escribir el archivo Excel con todas las hojas
write_xlsx(lista_dataframes, "Supplementary_Tables.xlsx")

cat("✅ Archivo Excel 'Supplementary_Tables.xlsx' creado con", length(lista_dataframes), "hojas.\n")


###Add MAG labels (from phylogenetic tree)

labels_tax <- read.delim("labels_tax.txt", header = TRUE, stringsAsFactors = FALSE)
View(labels_tax)

#merge labels
info_genomes <- info_genomes %>%
  left_join(select(labels_tax, NODE_ID, LABEL), by = c("user_genome" = "NODE_ID")) %>%
  relocate(LABEL, .after = user_genome)

write_xlsx(info_genomes, "info_genomes.xlsx")
