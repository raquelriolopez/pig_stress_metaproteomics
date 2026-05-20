library(dplyr)

# Crear dataframes con las columnas deseadas
MAGP_control <- MAGP_abundance[, c("MAG", "Control")]
MAGP_stress  <- MAGP_abundance[, c("MAG", "Stress")]

# Guardar como archivos .txt (tabulado, sin row.names)
write.table(MAGP_control, "MAGP_control.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(MAGP_stress,  "MAGP_stress.txt",  sep = "\t", row.names = FALSE, quote = FALSE)


taxonomy_MAG_C70C10_MetaP2 <- read_excel("C:/Users/RRIO/OneDrive - IRTA/Escritorio/DOCTORAT/SEGUNDO AÑO/Estancia/DATA for the stay/protein_legend, metadata, all_proteins abundance and MAG taxonomy/taxonomy_MAG_C70C10_MetaP2.xlsx")


labels_tax <- labels %>%
  left_join(
    taxonomy_MAG_C70C10_MetaP2 %>% select(user_genome, classification),
    by = c("NODE_ID" = "user_genome")
  )

View(labels_tax)

write.table(labels_tax,  "labels_tax.txt",  sep = "\t", row.names = FALSE, quote = FALSE)
