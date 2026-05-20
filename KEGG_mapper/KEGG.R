library(dplyr)
library(readr)
library(tidyr)
library(writexl)
library(devtools)

#Input

masterfile_genes_proteins_gutbrain <- read.delim("KEGG/masterfile_genes_proteins_gutbrain.txt")
View(masterfile_genes_proteins_gutbrain)
#Extract the KEGG data per MAG and past it into KEGG Mapper Reconstruct.

get_KEGG_mapper_MAG_clipboard <- function(mag) {
  t <- masterfile_genes_proteins_gutbrain %>%
    select(gene, Bin, ko_id) %>%
    filter(Bin == mag, !is.na(ko_id)) %>%
    select(-Bin)
  
  print(paste("Copied to clipboard. Protein count:", nrow(t)))
  write.table(t, paste0("clipboard-", 2^15), sep = "\t", row.names = FALSE, quote = FALSE) # ~32Mb limit (should be sufficient for one MAG)
}

#Then you write: get_KEGG_mapper_MAG_clipboard("the_name_of_your_mag_of_interest")

#Then you go to "KEGG Mapper reconstruct", paste the list and clic "Exec"
#https://www.kegg.jp/kegg/mapper/reconstruct.html


get_KEGG_mapper_MAG_clipboard("935_SemiBin_192") #s__Treponema_D sp018385415 (stress MAGG_PLSDA_clr; MAGG_volcano)
get_KEGG_mapper_MAG_clipboard("group6_bins.106") # g__RGIG464;s__ (stress MAGG_PLSDA_clr; MAGP_PLSDA_log; MAGP_volcano)
get_KEGG_mapper_MAG_clipboard("group9_bins.244") ############# s__Treponema_D sp016292805 (group9_bins.244) MAGG and MAGP
get_KEGG_mapper_MAG_clipboard("992_SemiBin_233") # Treponema_D sp002395155 Proteins, functions, KO_ID
get_KEGG_mapper_MAG_clipboard("MEGAHIT-MetaBAT2-S947.49") # Treponema_D sp018385315 MAGP
get_KEGG_mapper_MAG_clipboard("MEGAHIT-MetaBAT2-S948.142") #Treponema_D sp016293915

get_KEGG_mapper_MAG_clipboard("959_SemiBin_129") #g__Ruminococcus_E;s__
get_KEGG_mapper_MAG_clipboard("MEGAHIT-MetaBAT2-S938.43") # Ruminococcus E bovis
get_KEGG_mapper_MAG_clipboard("group7_bins.369") #Ruminococcaeae RUG12519
get_KEGG_mapper_MAG_clipboard("MEGAHIT-MetaBAT2-S972.52") ###############PeH17 sp004556165 (MAGP_PLSDA_log; MAGP_volcano)


head(final_table)
