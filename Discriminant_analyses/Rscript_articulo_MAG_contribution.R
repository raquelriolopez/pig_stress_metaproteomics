library(tidyverse)
library(KEGGREST)
library(ggplot2)
library(ggrepel)
library(ggtext)
install.packages("ggtext")

#MetaGenomics
#59 MAGs representation

mags59 <- PLSDA_MAGG_LOG_control_stress %>%
  rename(loadings = Loading_Value, Taxonomy = classification, MAG = Bin, Indicator = Association)

mags59$loadings <- mags59$loadings * -1

View(mags59)

mags59_clean <- mags59 %>%
  mutate(
    loadings = as.numeric(sub(",", ".", loadings)),
    last_taxon = sub(".*;(g__[^;]+).*", "\\1", Taxonomy),
    last_taxon = ifelse(grepl(";", last_taxon), 
                        sub(".*;(f__[^;]+).*", "\\1", Taxonomy),
                        last_taxon),
    last_taxon = ifelse(grepl(";", last_taxon), 
                        sub(".*;(o__[^;]+).*", "\\1", Taxonomy),
                        last_taxon),
    
    Species_label = ifelse(grepl(";s__[^;]+", Taxonomy),
                           sub(".*;s__([^;]+)", "\\1", Taxonomy),
                           paste0(last_taxon, " (NS)")),
    
    Species_short = ifelse(nchar(Species_label) > 40,
                           paste0(substr(Species_label, 1, 60), "..."),
                           Species_label),
    label_pos = ifelse(loadings > 0, 
                       loadings + 0.02 * diff(range(loadings)),
                       loadings - 0.02 * diff(range(loadings)))
  )

ggplot(mags59_clean, aes(x = reorder(MAG, loadings), y = loadings, fill = Indicator)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    aes(y = label_pos, label = Species_short,
        hjust = ifelse(loadings > 0, 0, 1)),
    size = 3, color = "black"
  ) +
  coord_flip() +
  scale_fill_manual(values = c("stress" = "red", "control" = "blue")) +
  labs(title = "MAG contribution to the PLS-DA", 
       x = "MAG", 
       y = "Loading") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    plot.caption = element_text(size = 8, color = "gray50")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.25)))


conteo_grupos <- mags55_clean %>%
  count(Indicator, name = "Total_MAGs") %>%
  mutate(Porcentaje = round(Total_MAGs / sum(Total_MAGs) * 100, 1))

print(conteo_grupos)

ggplot(conteo_grupos, aes(x = Indicator, y = Total_MAGs, fill = Indicator)) +
  geom_bar(stat = "identity") +
  geom_text(aes(la
                
head(significant_proteins_log)

negativos <- sum(significant_proteins_log$log2Fold_change < 0)

positivos <- sum(significant_proteins_log$log2Fold_change > 0)


cat("Número de proteínas con log2Fold_change negativo (control):", negativos, "\n")
cat("Número de proteínas con log2Fold_change positivo (stress):", positivos, "\n")
cat("Número de proteínas con log2Fold_change negativo (control):", negativos, "\n")
cat("Número de proteínas con log2Fold_change positivo (stress):", positivos, "\n")


# MetaProteomics

#39 MAGs representation
mags <- read.table("MAGP_LOG_association_control_stress_with_taxonomy.txt", header=T, sep="\t")
mags <- MAGP_LOG_association_control_stress_with_taxonomy
View(mags39)

mags39_clean <- mags39 %>%
  mutate(
    Loading_Value = as.numeric(sub(",", ".", Loading_Value)),
    last_taxon = sub(".*;(g__[^;]+).*", "\\1", classification),
    last_taxon = ifelse(grepl(";", last_taxon), 
                        sub(".*;(f__[^;]+).*", "\\1", classification),
                        last_taxon),
    last_taxon = ifelse(grepl(";", last_taxon), 
                        sub(".*;(o__[^;]+).*", "\\1", classification),
                        last_taxon),
    
    Species_label = ifelse(grepl(";s__[^;]+", classification),
                           sub(".*;s__([^;]+)", "\\1", classification),
                           paste0(last_taxon, " (NS)")),
    
    Species_short = ifelse(nchar(Species_label) > 50,
                           paste0(substr(Species_label, 1, 50), "..."),
                           Species_label),
    
    label_pos = ifelse(Loading_Value > 0, 
                       Loading_Value + 0.02 * diff(range(Loading_Value)),
                       Loading_Value - 0.02 * diff(range(Loading_Value))))

ggplot(mags39_clean, aes(x = reorder(Bin, Loading_Value), y = Loading_Value, fill = Association)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    aes(y = label_pos, label = Species_short,
        hjust = ifelse(Loading_Value > 0, 0, 0.95)),
    size = 3, color = "black"
  ) +
  coord_flip() +
  scale_fill_manual(values = c("stress" = "red", "control" = "blue")) +
  labs(title = "MAGP contribution to the PLS-DA", 
       x = "MAG", 
       y = "Loading") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    plot.margin = margin(1, 8, 1, 1, "cm"),
    plot.caption = element_text(size = 8, color = "gray50")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.3)))
