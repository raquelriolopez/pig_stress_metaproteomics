########Libraries
library(microeco)
library(magrittr)
library(tidyverse)
library(readxl)

########INPUT
MAGP_abundance <- read.table("protein_sums_per_bin_gutbrain.txt", header = TRUE, sep = "\t")
MAGG_abundance <- read.table("MAGG_abundance.tsv", header = TRUE, sep = "\t")

MAGP_abundance <- MAGP_abundance %>%  ######Transform LOG data into raw
  mutate(across(where(is.numeric), ~ 2^.x)) %>%
  rename(., "Genomes" = "Bin")

metadata <- read.table("metadataMetaP.txt", header = TRUE, sep = "\t") %>%
  rename(treatment = Tto) %>%    
  rename(pen = corral)
rownames(metadata) <- metadata$SampleID

taxonomy <- read_excel("C:/Users/RRIO/OneDrive - IRTA/Escritorio/DOCTORAT/SEGUNDO AÑO/Estancia/DATA for the stay/protein_legend, metadata, all_proteins abundance and MAG taxonomy/taxonomia_MetaP2.xlsx")

######### Prepare the data 
# Samples as rows, MAGs as columns
abundance_table <- MAGP_abundance %>% #change for MAGG_abundance or MAGP_abundance if needed
  column_to_rownames("Genomes") %>%
  as.data.frame()

# Format taxonomy
taxonomy_table <- taxonomy %>% 
  select(user_genome, classification) %>% 
  separate(classification, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% 
  column_to_rownames("user_genome")

# Microtable object
dataset <- microtable$new(
  sample_table = metadata,
  otu_table = abundance_table,
  tax_table = taxonomy_table
)

# Beta plots
# Recalculate beta diversity (Bray-Curtis, Jaccard, UniFrac)
dataset$cal_betadiv(unifrac = FALSE)
dataset$beta_diversity$bray

t_beta <- trans_beta$new(dataset = dataset, group = "treatment", measure = "bray")
t_beta$cal_ordination(method = "PCoA")

# Plot PCoA
t_beta$plot_ordination(plot_color = "treatment", plot_shape = "sex")

#########################################################################RESULTS - MAGG

####permANOVA based on adonis2, calculate r2 and p-value related to beta-diversity based on experimental groups

t1 <- trans_beta$new(dataset = dataset, group = "treatment", measure = "bray")

# manova for all groups when manova_all = TRUE

t1$cal_manova(manova_all = TRUE)
t1$res_manova
# > t1$res_manova
#           Df   SumOfSqs         R2         F Pr(>F) Significance
# treatment  1 0.08058831 0.01452311 0.8547543  0.553             
# Residual  58 5.46838095 0.98547689        NA     NA         <NA>
# Total     59 5.54896926 1.00000000        NA     NA         <NA>

  
t1$cal_manova(manova_all = FALSE)
t1$res_manova
# > t1$res_manova
#              Groups measure         F         R2 p.value p.adjusted Significance
# 1 control vs stress    bray 0.8547543 0.01452311   0.527      0.527  

t1$cal_manova(manova_all = FALSE, group = "treatment", by_group = "sex")
t1$res_manova
# > t1$res_manova
#   by_group            Groups measure         F         R2 p.value p.adjusted Significance
# 1        M control vs stress    bray 0.6924873 0.02413479   0.741      0.741             
# 2        F control vs stress    bray 0.6619519 0.02309514   0.697      0.697 

t1$cal_manova(manova_all = FALSE, group = "sex", by_group = "treatment")
t1$res_manova
# > t1$res_manova
#   by_group Groups measure         F         R2 p.value p.adjusted Significance
# 1  control M vs F    bray 0.9664021 0.03224952   0.436      0.436             
# 2   stress F vs M    bray 1.3720946 0.04836071   0.178      0.178  

t1$cal_manova(manova_all = FALSE, group = "pen", by_group = "treatment")
t1$res_manova
# > t1$res_manova
#    by_group    Groups measure         F         R2 p.value p.adjusted Significance
# 1   control  C5 vs C3    bray 3.3416072 0.19269305   0.006     0.0210            *
# 2   control  C5 vs C9    bray 2.3652646 0.15393582   0.033     0.0660             
# 3   control  C5 vs C7    bray 2.8911538 0.17116379   0.007     0.0210            *
# 4   control  C3 vs C9    bray 0.7789526 0.05653206   0.558     0.6696             
# 5   control  C3 vs C7    bray 0.6140153 0.04201551   0.827     0.8270             
# 6   control  C9 vs C7    bray 1.0512986 0.07481861   0.327     0.4905             
# 7    stress  C6 vs C4    bray 0.9467649 0.07312753   0.496     0.4960             
# 8    stress C6 vs C10    bray 3.0923262 0.20489394   0.005     0.0300            *
# 9    stress  C6 vs C8    bray 3.1902529 0.19704775   0.013     0.0390            *
# 10   stress C4 vs C10    bray 1.6421421 0.12037275   0.113     0.2260             
# 11   stress  C4 vs C8    bray 1.3617562 0.09481823   0.219     0.3285             
# 12   stress C10 vs C8    bray 1.0137277 0.07233819   0.398     0.4776             

t1$cal_manova(manova_set = "treatment + sex")
t1$res_manova
# > t1$res_manova
#           Df   SumOfSqs         R2         F Pr(>F) Significance
# treatment  1 0.08058831 0.01452311 0.8673613  0.512             
# sex        1 0.17239485 0.03106790 1.8554630  0.068             
# Residual  57 5.29598610 0.95440898        NA     NA         <NA>
# Total     59 5.54896926 1.00000000        NA     NA         <NA>



#########################################################################RESULTS - MAGP

####permANOVA based on adonis2, calculate r2 and p-value related to beta-diversity based on experimental groups

t1 <- trans_beta$new(dataset = dataset, group = "treatment", measure = "bray")

# manova for all groups when manova_all = TRUE

t1$cal_manova(manova_all = TRUE)
t1$res_manova
# > t1$res_manova
#           Df   SumOfSqs         R2         F Pr(>F) Significance
# treatment  1 0.04457293 0.01406738 0.8275492  0.581             
# Residual  58 3.12395882 0.98593262        NA     NA         <NA>
# Total     59 3.16853175 1.00000000        NA     NA         <NA>


t1$cal_manova(manova_all = FALSE)
t1$res_manova
# > t1$res_manova
#              Groups measure         F         R2 p.value p.adjusted Significance
# 1 control vs stress    bray 0.8275492 0.01406738   0.562      0.562  

t1$cal_manova(manova_all = FALSE, group = "treatment", by_group = "sex")
t1$res_manova
# > t1$res_manova
# by_group            Groups measure         F         R2 p.value p.adjusted Significance
# 1        M control vs stress    bray 0.7637882 0.02655381   0.689      0.689             
# 2        F control vs stress    bray 0.6129520 0.02142219   0.783      0.783

t1$cal_manova(manova_all = FALSE, group = "sex", by_group = "treatment")
t1$res_manova
# > t1$res_manova
#   by_group Groups measure        F         R2 p.value p.adjusted Significance
# 1  control M vs F    bray 0.597263 0.02017967   0.742      0.742             
# 2   stress F vs M    bray 1.368515 0.04824064   0.171      0.171 

t1$cal_manova(manova_all = FALSE, group = "pen", by_group = "treatment")
t1$res_manova
# > t1$res_manova
#    by_group    Groups measure        F         R2 p.value p.adjusted Significance
# 1   control  C5 vs C3    bray 3.188775 0.18551497   0.011     0.0330            *
# 2   control  C5 vs C9    bray 2.386688 0.15511380   0.056     0.0840             
# 3   control  C5 vs C7    bray 2.347294 0.14358915   0.008     0.0330            *
# 4   control  C3 vs C9    bray 1.980282 0.13219254   0.092     0.1104             
# 5   control  C3 vs C7    bray 1.033352 0.06873728   0.389     0.3890             
# 6   control  C9 vs C7    bray 2.749002 0.17455086   0.028     0.0560             
# 7    stress  C6 vs C4    bray 0.923046 0.07142635   0.486     0.4860             
# 8    stress C6 vs C10    bray 2.165415 0.15286634   0.052     0.1380             
# 9    stress  C6 vs C8    bray 2.458422 0.15903450   0.040     0.1380             
# 10   stress C4 vs C10    bray 1.539729 0.11371932   0.138     0.1656             
# 11   stress  C4 vs C8    bray 1.706801 0.11605524   0.069     0.1380             
# 12   stress C10 vs C8    bray 1.465140 0.10128767   0.135     0.1656             

t1$cal_manova(manova_set = "treatment + sex")
t1$res_manova
# > t1$res_manova
#           Df   SumOfSqs         R2         F Pr(>F) Significance
# treatment  1 0.04457293 0.01406738 0.8310767  0.565             
# sex        1 0.06689246 0.02111150 1.2472317  0.260             
# Residual  57 3.05706637 0.96482112        NA     NA         <NA>
# Total     59 3.16853175 1.00000000        NA     NA         <NA>


