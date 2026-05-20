# Decoding stress responses in the pig gut microbiome: a metaproteomic perspective for biomarker discovery

## Authors
Raquel Río-López \1, Judith Guitart-Matas \2, Magnus Ø. Arntzen \3, Adrià Clavell-Sansalvador \4, Ioanna-Theoni Vourlaki \4, Olga González-Rodríguez \4, L. Jesús García-Gil \5, Xavier Xifró \6, Yuliaxis Ramayo-Caldas\4\* and Antoni Dalmau\1\* 

\1 Animal Welfare Program, Institute of Agrifood Research and Technology (IRTA), 17121, Monells, Girona, Spain.
\2 Joint Research Unit IRTA-UAB in Animal Health, Animal Health Research Centre (CReSA), Autonomous University of Barcelona (UAB), Catalonia, Spain.
\3 Faculty of Chemistry, Biotechnology and Food Science, Norwegian University of Life Sciences (NMBU), 1433, Ås, Norway
\4 Animal Breeding and Genetics Program, Institute of Agrifood Research and Technology (IRTA), 08140, Caldes de Montbui, Barcelona, Spain.
\5 Digestive Diseases and Microbiota Group, Biomedical Research Institute of Girona (IDIBGI), 17190, Girona, Girona, Spain.
\6 New Therapeutic Targets Lab Research Group, Medical Sciences Department, Faculty of Medicine, University of Girona, 17071, Girona, Girona, Spain


\* Corresponding authors: yuliaxis.ramayo@irta.cat · antoni.dalmau@irta.cat

## Overview

This repository contains the bioinformatic workflows and analysis scripts associated with a genome-resolved metaproteomic study of the pig gut microbiome under chronic social stress. The study profiled 60 pig faecal samples from two conditions: a stress group (n=30; reduced space allowance + social mixing) and a control group (n=30; standard conditions). By combining shotgun metagenomics and nano LC-MS/MS metaproteomics, we linked expressed proteins directly to metagenome-assembled genomes (MAGs), enabling the identification of stress-responsive microbial taxa and functions with high discriminatory accuracy (99.51% via PLS-DA).

## Bioinformatic Workflow

### Metagenomics

- **Quality control & host decontamination**: nf-core/mag v3.3.0 (reads mapped 
  against *Sus scrofa* Sscrofa11.1 genome)
- **Assembly**: MEGAHIT v1.0.2 (individual and co-assemblies per pen)
- **Binning**: SemiBin2 + MetaBAT2, refined with DAS Tool
- **Dereplication**: dRep v3.5.0 (≥70% completeness, ≤10% contamination)
- **MAG quality**: CheckM2
- **Taxonomy**: GTDB-Tk v2.4.1
- **Functional annotation**: DRAM (PFAM-A, KOfam, dbCAN-V10)
- **Abundance**: CoverM
- **Metabolic potential**: MetQy v1.1.0 (KEGG module completion fractions)
- **Phylogenetic tree**: Phylophlan v3.0.60 + MAFFT v7.505 + RAxML v8.2.12 + 
  FastTree v2.1.11 (visualised in IToL v7.2.1)
  
### Metaproteomics

- **Database search**: FragPipe v19.0 + MSFragger (against MAG-derived protein 
  database, 911,040 entries; 1% FDR)
- **Quantification**: IonQuant (LFQ + match between runs)
- **Statistical analysis**: Perseus v2.1.3.0
- **Multivariate analysis**: PLS-DA via mixOmics v6.30.0 (R)
- **Differential abundance**: Student's t-test + volcano plots (ggplot2 v3.5.1)
- **Functional enrichment**: Hypergeometric test (phyper) + KEGG Mapper
- **Multi-omics integration**: Link-HD (RV coefficient)
- **Diversity**: Shannon index + Bray-Curtis / PerMANOVA (vegan v2.6.10)

## Contact
raquel.rio@irta.cat · yuliaxis.ramayo@irta.cat · antoni.dalmau@irta.cat