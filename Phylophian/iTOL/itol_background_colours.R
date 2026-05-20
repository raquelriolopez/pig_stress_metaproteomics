#itol magp tree background colours
# generate_tree_ranges_itol.R
# Genera TREE_COLORS (colored ranges) para iTOL, modo Label

input_file <- "labels.txt"
output_file <- "2. itol_phylum_background_Cover_Full"

# Labels
labels <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Norm cols
if (!("NODE_ID" %in% colnames(labels))) stop("Falta columna NODE_ID")
if (!("CLASS" %in% colnames(labels))) stop("Falta columna CLASS")

# Define colours
phylum_colors <- c(
  "Bacillota_A"        = "#a8ddb5", 
  "Bacillota_I"        = "#fc9272", 
  "Bacillota_B"        = "#9ebcda", 
  "Bacillota_C"        = "#fdae6b", 
  "Bacillota"          = "#c7e9c0", 
  "Bacteroidota"       = "#fff7bc", 
  "Spirochaetota"      = "#b3cde3",  
  "Methanobacteriota"  = "#d9d9d9", 
  "Elusimicrobiota"    = "#80b1d3", 
  "Pseudomonadota"     = "#fbb4ae", 
  "Myxococcota"        = "#decbe4", 
  "Verrucomicrobiota"  = "#ffffcc", 
  "Desulfobacterota"   = "#fccde5", 
  "Actinomycetota"     = "#ffd9b3",  
  "Fibrobacterota"     = "#d9f0a3",  
  "Planctomycetota"    = "#e0f3db", 
  "Patescibacteria"    = "#f0f0f0",  
  "Eremiobacterota"    = "#c5b0d5",  
  "Cyanobacteriota"    = "#D2DFFA",  
  "Thermoplasmatota"   = "#fff1a8"   
)


labels$COLOR <- phylum_colors[labels$CLASS]
labels$COLOR[is.na(labels$COLOR)] <- "#B0B0B0"

# TREE_COLORS
cat("TREE_COLORS\n", file = output_file)
cat("SEPARATOR TAB\n", file = output_file, append = TRUE)
cat("DATA\n", file = output_file, append = TRUE)
cat("#NODE_ID\tTYPE\tCOLOR\tLABEL\n", file = output_file, append = TRUE)

apply(labels, 1, function(x) {
  cat(paste0(x["NODE_ID"], "\trange\t", x["COLOR"], "\t", x["CLASS"], "\n"),
      file = output_file, append = TRUE)
})

cat(sprintf("✅ A'%s' file generated with %d entries (colored ranges).n", output_file, nrow(labels)))
