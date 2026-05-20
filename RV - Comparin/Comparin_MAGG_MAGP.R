devtools::install_github(repo="lauzingaretti/LinkHD")
library("LinkHD")
library(reshape2)
library("easyCODA")

#Input
MAGP_clr= read.table("MAGP_abundance_clr_New.txt", header=T, row.names=1, check.names=F)
MAGP_log = read.table("MAGP_abundance_log_New.txt", header=T, row.names=1, check.names=F)
MAGG_raw = read.table("MAGG_abundance.tsv", header=T, row.names=1, check.names=F)

#Adapt the format
MAGP_log_t <- as.data.frame(t(MAGP_log[-1]))
colnames(MAGP_log_t) <- MAGP_log$Bin  
MAGP_log <- MAGP_log_t

MAGP_clr_c <- as.data.frame(t(MAGP_clr[-1]))  
colnames(MAGP_clr_c) <- MAGP_clr$Bin
MAGP_clr <- MAGP_clr_c

MAGG_raw <- as.data.frame(t(MAGG_raw))
MAGG_log <- log2(MAGG_raw + 1)
MAGG_clr <- CLR(MAGG_raw + 0.0001, weight = TRUE)$LR

View(MAGP_clr)
View(MAGP_log)
View(MAGG_raw)
View(MAGG_log)
View(MAGG_clr)


#comprobar rownames
identical(rownames(MAGP_clr), rownames(MAGP_log)) && identical(rownames(MAGP_clr), rownames(MAGG_clr))

Datasets<- list(MAGP_clr, MAGP_log, MAGG_log, MAGG_clr)

testing = Datasets
names(testing)<-c("MAGP_clr", "MAGP_log", "MAGG_log", "MAGG_clr")


Output<-LinkData(testing,Distance=rep("euclidean", 4),Scale = FALSE,Center=FALSE) ###Me da error aquí

library(corrplot)
corrplot(Output@RV, method = 'number', type = 'upper', insig='blank', tl.col="black", addCoef.col ='black', number.cex = 0.8, order = 'original', diag=FALSE, number.digits=3)



############ Other way


verify_structure <- function(...) {
  datasets <- list(...)
  cat("\n=== Verificación de Estructura ===\n")
  cat("Número de datasets:", length(datasets), "\n")
  cat("Muestras comunes:", length(Reduce(intersect, lapply(datasets, rownames))), "\n\n")
  
  for(i in seq_along(datasets)) {
    cat("Dataset", i, ":\n")
    cat(" - Dimensiones:", dim(datasets[[i]]), "\n")
    cat(" - Clase:", class(datasets[[i]]), "\n")
    cat(" - NA values:", sum(is.na(datasets[[i]])), "\n")
    cat(" - Inf values:", sum(sapply(datasets[[i]], is.infinite)), "\n\n")
  }
}

verify_structure(MAGP_clr, MAGP_log, MAGG_log, MAGG_clr)


common_samples <- Reduce(intersect, list(rownames(MAGP_clr), 
                                         rownames(MAGP_log),
                                         rownames(MAGG_log),
                                         rownames(MAGG_clr)))

testing <- list(
  MAGP_clr = MAGP_clr[common_samples, ],
  MAGP_log = MAGP_log[common_samples, ],
  MAGG_log = MAGG_log[common_samples, ],
  MAGG_clr = MAGG_clr[common_samples, ]
)

try_LinkData <- function(data_list) {
  data_list <- lapply(data_list, function(x) as.matrix(as.data.frame(x)))
  data_list <- lapply(data_list, function(x) {
    x[is.na(x) | is.infinite(x)] <- 0
    return(x)
  })
  tryCatch({
    LinkData(data_list, 
             Distance = rep("euclidean", length(data_list)),
             Scale = FALSE,
             Center = FALSE)
  }, error = function(e) {
    message("Falló euclidean, intentando con distancia 'bray'")
    LinkData(data_list, 
             Distance = rep("bray", length(data_list)),
             Scale = FALSE,
             Center = FALSE)
  })
}

Output <- try_LinkData(testing)

if(exists("Output")) {
  corrplot(Output@RV, 
           method = 'number', 
           type = 'upper', 
           insig='blank',
           tl.col="black",
           addCoef.col= "black",
          number.cex=0.8,
          order="original",
          diag=FALSE,
          number.digits=3,
          mar=c(0,0,1,0),
          title="Matriz RV de Métodos de Normalización")} 
else {
  stop("No se pudo calcular la matriz RV. Revisa los datos de entrada.")
  }
