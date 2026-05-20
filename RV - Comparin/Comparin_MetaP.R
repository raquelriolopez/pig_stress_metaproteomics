
library("LinkHD")
library(reshape2)

MAGP_clr= read.table("MAGP_clr.txt", header=T, row.names=1, check.names=F)
MAGP_log = read.table("MAGP_log.txt", header=T, row.names=1, check.names=F)
MAG70_raw = read.table("MAG70_raw.txt", header=T, row.names=1, check.names=F)

#check rownames
identical(rownames(MAGP_clr), rownames(MAGP_log)) && identical(rownames(MAGP_log), rownames(MAG70_raw))

Datasets<- list(MAGP_clr, MAGP_log, MAG70_raw)


Normalization<-lapply(list(Datasets[[3]],Datasets[[3]]),function(x){DataProcessing(x,Method="Compositional")})

testing = list(Normalization[[1]], Datasets[[1]], Datasets[[2]])
names(testing)<-c("MAGG_clr","MAGP_clr", "MAGP_log")


Output<-LinkData(testing,Distance=rep("euclidean",3),Scale = FALSE,Center=FALSE)

library(corrplot)
corrplot(Output@RV, method = 'number', type = 'upper', insig='blank', tl.col="black", addCoef.col ='black', number.cex = 0.8, order = 'original', diag=FALSE, number.digits=3)

