#PLS-DA
library(mixOmics)
library("factoextra")
library("caTools")

#INPUTs
MAGP_abundance_raw <- read.table("MAGP_abundance_raw_New.txt", header=T, sep="\t") ##RAW MAG abundance from proteins
MAGP_abundance_raw <- as.data.frame(t(MAGP_abundance_raw))

MAGP_abundance_log <- read.table("MAGP_abundance_log_New.txt", header=T, sep="\t") ##LOG MAG abundance from proteins
MAGP_abundance_log <- as.data.frame(t(MAGP_abundance_log))

MAGP_abundance_clr <- read.table("MAGP_abundance_clr_New.txt", header=T, sep="\t") ##CLR MAG abundance from proteins
MAGP_abundance_clr <- as.data.frame(t(MAGP_abundance_clr))

View(MAGP_abundance_raw)
View(mag_abund_numeric)
#Transform only numeric values and keep bin_names RAW
bin_names <- MAGP_abundance_raw[1, ]

mag_abund_numeric <- MAGP_abundance_raw[-1, ]
row_names_num <- rownames(mag_abund_numeric)
mag_abund_numeric <- as.data.frame(lapply(mag_abund_numeric, as.numeric))
rownames(mag_abund_numeric) <- row_names_num
View(mag_abund_numeric)
MAGP_raw <- mag_abund_numeric

#Transform only numeric values
mag_abund_numeric_l <- MAGP_abundance_log[-1, ]
mag_abund_numeric_l <- as.data.frame(lapply(mag_abund_numeric_l, as.numeric))
rownames(mag_abund_numeric_l) <- rownames(mag_abund_numeric)
colnames(mag_abund_numeric_l) <- colnames(mag_abund_numeric)
View(mag_abund_numeric_l)
log_transformed <- mag_abund_numeric_l

#Transform only numeric values
mag_abund_numeric_c <- MAGP_abundance_clr[-1, ]
mag_abund_numeric_c <- as.data.frame(lapply(mag_abund_numeric_c, as.numeric))
rownames(mag_abund_numeric_c) <- rownames(mag_abund_numeric)
colnames(mag_abund_numeric_c) <- colnames(mag_abund_numeric)
View(mag_abund_numeric_c)
clr_transformed <- mag_abund_numeric_c


#Specify the input here (clr_transformed or log_transformed)
input_PLSDA <- MAGP_raw
input_PLSDA <- log_transformed
input_PLSDA <- clr_transformed

View(input_PLSDA)

#Rebuild with bin_names
input_w_names<- rbind(Bin = bin_names, input_PLSDA)

#write.table(input_w_names,"input_w_names.txt",row.names = T,col.names = T,quote = F,sep="\t")

#Metadata
metadata = read.table("metadataMetaP.txt", header=T, row.names=1, dec=",")
metadata = metadata[rownames(input_w_names)[-1], ]
Tto = as.factor(metadata$Tto)

#PLSDA
plsda <- plsda(input_PLSDA, Tto,ncomp=10,scale=F)

set.seed(30)
perf.pls <- perf(plsda, validation = "Mfold",criterion="all",folds = 4, progressBar = F,nrepeat = 100,)

#Optimal number of components

comp <- perf.pls[["choice.ncomp"]][6]

#BER for the number of components selected

err.n <- perf.pls[["error.rate"]]$BER[20+comp]
sd.n <- perf.pls[["error.rate.sd"]]$BER[20+comp]
err.total <- err.n + sd.n


#Variable important prediction (VIP). Contribution of each variable in the classification among lines/populations
vip <- data.frame(vip(plsda),stringsAsFactors = FALSE)

#A VIP higher than 1 was used as a threshold for selecting the variables with the highest contribution in the model

p <- 1
v.select <- vip[vip[,comp]>=p,]
v.ID <- row.names(v.select)

#Iterative process to obtain the best model with the lower BER

err <- 1
comp.cte <- 10

while (err.total<err) {
  
  err <- err.total
  sd.f <- sd.n
  comp.f <- comp
  
  vf<-v.ID
  vip <- data.frame(vip(plsda),stringsAsFactors = FALSE)
  v.select <- vip[vip[,comp] >= p,]
  
  v.ID <- row.names(v.select)
  
  if (length(v.ID) == 1) {
      break
  
    } 
  
  filter <- input_PLSDA[,colnames(input_PLSDA) %in% v.ID]
  
   if (length(v.ID)<comp.cte) {
    
    comp.cte <- length(v.ID)
  
  }
  

plsda <- plsda(filter,Tto,ncomp=comp.cte,scale=F)
  set.seed(30)
  perf.pls <- perf(plsda, validation = "Mfold",criterion="all",folds = 4,
                   progressBar = F,nrepeat = 100,)

  
  comp <- perf.pls[["choice.ncomp"]][6]
  err.n <- perf.pls[["error.rate"]]$BER[(comp.cte*2)+comp]
  sd.n <- perf.pls[["error.rate.sd"]]$BER[(comp.cte*2)+comp]
  err.total <- err.n + sd.n
}


#Balanced error rate of the model

err - sd.f
#RAW 0.319416          #LOG 0.2640489     # CLR 0.3036096

sd.f
#RAW 0.03659011      #LOG 0.04620047       # CLR 0.05153868

#Number of variables included in the model

length(vf) #RAW 5  #LOG 39 #CLR 13


#Save (CHOOSE LOG OR CLR)
write.table(vf,"relevant_MAGs_raw.txt",row.names = T,col.names = T,quote = F,sep="\t")
write.table(vf,"relevant_MAGs_log.txt",row.names = T,col.names = T,quote = F,sep="\t")
write.table(vf,"relevant_MAGs_clr.txt",row.names = T,col.names = T,quote = F,sep="\t")

#Add bin names to the saved tables (upload first the prior saved tables)
inverse <- as.data.frame(t(input_w_names))
leyenda_bins <- data.frame(codeBin = rownames(inverse), Bin = inverse$Bin)

relevant_MAGs_clr_names <- merge(relevant_MAGs_clr, leyenda_bins, by.x = "x", by.y = "codeBin", all.x = TRUE)
relevant_MAGs_log_names <- merge(relevant_MAGs_log, leyenda_bins, by.x = "x", by.y = "codeBin", all.x = TRUE)
write.table(relevant_MAGs_clr_names,"relevant_MAGs_clr.txt",row.names = T,col.names = T,quote = F,sep="\t")
write.table(relevant_MAGs_log_names,"relevant_MAGs_log.txt",row.names = T,col.names = T,quote = F,sep="\t")
write.table(leyenda_bins,"leyenda_bins.txt",row.names = T,col.names = T,quote = F,sep="\t")

#Final model
View(loadings)

filter <- input_PLSDA[,colnames(input_PLSDA) %in% vf]
View(filter)
Tto = as.factor(metadata$Tto)

#PLS-Plot of the final model

color <- c("#00008B","#EE0000","#2F2F28","#F0EFF0")

if (comp.f > 1) {
  
  plsda <- plsda(filter,Tto,ncomp=comp.f,scale=F)

  plotIndiv(plsda,ind.names = TRUE, legend=TRUE,style = "ggplot2",rep.space = "X-variate",
             ellipse = TRUE, centroid=FALSE,title = 'MAGP log',
             X.label = 'Comp 1', Y.label = 'Comp 2',col = color[1:2],abline=TRUE,cex = c(3,3),point.lwd = 1,
             size.title = rel(1.2), size.subtitle = rel(1.2), size.xlabel = rel(1),
             size.ylabel = rel(1.2), size.axis = rel(1), size.legend = rel(1),
             size.legend.title = rel(1),
             legend.title = "Class",alpha=1)+ theme_classic()
}


#PCA plot final model

#PCA-Plot of final model
pca <- prcomp(filter,scale=F)
fviz_pca_ind(pca,axes=c(1,2),geom = c("point","text"),col.ind=Tto,addEllipses = T,palette=color,ellipse.level = 0.95,pointsize = 1)

#Quality of the model

data.RF <- cbind(filter,Tto)
data.RF = as.data.frame(data.RF)
confusion.total <- matrix(ncol=2, nrow=2, 0)
x.total <- NULL
for (i in 1:10000) {
  
  sample = sample.split(data.RF$Tto, SplitRatio = .70)
  train = subset(data.RF, sample == TRUE)
  test  = subset(data.RF, sample == FALSE)
  dim(train)
  dim(test)
  
  x <- data.frame(table(test$Tto))
  x.total <- rbind(x,x.total)
  plsda.train <- plsda(train[, -ncol(train)], train$Tto,ncomp=comp.f,scale=F)
  test.predict <- predict(plsda.train, test[,-ncol(test)], dist = "mahalanobis.dist")
  prediction <- test.predict$class$mahalanobis.dist[,comp.f]
  
  confusion.mat <- get.confusion_matrix(truth = test$Tto, predicted = prediction)
  confusion.total <- confusion.total + confusion.mat
}

print(confusion.mat)
# LOG predicted.as.control predicted.as.stress
#control                    8                   1
#stress                     6                   3

#CLR predicted.as.1 predicted.as.2
#control                    4                   5
#stress                     2                   7

print(confusion.total)
# LOG predicted.as.control predicted.as.stress
#control                64436               25564
#stress                 23869               66131

#CLR  predicted.as.control predicted.as.stress
#control                58530               31470
#stress                 24576               65424

control <- sum(x.total$Freq[x.total$Var1 == "control"])  

stress <- sum(x.total$Freq[x.total$Var1 == "stress"])

confusion.total[1,] <- 100*confusion.total[1,] / control
confusion.total[2,] <- 100*confusion.total[2,] / stress

confusion.total


######## LOG  quality control VALUES

#             predicted.as.control predicted.as.stress
#control             71.59556            28.40444
#stress              26.52111            73.47889

#TN (Verdaderos Negativos) = 71.59556
#FP (Falsos Positivos) = 28.40444 
#FN (Falsos Negativos) = 26.52111
#TP (Verdaderos Positivos) = 73.47889

#Accuracy: 72.53722 %
#Sensitivity: 73.47889 %
#Specificity: 71.59556 %
#Precision: 72.12062 %
#F1 Score: 72.79342 %

######## CLR  quality control VALUES
#           predicted.as.control predicted.as.stress
#control             65.03333            34.96667
#stress              27.30667            72.69333

#True Negatives (TN): 65.03333  
#False Positives (FP): 34.96667
#False Negatives (FN): 27.30667 
#True Positives (TP): 72.69333

#Accuracy: 68.86333%
#Sensitivity: 72.69333%
#Specificity: 65.03333%
#Precision: 67.52956%
#F1 Score: 69.98789%


###########

# Extract values from confusion matrix
TN <- confusion.matrix[1, 1]
FP <- confusion.matrix[1, 2]
FN <- confusion.matrix[2, 1]
TP <- confusion.matrix[2, 2]

# Calculate metrics
accuracy <- (TP + TN) / (TP + TN + FP + FN)
sensitivity <- TP / (TP + FN)  # Recall
specificity <- TN / (TN + FP)
precision <- TP / (TP + FP)
f1_score <- 2 * ((precision * sensitivity) / (precision + sensitivity))

# Print results
cat("True Negatives (TN):", TN, "\n")
cat("False Positives (FP):", FP, "\n")
cat("False Negatives (FN):", FN, "\n")
cat("True Positives (TP):", TP, "\n")
cat("Accuracy:", accuracy * 100, "%\n")
cat("Sensitivity (Recall):", sensitivity * 100, "%\n")
cat("Specificity:", specificity * 100, "%\n")
cat("Precision:", precision * 100, "%\n")
cat("F1 Score:", f1_score * 100, "%\n")

##################################

final.splsda <- plsda(filter,Tto,ncomp=comp.f,scale=F)

plotVar(final.splsda, comp = c(1,2), cex = 3)


auc.splsda = auroc(final.splsda ,roc.comp = 1,plot=TRUE)
M<-auc.splsda$graph.Comp1 +
  scale_colour_manual(values =c("green","blue","red") ) + 
  ggtitle("") +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black")) +
  theme(axis.text.x=element_text(hjust=0.7))+
  theme(axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9))+
  theme(axis.text.y=element_text(hjust=0.7,size=9))+
  theme(legend.text=element_text(size=9))+
  theme(panel.background = element_rect(fill = "transparent"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.text.x = element_text(size=9),
        plot.background = element_rect(
          fill = "transparent"))+ theme(legend.title = element_blank()) 
		  
pdf("ROCPC1_MAGabundance_fromprot.pdf", width=6, height=6)
M
dev.off()


##
set.seed(40) # for reproducibility
perf.srbct <- perf(final.splsda, validation = "Mfold", folds = 4,
                   dist = 'max.dist', nrepeat = 100,
                   progressBar = FALSE) 
				   
pdf("perf.srbct_clr.pdf")
plot(perf.srbct, col = color.mixo(5))
dev.off()

#Selected variables per components

PC1 = selectVar(final.splsda, comp = 1)$value
PC2 = selectVar(final.splsda, comp = 2)$value

write.table(PC1, file="MAGs_PC1_clr_c70.txt", quote=F, sep="\t")
write.table(PC2, file="MAGs_PC2_clr_c70.txt", quote=F, sep="\t")

pdf("contibutionPC1_CLR_magfromprot.pdf", width=7, height=5)
plotLoadings(final.splsda, comp = 1, title = '', contrib = 'max', method = 'mean', legend.color=c("blue","red"), size.name=0.8, size.legend=0.8, ndisplay=30)
dev.off()


#heatmap

dose.col <- color.mixo(as.numeric(as.factor(metadata$Tto)))

dose.col <- gsub("#F68B33", "red", dose.col)
dose.col <- gsub("#388ECC", "blue", dose.col)

legend=list(legend = levels(Tto), # set of classes
            col = unique(dose.col), # set of colours
            title = "Class", # legend title
            cex = 0.7)


cim(final.splsda, comp=1, title ="Component 1", row.sideColors=dose.col)




################## ADD taxonomy to relevant_MAG
taxonomy <- taxonomia_MetaP2

relevant_MAGs_clr <- relevant_MAGs_clr %>% rename(user_genome = Bin)

relevant_MAGs_log <- relevant_MAGs_log %>% rename(user_genome = Bin)

relevant_56MAGs <- discriminantmags %>% rename(user_genome = MAG)
relevant_MAGs_clr_tax <- merge(relevant_MAGs_clr, taxonomy, by = "user_genome", all = FALSE)
relevant_MAGs_log_tax <- merge(relevant_MAGs_log, taxonomy, by = "user_genome", all = FALSE)

View(relevant_MAGs_log_tax)
relevant_MAGs_clr_tax <- relevant_MAGs_clr_tax %>%
  select(x, user_genome, classification, other_related_references)

relevant_MAGs_log_tax <- relevant_MAGs_log_tax %>% 
  select(x, user_genome, classification, other_related_references)

View(relevant_MAGs_clr_tax)
View(relevant_MAGs_log_tax)

write.table(relevant_MAGs_clr_tax, "relevant_MAGs_clr_tax.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(relevant_MAGs_log_tax, "relevant_MAGs_log_tax.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Filas completas que coinciden
filas_comunes_MAGP_MAGG <- merge(relevant_MAGs_clr_tax, relevant_56MAGs, 
                       by = "user_genome")
# Solo los user_genome comunes
comunes <- unique(filas_comunes$user_genome)
comunes_MAGGMAGP <- unique(filas_comunes_MAGP_MAGG$user_genome)

cat("Número de user_genome comunes:", length(comunes), "\n")
cat("Porcentaje en clr_tax:", mean(relevant_MAGs_clr_tax$user_genome %in% comunes)*100, "%\n")
#Porcentaje en clr_tax: 48.07692 %
cat("Porcentaje en log_tax:", mean(relevant_MAGs_log_tax$user_genome %in% comunes)*100, "%\n")
#Porcentaje en log_tax: 98.03922 %

#################### % of common MAGs in discriminant MAGG and MAGP (both clr)
cat("Porcentaje en 56MAGs:", mean(relevant_56MAGs$user_genome %in% comunes_MAGGMAGP)*100, "%\n")
#Porcentaje en 56MAGs: 33.92857 %
cat("Porcentaje en clr_tax:", mean(relevant_MAGs_clr_tax$user_genome %in% comunes_MAGGMAGP)*100, "%\n")
#Porcentaje en clr_tax: 18.26923 %

filas_comunes_MAGP_MAGG<- filas_comunes_MAGP_MAGG %>%
  select(x, user_genome, Indicator,loadings, Completeness, Contamination, classification, other_related_references)

write.table(filas_comunes_MAGP_MAGG, "filas_comunes_MAGP_MAGG.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
