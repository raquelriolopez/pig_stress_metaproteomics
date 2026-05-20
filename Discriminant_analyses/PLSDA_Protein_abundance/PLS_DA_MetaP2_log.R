#PLS-DA
library(mixOmics)
BiocManager::install("mixOmics")
library("factoextra")
library("caTools")
library(tidyverse)

metap= read.table("MetaP2_raw_filtfreq.tsv", header=T, check.names=F)
metap = t(metap)
### not needed they already have the names colnames(metap) <- paste0("P", 1:ncol(metap))

metadata = read.table("metadataMetaP.txt", header=T, row.names=1, dec=",")
metadata = metadata[rownames(metap),]

Tto = as.factor(metadata$Tto)

#Log transformation
metap = metap+0.0001
metap = log2(metap)
View(metap)
#PLSDA
plsda <- plsda(metap,Tto,ncomp=10,scale=F)


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
  
  filter <- metap[,colnames(metap) %in% v.ID]
  
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
0.002575083 #0.002736374


sd.f
0.006161753 #0.006734056


#Number of variables included in the model

length(vf)
726 #1001

View(metap)

#Save
write.table(vf,"relevant_metaboliteslog_1001_new.txt",row.names = T,col.names = T,quote = F,sep="\t")

#Final model

filter <- metap[,colnames(metap) %in% vf]
Tto = as.factor(metadata$Tto)

save(filter, Tto, comp.f, color, err, err.n, err.total, p, sd.f, sd.n, v.ID, vf, comp.cte, metadata, metap, perf.pls, plsda, v.select, vip, file = "recreate_plot.RData")

load("recreate_plot.RData")

dim(filter)  
table(Tto)   

View(Tto)

#PLS-Plot of the final model

color <- c("#00008B","#EE0000","#2F2F28","#F0EFF0")

if (comp.f > 1) {
  
  plsda <- plsda(filter,Tto,ncomp=comp.f,scale=F)
  png("PLSDA_plot.png", width = 800, height = 600)
  plotIndiv(plsda,ind.names = TRUE, legend=TRUE,style = "ggplot2",rep.space = "X-variate",
             ellipse = TRUE, centroid=FALSE,title = 'PLS-DA Proteins',
             X.label = 'Comp 1', Y.label = 'Comp 2',col = color[1:2],abline=TRUE,cex = c(3,3),point.lwd = 1,
             size.title = rel(1.2), size.subtitle = rel(1.2), size.xlabel = rel(1),
             size.ylabel = rel(1.2), size.axis = rel(1), size.legend = rel(1),
             size.legend.title = rel(1),
             legend.title = "Condition",alpha=1)+ theme_classic()
  dev.off()
  }


#PCA plot final model

#PCA-Plot of final model
pca <- prcomp(filter,scale=F)
fviz_pca_ind(pca,axes=c(1,2),geom = c("point","text"),col.ind=Tto,addEllipses = T,palette=color,ellipse.level = 0.95,pointsize = 1)


#Quality of the model

data.RF <- cbind(filter,Tto)
data.RF = as.data.frame(data.RF)
confusion.total <- matrix(ncol=2,nrow=2,0)
x.total <- NULL
for (i in 1:10000) {
  
  sample = sample.split(data.RF$Tto, SplitRatio = .70)
  train = subset(data.RF, sample == TRUE)
  test  = subset(data.RF, sample == FALSE)
  dim(train)
  dim(test)
  
  x <- data.frame(table(test$Tto))
  x.total <- rbind(x,x.total)
  plsda.train <- plsda(train[,-ncol(train)],train$Tto,ncomp=comp.f,scale=F)
  test.predict <- predict(plsda.train,test[,-ncol(test)],dist = "mahalanobis.dist")
  prediction <- test.predict$class$mahalanobis.dist[,comp.f]
  
  confusion.mat <- get.confusion_matrix(truth = test$Tto,predicted =prediction )
  confusion.total <- confusion.total + confusion.mat

  }


control <- sum(x.total$Freq[x.total$Var1 == 1])
stress <- sum(x.total$Freq[x.total$Var1 == 2])

confusion.total[1,] <- 100 * confusion.total[1,] / control
confusion.total[2,] <- 100 * confusion.total[2,] / stress

confusion.total

#REPETIDO 
predicted.as.1 predicted.as.2
1     99.9400000        0.06000
2      0.7988889       99.20111

TP (True Positives) = Predicciones correctas de "Control" (99.9400000)
TN (True Negatives) = Predicciones correctas de "Stress" (99.20111)
FP (False Positives) = Predicciones incorrectas de "Control" como "Stress" (0.06000)
FN (False Negatives) = Predicciones incorrectas de "Stress" como "Control" (0.7988889)

Accuracy: 99.57%
Sensitivity (Stress): 99.21%
Specificity (Control): 99.94%
Precision (Stress): 99.94%
F1 Score (Stress): 99.50%


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
		  
pdf("ROCPC1_mETAg2GENESlog.pdf", width=6, height=6)
M
dev.off()


##
set.seed(40) # for reproducibility
perf.srbct <- perf(final.splsda, validation = "Mfold", folds = 4,
                   dist = 'max.dist', nrepeat = 100,
                   progressBar = FALSE) 
				   
pdf("perf.srbct_log.pdf")
plot(perf.srbct, col = color.mixo(5))
dev.off()

#Selected variables per components

PC1 = selectVar(final.splsda, comp = 1)$value
PC2 = selectVar(final.splsda, comp = 2)$value

write.table(PC1, file="metapPC1_log_new.txt", quote=F, sep="\t")
write.table(PC2, file="metapPC2_log_new.txt", quote=F, sep="\t")

pdf("contibutionPC1_metap_log.pdf", width=7, height=5)
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

####Get names

Leyenda_MetaP2 <- read.table("Leyenda_MetaP2.txt", header = TRUE, sep = "\t", quote = "")
Leyenda_MetaP2$code_protein <- paste0("P", rownames(Leyenda_MetaP2))
  
View(Leyenda_MetaP2)
####Names
loadings <- final.splsda$loadings$X[, 1:comp.f] 
loadings_df <- as.data.frame(loadings)
loadings_df$code_protein <- rownames(loadings_df)
loadings_df <- merge(loadings_df, Leyenda_MetaP2, by = "code_protein", all.x = TRUE)

View(loadings_df)

# "comp1":
result_table <- loadings_df %>%
  mutate(
    Association = ifelse(comp1 < 0, "control", "stress"),  
    Loading_Value = comp1  
  ) %>%
  dplyr::select(code_protein, Protein, Association, Loading_Value, Description, comp1, comp2, comp3)


result_table <- result_table %>%
  arrange(desc(abs(Loading_Value)))


write.table(result_table, 
            "discriminant_proteins_control_stress_plsda_NEW.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

View(result_table)

result_table <- read.delim("discriminant_proteins_control_stress_plsda_CORRECT.txt", header = TRUE, sep ="\t")

sum(result_table$Association == "stress")





