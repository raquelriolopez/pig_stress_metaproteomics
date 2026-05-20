library(tidyverse)
library(mixOmics)
library("factoextra")
library("caTools")
library("easyCODA")


#input

summedLFQ_functions_filtfreq <- read.table("summedLFQ_functions_filtfreq.txt", header = T, sep = "\t") 

##########PLS-DA for functions

#First, we need to adapt the format to PLS-DA

wide_data <- summedLFQ_functions_filtfreq %>%
  dplyr::select(ko_id, SampleID, summedLFQ) %>%  #For some reason, select is masked by other packages that I'm using, but you can erase the explicit call
  tidyr::pivot_wider(
    names_from = SampleID,
    values_from = summedLFQ
  )

View(wide_data)

#Normalisation
#LOG
wide_data_log <- wide_data
wide_data_log[, -1] <- log2(wide_data[, -1] + 1) ####to avoid log(0) we add +1
write.table(wide_data_log, file = "functions_filtfreq_log.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
View(wide_data_log)
#CLR
wide_data_clr<- wide_data[, -1]
wide_data_clr <- CLR(wide_data_clr + 0.0001, weight = TRUE)$LR
wide_data_clr <- data.frame(ko_id = wide_data$ko_id, wide_data_clr)
View(wide_data_clr)
write.table(wide_data_clr, file = "functions_filtfreq_clr.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


#Extract the ko_id names
koid_names <- wide_data$ko_id

View(koid_names)

#transpose

input_PLSDA <- wide_data_log %>%
  dplyr::select(-ko_id) %>%          
  as.data.frame() %>%         
  t() %>%          
  as.data.frame()

input_plsda_w_names <- rbind(ko_id = koid_names, input_PLSDA)


View(input_PLSDA)
View(input_plsda_w_names)

#metadata

metadata = read.table("metadataMetaP.txt", header=T, row.names=1, dec=",")
metadata = metadata[rownames(input_plsda_w_names[-1, ]), ]
Tto = as.factor(metadata$Tto)


#PLSDA
plsda <- plsda(input_PLSDA, Tto, ncomp=10, scale=F)

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

err - sd.f
#RAW 0.3809733 #LOG 0.05134038 #CLR 0.05761402

sd.f
#RAW 0.05293757 #LOG 0.02665495 #CLR 0.01998884

#Number of variables included in the model

length(vf) #RAW 25 #LOG 251 #CLR 232

#Add bin names and save
inverse <- as.data.frame(t(input_plsda_w_names))
leyenda_koids <- data.frame(code_ko_id = rownames(inverse), ko_id = inverse$ko_id)
relevant_koid_log_names <- merge(vf, leyenda_koids, by.x = "x", by.y = "code_ko_id", all.x = TRUE)
relevant_koid_clr_names <- merge(vf, leyenda_koids, by.x = "x", by.y = "code_ko_id", all.x = TRUE)
write.table(relevant_koid_log_names,"relevant_koid_log_names.txt", row.names = T, col.names = T, quote = F, sep="\t")
write.table(relevant_koid_clr_names,"relevant_koid_clr_names.txt", row.names = T, col.names = T, quote = F, sep="\t")

#Final model

filter <- input_PLSDA[,colnames(input_PLSDA) %in% vf]
View(filter)
Tto = as.factor(metadata$Tto)

#PLS-Plot of the final model

color <- c("#00008B","#EE0000","#2F2F28","#F0EFF0")

if (comp.f > 1) {
  
  plsda <- plsda(filter,Tto,ncomp=comp.f,scale=F)
  
  plotIndiv(plsda,ind.names = TRUE, legend=TRUE,style = "ggplot2",rep.space = "X-variate",
            ellipse = TRUE, centroid=FALSE,title = 'KO_ID LOG',
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

#RAW        predicted.as.control predicted.as.stress
#control                    6                   3
#stress                     7                   2

#LOG predicted.as.control predicted.as.stress
#control                    8                   1
#stress                     0                   9

#CLR     predicted.as.control predicted.as.stress
#control                    8                   1
#stress                     0                   9

print(confusion.total)

#RAW     predicted.as.control predicted.as.stress
#control                51980               38020
#stress                 33628               56372

#LOG     predicted.as.control predicted.as.stress
#control                86115                3885
#stress                  7645               82355

#CLR     predicted.as.control predicted.as.stress
#control                87456                2544
#stress                  9162               80838

control <- sum(x.total$Freq[x.total$Var1 == "control"])  

stress <- sum(x.total$Freq[x.total$Var1 == "stress"])

confusion.total[1,] <- 100*confusion.total[1,] / control
confusion.total[2,] <- 100*confusion.total[2,] / stress

confusion.total
#RAW     predicted.as.control predicted.as.stress
#control             57.75556            42.24444
#stress              37.36444            62.63556

#LOG     predicted.as.control predicted.as.stress
#control            95.683333            4.316667
#stress              8.494444           91.505556

#CLR     predicted.as.control predicted.as.stress
#control             97.17333            2.826667
#stress              10.18000           89.820000


######## LOG  quality control VALUES

#TN (True Negatives)  = 95.683333
#FP (False Positives) = 4.316667
#FN (False Negatives) = 8.494444
#TP (True Positives)  = 91.505556

#Accuracy    = 93.59444%
#Sensitivity = 91.50556%
#Specificity = 95.68333%
#Precision   = 95.49444%
#F1 Score    = 93.468%

######## CLR quality control VALUES

#TN (True Negatives)  = 97.17333%
#FP (False Positives) = 2.826667%
#FN (False Negatives) = 10.18000%
#TP (True Positives)  = 89.82000%

#Accuracy    = 93.49667%
#Sensitivity = 89.82000%
#Specificity = 97.17333%
#Precision   = 96.94872%
#F1 Score    = 93.248%


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

pdf("ROCPC1_plsda_CLR.pdf", width=6, height=6)
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

write.table(PC1, file="koids_PC1_LOG.txt", quote=F, sep="\t")
write.table(PC2, file="koids_PC2_LOG.txt", quote=F, sep="\t")

write.table(PC1, file="koids_PC1_CLR.txt", quote=F, sep="\t")
write.table(PC2, file="koids_PC2_CLR.txt", quote=F, sep="\t")

plotLoadings(final.splsda, comp = 1, title = '', contrib = 'max', method = 'mean', legend.color=c("blue","red"), size.name=0.8, size.legend=0.8, ndisplay=30)
dev.off()
