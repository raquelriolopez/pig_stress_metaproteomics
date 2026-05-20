library(tidyverse)
library(mixOmics)
library("factoextra")
library("caTools")
library("easyCODA")


#input

MAGG_abundance <- read_tsv("MAGG_abundance.tsv") #Consider the file's format when reading it. (.txt, .tsv, .xlsx...)

MAGG_abundance_raw <- as.data.frame(t(MAGG_abundance)) #Adapt the format of your table so that it can be analysed by PLS-DA.
View(MAGG_abundance_raw)

#Transform only numeric values and keep bin_names RAW ###### This step is only needed if your values are not numeric and you have problems with transformations
bin_names <- MAGG_abundance_raw[1, ]
mag_abund_numeric <- MAGG_abundance_raw[-1, ]
row_names_num <- rownames(mag_abund_numeric)
mag_abund_numeric <- as.data.frame(lapply(mag_abund_numeric, as.numeric))
rownames(mag_abund_numeric) <- row_names_num
View(mag_abund_numeric)
MAGG_raw <- mag_abund_numeric
View(MAGG_raw)

#Log transformation
MAGG_log <- log2(MAGG_raw + 0.0001)
View(MAGG_log)
#CLR transformation   ### for microbiota abundances this is the best option (for compositional data)
MAGG_clr <- CLR(MAGG_raw + 0.0001, weight = TRUE)$LR
View(MAGG_clr)

#Specify the input here (clr_transformed or log_transformed). I prefer to create 'input_PLSDA' to standarise the code for future analyses
input_PLSDA <- MAGG_raw
input_PLSDA <- MAGG_log
input_PLSDA <- MAGG_clr

View(input_PLSDA)

#Rebuild with bin_names
input_w_names<- rbind(Bin = bin_names, input_PLSDA)

write.table(input_w_names,"MAGG_log_w_names.txt",row.names = T,col.names = T,quote = F,sep="\t")
write.table(input_w_names,"MAGG_clr_w_names.txt",row.names = T,col.names = T,quote = F,sep="\t")

#Metadata (You need a file with the information about your samples, with the factors that you want to consider, in our case the treatment 'Tto')
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
#RAW  0.3496162     #LOG 0.09695217 #CLR 0.09700779

sd.f
#RAW  0.04172795    #LOG 0.03395727 #CLR 0.03226612

#Number of variables included in the model

length(vf) #RAW 60 #LOG 59 #CLR 55


#Save (CHOOSE LOG OR CLR)
write.table(vf,"discriminant_MAGs_raw.txt",row.names = T,col.names = T,quote = F,sep="\t")
write.table(vf,"discriminant_MAGs_log.txt",row.names = T,col.names = T,quote = F,sep="\t")
write.table(vf,"discriminant_MAGs_clr.txt",row.names = T,col.names = T,quote = F,sep="\t")

#Add bin names to the saved tables (upload first the prior saved tables if you want)
inverse <- as.data.frame(t(input_w_names))
leyenda_bins <- data.frame(codeBin = rownames(inverse), Bin = inverse$Bin)

discriminant_MAGs_log <- vf
discriminant_MAGs_clr <- vf
discriminant_MAGs_clr_names <- merge(discriminant_MAGs_clr, leyenda_bins, by.x = "x", by.y = "codeBin", all.x = TRUE)
discriminant_MAGs_log_names <- merge(discriminant_MAGs_log, leyenda_bins, by.x = "x", by.y = "codeBin", all.x = TRUE)
write.table(discriminant_MAGs_clr_names,"discriminant_MAGs_clr.txt",row.names = T,col.names = T,quote = F,sep="\t")
write.table(discriminant_MAGs_log_names,"discriminant_MAGs_log.txt",row.names = T,col.names = T,quote = F,sep="\t")
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
            ellipse = TRUE, centroid=FALSE,title = 'MAGG clr',
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
#control                    7                   2
#stress                     0                   9

#CLR predicted.as.1 predicted.as.2
#control              7              2
#stress              1              9

print(confusion.total)
# LOG predicted.as.control predicted.as.stress
#control                79515               10485
#stress                  8589               81411

#CLR  predicted.as.control predicted.as.stress
#control          79569          10431
#stress           8692          81308

control <- sum(x.total$Freq[x.total$Var1 == "1"])  

stress <- sum(x.total$Freq[x.total$Var1 == "2"])

confusion.total[1,] <- 100*confusion.total[1,] / control
confusion.total[2,] <- 100*confusion.total[2,] / stress

confusion.total


######## LOG  quality control VALUES

#             predicted.as.control predicted.as.stress
#control            88.350000            11.65000
#stress              9.543333            90.45667

#TN (Verdaderos Negativos) = 88.350000
#FP (Falsos Positivos) = 11.65000
#FN (Falsos Negativos) = 9.543333
#TP (Verdaderos Positivos) = 90.45667


######## CLR  quality control VALUES
#           predicted.as.control predicted.as.stress
#control      88.410000       11.59000
#stress       9.657778       90.34222

#True Negatives (TN): 88.410000 
#False Positives (FP): 11.59000
#False Negatives (FN): 9.657778 
#True Positives (TP): 90.34222


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

write.table(PC1, file="MAGGs_PC1_log.txt", quote=F, sep="\t")
write.table(PC2, file="MAGGs_PC2_log.txt", quote=F, sep="\t")
write.table(PC1, file="MAGGs_PC1_clr.txt", quote=F, sep="\t")
write.table(PC2, file="MAGGs_PC2_clr.txt", quote=F, sep="\t")

pdf("contibutionPC1_CLR_magfromprot.pdf", width=7, height=5)
plotLoadings(final.splsda, comp = 1, title = '', contrib = 'max', method = 'mean', legend.color=c("blue","red"), size.name=0.8, size.legend=0.8, ndisplay=30)
dev.off()

###Check if your positive/negative loadings have been linked to the control or the stress condition, it depends on your data.
#log loading positivo stress loading negativo control
#clr loading positivo stress loading negativo control

#heatmap

dose.col <- color.mixo(as.numeric(as.factor(metadata$Tto)))

dose.col <- gsub("#F68B33", "red", dose.col)
dose.col <- gsub("#388ECC", "blue", dose.col)

legend=list(legend = levels(Tto), # set of classes
            col = unique(dose.col), # set of colours
            title = "Class", # legend title
            cex = 0.7)


cim(final.splsda, comp=1, title ="Component 1", row.sideColors=dose.col)


####Names
loadings <- final.splsda$loadings$X[, 1:comp.f]  # Usamos el número óptimo de componentes
loadings_df <- as.data.frame(loadings)
loadings_df$codeBin <- rownames(loadings_df)
loadings_df <- merge(loadings_df, leyenda_bins, by = "codeBin", all.x = TRUE)

head(leyenda_bins)
head(loadings_df)


# Asumiendo que el primer componente se llama "comp1":
result_table <- loadings_df %>%
  mutate(
    Association = ifelse(comp1 < 0, "control", "stress"),
    Loading_Value = comp1  # Guardamos el valor del loading para referencia
  ) %>%
  dplyr::select(codeBin, Association, Loading_Value, Bin, everything())

# Ordenar por magnitud del loading (importancia)
result_table <- result_table %>%
  arrange(desc(abs(Loading_Value)))

###AÑADIR TAXONOMÍA
library(readxl)
taxonomy_MAGG <- read_excel("taxonomy_MAG_C70C10_MetaP2.xlsx") %>%
  rename(
    Bin = user_genome
  )
View(taxonomy_MAGG)

result_table <- result_table %>%
  left_join(taxonomy_MAGG, by = "Bin") %>%
  dplyr::select(., Bin, Loading_Value, Association, classification, closest_genome_ani)

view(result_table)

# Guardar resultados
write.table(result_table, 
            "PLSDA_MAGG_LOG_control_stress.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

write.table(result_table, 
            "PLSDA_MAGG_CLR_control_stress.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
View(result_table)

sum(result_table$Association == "stress")  #21 LOG MAGS discriminate stress
                                           #18 CLR MAGs discriminate stress
sum(result_table$Association == "control") #38 LOG MAGs discriminate control
                                           #37 CLR MAGs discriminate control


MAGG_control <- subset(result_table, Association == "control")
MAGG_stress <- subset(result_table, Association == "stress")

write.table(MAGG_control, 
            "MAGG_LOG_control_plsda.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

write.table(MAGG_stress, 
            "MAGG_LOG_stress_plsda.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
###clr
write.table(MAGG_control, 
            "MAGG_CLR_control_plsda.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

write.table(MAGG_stress, 
            "MAGG_CLR_stress_plsda.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
