rm(list=ls())
out_dir <- paste("results/")
dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)


#### load libraries ##########
library(stringr)
library(SummarizedExperiment)
library(ggplot2)
library(Cairo)
library(limma)
library(edgeR)
library(reshape2)
library(tidyverse)
##############################

map_table <- read.table("../data/map_batch_1_2_3.txt",sep="\t",header=TRUE)
coldata <-map_table
rownames(coldata) <- coldata$Samples
coldata$names  <-  coldata$Sample 
coldata$Sample <- coldata$Sample %>% as.numeric
coldata <- coldata %>% filter(!grepl("RasG12D",Cell_line))
coldata <- coldata %>% mutate(Condition=ifelse(coldata$Condition=="bounce","bb" ,coldata$Condition)) %>%
    mutate(Cell_line=str_replace(coldata$Cell_line," ","_"))

for (i in 1:(dim(coldata)[2]))
{
    coldata[,i] <- as.factor(coldata[,i])
}
coldata$condition <- factor(coldata$Condition,levels=c("CTRL","AA","bb"))

########### to read from the Expression matrix instead of tximeta #########
gseLimma.path <- file.path("data/count_matrix.tsv")
gseLimma.table <- read.table(gseLimma.path,header=TRUE,sep="\t")
gseLimma.table <- gseLimma.table %>% column_to_rownames("gene.id")
gseLimma.mat <-  gseLimma.table %>% as.matrix
colnames(gseLimma.mat) <- gseLimma.mat %>% colnames %>% str_replace(.,"X","")
gseLimma.mat <-  gseLimma.mat[,order(as.numeric(colnames(gseLimma.mat)))]
gseLimma.mat <- gseLimma.mat[, colnames(gseLimma.mat) %in% (coldata$names %>% as.character) ]
gseLimma <- DGEList(gseLimma.mat,lib.size = colSums(gseLimma.mat))

combo <- paste(coldata$Condition,coldata$Cell_line,coldata$Week,sep="_") %>%
    str_replace(.,"CTRL","C") %>%
    str_replace(.,"AA","A") %>%
    str_replace(.,"bb","B")  %>%
    str_replace(.,"_$","")

design <- model.matrix( ~ 0 + combo + coldata$Batch)
colnames(design) <- gsub("combo","", colnames(design)) %>%  str_replace("coldata\\$Batch","Batch")



contrasts <- makeContrasts(
    mP4_p53KO_5= A_mP4_p53KO_5 - C_mP4_p53KO_1,
    mP4_p53R273H_5 = A_mP4_p53R273H_5 - C_mP4_p53R273H_1,
    mP4_5 = A_mP4_5 - C_mP4_1,
    mN10_5 = A_mN10_5 - C_mN10_1,
    ##
    mP4_p53KO_8= A_mP4_p53KO_8 - C_mP4_p53KO_1,
    mP4_p53R273H_8 = A_mP4_p53R273H_8 - C_mP4_p53R273H_1,
    mP4_8 = A_mP4_8 - C_mP4_1,
    mN10_8 = A_mN10_8 - C_mN10_1,
    ##
    mP4_p53KO_11= A_mP4_p53KO_11 - C_mP4_p53KO_1,
    mP4_p53R273H_11 = A_mP4_p53R273H_11 - C_mP4_p53R273H_1,
    mP4_11 = A_mP4_11 - C_mP4_1,
    mN10_11 = A_mN10_11 - C_mN10_1,
    ##
    mP4_p53KO_B= B_mP4_p53KO_11 - C_mP4_p53KO_1,
    mP4_p53R273H_B = B_mP4_p53R273H_11 - C_mP4_p53R273H_1,
    mP4_B = B_mP4_11 - C_mP4_1,
    mN10_B = B_mN10_11 - C_mN10_1,
    levels = design)


idx <- filterByExpr(gseLimma,min.count=55,min.total.counts=40 ,design) # filtering genes according to counts 
dds <- gseLimma[idx,]
dds <- calcNormFactors(dds) # this will be used by Limma 



####### voomwithqualityweights ##############
out <- file.path(out_dir,"voomWQW.pdf") #### VoomWQW and plot
CairoPDF(width = 11 , height = 7.5, file=out, pointsize=12, bg = "white")
v <- voomWithQualityWeights(dds,design,plot=TRUE)
dev.off()

fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrasts=contrasts)
efit <- eBayes(fit2)

out <- file.path(out_dir,"pval_histo_contrasts_voomWQW.png")  ##### plot pVal t-statistics
Cairo(width = 1700, height = 1100, file=out, type="png", pointsize=12, bg = "white", canvas = "white", units = "px", dpi = "auto")
par(mfrow=c(2,3))
colnames(contrasts) %>%  lapply(.,topTable,fit=efit,n=Inf) %>% lapply(.,dplyr::select,P.Value) %>% lapply(unlist)   %>% mapply(hist,.,main=colnames(contrasts))
dev.off()

deg.names <- colnames(contrasts) %>% # retrieving all the genes that have are significant at least once
    lapply(.,topTable,fit=efit,n=Inf,p.value=0.05) %>%
    lapply(rownames) %>%
    unlist %>%
    unique

deg.nb <- colnames(contrasts) %>% # nb of genes (not only the significant ones) (not the best way to code) 
    lapply(.,topTable,fit=efit,n=Inf) %>%
    lapply(nrow)  %>%
    unlist

deg.nb.c <- mapply(rep,contrasts %>% colnames,deg.nb,SIMPLIFY=FALSE) %>% # creating a column for the contrasts before binding all the rows
    unlist
names(deg.nb.c) <- NULL

deg.tab <- colnames(contrasts) %>% # binding the rows for all the contrasts
    lapply(.,topTable,fit=efit,n=Inf) %>%
    lapply(as.data.frame) %>% 
    lapply(rownames_to_column) %>% 
    do.call(rbind,.) %>%
    mutate(contrast=deg.nb.c) %>%
    dplyr::rename(gene.id=rowname) %>%
    dplyr::filter(gene.id %in% deg.names)

out.table <- file.path(out_dir,"deg_table_voomWQW.tsv")
write.table(deg.tab,out.table,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

