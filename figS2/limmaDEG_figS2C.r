rm(list=ls())
out_dir <- paste("../results_2/gene_expression_limma_time_series_BBvsW1")


dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)


#### load libraries ##########
library(stringr)
library(SummarizedExperiment)
library(ggplot2)
library(ggplus)
library(Cairo)
library(limma)
library(edgeR)
library(reshape2)
library(tidyverse)
##############################
## salmon.index <- "../../../shared_data/human_genome/salmon_index"
## makeLinkedTxome(indexDir=salmon.index,
##                 source="GENCODE",
##                 organism="Homo sapiens",
##                 release="39",
##                 genome="GRCh38",
##                 fasta="../../../shared_data/human_genome/gencode.v39.transcripts.fa",
##                 gtf="../../../shared_data/human_genome/gencode.v39.primary_assembly.annotation.gtf",
##                 write=FALSE)

quant_files.path <- "/isdata/sandelin/people/arnaud/pH_organoids/data/quant_files_path.txt"
quant_files <- read.table(quant_files.path,header=FALSE) %>% unlist



coldata <- NULL
for (files in quant_files)
{
    coldata <- rbind(coldata,data.frame(files, names=basename(dirname(files)), stringsAsFactors=FALSE))
}
map_table <- read.table("/isdata/sandelin/people/arnaud/pH_organoids/data/map_batch_1_2_3.txt",sep="\t",header=TRUE)
coldata$Sample <- coldata$names %>% as.numeric
coldata <- left_join(coldata,map_table) %>% arrange(as.numeric(names))
rownames(coldata) <- coldata$names
coldata <- coldata %>% filter(!grepl("RasG12D",Cell_line))#    %>% dplyr::filter(Week != 10) %>% dplyr::filter(Batch == 3)
coldata <- coldata %>% mutate(Condition=ifelse(coldata$Condition=="bounce","bb" ,coldata$Condition)) %>%
    mutate(Cell_line=str_replace(coldata$Cell_line," ","_"))

for (i in 1:(dim(coldata)[2]))
{
    coldata[,i] <- as.factor(coldata[,i])
}
coldata$condition <- factor(coldata$Condition,levels=c("CTRL","AA","bb"))


########### to read from the Expression matrix instead of tximeta #########
gseLimma.path <- file.path("../results_2/gene_expression_limma_time_series/Expression_matrix.tsv")
## gseLimma.path <- file.path("../data/counts_with_dots.tsv")
gseLimma.table <- read.table(gseLimma.path,header=TRUE,sep="\t")
rownames(gseLimma.table) <- gseLimma.table[,1]
gseLimma.mat <-  gseLimma.table %>% dplyr::select(-1) %>% as.matrix
colnames(gseLimma.mat) <- gseLimma.mat %>% colnames %>% str_replace(.,"X","")
gseLimma.mat <-  gseLimma.mat[,order(as.numeric(colnames(gseLimma.mat)))]
gseLimma.mat <- gseLimma.mat[, colnames(gseLimma.mat) %in% (coldata$names %>% as.character) ]
gseLimma <- DGEList(gseLimma.mat,lib.size = colSums(gseLimma.mat))

## ##### get list of expressed genes for each time point/geno ######
## coldata.mutate <- coldata %>% dplyr::select(Sample,Cell_line,Week,condition) %>%
##     mutate(Week=ifelse(Week==11 & condition=="bb","B",as.character(Week))) %>%
##     mutate(Cell_line.Week=paste(Cell_line,Week,sep=".")) %>%
##     dplyr::select(-condition)
## all_condi <- coldata.mutate$Cell_line.Week %>% unique
## genes.condi.df <- NULL
## list_condi <- list()
## for (condi in all_condi)
## {
##     list_condi[[condi]] <- coldata.mutate %>%
##         dplyr::filter(Cell_line.Week==condi) %>%
##         .$Sample %>%
##         as.character() %>%
##         as.numeric()
##     gseLimma.condi.geneList <- gseLimma[,colnames(gseLimma) %in% list_condi[[condi]]] %>%
##         .$counts %>% 
##         as.data.frame %>%
##         dplyr::filter(rowSums(. > 20) > 1) %>%
##         rownames()
##     genes.condi.df <- rbind.data.frame(genes.condi.df,cbind.data.frame(gseLimma.condi.geneList,condi))
## }
## colnames(genes.condi.df) <- c("gene.id","condi")
## genes.condi.df <- genes.condi.df %>% separate(condi,sep="\\.",into=c("Cell_line","Week"))

## write.table(genes.condi.df,file.path(out_dir,"Expressed_genes_condi.tsv"),sep="\t",quote=FALSE, row.names=FALSE,col.names=TRUE)

# ########## to read from tximeta (old code) ###############################
## se <- tximeta(coldata,skipSeqinfo=TRUE)
## gse <- summarizeToGene(se)

gseLimma <- DGEList(assays(gse)[["counts"]]) 
rownames(gseLimma) <- rownames(gseLimma) %>% str_replace(.,"\\..*","")

## gseLimma.abundance <- DGEList(assays(gse)[["abundance"]]) 
## rownames(gseLimma.abundance) <- rownames(gseLimma.abundance) %>% str_replace(.,"\\..*","")


combo <- paste(coldata$Condition,coldata$Cell_line,coldata$Week,sep="_") %>%
    str_replace(.,"CTRL","C") %>%
    str_replace(.,"AA","A") %>%
    str_replace(.,"bb","B")  %>%
    str_replace(.,"_$","")

design <- model.matrix( ~ 0 + combo + coldata$Batch)
colnames(design) <- gsub("combo","", colnames(design)) %>%  str_replace("coldata\\$Batch","Batch")

##### to write design matrix
design.df <-  design %>% as.data.frame %>% mutate(sample.id = coldata$names)
write.table(design.df,file.path(out_dir,"design_matrix.tsv"),sep="\t",quote=FALSE, row.names=FALSE)
#####




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

## to write contrast matrix
contrasts.df <- contrasts %>% as.data.frame %>% mutate(levels=rownames(.))
write.table(contrasts.df,file.path(out_dir,"contrasts_matrix.tsv"),sep="\t",quote=FALSE, row.names=FALSE)



idx <- filterByExpr(gseLimma,min.count=55,min.total.counts=40 ,design) # filtering genes according to counts 
dds <- gseLimma[idx,]
dds <- calcNormFactors(dds) # this will be used by Limma 


gseLimma$counts[idx,] %>% as.data.frame %>% mutate(gene.id=rownames(.)) %>% write.table(file.path(out_dir,"count_matrix.tsv"),quote=FALSE,sep="\t",row.names=FALSE, col.names=TRUE)

## gseLimma.abundance$counts[idx,] %>% as.data.frame %>% mutate(gene.id=rownames(.)) %>% write.table(file.path(out_dir,"TPM_matrix.tsv"),quote=FALSE,sep="\t",row.names=FALSE, col.names=TRUE)



####### voom ##############
out <- file.path(out_dir,"voom.pdf") #### Voom and plot
CairoPDF(width = 11 , height = 7.5, file=out, pointsize=12, bg = "white")
v <- voom(dds,design,plot=TRUE)
dev.off()

fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrasts=contrasts)
efit <- eBayes(fit2)

out <- file.path(out_dir,"pval_histo_contrasts_voom.png")  ##### plot pVal t-tests
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

out.table <- file.path(out_dir,"deg_table_voom.tsv")
write.table(deg.tab,out.table,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

voom_lists_dir <- file.path(out_dir,"voom_lists")
dir.create(voom_lists_dir,showWarnings=FALSE)

for (elt in colnames(contrasts))
    deg.tab %>%
        dplyr::filter(contrast==elt) %>%
        dplyr::filter(adj.P.Val < 0.05) %>%
        dplyr::select(gene.id) %>%
        write.table(paste(file.path(voom_lists_dir,elt),"_voom.list",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

####### voomwithqualityweights ##############



out <- file.path(out_dir,"voomWQW.pdf") #### VoomWQW and plot
CairoPDF(width = 11 , height = 7.5, file=out, pointsize=12, bg = "white")
v <- voomWithQualityWeights(dds,design,plot=TRUE)
dev.off()

fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrasts=contrasts)
efit <- eBayes(fit2)

out <- file.path(out_dir,"pval_histo_contrasts_voomWQW.png")  ##### plot pVal t-tests
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

out.table <- file.path(out_dir,"deg_table_voomWQW_all_batches_with_w10.tsv")
write.table(deg.tab,out.table,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

voomWQW_lists_dir <- file.path(out_dir,"voomWQW_lists_with")
dir.create(voomWQW_lists_dir,showWarnings=FALSE)


for (elt in colnames(contrasts))
    deg.tab %>%
        dplyr::filter(contrast==elt) %>%
        dplyr::filter(adj.P.Val < 0.05) %>%
        dplyr::select(gene.id) %>%
        write.table(paste(file.path(voomWQW_lists_dir,elt),"_voomWQW.list",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)


############ check stuff ###########
deg.table.2 <-read.table(file.path(out_dir,"deg_table_voomWQW_old.tsv"),header=TRUE,sep="\t")
## deg_wo10.table <-read.table(file.path(out_dir,"deg_table_voomWQW_wo10.tsv"),header=TRUE,sep="\t")

joined.tables <- left_join(deg.tab,deg.table.2 ,by=c("gene.id","contrast")) %>%
    dplyr::select(gene.id,logFC.x,logFC.y,adj.P.Val.x,adj.P.Val.y,contrast) %>%
    mutate(diff.logFC=logFC.y-logFC.x)
