rm(list=ls())
out_dir <- paste("results/")
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
library(magrittr)
library(gridExtra)
library(cowplot)
library(ggrepel)

##############################

rds.file <- readRDS("../data/20220516_gene_expr_organoids_drug_adaptation.rds")  ## RDS file containing SUmmarizedExperiment with Expression matrix (in TPM) and design matrix


#### Read count matrix/metadata  and filter

coldata <- rds.file %>% colData %>%
    as.data.frame %>%
    rownames_to_column("sampleID") %>% 
    mutate(Condition = Condition %>% str_replace("\\+","")) %>%
    dplyr::filter(Cell_line != "mI1_RasG12D" ) %>%
    mutate(Drug_adaptation=ifelse(Drug_adaptation=="No", "nDA", "DA")) %>% 
    mutate(Batch = as.character(Batch)) 

###################################
###### write metadata for the paper
coldata.sup_table <- rds.file %>%
    colData %>%
    as.data.frame %>%
    rownames_to_column("Library")
coldata.sup_table.batch4 <- coldata.sup_table %>%
    dplyr::filter(Batch==4) %>%
    mutate(libn = Library %>% str_replace("Res","")) %>%
    arrange(libn %>% as.numeric) %>%
    dplyr::select(-libn)
coldata.sup_table.no_batch4  <- coldata.sup_table %>%
    dplyr::filter(Batch!=4) %>%
    arrange(Library %>% as.numeric)
coldata1 <- read.table(file.path(out_dir,"coldata1.tsv"), sep="\t", header=TRUE)

coldata.sup_table.to_write <- left_join( rbind.data.frame(coldata.sup_table.no_batch4, coldata.sup_table.batch4),
                                        coldata1,
                                        by="Library") %>%
    mutate(group = paste(Condition,"week",Week)) %>%
    mutate(group = ifelse(Condition=="bb","AA->pH7.4",group)) %>%
    mutate(Condition = Condition %>% str_replace("EG","E+G"))%>%
    mutate(group=ifelse(grepl("E|G|E\\+G",Condition) ,
                 ifelse(Drug_adaptation=="No",paste0("Acute treatment : ", Condition),
                        paste0("Drug adaptation : ", Condition)),group))  %>%
    dplyr::select( -Condition, -Drug_adaptation)

write.table(coldata.sup_table.to_write, file.path(out_dir, "coldata_Sup_fig.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
###################################

## prepare design
combo <- coldata %$% paste(.$Condition, .$Cell_line, .$Week, .$Drug_adaptation, "batch",.$Batch ,sep="_") %>%
    str_replace(.,"CTRL","C") %>%
    str_replace(.,"AA","A") %>%
    str_replace(.,"bb","B")  %>%
    str_replace(.,"_$","") %>%
    str_replace(.,"mP4_p53","p53") %>%
    str_replace(.,"batch_","batch")
    


design <- model.matrix( ~ 0 + combo + coldata$Batch)
colnames(design) <- gsub("combo","", colnames(design)) %>%  str_replace("coldata\\$Batch","Batch")

##############

gseLimma <-  rds.file %>% assays() %>% .[["counts"]]  %>% .[,colnames(.) %in% coldata$sampleID]%>% DGEList
idx <- filterByExpr(gseLimma  ,min.count=55,min.total.counts=40
                   ,design) # filtering genes according to counts

# Log Expression with 1 as pseudocount
y <- gseLimma$counts
yprime <-log(y+1)
ylog <- removeBatchEffect(yprime,coldata$Batch)


##  lists of samples to study :

mN10_list <- combo %in% (combo[grepl("mN10",combo,fixed=TRUE)]  %>%  .[!grepl("_10_",.)]) %>% which(.)
p53KO_list <- combo %in% (combo[grepl("p53KO",combo,fixed=TRUE)] %>%  .[!grepl("_10_",.)]) %>% which(.)
p53R273H_list <- combo %in% (combo[grepl("p53R273H",combo,fixed=TRUE)] %>%  .[!grepl("_10_",.)]) %>% which(.)
mP4_list <- combo %in% (combo[grepl("mP4",combo,fixed=TRUE)] %>%  .[!grepl("_10_",.)]) %>% which(.)
mP4_Drug_list <- combo %in% (combo[grepl("mP4",combo,fixed=TRUE)] %>%  .[grepl("batch4",.)]) %>% which(.)
mP4_no_Drug_list <- mP4_list[!(mP4_list %in% mP4_Drug_list)]


mN10.df_list  <- coldata[mN10_list,] %>%
    mutate(color=paste(Condition,"week",Week)) %>%
    mutate(color=ifelse(Condition=="bb","AA->pH7.4",color)) %>%
    dplyr::select(sampleID, color) %>%
    mutate(color =   factor(color, levels = color  %>% unique %>% sort  %>% .[c(5,2,3,1,4)]))
p53KO.df_list  <- coldata[p53KO_list,] %>%
    mutate(color=paste(Condition,"week",Week)) %>%
    mutate(color=ifelse(Condition=="bb","AA->pH7.4",color)) %>%
    dplyr::select(sampleID, color) %>%
    mutate(color =   factor(color, levels = color  %>% unique %>% sort  %>% .[c(5,7,8,6,2,3,1,4)]))
p53R273H.df_list  <- coldata[p53R273H_list,] %>%
    mutate(color=paste(Condition,"week",Week)) %>%
    mutate(color=ifelse(Condition=="bb","AA->pH7.4",color)) %>%
    dplyr::select(sampleID, color)  %>%
    mutate(color =   factor(color, levels = color  %>% unique %>% sort  %>% .[c(5,2,3,1,4)]))
mP4.df_list  <- coldata[mP4_list,] %>%
    mutate(color = paste(Condition,"week",Week)) %>%
    mutate(color = ifelse(Condition=="bb","AA->pH7.4",color)) %>%
    mutate(Condition = Condition %>% str_replace("EG","E+G"))%>%
    mutate(color=ifelse(grepl("E|G|E\\+G",Condition) ,
                 ifelse(Drug_adaptation=="nDA",paste0("Acute treatment : ", Condition),
                        paste0("Drug adaptation : ", Condition)),color))  %>%
    mutate(color =   factor(color, levels = color  %>% unique %>% sort  %>% .[c(8,11,12,9,10,2,3,1,4,5,7,6,13,15,14)]))

mP4_noDrug.df_list  <- coldata[mP4_no_Drug_list,] %>%
    mutate(color = paste(Condition,"week",Week)) %>%
    mutate(color = ifelse(Condition=="bb","AA->pH7.4",color)) %>%
    mutate(Condition = Condition %>% str_replace("EG","E+G"))%>%
    mutate(color=ifelse(grepl("E|G|E\\+G",Condition) ,
                 ifelse(Drug_adaptation=="nDA",paste0("Acute treatment : ", Condition),
                        paste0("Drug adaptation : ", Condition)),color)) %>%
    mutate(color =   factor(color, levels = color  %>% unique %>% sort  %>% .[c(5,7,8,6,2,3,1,4)]))

mP4_Drug.df_list  <- coldata[mP4_Drug_list,] %>%
    mutate(color = paste(Condition,"week",Week)) %>%
    mutate(color = ifelse(Condition=="bb","AA->pH7.4",color)) %>%
    mutate(Condition = Condition %>% str_replace("EG","E+G"))%>%
    mutate(color=ifelse(grepl("E|G|E\\+G",Condition) ,
                 ifelse(Drug_adaptation=="nDA",paste0("Acute treatment : ", Condition),
                        paste0("Drug adaptation : ", Condition)),color)) %>%
    mutate(color =   factor(color, levels = color  %>% unique %>% sort  %>% .[c(4,5,1,3,2,6,8,7)]))


lists <- list(p53KO_list, mN10_list, p53R273H_list,mP4_no_Drug_list,mP4_Drug_list)
lists.df <- list( p53KO.df_list, mN10.df_list,  p53R273H.df_list,mP4_noDrug.df_list,mP4_Drug.df_list)
names(lists) <- c("p53KO", "mN10", "R273H", "mP4_Acid_adaptation",  "mP4_Drug_adaptation")

### defining function for legend

get_legend <- function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

## ###### vst #################
all_plots <- list() ## to store the plots
for (i in  1:length(lists))
{
    name_elt <- names(lists[i])
    elt <- lists[[i]]
    out <- file.path(out_dir,paste(name_elt,"_PCA.png",sep=""))
    zylog <- ylog[,elt]
    pcs = prcomp(t(zylog))
    axis_lab <- paste(as.character(colnames(summary(pcs)$importance))," ( ",as.character(round(summary(pcs)$importance[2,]*100)),"% )",sep="")
    tab <- as.data.frame(pcs$x) %>%
        rownames_to_column("sampleID") %>%
        left_join(lists.df[[i]]) 
    ## PCA Plot
    Cairo(width = 1500, height = 700, file=out, type="png", pointsize=12, bg = "white", canvas = "white", units = "px", dpi = "auto")
    p1 <- ggplot(tab) +
        theme_classic() +
        labs(color="Samples") +
        scale_colour_brewer(palette = "Set3") +
        geom_point(size=3,aes(x=PC1,y=PC2,color=color)) +
        ## geom_text(aes(x=PC1,y=PC2,label=sampleID),hjust=0, vjust=0 , fontface="bold") +
        geom_text_repel(aes(label=sampleID, x=PC1, y=PC2),size=2) +
        xlab(axis_lab[1]) +
        ylab(axis_lab[2]) +
        theme(text = element_text(size=9)) +
        theme(axis.ticks = element_line(size=0.1),axis.ticks.length = unit(0.04,"cm")) +
        theme(panel.border = element_blank()) +
        theme( legend.position="none") +
        ggtitle(name_elt)
    p2 <- ggplot(tab) +
        theme_classic() +
        labs(color="Samples") +
        scale_colour_brewer(palette = "Set3") +
        geom_point(size=3,aes(x=PC3,y=PC4,color=color)) +
        ## geom_text(aes(x=PC3,y=PC4,label=sampleID),hjust=0, vjust=0 , fontface="bold") +
        geom_text_repel(aes(label=sampleID, x=PC3, y=PC4), size=2) +
        xlab(axis_lab[3]) +
        ylab(axis_lab[4]) +
        theme(text = element_text(size=9)) +
        theme(axis.ticks = element_line(size=0.1),axis.ticks.length = unit(0.04,"cm")) +
        theme(panel.border = element_blank()) +
        theme( legend.text=element_text(size=8),legend.title=element_text(size=8)) +
        ggtitle("")
    legend <- get_legend(p2)
    p2 <- p2 + theme(legend.position="none")
    all_plots[[name_elt]] <- grid.arrange(p1, p2, legend , nrow = 1, ncol = 3,widths=c(1,1,0.3))
    grid.arrange(p1, p2, legend , nrow = 1, ncol = 3,widths=c(1,1,0.3)) 
    dev.off()
}



p <- plot_grid(all_plots[[1]],
               all_plots[[2]],
               all_plots[[3]],
               all_plots[[4]],
               all_plots[[5]],
               ncol=1)
all_plot.name <- file.path(out_dir,"PCA_all_geno.pdf")
ggpubr::ggexport(p, filename = all_plot.name,width = 8.3, height = 16)
