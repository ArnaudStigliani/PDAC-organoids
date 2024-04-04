rm(list=ls())        

######################
library(stringr)
library(ggplot2)
library(org.Mm.eg.db)
library(tidyverse)
library(clusterProfiler)
library(maditr)

########################
out_dir <- file.path("results")
dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)
##################

map_table <- read.table("../data/map_batch_1_2_3.txt",sep="\t",header=TRUE) 
coldata <-map_table
rownames(coldata) <- coldata$Samples
coldata$names  <-  coldata$Sample 
coldata$Sample <- coldata$Sample %>% as.numeric
coldata <- coldata %>% filter(!grepl("RasG12D",Cell_line))
coldata <- coldata %>% mutate(Condition=ifelse(coldata$Condition=="bounce","bb" ,coldata$Condition)) %>%
    mutate(Cell_line=str_replace(coldata$Cell_line," ","_")) %>%
    mutate(Cell_line = Cell_line %>%
               str_replace("mP4_p53KO", "p53KO") %>%
               str_replace("mP4_p53R273H", "p53R273H")) %>%
    dplyr::filter(!(Condition == "CTRL" & Week != 1)) %>% 
    mutate(combo = ifelse(Condition != "AA",
                          paste(Condition,Cell_line, sep="_"),
                          paste(paste0(Condition,".", Week), Cell_line, sep="_"))%>%
               str_replace(.,"CTRL","C") %>%
               str_replace(.,"AA","A") %>%
               str_replace(.,"bb","B")  %>%
               str_replace(.,"_$",""))  %>%
    dplyr::filter(Week != 10)
    
combo <- coldata %>%
    dplyr::select(sample.id=Sample, combo) 


########### to read from the Expression matrix instead of tximeta #########
gseLimma.path <- file.path( "data/TPM_batch_corrected.tsv" )
gseLimma.table <- read.table(gseLimma.path, header=TRUE) 

##### open logFC table
DEG.tab.name <- "results/deg_table_voomWQW.tsv"
DEG.tab <-  read.table(DEG.tab.name, header=TRUE,sep="\t") 



#### dir 
Expressed.genes <- gseLimma.table$gene.id 

gene.df <- bitr( Expressed.genes, fromType = "ENSEMBL",
                toType = c("SYMBOL", "ENTREZID"),
                OrgDb = org.Mm.eg.db) %>% dplyr::rename(gene.id=ENSEMBL,gene.name=SYMBOL) %>%
           mutate(gene.name = gene.name) %>% 
           distinct(gene.id, .keep_all = TRUE) 

abc.name <-"data/ABC_all.list"
abc.df <- read.table(abc.name, sep="\t", header=FALSE)
abc.genes <- abc.df[[1]]



gene.expr.df.test.no_filt <- gseLimma.table %>%
    dplyr::filter(gene.id %in% abc.genes, gene.id %in% Expressed.genes) %>%
    melt(id.vars="gene.id", value.name="counts",variable.name="sample.id") %>%
    mutate(sample.id = sample.id %>% str_replace("X","")) %>%
    right_join(combo %>% mutate(sample.id = sample.id %>% as.character)) %>%
    left_join(gene.df) %>% 
    group_by(combo, gene.name) %>%
    summarize(median.counts=median(counts), log.median.counts=log(1 + median.counts, 10)) %>%
    ungroup %>%
    separate(combo, into=c("condi","geno"), sep="_", remove=FALSE) %>%
    mutate(condi = factor(condi, c("C","A.5","A.8","A.11", "B")))



################### plot ABC genes


matrix.pre.plot <- gene.expr.df.test.no_filt  %>% 
    dcast(gene.name ~ combo, value.var = "log.median.counts", fill=0) %>% column_to_rownames("gene.name")
hclust.heatmap <- hclust(dist(matrix.pre.plot))  
levels.y <- hclust.heatmap$labels[hclust.heatmap$order %>% rev]


my_palette <- colorRampPalette(RColorBrewer::brewer.pal(n=8,name="RdYlBu"))(100) %>% rev
col_breaks <- c(seq(0.29, 1,length=25), 
                seq(1.01, 1.5,length=25),
                seq(1.51, 2.2,length=25),
                seq(2.21, 2.62,length=25))



g.svg <- ggplot(data = gene.expr.df.test.no_filt %>% mutate(gene.name = factor(gene.name, levels.y)), colour = "grey",  
               aes(y = gene.name , x = condi, fill = log.median.counts)) +  
    geom_tile() +
    scale_fill_gradientn(colours=my_palette,values=scales::rescale(col_breaks), breaks=c(1, 1.5, 2)) +
    ## scale_colour_gradientn(colours=my_palette,values=scales::rescale(col_breaks), guide=NULL) +
    ylab("") +
    xlab("") +
    theme(text = element_text(size=9)) +
    scale_y_discrete(expand=c(0,0)) +
    theme(panel.border = element_rect(size=0,fill=NA)) +
    theme( legend.text=element_text(size=8),legend.title=element_text(size=8)) +
    theme_classic() +
    theme(legend.key.size = unit(1, "cm"), legend.text=element_text(size=7),legend.title=element_text(size=7)) +
    theme(legend.margin = margin(0,3,0,-4),legend.box.margin = margin(0,3,0,-4)) +
    theme(strip.background=element_blank(),strip.text.x=element_text(size=8)) +
    theme(plot.margin = unit(c(5,-2,0,-10), "pt")) +
    theme(axis.line = element_line(colour = "white", size = NULL)) +
    theme(axis.text.x = element_text(size=8),axis.ticks.x=element_blank()) +
    facet_wrap(~geno, ncol=4)



out_heatmap_small_pdf <- file.path(out_dir,"heatmapt_small_abc.pdf")
ggsave(g.svg, filename = out_heatmap_small_pdf,width = 7, height = 7)
