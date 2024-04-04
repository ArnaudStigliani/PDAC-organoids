## rm(list=ls())
out_dir <- paste("/isdata/sandelin/people/arnaud/pH_organoids/results_2/gene_expression_limma_time_series_BBvsW1/figure2/heatmap_4")
Robj_dir <- paste("/isdata/sandelin/people/arnaud/pH_organoids/results_2/gene_expression_limma_time_series_BBvsW1/figure2/Robj")

dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(Robj_dir,showWarnings=FALSE,recursive=TRUE)


#### load libraries ##########

library(stringr)
library(ggplot2)
library(ggplus)
library(Cairo)
library(reshape2)
library(tidyverse)
library(pheatmap)
library(plyr)
library(cowplot)
##############################
## read tables ##


table_voomWQW.name <- file.path("/isdata/sandelin/people/arnaud/pH_organoids/results_2","gene_expression_limma_time_series_BBvsW1","deg_table_voomWQW.tsv")
table_voomWQW <- read.table(table_voomWQW.name,header=TRUE,sep="\t")

##############################
## heatmap voomWQW  ##
table_voomWQW.heatmap <- table_voomWQW %>% dplyr::select(gene.id,logFC,contrast) %>%
    mutate(contrast=str_replace(table_voomWQW$contrast,"mP4_p","mP4:p")) %>%
    separate(contrast,into=c("geno","contrast"),sep="_") %>%
    mutate(geno=str_replace(.$geno,"mP4:p","mP4_p")) %>%
    mutate(geno.contrast=paste(geno,contrast,sep=".")) %>%
    dplyr::select(-contrast,-geno) %>%
    dcast(gene.id ~ geno.contrast,value.var="logFC") ### wide format to cluster the rows later

rownames(table_voomWQW.heatmap) <- table_voomWQW.heatmap$gene.id
table_voomWQW.heatmap <- table_voomWQW.heatmap %>% dplyr::select(-gene.id)


######### color ########
### ggplot heatmap
my_palette <- colorRampPalette(RColorBrewer::brewer.pal(n=5,name="RdBu"))(100) %>% rev
col_breaks <- c(seq(-3,-0.5,length=25),
                seq(-0.49,-0.01,length=25),
                seq(0.01,0.49,length=25),
                seq(0.5,3,length=25))



matrix.pre.plot <- table_voomWQW.heatmap %>% filter(((abs(.) > 0.5)  %>% rowSums) > 1) ## filter genes that don't have |logFC| > 0.5 in at least 2 conditions
hclust.heatmap <- hclust(dist(matrix.pre.plot))  
levels.y <- hclust.heatmap$labels[hclust.heatmap$order %>% rev] ## get order of hclust to order the rows of the heatmap




heatmap.pre.plot <- table_voomWQW.heatmap %>% filter(((abs(.) > 0.5)  %>% rowSums) > 1) %>%
    rownames_to_column("gene.id") %>% 
    melt(varnames=c(rownames(.) %>% .[-1]),variable.name="geno.time",value.name="logFC") %>%
    separate(geno.time,c("geno","time"),sep="\\.") %>%
    mutate(time=ifelse(time=="B","R",time)) %>%
    mutate(time=factor(time,levels=c("5","8","11","R"))) %>%
    mutate(geno=str_replace(.$geno,"_","\n")) %>%
    mutate(gene.id=factor(gene.id,levels=levels.y)) %>%  ## use th orders of hclust to order the rows  (by assigning factor levels))
    mutate(logFC=ifelse(abs(logFC)<3,logFC,3*sign(logFC))) 


g.svg <- ggplot(data = heatmap.pre.plot, 
  , aes(y = gene.id , x = time , color = logFC,fill=logFC)) +  
    geom_tile() +
    scale_fill_gradientn(colours=my_palette,values=scales::rescale(col_breaks)) +
    scale_colour_gradientn(colours=my_palette,values=scales::rescale(col_breaks)) +
    ylab("") +
    xlab("") +
    theme(text = element_text(size=9)) +
    scale_x_discrete(drop=TRUE) +
    theme(panel.border = element_rect(size=0,fill=NA)) +
    labs(size="-log(q.val)")  + 
    theme( legend.text=element_text(size=8),legend.title=element_text(size=8)) +
    facet_grid(cols=vars(geno),scales="free_x") +
    theme_classic() +
    theme(legend.key.size = unit(0.35, "cm"), legend.text=element_text(size=7),legend.title=element_text(size=7)) +
    theme(legend.margin = margin(0,3,0,-4),legend.box.margin = margin(0,3,0,-4)) +
    theme(strip.background=element_blank(),strip.text.x=element_text(size=8)) +
    theme(plot.margin = unit(c(5,-2,1,-10), "pt")) +
    theme(axis.line = element_line(colour = "white", size = NULL)) +
    theme(axis.text.x = element_text(size=8),axis.ticks.x=element_blank()) +
    theme(axis.text.y = element_text(size=0),axis.ticks.y=element_blank()) 


saveRDS(g.svg,file=file.path(Robj_dir,"heatmap.rdata"))
out_heatmap_small_svg <- file.path(out_dir,"heatmapt_small.svg")
out_heatmap_small_pdf <- file.path(out_dir,"heatmapt_small.pdf")
ggsave(g.svg, filename = out_heatmap_small_svg,width = 7, height = 7)
ggsave(g.svg, filename = out_heatmap_small_pdf,width = 7, height = 7)

