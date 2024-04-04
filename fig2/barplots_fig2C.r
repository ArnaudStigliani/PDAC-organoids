
## rm(list=ls())
args = c("5","8","11","B")
out_dir <- file.path("./resutls/")
dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)

#### load libraries ##########
library(stringr)
library(ggplot2)
library(ggplus)
library(Cairo)
library(reshape2)
library(pheatmap)
library(ggVennDiagram)
library(ComplexUpset)
library(gprofiler2)
library(org.Mm.eg.db)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(parallel)
library(tidyverse)
library(Matrix)
library(ggtext)

##########################


DEG.tab.name <- "../data/deg_table_voomWQW.tsv"

DEG.tab <-  read.table(DEG.tab.name, header=TRUE,sep="\t")

genes.up <- list()
genes.down <- list()
genes.all <- list()
for(week.args in args)
{
    DEG.tab.week.args  <-  DEG.tab %>%
        mutate(contrast=str_replace(DEG.tab$contrast,"mP4_p","mP4:p")) %>%
        separate(contrast,into=c("geno","contrast"),sep="_") %>%
        mutate(geno=str_replace(.$geno,"mP4:p","mP4_p")) %>%
        mutate(contrast.geno=paste(contrast,geno,sep=".")) %>%
        dplyr::filter(contrast==week.args) 
####################
    DEG.tab.dcast <-  DEG.tab.week.args %>%  mutate(logFC = ifelse(adj.P.Val < 0.05,logFC,0)) %>%
        mutate(logFC = ifelse(abs(logFC) > 0.5  ,logFC,0)) %>% 
        dcast(gene.id ~ geno,value.var="logFC")

##### making lists

    DEG.tab.dcast.logFC <- DEG.tab.dcast[,c(2:5)]
    rownames(DEG.tab.dcast.logFC) <- DEG.tab.dcast %>% .$gene.id

    genes.up[[week.args]] <- data.table::between(DEG.tab.dcast.logFC, 0, Inf,incbounds=FALSE)
    genes.down[[week.args]] <- data.table::between(DEG.tab.dcast.logFC, -Inf, -0,incbounds=FALSE)
    genes.all[[week.args]] <- rbind(genes.down[[week.args]],genes.up[[week.args]])
    colnames(genes.up[[week.args]]) <- paste(colnames(genes.up[[week.args]]),week.args,sep=".")
    colnames(genes.down[[week.args]]) <- paste(colnames(genes.down[[week.args]]),week.args,sep=".")
    colnames(genes.all[[week.args]]) <- paste(colnames(genes.all[[week.args]]),week.args,sep=".")
}

colnames.mats <- lapply(genes.up,colnames) %>% unlist %>% str_replace("_"," ")
genes.up.mat <- bdiag(genes.up) %>%  as.matrix %>% as.data.frame %>% setNames(colnames.mats)
genes.down.mat <- bdiag(genes.down) %>%  as.matrix %>% as.data.frame %>% setNames(colnames.mats)
genes.all.mat <- bdiag(genes.all) %>%  as.matrix %>% as.data.frame %>% setNames(colnames.mats)


########### preparing categories
b <- 1:4 %>% lapply(combn,x=4) %>% lapply(as.data.frame) %>% lapply(as.list)

c <- list()
for (elt in b)
{
    for (elt2 in elt)
    {
        c[[(length(c)+1)]] <- elt2
    }
}

a <- list()
j <- 0
for (i in c(0,4,8,12))
{
    for (elt in c)
    {
        j <- j+1
        a[[j]] <- elt + i 
    }
}



##### exclusive lists
lists.up <- list()
for (elt in a)
{
    lists.up[[colnames(genes.up.mat)[elt] %>% paste(collapse=":")]] <- rownames(genes.up.mat)[rowSums(genes.up.mat[,elt] %>% as.matrix) == rowSums(genes.up.mat) & rowSums(genes.up.mat[,elt] %>% as.matrix) == length(elt)]
}

lists.down <- list()
for (elt in a)
{
    lists.down[[colnames(genes.down.mat)[elt] %>% paste(collapse=":")]] <- rownames(genes.down.mat)[rowSums(genes.down.mat[,elt] %>% as.matrix) == rowSums(genes.down.mat) & rowSums(genes.down.mat[,elt] %>% as.matrix)== length(elt)]
}

lists.all <- list()
for (elt in a)
{
    lists.all[[colnames(genes.all.mat)[elt] %>% paste(collapse=":")]] <- rownames(genes.all.mat)[rowSums(genes.all.mat[,elt] %>% as.matrix) == rowSums(genes.all.mat) & rowSums(genes.all.mat[,elt] %>% as.matrix)== length(elt)]
}


### prepare format list

pattern.mat <- matrix(0,nrow=length(lists.up),ncol=16) 
colnames(pattern.mat) <- colnames(genes.up.mat)
rownames(pattern.mat) <- names(lists.up)

for (i in 1:length(lists.up))
{
    pattern.mat[i,a[[i]]] <- 1
}


palette <- c("white", colorRampPalette(RColorBrewer::brewer.pal(n=9,"Blues")[5:8])(100)[c(25,50,75,100)]) %>% as.data.frame
palette <- palette %>% mutate(col.levels=c("0","5","8","11","B")) %>%
    dplyr::rename(col=".")


palette <- c("white", RColorBrewer::brewer.pal(n=5,"Blues")[-1]) %>% as.data.frame
palette <- palette %>% mutate(col.levels=c("0","5","8","11","B")) %>%
    dplyr::rename(col=".")

## palette <- c("white",RColorBrewer::brewer.pal(n=4,"Set2")) %>% as.data.frame
## palette <- palette %>% mutate(col.levels=c("0","5","8","11","B")) %>%
##     dplyr::rename(col=".")


pattern.df <- pattern.mat %>% as.data.frame  %>% rownames_to_column("group") %>%
    mutate(group = group %>% str_replace_all(":",".") %>% str_replace_all(" ",".")) %>% 
    mutate(group=factor(group,levels=group)) %>%
    melt(variable.name="geno.contrast") %>%
    mutate(geno.contrast=factor(geno.contrast,levels=colnames(pattern.mat) %>% rev)) %>%
    mutate(contrast=str_replace(.$group,".*\\.","")) %>%
    mutate(col.levels=ifelse(value==1,contrast,"0")) %>%
    mutate(col.levels=factor(col.levels,levels=c("0","5","8","11","B"))) %>%
    left_join(palette,by="col.levels") %>%
    mutate(geno=str_replace(.$geno.contrast,"\\..*","")) %>%
    mutate(geno=factor(geno,levels=c("mN10","mP4", "mP4 p53KO", "mP4 p53R273H") %>% rev )) %>%
    dplyr::filter(value==1)

levels.inter <- matrix(pattern.df$group %>% levels,ncol=4) %>% t %>% c  

pattern.df <- pattern.df %>%
    dplyr::mutate(group=factor(group,levels=levels.inter))

####### format lists
lists.up.df <- lists.up %>% 
    lapply(length) %>%
    as.data.frame %>%
    t %>%
    as.data.frame %>% 
    rownames_to_column("group") %>%
    dplyr::select(counts.up=V1,group) %>%
    mutate(group=factor(group,levels=group)) %>%
    mutate(group=factor(group,levels=matrix(group %>% levels,ncol=4) %>% t %>% c  )) %>%
    mutate(col.levels=str_replace(.$group,".*\\.","")) %>%
    left_join(palette,by="col.levels")

lists.down.df <- lists.down %>% 
    lapply(length) %>%
    as.data.frame %>%
    t %>%
    as.data.frame %>% 
    rownames_to_column("group") %>%
    dplyr::select(counts.down=V1,group) %>%
    mutate(group=factor(group,levels=group)) %>%
    mutate(group=factor(group,levels=matrix(group %>% levels,ncol=4) %>% t %>% c  )) %>%
    mutate(col.levels=str_replace(.$group,".*\\.","")) %>%
    left_join(palette,by="col.levels")

by.join <- colnames(lists.down.df)[-1]

lists.concat.df <- full_join(lists.down.df,lists.up.df,by=by.join) %>%
    dplyr::filter(counts.up > 100, counts.down > 100)

#### prepare ggplot objects
pattern.plot <- ggplot(pattern.df %>% dplyr::filter(group %in% lists.concat.df$group)) +
    geom_point(aes(x=group,y=geno,color=col),size=1.5) +
    geom_line( aes(x=group,y=geno,group=group,color=col),size=0.5) +
    scale_color_identity() +
    theme_bw() +
    theme(plot.margin = unit(c(0,10,0,0), "pt"), axis.title=element_blank(),axis.text.y=element_text(size=8)) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
    theme( axis.text.x = element_blank(),axis.ticks.x=element_blank())
out_pdf <- file.path(out_dir,"pattern.pdf")
ggpubr::ggexport(pattern.plot, filename = out_pdf,height=1,width=7) 


bar.up <- ggplot(lists.concat.df) +
    geom_col(aes(x=group,y=counts.up,fill=col)) +
    scale_fill_identity()  +
    theme_bw() +    
    ylab("Up regulated \n genes") + 
    theme(plot.margin = unit(c(5,10,0,0), "pt"),axis.title.x=element_blank(), axis.title.y=element_text(size=8)) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
    theme( axis.text.x = element_blank(),axis.ticks.x=element_blank()) +
    ylim(c(0,max(lists.up.df$counts) + 150)) +
    geom_text(size=2.2, aes(y=counts.up,x=group,label=counts.up),vjust=-0.5)

out_pdf <- file.path(out_dir,"bar_up.pdf")
ggpubr::ggexport(bar.up, filename = out_pdf,height=1,width=7)


bar.down <- ggplot(lists.concat.df) +
    geom_col(aes(x=group,y=counts.down,fill=col)) +
    scale_fill_identity()  +
    theme_bw() +    
    ylab("Down regulated \n genes") + 
    theme(plot.margin = unit(c(5,10,0,0), "pt"),axis.title.x=element_blank(), axis.title.y=element_text(size=8)) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
    theme( axis.text.x = element_blank(),axis.ticks.x=element_blank()) +
    ylim(c(0,max(lists.down.df$counts) + 150)) +
    geom_text(size=2.2, aes(y=counts.down,x=group,label=counts.down),vjust=-0.5)



lists.down.df$col.levels[lists.down.df$col.levels=="B"] <- "AA->7.4"
lists.down.df$col.levels[lists.down.df$col.levels=="5"] <- "AAweek5"
lists.down.df$col.levels[lists.down.df$col.levels=="8"] <- "AAweek8"
lists.down.df$col.levels[lists.down.df$col.levels=="11"] <- "AAweek11"

cols.names <- lists.down.df$col.levels %>% unique
cols <- lists.down.df$col %>% unique
names(cols) <- cols.names

pre.legend <- ggplot(lists.down.df %>% mutate(col.levels=factor(col.levels,levels=names(cols)))) +
    geom_col(aes(x=group,y=counts.down,fill=col.levels)) +
    scale_fill_manual(values=cols) +
    guides(fill=guide_legend(title="Time\npoint")) +
    theme(plot.margin = unit(c(0,0,0,2), "pt")) +
    theme(legend.margin = margin(0,0,0,0),legend.box.margin = margin(0,0,0,0)) + 
    theme(legend.key.size = unit(0.40, "cm"), legend.text=element_text(size=8),legend.title=element_text(size=8))

legend <- as_ggplot(get_legend(pre.legend))


out_pdf <- file.path(out_dir,"legend.pdf")
ggpubr::ggexport(legend, filename = out_pdf,height=10,width=15)


## plot stuff 
gg.up <- ggdraw() +
    draw_plot(bar.up , x = 0 , y = 0.6 , width = 1, heigh = 0.4 ) +
    draw_plot(pattern.plot, x = 0 , y = 0 , width = 1, heigh = 0.6 ) 

gg.down <- ggdraw() +
    draw_plot(bar.down , x = 0 , y = 0.6 , width = 1, heigh = 0.4 ) +
    draw_plot(pattern.plot, x = 0 , y = 0 , width = 1, heigh = 0.6 ) 

out_pdf <- file.path(out_dir,"Upset_up.pdf")
ggpubr::ggexport(gg.up, filename = out_pdf,height= 5 ,width=15)

out_pdf <- file.path(out_dir,"Upset_down.pdf")
ggpubr::ggexport(gg.down, filename = out_pdf,height= 5 ,width=15)

gg.combo <- ggdraw() +
    draw_plot(bar.up , x = 0.0355 , y = 0.6 , width = 0.88, heigh = 0.4 ) +
    draw_plot(legend , x = 0.94 , y = 0.4 , width = 0.016, heigh = 0.15 ) +
    draw_plot(bar.down , x = 0.0355 , y = 0.2 , width = 0.88, heigh = 0.4 ) +
    draw_plot(pattern.plot, x = 0.010 , y = 0 , width = 0.905, heigh = 0.2 ) 

out_pdf <- file.path(out_dir,"Upset_all.pdf")
ggpubr::ggexport(gg.combo, filename = out_pdf,height= 5  ,width=10)


