rm(list=ls())
out_dir <- file.path("./results")
Robj_dir <- paste("./results/Robj")

dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(Robj_dir,showWarnings=FALSE,recursive=TRUE)

#### load libraries ##########
library(tidyverse)
library(reshape2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(parallel)
library(stringr)
library(ggplot2)
library(ggplus)
library(Cairo)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(msigdbr)
library(svglite)
#########################


DEG.tab.name <- "../data/deg_table_voomWQW.tsv"
DEG.tab <-  read.table(DEG.tab.name, header=TRUE,sep="\t")

DEG.tab.week.args  <-  DEG.tab %>%
    mutate(contrast=str_replace(DEG.tab$contrast,"mP4_p","mP4:p")) %>%
    separate(contrast,into=c("geno","contrast"),sep="_") %>%
    mutate(geno=str_replace(.$geno,"mP4:p","mP4_p")) %>%
    mutate(contrast.geno=paste(contrast,geno,sep=".")) %>% 
    mutate(geno.contrast=paste(geno,contrast,sep="."))

Expressed.genes <- DEG.tab.week.args$gene.id %>% unique
gene.df <- bitr( Expressed.genes, fromType = "ENSEMBL",
                toType = c("ENTREZID" ,"SYMBOL"),
                OrgDb = org.Mm.eg.db) %>% dplyr::rename(gene.id=ENSEMBL)




DEG.tab.week.args.names <- DEG.tab.week.args %>%
    left_join(gene.df) %>%
    drop_na() %>%
    distinct(gene.id,geno.contrast,.keep_all=TRUE)
    
    

list.contrasts <- list()
for (group in DEG.tab.week.args.names$geno.contrast %>% unique)
{
    list.contrasts[[group]] <- DEG.tab.week.args.names %>%
        dplyr::filter(geno.contrast==group) %>%
        arrange(desc(logFC)) %>%
        .$logFC
    names(list.contrasts[[group]]) <- DEG.tab.week.args.names %>%
        dplyr::filter(geno.contrast==group) %>%
        arrange(desc(logFC)) %>%
        .$gene.id
}

list.contrasts.kegg <- list()
for (group in DEG.tab.week.args.names$geno.contrast %>% unique)
{
    list.contrasts.kegg[[group]] <- DEG.tab.week.args.names %>%
        dplyr::filter(geno.contrast==group) %>%
        arrange(desc(logFC)) %>%
        .$logFC
    names(list.contrasts.kegg[[group]]) <- DEG.tab.week.args.names %>%
        dplyr::filter(geno.contrast==group) %>%
        arrange(desc(logFC)) %>%
        .$ENTREZID
}



gse.list <- mclapply(list.contrasts,
                     gseGO,
                     ont ="ALL",
                     keyType = "ENSEMBL",
                     minGSSize = 3,
                     maxGSSize = 800,
                     pvalueCutoff = 1,
                     verbose = TRUE,
                     eps=0,
                     OrgDb = org.Mm.eg.db,
                     mc.cores = 16)


gsekegg.list <- mclapply(list.contrasts.kegg,
                       gseKEGG,
                       organism='mmu',
                       minGSSize = 3,
                       maxGSSize = 800,
                       pvalueCutoff = 1,
                       eps=0,
                       verbose = TRUE  ,
                       mc.cores=16)


gse.list.df <- gse.list  %>%  mapply(mutate,.,group=names(.),SIMPLIFY=FALSE)  %>%
    lapply(as_tibble) %>% 
    do.call(what=rbind.data.frame,arg=.) %>%
    mutate(group=factor(group,levels=str_sort(group %>% unique, numeric=TRUE))) %>% 
    dplyr::select( pvalue, p.adjust, ID, Description, setSize, NES, enrichmentScore, group) %>%
    mutate(NES=ifelse(abs(NES) > 2, sign(NES) * 2, NES)) 


gsekegg.list.df <- gsekegg.list  %>%  mapply(mutate,.,group=names(.),SIMPLIFY=FALSE)  %>%
    lapply(as_tibble) %>% 
    do.call(what=rbind.data.frame,arg=.) %>%
    mutate(group=factor(group,levels=str_sort(group %>% unique, numeric=TRUE))) %>% 
    dplyr::select( pvalue, p.adjust, ID, Description, setSize, NES, enrichmentScore, group) %>%
    mutate(NES=ifelse(abs(NES) > 2, sign(NES) * 2, NES)) 



###### MSIG and df drug resistance genes analysis
df.newlist.name <- file.path("../../data/","drug_resistance_genes_carnonical.mouse.rds")
df.newlist <- rbind.data.frame(readRDS(df.newlist.name) %>%
                                  dplyr::mutate(gs_description=paste("df",Drug,Directionality,sep=".")),
                                  readRDS(df.newlist.name) %>%
                                  dplyr::mutate(gs_description="df.all"),
                                  readRDS(df.newlist.name) %>%
                                  dplyr::mutate(gs_description=paste("df",Directionality,sep="."))) %>% 
    dplyr::select(gs_description,ensembl_gene=gene_id)



## df.genes.df.name <- "../data/drug_resistance_info.mouse.csv"
## df.genes.df <- read.table(df.genes.df.name,header=TRUE,sep=",") %>%
##     dplyr::mutate(gs_description=paste("df",Drug,Directionality,sep=".")) %>%
##     dplyr::select(gs_description,ensembl_gene=gene_id)

gene2msig.prefilt <- msigdbr(species = "Mus musculus", category = "C2") %>%
    dplyr::filter(gs_description!="")

gene2msig <- gene2msig.prefilt %>% 
    dplyr::select(gs_description, ensembl_gene) %>%
    rbind.data.frame(df.newlist)




gsemsig.list <- mclapply(list.contrasts,
                     GSEA,
                     TERM2GENE=gene2msig,
                     minGSSize = 3,
                     maxGSSize = 800,
                     pvalueCutoff = 1,
                     verbose = TRUE,
                     eps=0,
                     mc.cores=16)

gsemsig.list.df <- gsemsig.list  %>%  mapply(mutate,.,group=names(.),SIMPLIFY=FALSE)  %>%
    lapply(as_tibble) %>% 
    do.call(what=rbind.data.frame,arg=.) %>%
    mutate(group=factor(group,levels=str_sort(group %>% unique, numeric=TRUE))) %>%
    ## mutate(p.adjust=ifelse(grepl("df",Description),pvalue,p.adjust)) %>%
    dplyr::select( pvalue, p.adjust, ID, Description, setSize, NES, enrichmentScore, group) %>%
    mutate(NES=ifelse(abs(NES) > 2, sign(NES) * 2, NES)) %>%
    mutate(pvalue = 10**(floor(abs(log(pvalue,10)) * sign(log(pvalue,10)))))  
gsemsig.list.df[grepl("df",gsemsig.list.df$Description),]$p.adjust <- p.adjust(gsemsig.list.df[grepl("df",gsemsig.list.df$Description),]$pvalue,method="fdr") 
gsemsig.list.df <- gsemsig.list.df %>%  mutate( Description = ifelse(Description=="df.Erlotinib.resistance", "Erlotinib resistance",Description)) %>%
    dplyr::filter(! grepl("df",Description))


set.names <- c( "Genes up-regulated in MCF7 cells (breast cancer) at 24 h of",
               "Genes up-regulated by hypoxia in TRAMP-C cells (prostatic ca",
               "Genes down-regulated in MCF7 cells (breast cancer) after kno",
               ## "Genes down-regulated in H1975 cells (non-small cell lung can",
               "The \'Cervical Cancer Proliferation Cluster\' (CCPC): genes wh")

reset.names <- c("Up-regulated in breast cancer after estradiol treatment",
                "Up-regulated by hypoxia in prostatic cancer",
                "Down-regulated in breast cancer after HIF1A KO",
                ## "Down-regulated in lung cancer after gefitinib treatment",
                "Cervical Cancer Proliferation Cluster ")


tab <- NULL
for (i in (1:length(set.names)))
{
    gsemsig.list.df %>% dplyr::filter(grepl(set.names[i],Description,fixed=TRUE)) %>% dim %>% print
    print(set.names[i])
    tab <- rbind(tab, gsemsig.list.df %>% dplyr::filter(grepl(set.names[i],Description,fixed=TRUE)) %>%
                      mutate(Description.2=reset.names[i]))
}
tab.2 <- tab %>%
    mutate(Description.2= ifelse( ! grepl("HIF2A",Description),Description.2,"Down-regulated in breast cancer after HIF1A and HIF2A KO"))  %>%
    dplyr::filter(!grepl("6h",Description)) %>%
    dplyr::filter(!grepl("FOXA2 ",Description))  %>%
    dplyr::select(-Description) %>% 
    dplyr::rename(Description=Description.2)  %>%
    bind_rows(gsemsig.list.df %>% dplyr::filter(grepl("Erlotinib",Description)))



##### Gene sets informations (paper writing)
gsemsig.list.df.info <- gsemsig.list  %>%  mapply(mutate,.,group=names(.),SIMPLIFY=FALSE)  %>%
    lapply(as_tibble) %>% 
    do.call(what=rbind.data.frame,arg=.)
tab.info <- NULL
for (i in (1:length(set.names)))
{
    tab.info <- rbind(tab.info, gsemsig.list.df.info %>% dplyr::filter(grepl(set.names[i],Description,fixed=TRUE)) %>%
                      mutate(Description.2=reset.names[i]))
}
tab.df <- tab.info %>%
    mutate(Description.2= ifelse( ! grepl("HIF2A",Description),Description.2,"Down-regulated in breast cancer after HIF1A and HIF2A KO"))  %>%
    dplyr::filter(!grepl("6h",Description)) %>%
    dplyr::filter(!grepl("FOXA2 ",Description))   %>%
    distinct(Description,.keep_all=TRUE)  %>%
    dplyr::select( Description, Description.2)
write.table(tab.df,file.path(out_dir,"cancer_related_terms.tsv"),sep="\t",row.names=FALSE, col.names=FALSE)
#############################################



gse.list.df.all <- rbind(gse.list.df,gsekegg.list.df) %>% mutate(GO_class="") %>%
    rbind(.,tab.2 %>% mutate(GO_class="Cancer related")) %>%
    separate(group , into=c("geno","time"),sep="\\.",remove=FALSE) %>%
    mutate(geno=str_replace(.$geno,"_","\n"))





####### to get the real names
gse.list.df.all  %>%  dplyr::filter( grepl("phosphorylation",Description,ignore.case=TRUE )) %>% arrange(pvalue) %>% .$Description %>% head   # %>% unique
gse.list.df.all  %>%  dplyr::filter( grepl("autophagy",Description,ignore.case=TRUE )) %>% arrange(pvalue) %>% .$Description %>% head   # %>% unique


gseall.list.df.to_plot <- gse.list.df.all %>% dplyr::filter(Description=="positive regulation of cell cycle" |
                                                            Description=="cervical cancer proliferation cluser" |
                                                            Description=="mitotic cell cycle" |
                                                            Description=="DNA replication" |
                                                            ## Description=="ribosome biogenesis" |
                                                            Description=="positive regulation of protein phosphorylation" |
                                                            Description=="chemotaxis" |
                                                            Description=="cell adhesion" |
                                                            Description=="extracellular space" |
                                                            Description=="extracellular matrix" |
                                                            Description=="regulation of locomotion" |
                                                            Description=="cell surface" |
                                                            Description=="regulation of cell motility" |
                                                            Description=="epithelial to mesenchymal transition" |
                                                            Description=="circulatory system development" |
                                                            Description=="Down-regulated in breast cancer after HIf1A K0" |
                                                            Description=="Down-regulated in breast cancer after HIF1A and HIF2A KO" |
                                                            Description=="Upregulated by hypoxia in prostate cancer" |
                                                            Description=="angiogenesis" |
                                                            Description=="response to hypoxia" |
                                                            Description=="vascular process in circulatory system" |
                                                            Description=="immune response" |
                                                            Description=="inflammatory response" |
                                                            Description=="regulation of cytokine production" |
                                                            Description=="secretion by cell"|
                                                            Description=="cilium assembly" |
                                                            Description=="cilium movement" |
                                                            Description=="axoneme" |
                                                            Description=="ciliary plasm" |
                                                            Description=="ciliary basal body" |
                                                            Description=="microtubule bundle formation" |
                                                            Description=="mitochondrion organization" |
                                                            Description=="Carbon metabolism" |
                                                            ## Description=="positive regulation of phosphate metabolic process" |
                                                            Description=="Platinum drug resistance" |
                                                            Description=="EGFR tyrosine kinase inhibitor resistance" |
                                                            Description=="Up-regulated in breast cancer after estradiol treatment" |
                                                            Description=="glutathione metabolic process" |
                                                            Description=="autophagy" |
                                                            Description=="Cellular senescence - Mus musculus (house mouse)" |
                                                            Description=="Metabolism of xenobiotics by cytochrome P450" |
                                                            GO_class=='Cancer related') %>%
    mutate(Description= Description %>% str_replace(" - Mus musculus.*", "")) %>% 
    ## mutate(pvalue.adj = p.adjust(pvalue, method="BH"))  %>%
    mutate(pvalue.adj = p.adjust)  %>%
    dplyr::filter(pvalue.adj<0.05) %>% 
    mutate(pvalue.adj = ifelse(pvalue.adj< 10^-10,10^-10,pvalue.adj))  %>%
    mutate(pvalue.adj = 10**(floor(abs(log(pvalue.adj,10)) * sign(log(pvalue.adj,10))))) %>%
    mutate(GO_class=ifelse(Description=="positive regulation of cell cycle","Cell\nproliferation\nand growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="cervical cancer proliferation cluser","Cell\nproliferation\nand growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="mitotic cell cycle","Cell\nproliferation\nand growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="DNA replication","Cell\nproliferation\nand growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="ribosome biogenesis","Cell\nproliferation\nand growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="chemotaxis","EMT, extracellular\nmatrix and growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="cell adhesion","EMT, extracellular\nmatrix and growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="extracellular space","EMT, extracellular\nmatrix and growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="extracellular matrix","EMT, extracellular\nmatrix and growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="regulation of locomotion","EMT, extracellular\nmatrix and growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="cell surface","EMT, extracellular\nmatrix and growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="regulation of cell motility","EMT, extracellular\nmatrix and growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="epithelial to mesenchymal transition","EMT, extracellular\nmatrix and growth",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="circulatory system development","Hypoxia\nresponse and\nangiogenesis",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="Down-regulated in breast cancer after HIF1A KO","Hypoxia\nresponse and\nangiogenesis",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="Down-regulated in breast cancer after HIF1A and HIF2A KO","Hypoxia\nresponse and\nangiogenesis",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="Upregulated by hypoxia in prostate cancer","Hypoxia\nresponse and\nangiogenesis",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="angiogenesis","Hypoxia\nresponse and\nangiogenesis",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="response to hypoxia","Hypoxia\nresponse and\nangiogenesis",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="vascular process in circulatory system","Hypoxia\nresponse and\nangiogenesis",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="immune response","Cell stress\nand cytokine",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="inflammatory response","Cell stress\nand cytokine",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="regulation of cytokine production","Cell stress\nand cytokine",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="secretion by cell","Cell stress\nand cytokine",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="cilium assembly","Cilium",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="cilium movement","Cilium",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="axoneme","Cilium",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="ciliary plasm","Cilium",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="ciliary basal body","Cilium",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="microtubule bundle formation","Cilium",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="mitochondrion organization","Metab.",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="Carbon metabolism","Metab.",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="positive regulation of protein phosphorylation","Metab.",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="autophagy", "Autophagy\nand senescence", GO_class)) %>%
    mutate(GO_class=ifelse(Description=="Cellular senescence", "Autophagy\nand senescence",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="Platinum drug resistance","Response\nto drugs",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="glutathione metabolic process","Response\nto drugs",GO_class)) %>%
    mutate(GO_class=ifelse(grepl("Erlotinib",Description),"Response\nto drugs",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="EGFR tyrosine kinase inhibitor resistance","Response\nto drugs",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="Up-regulated in breast cancer after estradiol treatment","Response\nto drugs",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="Down-regulated in lung cancer after gefitinib treatment","Response\nto drugs",GO_class)) %>%
    mutate(GO_class=ifelse(Description=="Metabolism of xenobiotics by cytochrome P450","Response\nto drugs",GO_class)) %>% 
    mutate(time=ifelse(time=="B","R",time)) %>%
    mutate(time=factor(time,levels=c("5","8","11","R"))) %>%
    mutate(Description=str_to_title(Description)) %>%
    dplyr::filter(GO_class!="Cancer related")  %>%
    mutate(GO_class=factor(GO_class, levels = GO_class %>% unique %>% .[c(1,2,3,6,4,7,5,8)]))


##### small pdf file ... 

g.svg <- ggplot(data = gseall.list.df.to_plot, aes(y = Description , x = time , color = NES, size = -log(pvalue.adj,10))) +  
    geom_point() +
    scale_color_gradientn(colours=c(scales::muted("blue"),"white" ,"white", scales::muted("red")),values=scales::rescale(c(-2,-0.3,0.3,2)),breaks=c(-2,-0.3,0.3,2),limits=c(-2,2)) +
    theme_bw() +
    ylab("") +
    xlab("") +
    ## scale_x_discrete(breaks=c("5","8","11","R"),limits=c("5","8","11","R")) +
    scale_size_continuous(range=c(1.5,3.5),breaks= c(2,5,8),labels=c("2","5","> 8")) + 
    theme(text = element_text(size=9)) +
    theme(axis.ticks = element_line(size=0.1),axis.ticks.length = unit(0.04,"cm")) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.minor = element_line(size = 0.2), panel.grid.major = element_line(size = 0.2)) +
    labs(size="-log(q.val)")  + 
    theme(legend.margin = margin(0,3,0,0),legend.box.margin = margin(0,3,0,0)) +
    theme( legend.text=element_text(size=8),legend.title=element_text(size=8)) +
    theme(plot.margin = unit(c(3,0,-2,0), "pt")) +
    facet_grid(vars(GO_class),vars(geno),scales="free_y",space="free_y") +
    theme(strip.placement="outside") +
    theme( strip.text.y = element_text(size=7),strip.background=element_blank()) +
    theme(legend.key.size = unit(0.35, "cm"), legend.text=element_text(size=7),legend.title=element_text(size=7))



out_dotplot_small_pdf<- file.path(out_dir,"dotplot_small.pdf")
ggpubr::ggexport(g.svg, filename = out_dotplot_small_pdf,width = 7, height = 5 +  nrow(gseall.list.df.to_plot)*(1/100))


out_dotplot_small_svg <- file.path(out_dir,"dotplot_small.svg")
ggsave(g.svg, filename = out_dotplot_small_svg,width = 7, height = 2 +  nrow(gseall.list.df.to_plot)*(1/100))

saveRDS(g.svg,file=file.path(Robj_dir,"gsea.rdata"))

