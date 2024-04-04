rm(list=ls())
library(ggplot2)
library(ggh4x)
library(ggpubr)
library(Cairo)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(tidyverse)
library(rstatix)
library(stats)
library(readxl)
library(reshape2)
library(olsrr)
library(MASS)
library(ordinal)
library(emmeans)
library(stats)
library(contrast)
library(multcomp)
library(nlme)
library(lme4)

out_dir <- file.path("./results")
dir.create(out_dir,showWarnings=FALSE, recursive=TRUE)

##### read df data


df_new_combined.rds = readRDS("../data/20240308_merged_data_of_viability_and_organoids_size.rds")
tab.combined <- df_new_combined.rds %>%
    dplyr::select(pictures_included, drug, acidosis.condition, geno, batch, pair, log_mean_total_sphere_surface_area, logv) %>%
    mutate(drug = drug %>% str_replace("E\\+G","E_G")) %>% 
    mutate(drug.condi.quant =  drug  %>% str_replace("10",".10") %>% str_replace("25",".25") %>%
               factor(levels=c("Untreated", "E.10", "G.10","E_G.10", "E.25", "G.25", "E_G.25"))) %>%
        dplyr::filter(! grepl("25", drug.condi.quant)) %>%
    mutate(drug.condi.quant = drug.condi.quant %>% factor(levels=c("Untreated", "E.10", "G.10","E_G.10"))) %>% 
    dplyr::select(-drug) %>% 
    dplyr::rename(pairv = pair, pairs = batch,   logs = log_mean_total_sphere_surface_area) %>%
    mutate(acidosis.condition  = acidosis.condition %>%
               str_replace("Ctrl","ctrl") %>% str_replace("AA->pH7.4","bb") %>% str_replace("AA","aa") %>%
               factor(c("ctrl","aa","bb"))) %>%
    mutate(geno = ifelse(geno =="p53R273H", "R273H",geno) %>% factor) %>%
    mutate(pairs = pairs %>% factor %>% as.numeric %>% as.factor, pairv = pairv %>% factor %>% as.numeric %>% as.factor) 





df.viability <- tab.combined %>%
    dplyr::select(-logs, -pairs, -pictures_included) %>%
    dplyr::rename(pair = pairv) %>%
    drop_na() %>%
    group_by(geno, pair, drug.condi.quant) %>%
    add_count %>%
    ungroup %>%
    dplyr::filter(n==3)


contrasts.drug <- list()
contrasts.df.list <- list()
ano.df.list <- list()
coeffs <- list()
df.temp.list <- list()

##
for (geno.sel in df.viability$geno %>% levels)
{
    lm.logv <- lm(logv ~       acidosis.condition*drug.condi.quant + pair,
                  data=df.viability %>% dplyr::filter(geno==geno.sel))
    lm.logv.test <- lm(logv ~       acidosis.condition*drug.condi.quant + pair,
                      data=df.viability %>% dplyr::filter(geno==geno.sel))
    a <- predict(lm.logv.test, df.viability %>% dplyr::filter(geno==geno.sel))
    residuals  <-  df.viability %>% dplyr::filter(geno==geno.sel) %>% .$logv -a 
    lm.logv.test$coefficients[grepl("pair", lm.logv.test$coefficients %>% names)]  <- 0
    lm.logv.test$coefficients[1] <- 11
    a <- predict(lm.logv.test, df.viability %>% dplyr::filter(geno==geno.sel)) + residuals
    df.temp.list[[geno.sel]]  <- df.viability %>% dplyr::filter(geno==geno.sel) %>% mutate(logv.bc=a)
    coeffs[[geno.sel]] <- lm.logv %>% summary() %>% .$coeff %>% 
        as.data.frame() %>%
        setNames(c("Estimate","ste","tval","p")) %>% 
        rownames_to_column("coeff") %>%
        dplyr::select(coeff, Estimate, p) %>%
        dplyr::filter(! grepl("pair|Intercept", coeff)) %>%
        mutate(coeff = coeff %>% str_replace("acidosis.condition","") %>% str_replace("drug.condi.quant","")) %>%
        mutate(geno = geno.sel) %>%
        mutate(coeff.dcast = paste0(ifelse(grepl("aa",coeff),  "AAw11", ifelse(grepl("bb",coeff),"AA->pH7.4", "")),
                                    ifelse(grepl(":", coeff)," x ",""),
                                    ifelse(grepl("10", coeff),"D.10",""))) %>%
        mutate(drug = coeff %>% str_replace("\\..*", "") %>% str_replace("(aa|bb):*","") %>% str_replace("_","+")) %>%
        mutate(condi = ifelse(drug =="", "Untreated", "Treated"))               
    ##
    cont.bb.aa <- list()
    cont.c.bb <- list()
    cont.aa.c <- list()
    contrasts.df.list[[geno.sel]] <- list()
    pair.sel  <-  df.viability %>% dplyr::filter(geno == geno.sel) %>% .$pair %>% .[1] %>% as.character
    for (elt in df.viability$drug.condi.quant %>% levels %>% .[-1])
    {
        h.init <- lm.logv$coefficients %>% names
        h.c <- grepl(paste0("^drug\\..*","quant", elt), h.init) %>% as.numeric %>% matrix(nrow=1)
        h.aa <- grepl(paste0("^((?!.*bb.*).*)","quant",elt), h.init, perl=TRUE) %>% as.numeric  %>% matrix(nrow=1)
        h.bb <- grepl(paste0("^((?!.*aa.*).*)","quant",elt), h.init, perl=TRUE) %>% as.numeric  %>% matrix(nrow=1)
        cont.bb.aa[[elt]] <- glht(lm.logv, linfct = h.bb - h.aa) %>% tidy
        cont.c.bb[[elt]] <- glht(lm.logv, linfct = h.c - h.bb) %>% tidy
        cont.aa.c[[elt]] <- glht(lm.logv, linfct = h.aa - h.c) %>% tidy
        contrasts.df.list[[geno.sel]][[elt]] <- data.frame(
            c("bb" ,"aa" , cont.bb.aa[[elt]]$adj.p.value , cont.bb.aa[[elt]]$estimate, elt),
            c("ctrl", "bb"  , cont.c.bb[[elt]]$adj.p.value , cont.c.bb[[elt]]$estimate, elt),
            c("aa" , "ctrl" , cont.aa.c[[elt]]$adj.p.value , cont.aa.c[[elt]]$estimate, elt))  %>%
            t %>%
            as.data.frame %>% 
            remove_rownames %>%
            setNames(c("group1", "group2", "p", "Contrast", "drug")) %>%
            mutate(Contrast = Contrast %>% as.numeric %>% round(3)) %>%
            mutate(p = p %>% as.numeric %>% round(3))
    }
    contrasts.df.list[[geno.sel]] <- contrasts.df.list[[geno.sel]] %>%
        do.call(rbind.data.frame, .) %>%
        remove_rownames  %>%
        dplyr::select(-Contrast)
    ##
#### get values to plot
    contrasts.drug[[geno.sel]] <- list()
    for (elt in df.viability$acidosis.condition %>% levels)
    {
        Unt <-   list(acidosis.condition = elt , pair = pair.sel, drug.condi.quant = "Untreated")
        E10 <-   list(acidosis.condition = elt , pair = pair.sel, drug.condi.quant = "E.10")
        G10 <-   list(acidosis.condition = elt , pair = pair.sel, drug.condi.quant =  "G.10")
        EG10 <-  list(acidosis.condition = elt , pair = pair.sel, drug.condi.quant = "E_G.10")
        ##
        cont.E10 <- contrast(lm.logv, E10, Unt)
        cont.G10 <- contrast(lm.logv, G10, Unt)
        cont.EG10 <- contrast(lm.logv, EG10, Unt)
        ##
        contrasts.drug[[geno.sel]][[elt]] <- data.frame(
            c("E.10" , cont.E10$Contrast, cont.E10$Lower, cont.E10$Upper, elt),
            c("G.10" , cont.G10$Contrast, cont.G10$Lower, cont.G10$Upper, elt),
            c("E_G.10" , cont.EG10$Contrast, cont.EG10$Lower, cont.EG10$Upper,  elt)) %>% 
            t %>%
            as.data.frame %>% 
            remove_rownames %>%
            setNames(c("drug", "Contrast", "Lower.CI", "Upper.CI","pH")) %>%
            mutate(Contrast = Contrast %>% as.numeric %>% round(3)) %>%
            mutate(Lower.CI = Lower.CI %>% as.numeric %>% round(3)) %>%
            mutate(Upper.CI = Upper.CI %>% as.numeric %>% round(3))
    }
    contrasts.drug[[geno.sel]] <- contrasts.drug[[geno.sel]] %>%
        do.call(rbind.data.frame, .) %>% remove_rownames
    ##
    ##
    ano.df.list[[geno.sel]] <- list()
    for (elt in df.viability$acidosis.condition %>% levels)
    {
        dosage ="10"
        df.viability.null <- df.viability %>%
            mutate(drug.condi.quant = drug.condi.quant %>%
                       as.character %>% 
                       ifelse(acidosis.condition == elt & (grepl(dosage, drug.condi.quant) | drug.condi.quant == "Untreated" ), "null", .) %>%
                       factor(c("Untreated", "E.10", "G.10","E_G.10", "null")))
        lm.logv.null <- lm(logv ~  acidosis.condition*drug.condi.quant + pair,
                           data=df.viability.null %>% dplyr::filter(geno==geno.sel))
        elt.dosage  <-  paste0(elt,".", dosage)
        ano.df.list[[geno.sel]][[elt.dosage]] <- anova(lm.logv, lm.logv.null)
    }
    ##
    ano.df.list[[geno.sel]]  <- ano.df.list[[geno.sel]] %>%
        lapply(as.data.frame) %>%
        lapply(drop_na) %>%
        lapply("[", "Pr(>F)") %>%
        lapply(setNames, "p") %>%
        mapply(mutate, . , ph.dosage =names(.), SIMPLIFY=FALSE) %>%
        do.call(rbind.data.frame, .) %>%
        separate(ph.dosage, into=c("ph","dosage"), sep="\\.") %>% 
        remove_rownames
}

df_logv.bc <- df.temp.list %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames

ano.df  <-  ano.df.list %>%
    mapply(mutate, ., geno = names(.), SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames  %>%
    mutate(p = p %>% formatC(format = "e", digit = 2))

contrasts.drug.df  <-  contrasts.drug %>%
    mapply(mutate, ., geno = names(.), SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames

contrasts.df  <-  contrasts.df.list %>%
    mapply(mutate, ., geno = names(.), SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames

stats.df.barplots <-  contrasts.drug.df %>%
    mutate(ypos.temp=ifelse(Contrast > 0, Upper.CI, 0)) %>%
    dplyr::select(-Contrast, -Lower.CI, - Upper.CI) %>%
    mutate(group1=pH) %>% 
    left_join(contrasts.df, by=c("geno", "drug", "group1"), multiple = "all") %>%
    arrange(geno, drug, desc(ypos.temp))  %>%
    group_by(geno, drug) %>%
    mutate(ypos = max(ypos.temp))  %>%
    ungroup %>%               
    dplyr::filter(p<=0.1) %>%
    mutate(p.numeric = p , p = ifelse(p  <= 0.05, paste0(p, "*"), p)) %>% 
    group_by(geno, drug) %>%
    mutate(id = row_number()) %>%
    mutate(ypos= ypos + 0.1 + (id-1) * 0.6) %>%
    ungroup %>%
    dplyr::select(drug, pH, geno, y.position=ypos, group1, group2, p, p.numeric)   %>%
    dplyr::filter(group2 != "aa") %>%
    mutate(drug = factor(drug, levels =c("E.10", "G.10", "E_G.10" )))


floor_pval <- function(x)
{
    a <- abs(log(x,10))
    sign.x <- sign(log(x,10))
    y <- ifelse(a %% 1 < -(log(0.05, 10) +1),  10**ceiling((a*sign.x)), 10**floor((a*sign.x)))
    return(y)
}

coeffs.df <- coeffs %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames %>%
    mutate(log.p = floor_pval(p) %>% log(.,10) %>% ifelse(. < -3, -3, .)) %>%
    mutate(condi = factor(condi, c("Untreated","Treated"))) %>%
    mutate(drug = ifelse(drug =="", "Untreated", drug) %>% factor(., c("Untreated", "G", "E", "E+G"))) %>%
    mutate(Estimate = round(Estimate, 1))

#### fig S3 dotplot
##################
g <- ggplot(coeffs.df) +
    geom_point( aes(x=drug, y = coeff.dcast, size = -log.p, colour = Estimate)) +
    scale_size_continuous(range=c(0.5,8),breaks= c(0, 1,2,3),labels=c("0", "1","2","> 3")) +
    facet_grid2(vars(condi), vars(geno), scales= "free", independent="x", space = "free_y") + 
    scale_color_gradientn(colours=c(scales::muted("blue"),"white" ,"white", scales::muted("red")),
                          values=scales::rescale(c(-2,-0.2,0.2,2)),breaks=c(-2,-0.2,0.2,2),limits=c(-2,2)) +
    theme_bw() +
    ggnewscale::new_scale_colour() + 
    geom_text(data = coeffs.df %>% dplyr::filter(p<0.05),
              aes(colour = ifelse(abs(Estimate)>1.3, "white", "black"), label = Estimate, x = drug, y = coeff.dcast),
              size = 3 ) +
    scale_colour_manual(values = c("black", "white")) +
    guides(colour="none") +
    ggtitle("Per genotype, log10(viability) ~ pH * drug + batch ")
ggpubr::ggexport(g, filename = file.path(out_dir, "dotplot_logv_fig3.pdf"), width = 8, height = 6)


###################
#### fig 3 boxplots logv
stat.test.10.logv <- ano.df %>%
    dplyr::rename(acidosis.condition = ph) %>%
    mutate(group1 = "Untreated", group2 = "E_G.10") %>%
    mutate(y.position = max(df.viability[, "logv"])* 1.05) %>%
    mutate(p.numeric = as.numeric(p)) %>%
    mutate(p = ifelse(p.numeric < 0.05, paste0(p,"*"), p))

## boxplot
g <- ggboxplot(df.viability %>% ## boxplot
               mutate(drug.condi.quant = drug.condi.quant %>% factor(c("Untreated", "E.10","G.10", "E_G.10"))),
               x = "acidosis.condition", y = "logv", color = "drug.condi.quant", add.params = list(size = 1, shape = 1),
               add = "jitter",
               facet.by="geno", nrow=1, ncol = 4) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, 
                                     vjust = 1, 
                                     hjust=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_color_manual(values=c("#F3AF7F","#00A087", "#61ACCE","#905EA9")) +
    stat_pvalue_manual(stat.test.10.logv %>% dplyr::filter(p.numeric < 0.05), 
                       x = "acidosis.condition", size = 3,label = "p",
                       fontface = "bold") +
    stat_pvalue_manual(stat.test.10.logv %>% dplyr::filter(p.numeric >= 0.05), 
                       x = "acidosis.condition", size = 3,label = "p")
ggpubr::ggexport(g, filename = file.path(out_dir, "boxplots_logv_fig3.pdf"), width = 10, height = 4)

g <- ggboxplot(df_logv.bc %>%   # boxplot
               mutate(drug.condi.quant = drug.condi.quant %>% factor(c("Untreated", "E.10","G.10", "E_G.10"))),
               x = "acidosis.condition", y = "logv.bc", color = "drug.condi.quant", add.params = list(size = 1, shape = 1),
               add = "jitter",
               facet.by="geno", nrow=1, ncol = 4) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, 
                                     vjust = 1, 
                                     hjust=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_color_manual(values=c("#F3AF7F","#00A087", "#61ACCE","#905EA9")) +
    stat_pvalue_manual(stat.test.10.logv %>% dplyr::filter(p.numeric < 0.05), 
                       x = "acidosis.condition", size = 3,label = "p",
                       fontface = "bold") +
    stat_pvalue_manual(stat.test.10.logv %>% dplyr::filter(p.numeric >= 0.05), 
                       x = "acidosis.condition", size = 3,label = "p")
ggpubr::ggexport(g, filename = file.path(out_dir, "boxplots_logvBC_fig3.pdf"), width = 10, height = 4)


## barplots fig3 logv
bar.df.fig3.logv <- contrasts.drug.df %>%
    mutate(drug = drug %>% factor(c("E.10","G.10","E_G.10"))) %>%
    mutate(pH = pH %>% factor(c("ctrl", "aa", "bb")))  

g  <-   ggplot(bar.df.fig3.logv,  aes(x = pH, y = Contrast)) +
    geom_errorbar(aes(ymin = ifelse(Contrast < 0, Lower.CI, Contrast),
                      ymax = ifelse(Contrast > 0, Upper.CI, Contrast)), width = 0.25) +
    geom_bar(aes(fill = pH), stat="identity", color="black", position=position_dodge(), width = 0.6) +
    scale_fill_manual(values=c("#8A8F9A", "#FFB700", "#E8217B")) +
    geom_hline(yintercept = 0, color = "black")+
    theme_bw() +
    facet_grid(fct_relevel(drug,c("E.10","G.10","E_G.10")) ~ geno, scales = "free_x") +
    stat_pvalue_manual(stats.df.barplots %>% dplyr::filter(!grepl("25", drug), p.numeric >= 0.05), label = "p") +
    stat_pvalue_manual(stats.df.barplots %>% dplyr::filter(!grepl("25", drug), p.numeric < 0.05), label = "p",
                       fontface="bold") +
    ylim(c(NA, 1))
ggpubr::ggexport(g, filename = file.path(out_dir, "barplots_logv_fig3.pdf"), width = 6, height = 5)



##################
#### same for area
##################

tab.surface <- tab.combined %>%
    dplyr::select(-logv, -pairv) %>%
    dplyr::rename(pair = pairs) %>%
    drop_na() %>%
    group_by(geno, pair, drug.condi.quant) %>%
    add_count() %>%
    ungroup %>% 
    dplyr::filter(n==3)


    
contrasts.drug <- list()
contrasts.df.list <- list()
ano.df.list <- list()
coeffs <- list()
df.temp.list <- list()
##
for (geno.sel in tab.surface$geno %>% levels)
{
    lm.logs <- lm(logs ~       acidosis.condition*drug.condi.quant + pair,
                  data=tab.surface %>% dplyr::filter(geno==geno.sel))
    lm.logs.test <- lm(logs ~       acidosis.condition*drug.condi.quant + pair,
                      data=tab.surface %>% dplyr::filter(geno==geno.sel))
    a <- predict(lm.logs.test, tab.surface %>% dplyr::filter(geno==geno.sel))
    residuals  <-  tab.surface %>% dplyr::filter(geno==geno.sel) %>% .$logs -a 
    lm.logs.test$coefficients[grepl("pair", lm.logs.test$coefficients %>% names)]  <- 0
    lm.logs.test$coefficients[0] <- 7
    a <- predict(lm.logs.test, tab.surface %>% dplyr::filter(geno==geno.sel)) + residuals
    df.temp.list[[geno.sel]]  <- tab.surface %>% dplyr::filter(geno==geno.sel) %>% mutate(logs.bc=a)
    coeffs[[geno.sel]] <- lm.logs %>% summary() %>% .$coeff %>% 
        as.data.frame() %>%
        setNames(c("Estimate","ste","tval","p")) %>% 
        rownames_to_column("coeff") %>%
        dplyr::select(coeff, Estimate, p) %>%
        dplyr::filter(! grepl("pair|Intercept", coeff)) %>%
        mutate(coeff = coeff %>% str_replace("acidosis.condition","") %>% str_replace("drug.condi.quant","")) %>%
        mutate(geno = geno.sel) %>%
        mutate(coeff.dcast = paste0(ifelse(grepl("aa",coeff),  "AAw11", ifelse(grepl("bb",coeff),"AA->pH7.4", "")),
                                    ifelse(grepl(":", coeff)," x ",""),
                                    ifelse(grepl("10", coeff),"D.10",""))) %>%
        mutate(drug = coeff %>% str_replace("\\..*", "") %>% str_replace("(aa|bb):*","") %>% str_replace("_","+")) %>%
        mutate(condi = ifelse(drug =="", "Untreated", "Treated"))                   
    ##
    cont.bb.aa <- list()
    cont.c.bb <- list()
    cont.aa.c <- list()
    contrasts.df.list[[geno.sel]] <- list()
    pair.sel  <-  tab.surface %>% dplyr::filter(geno == geno.sel) %>% .$pair %>% .[1] %>% as.character
    for (elt in tab.surface$drug.condi.quant %>% levels %>% .[-1])
    {
        h.init <- lm.logs$coefficients %>% names
        h.c <- grepl(paste0("^drug\\..*","quant",elt), h.init) %>% as.numeric %>% matrix(nrow=1)
        h.aa <- grepl(paste0("^((?!.*bb.*).*)","quant",elt), h.init, perl=TRUE) %>% as.numeric  %>% matrix(nrow=1)
        h.bb <- grepl(paste0("^((?!.*aa.*).*)","quant",elt), h.init, perl=TRUE) %>% as.numeric  %>% matrix(nrow=1)
        cont.bb.aa[[elt]] <- glht(lm.logs, linfct = h.bb - h.aa) %>% tidy
        cont.c.bb[[elt]] <- glht(lm.logs, linfct = h.c - h.bb) %>% tidy
        cont.aa.c[[elt]] <- glht(lm.logs, linfct = h.aa - h.c) %>% tidy
        contrasts.df.list[[geno.sel]][[elt]] <- data.frame(
            c("bb" ,"aa" , cont.bb.aa[[elt]]$adj.p.value , cont.bb.aa[[elt]]$estimate, elt),
            c("ctrl", "bb"  , cont.c.bb[[elt]]$adj.p.value , cont.c.bb[[elt]]$estimate, elt),
            c("aa" , "ctrl" , cont.aa.c[[elt]]$adj.p.value , cont.aa.c[[elt]]$estimate, elt))  %>%
            t %>%
            as.data.frame %>% 
            remove_rownames %>%
            setNames(c("group1", "group2", "p", "Contrast", "drug")) %>%
            mutate(Contrast = Contrast %>% as.numeric %>% round(3)) %>% 
            mutate(p = p %>% as.numeric %>% signif(2))
    }
    contrasts.df.list[[geno.sel]] <- contrasts.df.list[[geno.sel]] %>%
        do.call(rbind.data.frame, .) %>%
        remove_rownames %>%
        dplyr::select(-Contrast)
    ##
#### get values to plot
    contrasts.drug[[geno.sel]] <- list()
    for (elt in tab.surface$acidosis.condition %>% levels)
    {
        Unt <-   list(acidosis.condition = elt ,  pair = pair.sel, drug.condi.quant = "Untreated")
        E10 <-   list(acidosis.condition = elt ,  pair = pair.sel, drug.condi.quant = "E.10")
        G10 <-   list(acidosis.condition = elt ,  pair = pair.sel, drug.condi.quant =  "G.10")
        EG10 <-  list(acidosis.condition = elt ,  pair = pair.sel, drug.condi.quant = "E_G.10")
        ##
        cont.E10 <- contrast(lm.logs, E10, Unt)
        cont.G10 <- contrast(lm.logs, G10, Unt)
        cont.EG10 <- contrast(lm.logs, EG10, Unt)
        ##
        contrasts.drug[[geno.sel]][[elt]] <- data.frame(
            c("E.10" , cont.E10$Contrast, cont.E10$Lower, cont.E10$Upper, elt),
            c("G.10" , cont.G10$Contrast, cont.G10$Lower, cont.G10$Upper, elt),
            c("E_G.10" , cont.EG10$Contrast, cont.EG10$Lower, cont.EG10$Upper,  elt)) %>% 
            t %>%
            as.data.frame %>% 
            remove_rownames %>%
            setNames(c("drug", "Contrast", "Lower.CI", "Upper.CI","pH")) %>%
            mutate(Contrast = Contrast %>% as.numeric %>% round(3)) %>%
            mutate(Lower.CI = Lower.CI %>% as.numeric %>% round(3)) %>%
            mutate(Upper.CI = Upper.CI %>% as.numeric %>% round(3))
    }
    contrasts.drug[[geno.sel]] <- contrasts.drug[[geno.sel]] %>%
        do.call(rbind.data.frame, .) %>% remove_rownames
    ##
    ##
    ano.df.list[[geno.sel]] <- list()
    for (elt in tab.surface$acidosis.condition %>% levels)
    {
        for (dosage in c("10"))
        {
            tab.surface.null <- tab.surface %>%
                mutate(drug.condi.quant = drug.condi.quant %>%
                           as.character %>% 
                           ifelse(acidosis.condition == elt & (grepl(dosage, drug.condi.quant) | drug.condi.quant == "Untreated" ), "null", .) %>%
                           factor(c("Untreated", "E.10", "G.10","E_G.10", "null")))
            lm.logs.null <- lm(logs ~  acidosis.condition*drug.condi.quant + pair,
                               data=tab.surface.null %>% dplyr::filter(geno==geno.sel))
            elt.dosage  <-  paste0(elt,".", dosage)
            ano.df.list[[geno.sel]][[elt.dosage]] <- anova(lm.logs, lm.logs.null)
        }
    }
    ##
    ano.df.list[[geno.sel]]  <- ano.df.list[[geno.sel]] %>%
        lapply(as.data.frame) %>%
        lapply(drop_na) %>%
        lapply("[", "Pr(>F)") %>%
        lapply(setNames, "p") %>%
        mapply(mutate, . , ph.dosage =names(.), SIMPLIFY=FALSE) %>%
        do.call(rbind.data.frame, .) %>%
        separate(ph.dosage, into=c("ph","dosage"), sep="\\.") %>% 
        remove_rownames
}

df_area.bc <- df.temp.list %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames

ano.df  <-  ano.df.list %>%
    mapply(mutate, ., geno = names(.), SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames  %>%
    mutate(p = p %>% formatC(format = "e", digit = 2))

contrasts.drug.df  <-  contrasts.drug %>%
    mapply(mutate, ., geno = names(.), SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames

contrasts.df  <-  contrasts.df.list %>%
    mapply(mutate, ., geno = names(.), SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames 

stats.df.barplots <-  contrasts.drug.df %>%
    mutate(ypos.temp=ifelse(Contrast > 0, Upper.CI, 0)) %>%
    dplyr::select(-Contrast, -Lower.CI, - Upper.CI) %>%
    mutate(group1=pH) %>% 
    left_join(contrasts.df, by=c("geno", "drug", "group1"), multiple = "all") %>%
    arrange(geno, drug, desc(ypos.temp))  %>%
    group_by(geno, drug) %>%
    mutate(ypos = max(ypos.temp))  %>%
    ungroup %>%               
    dplyr::filter(p < 0.1) %>% 
    mutate(p.numeric = p , p = ifelse(p  <= 0.05, paste0(p, "*"), p)) %>% 
    group_by(geno, drug) %>%
    mutate(id = row_number()) %>%
    mutate(ypos= ypos + 0.1 + (id-1) * 0.3) %>%
    ungroup %>%
    dplyr::select(drug, pH, geno, y.position=ypos, group1, group2, p, p.numeric)  %>%
    dplyr::filter(group2 != "aa")  %>% 
    mutate(drug = factor(drug, levels =c("E.10", "G.10", "E_G.10" )))



coeffs.df <- coeffs %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames %>%
    mutate(log.p = floor_pval(p) %>% log(.,10) %>% ifelse(. < -3, -3, .)) %>%
    mutate(condi = factor(condi, c("Untreated","Treated"))) %>%
    mutate(drug = ifelse(drug =="", "Untreated", drug) %>% factor(., c("Untreated", "G", "E", "E+G"))) %>%
    mutate(Estimate = round(Estimate, 1))

#### fig S3 dotplot
##################
g <- ggplot(coeffs.df) +
    geom_point( aes(x=drug, y = coeff.dcast, size = -log.p, colour = Estimate)) +
    scale_size_continuous(range=c(0.5,8),breaks= c(0, 1,2,3),labels=c("0", "1","2","> 3")) +
    facet_grid2(vars(condi), vars(geno), scales= "free", independent="x", space = "free_y") + 
    scale_color_gradientn(colours=c(scales::muted("blue"),"white" ,"white", scales::muted("red")),
                          values=scales::rescale(c(-1,-0.1,0.1,1)),breaks=c(-1,-0.1,0.1,1),limits=c(-1,1)) +
    theme_bw() +
    ggnewscale::new_scale_colour() + 
    geom_text(data = coeffs.df %>% dplyr::filter(p<0.05),
              aes(colour = ifelse(abs(Estimate)>0.8, "white", "black"), label = Estimate, x = drug, y = coeff.dcast),
              size = 3 ) +
    scale_colour_manual(values = c("black", "white")) +
    guides(colour="none") +
    ggtitle("Per genotype, log10(surface) ~ pH * drug + batch ")
ggpubr::ggexport(g, filename = file.path(out_dir, "dotplot_logs_fig3.pdf"), width = 8, height = 6)

###################
#### fig 3 boxplots logs
stat.test.10.logs <- ano.df %>%
    dplyr::rename(acidosis.condition = ph) %>%
    mutate(group1 = "Untreated", group2 = "E_G.10") %>%
    mutate(y.position = max(tab.surface[, "logs"])* 1.05) %>%
    mutate(p.numeric = as.numeric(p)) %>%
    mutate(p = ifelse(p.numeric < 0.05, paste0(p,"*"), p))

## boxplot
g <- ggboxplot(tab.surface %>%   # boxplot
               mutate(drug.condi.quant = drug.condi.quant %>% factor(c("Untreated", "E.10","G.10", "E_G.10"))),
               x = "acidosis.condition", y = "logs", color = "drug.condi.quant", add.params = list(size = 1, shape = 1),
               add = "jitter",
               facet.by="geno", nrow=1, ncol = 4) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, 
                                     vjust = 1, 
                                     hjust=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_color_manual(values=c("#F3AF7F","#00A087", "#61ACCE","#905EA9")) +
    stat_pvalue_manual(stat.test.10.logs %>% dplyr::filter(p.numeric < 0.05), 
                       x = "acidosis.condition", size = 3,label = "p",
                       fontface = "bold") +
    stat_pvalue_manual(stat.test.10.logs %>% dplyr::filter(p.numeric >= 0.05), 
                       x = "acidosis.condition", size = 3,label = "p")
ggpubr::ggexport(g, filename = file.path(out_dir, "boxplots_logs_fig3.pdf"), width = 10, height = 4)

g <- ggboxplot(df_area.bc %>%   # boxplot
               mutate(drug.condi.quant = drug.condi.quant %>% factor(c("Untreated", "E.10","G.10", "E_G.10"))),
               x = "acidosis.condition", y = "logs.bc", color = "drug.condi.quant", add.params = list(size = 1, shape = 1),
               add = "jitter",
               facet.by="geno", nrow=1, ncol = 4) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, 
                                     vjust = 1, 
                                     hjust=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_color_manual(values=c("#F3AF7F","#00A087", "#61ACCE","#905EA9")) +
    stat_pvalue_manual(stat.test.10.logs %>% dplyr::filter(p.numeric < 0.05), 
                       x = "acidosis.condition", size = 3,label = "p",
                       fontface = "bold") +
    stat_pvalue_manual(stat.test.10.logs %>% dplyr::filter(p.numeric >= 0.05), 
                       x = "acidosis.condition", size = 3,label = "p")
ggpubr::ggexport(g, filename = file.path(out_dir, "boxplots_logsBC_fig3.pdf"), width = 10, height = 4)


## barplots fig3 logs
bar.df.fig3.logs <- contrasts.drug.df %>%
    mutate(drug = drug %>% factor(c("E.10","G.10","E_G.10"))) %>%
    mutate(pH = pH %>% factor(c("ctrl", "aa", "bb")))  

g  <-   ggplot(bar.df.fig3.logs,  aes(x = pH, y = Contrast)) +
    geom_errorbar(aes(ymin = ifelse(Contrast < 0, Lower.CI, Contrast),
                      ymax = ifelse(Contrast > 0, Upper.CI, Contrast)), width = 0.25) +
    geom_bar(aes(fill = pH), stat="identity", color="black", position=position_dodge(), width = 0.6) +
    scale_fill_manual(values=c("#8A8F9A", "#FFB700", "#E8217B")) +
    geom_hline(yintercept = 0, color = "black")+
    theme_bw() +
    facet_grid(fct_relevel(drug,c("E.10","G.10","E_G.10")) ~ geno, scales = "free_x") +
    stat_pvalue_manual(stats.df.barplots %>% dplyr::filter(!grepl("25", drug), p.numeric >= 0.05), label = "p") +
    stat_pvalue_manual(stats.df.barplots %>% dplyr::filter(!grepl("25", drug), p.numeric < 0.05), label = "p",
                       fontface="bold") +
    ylim(c(NA, 1))
ggpubr::ggexport(g, filename = file.path(out_dir, "barplots_logs_fig3.pdf"), width = 6, height = 5)



