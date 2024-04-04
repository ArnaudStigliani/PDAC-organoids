rm(list=ls())
library(ggplot2)
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
library(sva)
#### test new table
out_dir <- file.path("./results")
dir.create(out_dir,showWarnings=FALSE, recursive=TRUE)
df_new_combined.rds = readRDS("../data/20240308_merged_data_of_viability_and_organoids_size.rds")
tab.combined <- df_new_combined.rds %>%
    dplyr::select(pictures_included, drug, acidosis.condition, geno, batch, pair, log_mean_total_sphere_surface_area, logv) %>%
    mutate(drug = drug %>% str_replace("E\\+G","E_G")) %>% 
    mutate(drug.condi.quant =  drug  %>% str_replace("10",".10") %>% str_replace("25",".25") %>%
               factor(levels=c("Untreated", "E.10", "G.10","E_G.10", "E.25", "G.25", "E_G.25"))) %>%
    dplyr::select(-drug) %>% 
    dplyr::rename(pairv = pair, pairs = batch,   logs = log_mean_total_sphere_surface_area) %>%
    mutate(acidosis.condition  = acidosis.condition %>%
               str_replace("Ctrl","ctrl") %>% str_replace("AA->pH7.4","bb") %>% str_replace("AA","aa") %>%
               factor(c("ctrl","aa","bb"))) %>%
    mutate(geno = ifelse(geno =="p53R273H", "R273H",geno) %>% factor) %>%
    mutate(pairs = pairs %>% factor %>% as.numeric %>% as.factor,
           pairv = pairv %>% factor %>% as.numeric %>% as.factor)  %>%
    dplyr::filter(drug.condi.quant == "Untreated")


tab.combined.v  <-  tab.combined %>%
    dplyr::select(-pairs, -logs, -pictures_included) %>%
    drop_na 


lm.logv <- lm(logv ~        pairv + geno*acidosis.condition,
              data=tab.combined.v )


a1  <-   list(geno="mN10",  acidosis.condition = "ctrl", pairv = "1")
b1  <-   list(geno="mP4",  acidosis.condition = "ctrl", pairv = "1")
c1  <-   list(geno="p53KO",  acidosis.condition = "ctrl", pairv = "1")
d1  <-   list(geno="R273H",  acidosis.condition = "ctrl", pairv = "1")

print(contrast(lm.logv, a1, d1 ))
print(contrast(lm.logv, b1, d1 ))
print(contrast(lm.logv, c1, d1 ))
print(contrast(lm.logv, b1, c1 ))
print(contrast(lm.logv, b1, d1 ))
print(contrast(lm.logv, c1, d1 ))



##### remove batch effect
lm.logv.test <- lm(logv ~   pairv + geno*acidosis.condition,
                   data=tab.combined.v)
a <- predict(lm.logv.test, tab.combined.v)
residuals  <-  tab.combined.v$logv -a 
lm.logv.test$coefficients[grepl("pair", lm.logv.test$coefficients %>% names)]  <- 0
lm.logv.test$coefficients[1] <- 11
a <- predict(lm.logv.test,tab.combined.v) + residuals
tab.combined.v <- tab.combined.v %>% mutate(logv.nobatch = a)
#######################




c.mN10 <-   list(geno="mN10",  acidosis.condition = "ctrl", pairv = "1")
aa.mN10 <-   list(geno="mN10",  acidosis.condition = "aa", pairv = "1")
bb.mN10 <-   list(geno="mN10",  acidosis.condition = "bb", pairv = "1")
c.mP4 <-   list(geno="mP4",  acidosis.condition = "ctrl", pairv = "1")
aa.mP4 <-   list(geno="mP4",  acidosis.condition = "aa", pairv = "1")
bb.mP4 <-   list(geno="mP4",  acidosis.condition = "bb", pairv = "1")
c.p53KO <-   list(geno="p53KO",  acidosis.condition = "ctrl", pairv = "1")
aa.p53KO <-   list(geno="p53KO",  acidosis.condition = "aa", pairv = "1")
bb.p53KO <-   list(geno="p53KO",  acidosis.condition = "bb", pairv = "1")
c.R273H <-   list(geno="R273H",  acidosis.condition = "ctrl", pairv = "1")
aa.R273H <-   list(geno="R273H",  acidosis.condition = "aa", pairv = "1")
bb.R273H <-   list(geno="R273H",  acidosis.condition = "bb", pairv = "1")



## contrasts vs ctrl
cont.logv.aa.mN10  <- contrast(lm.logv, aa.mN10, c.mN10)
print(cont.logv.aa.mN10)
cont.logv.bb.mN10  <- contrast(lm.logv, bb.mN10, c.mN10)
print(cont.logv.bb.mN10)
cont.logv.aa.mP4  <- contrast(lm.logv, aa.mP4, c.mP4)
print(cont.logv.aa.mP4)
cont.logv.bb.mP4  <- contrast(lm.logv, bb.mP4, c.mP4)
print(cont.logv.bb.mP4)
cont.logv.aa.p53KO  <- contrast(lm.logv, aa.p53KO, c.p53KO)
print(cont.logv.aa.p53KO)
cont.logv.bb.p53KO  <- contrast(lm.logv, bb.p53KO, c.p53KO)
print(cont.logv.bb.p53KO)
cont.logv.aa.R273H  <- contrast(lm.logv, aa.R273H, c.R273H)
print(cont.logv.aa.R273H)
cont.logv.bb.R273H  <- contrast(lm.logv, bb.R273H, c.R273H)
print(cont.logv.bb.R273H)

## contrasts bb vs aa
cont.logv.aabb.mN10  <- contrast(lm.logv, aa.mN10, bb.mN10)
print(cont.logv.aabb.mN10)
cont.logv.aabb.mP4  <- contrast(lm.logv, aa.mP4, bb.mP4)
print(cont.logv.aabb.mP4)
cont.logv.aabb.p53KO  <- contrast(lm.logv, aa.p53KO, bb.p53KO)
print(cont.logv.aabb.p53KO)
cont.logv.aabb.R273H  <- contrast(lm.logv, aa.R273H,bb.R273H)
print(cont.logv.aabb.R273H)

## overall pvalues
tab.combined.v2 <- tab.combined.v %>% 
    mutate(acidosis.condition.mN10null = acidosis.condition %>%
               as.character %>%
               ifelse(geno=="mN10", "null", .) %>%
               factor(levels=c("ctrl","aa","bb","null"))) %>%
    mutate(acidosis.condition.mP4null = acidosis.condition %>%
               as.character %>%
               ifelse(geno=="mP4", "null", .) %>%
               factor(levels=c("ctrl","aa","bb","null"))) %>% 
    mutate(acidosis.condition.p53KOnull = acidosis.condition %>%
               as.character %>%
               ifelse(geno=="p53KO", "null", .) %>%
               factor(levels=c("ctrl","aa","bb","null"))) %>% 
    mutate(acidosis.condition.R273Hnull = acidosis.condition %>%
               as.character %>%
               ifelse(geno=="R273H", "null", .) %>%
               factor(levels=c("ctrl","aa","bb","null")))



lm.logv <- lm(logv ~        pairv + geno*acidosis.condition,
              data=tab.combined.v )

lm.logv.mN10null <- lm(logv ~     pairv + geno*acidosis.condition.mN10null,
              data=tab.combined.v2 )
lm.logv.mP4null <- lm(logv ~    pairv + geno*acidosis.condition.mP4null,
              data=tab.combined.v2 )
lm.logv.p53KOnull <- lm(logv ~     pairv + geno*acidosis.condition.p53KOnull,
              data=tab.combined.v2 )
lm.logv.R273Hnull <- lm(logv ~     pairv + geno*acidosis.condition.R273Hnull,
              data=tab.combined.v2 )


anova(lm.logv, lm.logv.mN10null)
anova(lm.logv, lm.logv.mP4null)
anova(lm.logv, lm.logv.p53KOnull)
anova(lm.logv, lm.logv.R273Hnull)


a <- anova(lm.logv, lm.logv.mN10null)
b <- anova(lm.logv, lm.logv.mP4null)
c <- anova(lm.logv, lm.logv.p53KOnull)
d <- anova(lm.logv, lm.logv.R273Hnull)


### plot stuff
stat.test.cbb <- tab.combined.v %>%
    group_by( geno) %>%
    dplyr::filter(acidosis.condition!="aa") %>% 
    wilcox_test(data= . , logv ~ acidosis.condition) 
stat.test.caa <- tab.combined.v %>%
    group_by( geno) %>%
    dplyr::filter(acidosis.condition!="bb") %>% 
    wilcox_test(data= . , logv ~ acidosis.condition) 
stat.test.aabb <- tab.combined.v %>%
    group_by( geno) %>%
    dplyr::filter(acidosis.condition!="ctrl") %>% 
    wilcox_test(data= . , logv ~ acidosis.condition) 
stats <- rbind.data.frame(stat.test.aabb, stat.test.cbb, stat.test.caa)


pvals <- list(cont.logv.aa.mN10, cont.logv.bb.mN10, cont.logv.aa.mP4, cont.logv.bb.mP4, cont.logv.aa.p53KO, cont.logv.bb.p53KO, cont.logv.aa.R273H, cont.logv.bb.R273H, cont.logv.aabb.mN10, cont.logv.aabb.mP4, cont.logv.aabb.p53KO, cont.logv.aabb.R273H) %>%
    lapply(getElement, "Pvalue") %>%
    do.call(rbind.data.frame, .) %>%
    setNames("p") %>%
    mutate(p = formatC(p, format ="e", digit=2)) %>% 
    mutate(geno=c("mN10","mN10","mP4","mP4","p53KO","p53KO","R273H","R273H","mN10","mP4","p53KO","R273H" )) %>%
    mutate(group2=c("aa","bb","aa","bb","aa","bb","aa","bb", "bb", "bb", "bb", "bb"),
           group1=c("ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl", "aa", "aa", "aa", "aa"))


tab.anova = list(a,b,c,d) %>%
    lapply(getElement, "Pr(>F)") %>%
    sapply("[", 2) %>%
    data.frame("pval"=.) %>%
    mutate(geno=c("mN10","mP4","p53KO","R273H")) %>%
    mutate(acidosis.condition="aa", logv=0) %>%
    mutate(p.numeric = as.numeric(pval)) %>% 
    mutate(pval=paste0("pval anova = ", formatC(pval, format ="e", digit=2))) %>%
    mutate(pval=ifelse(p.numeric < 0.05, paste0(pval,"*"), pval))

stat.test <-  stats %>%
    dplyr::select(-p) %>%
    left_join(pvals) %>%
    mutate(p.numeric = as.numeric(p)) %>%
    mutate(p = ifelse(p.numeric < 0.05, paste0(p,"*"), p)) %>% 
    group_by(geno) %>%
    arrange(geno, p.numeric) %>% 
    mutate(y.position = 12.5 + 0.3 * c(0:2)) %>%
    ungroup %>%
    dplyr::filter(p.numeric < 0.1)


### 2nd gen plot
g <- ggboxplot(tab.combined.v,x = "acidosis.condition", y="logv", color = "acidosis.condition", 
               add.params = list(size = 1, shape = 1),
               add = "jitter",
               facet.by="geno", nrow=1, ncol = 4) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, 
                                     vjust = 1, 
                                     hjust=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
        ## stat_summary(fun.y=mean, geom="crossbar", linetype ="dotted", width = 0.75) +
    scale_color_manual(values=c("#8A8F9A", "#FFB700", "#E8217B")) +
    stat_pvalue_manual(stat.test %>% dplyr::filter(p.numeric < 0.05 ),
                       size = 3,label = "p", tip.length=0,
                       fontface =  "bold")  +
    stat_pvalue_manual(stat.test %>% dplyr::filter(p.numeric > 0.05 ),
                       size = 3,label = "p", tip.length=0) +
    geom_text(data=tab.anova, aes(label=pval, x=1.7, y=8.7, fontface = ifelse(p.numeric < 0.05, "bold", "plain")), size=3)
ggpubr::ggexport(g, filename = file.path(out_dir, "boxplots_logv_fig1_2ndGen.pdf"), width = 10, height = 4)


### 2nd gen plot
g <- ggboxplot(tab.combined.v,x = "acidosis.condition", y="logv.nobatch", color = "acidosis.condition", 
               add.params = list(size = 1, shape = 1),
               add = "jitter",
               facet.by="geno", nrow=1, ncol = 4) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, 
                                     vjust = 1, 
                                     hjust=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_color_manual(values=c("#8A8F9A", "#FFB700", "#E8217B")) +
    stat_pvalue_manual(stat.test %>% dplyr::filter(p.numeric < 0.05 ),
                       size = 3,label = "p",
                       fontface =  "bold")  +
    stat_pvalue_manual(stat.test %>% dplyr::filter(p.numeric > 0.05 ),
                       size = 3,label = "p") +
    geom_text(data=tab.anova, aes(label=pval, x=1.7, y=8.7, fontface = ifelse(p.numeric < 0.05, "bold", "plain")), size=3)
ggpubr::ggexport(g, filename = file.path(out_dir, "boxplots_logvBC_fig1_2ndGen.pdf"), width = 10, height = 4)



##################
#### same for area
##################


tab.combined.s <- tab.combined %>% dplyr::select(-pairv, -logv) %>%
    drop_na

lm.logs <- lm(logs ~     pairs +  geno*acidosis.condition,
              data=tab.combined.s )

c.mN10 <-   list(geno="mN10",  acidosis.condition = "ctrl", pairs = "1")
aa.mN10 <-   list(geno="mN10",  acidosis.condition = "aa", pairs = "1")
bb.mN10 <-   list(geno="mN10",  acidosis.condition = "bb", pairs = "1")
c.mP4 <-   list(geno="mP4",  acidosis.condition = "ctrl", pairs = "1")
aa.mP4 <-   list(geno="mP4",  acidosis.condition = "aa", pairs = "1")
bb.mP4 <-   list(geno="mP4",  acidosis.condition = "bb", pairs = "1")
c.p53KO <-   list(geno="p53KO",  acidosis.condition = "ctrl", pairs = "1")
aa.p53KO <-   list(geno="p53KO",  acidosis.condition = "aa", pairs = "1")
bb.p53KO <-   list(geno="p53KO",  acidosis.condition = "bb", pairs = "1")
c.R273H <-   list(geno="R273H",  acidosis.condition = "ctrl", pairs = "1")
aa.R273H <-   list(geno="R273H",  acidosis.condition = "aa", pairs = "1")
bb.R273H <-   list(geno="R273H",  acidosis.condition = "bb", pairs = "1")



##### remove batch effect
lm.logs.test <- lm(logs ~   pairs + geno*acidosis.condition,
                   data=tab.combined.s)
a <- predict(lm.logs.test, tab.combined.s)
residuals  <-  tab.combined.s$logs -a 
lm.logs.test$coefficients[grepl("pair", lm.logs.test$coefficients %>% names)]  <- 0
lm.logs.test$coefficients[1] <- 6.5
a <- predict(lm.logs.test,tab.combined.s) + residuals
tab.combined.s <- tab.combined.s %>% mutate(logs.nobatch = a)
#########################"


a1  <-   list(geno="mN10",  acidosis.condition = "ctrl", pairs = "1")
b1  <-   list(geno="mP4",  acidosis.condition = "ctrl", pairs = "1")
c1  <-   list(geno="p53KO",  acidosis.condition = "ctrl", pairs = "1")
d1  <-   list(geno="R273H",  acidosis.condition = "ctrl", pairs = "1")

print(contrast(lm.logs, a1, d1 ))
print(contrast(lm.logs, b1, d1 ))
print(contrast(lm.logs, c1, d1 ))
print(contrast(lm.logs, b1, c1 ))
print(contrast(lm.logs, b1, d1 ))
print(contrast(lm.logs, c1, d1 ))





## contrasts vs ctrl
cont.logs.aa.mN10  <- contrast(lm.logs, aa.mN10, c.mN10)
print(cont.logs.aa.mN10)
cont.logs.bb.mN10  <- contrast(lm.logs, bb.mN10, c.mN10)
print(cont.logs.bb.mN10)
cont.logs.aa.mP4  <- contrast(lm.logs, aa.mP4, c.mP4)
print(cont.logs.aa.mP4)
cont.logs.bb.mP4  <- contrast(lm.logs, bb.mP4, c.mP4)
print(cont.logs.bb.mP4)
cont.logs.aa.p53KO  <- contrast(lm.logs, aa.p53KO, c.p53KO)
print(cont.logs.aa.p53KO)
cont.logs.bb.p53KO  <- contrast(lm.logs, bb.p53KO, c.p53KO)
print(cont.logs.bb.p53KO)
cont.logs.aa.R273H  <- contrast(lm.logs, aa.R273H, c.R273H)
print(cont.logs.aa.R273H)
cont.logs.bb.R273H  <- contrast(lm.logs, bb.R273H, c.R273H)
print(cont.logs.bb.R273H)

## contrasts bb vs aa
cont.logs.aabb.mN10  <- contrast(lm.logs, aa.mN10, bb.mN10)
print(cont.logs.aabb.mN10)
cont.logs.aabb.mP4  <- contrast(lm.logs, aa.mP4, bb.mP4)
print(cont.logs.aabb.mP4)
cont.logs.aabb.p53KO  <- contrast(lm.logs, aa.p53KO, bb.p53KO)
print(cont.logs.aabb.p53KO)
cont.logs.aabb.R273H  <- contrast(lm.logs, aa.R273H,bb.R273H)
print(cont.logs.aabb.R273H)

## overall pvalues
tab.combined.s2 <- tab.combined.s %>% 
    mutate(acidosis.condition.mN10null = acidosis.condition %>%
               as.character %>%
               ifelse(geno=="mN10", "null", .) %>%
               factor(levels=c("ctrl","aa","bb","null"))) %>%
    mutate(acidosis.condition.mP4null = acidosis.condition %>%
               as.character %>%
               ifelse(geno=="mP4", "null", .) %>%
               factor(levels=c("ctrl","aa","bb","null"))) %>% 
    mutate(acidosis.condition.p53KOnull = acidosis.condition %>%
               as.character %>%
               ifelse(geno=="p53KO", "null", .) %>%
               factor(levels=c("ctrl","aa","bb","null"))) %>% 
    mutate(acidosis.condition.R273Hnull = acidosis.condition %>%
               as.character %>%
               ifelse(geno=="R273H", "null", .) %>%
               factor(levels=c("ctrl","aa","bb","null")))



lm.logs.mN10null <- lm(logs ~  pairs +  geno*acidosis.condition.mN10null,
              data=tab.combined.s2 )
lm.logs.mP4null <- lm(logs ~  pairs +  geno*acidosis.condition.mP4null,
              data=tab.combined.s2 )
lm.logs.p53KOnull <- lm(logs ~  pairs +  geno*acidosis.condition.p53KOnull,
              data=tab.combined.s2 )
lm.logs.R273Hnull <- lm(logs ~  pairs +  geno*acidosis.condition.R273Hnull,
              data=tab.combined.s2 )


anova(lm.logs, lm.logs.mN10null)
anova(lm.logs, lm.logs.mP4null)
anova(lm.logs, lm.logs.p53KOnull)
anova(lm.logs, lm.logs.R273Hnull)




a <- anova(lm.logs, lm.logs.mN10null)
b <- anova(lm.logs, lm.logs.mP4null)
c <- anova(lm.logs, lm.logs.p53KOnull)
d <- anova(lm.logs, lm.logs.R273Hnull)

### plot stuff
### compile pvals

pvals <- list(cont.logs.aa.mN10, cont.logs.bb.mN10, cont.logs.aa.mP4, cont.logs.bb.mP4, cont.logs.aa.p53KO, cont.logs.bb.p53KO, cont.logs.aa.R273H, cont.logs.bb.R273H, cont.logs.aabb.mN10, cont.logs.aabb.mP4, cont.logs.aabb.p53KO, cont.logs.aabb.R273H) %>%
    lapply(getElement, "Pvalue") %>%
    do.call(rbind.data.frame, .) %>%
    setNames("p") %>%
    mutate(p = formatC(p, format ="e", digit=2)) %>% 
    mutate(geno=c("mN10","mN10","mP4","mP4","p53KO","p53KO","R273H","R273H","mN10","mP4","p53KO","R273H" )) %>%
    mutate(group2=c("aa","bb","aa","bb","aa","bb","aa","bb", "bb", "bb", "bb", "bb"),
           group1=c("ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl", "aa", "aa", "aa", "aa"))


tab.anova = list(a,b,c,d) %>%
    lapply(getElement, "Pr(>F)") %>%
    sapply("[", 2) %>%
    data.frame("pval"=.) %>%
    mutate(geno=c("mN10","mP4","p53KO","R273H")) %>%
    mutate(acidosis.condition="aa", logs=0) %>%
    mutate(p.numeric = as.numeric(pval)) %>% 
    mutate(pval=paste0("pval anova = ", formatC(pval, format ="e", digit=2))) %>%
    mutate(pval=ifelse(p.numeric < 0.05, paste0(pval,"*"), pval))
    

stat.test <-  stats %>%
    dplyr::select(-p) %>%
    left_join(pvals) %>%
    mutate(p.numeric = as.numeric(p)) %>%
    mutate(p = ifelse(p.numeric < 0.05, paste0(p,"*"), p)) %>% 
    group_by(geno) %>%
    arrange(geno, p.numeric) %>% 
    mutate(y.position = 7.3 + 0.15 * c(0:2)) %>%
    ungroup %>%
    dplyr::filter(p.numeric < 0.1)


### 2nd gen plot
g <- ggboxplot(tab.combined.s,x = "acidosis.condition", y="logs", color = "acidosis.condition", 
               add.params = list(size = 1, shape = 1),
               add = "jitter",
               facet.by="geno", nrow=1, ncol = 4) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, 
                                     vjust = 1, 
                                     hjust=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_color_manual(values=c("#8A8F9A", "#FFB700", "#E8217B")) +
    stat_pvalue_manual(stat.test %>% dplyr::filter(p.numeric < 0.05 ),
                       size = 3,label = "p", tip.length=0,
                       fontface =  "bold")  +
    stat_pvalue_manual(stat.test %>% dplyr::filter(p.numeric > 0.05 ),
                       size = 3,label = "p", tip.length=0) +
    geom_text(data=tab.anova, aes(label=pval, x=1.7, y=5.5, fontface = ifelse(p.numeric < 0.05, "bold", "plain")), size=3)
ggpubr::ggexport(g, filename = file.path(out_dir, "boxplots_logs_fig1_2ndGen.pdf"), width = 10, height = 4)


### 2nd gen plot
g <- ggboxplot(tab.combined.s,x = "acidosis.condition", y="logs.nobatch", color = "acidosis.condition", 
               add.params = list(size = 1, shape = 1),
               add = "jitter",
               facet.by="geno", nrow=1, ncol = 4) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, 
                                     vjust = 1, 
                                     hjust=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_color_manual(values=c("#8A8F9A", "#FFB700", "#E8217B")) +
    stat_pvalue_manual(stat.test %>% dplyr::filter(p.numeric < 0.05 ),
                       size = 3,label = "p", tip.length=0,
                       fontface =  "bold")  +
    stat_pvalue_manual(stat.test %>% dplyr::filter(p.numeric > 0.05 ),
                       size = 3,label = "p", tip.length=0) +
    geom_text(data=tab.anova, aes(label=pval, x=1.7, y=5.5, fontface = ifelse(p.numeric < 0.05, "bold", "plain")), size=3)
ggpubr::ggexport(g, filename = file.path(out_dir, "boxplots_logsBC_fig1_2ndGen.pdf"), width = 10, height = 4)


