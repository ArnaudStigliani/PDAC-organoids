library(tidyverse)
library(scales)
DEG_ko_aa <- DEG_all %>% filter(contrast == "AA_KO_11") %>% select(gene_id, logFC) %>% rename(w11_aa = logFC) %>% remove_rownames()
DEG_ko_bb <- DEG_all %>% filter(contrast == "BB_KO_11") %>% select(gene_id, logFC) %>% rename(w11_bb = logFC) %>% remove_rownames()

DEG_aa <- DEG_all %>% filter(contrast == "AA_11") %>% select(gene_id, logFC) %>% rename(w11_aa = logFC) %>% remove_rownames()
DEG_bb <- DEG_all %>% filter(contrast == "BB_11") %>% select(gene_id, logFC) %>% rename(w11_bb = logFC) %>% remove_rownames()
DEG_g <- DEG_all %>% filter(contrast == "DA_G_12") %>% select(gene_id, logFC) %>% rename(G = logFC) %>% remove_rownames()
DEG_e <- DEG_all %>% filter(contrast == "DA_E_12") %>% select(gene_id, logFC) %>% rename(E = logFC) %>% remove_rownames()
DEG_ge <- DEG_all %>% filter(contrast == "DA_EG_12") %>% select(gene_id, logFC) %>% rename(EG = logFC) %>% remove_rownames()

DEG_acuteg <- DEG_all %>% filter(contrast == "nDA_G_12") %>% select(gene_id, logFC) %>% rename(G_acute = logFC) %>% remove_rownames()
DEG_acutee <- DEG_all %>% filter(contrast == "nDA_E_12") %>% select(gene_id, logFC) %>% rename(E_acute = logFC) %>% remove_rownames()
DEG_acutege <- DEG_all %>% filter(contrast == "nDA_EG_12") %>% select(gene_id, logFC) %>% rename(EG_acute = logFC) %>% remove_rownames()




aa.g

#mp4
aa.g.1 <- merge(DEG_aa, DEG_g, by = "gene_id")
aa.e.1 <- merge(DEG_aa, DEG_e, by = "gene_id")
aa.ge.1 <- merge(DEG_aa, DEG_ge, by = "gene_id")
bb.g.1 <- merge(DEG_bb, DEG_g, by = "gene_id")
bb.e.1 <- merge(DEG_bb, DEG_e, by = "gene_id")
bb.ge.1 <- merge(DEG_bb, DEG_ge, by = "gene_id")


aa.acuteg.1 <- merge(DEG_aa, DEG_acuteg, by = "gene_id")
aa.acutee.1 <- merge(DEG_aa, DEG_acutee, by = "gene_id")
aa.acutege.1 <- merge(DEG_aa, DEG_acutege, by = "gene_id")
bb.acuteg.1 <- merge(DEG_bb, DEG_acuteg, by = "gene_id")
bb.acutee.1 <- merge(DEG_bb, DEG_acutee, by = "gene_id")
bb.acutege.1 <- merge(DEG_bb, DEG_acutege, by = "gene_id")




# p53ko
# ko.aa.g.1 <- merge(DEG_ko_aa, DEG_g, by = "gene_id")
# ko.aa.e.1 <- merge(DEG_ko_aa, DEG_e, by = "gene_id")
# ko.aa.ge.1 <- merge(DEG_ko_aa, DEG_ge, by = "gene_id")
# ko.bb.g.1 <- merge(DEG_ko_bb, DEG_g, by = "gene_id")
# ko.bb.e.1 <- merge(DEG_ko_bb, DEG_e, by = "gene_id")
# ko.bb.ge.1 <- merge(DEG_ko_bb, DEG_ge, by = "gene_id")
# 
# ko.aa.acuteg.1 <- merge(DEG_ko_aa, DEG_acuteg, by = "gene_id")
# ko.aa.acutee.1 <- merge(DEG_ko_aa, DEG_acutee, by = "gene_id")
# ko.aa.acutege.1 <- merge(DEG_ko_aa, DEG_acutege, by = "gene_id")
# ko.bb.acuteg.1 <- merge(DEG_ko_bb, DEG_acuteg, by = "gene_id")
# ko.bb.acutee.1 <- merge(DEG_ko_bb, DEG_acutee, by = "gene_id")
# ko.bb.acutege.1 <- merge(DEG_ko_bb, DEG_acutege, by = "gene_id")
# 
# aa.g.1 <- ko.aa.g.1
# aa.e.1 <- ko.aa.e.1
# aa.ge.1 <- ko.aa.ge.1
# bb.g.1 <- ko.bb.g.1
# bb.e.1 <- ko.bb.e.1
# bb.ge.1 <- ko.bb.ge.1
# 
# aa.acuteg.1 <- ko.aa.acuteg.1
# aa.acutee.1 <- ko.aa.acutee.1
# aa.acutege.1 <- ko.aa.acutege.1
# bb.acuteg.1 <- ko.bb.acuteg.1
# bb.acutee.1 <- ko.bb.acutee.1
# bb.acutege.1 <- ko.bb.acutege.1




write.xlsx(as.data.frame(aa.g.1 %>% filter(w11_aa < -0.25, G < -0.25)), file = "ko-aa-g-shared-down.xlsx", row.names = FALSE)
write.xlsx(as.data.frame(aa.e.1 %>% filter(w11_aa < -0.25, E < -0.25)), file = "mp4-aa-e-shared-down.xlsx", row.names = FALSE)
write.xlsx(as.data.frame(aa.ge.1 %>% filter(w11_aa < -0.25, EG < -0.25)), file = "mp4-aa-ge-shared-down.xlsx", row.names = FALSE)
write.xlsx(as.data.frame(bb.g.1 %>% filter(w11_bb < -0.25, G < -0.25)), file = "mp4-bb-g-shared-down.xlsx", row.names = FALSE)
write.xlsx(as.data.frame(bb.e.1 %>% filter(w11_bb < -0.25, E < -0.25)), file = "mp4-bb-e-shared-down.xlsx", row.names = FALSE)
write.xlsx(as.data.frame(bb.ge.1 %>% filter(w11_bb < -0.25, EG < -0.25)), file = "mp4-bb-ge-shared-down.xlsx", row.names = FALSE)
write.xlsx(as.data.frame(aa.acuteg.1 %>% filter(w11_aa < -0.25, G_acute < -0.25)), file = "mp4-aa-g-acute-shared-down.xlsx", row.names = FALSE)
write.xlsx(as.data.frame(aa.acutee.1 %>% filter(w11_aa < -0.25, E_acute < -0.25)), file = "mp4-aa-e-acute-shared-down.xlsx", row.names = FALSE)
write.xlsx(as.data.frame(aa.acutege.1 %>% filter(w11_aa < -0.25, EG_acute < -0.25)), file = "mp4-aa-ge-acute-shared-down.xlsx", row.names = FALSE)
write.xlsx(as.data.frame(bb.acuteg.1 %>% filter(w11_bb < -0.25, G_acute < -0.25)), file = "mp4-bb-g-acute-shared-down.xlsx", row.names = FALSE)
write.xlsx(as.data.frame(bb.acutee.1 %>% filter(w11_bb < -0.25, E_acute < -0.25)), file = "mp4-bb-e-acute-shared-down.xlsx", row.names = FALSE)
write.xlsx(as.data.frame(bb.acutege.1 %>% filter(w11_bb < -0.25, EG_acute < -0.25)), file = "mp4-bb-ge-acute-shared-down.xlsx", row.names = FALSE)





#permutation
perm = 5000

aa.g.up.ratios <- c()
aa.e.up.ratios <- c()
aa.ge.up.ratios <- c()
bb.g.up.ratios <- c()
bb.e.up.ratios <- c()
bb.ge.up.ratios <- c()

aa.g.acute.up.ratios <- c()
aa.e.acute.up.ratios <- c()
aa.ge.acute.up.ratios <- c()
bb.g.acute.up.ratios <- c()
bb.e.acute.up.ratios <- c()
bb.ge.acute.up.ratios <- c()

ngene <- nrow(aa.g.1)

for (i in 1:perm){
  #adaptation
  aa.g.p1 <- aa.g.1 %>% mutate(w11_aa_1 = sample(w11_aa), G_1 = sample(G))
  nratio <- (aa.g.p1 %>% filter(w11_aa_1 > 0.25, G_1 > 0.25) %>% nrow())/ngene
  aa.g.up.ratios <- append(aa.g.up.ratios, nratio)
}

#hist(aa.g.up.ratios)

for (i in 1:perm){  
  aa.e.p1 <- aa.e.1 %>% mutate(w11_aa_1 = sample(w11_aa), E_1 = sample(E))
  nratio <- (aa.e.p1 %>% filter(w11_aa_1 > 0.25, E_1 > 0.25) %>% nrow())/ngene
  aa.e.up.ratios <- append(aa.e.up.ratios, nratio)
}

#hist(aa.e.up.ratios)

for (i in 1:perm){  
  aa.ge.p1 <- aa.ge.1 %>% mutate(w11_aa_1 = sample(w11_aa), EG_1 = sample(EG))
  nratio <- (aa.ge.p1 %>% filter(w11_aa_1 > 0.25, EG_1 > 0.25) %>% nrow())/ngene
  aa.ge.up.ratios <- append(aa.ge.up.ratios, nratio)
}

hist(aa.ge.up.ratios)


for (i in 1:perm){
  bb.g.p1 <- bb.g.1 %>% mutate(w11_bb_1 = sample(w11_bb), G_1 = sample(G))
  nratio <- (bb.g.p1 %>% filter(w11_bb_1 > 0.25, G_1 > 0.25) %>% nrow())/ngene
  bb.g.up.ratios <- append(bb.g.up.ratios, nratio)
}

for (i in 1:perm){
  bb.e.p1 <- bb.e.1 %>% mutate(w11_bb_1 = sample(w11_bb), E_1 = sample(E))
  nratio <- (bb.e.p1 %>% filter(w11_bb_1 > 0.25, E_1 > 0.25) %>% nrow())/ngene
  bb.e.up.ratios <- append(bb.e.up.ratios, nratio)
}

for (i in 1:perm){
  bb.ge.p1 <- bb.ge.1 %>% mutate(w11_bb_1 = sample(w11_bb), EG_1 = sample(EG))
  nratio <- (bb.ge.p1 %>% filter(w11_bb_1 > 0.25, EG_1 > 0.25) %>% nrow())/ngene
  bb.ge.up.ratios <- append(bb.ge.up.ratios, nratio)
}

#acute
for (i in 1:perm){
  aa.g.acute.p1 <- aa.acuteg.1 %>% mutate(w11_aa_1 = sample(w11_aa), G_1 = sample(G_acute))
  nratio <- (aa.g.acute.p1 %>% filter(w11_aa_1 > 0.25, G_1 > 0.25) %>% nrow())/ngene
  aa.g.acute.up.ratios <- append(aa.g.acute.up.ratios, nratio)
}

for (i in 1:perm){
  aa.e.acute.p1 <- aa.acutee.1 %>% mutate(w11_aa_1 = sample(w11_aa), E_1 = sample(E_acute))
  nratio <- (aa.e.acute.p1 %>% filter(w11_aa_1 > 0.25, E_1 > 0.25) %>% nrow())/ngene
  aa.e.acute.up.ratios <- append(aa.e.acute.up.ratios, nratio)
}

for (i in 1:perm){
  aa.ge.acute.p1 <- aa.acutege.1 %>% mutate(w11_aa_1 = sample(w11_aa), EG_1 = sample(EG_acute))
  nratio <- (aa.ge.acute.p1 %>% filter(w11_aa_1 > 0.25, EG_1 > 0.25) %>% nrow())/ngene
  aa.ge.acute.up.ratios <- append(aa.ge.acute.up.ratios, nratio)
}

for (i in 1:perm){
  bb.g.acute.p1 <- bb.acuteg.1 %>% mutate(w11_bb_1 = sample(w11_bb), G_1 = sample(G_acute))
  nratio <- (bb.g.acute.p1 %>% filter(w11_bb_1 > 0.25, G_1 > 0.25) %>% nrow())/ngene
  bb.g.acute.up.ratios <- append(bb.g.acute.up.ratios, nratio)
}

for (i in 1:perm){
  bb.e.acute.p1 <- bb.acutee.1 %>% mutate(w11_bb_1 = sample(w11_bb), E_1 = sample(E_acute))
  nratio <- (bb.e.acute.p1 %>% filter(w11_bb_1 > 0.25, E_1 > 0.25) %>% nrow())/ngene
  bb.e.acute.up.ratios <- append(bb.e.acute.up.ratios, nratio)
}

for (i in 1:perm){ 
  bb.ge.acute.p1 <- bb.acutege.1 %>% mutate(w11_bb_1 = sample(w11_bb), EG_1 = sample(EG_acute))
  nratio <- (bb.ge.acute.p1 %>% filter(w11_bb_1 > 0.25, EG_1 > 0.25) %>% nrow())/ngene
  bb.ge.acute.up.ratios <- append(bb.ge.acute.up.ratios, nratio)
}

aa.g.up.ob.ratio <- (aa.g.1 %>% filter(w11_aa > 0.25, G > 0.25) %>% nrow())/ngene
aa.e.up.ob.ratio <- (aa.e.1 %>% filter(w11_aa > 0.25, E > 0.25) %>% nrow())/ngene
aa.ge.up.ob.ratio <- (aa.ge.1 %>% filter(w11_aa > 0.25, EG > 0.25) %>% nrow())/ngene
bb.g.up.ob.ratio <- (bb.g.1 %>% filter(w11_bb > 0.25, G > 0.25) %>% nrow())/ngene
bb.e.up.ob.ratio <- (bb.e.1 %>% filter(w11_bb > 0.25, E > 0.25) %>% nrow())/ngene
bb.ge.up.ob.ratio <- (bb.ge.1 %>% filter(w11_bb > 0.25, EG > 0.25) %>% nrow())/ngene

aa.g.acute.up.ob.ratio <- (aa.acuteg.1 %>% filter(w11_aa > 0.25, G_acute > 0.25) %>% nrow())/ngene
aa.e.acute.up.ob.ratio <- (aa.acutee.1 %>% filter(w11_aa > 0.25, E_acute > 0.25) %>% nrow())/ngene
aa.ge.acute.up.ob.ratio <- (aa.acutege.1 %>% filter(w11_aa > 0.25, EG_acute > 0.25) %>% nrow())/ngene
bb.g.acute.up.ob.ratio <- (bb.acuteg.1 %>% filter(w11_bb > 0.25, G_acute > 0.25) %>% nrow())/ngene
bb.e.acute.up.ob.ratio <- (bb.acutee.1 %>% filter(w11_bb > 0.25, E_acute > 0.25) %>% nrow())/ngene
bb.ge.acute.up.ob.ratio <- (bb.acutege.1 %>% filter(w11_bb > 0.25, EG_acute > 0.25) %>% nrow())/ngene


write.xlsx(as.data.frame(bb.acutege.1 %>% filter(w11_bb > 0.25, EG_acute > 0.25)), file = "mp4-bb-ge-acute-shared-up.xlsx", row.names = FALSE)


#----------------------down------------------

aa.g.down.ratios <- c()
aa.e.down.ratios <- c()
aa.ge.down.ratios <- c()
bb.g.down.ratios <- c()
bb.e.down.ratios <- c()
bb.ge.down.ratios <- c()

aa.g.acute.down.ratios <- c()
aa.e.acute.down.ratios <- c()
aa.ge.acute.down.ratios <- c()
bb.g.acute.down.ratios <- c()
bb.e.acute.down.ratios <- c()
bb.ge.acute.down.ratios <- c()

ngene <- nrow(aa.g.1)

for (i in 1:perm){
  #adaptation
  aa.g.p1 <- aa.g.1 %>% mutate(w11_aa_1 = sample(w11_aa), G_1 = sample(G))
  nratio <- (aa.g.p1 %>% filter(w11_aa_1 < -0.25, G_1 < -0.25) %>% nrow())/ngene
  aa.g.down.ratios <- append(aa.g.down.ratios, nratio)
}

#hist(aa.g.down.ratios)

for (i in 1:perm){  
  aa.e.p1 <- aa.e.1 %>% mutate(w11_aa_1 = sample(w11_aa), E_1 = sample(E))
  nratio <- (aa.e.p1 %>% filter(w11_aa_1 < -0.25, E_1 < -0.25) %>% nrow())/ngene
  aa.e.down.ratios <- append(aa.e.down.ratios, nratio)
}

#hist(aa.e.down.ratios)

for (i in 1:perm){  
  aa.ge.p1 <- aa.ge.1 %>% mutate(w11_aa_1 = sample(w11_aa), EG_1 = sample(EG))
  nratio <- (aa.ge.p1 %>% filter(w11_aa_1 < -0.25, EG_1 < -0.25) %>% nrow())/ngene
  aa.ge.down.ratios <- append(aa.ge.down.ratios, nratio)
}

hist(aa.ge.down.ratios)


for (i in 1:perm){
  bb.g.p1 <- bb.g.1 %>% mutate(w11_bb_1 = sample(w11_bb), G_1 = sample(G))
  nratio <- (bb.g.p1 %>% filter(w11_bb_1 < -0.25, G_1 < -0.25) %>% nrow())/ngene
  bb.g.down.ratios <- append(bb.g.down.ratios, nratio)
}

for (i in 1:perm){
  bb.e.p1 <- bb.e.1 %>% mutate(w11_bb_1 = sample(w11_bb), E_1 = sample(E))
  nratio <- (bb.e.p1 %>% filter(w11_bb_1 < -0.25, E_1 < -0.25) %>% nrow())/ngene
  bb.e.down.ratios <- append(bb.e.down.ratios, nratio)
}

for (i in 1:perm){
  bb.ge.p1 <- bb.ge.1 %>% mutate(w11_bb_1 = sample(w11_bb), EG_1 = sample(EG))
  nratio <- (bb.ge.p1 %>% filter(w11_bb_1 < -0.25, EG_1 < -0.25) %>% nrow())/ngene
  bb.ge.down.ratios <- append(bb.ge.down.ratios, nratio)
}

#acute
for (i in 1:perm){
  aa.g.acute.p1 <- aa.acuteg.1 %>% mutate(w11_aa_1 = sample(w11_aa), G_1 = sample(G_acute))
  nratio <- (aa.g.acute.p1 %>% filter(w11_aa_1 < -0.25, G_1 < -0.25) %>% nrow())/ngene
  aa.g.acute.down.ratios <- append(aa.g.acute.down.ratios, nratio)
}

for (i in 1:perm){
  aa.e.acute.p1 <- aa.acutee.1 %>% mutate(w11_aa_1 = sample(w11_aa), E_1 = sample(E_acute))
  nratio <- (aa.e.acute.p1 %>% filter(w11_aa_1 < -0.25, E_1 < -0.25) %>% nrow())/ngene
  aa.e.acute.down.ratios <- append(aa.e.acute.down.ratios, nratio)
}

for (i in 1:perm){
  aa.ge.acute.p1 <- aa.acutege.1 %>% mutate(w11_aa_1 = sample(w11_aa), EG_1 = sample(EG_acute))
  nratio <- (aa.ge.acute.p1 %>% filter(w11_aa_1 < -0.25, EG_1 < -0.25) %>% nrow())/ngene
  aa.ge.acute.down.ratios <- append(aa.ge.acute.down.ratios, nratio)
}

for (i in 1:perm){
  bb.g.acute.p1 <- bb.acuteg.1 %>% mutate(w11_bb_1 = sample(w11_bb), G_1 = sample(G_acute))
  nratio <- (bb.g.acute.p1 %>% filter(w11_bb_1 < -0.25, G_1 < -0.25) %>% nrow())/ngene
  bb.g.acute.down.ratios <- append(bb.g.acute.down.ratios, nratio)
}

for (i in 1:perm){
  bb.e.acute.p1 <- bb.acutee.1 %>% mutate(w11_bb_1 = sample(w11_bb), E_1 = sample(E_acute))
  nratio <- (bb.e.acute.p1 %>% filter(w11_bb_1 < -0.25, E_1 < -0.25) %>% nrow())/ngene
  bb.e.acute.down.ratios <- append(bb.e.acute.down.ratios, nratio)
}

for (i in 1:perm){ 
  bb.ge.acute.p1 <- bb.acutege.1 %>% mutate(w11_bb_1 = sample(w11_bb), EG_1 = sample(EG_acute))
  nratio <- (bb.ge.acute.p1 %>% filter(w11_bb_1 < -0.25, EG_1 < -0.25) %>% nrow())/ngene
  bb.ge.acute.down.ratios <- append(bb.ge.acute.down.ratios, nratio)
}

aa.g.down.ob.ratio <- (aa.g.1 %>% filter(w11_aa < -0.25, G < -0.25) %>% nrow())/ngene
aa.e.down.ob.ratio <- (aa.e.1 %>% filter(w11_aa < -0.25, E < -0.25) %>% nrow())/ngene
aa.ge.down.ob.ratio <- (aa.ge.1 %>% filter(w11_aa < -0.25, EG < -0.25) %>% nrow())/ngene
bb.g.down.ob.ratio <- (bb.g.1 %>% filter(w11_bb < -0.25, G < -0.25) %>% nrow())/ngene
bb.e.down.ob.ratio <- (bb.e.1 %>% filter(w11_bb < -0.25, E < -0.25) %>% nrow())/ngene
bb.ge.down.ob.ratio <- (bb.ge.1 %>% filter(w11_bb < -0.25, EG < -0.25) %>% nrow())/ngene

aa.g.acute.down.ob.ratio <- (aa.acuteg.1 %>% filter(w11_aa < -0.25, G_acute < -0.25) %>% nrow())/ngene
aa.e.acute.down.ob.ratio <- (aa.acutee.1 %>% filter(w11_aa < -0.25, E_acute < -0.25) %>% nrow())/ngene
aa.ge.acute.down.ob.ratio <- (aa.acutege.1 %>% filter(w11_aa < -0.25, EG_acute < -0.25) %>% nrow())/ngene
bb.g.acute.down.ob.ratio <- (bb.acuteg.1 %>% filter(w11_bb < -0.25, G_acute < -0.25) %>% nrow())/ngene
bb.e.acute.down.ob.ratio <- (bb.acutee.1 %>% filter(w11_bb < -0.25, E_acute < -0.25) %>% nrow())/ngene
bb.ge.acute.down.ob.ratio <- (bb.acutege.1 %>% filter(w11_bb < -0.25, EG_acute < -0.25) %>% nrow())/ngene


table(aa.e.up.ob.ratio > aa.e.up.ratios)

#----------summary
perm.res <- tribble(
  ~treatment, ~drug, ~acid, ~shared, ~ratio, ~perm_mean, ~perm_std,
  "adaptation", "G", "aa", "up", aa.g.up.ob.ratio, mean(aa.g.up.ratios), sd(aa.g.up.ratios),
  "adaptation", "E", "aa", "up", aa.e.up.ob.ratio, mean(aa.e.up.ratios), sd(aa.e.up.ratios),
  "adaptation", "GE", "aa", "up", aa.ge.up.ob.ratio, mean(aa.ge.up.ratios), sd(aa.ge.up.ratios), 
  "adaptation", "G", "bb", "up", bb.g.up.ob.ratio, mean(bb.g.up.ratios), sd(bb.g.up.ratios),
  "adaptation", "E", "bb", "up", bb.e.up.ob.ratio, mean(bb.e.up.ratios), sd(bb.e.up.ratios),
  "adaptation", "GE", "bb", "up", bb.ge.up.ob.ratio, mean(bb.ge.up.ratios), sd(bb.ge.up.ratios), 
  "acute", "G", "aa", "up", aa.g.acute.up.ob.ratio, mean(aa.g.acute.up.ratios), sd(aa.g.acute.up.ratios),
  "acute", "E", "aa", "up", aa.e.acute.up.ob.ratio, mean(aa.e.acute.up.ratios), sd(aa.e.acute.up.ratios),
  "acute", "GE", "aa", "up", aa.ge.acute.up.ob.ratio, mean(aa.ge.acute.up.ratios), sd(aa.ge.acute.up.ratios),
  "acute", "G", "bb", "up", bb.g.acute.up.ob.ratio, mean(bb.g.acute.up.ratios), sd(bb.g.acute.up.ratios),
  "acute", "E", "bb", "up", bb.e.acute.up.ob.ratio, mean(bb.e.acute.up.ratios), sd(bb.e.acute.up.ratios),
  "acute", "GE", "bb", "up", bb.ge.acute.up.ob.ratio, mean(bb.ge.acute.up.ratios), sd(bb.ge.acute.up.ratios),
  "adaptation", "G", "aa", "down", aa.g.down.ob.ratio, mean(aa.g.down.ratios), sd(aa.g.down.ratios),
  "adaptation", "E", "aa", "down", aa.e.down.ob.ratio, mean(aa.e.down.ratios), sd(aa.e.down.ratios),
  "adaptation", "GE", "aa", "down", aa.ge.down.ob.ratio, mean(aa.ge.down.ratios), sd(aa.ge.down.ratios), 
  "adaptation", "G", "bb", "down", bb.g.down.ob.ratio, mean(bb.g.down.ratios), sd(bb.g.down.ratios),
  "adaptation", "E", "bb", "down", bb.e.down.ob.ratio, mean(bb.e.down.ratios), sd(bb.e.down.ratios),
  "adaptation", "GE", "bb", "down", bb.ge.down.ob.ratio, mean(bb.ge.down.ratios), sd(bb.ge.down.ratios), 
  "acute", "G", "aa", "down", aa.g.acute.down.ob.ratio, mean(aa.g.acute.down.ratios), sd(aa.g.acute.down.ratios),
  "acute", "E", "aa", "down", aa.e.acute.down.ob.ratio, mean(aa.e.acute.down.ratios), sd(aa.e.acute.down.ratios),
  "acute", "GE", "aa", "down", aa.ge.acute.down.ob.ratio, mean(aa.ge.acute.down.ratios), sd(aa.ge.acute.down.ratios),
  "acute", "G", "bb", "down", bb.g.acute.down.ob.ratio, mean(bb.g.acute.down.ratios), sd(bb.g.acute.down.ratios),
  "acute", "E", "bb", "down", bb.e.acute.down.ob.ratio, mean(bb.e.acute.down.ratios), sd(bb.e.acute.down.ratios),
  "acute", "GE", "bb", "down", bb.ge.acute.down.ob.ratio, mean(bb.ge.acute.down.ratios), sd(bb.ge.acute.down.ratios),
)



perm.res <- perm.res %>% mutate(zscore = round((ratio-perm_mean)/perm_std,2) )
perm.res <- perm.res %>% mutate(log2fc=round(log2(ratio/perm_mean),2) )

write.xlsx(as.data.frame(perm.res), file = "perm_res_mp4.xlsx", row.names = FALSE)

perm.res$treatment <- factor(perm.res$treatment, levels = c("adaptation", "acute"))
perm.res$shared <- factor(perm.res$shared, levels = c("up", "down"))


breaks <- c(min(perm.res$log2fc), 0)
breaks2 <- c(0, max(perm.res$log2fc))

ggplot(perm.res) +
  geom_tile(aes(x = acid, y = drug, fill = log2fc),
            color = "white",
            lwd = 1.5,
            linetype = 1) +
  geom_text(aes(x = acid, y = drug, label = log2fc), color = "black", size = 4) +
  facet_grid(treatment~shared)+
  scale_color_gradient2(limits = c(-0.6, 0.6),
                        low = muted("green"),
                        mid = "white",
                        high = muted("blue"),
                        midpoint = 0,
                        space = "Lab",
                        guide = "colourbar",
                        aesthetics = "fill") +
  coord_fixed() + theme_bw()




