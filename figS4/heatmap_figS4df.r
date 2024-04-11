aa.g.ups <- aa.g %>% filter(w11_aa > 0.25, G > 0.25)
aa.g.ups.genes <- aa.g.ups$gene_id
aa.e.ups <- aa.e %>% filter(w11_aa > 0.25, E > 0.25)
aa.e.ups.genes <- aa.e.ups$gene_id
aa.ge.ups <- aa.ge %>% filter(w11_aa > 0.25, EG > 0.25)
aa.ge.ups.genes <- aa.ge.ups$gene_id

bb.g.ups <- bb.g %>% filter(w11_bb > 0.25, G > 0.25)
bb.g.ups.genes <- bb.g.ups$gene_id
bb.e.ups <- bb.e %>% filter(w11_bb > 0.25, E > 0.25)
bb.e.ups.genes <- bb.e.ups$gene_id
bb.ge.ups <- bb.ge %>% filter(w11_bb > 0.25, EG > 0.25)
bb.ge.ups.genes <- bb.ge.ups$gene_id


up_at_least_two <- c(Reduce(intersect, list(aa.g.ups.genes, aa.e.ups.genes)),
                        Reduce(intersect, list(aa.g.ups.genes, aa.ge.ups.genes)),
                        Reduce(intersect, list(aa.e.ups.genes, aa.ge.ups.genes)))


aa.acuteg.ups <- aa.acuteg %>% filter(w11_aa > 0.25, G_acute > 0.25)
aa.acuteg.ups.genes <- aa.acuteg.ups$gene_id
aa.acutee.ups <- aa.acutee %>% filter(w11_aa > 0.25, E_acute > 0.25)
aa.acutee.ups.genes <- aa.acutee.ups$gene_id
aa.acutege.ups <- aa.acutege %>% filter(w11_aa > 0.25, EG_acute > 0.25)
aa.acutege.ups.genes <- aa.acutege.ups$gene_id


bb.acuteg.ups <- bb.acuteg %>% filter(w11_bb > 0.25, G_acute > 0.25)
bb.acuteg.ups.genes <- bb.acuteg.ups$gene_id
bb.acutee.ups <- bb.acutee %>% filter(w11_bb > 0.25, E_acute > 0.25)
bb.acutee.ups.genes <- bb.acutee.ups$gene_id
bb.acutege.ups <- bb.acutege %>% filter(w11_bb > 0.25, EG_acute > 0.25)
bb.acutege.ups.genes <- bb.acutege.ups$gene_id



up_at_least_two_acute <- c(Reduce(intersect, list(aa.acuteg.ups.genes, aa.acutee.ups.genes)),
                     Reduce(intersect, list(aa.acuteg.ups.genes, aa.acutege.ups.genes)),
                     Reduce(intersect, list(aa.acutee.ups.genes, aa.acutege.ups.genes)))



selected_up <- c(up_at_least_two, up_at_least_two_acute)


# 
# DA_up_genes_1 <- list(
#   w11A_E = aa.acutee.ups.genes,
#   w11A_G = aa.acuteg.ups.genes,
#   w11A_GE = aa.acutege.ups.genes
# )
# 
# 
# upset(fromList(DA_up_genes_1))
# 
# up_all <- Reduce(intersect, DA_up_genes_1)
# up_at_least_two_fc <- all_fc %>% filter(gene_id %in% up_at_least_two_acute)
# up_at_least_two_fc_heatmap <- up_at_least_two_fc %>% select(-gene_id, -E, -G, -EG) %>% 
#   drop_na() %>% column_to_rownames("entrez_name") 


selected_up_fc <- all_fc %>% filter(gene_id %in% selected_up)
selected_up_fc_heatmap <- selected_up_fc %>% select(-gene_id, -w11_aa, -w11_bb) %>% 
  drop_na() %>% column_to_rownames("entrez_name") 


# myBreaks <- c(seq(min(selected_up_fc_heatmap), 0, length.out=ceiling(paletteLength/2)+1), 
#               seq(0.01, max(selected_up_fc_heatmap), length.out=floor(paletteLength/2)))


my_palette <- colorRampPalette(RColorBrewer::brewer.pal(n=7,name="RdBu"))(100) %>% rev
col_breaks <- c(seq(-3,-0.5,length=30),
                seq(-0.49,-0.01,length=20),
                seq(0.01,0.49,length=20),
                seq(0.5,4,length=30))


pheatmap(selected_up_fc_heatmap, cluster_rows = TRUE, cluster_cols = FALSE, color = my_palette, breaks = col_breaks, show_rownames = FALSE, cutree_rows = 2)
rowMeta_df <- data.frame(Sig = rep("", 930), 
                         stringsAsFactors = F,
                         row.names = rownames(selected_up_fc_heatmap))
up_pathway_genes <- Book2$gene_up
for (gene_v in up_pathway_genes) rowMeta_df[rownames(rowMeta_df) == gene_v, "Sig"] <- gene_v

upheat <- pheatmap(selected_up_fc_heatmap, cluster_rows = TRUE, 
         cluster_cols = FALSE, color = myColor, breaks = myBreaks, show_rownames = TRUE, labels_row = rowMeta_df$Sig, fontsize_row = 3)


selected_up_fc_heatmap_acid <- selected_up_fc %>% select(w11_aa, w11_bb, entrez_name) %>% 
  drop_na() %>% column_to_rownames("entrez_name") 

row.order <- upheat$tree_row$order
selected_up_fc_heatmap_acid_1 <- selected_up_fc_heatmap_acid[row.order, ]

acid_down_heatmap <- pheatmap(selected_up_fc_heatmap_acid_1, cluster_rows = FALSE, cluster_cols = FALSE, 
                              color = my_palette, breaks = col_breaks, show_rownames = FALSE)







upheat <- pheatmap(selected_up_fc_heatmap, cluster_rows = TRUE, cluster_cols = FALSE, color = myColor, breaks = myBreaks, show_rownames = FALSE, cutree_rows = 2)
cuts <- cutree(upheat$tree_row, 3)
upcuts <- data.frame(cuts)
upcuts$genes <- rownames(upcuts)
rownames(upcuts) <- NULL
write.xlsx(as.data.frame(upcuts), file = "upcuts.xlsx", row.names = FALSE)


write.xlsx(as.data.frame(up_at_least_two_fc), file = "acute_up_0.25.xlsx", row.names = FALSE)

counts <- tribble(
  ~treatment, ~upaa, ~upbb, ~downaa, ~downbb,
  "adaptation-E",  nrow(aa.g.ups), nrow(bb.g.ups), nrow(aa.g.downs), nrow(bb.g.downs), 
  "adaptation-G",  nrow(aa.e.ups), nrow(bb.e.ups), nrow(aa.e.downs), nrow(bb.e.downs), 
  "adaptation-EG",  nrow(aa.ge.ups), nrow(bb.ge.ups), nrow(aa.ge.downs), nrow(bb.ge.downs), 
  "acute-E", nrow(aa.acuteg.ups), nrow(bb.acuteg.ups), nrow(aa.acuteg.downs), nrow(bb.acuteg.downs),
  "acute-G",  nrow(aa.acutee.ups), nrow(bb.acutee.ups), nrow(aa.acutee.downs), nrow(bb.acutee.downs),
  "acute-EG",  nrow(aa.acutege.ups), nrow(bb.acutege.ups), nrow(aa.acutege.downs), nrow(bb.acutege.downs)
)


ratio <- counts %>% pivot_longer(-treatment, names_to = "con", values_to = "ratio" ) %>% mutate(ratio = round(ratio/nrow(aa.g),3))
ratio$treatment <- factor(ratio$treatment, levels = c( "acute-EG",  "acute-G","acute-E", "adaptation-EG", "adaptation-G", "adaptation-E"))
ratio$con <- factor(ratio$con, levels = c( "upaa",  "upbb", "downaa", "downbb"))


ggplot(ratio, aes(x = con, y = treatment, fill = ratio)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  geom_text(aes(label = ratio), color = "white", size = 4) +
  scale_fill_continuous(high = "black", low = "grey")
  coord_fixed()


