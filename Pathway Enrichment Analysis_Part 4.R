# Step 4: Pathway Enrichment Analysis

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
install.packages("GEOquery")
library(GEOquery)
# Map probes to Entrez IDs
gpl <- getGEO("GPL570", destdir = ".")
probe_to_entrez <- Table(gpl)[, c("ID", "ENTREZ_GENE_ID")]
probe_to_entrez <- probe_to_entrez[probe_to_entrez$ENTREZ_GENE_ID != "", ]
significant_genes_df <- as.data.frame(significant_genes)
significant_genes_df$probe_id <- rownames(significant_genes_df)
significant_genes_df <- merge(significant_genes_df, probe_to_entrez, by.x = "probe_id", by.y = "ID")
entrez_ids <- significant_genes_df$ENTREZ_GENE_ID

# GO enrichment analysis
go_enrichment <- enrichGO(gene = entrez_ids,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
summary(go_enrichment)
dotplot(go_enrichment, showCategory = 8, title = "GO Enrichment Analysis (Biological Process)") +
  theme(axis.text.y = element_text(size = 10))

# Save GO results
write.csv(as.data.frame(go_enrichment), file = file.path("go_enrichment_results.csv"))


##################  WILL FIX THIS LATER - NOT WORKING ENTRZ IDS NOT LABELLING #################

# KEGG enrichment analysis
# kegg_mapping <- bitr(entrez_ids, 
#                      fromType = "ENTREZID", 
#                      toType = "PATH", 
#                      OrgDb = org.Hs.eg.db)
# valid_kegg_ids <- kegg_mapping$PATH[!is.na(kegg_mapping$PATH)]
# cat("Number of valid KEGG IDs:", length(valid_kegg_ids), "\n")
# if (length(valid_kegg_ids) == 0) {
#   stop("No valid KEGG IDs found.")
# }
# kegg_enrichment <- enrichKEGG(
#   gene = valid_kegg_ids,
#   organism = "hsa",
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.5,
#   qvalueCutoff = 0.5,
# )
# if (nrow(as.data.frame(kegg_enrichment)) > 0) {
#   print(summary(kegg_enrichment))
#   dotplot(kegg_enrichment, 
#           showCategory = 15, 
#           title = "KEGG Pathway Enrichment Analysis") +
#     theme(axis.text.y = element_text(size = 10))
#   write.csv(as.data.frame(kegg_enrichment), file = "kegg_enrichment_results.csv")
# } else {
#   cat("No significant KEGG pathways found.\n")
# }



