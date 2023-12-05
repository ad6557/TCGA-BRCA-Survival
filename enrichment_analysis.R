library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

load("preprocessing+genefinding.RData")

gene_list_symbol <- select1

### GO Enrichment Analysis ###

ego <- enrichGO(gene = gene_list_symbol,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP")

# GO enrichment plot
# View(ego@result)
plot(barplot(ego,showCategory = 20))
goplot(ego)
# Gene-Concept Network
ego2 <- simplify(ego)
cnetplot(ego2, circular = TRUE)
# Enrichment Map


### KEGG Enrichment Analysis ###

# Convert to Entrez IDs
gene_list_entrez <- mapIds(org.Hs.eg.db, 
                           keys = gene_list_symbol, 
                           column = "ENTREZID", 
                           keytype = "SYMBOL", 
                           multiVals = "first")
kegg <- enrichKEGG(gene = gene_list_entrez,
                   organism = 'hsa',  # 'hsa' for Homo sapiens
                   keyType = "kegg")

# KEGG enrichment plot
# View(kegg@result)
plot(barplot(kegg,showCategory = 20))
