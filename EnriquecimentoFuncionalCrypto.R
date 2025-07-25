# Instalar pacotes necessários
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "KEGGREST", "pathview", "enrichplot", "GOSemSim", "dplyr"))

# Carregar bibliotecas
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(GOSemSim)
library(dplyr)
library(pathview)

# Carregar lista de genes (coluna única)
gene_symbols = read.csv('listaCryptococcus.csv', sep = ';', dec = '.', header = TRUE)
genesCrypto = gene_symbols$GeneID

# Enriquecimento KEGG para Cryptococcus neoformans (cne)
kegg_resultCrypto <- enrichKEGG(
  gene          = genesCrypto,
  organism      = "cng",
  keyType       = "kegg",
  pvalueCutoff  = 0.05
)

write.csv(kegg_resultCrypto@result, "kegg_enrichment_resultsCrypto.csv", row.names = FALSE)

head(kegg_resultCrypto)
kegg_resultCrypto

kegg_resultCrypto@result$Description


viasKeggSignificantes <- subset(kegg_result@result, p.adjust < 0.05)[, c("ID", "Description", "Count")]
print(viasKeggSignificantes)

# Visualização KEGG (substitua pela via de interesse)
pathview(
  gene.data   = genes,
  pathway.id  = "cne00500", 
  species     = "cne",
  gene.idtype = "KEGG"
)

# -------- ENRIQUECIMENTO GO MANUAL --------
# É necessário fornecer os arquivos:
# - go_term2gene.csv: colunas GO_ID, GeneID
# - go_term2name.csv: colunas GO_ID, TermDescription

go_gene_map <- read.csv("go_term2gene.csv")
go_term_names <- read.csv("go_term2name.csv")

ego_GO <- enricher(
  gene        = genes,
  TERM2GENE   = go_gene_map,
  TERM2NAME   = go_term_names,
  pvalueCutoff = 0.05
)

# Gráfico de barras
ego_bp <- ego_GO@result %>%
  top_n(20, -pvalue) %>%
  arrange(pvalue) %>%
  select(Description, Count, pvalue)

barPlot_BP <- ggplot(ego_bp, aes(x = reorder(Description, pvalue), y = Count, fill = -log10(pvalue))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 20 Termos GO Enriquecidos", x = "Termos GO", y = "Número de Genes", fill = "-log10(p-value)") +
  theme_minimal()

ggsave("barPlot_Crypto_GO.jpg", plot = barPlot_BP, width = 10, height = 6, dpi = 600)

# Cnetplot
cnetplot_BP <- cnetplot(ego_GO, showCategory = 10, categorySize = "pvalue")
ggsave("cnetplot_Crypto_GO.jpg", plot = cnetplot_BP, width = 8, height = 6, dpi = 600)

# Mapa de Enriquecimento
ego_GO <- pairwise_termsim(ego_GO)
emapplot(ego_GO, showCategory = 20)
