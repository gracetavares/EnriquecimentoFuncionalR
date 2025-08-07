# Instalar pacotes necessários
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
if (!requireNamespace("here", quietly = TRUE)){
  install.packages("here")}
if (!requireNamespace (c("clusterProfiler", "KEGGREST", "pathview", "enrichplot", "GOSemSim", "dplyr",  quietly = TRUE))){
  BiocManager::install (c("clusterProfiler", "KEGGREST", "pathview", "enrichplot", "GOSemSim", "dplyr"))}
if (!requireNamespace("AnnotationForge", quietly = TRUE)) {
  BiocManager::install("AnnotationForge")}
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi")}
if (!requireNamespace("RSQLite", quietly = TRUE)) {
  install.packages("RSQLite")}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")}


# Carregar bibliotecas
library(clusterProfiler) #enriquecimento funcional com KEGG
library(enrichplot) 
library(ggplot2)
library(GOSemSim)
library(dplyr) #manipulação de dados
library(pathview)
library(here) #manipulação de path
library(tidyr)
library(stringr) #para busca por padrão na string
library(AnnotationForge) #fazer o db
library(AnnotationDbi) #fazer o db
library(RSQLite) #fazer o db
library(devtools) #instalação de pacote local

# Cria a pasta "resultados" se não existir
if (!dir.exists(here("data/resultados"))) dir.create(here("data/resultados"))

setwd("~/UFMG/ArtigoSplicing/data")

# Carregar lista de genes (coluna única)
gene_symbols = read.csv('listaCryptococcusEntrezID.csv', sep = ';', dec = '.', header = TRUE)
genesCrypto = gene_symbols$GeneID

# Enriquecimento KEGG para Cryptococcus neoformans (cne)
kegg_resultCrypto <- enrichKEGG(
  gene          = genesCrypto,
  organism      = "cng",
  keyType       = "kegg",
  pvalueCutoff  = 0.25
  #pAdjustMethod = 0.05
)



write.csv(kegg_resultCrypto@result, here("data/resultados", "keggEnrichmentResultsCryptoPValue0.25csv"), row.names = FALSE)

head(kegg_resultCrypto)
kegg_resultCrypto

kegg_resultCrypto@result$Description


viasKeggSignificantes <- subset(kegg_resultCrypto@result, p.adjust < 0.25)[, c("ID", "Description", "Count")]
print(viasKeggSignificantes)

# Visualização KEGG (substitua pela via de interesse)
setwd(here("data/resultados"))

pathview(
  gene.data   = genesCrypto,
  pathway.id  = "cng04144",   
  species     = "cng",
  gene.idtype = "KEGG"
)

# -------- ENRIQUECIMENTO GO --------

#Voltando para o diretório anterior de trabalho
setwd("~/UFMG/ArtigoSplicing/data")

#lista de genes do EntrezID
geneEntrezID = gene_symbols$EntrezID


if (!dir.exists(here("data/org.Cneoformans_var_grubii_H99.eg.db"))) {
  #fazendo o orgDB
  # Extrair as colunas importantes para o gene_info
  org = read.csv("data/uniprotCrypto3.csv", sep = ";", dec = ".", header = TRUE)
  
  gene_info = org %>%
    transmute(
      GeneID = Gene_Names,
      Symbol = word(Gene_Names, 1),    # pega o primeiro nome do campo Gene names
      Description = Protein.names
    ) %>%
    dplyr::rename(GID = GeneID) %>%
    #dplyr::select(GID) %>%
    dplyr::filter(str_detect(GID, "^CNAG_"))
  gene_info = distinct(gene_info)
  
  colnames(gene_info)
  
  dfLong = org %>%
    tidyr::separate_rows(Gene_Ontology_IDs, sep = "; \\s*") %>%
    dplyr::distinct(Gene_Names, Gene_Ontology_IDs, .keep_all = TRUE) %>%
    dplyr::rename(GID = Gene_Names, GO = Gene_Ontology_IDs ) %>%
    dplyr::select(GID, GO) %>%
    dplyr::filter(str_detect(GID, "^CNAG_"))
  
  
  dfLong$EVIDENCE <- "IEA"
  
  # Reordenar para garantir a ordem certa
  dfLong <- dfLong[, c("GID", "GO", "EVIDENCE")]
  head(dfLong)
  
  dfLong = as.data.frame(dfLong)
  
  # Confirmar estrutura
  str(dfLong)
  
  # Remove espaços, forçando padrão GO:XXXXXXX
  dfLong$GO = gsub("\\s", "", dfLong$GO)
  #Remove linhas duplicadas
  duplicated_rows = dfLong[duplicated(dfLong), ]
  dfLong = dfLong[!duplicated(dfLong), ]
  table(nchar(dfLong$GO))
  
  #Retira os valores com < 10 caracteres no GO
  head (dfLong)
  dfLong <- dfLong[nchar(dfLong$GO) == 10 & grepl("^GO:\\d{7}$", dfLong$GO), ]
  
  
  # Passo 2: rodar a função correta
  makeOrgPackage(
    gene_info = gene_info,
    go = dfLong,
    version = "0.1",
    maintainer = "Grace Avelar <grace.biologia@gmail.com>",
    author = "Grace Avelar",
    outputDir = "data",
    tax_id = "235443",
    genus = "Cryptococcus",
    species = "neoformans.var.grubii.H99",
    goTable = "go"
  )
  
  
  if (!requireNamespace("org.Cneoformans.var.grubii.H99.eg.db", quietly = TRUE)) {
    devtools::install(normalizePath(here("data/org.Cneoformans.var.grubii.H99.eg.db")),
                     repos = NULL, type = "source")
  }
  
}

library(org.Cneoformans.var.grubii.H99.eg.db)

ego = enrichGO(
  gene  = geneList, #vetor com os genes de interesse
  OrgDb        = org.Cneoformans.var.grubii.H99.eg.db,
  keyType      = "GID",
  ont          = "CC", # "BP" = Biological Process (processos biológicos), 
                       # "MF" = Molecular Function (funções moleculares)
                       # "CC" = Cellular Component (componentes celulares)
  pAdjustMethod = "bonferroni",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  #pvalueCutoff = 0.05, # Limite de valor de p bruto para considerar significância (ex: 0.05).
  #qvalueCutoff = 0.2, #Limite de q-valor (FDR ajustado) para filtrar os resultados (ex: 0.2).
  readable = FALSE
  )

library(dplyr)


# Gráfico de barras
ego_cc <- ego@result%>%
  top_n(20)%>%
  arrange(pvalue) %>%
  dplyr::select(Description, Count, pvalue)

barPlot_CC <- ggplot(ego_cc, aes(x = reorder(Description, pvalue), 
                                 y = Count, fill = -log10(pvalue))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 20 Termos GO Enriquecidos", x = "Termos GO", 
       y = "Número de Genes", fill = "-log10(p-value)") +
  theme_minimal()

ggsave("resultados/barPlot_Crypto_GO_CC2.jpg", plot = barPlot_BP, width = 10, height = 6, dpi = 600)

# Cnetplot
cnetplot_CC = cnetplot(ego, showCategory = 10, categorySize = "pvalue")
ggsave("resultados/cnetplot_Crypto_GOBP.jpg", plot = cnetplot_BP, width = 8, height = 6, dpi = 600)

# Mapa de Enriquecimento
ego_GO <- pairwise_termsim(ego_GO)
emapplot(ego_GO, showCategory = 20)
