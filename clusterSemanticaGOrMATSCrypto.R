# Instalar as bibliotecas necessárias (lista de bibliotecas do CRAN)
pacotes = c ("BiocManager", "dplyr", "readr", "here", "circlize", "tibble" )
instalar = pacotes[!(pacotes %in% installed.packages()[,"Package"])]
if (any(!(pacotes %in% installed.packages()[,"Package"]))) {
  install.packages(instalar)
}

# Instalar as bibliotecas necessárias (Lista de bibliotecas do Bioconductor)
biocPacotes = c("ComplexHeatmap", "clusterProfiler", "simplifyEnrichment", "org.Mm.eg.db",
                "GenomicFeatures", "DESeq2", "clusterProfiler")
instalarBioc = biocPacotes[!(biocPacotes %in% installed.packages()[,"Package"])]
if (any(!(biocPacotes %in% installed.packages()[,"Package"]))) {
  BiocManager::install(instalarBioc)
}

#Carregar as bibibliotecas
library(dplyr)
library(tidyr)
library(readr)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db)   # trocar pelo organismo correto
library(simplifyEnrichment)
library(circlize)
library(tibble)
library(here)


#Diretório de trabalho
setwd("~/UFMG/ArtigoSplicing/EnriquecimentoFuncionalR/data")

# Criar a pasta "resultados" se não existir
if (!dir.exists(here("data/resultados"))) dir.create(here("data/resultados"))
setwd ("D:/UFMG/splicing/rmats_results_Tamires/Cneoformans_h99/ypdXbal10d/results2filtered")

# Carregar lista de genes (coluna única)
rmatsTableRI = read_tsv('RI.MATS.JC.txt') %>%
  mutate(Event = "RI")
rmatsTableA3SS = read_tsv('A3SS.MATS.JC.txt') %>%
  mutate(Event = "A3SS")
rmatsTableA5SS = read_tsv('A5SS.MATS.JC.txt') %>%
  mutate(Event = "A5SS")
rmatsTableMXE = read_tsv('MXE.MATS.JC.txt') %>%
  mutate(Event = "MXE")
rmatsTableSE = read_tsv('SE.MATS.JC.txt') %>%
  mutate(Event = "SE")


# Juntar as duas tabelas
rmatsTableGORI = rmatsTableRI %>%
  dplyr::select(GeneID,Event, FDR, IncLevelDifference) # acrescenta entrezID

rmatsTableGOA3SS = rmatsTableA3SS %>%
  dplyr::select(GeneID, Event, FDR, IncLevelDifference) 

rmatsTableGOA5SS = rmatsTableA5SS %>%
  dplyr::select(GeneID, Event, FDR, IncLevelDifference)

rmatsTableGOMXE = rmatsTableMXE %>%
  dplyr::select(GeneID, Event, FDR, IncLevelDifference) %>%
  mutate(across(c(FDR, IncLevelDifference), as.numeric))  

rmatsTableGOSE = rmatsTableSE %>%
  dplyr::select(GeneID, Event, FDR, IncLevelDifference) %>%
  mutate(across(c(FDR, IncLevelDifference), as.numeric))
  

rmatsTableAll = bind_rows(rmatsTableGORI, rmatsTableGOA3SS,
                          rmatsTableGOA5SS, rmatsTableGOMXE, 
                          rmatsTableGOSE)

# filtrar eventos significativos
EventosSignif = rmatsTableAll %>%
  filter(FDR <= 0.05, abs(IncLevelDifference) >= 0.1) %>%
  pull(GeneID) 

heatmapMatrix2 = rmatsTableAll %>%
  as.data.frame() %>%
  filter(GeneID %in% EventosSignif)


############Análise de enriquecimento funcional####################

#lista de genes do EntrezID
geneEntrezID = heatmapMatrix2$EntrezID
geneList = heatmapMatrix2$GeneID

if (!dir.exists(here("data/org.Cneoformans.var.grubii.H99.eg.db"))) {
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
  pAdjustMethod = "BH",
  #pvalueCutoff = 1,
  #qvalueCutoff = 1,
  pvalueCutoff = 0.05, # Limite de valor de p bruto para considerar significância (ex: 0.05).
  qvalueCutoff = 0.1, #Limite de q-valor (FDR ajustado) para filtrar os resultados (ex: 0.2).
  readable = FALSE #TODO: Verificar pq o TRUE não funciona!
)

#########Cluster semântico ######################


goClusters = simplifyGO(ego@result$ID,
                        plot = TRUE)#,
                        column_title = ("7 GO terms clustered by 'binary_cut'- Molecular Function"))

