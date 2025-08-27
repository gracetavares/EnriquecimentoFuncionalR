# Meu Projeto 🚀

Este projeto teve como objetivo realizar a análise de enriquecimento funcional em dados genômicos de Cryptococcus neoformans H99 e Mus musculus. Para a análise foram utilizadas as ontologias do Gene ontology e as redes metabólicas do KEGG.  

## Sumário
- [Descrição do projeto](#descricao-do-projeto)
- [Sobre os scripts](#sobre-os-scripts)
    - [EnriquecimentoFuncionalCrypto.R](#enriquecimentofuncionalcryptor)
        - [Dependências](#-dependências)
    - [EnriquecimentoFuncional.R](#enriquecimentofuncionalr)
        - [Dependências](#-dependências-1)
    - [clusterSemanticaGOrMATSCrypto.R](#clustersemanticagormatscryptor)
        - [Dependências](#-dependências-2)
- [Status do projeto](#status-do-projeto)

## Descrição do projeto
Os scripts foram construídos em linguagem R e têm como objetivo realizar a análise de enriquecimento funcional em dados genômicos de Cryptococcus neoformans H99 e Mus musculus. Para a análise foram utilizadas as ontologias do Gene ontology (GO) e as redes metabólicas do KEGG. 
Os scripts que possuem o termo *semantica* realizam a análise com base nos clusters de semântica das ontologias do Gene Ontology. 
Os scripts podem ser utilizados de forma individual.

## Sobre os scripts
### EnriquecimentoFuncionalCrypto.R
Esse _script_ realiza a análise de enriquecimento funcional contra as ontologias do GO e as redes metabólicas do KEGG utilizando um banco de dados de _Cryptococcus neoformans var. H99_ construído manualmente. 
#### 📦 Dependências
Esse _script_ utiliza a **versão 4.5.1** do R e os seguintes pacotes do R:
- clusterProfiler
- enrichplot
- ggplot2
- GOSemSim
- dplyr
- pathview
- here 
- tidyr
- stringr
- AnnotationForge
- AnnotationDbi
- RSQLite
- devtools

_Para instalar:_
Todos os pacotes estão listados para instalação automática, caso necessário, no início do _script_. Selecione as linhas de instalação e execute-as.

### EnriquecimentoFuncional.R
Esse _script_ realiza a análise de enriquecimento funcional contra as ontologias do GO e as redes metabólicas do KEGG utilizando um banco de dados de _Mus musculus_ obtido da bioblioteca org.Mm.eg.db. 
#### 📦 Dependências
Esse _script_ utiliza a **versão 4.5.1** do R e os seguintes pacotes do R:
- clusterProfiler
- org.Mm.eg.db
- enrichplot
- ggplot2
- GOSemSim
- dplyr

_Para instalar:_
Todos os pacotes estão listados para instalação automática, caso necessário, no início do _script_. Selecione as linhas de instalação e execute-as.

### clusterSemanticaGOrMATSCrypto.R
Esse _script_ realiza a análise de enriquecimento funcional do resultado provenente do rMATS contra as ontologias do GO  utilizando um banco de dados de _Cryptococcus neoformans var. H99_ construído manualmente. A análise é feita baseando-se em _clusters_ semânticos das ontologias do Gene Ontology identificados e ao final gera um heatmap com a nuvem de palavras das ontologias mais representativas em cada _cluster_.

#### 📦 Dependências
Esse _script_ utiliza a **versão 4.5.1** do R e os seguintes pacotes do R:
- dplyr
- tidyr
- readr
- ComplexHeatmap
- clusterProfiler
- simplifyEnrichment
- circlize
- tibble
- here

_Para instalar:_
Todos os pacotes estão listados para instalação automática, caso necessário, no início do _script_. Selecione as linhas de instalação e execute-as.

## Status do projeto
![Badge em Desenvolvimento](http://img.shields.io/static/v1?label=STATUS&message=EM%20DESENVOLVIMENTO&color=GREEN&style=for-the-badge)

