# Se voce nao tem o qiime2R roda esse codigo:
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

# Carrefa as librarias necessarias
pacman::p_load("qiime2R", "dplyr", "tibble", "phyloseq", "stringr")


metadata <- read.csv("metadata_table2019.csv") #Caminho do seu arquivo metadata
taxonomy_qza <- read_qza("pnat2019taxonomy_merge.qza") #QZA com a taxonomia
seqtab_qza <- read_qza("pnat2019table-arch-bact.qza") #QZA com as sequencias

# Esta funcao foi criada para colocar os nomes da terceira coluna da metadata
# nos nomes das amostras do QZA e criar um objeto phyloseq com os 3 arquivos fornecidos

qiime_to_physeq <- function(metadata, taxonomy, seqtab)  {# Carregamos a metadata
  metadata <- (metadata)
  
  # Renomeamos as linhas da metadata com o nome correto das amostras
  metadata <- (data.frame(metadata, row.names = metadata[,3])) #Muda row.names para
  # o numero que vc 
  
  # Renomeamos as linhas da metadata
  
  ## Carregar o arquivo da taxonomia
  
  taxonomy <- parse_taxonomy(taxonomy$data)
  
  ## Carregar o arquivo com as sequencias
  
  #Fazemos transposicao para renomear as lihas
  seqtab <- t(seqtab_qza$data)
  
  # Convertimos em dataframe para facilitar a manipulacao com o pacote "dplyr"
  seqtab <- as.data.frame(seqtab)
  
  # Adicionamos uma coluna com os nomes atuais das linhas 
  seqtab <- seqtab %>%
    mutate (X = row.names(seqtab))
  
  # Permite que adicionemos os nomes que de fato queremos que tenham nossas
  # amostras
  
  seqtab_merge <- merge(
    seqtab, #Nome do arquivo 1 com os ASVs
    metadata[c("X", "Codigo")], #Nome do arquivo metadata com a coluna 
    # referencia ("samples id x")e a coluna com os nomes corretos das amostras
    by = "X") # Nome das colunas de referencias para juntar  os dois dataframes,
  # se tiverem nomes diferentes ussar "by.x" e "by.y"
  
  # Renomear as linhas com os nomes das amostras
  rownames(seqtab_merge) <- seqtab_merge$Codigo
  
  # Excluir as colunas adicionadas "Codigo" e "X" para so ter a abudnacia de ASVs
  seqtab_merge <- seqtab_merge[, !colnames(seqtab_merge) %in% c("Codigo", "X")]
  
  # Convertemos em matrix para poder ussar o phyloseq
  taxonomy <- as.matrix(taxonomy)
  seqtab_merge <- as.matrix(seqtab_merge)
  
  ps <- phyloseq(otu_table(seqtab_merge, taxa_are_rows=FALSE), 
                  sample_data(metadata), 
                  tax_table(taxonomy))
  
  return(ps)}