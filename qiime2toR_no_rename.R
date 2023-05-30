# Se voce nao tem o qiime2R roda esse codigo:
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

# Carrefa as librarias necessarias
pacman::p_load("qiime2R", "dplyr", "tibble", "phyloseq", "stringr")

# Esta funcao foi criada para colocar os nomes da terceira coluna da metadata
# nos nomes das amostras do QZA e criar um objeto phyloseq com os 3 arquivos fornecidos

qiime_to_physeq <- function(metadata, taxonomy, seqtab)  {# Carregamos a metadata
  metadata <- (metadata)
  taxonomy <- parse_taxonomy(taxonomy$data)
  seqtab <- t(seqtab$data)
    seqtab <- as.data.frame(seqtab)
  
  taxonomy <- as.matrix(taxonomy)
  seqtab_merge <- as.matrix(seqtab)
  
  ps <- phyloseq(otu_table(seqtab_merge, taxa_are_rows=FALSE), 
                  sample_data(metadata), 
                  tax_table(taxonomy))
  
  return(ps)}

# Exemplos 
metadata_rfi
asv_table <- read_qza("~/melhoramento_animal/priscilla/Asv_table.qza")
taxonomy <- read_qza("~/melhoramento_animal/priscilla/29_05_23_OUTPUT.qza")
  
ps <- qiime_to_physeq(metadata =metadata_rfi, taxonomy = taxonomy, seqtab = asv_table)

