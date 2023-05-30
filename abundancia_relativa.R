library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble") 

#CONTANDO OS GENEROS, FAMILIAS, FILOS

# Essa funcao permite voce contar quantas familias/generos/especies unicos tem

contador_tax <- function (ps_object, tax_level){
# Create a factor corresponding to the Genera
genfac = factor(tax_table(ps_object)[, tax_level])
# Tabulate the counts for each genera in each sample
gentab = apply(tax_table(ps_object), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
return(gentab)
}

# Exemplo
taxa_count <- contador_tax (ps, "Family")
#Quantidade de reads para cada nivel taxonomico
write.csv(contador_tax, file='Family_tab.csv')


#Imprimir a abundancia relativa
abundancia_relativa.tudo <- (ps %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
                               arrange(OTU) %>% rename(ASV = OTU) %>% 
                               select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species, Sample, Abundance) %>%
                               spread(Sample, Abundance))

# Por nivel taxonomico
psr <- tax_glom(ps, "Phylum") #Reemplazar pelo nivel desejado
ps0 <- transform_sample_counts(psr, function(x) x / sum(x))
ps1 <- merge_samples(ps0)

write.table(ps1 %>% tax_glom(taxrank = "Phylum") %>% 
              transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>%
              select(Phylum, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "ps.relative_abundance.phylum.csv", sep = "\t", quote = F, row.names = F, col.names = T)

# Graficos de abundancia

#com filtro
plot_graphic <- function(ps_object, tax_level, condition){
psr <- tax_glom(ps_object, tax_level)
ps0 <- transform_sample_counts(psr, function(x) x / sum(x))
ps1 <- merge_samples(ps0, condition)
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill=tax_level) 
}

# Exemplo

plot_graphic(ps, "Phylum", "Ano")