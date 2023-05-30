# Plota a diversidade alpha de shannon

#MUDAR ANO PARA O NOME DA COLUNA NA METADATA QUE VC QUER PLOTAR (Ex periodo)
plot_richness(ps, x="Ano", 'Ano', measures=c("Observed", "Shannon"))

# estimate alpha-diversity: table/data.frame
ab_diversity = estimate_richness(physeq = ps, measures = c("Shannon", "Observed", "Simpson")) 
## add metadata columns to alpha-diversity data frame 
metadata <- sample_data(object = ps) %>% 
  data.frame(.) # convert to a data.frame
# check if all the rownames between metadata and alpha-diversity match to join them blindly
# add metadata to alpha-div
alpha_div_metadata <- cbind(ab_diversity, metadata)
# determine the median by groups: change the group 'SampleType' to the sample_data column 
#name group that you want to use
alpha_div_metadata_median <- alpha_div_metadata %>% 
  group_by(Ano) %>% 
  summarise_at(vars(Observed:Shannon), mean) %>% 
  as.data.frame(.)

# print results
View(alpha_div_metadata_median)
View(alpha_div_metadata)

#NMDS

ps <- readRDS("ps.rds")
ps.ord <- ordinate(ps, "NMDS", "bray")
p1 = plot_ordination(ps, ps.ord, type="samples", title="Anos", color = "Ano")
print(p1)


#PERMANOVA

# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(ps), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(phyloseq::otu_table(ps))
set.seed(100)
variaveis <- as(sample_data(ps), "data.frame")
variaveis <- variaveis %>%
  select(2) %>% 
  mutate (Ano = as.factor(Ano))

library(vegan)
library(permute)
library(lattice)
###OTUdf.log<- log10(OTUdf+1) #transforma para log
dist <- vegdist(OTUdf, method = "bray") #cria um  objeto dist
dist


ad.csv <- adonis2(dist ~ variaveis, method ="bray",perm=9999, data= variaveis) #so a invasao
ad.csv


# Calculate bray curtis distance matrix
ps_bray <- phyloseq::distance(ps, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps))

# Adonis test
a <- adonis2(ps_bray ~ Ano, data = sampledf, by = NULL)
a

#betadisp

dispersion <- betadisper(ps_bray, group=sampledf$Ano)
anova(dispersion)
permutest(dispersion)
permutest(dispersion, pairwise = TRUE)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse

