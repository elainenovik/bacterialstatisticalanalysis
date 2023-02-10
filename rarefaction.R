library(mirlyn)
library(readxl)
library(biomformat)
library(phyloseq)
library(Biostrings)

otu <- read.table("otu_table.tsv", row.names = 1, header = TRUE, sep = "\t")

otu.table <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
View(otu.table)
tax_table(otu.table) <- otu
rep.seqs <- Biostrings::readDNAStringSet("dna-sequences.fasta", format = "fasta")
sampledata<- read.table("sampledata.txt")
sample.data<-phyloseq::sample_names(sampledata)

expt <- phyloseq::phyloseq(otu.table, rep.seqs, sample.data)
expt


rarefy2<- rarefy_even_depth(expt)
View(rarefy2)

library(vegan)
rarecurve(otu_table(expt), step = 50, cex = 0.5)


plot_bar(expt, fill = "F")
View(expt)
rarefy1 <- alphacone(expt, rep = 1000)
exptmirl<- mirl(expt)
div <- alphadivDF(exptmirl, diversity = "shannon")

ps = import_biom("new_otu_table2.biom")

library(qiime2R)

phyloseqob <- qza_to_phyloseq(features = "table.qza", taxonomy = "silva-taxonomy.qza", metadata = "metadata-sample.tsv", tmp = "temp")
View(phyloseqob)

#Remove BANMAX
newphyloseqob = subset_samples(phyloseqob, sample_names(phyloseqob) != "ElaineNovik-BAN-MAX-16S")
View(newphyloseqob)

replacerare<-mirl(newphyloseqob, rep=1000, replace = TRUE)
View(replacerare)

alphadiv1<-alphacone(newphyloseqob, rep = 1000,
                     steps = seq(from = 0.001, to = 1, by = 0.01),
                     diversity = "shannon", replace = "FALSE") 
alphadiv2<-alphacone(newphyloseqob, rep = 100,
                     steps = seq(from = 0.001, to = 1, by = 0.01),
                     diversity = "shannon", replace = "FALSE")


alphadiv3<- alphadivDF(replacerare)
betadiv<-betamatPCA(replacerare)
View(betadiv)

mds <- metaMDS(betadiv)

names(alphadiv3)[names(alphadiv3) == "location"] <- "Alga_Sample"
View(alphadiv3)
hist(alphadiv3$DiversityIndex)
divplot<-ggplot(alphadiv3, aes(Alga_Sample, DiversityIndex))+ geom_jitter(aes(colour = Alga_Sample), width = 0.25, height = 0.25)
divplot

divbox<-ggplot(alphadiv3, aes(x = location, y = DiversityIndex)) + 
  geom_boxplot()
divbox


anova_test<-aov(DiversityIndex ~ Alga_Sample, data = alphadiv3)
summary(anova_test)
plot(anova_test)

resid<- residuals(object=anova_test)
shapiro.test(resid)

tukey<-TukeyHSD(anova_test)
plot(tukey)
