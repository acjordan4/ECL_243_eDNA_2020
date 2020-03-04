#ECL 243 project pipeline
library(dada2)
path <- "~/Desktop/ECL243/Samples"
#path <- "C:/Users/Awnna/Desktop/ECL243/Samples/"
list.files(path)

#read in the names of the fastq files and list into vector
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), '[',1)
sample.names


#inspect read quality of forward and reverse reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#filter and trim the reads. trimming down to 200bp
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

output<- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(200,200), maxN = 0, maxEE = c(2,5), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
head(output)

#learn error rates and visualize
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

#sample inference
#can pool samples together to infer sample inference or do psuedo pooling
dadaFs <- dada(filtFs, err = errF)
dadaRs <- dada(filtRs, err = errR)

#DADA-Class output from above two lines,
#see help("dada-class") for full informaiton on the output
dadaFs[[1]]

#merge paired reads
#merged sequences only output if forward and reverse reads overlap by at least 12 bases and are identical in the overlap
#if most of your reads do not merge, then change truncation requirements
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
head(mergers[[1]])
mergers.check <- mergers[[1]]
write.csv(mergers.check, "output/checks.merges.csv")

#construct sequence table
#this sequence table lists each sequence as a column name and then the sample number is the row number
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#look at the distribution of sequence lengths
seq.lengths <- table(nchar(getSequences(seqtab)))
write.csv(seq.lengths, "output/distribution_sequence_lengths.csv")

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

head(seqtab.nochim)


#track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(output, sapply(dadaFs, getN), sapply(dadaRs,getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
rownames(track) <- sample.names
head(track)
write.csv(track, "output/track_reads_pipeline.csv")

#Assign Taxonomy
#taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Awnna/Desktop/ECL243/Samples/silva_nr_v132_train_set.fa")
#taxa <- addSpecies(taxa, "C:/Users/Awnna/Desktop/ECL243/Samples/silva_species_assignment_v132.fa")
taxa.gg <- assignTaxonomy(seqtab.nochim, "C:/Users/Awnna/Desktop/ECL243/Samples/rdp_train_set_14.fa.gz")
taxa.gg <- addSpecies(taxa.gg, "C:/Users/Awnna/Desktop/ECL243/Samples/rdp_species_assignment_14.fa.gz")

#and inspect
taxa.print <- taxa
rownames(taxa.print)<- NULL
head(taxa.print)
taxa.print[1,1:7]
write.csv(taxa.print,"output/frogs.txt")
names<- read.csv("output/frogs.txt", header =TRUE)
colnames(names)
names$name <- paste(names$Kingdom, names$Phylum, names$Class, names$Order, names$Family, names$Genus, names$Species)
names$name

#writing fasta file to feed into MSA for alignment
library(tidyverse)
library(seqinr)

#old code where names are a sequence of numbers
#names <- seq(from = 0, to = length(seqtab.nochim))
#seqs<-colnames(seqtab.nochim)
#write.fasta(as.list(colnames(seqtab.nochim)), names = names, "seqtab_nochim.txt", open = "w", as.string = FALSE)
#write.fasta(as.list(rownames(taxa)), names = taxa.print, "taxa.txt", open = "w", as.string = FALSE)
write.fasta(as.list(rownames(taxa)), names = names$name, "taxa.txt", open = "w", as.string = FALSE)

#Align Sequences
library(msa)
seqs.fasta <- readDNAStringSet("taxa.txt")
start_time <- Sys.time()
seqs.aligned <- msa(seqs.fasta, method="ClustalOmega")
end_time <- Sys.time()
end_time - start_time
seqs.aligned
writeXStringSet(unmasked(seqs.aligned), file="seqs_aligned.fasta")

#for this metadata function, I removed the other info in the readme file. will add to the github repository
#metadata <- read.delim('~/ECL243/README_for_amphibian_skin_microbiome1.txt', sep = "\t", col.names=c("sampleIDs","species","life_stage","site","year","bd_test"), row.names = NULL)
metadata <- read.delim('code/README_for_amphibian_skin_microbiome1.txt', sep = "\t", col.names=c("sampleIDs","species","life_stage","site","year","bd_test"), row.names = NULL)
dim(metadata)
metadata$sampleIDs

#Phyloseq statistical analysis and data visualization
library(phyloseq)
library(Biostrings)
library(ggplot2)

theme_set(theme_bw())
samples.out <- rownames(seqtab.nochim)
frog <- metadata$sampleIDs
species <- metadata$species
lifestage <- metadata$life_stage
bd_test <- metadata$bd_test
samdf <- data.frame(Frog = frog, frogspecies = species, lifestage = lifestage, bd_test = bd_test)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna) #add phylo tree here merge_phyloseq(ps,dna,tree)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#Filter the taxa based on mito and chloro
#taxa.filter <- remove_taxa(c("Mitochondria", "Chloroplast"), ps)

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna) #add phylo tree here merge_phyloseq(ps,dna,tree)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

ps <- subset_taxa(ps, Order!="Chloroplast" |Family !="Mitochondria" | is.na(Class))

library("microbiome")
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:10]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
#ps.top10 <- aggregate_taxa(ps, level= 'Family', top = 10)
p <- plot_bar(ps.top20, x="species", fill="Family", facet_grid=~samdf$species)
#p <- plot_bar(ps, fill = "Family", facet_grid=~samdf$species)
p + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")
ggsave("p.pdf", height = 25 , width = 30)


pdf("samples_families.pdf")
plot_bar(ps, fill = "Phylum")
dev.off()

plot_richness(ps, x="bd_test", measures=c("Shannon", "Simpson"), color="Species")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Species", title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Species:bd_test", fill="Family")

#heatmap function for abundance of each organism by phylum from https://joey711.github.io/phyloseq/import-data.html
#plot_heatmap(ps, taxa.label="Phylum")
