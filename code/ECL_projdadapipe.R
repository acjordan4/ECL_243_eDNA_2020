#ECL 243 project pipeline
library(dada2)
path <- "C:/Users/MariyK/Desktop/ECL243/Samples"
list.files(path)

#read in the names of the fastq files and list into vector
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), '[',1)
sample.names

#inspect read quality of forward reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#filter and trim the reads
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

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#look at the distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

head(seqtab.nochim)

#Trying to write the sequences without chimeras to a text file to import into mafft or fast tree2 to
#isualize the phylogenies without prior taxonomy assignments
#so this doesn't work,, there's a couple reasons but mainly the format of the seqtab.nochim data isn't suitable/ agreeable to
#the required form for the write.fasta function

library(seqinr)
names <- seq(from = 0, to = length(seqtab.nochim))
sequences <-write.table(seqtab.nochim, file = "seqtab_nochim1.txt", sep = '\t')
write.fasta(sequences = sequences, names = names, "seqtab_nochim.txt", open = "w", as.string = FALSE)


       
#track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(output, sapply(dadaFs, getN), sapply(dadaRs,getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
ronames(track) <- sample.names
head(track)

#Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/MariyK/Desktop/ECL243/Samples/silva_nr_v132_train_set.fa.gz")
taxa <- addSpecies(taxa, "C:/Users/MariyK/Desktop/ECL243/Samples/silva_species_assignment_v132.fa.gz")
#and inspect
taxa.print <- taxa
rownames(taxa.print)
head(taxa.print)

#can check the accuracy of the assignment with additional code
#code not listed here butcan be found on dada2 documentation
#would also need the sample data to run this accuracy check


#Phyloseq statistical analysis and data visualization
library(phyloseq)
library(Biostrings)
library(ggplot2)

#to make the phyloseq table, need to look at example data in the dada2 tutorial
theme_set(theme_bw())
samples.out <- rownames(seqtab.nochim)
frog <- sapply(strsplit(samples.out, "D"), '[' , 1)

#This one doesn't work right now but should be inline with the rest of the code above
bact <- substr(subject, 1,1)



