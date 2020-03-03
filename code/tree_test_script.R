
###Ape Example

library(ape)
library(phangorn)
library(seqinr)

frogs_phyDat<- msaConvert(seqs.alligned, type= "phangorn::phyDat")
frogs_DnaBin <- msaConvert(seqs.alligned, type= "ape::DNAbin")

frogs10 <- subset(frogs_phyDat, 1:10)
frogs10 <- subset(frogs_DnaBin, 1:10)



## Other Method
dna_dist <- dist.dna(frogs_DnaBin, model="TN93")
tre <- njs(dna_dist)
class(tre)

##    Neighbor Joining, UPGMA, Maximum Parsimony                        
mt <- modelTest(frogs10)
print(mt)

                           
dna_dist <- dist.dna(frogs_DnaBin, model="TN93")

frogs_UPGMA <- upgma(dna_dist)
frogs_NJ <- NJ(dna_dist)
plot(frogs_UPGMA, main = "UPGMA")
plot(frogs_NJ, main = "Neighbor Joining")

parsimony(frogs_UPGMA, frogs10)
parsimony(frogs_NJ, frogs10)
frog_optim <- optim.parsimony(frogs_NJ, frogs_phyDat)
frogs_pratchet <- pratchet(frogs_phyDat)
plot(frog_optim)
plot(frogs_pratchet)

#Maximum LIkelihood nad Bootstrapping
fit <- pml(frogs_NJ, frogs_phyDat)
print(fit)
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=T, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")

write.tree(bs, file="bootstrap_samples.tsv")

                           