### Ying's code to run pca

# download packages
source("http://bioconductor.org/biocLite.R")
biocLite("gdsfmt")
biocLite("SNPRelate")

#load R packages
library(gdsfmt)
library(SNPRelate)

# working directory
rundate=20170914
setwd(paste("/Users/annabelbeichman/Documents/UCLA/Otters/CaptureDataProcessing_Analysis/CaptureResults/",sep=""))

########################## format vcf to gds #####################
# #read vcf, and reformat to gds
vcf.fn = paste("VCF/",rundate,"/allPops/elut.raw_variants.",rundate,".HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.vcf.gz",sep="")
gdsFile <- paste("PCA/",rundate,"/allPops/elut.raw_variants.",rundate,".HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.gds",sep="") # what you're going to name output gds
snpgdsVCF2GDS(vcf.fn, gdsFile, method="biallelic.only")
# VCF Format ==> SNP GDS Format
# Method: exacting biallelic SNPs
# Number of samples: 26
# Parsing "allPops/elut.raw_variants.20170914.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.vcf.gz" ...
# import 50881 variants.
# + genotype   { Bit2 26x50881, 323.0K } *
#   Optimize the access efficiency ...
# Clean up the fragments of GDS file:
#   open the file 'example.gds' (751.9K)
# # of fragments: 48
# save to 'example.gds.tmp'
# rename 'example.gds.tmp' (751.5K, reduced: 336B)
# # of fragments: 20

# get summary:
snpgdsSummary(paste("PCA/",rundate,"/allPops/elut.raw_variants.",rundate,".HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.gds",sep=""))
#
#The file name: /Users/annabelbeichman/Documents/UCLA/Otters/CaptureDataProcessing_Analysis/CaptureResults/VCF/20170914/allPops/elut.raw_variants.20170914.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.gds 
#The total number of samples: 26 
#The total number of SNPs: 50881 
#SNP genotypes are stored in SNP-major mode (Sample X SNP).

############# open gds file  ##############
genofile <- snpgdsOpen(gdsFile)
head(genofile)

############## pca #####################
pca <- snpgdsPCA(genofile,autosome.only = F,maf=0.06,missing.rate = 0.2) # where did these numbers come from?
# are these cutoffs?
# Excluding 20,571 SNPs (monomorphic: TRUE) 
# are these monomorphic alternate? I guess??
# Principal Component Analysis (PCA) on genotypes:
# Excluding 20,571 SNPs (monomorphic: TRUE, MAF: 0.06, missing rate: 0.2)
# Working space: 26 samples, 30,310 SNPs
# using 1 (CPU) core
# PCA:    the sum of all selected genotypes (0,1,2) = 1054622
# CPU capabilities: Double-Precision SSE2
# Tue Sep 19 16:31:07 2017    (internal increment: 18744)
# [==================================================] 100%, completed in 0s
# Tue Sep 19 16:31:07 2017    Begin (eigenvalues and eigenvectors)
# Tue Sep 19 16:31:07 2017    Done.

# variance proportion (%)
#variance proportion (%)
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc

# population information
popmap <- read.table(paste("PCA/",rundate,"/popMap.",rundate,".txt",sep=""),header=T)
sample.id = as.character(popmap$SAMPLE)
pop_code = as.character(popmap$POPULATION)                     
cap_code = as.character(popmap$CAPTURE)

#make a data.frame
tab <- data.frame(sample.id = pca$sample.id, pop = factor(pop_code)[match(pca$sample.id, sample.id)], cap=factor(cap_code)[match(pca$sample.id,sample.id)], EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
head(tab)

#plot first 2 pc
#plot(tab$EV1, tab$EV2, col=as.integer(tab$pop), pch=16,ylab=paste("PC2 (", pc[2],"%)"), xlab=paste("PC1 (", pc[1],"%)"))
#legend("top", legend=levels(tab$pop), pch=16, col=1:nlevels(tab$pop))

require(ggplot2)
ggplot(tab,aes(x=EV1,y=EV2,color=interaction(pop,cap)))+
  geom_point()+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))

# other PCs (confused)
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digit=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=interaction(tab$pop,tab$cap), labels=lbls)

#########choose a subset of samples, eg, remove some samples ##########
CA_only <- tab[tab$pop != "KUR",]$sample.id
pca <- snpgdsPCA(genofile, autosome.only = F, sample.id=CA_only, maf=0.06,missing.rate = 0.2)
tab <- data.frame(sample.id = pca$sample.id, pop = factor(pop_code)[match(pca$sample.id, sample.id)],cap=factor(cap_code)[match(pca$sample.id,sample.id)], EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc
ggplot(tab,aes(x=EV1,y=EV2,color=interaction(pop,cap)))+
  geom_point()+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))
# other PCs (confused)
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digit=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=interaction(tab$pop,tab$cap), labels=lbls)

