#load R packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")
library(gdsfmt)
library(SNPRelate)
require(ggplot2)
require(ggally)
require(ggrepel)

########### set up colors ############
pops=c("CA","AK","AL","COM","KUR")  # your populations
require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=length(pops),name = "Dark2")
colors=list(California=colorPal[1],Alaska=colorPal[2],Aleutian=colorPal[3],Commander=colorPal[4],Kuril=colorPal[5]) 


calldate=20180806 # date gt's were called in format YYYYMMDD


todaysdate=format(Sys.Date(),format="%Y%m%d")
# file locations:

indir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/snps_gds/",calldate,"/",sep="") # this is where your snp vcf file is and where you will save your gds file (downloaded from hoffman)
plotoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/PCA/",calldate,"/",sep="")
fileoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PCA/",calldate,"/",sep="")
dir.create(plotoutdir,recursive = T)
dir.create(fileoutdir,recursive = T)

#open the gds file
# before removing admixed/relatives:
#genofile <- snpgdsOpen(paste(indir,"snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.gds",sep=""))
# after removing admixed/relatives:
genofile <- snpgdsOpen(paste(indir,"downsampled.COM.snp_8_rmRelativesAdmixed_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.gds",sep=""))

# 20180802: adding LD snp pruning: (1min); r2 threshold : 0.2; recommended by SNPRelate tutorial
# https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#ld-based-snp-pruning
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F)
# I double checked, and it doesn't exclude any chromosomes due to any 22 chr cutoffs. cool
head(snpset)
# Get all selected snp id
snpset.id <- unlist(snpset)
head(snpset.id)
#pca (fast)
sink(paste(fileoutdir,"/PCA.summary.postAdmixedRelativeRemoval",todaysdate,".txt",sep=""))
pca <- snpgdsPCA(genofile,snp.id=snpset.id,autosome.only = F, maf=0.06, missing.rate=0.2)
sink()
#variance proportion (%)
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc

#population information (doesn't matter if file is downsampled or not; has the record for every sequenced individual in it)
popmap = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/allSamples.Populations.txt",header=T) # this includes the RWAB samples
head(popmap)

sample.id = as.character(popmap$Sample)
pop1_code = as.character(popmap$PrimaryPop)
pop2_code = as.character(popmap$SecondaryPop)
seq_code = as.character(popmap$sequencer)
# can add as many of these as you want
#make a data.frame
tab <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2=factor(pop2_code)[match(pca$sample.id,sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
head(tab)

#plot first 2 pc coloring by primary population
require(ggplot2)
p1 <- ggplot(tab,aes(x=EV1,y=EV2,color=pop1,shape=sequencer))+
  geom_point()+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_color_manual(values=unlist(colors))+
  scale_shape_manual(values=c(1,16))
p1
ggsave(paste(plotoutdir,"/PCA.inclCA.LDPruned.postRelativeRemoval.AdmixedRemoval.downsampledCOM.",todaysdate,".pdf",sep=""),p1,device="pdf",width = 8,height=5)



#plot pc pairs for the first four pc 
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digit=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop1,labels=lbls)

#choose a subset of samples, eg exclude California & AK
sub <- tab$sample.id[tab$pop1!="California"]
sink(paste(fileoutdir,"/PCA.record.excl.CA.downsampledCOM.",todaysdate,".txt",sep=""))
pca <- snpgdsPCA(genofile, snp.id=snpset.id,autosome.only = F, sample.id=sub, maf=0.06)
sink()
# check output! 
tab2 <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2 = factor(pop2_code)[match(pca$sample.id, sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)

#plot first 2 pc coloring by primary population, excluding California
p2 <- ggplot(tab2,aes(x=EV1,y=EV2,color=pop1,shape=sequencer,label=sample.id))+
  geom_point()+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  scale_color_manual(values=unlist(colors))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,16))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))
p2
ggsave(paste(plotoutdir,"/PCA.excludeCA.LDPruned.postRelativeRemoval.AdmixedRemoval.downsampledCOM.",todaysdate,".pdf",sep=""),p2,device="pdf",width = 8,height=5)


#plot first 2 pc coloring by primary population, excluding California & Alaska
p3 <- ggplot(tab2,aes(x=EV1,y=EV2,color=pop1,shape=sequencer,label=sample.id))+
  geom_point()+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,16))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))+
  geom_text_repel(aes(label=sample.id))
p3
ggsave(paste(plotoutdir,"/PCA.excludeCA.LDPruned.postRelativeRemoval.AdmixedRemoval.withNames.downsampledCOM.",todaysdate,".pdf",sep=""),p3,device="pdf",width = 16,height=10)

#close gds file
snpgdsClose(genofile)
