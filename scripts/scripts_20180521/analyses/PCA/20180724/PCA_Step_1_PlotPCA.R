#load R packages
library(gdsfmt)
library(SNPRelate)


calldate=20180724 # date gt's were called in format YYYYMMDD
todaysdate=format(Sys.Date(),format="%Y%m%d")
# file locations:

indir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/snps_gds/",calldate,"/",sep="") # this is where your snp vcf file is and where you will save your gds file (downloaded from hoffman)
plotoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/PCA/",calldate,sep="")
fileoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PCA/",calldate,sep="")
dir.create(plotoutdir,recursive = T)
dir.create(fileoutdir,recursive = T)

#open the gds file
genofile <- snpgdsOpen(paste(indir,"/snp_5_passingAllFilters_postMerge_raw_variants.gds",sep=""))

# 20180802: adding LD snp pruning: (1min); r2 threshold : 0.2; recommended by SNPRelate tutorial
# https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#ld-based-snp-pruning
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F)
head(snpset)
# Get all selected snp id
snpset.id <- unlist(snpset)
head(snpset.id)
#pca (fast)
sink(paste(fileoutdir,"/PCA.summary.",todaysdate,".txt",sep=""))
pca <- snpgdsPCA(genofile,snp.id=snpset.id,autosome.only = F, maf=0.06, missing.rate=0.2)
sink()
#variance proportion (%)
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc

#population information
popmap = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/allSamples.Populations.txt",header=T)
head(popmap)
sample.id = as.character(popmap$Sample)
pop1_code = as.character(popmap$PrimaryPop)
pop2_code = as.character(popmap$SecondaryPop)

# can add as many of these as you want
#make a data.frame
tab <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2=factor(pop2_code)[match(pca$sample.id,sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
head(tab)

#plot first 2 pc coloring by primary population
require(ggplot2)
p1 <- ggplot(tab,aes(x=EV1,y=EV2,color=pop1))+
  geom_point()+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))
p1
ggsave(paste(plotoutdir,"/PCA.inclCA.LDPruned.",todaysdate,".pdf",sep=""),p1,device="pdf",width = 8,height=5)



#plot pc pairs for the first four pc 
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digit=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop1, labels=lbls)

#choose a subset of samples, eg exclude California
sub <- tab$sample.id[tab$pop1!="California"]
sink(paste(fileoutdir,"/PCA.record.excl.CA.",todaysdate,".txt",sep=""))
pca <- snpgdsPCA(genofile, snp.id=snpset.id,autosome.only = F, sample.id=sub, maf=0.06)
sink()
# check output! 
tab2 <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2 = factor(pop2_code)[match(pca$sample.id, sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)

#plot first 2 pc coloring by primary population, excluding California
p2 <- ggplot(tab2,aes(x=EV1,y=EV2,color=pop1))+
  geom_point()+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))
p2
ggsave(paste(plotoutdir,"/PCA.excludeCA.LDPruned.",todaysdate,".pdf",sep=""),p2,device="pdf",width = 8,height=5)

#close gds file
snpgdsClose(genofile)
