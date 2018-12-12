#load R packages
library(gdsfmt)
library(SNPRelate)
require(ggplot2)
require(GGally)
require(ggrepel)
############### set up your colors -- keep this consistent across all plots ######
pops=c("CA","BAJ","AK","AL","COM","KUR")  # your populations

require(RColorBrewer)
display.brewer.pal("Dark2",n=8)
colorPal=RColorBrewer::brewer.pal(n=8,name = "Dark2")
# make Baja brown
colors=list(California=colorPal[1],Baja=colorPal[7],Alaska=colorPal[2],Aleutian=colorPal[3],Commander=colorPal[4],Kuril=colorPal[5]) # your population colors
##################################################################

calldate=20181119 # date gt's were called in format YYYYMMDD
todaysdate=format(Sys.Date(),format="%Y%m%d")
# file locations:

indir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/snps_gds/",calldate,"/QD8_containsLowCovInds/",sep="") # this is where your snp vcf file is and where you will save your gds file (downloaded from hoffman)
plotoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/PCA/",calldate,"/QD8_containsLowCovIndsAndRWAB/",sep="")
fileoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PCA/",calldate,"/",sep="")
dir.create(plotoutdir,recursive = T)
dir.create(fileoutdir,recursive = T)

#open the gds file
genofile <- snpgdsOpen(paste(indir,"snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.gds",sep=""))

# 20180802: adding LD snp pruning: (1min); r2 threshold : 0.2; recommended by SNPRelate tutorial
# https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#ld-based-snp-pruning
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F)
# I double checked, and it doesn't exclude any chromosomes due to any 22 chr cutoffs. cool
head(snpset)
# Get all selected snp id
snpset.id <- unlist(snpset)
head(snpset.id)
#pca (fast)
#sink(paste(fileoutdir,"/PCA.summary.",todaysdate,".txt",sep=""))


########### pop map ########
#population information
popmap = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/allSamples.Populations.PCA.txt",header=T) # this includes the RWAB samples
sample.id = as.character(popmap$Sample)
pop1_code = as.character(popmap$PrimaryPop)
pop2_code = as.character(popmap$SecondaryPop)
seq_code = as.character(popmap$sequencer)
note_code= as.character(popmap$notes) 
############################## First PCA: maf 0.06; all individuals (except low coverage which were removed before) ###########
# no maf filter!
samplesInGenofile=read.gdsn(index.gdsn(genofile, "sample.id"))
subCA <- samplesInGenofile[grep("CA",samplesInGenofile)]
subBAJ <- samplesInGenofile[grep("BAJ",samplesInGenofile)]
sub <- c(subCA, subBAJ)
# get rid of samples that are relatives, have high-nocall or are pca outliers:
sub2 <- sub[sub %in% popmap[popmap$notes=="good" ,]$Sample]
pca <- snpgdsPCA(genofile,snp.id=snpset.id,autosome.only = F, maf=0.06, missing.rate=0.2,sample.id = sub2)
#sink()
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc

tab1a <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2=factor(pop2_code)[match(pca$sample.id,sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],note=factor(note_code)[match(pca$sample.id,sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
head(tab1a)

#plot first 2 pc coloring by primary population
p0 <- ggplot(tab1a,aes(x=EV1,y=EV2,color=pop1,shape=note))+
  geom_point(size=3)+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  #scale_shape_manual(values=c(1,16,8))+
  scale_color_manual(values=unlist(colors))+
  geom_text_repel(aes(label=sample.id))
  
p0
