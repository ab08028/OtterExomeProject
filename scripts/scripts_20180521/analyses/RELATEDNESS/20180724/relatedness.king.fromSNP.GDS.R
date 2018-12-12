#load R packages
library(gdsfmt)
library(SNPRelate)
# guide: https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#principal-component-analysis-pca
# For relatedness analysis, identity-by-descent (IBD) estimation in SNPRelate can be done by either the method of moments (MoM) (Purcell et al., 2007) or maximum likelihood estimation (MLE) (Milligan, 2003; Choi et al., 2009). For both of these methods it is preffered to use a LD pruned SNP set.

calldate=20180724 # date gt's were called in format YYYYMMDD
todaysdate=format(Sys.Date(),format="%Y%m%d")
# file locations:

indir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/snps_gds/",calldate,"/",sep="") # this is where your snp vcf file is and where you will save your gds file (downloaded from hoffman)
plotoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/Relatedness/",calldate,sep="")
fileoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/Relatedness/",calldate,sep="")
dir.create(plotoutdir,recursive = T)
dir.create(fileoutdir,recursive = T)

#open the gds file
genofile <- snpgdsOpen(paste(indir,"/snp_5_passingAllFilters_postMerge_raw_variants.gds",sep=""))

# LD snp pruning: (1min); r2 threshold : 0.2; recommended by SNPRelate tutorial
# https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#ld-based-snp-pruning
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F)
head(snpset)
# Get all selected snp id
snpset.id <- unlist(snpset)
head(snpset.id)

#population information
popmap = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/allSamples.Populations.txt",header=T)
head(popmap)
sample.id = as.character(popmap$Sample)
pop1_code = as.character(popmap$PrimaryPop)
pop2_code = as.character(popmap$SecondaryPop)

# separate by population:
AK.id <- popmap[popmap$PrimaryPop == "Alaska",]$Sample
CA.id <- popmap[popmap$PrimaryPop == "California",]$Sample
KUR.id <- popmap[popmap$PrimaryPop == "Kuril",]$Sample
COM.id <- popmap[popmap$PrimaryPop == "Commander",]$Sample
AL.id <- popmap[popmap$PrimaryPop == "Aleutian",]$Sample


########### Relationship inference Using KING method of moments (all pops together) ########
# sink()
#sink(paste(fileoutdir,"/KING.relatedness.summary.",todaysdate,".txt",sep=""))
ibd.robust <- snpgdsIBDKING(genofile,autosome.only = F,remove.monosnp = T,missing.rate = 0.2, maf=0.05,snp.id=snpset.id) # removing monomorphing snps; adding ld pruning; this tightens it up a lot
#sink()
names(ibd.robust)
# close genofile
snpgdsClose(genofile)
# Excluding 12,655 SNPs (monomorphic: TRUE, MAF: 0.05, missing rate: 0.2)
# Working space: 118 samples, 5,697 SNPs
# Pairs of individuals
dat <- snpgdsIBDSelection(ibd.robust)
head(dat)

##################### Prepare data for plotting #################
# check that no animals from diff popluations have high kinship (would be sign of sample mixup)
# want to color by populations? or ones that are in same population? 
dat2 <- merge(dat,popmap,by.x="ID1",by.y="Sample")
dat3 <- merge(dat2,popmap,by.x = "ID2",by.y="Sample",suffixes = c("",".2"))
head(dat3)
dat3$match <- "different pops" # are ind 1 and ind 2 from same population, color green if so
dat3[dat3$PrimaryPop==dat3$PrimaryPop.2,]$match <- "same pop"
# make a color column
dat3$colorLabel <- "different pops"
dat3[dat3$match=="same pop",]$colorLabel <- as.character(dat3[dat3$match=="same pop",]$PrimaryPop)


##################### Plot with all comparisons (including across pops) ################

p3a <- ggplot(data=dat3,aes(x=IBS0, y=kinship))+
  geom_point(aes(color=colorLabel),alpha=0.3)+
  xlab("Proportion of Zero IBS")+
  ylab("Estimated Kinship Coefficient (KING-robust)")+
  ggtitle(paste("King Relatedness Results based on ",as.character(length(ibd.robust$snp.id))," LD Pruned SNPs, maf cutoff (0.05), no monomorphic",sep=""))+
  theme_bw()+
  theme(legend.title = element_blank())

p3a
ggsave(paste(plotoutdir,"/relatedness.King.coloredByPop",todaysdate,".pdf",sep=""),p3a,device="pdf",width = 8,height=5)
########################## Manichaikul boxes ################
# from the Manichaikul manuscript: cutoffs (see my Manichaikul presentation)
# values: 
a = 1/ (2^(3/2)) # 0.3535534
b = 1/ (2^(5/2)) # 0.1767767
c = 1/ (2^(7/2)) # 0.08838835 
d = 1/ (2^(9/2)) # 0.04419417
unrelatedBox=data.frame(ymin=0,ymax=d,xmin=(1-b),xmax=1)
thirdDegreeBox=data.frame(ymin=d,ymax=c,xmin=1-a,xmax=1-b)
secondDegreeBox=data.frame(ymin=c,ymax=b,xmin=0.365,xmax=a)
fullSibBox=data.frame(ymin=b,ymax=a,xmin=0.1,xmax=0.365)
ParentBox=data.frame(ymin=b,ymax=a,xmin=0,xmax=0.1)
MonozygoticBox=data.frame(ymin=a,ymax=Inf,xmin=0,xmax=0.1)
# these values are from manichaikul 2010 table 1 
# when using geom rect set aes for each thing separately


##################### Plot with only within pop comparisons ################

p4a <- ggplot(data=dat3[dat3$match=="same pop",],aes(x=IBS0, y=kinship))+
  geom_point(aes(color=colorLabel),alpha=0.8)+
  xlab("Proportion of Zero IBS")+
  ylab("Estimated Kinship Coefficient (KING-robust)")+
  ggtitle(paste("King Relatedness Results based on ",as.character(length(ibd.robust$snp.id))," LD Pruned SNPs, maf cutoff (0.05), no monomorphic",sep=""))+
  theme_bw()+
  theme(legend.title = element_blank())
p4a
ggsave(paste(plotoutdir,"/relatedness.King.coloredByPop.OnlyInPopComparisons.",todaysdate,".pdf",sep=""),p4a,device="pdf",width = 8,height=5)


##################### Plot with all comparisons (across all pops) with Manichaikul boxes #################

p3b <- ggplot()+
  geom_rect(data=unrelatedBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=unrelatedBox$xmin+0.1,y=unrelatedBox$ymin+0.02,label="unrelated"))+
  geom_rect(data=thirdDegreeBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=thirdDegreeBox$xmin+0.1,y=thirdDegreeBox$ymin+0.02,label="third degree"))+
  geom_rect(data=secondDegreeBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=secondDegreeBox$xmin+0.05,y=secondDegreeBox$ymin+0.02,label="second degree"))+
  geom_rect(data=fullSibBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=fullSibBox$xmin+0.1,y=fullSibBox$ymin+0.1,label="full sib"))+
  geom_rect(data=ParentBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3,fill="blue")+
  geom_text(aes(x=ParentBox$xmin+.1,y=ParentBox$ymin+.15,label="parent-offspring"))+
  geom_rect(data=MonozygoticBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3,fill="red")+
  geom_text(aes(x=MonozygoticBox$xmin + 0.05,y=MonozygoticBox$ymin + 0.15,label="monozygotic"))+
  geom_point(data=dat3,aes(x=IBS0, y=kinship,color=colorLabel),alpha=0.3)+
  xlab("Proportion of Zero IBS")+
  ylab("Estimated Kinship Coefficient (KING-robust)")+
  theme_bw()+
  ggtitle(paste("King Relatedness Results based on ",as.character(length(ibd.robust$snp.id))," LD Pruned SNPs, maf cutoff (0.05), no monomorphic",sep=""))+
  scale_x_continuous(limits = c(0,1))+
  theme(legend.title = element_blank())
p3b
ggsave(paste(plotoutdir,"/relatedness.King.Boxes.",todaysdate,".pdf",sep=""),p3b,device="pdf",width = 8,height=5)


##################### Plot with all just within pop comparisons #################
p4b <- ggplot()+
  geom_rect(data=unrelatedBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=unrelatedBox$xmin+0.1,y=unrelatedBox$ymin+0.02,label="unrelated"))+
  geom_rect(data=thirdDegreeBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=thirdDegreeBox$xmin+0.1,y=thirdDegreeBox$ymin+0.02,label="third degree"))+
  geom_rect(data=secondDegreeBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=secondDegreeBox$xmin+0.05,y=secondDegreeBox$ymin+0.02,label="second degree"))+
  geom_rect(data=fullSibBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=fullSibBox$xmin+0.1,y=fullSibBox$ymin+0.1,label="full sib"))+
  geom_rect(data=ParentBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3,fill="blue")+
  geom_text(aes(x=ParentBox$xmin+.1,y=ParentBox$ymin+.15,label="parent-offspring"))+
  geom_rect(data=MonozygoticBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3,fill="red")+
  geom_text(aes(x=MonozygoticBox$xmin + 0.05,y=MonozygoticBox$ymin + 0.15,label="monozygotic"))+
  geom_point(data=dat3[dat3$match=="same pop",],aes(x=IBS0, y=kinship,color=colorLabel),alpha=0.3)+
  xlab("Proportion of Zero IBS")+
  ylab("Estimated Kinship Coefficient (KING-robust)")+
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle(paste("King Relatedness Results based on ",as.character(length(ibd.robust$snp.id))," LD Pruned SNPs, maf cutoff (0.05), no monomorphic",sep=""))+
  scale_x_continuous(limits = c(0,1))
p4b
ggsave(paste(plotoutdir,"/relatedness.King.Boxes.onlyWithinPops.",todaysdate,".pdf",sep=""),p4b,device="pdf",width = 8,height=5)

# how does fixed differences affect this? if both are AA AA due to fixed diff, it might count?
# but I think monomorphic sites are excluded right? 

# looks like a few samples have almost zero prop of zero IBS (share a lot)
# and kinsihp of 0.5 (monozygotic twins or sample sample twice)
# which samples are those 
relatives <- dat[dat$kinship>0.4,]
# want to pick one of each one
relatives

############################# Plot average kinship coeff per population ##############
head(dat3)
require(dplyr)
meanKinship <- dat3[dat3$match=="same pop",] %>%
  group_by(PrimaryPop) %>%
  summarize(MeanKinship=mean(kinship))
p5 <- ggplot(data=meanKinship,aes(x=reorder(PrimaryPop,MeanKinship),y=MeanKinship,fill=PrimaryPop))+
  geom_bar(stat="identity")+
  theme_bw()+
  theme(legend.title = element_blank())+
  ylab("Mean Kinship (KING)")+
  xlab("Population")+
  ggtitle("Mean Kinship")
p5
ggsave(paste(plotoutdir,"/relatedness.meanKinship.King.",todaysdate,".pdf",sep=""),p5,device="pdf",width = 8,height=5)

