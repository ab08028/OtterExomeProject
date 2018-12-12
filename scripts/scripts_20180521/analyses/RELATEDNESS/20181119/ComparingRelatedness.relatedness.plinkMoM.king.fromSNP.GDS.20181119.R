#load R packages
library(gdsfmt)
library(SNPRelate)
require(ggplot2)
# guide: https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#principal-component-analysis-pca
# For relatedness analysis, identity-by-descent (IBD) estimation in SNPRelate can be done by either the method of moments (MoM) (Purcell et al., 2007) or maximum likelihood estimation (MLE) (Milligan, 2003; Choi et al., 2009). For both of these methods it is preffered to use a LD pruned SNP set.

calldate=20181119 # date gt's were called in format YYYYMMDD
todaysdate=format(Sys.Date(),format="%Y%m%d")
# file locations:
pops=c("CA","BAJ","AK","AL","COM","KUR")
indir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/snps_gds/",calldate,"/",sep="") # this is where your snp vcf file is and where you will save your gds file (downloaded from hoffman)
plotoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/Relatedness/",calldate,sep="")
fileoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/Relatedness/",calldate,sep="")
dir.create(plotoutdir,recursive = T)
dir.create(fileoutdir,recursive = T)

#open the gds file
genofile <- snpgdsOpen(paste(indir,"/snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.gds",sep=""))
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
# LD snp pruning: (1min); r2 threshold : 0.2; recommended by SNPRelate tutorial
# https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#ld-based-snp-pruning
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F)
head(snpset)
# Get all selected snp id
snpset.id <- unlist(snpset)
head(snpset.id)

#population information
popmap = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/allSamples.Populations.PCA.txt",header=T)
head(popmap)
sample.id = as.character(popmap$Sample)
pop1_code = as.character(popmap$PrimaryPop)
pop2_code = as.character(popmap$SecondaryPop)

# separate by population and make sure the samples are in dataset (removed some low cov ind)
AK.id <- popmap[popmap$PrimaryPop == "Alaska" & popmap$Sample %in% samp.id,]$Sample
CA.id <- popmap[popmap$PrimaryPop == "California" & popmap$Sample %in% samp.id,]$Sample
KUR.id <- popmap[popmap$PrimaryPop == "Kuril" & popmap$Sample %in% samp.id,]$Sample
COM.id <- popmap[popmap$PrimaryPop == "Commander" & popmap$Sample %in% samp.id,]$Sample
AL.id <- popmap[popmap$PrimaryPop == "Aleutian" & popmap$Sample %in% samp.id,]$Sample
BAJ.id <- popmap[popmap$PrimaryPop == "Baja" & popmap$Sample %in% samp.id,]$Sample
########## RElationship inference using Plink IBD #########
# each population separately
# Estimate IBD coefficients
# removing monosnps (20181211)
AK.ibd <- snpgdsIBDMoM(genofile, sample.id=AK.id, snp.id=snpset.id, remove.monosnp = T,
                    maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
CA.ibd <- snpgdsIBDMoM(genofile, sample.id=CA.id, snp.id=snpset.id,remove.monosnp = T,
                       maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
KUR.ibd <- snpgdsIBDMoM(genofile, sample.id=KUR.id, snp.id=snpset.id,remove.monosnp = T,
                       maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
COM.ibd <- snpgdsIBDMoM(genofile, sample.id=COM.id, snp.id=snpset.id,remove.monosnp = T,
                       maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
AL.ibd <- snpgdsIBDMoM(genofile, sample.id=AL.id, snp.id=snpset.id,remove.monosnp = T,
                       maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
BAJ.ibd <- snpgdsIBDMoM(genofile, sample.id=BAJ.id, snp.id=snpset.id,remove.monosnp = T,
                        maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
# baja + CA:
CA.BAJ.ibd <- snpgdsIBDMoM(genofile, sample.id=c(as.character(BAJ.id),as.character(CA.id)), snp.id=snpset.id,
                           maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
# baja + CA + AK:
CA.AK.BAJ.ibd <- snpgdsIBDMoM(genofile, sample.id=c(as.character(BAJ.id),as.character(CA.id),as.character(AK.id)), snp.id=snpset.id,
                              maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
# Make a data.frame of them all
AK.ibd.coeff <- snpgdsIBDSelection(AK.ibd)
AK.ibd.coeff$population <- "AK"
CA.ibd.coeff <- snpgdsIBDSelection(CA.ibd)
CA.ibd.coeff$population <- "CA"
KUR.ibd.coeff <- snpgdsIBDSelection(KUR.ibd)
KUR.ibd.coeff$population <- "KUR"
COM.ibd.coeff <- snpgdsIBDSelection(COM.ibd)
COM.ibd.coeff$population <- "COM"
AL.ibd.coeff <- snpgdsIBDSelection(AL.ibd)
AL.ibd.coeff$population <- "AL"
BAJ.ibd.coeff <- snpgdsIBDSelection(BAJ.ibd)
BAJ.ibd.coeff$population <- "BAJ"
# put into one dataframe
all.ibd.coeff <- rbind(CA.ibd.coeff,AK.ibd.coeff,AL.ibd.coeff,COM.ibd.coeff,KUR.ibd.coeff,BAJ.ibd.coeff)
# color column
all.ibd.coeff$KinshipGreaterThan0.2 <- "no"
all.ibd.coeff[all.ibd.coeff$kinship > 0.2,]$KinshipGreaterThan0.2 <- "yes"
all.ibd.coeff$pop1 <- sapply(strsplit(all.ibd.coeff$ID1,"_"),"[",3)
all.ibd.coeff$pop2 <- sapply(strsplit(all.ibd.coeff$ID2,"_"),"[",3)
# fix BER/MED to COM
all.ibd.coeff[all.ibd.coeff$pop1=="BER" | all.ibd.coeff$pop1=="MED",]$pop1 <- "COM"
all.ibd.coeff[all.ibd.coeff$pop2=="BER" | all.ibd.coeff$pop2=="MED",]$pop2 <- "COM"
all.ibd.coeff$pop1==all.ibd.coeff$pop2 # should all be true now (no cross population comparisons)
p0 <- ggplot(all.ibd.coeff,aes(x=k0,y=k1,color=KinshipGreaterThan0.2))+
  geom_point()+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits=c(0,1))+
  geom_abline(intercept = 1,slope = -1,color="red",linetype="dashed")+
  facet_wrap(~population)+
  theme_bw()+
  scale_color_manual(values=c("dodgerblue","orange"))
p0
ggsave(paste(plotoutdir,"/relatedness.PlinkMoM.allPops.",todaysdate,".pdf",sep=""),p0,device="pdf",width = 8,height=5)

write.table(all.ibd.coeff[all.ibd.coeff$KinshipGreaterThan0.2=="yes",],paste(fileoutdir,"/relatedness.plinkMoM.",todaysdate,".individualsKinshipGT0.2.txt",sep=""),quote=F,row.names = F)

################### plot all kinship coefficients in a heatmap ##########
pops=c("CA" , "BAJ", "AK" , "AL" ,  "KUR")

for(pop in pops){
  # note: pop1 and pop2 should be the same, but this is defensive coding in case they aren't
  info = all.ibd.coeff[all.ibd.coeff$pop1==pop & all.ibd.coeff$pop2==pop,]
  p0b <- ggplot(info,aes(x=ID1,y=ID2,fill=kinship))+
    geom_tile()+
    scale_fill_gradientn(colours=c("white","gray","pink"))+
    geom_text(aes(label = round(kinship,2)))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(paste(pop," kinship coefficients from Plink\ncontains relatives/admixed/outliers"))
  p0b
  ggsave(paste(plotoutdir,"/",pop,".kinshipHeatMap.PlinkMoM.",todaysdate,".pdf",sep=""),p0b,device="pdf",width = 10,height=8)
}
# do Com separately because it's so big: 
pops="COM"
for(pop in pops){
  # note: pop1 and pop2 should be the same, but this is defensive coding in case they aren't
  info = all.ibd.coeff[all.ibd.coeff$pop1==pop & all.ibd.coeff$pop2==pop,]
  p0b <- ggplot(info,aes(x=ID1,y=ID2,fill=kinship))+
    geom_tile()+
    scale_fill_gradientn(colours=c("white","gray","pink"))+
    geom_text(aes(label = round(kinship,2)))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(paste(pop," kinship coefficients from Plink\ncontains relatives/admixed/outliers"))
  p0b
  ggsave(paste(plotoutdir,"/",pop,".kinshipHeatMap.PlinkMoM.",todaysdate,".pdf",sep=""),p0b,device="pdf",width = 20,height=16) # much bigger
}

############# California and Baja together ###########
CA.BAJ.ibd.coeff <- snpgdsIBDSelection(CA.BAJ.ibd)
CA.BAJ.ibd.coeff$pop1 <- sapply(strsplit(CA.BAJ.ibd.coeff$ID1,"_"),"[",3)
CA.BAJ.ibd.coeff$pop2 <- sapply(strsplit(CA.BAJ.ibd.coeff$ID2,"_"),"[",3)
CA.BAJ.ibd.coeff$KinshipGreaterThan0.2 <- "no"
CA.BAJ.ibd.coeff[CA.BAJ.ibd.coeff$kinship > 0.2,]$KinshipGreaterThan0.2 <- "yes"

p1a <- ggplot(CA.BAJ.ibd.coeff,aes(x=k0,y=k1,color=interaction(pop1,pop2)))+
  geom_point()+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits=c(0,1))+
  geom_abline(intercept = 1,slope = -1,color="red",linetype="dashed")+
  theme_bw()
  #scale_color_manual(values=c("dodgerblue","orange"))
p1a
ggsave(paste(plotoutdir,"/relatedness.PlinkMoM.CA.Baja.",todaysdate,".pdf",sep=""),p1a,device="pdf",width = 8,height=5)


############## CA + Baja plot kinship coefficients: heat map ################
p1b <- ggplot(CA.BAJ.ibd.coeff,aes(x=ID1,y=ID2,fill=kinship))+
  geom_tile()+
  scale_fill_gradientn(colours=c("white","gray","pink"))+
  geom_text(aes(label = round(kinship,2)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("CA-BAJ Kinship coefficients from Plink")
p1b
ggsave(paste(plotoutdir,"/CA.BAJ.kinshipHeatmap.PlinkMoM.CA.Baja.",todaysdate,".pdf",sep=""),p1b,device="pdf",width = 10,height=8)
############# California and Baja and Alaska together ###########
# CA.AK.BAJ.ibd.coeff <- snpgdsIBDSelection(CA.AK.BAJ.ibd)
# CA.AK.BAJ.ibd.coeff$pop1 <- sapply(strsplit(CA.AK.BAJ.ibd.coeff$ID1,"_"),"[",3)
# CA.AK.BAJ.ibd.coeff$pop2 <- sapply(strsplit(CA.AK.BAJ.ibd.coeff$ID2,"_"),"[",3)
# CA.AK.BAJ.ibd.coeff$KinshipGreaterThan0.2 <- "no"
# CA.AK.BAJ.ibd.coeff[CA.AK.BAJ.ibd.coeff$kinship > 0.2,]$KinshipGreaterThan0.2 <- "yes"
# 
# p2 <- ggplot(CA.AK.BAJ.ibd.coeff,aes(x=k0,y=k1,color=interaction(pop1,pop2)))+
#   geom_point()+
#   scale_x_continuous(limits = c(0,1))+
#   scale_y_continuous(limits=c(0,1))+
#   geom_abline(intercept = 1,slope = -1,color="red",linetype="dashed")+
#   theme_bw()
#   #scale_color_manual(values=c("dodgerblue","orange"))
# p2
# ggsave(paste(plotoutdir,"/relatedness.PlinkMoM.CA.Baja.AK.",todaysdate,".pdf",sep=""),p1,device="pdf",width = 8,height=5)


################## less interpetable: KING (skipping) ###############################
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

