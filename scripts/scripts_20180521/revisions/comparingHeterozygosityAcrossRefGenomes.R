############ Plotting individual heterozygosity 
require(ggplot2)
require(dplyr)
require(RColorBrewer)
wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/hetPerIndWithDifferentRefGenomes/"
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/hetPerIndWithDifferentRefGenomes/plots/"
colorPal=RColorBrewer::brewer.pal(n=8,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[7],BC=colorPal[7],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors
refs=c("MFUR","SSO","NSO")
################ All sites, no missingness filter, no admixed/relatives #############
allResults <- data.frame()
for(ref in refs){
file=read.table(paste(wd,ref,"/all_9_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz.perIndHet.",ref,".txt",sep=""),header=T)
file$ref <- ref
allResults <- rbind(allResults,file)
}

require(ggplot2)
p1 <- ggplot(allResults,aes(x=ref,y=HeterozygosityPerSite))+
  geom_boxplot()
p1
means <- allResults %>%
  group_by(ref) %>%
  summarise(meanHet = mean(HeterozygosityPerSite))
means

allResults$pop <- unlist(lapply(strsplit(as.character(allResults$Sample),"_"),"[",3))
# get overall pop (ber and med combine into COM)
allResults$pop2 <- allResults$pop
allResults[allResults$pop=="BER"|allResults$pop=="MED",]$pop2 <- "COM"
allResults[allResults$pop=="BAJ",]$pop2 <- "BC"

p2 <- ggplot(allResults, aes(x=pop,y=HeterozygosityPerSite,fill=ref))+
  geom_boxplot()

p2

# want to plot with line between them 
allResults$ref  <- factor(allResults$ref,levels=c("MFUR","SSO","NSO"))
p3 <- ggplot(allResults,aes(x=ref,y=HeterozygosityPerSite,group=Sample,color=ref))+ 
  geom_path(color="gray")+
  geom_point(alpha=0.8)
  
p3

## want a scatter plot of each set of genomes
require(reshape2)
#?dcast
# http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/


############## All sites (no miss.) Plots for manuscript ############
allResults_HetPerRef <- dcast(allResults,pop2+Sample~ref,value.var="HeterozygosityPerSite") # woah where has THIS been my whole life. 
### MFUR vs SSO
p4a <- ggplot(allResults_HetPerRef,aes(x=MFUR,y=SSO,color=pop2))+
  geom_point(size=3,alpha=0.75)+
  geom_abline(slope=1)+
  scale_color_manual(values=unlist(colors))+
  theme_bw()+  
  xlim(c(7e-5,2.1e-4))+
  ylim(c(7e-5,2.1e-4))+
  ggtitle("All-sites per individual het.\n(no missingness filter)")+
  theme(legend.title=element_blank())
p4a
ggsave(paste(plot.dir,"all9.MFUR.vs.SSO.HetPerInd.nomissingfilter.pdf",sep=""),p4a,device="pdf",height=3,width=4)

## mfur vs nso
p4b <- ggplot(allResults_HetPerRef,aes(x=MFUR,y=NSO,color=pop2))+
  geom_point(size=3,alpha=0.75)+
  geom_abline(slope=1)+
  scale_color_manual(values=unlist(colors))+
  theme_bw()+  
  xlim(c(7e-5,2.1e-4))+
  ylim(c(7e-5,2.1e-4))+
  ggtitle("All-sites per individual het.\n(no missingness filter)")+
  theme(legend.title=element_blank())
p4b
ggsave(paste(plot.dir,"all9.MFUR.vs.NSO.HetPerInd.nomissingfilter.pdf",sep=""),p4b,device="pdf",height=3,width=4)
## nso vs sso
p4c <- ggplot(allResults_HetPerRef,aes(x=SSO,y=NSO,color=pop2))+
  geom_point(size=3,alpha=0.75)+
  geom_abline(slope=1)+
  scale_color_manual(values=unlist(colors))+
  theme_bw()+  
  xlim(c(7e-5,2.1e-4))+
  ylim(c(7e-5,2.1e-4))+
  ggtitle("All-sites per individual het.\n(no missingness filter)")+
  theme(legend.title=element_blank())
p4c
ggsave(paste(plot.dir,"all9.SSO.vs.NSO.HetPerInd.nomissingfilter.pdf",sep=""),p4c,device="pdf",height=3,width=4)

############# sandboxing #############
# want to find ones where it goes down -- only 5 out of 84.
samples <- unlist(unique(allResults$Sample))
samples
HetGetsBiggerInSSO <- list()
HetGetsLowerOrSameInSSO <- list()
for(sample in samples){
  MFUR_SSO_Diff <- allResults[allResults$Sample==sample & allResults$ref=="MFUR",]$HeterozygosityPerSite - allResults[allResults$Sample==sample & allResults$ref=="SSO",]$HeterozygosityPerSite
  # if MFUR - SSO is negative, that means SSO is bigger than mfur:
  if(MFUR_SSO_Diff < 0){
    HetGetsBiggerInSSO <- c(HetGetsBiggerInSSO,sample)
  } else if(MFUR_SSO_Diff >= 0){
    HetGetsLowerOrSameInSSO <- c(HetGetsLowerOrSameInSSO,sample)
  }
}

############## Neutral Region Het #############

wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/hetPerIndWithDifferentRefGenomes/"
refs=c("MFUR","SSO") # didn't do NSO
# concat files
neutralResults <- data.frame()
for(ref in refs){
  file=read.table(paste(wd,ref,"/neutral.all_9_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz.perIndHet.",ref,".txt",sep=""),header=T)
  file$ref <- ref
  neutralResults <- rbind(neutralResults,file)
}
neutralResults$pop <- unlist(lapply(strsplit(as.character(neutralResults$Sample),"_"),"[",3))
neutralResults$pop2 <- neutralResults$pop
neutralResults[neutralResults$pop=="BER"|neutralResults$pop=="MED",]$pop2 <- "COM"
neutralResults[neutralResults$pop=="BAJ",]$pop2 <- "BC"
require(dplyr)
require(ggplot2)
NeutralMeans <- neutralResults %>%
  group_by(ref) %>%
  summarise(meanHet = mean(HeterozygosityPerSite))
NeutralMeans

p1 <- ggplot(neutralResults,aes(x=ref,y=HeterozygosityPerSite))+
  geom_boxplot()
p1
## is this evidence for the mfur regions being MORE neutral? argues in favor of our dropouts I think.
neutralResults$ref  <- factor(neutralResults$ref,levels=c("MFUR","SSO"))
p3 <- ggplot(neutralResults,aes(x=ref,y=HeterozygosityPerSite,group=Sample,color=ref))+ 
  geom_path(color="gray")+
  geom_point(alpha=0.8)

p3

############## Neutral (no miss.) Plots for manuscript ############
neutralResults_HetPerRef <- dcast(neutralResults,pop2+Sample~ref,value.var="HeterozygosityPerSite") # woah where has THIS been my whole life. 
### MFUR vs SSO
p6a <- ggplot(neutralResults_HetPerRef,aes(x=MFUR,y=SSO,color=pop2))+
  geom_point(size=3,alpha=0.75)+
  geom_abline(slope=1)+
  scale_color_manual(values=unlist(colors))+
  theme_bw()+  
  xlim(c(7e-5,2.1e-4))+
  ylim(c(7e-5,2.1e-4))+
  ggtitle("Neutral per individual het.\n(no missingness filter)")+
  theme(legend.title=element_blank())
p6a
ggsave(paste(plot.dir,"neutral.MFUR.vs.SSO.HetPerInd.nomissingfilter.pdf",sep=""),p6a,device="pdf",height=3,width=4)

## don't have nso neutral regions.


############## all_10 -- this has admixed individuals included (and relatives) adn 2 pca outliers ; it also has 20% missingness filter #############
wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/hetPerIndWithDifferentRefGenomes/"
refs=c("MFUR","SSO") # didn't do NSO
# concat files
all10Results <- data.frame()
for(ref in refs){
  file=read.table(paste(wd,ref,"/all_10_maxHetFilter_0.75_ForRevPerIndHet_keepRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.2_rmBadIndividuals_passingFilters_raw_variants.vcf.gz.perIndHet.",ref,".txt",sep=""),header=T)
  file$ref <- ref
  all10Results <- rbind(all10Results,file)
}
all10Results$pop <- unlist(lapply(strsplit(as.character(all10Results$Sample),"_"),"[",3))
# need pop2 and label as admixed/relative/outlier:
all10Results$pop2 <- all10Results$pop
all10Results[all10Results$pop=="BER"|all10Results$pop=="MED",]$pop2 <- "COM"
all10Results[all10Results$pop=="BAJ",]$pop2 <- "BC"
admixedOutliersLabels <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/hetPerIndWithDifferentRefGenomes/ListOfRemovedIndividuals.Dups.Rels.Admixed.txt",header=T)

all10Results$label <- "Kept in SFS"
all10Results[all10Results$Sample %in% admixedOutliersLabels[admixedOutliersLabels$Label=="Outlier",]$Sample,]$label <- "PCA Outlier"
all10Results[all10Results$Sample %in% admixedOutliersLabels[admixedOutliersLabels$Label=="Admixed",]$Sample,]$label <- "Admixed"
all10Results[all10Results$Sample %in% admixedOutliersLabels[admixedOutliersLabels$Label=="Relative",]$Sample,]$label <- "Relative/Duplicate"
all10Results[all10Results$Sample %in% admixedOutliersLabels[admixedOutliersLabels$Label=="Duplicate",]$Sample,]$label <- "Relative/Duplicate"

require(dplyr)
require(ggplot2)
all10Means <- all10Results %>%
  group_by(ref) %>%
  summarise(meanHet = mean(HeterozygosityPerSite),meanCalledSites=mean(TotalCallCount))
all10Means

p1 <- ggplot(all10Results,aes(x=ref,y=HeterozygosityPerSite))+
  geom_boxplot()
p1
## is this evidence for the mfur regions being MORE neutral? argues in favor of our dropouts I think.
all10Results$ref  <- factor(all10Results$ref,levels=c("MFUR","SSO"))
p3 <- ggplot(all10Results,aes(x=ref,y=HeterozygosityPerSite,group=Sample,color=ref))+ 
  geom_path(color="gray")+
  geom_point(alpha=0.8)

p3
# scatter
all10Results_HetPerRef <- dcast(all10Results,pop2+Sample+label~ref,value.var="HeterozygosityPerSite") # woah where has THIS been my whole life. 
######################### all10 with admixed and missingness Plots for manuscript ##########
p7a <- ggplot(all10Results_HetPerRef,aes(x=MFUR,y=SSO,color=pop2))+
  geom_point(size=3,alpha=0.75)+
  geom_abline(slope=1)+
  scale_color_manual(values=unlist(colors))+
  theme_bw()+  
  xlim(c(7e-5,2.1e-4))+
  ylim(c(7e-5,2.1e-4))+
  ggtitle("All-sites per individual het.\n(20% missingness filter w/ admx/rels)")+
  theme(legend.title=element_blank())
  #scale_shape_manual(values=c(8,1,4,7))
p7a
ggsave(paste(plot.dir,"all10.MFUR.vs.SSO.HetPerInd.wmissingfilter.wadmixed.pdf",sep=""),p7a,device="pdf",height=3.2,width=4.2)


############ neutral all_10 ##################
wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/hetPerIndWithDifferentRefGenomes/"
refs=c("MFUR","SSO") # didn't do NSO
# concat files
neutral10Results <- data.frame()
for(ref in refs){
  file=read.table(paste(wd,ref,"/neutral.all_10_maxHetFilter_0.75_ForRevPerIndHet_keepRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.2_rmBadIndividuals_passingFilters_raw_variants.vcf.gz.perIndHet.",ref,".txt",sep=""),header=T)
  file$ref <- ref
  neutral10Results <- rbind(neutral10Results,file)
}
neutral10Results$pop <- unlist(lapply(strsplit(as.character(neutral10Results$Sample),"_"),"[",3))
neutral10Results$pop2 <- neutral10Results$pop
neutral10Results[neutral10Results$pop=="BER"|neutral10Results$pop=="MED",]$pop2 <- "COM"
neutral10Results[neutral10Results$pop=="BAJ",]$pop2 <- "BC"
require(dplyr)
require(ggplot2)
neutral10Means <- neutral10Results %>%
  group_by(ref) %>%
  summarise(meanHet = mean(HeterozygosityPerSite),meanCalledSites=mean(TotalCallCount))
neutral10Means

p1 <- ggplot(neutral10Results,aes(x=ref,y=HeterozygosityPerSite))+
  geom_boxplot()
p1

neutral10Results$ref  <- factor(neutral10Results$ref,levels=c("MFUR","SSO"))
p3 <- ggplot(neutral10Results,aes(x=ref,y=HeterozygosityPerSite,group=Sample,color=ref))+ 
  geom_path(color="gray")+
  geom_point(alpha=0.8)

p3


neutral10Results_HetPerRef <- dcast(neutral10Results,pop2+Sample~ref,value.var="HeterozygosityPerSite") # woah where has THIS been my whole life. 
### MFUR vs SSO
p8a <- ggplot(neutral10Results_HetPerRef,aes(x=MFUR,y=SSO,color=pop2))+
  geom_point(size=3,alpha=0.75)+
  geom_abline(slope=1)+
  scale_color_manual(values=unlist(colors))+
  theme_bw()+  
  xlim(c(7e-5,2.1e-4))+
  ylim(c(7e-5,2.1e-4))+
  ggtitle("Neutral per individual het.\n(20% missingness filter w/ admx/rels)")+
  theme(legend.title=element_blank())
p8a
ggsave(paste(plot.dir,"neutral10.MFUR.vs.SSO.HetPerInd.wmissingfilter.wadmixed.pdf",sep=""),p8a,device="pdf",height=3.2,width=4.2)

###################### Plotting admixed heterozygosity in neutral regions #########
#ggplot(all10Results_HetPerRef,aes(x=MFUR,y=SSO,color=pop2,shape=label))+
  geom_point()

neutral10Results_HetPerRef$label <- "Kept in SFS"
neutral10Results_HetPerRef[neutral10Results_HetPerRef$Sample %in% admixedOutliersLabels[admixedOutliersLabels$Label=="Outlier",]$Sample,]$label <- "PCA Outlier"
neutral10Results_HetPerRef[neutral10Results_HetPerRef$Sample %in% admixedOutliersLabels[admixedOutliersLabels$Label=="Admixed",]$Sample,]$label <- "Admixed"
neutral10Results_HetPerRef[neutral10Results_HetPerRef$Sample %in% admixedOutliersLabels[admixedOutliersLabels$Label=="Relative",]$Sample,]$label <- "Relative/Duplicate"
neutral10Results_HetPerRef[neutral10Results_HetPerRef$Sample %in% admixedOutliersLabels[admixedOutliersLabels$Label=="Duplicate",]$Sample,]$label <- "Relative/Duplicate"
head(neutral10Results_HetPerRef)
write.table(neutral10Results_HetPerRef,"/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/hetPerIndWithDifferentRefGenomes/HetPerIndividual.neutral10.includesAdmixedRelatives.SSO.MFUR.20percMissingnessFilter.NeutralRegions.txt",col.names=T,sep="\t",quote=F,row.names=F)
# get average per ind het among admixed? 
#### you are here -- need to figure out way to deal with this ####
neutral10Results_HetPerRef %>% 
  group_by(label) %>%
  summarise(meanMFUR=mean(MFUR),meanSSO=mean(SSO))
# A tibble: 4 x 3
#label              meanMFUR  meanSSO
#<chr>                 <dbl>    <dbl>
#1 Admixed            0.000143 0.000158
#2 Kept in SFS        0.000133 0.000146
#3 PCA Outlier        0.000126 0.000136
#4 Relative/Duplicate 0.000143 0.000154
#################### Plotting het in the 2.3Mb that were excluded from ferret #########
excludedFerretHet_fromall9 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/hetPerIndWithDifferentRefGenomes/MFUR/SitesThatWereExcludedDueToRepeatProximity/SitesThatFailedNeutralFilt.all_9_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz.perIndHet.MFUR.txt",header=T)

excludedFerretHet_fromall10 <-read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/hetPerIndWithDifferentRefGenomes/MFUR/SitesThatWereExcludedDueToRepeatProximity/SitesThatFailedNeutralFilt.all_10_maxHetFilter_0.75_ForRevPerIndHet_keepRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_0.2_rmBadIndividuals_passingFilters_raw_variants.vcf.gz.perIndHet.MFUR.txt",header=T) # has missingness filter
 ### want to compare to neutral het from 9 and 10
combine_includedExcludedMfur_fromall9 <- merge(excludedFerretHet_fromall9,neutralResults[neutralResults$ref=="MFUR",],by="Sample",suffixes = c(".excluded",".included"))
head(combine_includedExcludedMfur_fromall9)
p9a <- ggplot(combine_includedExcludedMfur_fromall9,aes(x=HeterozygosityPerSite.included,y=HeterozygosityPerSite.excluded,color=pop2))+
  geom_point(size=3,alpha=0.75)+
  geom_abline(slope=1)+
  scale_color_manual(values=unlist(colors))+
  theme_bw()+  
  xlim(c(7e-5,2.1e-4))+
  ylim(c(7e-5,2.1e-4))+
  ggtitle("Neutral per individual het.\n(no missingness filter w/out admx/rels)")+
  theme(legend.title=element_blank())
p9a
ggsave(paste(plot.dir,"included.vs.excluded.neutralSites9.MFUR.Only.HetPerInd.NOmissingfilter.NOadmixed.pdf",sep=""),p9a,device="pdf",height=3.2,width=4.2)

# get means:
mean(excludedFerretHet_fromall9$HeterozygosityPerSite) # 0.0001447941 in excluded regions 
mean(neutralResults[neutralResults$ref=="MFUR",]$HeterozygosityPerSite) #  0.0001641305 in included regions

# okay with missingness filter applied:

# okay so based on 
#### 
combine_includedExcludedMfur_fromall10 <- merge(excludedFerretHet_fromall10,neutral10Results[neutral10Results$ref=="MFUR",],by="Sample",suffixes = c(".excluded",".included"))
head(combine_includedExcludedMfur_fromall10)
p9b <- ggplot(combine_includedExcludedMfur_fromall10,aes(x=HeterozygosityPerSite.included,y=HeterozygosityPerSite.excluded,color=pop2))+
  geom_point(size=3,alpha=0.75)+
  geom_abline(slope=1)+
  scale_color_manual(values=unlist(colors))+
  theme_bw()+  
  xlim(c(7e-5,2.1e-4))+
  ylim(c(7e-5,2.1e-4))+
  ggtitle("Neutral per individual het.\n(20% missingness filter w admx/rels)")+
  theme(legend.title=element_blank())
p9b
ggsave(paste(plot.dir,"included.vs.excluded.neutralSites10.MFUR.Only.HetPerInd.wmissingfilter.wadmixed.pdf",sep=""),p9b,device="pdf",height=3.2,width=4.2)
mean(excludedFerretHet_fromall10$HeterozygosityPerSite) # 0.0001223757 in excluded regions 
mean(neutral10Results[neutral10Results$ref=="MFUR",]$HeterozygosityPerSite) #  0.0001340966 
################## What are summaries I want to put in paper ###########
