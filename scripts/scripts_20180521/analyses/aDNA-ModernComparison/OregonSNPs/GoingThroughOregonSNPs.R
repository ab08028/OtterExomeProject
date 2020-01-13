require(ggplot2)
require(clipr)
require(reshape2)
############### Got oregon pseudohaploids ##############
CAfreqs <- read.table('/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/OregonCompareToModernAlleleFreqs/modern_alleleFreqsInCA_AK_atOregonSites/CA.AlleleFreqs.frq',header=F,skip=1,sep="\t")
colnames(CAfreqs) <- c("CHROM",	"POS",	"N_ALLELES",	"N_CHR","Allele1Info","Allele2Info")
AKfreqs <- read.table('/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/OregonCompareToModernAlleleFreqs/modern_alleleFreqsInCA_AK_atOregonSites/AK.AlleleFreqs.frq',header=F,skip=1,sep="\t")
colnames(AKfreqs) <- c("CHROM",	"POS",	"N_ALLELES",	"N_CHR","Allele1Info","Allele2Info")
######### get alelle freq info:
AKfreqs$Allele1 <- unlist(lapply(strsplit(as.character(AKfreqs$Allele1Info),":"),"[",1))
AKfreqs$Allele1Freq <- as.numeric(unlist(lapply(strsplit(as.character(AKfreqs$Allele1Info),":"),"[",2)))
AKfreqs$Allele2 <- unlist(lapply(strsplit(as.character(AKfreqs$Allele2Info),":"),"[",1))
AKfreqs$Allele2Freq <- as.numeric(unlist(lapply(strsplit(as.character(AKfreqs$Allele2Info),":"),"[",2)))
AKfreqs$population <- "AK"
### CA:
CAfreqs$Allele1 <- unlist(lapply(strsplit(as.character(CAfreqs$Allele1Info),":"),"[",1))
CAfreqs$Allele1Freq <- as.numeric(unlist(lapply(strsplit(as.character(CAfreqs$Allele1Info),":"),"[",2)))
CAfreqs$Allele2 <- unlist(lapply(strsplit(as.character(CAfreqs$Allele2Info),":"),"[",1))
CAfreqs$Allele2Freq <- as.numeric(unlist(lapply(strsplit(as.character(CAfreqs$Allele2Info),":"),"[",2)))
CAfreqs$population <- "CA"

CA_AK_Freq <- merge(CAfreqs,AKfreqs,by=c("CHROM","POS"), suffixes = c(".CA",".AK"))
#CA_AK_Freq <- rbind(CAfreqs,AKfreqs)
# want to omit any NaN sites where I don't have data in both pops
CA_AK_Freq_noMissing <- na.omit(CA_AK_Freq)
CA_AK_Freq_noMissing$siteID <- paste(CA_AK_Freq_noMissing$CHROM,"_",CA_AK_Freq_noMissing$POS,sep="")


########### add in oregon pseudohaps at those sites #########

oregonPseudoHaps <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/OregonCompareToModernAlleleFreqs/20191216-ORs-pseudohaps/intersectWithModernData/Oregon.Modern.mfur.intersection.superfile.bed",header=T,comment.char = "&") # intersected with modern data (making skip char something other than "#" because I want the header)
oregonPseudoHaps_usefulInfo <- oregonPseudoHaps[,c("chr","pos","ind0","ind1")]
colnames(oregonPseudoHaps_usefulInfo) <- c("CHROM","POS","OR1","OR2")
CA_AK_OR <- merge(CA_AK_Freq_noMissing,oregonPseudoHaps_usefulInfo,by=c("CHROM","POS"))
head(CA_AK_OR)

##################### make sure alleles are the same ##############
CA_AK_OR[CA_AK_OR$Allele1.CA!=CA_AK_OR$Allele1.AK,] # okay cool they are
CA_AK_OR[CA_AK_OR$Allele2.CA!=CA_AK_OR$Allele2.AK,]
######## allele 1 and 2 labels are the same for CA and AK

# so can subtract allele frequencies 
#CA_AK_OR$AK_CA_Allele1FreqDiff <- abs(CA_AK_OR$Allele1Freq.AK - CA_AK_OR$Allele1Freq.CA)
# remove alelles where allele freq dif is 0:
#CA_AK_OR_variable <- CA_AK_OR[CA_AK_OR$AK_CA_Allele1FreqDiff!=0,]

#require(ggplot2)
#CA_AK_OR$siteID <- paste(CA_AK_OR$CHROM,"_",CA_AK_OR$POS,sep="")

# pie chart plot?

#ggplot(CA_AK_OR,aes(x=siteID,y=Allele1Freq.CA,color=OR1))+
 # geom_point()

# make a nice table:

# label transitions/transversions/triallelic?
CA_AK_OR$allelePair <- paste(CA_AK_OR$Allele1.CA,CA_AK_OR$Allele2.CA,sep=',')
transitions <- c("C,T","T,C","G,A","A,G")
CA_AK_OR$TiTv <- "Tv"
CA_AK_OR[CA_AK_OR$allelePair %in% transitions,]$TiTv <- "Ti"
CA_AK_OR$diff <- abs(CA_AK_OR$Allele1Freq.AK-CA_AK_OR$Allele1Freq.CA)
write_clip(CA_AK_OR[,c("siteID","Allele1.CA","Allele1Freq.CA","Allele1Freq.AK","Allele2.CA","Allele2Freq.CA","Allele2Freq.AK","OR1","OR2","TiTv")])


###### copy into excel ########
