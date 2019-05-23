require(reshape2)
require(ggplot2)
require(dplyr)
################### comparing GL (geno lhood), GP (geno posterior), and MAF ############
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/Heterozygosity/experiments/compareGL-GP-MAF/"
sampleListDir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/" # BE VERY CAREFUL WHAT LIST YOU USE; IT MUST MATCH THE LIST USED IN ANGSD!
rundate="20190522"

######## dataset 1: high coverage (downsampled) modern + aDNA just for scaffold GL896898.1 ############
state="highcov"
GLs=read.table(paste(data.dir,rundate,"-experiment-GLvsGP-",state,"/posteriorProbabilities/angsdOut.mappedTomfur.OrlandoSettings.beagle.gz",sep=""),header=T)
GPs=read.table(paste(data.dir,rundate,"-experiment-GLvsGP-",state,"/posteriorProbabilities/angsdOut.mappedTomfur.OrlandoSettings.beagle.gprobs.gz",sep=""),header=T)
#MAFs=read.table(paste(data.dir,rundate,"-experiment-GLvsGP-",state,"/posteriorProbabilities/angsdOut.mappedTomfur.OrlandoSettings.mafs.gz",sep=""),header=T)

######### melt each one #########
GLs_melt <- melt(GLs,id.vars = c("marker"  ,"allele1", "allele2")) 
colnames(GLs_melt) <- c("marker"  ,"allele1", "allele2","Ind_GT","GL")

GPs_melt <- melt(GPs,id.vars = c("marker"  ,"allele1", "allele2")) 
colnames(GPs_melt) <- c("marker"  ,"allele1", "allele2","Ind_GT","GP")
# try just two individuals to start with?
#inds=c("A30_Elut_CA_SM_35_SN1_CAP","116_Elut_CA_307_downsamp") # these are Ind2 and Ind3

# want to get mean across ancient and modern individuals 
#Ind0 Ind1 Ind2 are ancient; Ind3 - Ind 8 are modern#
GPs_melt$label <- NA

GPs_melt[grepl("Ind[0-2]",GPs_melt$Ind_GT),]$label <- "ancient" # label the ancient samples
GPs_melt[grepl("Ind[3-8]",GPs_melt$Ind_GT),]$label <- "modern" # label the modern samples

# label the different genotypes
GPs_melt$GTLabel <- NA
GPs_melt[!grepl("\\.",GPs_melt$Ind_GT),]$GTLabel <- "0/0" # Ind with no "." are hom ref (originally the first column for each individual)
GPs_melt[grepl("\\.1",GPs_melt$Ind_GT),]$GTLabel <- "0/1" # X.1 are heterozygotes (originally the second column for each individual)
GPs_melt[grepl("\\.2",GPs_melt$Ind_GT),]$GTLabel <- "1/1" # X.2 are hom alt (originally the third column for each individual)

# get the mean across ancient and modern invididuals for each site (marker) grouped by ancient modern and by genotype
meanGP_modernAncient <- GPs_melt %>%
  group_by(label,marker,GTLabel) %>%
  summarise(meanGP = mean(GP))

p1 <- ggplot(meanGP_modernAncient,aes(x=meanGP,color=GTLabel))+
  geom_density()+
  facet_wrap(~label)
# then do the same for GLs
# want to get mean across ancient and modern individuals 
#Ind0 Ind1 Ind2 are ancient; Ind3 - Ind 8 are modern#
GLs_melt$label <- NA

GLs_melt[grepl("Ind[0-2]",GLs_melt$Ind_GT),]$label <- "ancient" # label the ancient samples
GLs_melt[grepl("Ind[3-8]",GLs_melt$Ind_GT),]$label <- "modern" # label the modern samples

# label the different genotypes
GLs_melt$GTLabel <- NA
GLs_melt[!grepl("\\.",GLs_melt$Ind_GT),]$GTLabel <- "0/0" # Ind with no "." are hom ref (originally the first column for each individual)
GLs_melt[grepl("\\.1",GLs_melt$Ind_GT),]$GTLabel <- "0/1" # X.1 are heterozygotes (originally the second column for each individual)
GLs_melt[grepl("\\.2",GLs_melt$Ind_GT),]$GTLabel <- "1/1" # X.2 are hom alt (originally the third column for each individual)

meanGL_modernAncient <- GLs_melt %>%
  group_by(label,marker,GTLabel) %>%
  summarise(meanGL = mean(GL))
# get the mean across ancient and modern invididuals for each site (marker) grouped by ancient modern and by genotype
#p1 <- ggplot(meanGL_modernAncient,aes(x=meanGL,color=GTLabel))+
#  geom_density()+
#  facet_wrap(~label)

#### combine GL and GP means: 
meanGL_GP = inner_join(meanGL_modernAncient,meanGP_modernAncient)
# pretty fast
###### plot: ########
#GLvsGP_plot <- ggplot(head(meanGL_GP),aes(x=meanGL,y=meanGP,color=GTLabel))+
#  geom_point()+
#  facet_wrap(~label)
#GLvsGP_plot
#ggsave(paste(data.dir,"meanGL.vs.meanGP.OneScaffold.png",sep=""),GLvsGP_plot,width=10,height=7)

GLvsGP_plot2 <- ggplot(meanGL_GP,aes(x=meanGL,y=meanGP,color=label))+
  geom_point()+
  facet_wrap(~GTLabel)+
  geom_abline(slope=1,intercept=0)+
  ggtitle(paste("Based on ",state," dataset"))
GLvsGP_plot2
ggsave(paste(data.dir,state,".meanGL.vs.meanGP.OneScaffold.2.png",sep=""),GLvsGP_plot2,width=10,height=7)

############## incorporate MAFs? ###############

