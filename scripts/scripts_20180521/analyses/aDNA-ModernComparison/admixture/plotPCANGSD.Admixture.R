library(RcppCNPy)
require(ggplot2)
require(reshape2)
######### testing admixture
# binary npy file needs to be read as:
minMafs=c("0.05","0.12","0.2")
refs=c("elut","mfur")
#################################### high coverage #######################################
date="20190701-highcov-AFprior-MajorMinor4"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/admixture/",date,"/",sep="")
# be sure to pick right sample list!!! #
sampleList=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt")
for(ref in refs){
  for(minMaf in minMafs){
    state="1e-06.snpsOnly.TransvOnly"
    inputfilename <- paste(data.dir,"pcAngsd.",ref,".",state,".minMaf.",minMaf,".admix.Q.npy",sep="")
    admixQ  <- npyLoad(inputfilename) # ready from binary ## NOTE; this test dataset is only based on a tiny part of a chromosome. Only 4 snps! I think this test was mapped to elut. let's see
    # plot loadings
    admixQdf <- data.frame(admixQ)
    colnames(sampleList) <- "sampleID"
    admixQdf <- cbind(admixQdf,sampleList)
    head(admixQdf)
    #
    # order individuals by population
    admixQdf$population <- unlist(lapply(strsplit(as.character(admixQdf$sampleID),"_"),"[",3))
    admixQdf[grep("^A",admixQdf$sampleID),]$population <- "Ancient-CA"
    admixQdf_melt <- melt(admixQdf)
    admixQdf_melt$sampleID <- factor(admixQdf_melt$sampleID, levels=c("129_Elut_AK_AL4660","126_Elut_AK_AF3394","55_Elut_AK_AF3736","116_Elut_CA_307","141_Elut_CA_419","140_Elut_CA_403","A29_Elut_CA_SM_30_SN2_CAP","A13_Elut_CA_AN_388_SN1_2CAP_screen","A30_Elut_CA_SM_35_SN1_CAP"))
    admixQdf_melt$population <- factor(admixQdf_melt$population,levels=c("AK","CA","Ancient-CA"))
    p1 <- ggplot(admixQdf_melt,aes(x=sampleID,y=value,fill=variable))+
      geom_col(position="stack")+
      theme_bw()+
      theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank())+
      ggtitle(paste("PCAngsd Admixture\n",date,"\nref = ",ref,"; ",state,"; minMaf: ",minMaf,sep=""))+
      ylab("")+
      xlab("")+
      facet_wrap(~population,nrow=1,scales="free_x")
    p1
    ggsave(paste(data.dir,"admixturePlot.ref.",ref,".minMaf.",minMaf,".pdf",sep=""),height=5,width=9)
    # Number of sites after MAF filtering (0.12): 14499 <--- how to get this number efficiently? this is output to screen during pcangsd 
    # want to deal with missing data better -- maybe pre-filtering beagle file to have (some? no? missing data)
  }}

#################################### high coverage -- no missing data #######################################
date="20191210-highcov-AFprior-MajorMinor4-noMissing"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/admixture/",date,"/",sep="")
# be sure to pick right sample list!!! #
sampleList=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt")
for(ref in refs){
  for(minMaf in minMafs){
    state="1e-06.snpsOnly.TransvOnly.noMiss"
    inputfilename <- paste(data.dir,"pcAngsd.",ref,".",state,".minMaf.",minMaf,".admix.Q.npy",sep="")
    admixQ  <- npyLoad(inputfilename) # ready from binary ## NOTE; this test dataset is only based on a tiny part of a chromosome. Only 4 snps! I think this test was mapped to elut. let's see
    # plot loadings
    admixQdf <- data.frame(admixQ)
    colnames(sampleList) <- "sampleID"
    admixQdf <- cbind(admixQdf,sampleList)
    head(admixQdf)
    #
    # order individuals by population
    admixQdf$population <- unlist(lapply(strsplit(as.character(admixQdf$sampleID),"_"),"[",3))
    admixQdf[grep("^A",admixQdf$sampleID),]$population <- "Ancient-CA"
    admixQdf_melt <- melt(admixQdf)
    admixQdf_melt$sampleID <- factor(admixQdf_melt$sampleID, levels=c("129_Elut_AK_AL4660","126_Elut_AK_AF3394","55_Elut_AK_AF3736","116_Elut_CA_307","141_Elut_CA_419","140_Elut_CA_403","A29_Elut_CA_SM_30_SN2_CAP","A13_Elut_CA_AN_388_SN1_2CAP_screen","A30_Elut_CA_SM_35_SN1_CAP"))
    admixQdf_melt$population <- factor(admixQdf_melt$population,levels=c("AK","CA","Ancient-CA"))
    p1 <- ggplot(admixQdf_melt,aes(x=sampleID,y=value,fill=variable))+
      geom_col(position="stack")+
      theme_bw()+
      theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank())+
      ggtitle(paste("PCAngsd Admixture\n",date,"\nref = ",ref,"; ",state,"; minMaf: ",minMaf,sep=""))+
      ylab("")+
      xlab("")+
      facet_wrap(~population,nrow=1,scales="free_x")
    p1
    ggsave(paste(data.dir,"admixturePlot.ref.",ref,".minMaf.",minMaf,".pdf",sep=""),height=5,width=9)
    # Number of sites after MAF filtering (0.12): 14499 <--- how to get this number efficiently? this is output to screen during pcangsd 
    # want to deal with missing data better -- maybe pre-filtering beagle file to have (some? no? missing data)
  }}

#################################### low coverage #######################################
date="20190701-lowcov-AFprior-MajorMinor4"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/admixture/",date,"/",sep="")
# be sure to pick right sample list!!
sampleList=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.LowCoverageOnly.BeCarefulOfOrder.txt")
for(ref in refs){
  for(minMaf in minMafs){
    state="1e-06.snpsOnly.TransvOnly"
    inputfilename <- paste(data.dir,"pcAngsd.",ref,".",state,".minMaf.",minMaf,".admix.Q.npy",sep="")
    admixQ  <- npyLoad(inputfilename) # ready from binary ## NOTE; this test dataset is only based on a tiny part of a chromosome. Only 4 snps! I think this test was mapped to elut. let's see
    # plot loadings
    admixQdf <- data.frame(admixQ)
    colnames(sampleList) <- "sampleID"
    admixQdf <- cbind(admixQdf,sampleList)
    head(admixQdf)
    #
    # order individuals by population
    admixQdf$population <- unlist(lapply(strsplit(as.character(admixQdf$sampleID),"_"),"[",3))
    admixQdf[grep("^A",admixQdf$sampleID),]$population <- "Ancient-CA"
    admixQdf_melt <- melt(admixQdf)
    admixQdf_melt$sampleID <- factor(admixQdf_melt$sampleID, levels=c("129_Elut_AK_AL4660_downsamp","126_Elut_AK_AF3394_downsamp","55_Elut_AK_AF3736_downsamp","116_Elut_CA_307_downsamp","141_Elut_CA_419_downsamp","140_Elut_CA_403_downsamp","A29_Elut_CA_SM_30_SN2_CAP","A13_Elut_CA_AN_388_SN1_2CAP_screen","A30_Elut_CA_SM_35_SN1_CAP"))
    admixQdf_melt$population <- factor(admixQdf_melt$population,levels=c("AK","CA","Ancient-CA"))
    p1 <- ggplot(admixQdf_melt,aes(x=sampleID,y=value,fill=variable))+
      geom_col(position="stack")+
      theme_bw()+
      theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank())+
      ggtitle(paste("PCAngsd Admixture\n",date,"\nref = ",ref,"; ",state,"; minMaf: ",minMaf,sep=""))+
      ylab("")+
      xlab("")+
      facet_wrap(~population,nrow=1,scales="free_x")
    p1
    ggsave(paste(data.dir,"admixturePlot.ref.",ref,".minMaf.",minMaf,".pdf",sep=""),height=5,width=9)
    # Number of sites after MAF filtering (0.12): 14499 <--- how to get this number efficiently? this is output to screen during pcangsd 
    # want to deal with missing data better -- maybe pre-filtering beagle file to have (some? no? missing data)
  }}


#################################### low coverage -- no missing data #######################################
date="20191210-lowcov-AFprior-MajorMinor4-noMissing"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/admixture/",date,"/",sep="")
# be sure to pick right sample list!!
sampleList=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.LowCoverageOnly.BeCarefulOfOrder.txt")
for(ref in refs){
  for(minMaf in minMafs){
    state="1e-06.snpsOnly.TransvOnly.noMiss"
    inputfilename <- paste(data.dir,"pcAngsd.",ref,".",state,".minMaf.",minMaf,".admix.Q.npy",sep="")
    admixQ  <- npyLoad(inputfilename) # ready from binary ## NOTE; this test dataset is only based on a tiny part of a chromosome. Only 4 snps! I think this test was mapped to elut. let's see
    # plot loadings
    admixQdf <- data.frame(admixQ)
    colnames(sampleList) <- "sampleID"
    admixQdf <- cbind(admixQdf,sampleList)
    head(admixQdf)
    #
    # order individuals by population
    admixQdf$population <- unlist(lapply(strsplit(as.character(admixQdf$sampleID),"_"),"[",3))
    admixQdf[grep("^A",admixQdf$sampleID),]$population <- "Ancient-CA"
    admixQdf_melt <- melt(admixQdf)
    admixQdf_melt$sampleID <- factor(admixQdf_melt$sampleID, levels=c("129_Elut_AK_AL4660_downsamp","126_Elut_AK_AF3394_downsamp","55_Elut_AK_AF3736_downsamp","116_Elut_CA_307_downsamp","141_Elut_CA_419_downsamp","140_Elut_CA_403_downsamp","A29_Elut_CA_SM_30_SN2_CAP","A13_Elut_CA_AN_388_SN1_2CAP_screen","A30_Elut_CA_SM_35_SN1_CAP"))
    admixQdf_melt$population <- factor(admixQdf_melt$population,levels=c("AK","CA","Ancient-CA"))
    p1 <- ggplot(admixQdf_melt,aes(x=sampleID,y=value,fill=variable))+
      geom_col(position="stack")+
      theme_bw()+
      theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank())+
      ggtitle(paste("PCAngsd Admixture\n",date,"\nref = ",ref,"; ",state,"; minMaf: ",minMaf,sep=""))+
      ylab("")+
      xlab("")+
      facet_wrap(~population,nrow=1,scales="free_x")
    p1
    ggsave(paste(data.dir,"admixturePlot.ref.",ref,".minMaf.",minMaf,".pdf",sep=""),height=5,width=9)
    # Number of sites after MAF filtering (0.12): 14499 <--- how to get this number efficiently? this is output to screen during pcangsd 
    # want to deal with missing data better -- maybe pre-filtering beagle file to have (some? no? missing data)
  }}



#################################### experiment (too few data!): with extra individuals (elut only) ###########
refs="elut" # only for now
date="20191205-highcov-AFprior-MajorMinor4-IncludesMoreInds-notAllIndelRealignedYet-noA7"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/admixture/",date,"/",sep="")
sampleList=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/angsd.bamList.mappedtoElutfullpaths.HighCovPlusADNAOnly.PLUSEXTRAADNA.notIndelRealYet.IDsOnly.noA7.txt")
for(ref in refs){
  for(minMaf in minMafs){
    state="1e-06.snpsOnly.TransvOnly"
    inputfilename <- paste(data.dir,"pcAngsd.",ref,".",state,".minMaf.",minMaf,".admix.Q.npy",sep="")
    admixQ  <- npyLoad(inputfilename) # ready from binary ## NOTE; this test dataset is only based on a tiny part of a chromosome. Only 4 snps! I think this test was mapped to elut. let's see
    # plot loadings
    admixQdf <- data.frame(admixQ)
    colnames(sampleList) <- "sampleID"
    admixQdf <- cbind(admixQdf,sampleList)
    head(admixQdf)
    #
    # order individuals by population
    admixQdf$population <- unlist(lapply(strsplit(as.character(admixQdf$sampleID),"_"),"[",3))
    #admixQdf[grep("^A",admixQdf$sampleID),]$population <- "Ancient-CA"
    admixQdf_melt <- melt(admixQdf)
    #admixQdf_melt$sampleID <- factor(admixQdf_melt$sampleID, levels=c("129_Elut_AK_AL4660","126_Elut_AK_AF3394","55_Elut_AK_AF3736","116_Elut_CA_307","141_Elut_CA_419","140_Elut_CA_403","A29_Elut_CA_SM_30_SN2_CAP","A13_Elut_CA_AN_388_SN1_2CAP_screen","A30_Elut_CA_SM_35_SN1_CAP"))
    #admixQdf_melt$population <- factor(admixQdf_melt$population,levels=c("AK","CA","Ancient-CA"))
    p1 <- ggplot(admixQdf_melt,aes(x=sampleID,y=value,fill=variable))+
      geom_col(position="stack")+
      theme_bw()+
      theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank())+
      ggtitle(paste("PCAngsd Admixture\n",date,"\nref = ",ref,"; ",state,"; minMaf: ",minMaf,sep=""))+
      ylab("")+
      xlab("")+
      facet_wrap(~population,nrow=1,scales="free_x")
    p1
    ggsave(paste(data.dir,"admixturePlot.ref.",ref,".minMaf.",minMaf,".pdf",sep=""),height=3,width=7)
  }
}

