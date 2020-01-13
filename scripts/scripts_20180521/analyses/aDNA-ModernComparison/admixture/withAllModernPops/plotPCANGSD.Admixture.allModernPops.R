library(RcppCNPy)
require(ggplot2)
require(reshape2)
require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[1],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5])
######### testing admixture
# binary npy file needs to be read as:
minMafs=c("0.025","0.05","0.12","0.2")
refs=c("elut","mfur")
eStates=c("") # skipping e5 (did previously)
# to do: fix order of individuals
#
#################################### high coverage #######################################
date="20191212-highcov-AFprior-MajorMinor4-plusCOM-KUR-AL"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/admixture/",date,"/",sep="")
# be sure to pick right sample list!!! #
sampleList=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/angsd.bamList.mappedtoMfurfullpaths.HighCovPlusADNA.PlusCOM.KUR.AL.IDsOnly.txt") # same list of IDs for mapped to Mfur or Elut (Same order) so just using Mfur version (same)
for(ref in refs){
  for(minMaf in minMafs){
    for(eState in eStates){
    state="1e-06.snpsOnly.TransvOnly"
    inputfilename <- paste(data.dir,"pcAngsd.",ref,".",state,".minMaf.",minMaf,eState,".admix.Q.npy",sep="")
    admixQ  <- npyLoad(inputfilename) # ready from binary ## NOTE; this test dataset is only based on a tiny part of a chromosome. Only 4 snps! I think this test was mapped to elut. let's see
    # plot loadings
    admixQdf <- data.frame(admixQ)
    colnames(sampleList) <- "sampleID"
    admixQdf <- cbind(admixQdf,sampleList)
    head(admixQdf)
    #
    # order individuals by population
    admixQdf$population <- unlist(lapply(strsplit(as.character(admixQdf$sampleID),"_"),"[",3))
    admixQdf[grep("^A",admixQdf$sampleID),]$population <- "Anc-CA"
    admixQdf[grep("BER",admixQdf$sampleID),]$population <- "COM"
    admixQdf_melt <- melt(admixQdf)
    #FIX LEVELS:
    admixQdf_melt$sampleID <- factor(admixQdf_melt$sampleID, levels=c("A29_Elut_CA_SM_30_SN2_CAP","A13_Elut_CA_AN_388_SN1_2CAP_screen","A30_Elut_CA_SM_35_SN1_CAP","116_Elut_CA_307","141_Elut_CA_419","140_Elut_CA_403","129_Elut_AK_AL4660","126_Elut_AK_AF3394","55_Elut_AK_AF3736","111_Elut_AL_AT_GE91133","121_Elut_AL_AT_GE91135","124_Elut_AL_AT_GE91143","136_Elut_BER_46","137_Elut_BER_88","50_Elut_BER_100","100_Elut_KUR_24","101_Elut_KUR_3","102_Elut_KUR_4"))
   
    admixQdf_melt$population <- factor(admixQdf_melt$population,levels=c("Anc-CA","CA","AK","AL","COM","KUR"))
    p1 <- ggplot(admixQdf_melt,aes(x=sampleID,y=value,fill=variable))+
      geom_col(position="stack")+
      theme_bw()+
      theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank())+
      ggtitle(paste("PCAngsd Admixture\n",date,"\nref = ",ref,"; ",state,"; minMaf: ",minMaf,sep=""))+
      ylab("")+
      xlab("")+
      facet_wrap(~population,nrow=1,scales="free_x")+
      theme(panel.border = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.spacing.x = unit(0.1,"line"),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.background = element_rect("transparent"),
            strip.text = element_text(size=11))
      
      #+
      #scale_fill_manual(values=c(V1=colorPal[1],V4=colorPal[2],V5=colorPal[3],V3=colorPal[4],V2=colorPal[5])) ### adjust
    p1
    ggsave(paste(data.dir,"admixturePlot.ref.",ref,".minMaf.",minMaf,eState,".pdf",sep=""),p1,height=4,width=5)
    # Number of sites after MAF filtering (0.12): 14499 <--- how to get this number efficiently? this is output to screen during pcangsd 
    # want to deal with missing data better -- maybe pre-filtering beagle file to have (some? no? missing data)
    
    
    ######### formattted for manuscript ##########
    p1b <- ggplot(admixQdf_melt,aes(x=sampleID,y=value,fill=variable))+
      geom_col(position="stack")+
      theme_bw()+
      theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank())+
      #ggtitle(paste("PCAngsd Admixture\n",date,"\nref = ",ref,"; ",state,"; minMaf: ",minMaf,sep=""))+
      ylab("")+
      xlab("")+
      facet_wrap(~population,nrow=1,scales="free_x")+
      theme(panel.border = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.spacing.x = unit(0.1,"line"),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.background = element_rect("transparent"),
            strip.text = element_text(size=12),
            axis.text.y=element_text(size=14))
    
    #+
    #scale_fill_manual(values=c(V1=colorPal[1],V4=colorPal[2],V5=colorPal[3],V3=colorPal[4],V2=colorPal[5])) ### adjust
    p1b
    ggsave(paste(data.dir,"admixturePlot.ref.",ref,".minMaf.",minMaf,eState,".formattedForManuscript.pdf",sep=""),p1b,height=3.6,width=4.7)
    
  }}}

######### for manuscript ##########

