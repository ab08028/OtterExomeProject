########## plot bootstraps and point estimates based on good windows that passed filters ##########
require(ggplot2)
require(dplyr)
dates=c("20190701-lowcov-AFprior-MajorMinor4","20190701-highcov-AFprior-MajorMinor4")
minGP=0.95
minDepth=2
minInd=1 # note this is min ind with coverage at a *site*; when windows were filtered it required all 9 inds to have data wihtin the window
for(angsdDate in dates){
  data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/compareMisSynDists_withBootstraps/AK-CA_separate/",angsdDate,"/PointEstsPlusBootstraps/",sep="")
  ## need to update infile name to have all the filter info maybe?
  ######## point estimates are from summing up GPs across rescaled good bins
  pointEst=read.table(paste(data.dir,"angsdOut.mappedTomfur.Bootstraps.GPs.ProbCutoff.",minGP,".DepthCutoff.",minDepth,".minInd.",minInd,".",angsdDate,".Modern.Ancient.PointEstimatesBasedonGoodBins.txt",sep=""),header=T)
  head(pointEst)     
  
  boots = read.table(paste(data.dir,"angsdOut.mappedTomfur.Bootstraps.GPs.ProbCutoff.",minGP,".DepthCutoff.",minDepth,".minInd.",minInd,".",angsdDate,".Modern.Ancient.allBoots.txt",sep=""),header=T)
  
  # in step 7bii already did a final rescaling of the derived allele counts in the bootstraps to account for hte fact that slightly different numbers of sites are drawn per bootstrap depending on the windows chosen (all are very close)
  # but also want to do this for homAlt: (ideally do in previous step, this is just a lazy shortcut, but it's the same answer and faster to do now:)
  boots$homAlt_Rescaled <- (boots$homAlt/boots$totalCDSSitesInBoot) * boots$SitesToDraw # rescale by homAlt to get spot on estimates of exaclty the homAlt amount per bootstrap if they drew the exact same amount of sites. Only changes numbers by teeny amnt eg 3742 --> 3740
  
  # order factors:
  boots$Consequence_BroadName <- factor(boots$Consequence_BroadName,levels=c("synonymous_variant","missense_variant","stop_gained"))
  pointEst$Consequence_BroadName <- factor(pointEst$Consequence_BroadName,levels=c("synonymous_variant","missense_variant","stop_gained"))
  
  ###### plot ######
  p1 <- ggplot(boots,aes(x=group,y=derivedAlleles_Rescaled))+
    geom_violin()+
    geom_point(data=pointEst,aes(x=group,y=derivedAlleles))+
    theme_bw()+
    facet_wrap(~sites~Consequence_BroadName,scales="free")+
    ggtitle("Derived Alleles -- calculated from windows of genome passing filters\n1000 bootstraps")
  p1
  ggsave(paste(data.dir,"derivedAlleles.PointPlusBoots.pdf",sep=""),p1,device="pdf",height=5,width=7)
  
  p2 <-  ggplot(boots,aes(x=group,y=homAlt_Rescaled))+
    geom_violin()+
    geom_point(data=pointEst,aes(x=group,y=homAlt))+
    theme_bw()+
    facet_wrap(~sites~Consequence_BroadName,scales="free")+
    ggtitle("Homozygous alternate genotypes -- calculated from windows of genome passing filters\n1000 bootstraps")      
  p2
  ggsave(paste(data.dir,"homAlt.PointPlusBoots.pdf",sep=""),p2,device="pdf",height=5,width=7)
  ##### just plot CA ancient and modern ####
  ###### plot ######
  p3 <- ggplot(boots[boots$group!="Modern-AK",],aes(x=group,y=derivedAlleles_Rescaled))+
    geom_violin()+
    geom_point(data=pointEst[pointEst$group!="Modern-AK",],aes(x=group,y=derivedAlleles))+
    theme_bw()+
    facet_wrap(~sites~Consequence_BroadName,scales="free")+
    ggtitle("Derived Alleles -- calculated from windows of genome passing filters\n1000 bootstraps")
  p3
  ggsave(paste(data.dir,"derivedAlleles.PointPlusBoots.noAK.pdf",sep=""),p3,device="pdf",height=5,width=7)
  
  p4 <- ggplot(boots[boots$group!="Modern-AK",],aes(x=group,y=homAlt_Rescaled))+
    geom_violin()+
    geom_point(data=pointEst[pointEst$group!="Modern-AK",],aes(x=group,y=homAlt))+
    theme_bw()+
    facet_wrap(~sites~Consequence_BroadName,scales="free")+
    ggtitle("Homozygous alternate genotypes -- calculated from windows of genome passing filters\n1000 bootstraps")
  p4
  ggsave(paste(data.dir,"homAlt.PointPlusBoots.noAK.pdf",sep=""),p4,device="pdf",height=5,width=7)
  
  ################# carry out Z test --- FINISH CODING THIS!!! ##############
  ######### derived alleles:
  ######### Significance testing ##########
  # Calculate a z score for bootstrap
  # z = [(point estimate (ancient) - point estimate (modern) )- (truemean1 - truemean2)]/sqrt(variance (from bootstrap) 1+variance (from bootstrap) 2) 
  # (truemean1 - truemean2) <-- this part gets set to 0 (assume true mean1=truemean2)
  # and you don't divide variance by n. 
  ### do for each category (syn, mis, sg)
  # and do for Ti and Tv
  # get variance:
  # testing dplyr grouping: 
  #dim(boots[boots$group=="Ancient" & boots$sites=="TvOnly" & boots$Consequence_BroadName=="synonymous_variant",]) # should be 1000; yes that's right
  #var(boots[boots$group=="Ancient" & boots$sites=="TvOnly" & boots$Consequence_BroadName=="synonymous_variant",]$derivedAlleles_Rescaled)
  
  variances <- boots %>% 
    group_by(group,sites,Consequence_BroadName) %>%
    summarise(derivedAlleles_Var=var(derivedAlleles_Rescaled),homAlt_Var=var(homAlt))
  variances
  
  ### make a new df that has all the info: #need to sep mod AK and CA out:
  variance_anc <- variances[variances$group=="Ancient",]
  variance_mod_CA <- variances[variances$group=="Modern-CA",]
  variance_mod_AK <- variances[variances$group=="Modern-AK",]
  
  pointEst_anc <- pointEst[pointEst$group=="Ancient",c("sites","Consequence_BroadName","homAlt","derivedAlleles")]
  pointEst_mod_CA <- pointEst[pointEst$group=="Modern-CA",c("sites","Consequence_BroadName","homAlt","derivedAlleles")]
  pointEst_mod_AK <- pointEst[pointEst$group=="Modern-AK",c("sites","Consequence_BroadName","homAlt","derivedAlleles")]
  
  # put together point ests and variance for modern/ancient:
  modCombo_CA <- merge(pointEst_mod_CA,variance_mod_CA,by=c("sites","Consequence_BroadName"))
  modCombo_AK <- merge(pointEst_mod_AK,variance_mod_AK,by=c("sites","Consequence_BroadName"))
  ancCombo <- merge(pointEst_anc,variance_anc,by=c("sites","Consequence_BroadName"))
  
  # then need to merge the two combos:
  # this now has ancient and modern data combined with point ests from bins and bootstraps:
  # maybe separate Ca and AK combos? Do I even want to test ak? or just CA?
  allCombo1 <- merge(modCombo_CA,ancCombo,by=c("sites","Consequence_BroadName"),suffixes=c(".modCA",".anc"))
  
  ################## YOU ARE STUCK HERE!!!!!!!!!!!! NEED TO REDO Z TEST #######
  allCombo2 <- merge(allCombo1,modCombo_AK,by=c("sites","Consequence_BroadName"),suffixes=c("",".modAK"))
  
  head(allCombo)
  head(ancCombo)
  # this worked!
  
  # now calculate Z:
  # first the diff in point estimates for derived alleles and then hom Alt (anc - mod )
  ### NUMERATORs OF Z score: ####
  allCombo$derivedAlleles_xDiff <- allCombo$derivedAlleles.anc - allCombo$derivedAlleles.mod  # diff in point estimates 
  allCombo$homAlt_xDiff <- allCombo$homAlt.anc - allCombo$homAlt.mod  # diff in point estimates 
  ##### DENOMINATORS OF Z SCORE: ####
  # sqrt(VARIANCE1+VARIANCE2) (recall that variance is written as sigma squared so you don't need to additionally square the variance.)
  allCombo$derivedAlleles_sqVarSum <- sqrt(allCombo$derivedAlleles_Var.mod + allCombo$derivedAlleles_Var.anc)
  allCombo$homAlt_sqVarSum <- sqrt(allCombo$homAlt_Var.mod + allCombo$homAlt_Var.anc)
  
  ########## Z score and two sided p values: ##########
  allCombo$derivedAlleles_ZScore <- allCombo$derivedAlleles_xDiff / allCombo$derivedAlleles_sqVarSum
  allCombo$derivedAlleles_PValue2Sided <- 2*pnorm(-abs(allCombo$derivedAlleles_ZScore))
  
  allCombo$homAlt_ZScore <- allCombo$homAlt_xDiff / allCombo$homAlt_sqVarSum
  allCombo$homAlt_PValue2Sided <- 2*pnorm(-abs(allCombo$homAlt_ZScore))
  
  ####### write out table #######
  write.table(allCombo,paste(data.dir,"ZtestCalculations.",angsdDate,".txt",sep=""),quote=F,row.names = F,sep="\t")
}

