######## you are here ######

########## plot bootstraps and point estimates based on good windows that passed filters ##########
require(ggplot2)
require(dplyr)
require(RColorBrewer)
# not loving the color scheme.
colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[7],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors
dates=c("20190701-lowcov-AFprior-MajorMinor4","20190701-highcov-AFprior-MajorMinor4")
#dates=c("20190701-lowcov-AFprior-MajorMinor4","20190701-highcov-AFprior-MajorMinor4")
minGP=0.95
minDepth=2
minInd=1 # note this is min ind with coverage at a *site*; when windows were filtered it required all 9 inds to have data wihtin the window
for(angsdDate in dates){
  data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/compareMisSynDists_withBootstraps/AK-CA-separate/",angsdDate,"/PointEstsPlusBootstraps/",sep="")
  ## need to update infile name to have all the filter info maybe?
  ######## point estimates are from summing up GPs across rescaled good bins
  pointEst=read.table(paste(data.dir,"angsdOut.mappedTomfur.Bootstraps.GPs.ProbCutoff.",minGP,".DepthCutoff.",minDepth,".minInd.",minInd,".",angsdDate,".Modern.Ancient.PointEstimatesBasedonGoodBins.txt",sep=""),header=T)
  head(pointEst)     
  ### need to get avg per group between 3 inds per group:
  averagePointEsts <- pointEst %>%
    group_by(group,sites,Consequence_BroadName) %>%
    summarise(meanPtEstHomRef=mean(homRef),meanPtEstHet=mean(het),meanPtEstHomAlt=mean(homAlt),meanPtEstDerivedAlleles=mean(derivedAlleles))
  boots = read.table(paste(data.dir,"angsdOut.mappedTomfur.Bootstraps.GPs.ProbCutoff.",minGP,".DepthCutoff.",minDepth,".minInd.",minInd,".",angsdDate,".perGroupAveragesBootstraps.useThis.txt",sep=""),header=T)
  
  # in step 7bii already did a final rescaling of the derived allele counts in the bootstraps to account for hte fact that slightly different numbers of sites are drawn per bootstrap depending on the windows chosen (all are very close)
  # but also want to do this for homAlt: (ideally do in previous step, this is just a lazy shortcut, but it's the same answer and faster to do now:)
  boots$homAlt_Rescaled <- (boots$homAlt/boots$totalCDSSitesInBoot) * boots$SitesToDraw # rescale by homAlt to get spot on estimates of exaclty the homAlt amount per bootstrap if they drew the exact same amount of sites. Only changes numbers by teeny amnt eg 3742 --> 3740
  
  # order factors:
  boots$Consequence_BroadName <- factor(boots$Consequence_BroadName,levels=c("synonymous_variant","missense_variant","stop_gained"))
  averagePointEsts$Consequence_BroadName <- factor(averagePointEsts$Consequence_BroadName,levels=c("synonymous_variant","missense_variant","stop_gained"))
  
  ###### plot ######
  p1 <- ggplot(boots,aes(x=group,y=derivedAlleles_Rescaled))+
    geom_violin(alpha=0.7,aes(fill=group))+
    geom_point(data=averagePointEsts,aes(x=group,y=meanPtEstDerivedAlleles))+
    theme_bw()+
    scale_fill_manual(values=c("#1B9E77","#D95F02", "#1B9E77"))+
    theme(legend.position = "none")+
    facet_wrap(~sites~Consequence_BroadName,scales="free")+
    ggtitle("Derived Alleles -- calculated from windows of genome passing filters\n1000 bootstraps\nNew approach: averaging across individuals per bootstrap not per window")
  p1
  ggsave(paste(data.dir,"derivedAlleles.PointPlusBoots.pdf",sep=""),p1,device="pdf",height=6,width=9)
  
  p2 <-  ggplot(boots,aes(x=group,y=homAlt_Rescaled))+
    geom_violin(alpha=0.7,aes(fill=group))+
    geom_point(data=averagePointEsts,aes(x=group,y=meanPtEstHomAlt))+
    theme_bw()+
    scale_fill_manual(values=c("#1B9E77","#D95F02", "#1B9E77"))+
    theme(legend.position = "none")+
    facet_wrap(~sites~Consequence_BroadName,scales="free")+
    ggtitle("Homozygous alternate genotypes -- calculated from windows of genome passing filters\n1000 bootstraps\nNew approach: averaging across individuals per bootstrap not per window")      
  p2
  ggsave(paste(data.dir,"homAlt.PointPlusBoots.pdf",sep=""),p2,device="pdf",height=6,width=9)
  
  ################# carry out Z test ##############
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
    summarise(derivedAlleles_Var=var(derivedAlleles_Rescaled),homAlt_Var=var(homAlt_Rescaled))
  variances
  
  ### make a new df that has all the info:
  variance_anc <- variances[variances$group=="Ancient",]
  variance_modAK <- variances[variances$group=="Modern-AK",]
  variance_modCA <- variances[variances$group=="Modern-CA",]
  
  avg_pointEst_anc <- averagePointEsts[averagePointEsts$group=="Ancient",c("sites","Consequence_BroadName","meanPtEstHomAlt","meanPtEstDerivedAlleles")]
  avg_pointEst_modAK <- averagePointEsts[averagePointEsts$group=="Modern-AK",c("sites","Consequence_BroadName","meanPtEstHomAlt","meanPtEstDerivedAlleles")]
  avg_pointEst_modCA <- averagePointEsts[averagePointEsts$group=="Modern-CA",c("sites","Consequence_BroadName","meanPtEstHomAlt","meanPtEstDerivedAlleles")]
  # put together point ests and variance for modern/ancient:
  modComboCA <- merge(avg_pointEst_modCA,variance_modCA,by=c("sites","Consequence_BroadName"))
  modComboAK <- merge(avg_pointEst_modAK,variance_modAK,by=c("sites","Consequence_BroadName"))
  ancCombo <- merge(avg_pointEst_anc,variance_anc,by=c("sites","Consequence_BroadName"))
  
  # then need to merge the two combos:
  # this now has ancient and modern data combined with point ests from bins and bootstraps:
  CACombo <- merge(modComboCA,ancCombo,by=c("sites","Consequence_BroadName"),suffixes=c(".1",".2"))
  AKCombo <-  merge(modComboAK,ancCombo,by=c("sites","Consequence_BroadName"),suffixes=c(".1",".2"))
  ModCombo <-  merge(modComboCA,modComboAK,by=c("sites","Consequence_BroadName"),suffixes=c(".1",".2"))
  #allCombo <- merge(modCombo,ancCombo,by=c("sites","Consequence_BroadName"),suffixes=c(".mod",".anc"))
  combos = list(CACombo,AKCombo,ModCombo)
  # now calculate Z:
  # first the diff in point estimates for derived alleles and then hom Alt (anc - mod )
  ##### NUMERATORs OF Z score: ######
  for(allCombo in combos){
  allCombo= data.frame(allCombo)
  group1 = unlist(unique(allCombo$group.1))
  group2 = unlist(unique(allCombo$group.2))
  allCombo$derivedAlleles_xDiff <- allCombo$meanPtEstDerivedAlleles.2 - allCombo$meanPtEstDerivedAlleles.1  # diff in point estimates 
  allCombo$homAlt_xDiff <- allCombo$meanPtEstHomAlt.2 - allCombo$meanPtEstHomAlt.1  # diff in point estimates 
  ##### DENOMINATORS OF Z SCORE: ####
  # sqrt(VARIANCE1+VARIANCE2) (recall that variance is written as sigma squared so you don't need to additionally square the variance.)
  allCombo$derivedAlleles_sqVarSum <- sqrt(allCombo$derivedAlleles_Var.2 + allCombo$derivedAlleles_Var.1)
  allCombo$homAlt_sqVarSum <- sqrt(allCombo$homAlt_Var.2 + allCombo$homAlt_Var.1)
  
  ########## Z score and two sided p values: ##########
  allCombo$derivedAlleles_ZScore <- allCombo$derivedAlleles_xDiff / allCombo$derivedAlleles_sqVarSum
  allCombo$derivedAlleles_PValue2Sided <- 2*pnorm(-abs(allCombo$derivedAlleles_ZScore))
  
  allCombo$homAlt_ZScore <- allCombo$homAlt_xDiff / allCombo$homAlt_sqVarSum
  allCombo$homAlt_PValue2Sided <- 2*pnorm(-abs(allCombo$homAlt_ZScore))
  
  ####### write out table #######
  write.table(allCombo,paste(data.dir,group1,".vs.",group2,".ZtestCalculations.",angsdDate,".txt",sep=""),quote=F,row.names = F,sep="\t")
}}
