########### Plot scatter plots of DP, QD and QUAL for variant sites (invariant sites don't have QD)
########################## Plot DP #########################
genotypeDate="20181119"
scaffold="GL896899.1"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/DP_QD_DIST/",genotypeDate,"/jointDP_QD_QUAL/",sep="")
# this is the result of extract_QD_DP_QUAL.py for one scaffold
# before the min500dp filter:
jointDist_preDPFilter = read.table(paste(data.dir,scaffold,".raw_variants.joint.QD.DP.QUAL.dist.2.txt.gz",sep=""),sep="\t",strip.white = T,header=T)
jointDist_postDPFilter = read.table(paste(data.dir,scaffold,".all_1_TrimAlt_raw_variants.joint.QD.DP.QUAL.dist.txt.gz",sep=""),sep="\t",strip.white = T,header=T)
jointDist_postAllFilters = read.table(paste(data.dir,scaffold,".all_5_passingFilters_raw_variants.joint.QD.DP.QUAL.dist.txt.gz",sep=""),sep="\t",strip.white = T,header=T)
####################### prefiltering ####################################
######## QUAL vs DP ########
p0 <- ggplot(jointDist_preDPFilter,aes(x=QUAL,y=DP,color=QD))+
  geom_point(alpha=0.6)+
  theme_bw()+
  geom_hline(yintercept = 500)
p0
ggsave(paste(data.dir,scaffold,".DP.vs.QUAL.coloredByQD.preDP500.png",sep=""),p0,width=7,height=5)
######## QUAL vs QD ########
p1 <- ggplot(jointDist_preDPFilter,aes(x=QUAL,y=QD,color=DP))+
  geom_point(alpha=0.6)+
  theme_bw()+
  geom_hline(yintercept = 2)
p1
ggsave(paste(data.dir,scaffold,".QUAL.vs.QD.coloredByDP.preDP500.png",sep=""),p1,width=7,height=5)

######## QD vs DP ########
p2 <- ggplot(jointDist_preDPFilter,aes(x=DP,y=QD,color=QUAL))+
  geom_point(alpha=0.6)+
  theme_bw()+
  geom_hline(yintercept = 2)
p2
ggsave(paste(data.dir,scaffold,".DP.vs.QD.coloredByQUAL.preDP500.png",sep=""),p2,width=7,height=5)

####################### postfiltering ####################################
######## QUAL vs DP ########
p0 <- ggplot(jointDist_postDPFilter,aes(x=QUAL,y=DP,color=QD))+
  geom_point(alpha=0.6)+
  theme_bw()+
  geom_hline(yintercept = 500)
p0
ggsave(paste(data.dir,scaffold,".DP.vs.QUAL.coloredByQD.postDP500.png",sep=""),p0,width=7,height=5)
######## QUAL vs QD ########
p1 <- ggplot(jointDist_postDPFilter,aes(x=QUAL,y=QD,color=DP))+
  geom_point(alpha=0.6)+
  theme_bw()+
  geom_hline(yintercept = 2)
p1
ggsave(paste(data.dir,scaffold,".QUAL.vs.QD.coloredByDP.postDP500.png",sep=""),p1,width=7,height=5)

######## QD vs DP ########
p2 <- ggplot(jointDist_postDPFilter,aes(x=DP,y=QD,color=QUAL))+
  geom_point(alpha=0.6)+
  theme_bw()+
  geom_hline(yintercept = 2)
p2
ggsave(paste(data.dir,scaffold,".DP.vs.QD.coloredByQUAL.postDP500.png",sep=""),p2,width=7,height=5)

########################## After all filtering (GATKHF, QD < 2, clustering, etc.) ######################
p3 <- ggplot(jointDist_postAllFilters,aes(x=QUAL,y=DP,color=QD))+
  geom_point(alpha=0.6)+
  theme_bw()
p3
ggsave(paste(data.dir,scaffold,".DP.vs.QUAL.all_5_filteredSites.png",sep=""),p3,width=7,height=5)

######################## Experiments ####################
# only QD < 2
p3 <- ggplot(jointDist_postDPFilter[jointDist_postDPFilter$QD < 2,],aes(x=QUAL,y=DP,color=QD))+
  geom_point(alpha=0.6)+
  theme_bw()+
  ggtitle("QD < 2 sites only")
p3
# used to be QD < 5; talked to Kirk and making it < 8
# exclude QD < 8 
p5a <- ggplot(jointDist_postAllFilters[jointDist_postAllFilters$QD >= 8,],aes(x=QUAL,y=DP,color=QD))+
  geom_point(alpha=0.75)+
  theme_bw()+
  ggtitle("QD > 8 sites only, post all filters")
p5a
ggsave(paste(data.dir,scaffold,".DP.vs.QUAL.all_5_filteredSites.noQDlt8.png",sep=""),p5a,width=7,height=5)
# exclude QD < 5 and show in red
p5b <- ggplot(jointDist_postAllFilters[jointDist_postAllFilters$QD >= 8,],aes(x=QUAL,y=DP,color=QD))+
  geom_point(alpha=0.75)+
  geom_point(data=jointDist_postAllFilters[jointDist_postAllFilters$QD < 8,],aes(x=QUAL,y=DP),color="red",alpha=0.2)+
  theme_bw()+
  ggtitle("QD > 8 sites comparison")
p5b
ggsave(paste(data.dir,scaffold,".DP.vs.QUAL.all_5_filteredSites.QDlt8inred.png",sep=""),p5b,width=7,height=5)


