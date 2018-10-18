########## Want to calculate pi per site for each population
########## and watterson's theta for each population (so need up to date sample size)
# so I need the total snp and the pi files
# for pi, if it was calculated from the callable neutral site regions for the population
# note that the coords are 0-based and half open so you can just subtract end - start to get total sites in the region
# so you sum up pi across all regions and you sum up all callable sites and divide sum(pi)/sum(callable sites)
# note that if you used an overall set of neutral regions, rather than the specifically callable sites WITHIN neutral regions, then the region start and end won't give you the total of callable sites <-- caution!!
# for watterson's theta, you need the total snps across all regions divided by the total callable sites and the harmonic number for the number of total chromosomes (2N = n ). 
require(ggplot2)
require(dplyr)
require(gridExtra)# Sample sizes for per-population files (removed admixed, relatives, etc.) this is the same that went into SFS:
pops=c("CA","AK","AL","COM","KUR")  # your populations
rundate=20180806 # genotype calling date
sfsdate=20181009 # date sfses were made
prefix="all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled"
############### set up your colors -- keep this consistent across all plots ######
require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=length(pops),name = "Dark2")
colors=list(CA=colorPal[1],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors

################### set up dir structure  ####################
todaysdate=format(Sys.Date(),format="%Y%m%d") # date you make plots

data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PI_THETA/",rundate,"/",sep="")
sfs.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",rundate,"/neutralSFS/",sep="")
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/PI_THETA/",rundate,"/",sep="")
analysis.dir=
dir.create(plot.dir,recursive = T,showWarnings = F)

################ Get total neutral callable sites from table ############
totalNeut <- read.table(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TotalCallableNeutralSites/",rundate,"/summary.neutralCallableSites.perPop.txt",sep=""),header = T)
################ Get Sample Size from SFS ###############
########### THE SFS MUST BE FOLDED FOR THIS TO WORK! 
# get sample sizes from the SFS; you calculate pi and theta from the same pop-specific vcf files
allSFS_list=list()
allSFS_length_list=list()
# this is a nice for-loop; it generates and writes out a plot for each population
for(i in (1:length(pops))){
  pop=pops[i]
  sfs <- read.table(paste(sfs.dir,pop,"_",prefix,".filtered.sfs.",sfsdate,".out",sep=""),header=T,stringsAsFactors = F) 
  sfs$population <- pop
  popcolor=colors[pop]
  allSFS_list[[i]] <- sfs # add the SFS to your list
}
allSFS <- bind_rows(allSFS_list)

ss <- allSFS %>%
  group_by(population) %>%
  tally()
colnames(ss) <- c("population","numDiploidIndividuals")
#################### calculate pi per site #####################
# based on output of tanya's script
# which gives you pi (not pi per site) for each region of called NEUTRAL sites
#### !!! Note: the callable sites calculation ONLY WORKS if the region file you gave Tanya's script was the neutral CALLED sites for each population, rather than some abstract set of regions where sites may or may not be called. In my analysis, I use the per-population called neutral sites, so it is okay to calculated callable sites from the output of Tanya's script.!!! ####
#######***** SOMETHING WEIRD WITH PI SCRIPT, I THINK BECAUSE OF LACK OF CHROMOSOME DESIGNATION ***** deal with this! Aha -- figured out why regions differ; without chromosomes, there is definitely some condensing of regions that have exact same start/stop points. Need to reclaculate pi!!! #######
allPi <- data.frame()
for(i in (1:length(pops))){
  pop=pops[i]
  pi <- read.table(paste(data.dir,pop,"_",prefix,".pi.perPopCallableSiteNeutralRegions.out",sep=""),header=F,stringsAsFactors = F)
  colnames(pi) <- c("region_start0based","region_end_halfOpen","pi_per_region")
  pi$siteCount <- pi$region_end_halfOpen - pi$region_start0based # note that this subtraction works because the coordinates are zero based and half open
  piSumAllRegions=sum(pi$pi_per_region)
  #callableSitesTotal=sum(pi$siteCount) # this only works because I used called sites in neutral regions per pop as my regions; too unreliable. using a special table instead:
  callableSitesTotal <- totalNeut[totalNeut$pop==pop,]$totalCalledNeutralSites
  piPerSite=(piSumAllRegions/callableSitesTotal)
  allPi <- rbind.data.frame(allPi,(cbind.data.frame(pop,piPerSite,piSumAllRegions,callableSitesTotal)))
  # add the pi per site to your dataframe
}
write.table(allPi,paste(data.dir,"allPops.piPerSitePerPopulation.",todaysdate,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
# temp plot (make a separate real plotting script)
piPlot <- ggplot(allPi,aes(x=pop,y=piPerSite,fill=pop))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=unlist(colors))+
  theme_bw()+
  theme(legend.title = element_blank())+
  xlab("Population")+
  ggtitle("Pi Per Neutral Site")
ggsave(paste(plot.dir,"allPops.piPerSitePerPopulation.",todaysdate,".pdf",sep=""),piPlot,device = "pdf",height=5,width=7)
#################### calculate watterson's theta #####################
# note harmonic number is sum from 1 to n of 1/n. n in this case is the number of choromosomes, not individuals, so it's 2*N if N is number of diploid inds.
# define a function to calc harm number:

#### HARMONIC NUMBER

# An = SUM(i=1-->n-1)(1/i)
# Caclulate harmonic number for 2n - 1 : 
# this function will calculate the harmonic number for a 
harmonicNumber_2nMinus1 <- function(num_indv){
  n_dip=num_indv # the numver of diploid individuals
  n_hap=2*n_dip # number of haploids
  n_harm=n_hap - 1 # number of haps minus 1 for harmonic number
  harm_number=0
  for(i in 1:(n_hap-1)){
    harm_number <- harm_number + 1/i
  }
  return(harm_number)
}

# example for 26 individuals:
#harmonicNumber_2nMinus1(26) # is the harmonic # for 2*26-1 = 51 -->  4.518813 # checked it here for 51: https://www.math.utah.edu/~carlson/teaching/calculus/harmonic.html

#### WATTERSON'S THETA (uses harmonic number function)
# this function will calculate watterson's theta using the number of diploid inds, 
# number of SNPs, and the number of callable (HQ) sites used for the calculation of # of SNPs. 
wattersons_theta <- function(num_indv,numSNPs,callablesites){
  harm_num <- harmonicNumber_2nMinus1(num_indv)
  theta <- numSNPs / harm_num
  theta_divbyL <- theta/callablesites
  return(theta_divbyL)
}


allTheta <- data.frame()
for(i in (1:length(pops))){
  pop=pops[i]
  # isntead of using Tanya's "Total SNPs" whcih seems to be having problems,
  #just get total snps from the sfs!
  #totalSNPs <- read.table(paste(data.dir,pop,"_",prefix,".totalSNPs.perPopCallableSiteNeutralRegions.out",sep=""),header=F,stringsAsFactors = F)
  #colnames(totalSNPs) <- c("region_start0based","region_end_halfOpen","SNPCount")
  #totalSNPs$siteCount <- totalSNPs$region_end_halfOpen - totalSNPs$region_start0based # note that this subtraction works because the coordinates are zero based and half open
  totalSNPsSum=sum(allSFS[allSFS$population==pop,]$num_variants)
  # get from Pi thing for now, while Total SNPs is being wonky:
  #callableSitesTotal=sum(totalSNPs$siteCount) # this only works because I used called sites in neutral regions per pop as my regions
  callableSitesTotal=allPi[allPi$pop==pop,]$callableSitesTotal
  # sample size in chromosomes for this population (which is 2x number of individuals, calculated earlier in the script)
  SampleSize=ss[ss$population==pop,]$numDiploidIndividuals # this is number of individuals NOT the number of chromosomes, because my harmonic number script multiplies it by 2 to get the number of chromosomes, so you don't need to do that yourself.
  # this uses my watterson's theta funciton; note that even though you input number of individuals, within the function it converts that to 2*numInd which is num of chromosomes.
  ThetaW <- wattersons_theta(num_indv = SampleSize,numSNPs = totalSNPsSum,callablesites = callableSitesTotal)
  allTheta <- rbind.data.frame(allTheta,(cbind.data.frame(pop,ThetaW,SampleSize,totalSNPsSum,callableSitesTotal))) # add the pi per site to your dataframe
}

thetaWPlot <- ggplot(allTheta,aes(x=pop,y=ThetaW,fill=pop))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=unlist(colors))+
  theme_bw()+
  theme(legend.title = element_blank())+
  xlab("Population")+
  ggtitle("Watterson's Theta Per Neutral Site")
thetaWPlot
write.table(allTheta,paste(data.dir,"allPops.WattersonsThetaPerSitePerPopulation.",todaysdate,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
ggsave(paste(plot.dir,"allPops.WattersonsTheta.PerSitePerPopulation.",todaysdate,".pdf",sep=""),thetaWPlot,device = "pdf",height=5,width=7)


############ Write out population information (work in progress) ########## 
for(pop in pops){
  sink(paste(data.dir,pop,".populationInformation.pi.theta.callableSites.txt",sep=""))
  cat(paste(todaysdate,pop,"\n## This is information about each population based on Tanya's popGen tools script that calcualtes pi across supplied regions.\nThe regions I used were putatively neutral regions intersected with the called sites from the filtered population-specific VCF file.\n##Therefore, the total called sites can be calculated from the regions (this isn't the case if you supply generic genomic regions in which case the sites may or may not be called inside them)\n\n"))
  cat(paste("Population:",pop,"\n"))
  cat(paste("Sample Size (based on SFS):",ss[ss$population==pop,]$numDiploidIndividuals,"\n"))
  cat(paste("Chromosome Count (2x sample size):",2*ss[ss$population==pop,]$numDiploidIndividuals,"\n"))
  cat(paste("Callable Neutral Sites:",allPi[allPi$pop==pop,]$callableSitesTotal,"\n"))
  cat(paste("pi (per site; neutral regions only):",round(allPi[allPi$pop==pop,]$piPerSite,5),"\n"))
  cat(paste("Watterson's Theta:",round(allTheta[allTheta$pop==pop,]$ThetaW,5),"\n"))
  sink()
}
