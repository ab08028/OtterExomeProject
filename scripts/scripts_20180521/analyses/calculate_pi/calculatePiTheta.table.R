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
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/PI_THETA/",rundate,sep="")
analysis.dir=
dir.create(plot.dir,recursive = T,showWarnings = F)
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
  tally() %>% 
  mutate(chromCount=2*n) # multiply by 

#################### calculate pi per site #####################
# based on output of tanya's script
# which gives you pi (not pi per site) for each region of called NEUTRAL sites


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

# for each population, carry out wattersons_theta calculation (maybe put all together and use dplyr?)
############ Write out population information (work in progress) ########## 
pop="CA"
sink(paste(data.dir,pop,".populationInformation.pi.theta.callableSites.txt",sep=""))
cat(paste("Population:",pop))
cat(paste("Sample Size (based on SFS):",ss[ss$population==pop,]$n))
cat(paste("Chromosome Count (2x sample size):",ss[ss$population==pop,]$chromCount))
cat(paste("Callable Neutral Sites:",callneut))
cat(paste("pi per site:",pi))
sink()