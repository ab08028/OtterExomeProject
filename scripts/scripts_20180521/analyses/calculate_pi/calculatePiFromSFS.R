############ This calculates S, pi and Theta from SFS that is in fsc2 format from easy sfs
# it is a slightly odd format because it is in fact *FOLDED* but it has zeroes in half the bins 
# Want to get rid of those zeroes and convert to R format
######## be very careful not to use this script with an SFS in any other input format except from EASY SFS 
genotypeDate="20181119"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/fastsimcoal2-plusMonomorphic/",sep="")
out.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PI_THETA/empirical/",genotypeDate,"/",sep="")
dir.create(out.dir,showWarnings = T,recursive = T)
#populations=c("AK","CA","AL","BER","MED","KUR")
populations=c("AK","CA","AL","COM","KUR")
################# functions #################
############# convert sfs from fsc2 to R format ##############
convertToRFormat <- function(sfs_fscFormat){
  sfs=data.frame(t(sfs_fscFormat))
  colnames(sfs) <- c("count")
  # calculate pi: from sfs in fastsimcoal format (or dadi format?) or R format
  # input: frequency 
  # then make a function
  sfs$frequency <- as.numeric(lapply(strsplit(as.character(rownames(sfs)),"_"),"[",2))
  return(sfs)
}
###################### fold sfs function ############
foldSFS <- function(sfs){
  foldedSFS <- data.frame()
  ss=length(sfs$frequency) - 1 # this is the ss in chromosomes
  foldedBin=ss/2  # half as many as ss ; do seq from 0 --> foldedBins to include monomorphics
  # note you have to start seq at 0 to include monomorphic bin 
  for(i in seq(0,foldedBin)){
    # if it's not the center point (which doesn't get added together)
    # see wakeley coalescent eq 1.2
    if(i==ss-i){
      foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count)/2)) # if it's add midpoint (foldedLen, including 0 monomorphic bin), add together and divide by two (equivalent of not adding together)
    }
    else{
      foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count))) # if not at mid point, just add together like normal
    }
  }
  return(foldedSFS)
}
####################### calculate pi from folded SFS function #######################
# write a calculation function that takes in n, and a vector of frequencies (excluding monomorphic)
# then do the processing upstream
# unfoldedHaploidSampleSize is n
calculatePiFromSFS_empData <- function(unfoldedHaploidSampleSize,totalSimulatedSites,sfsFoldedExcludingMonomorphic){
  n=as.numeric(unfoldedHaploidSampleSize)
  #print(n)
  totalSites=totalSimulatedSites
  sumTotal=0
  # need to skip over monomorphic so add 1 to each thing in the sequence
  # make sure monomorphic is excluded: (should already be, but this makes sure)
  sfsExclMono = sfsFoldedExcludingMonomorphic[sfsFoldedExcludingMonomorphic$frequency!=0,]
  #print(sfsExclMono)
  for(i in seq(1,n/2)){
    # greek letter eta (pg 16 of Wakeley) is the folded count
    eta = sfsExclMono[i,2] # second column is count 
    #print(eta)
    # from the equation 
    sumTotal=sumTotal+(i*(n-i)*eta)
    #print(sumTotal)
  }
  # then multiply sumTotal by 1/(n choose 2)
  pi_overall = (1/(choose(n,2))) * sumTotal
  #print(pi_overall)
  # then divide by total sites to get pi per site (note that monomorphic sites were also projected)
  pi_perSite = pi_overall/totalSites
  return(pi_perSite)
}

####################### functions to calculate S and Wattersons Theta ############
# watterson's theta uses n-1 haploid individuals for harmonic number
# my harmonicNumber_2nMinus1 function takes in the number of *diploid* individuals
# and converts it to haploid 
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

####### read in SFS, convert to R format, 'fold' to get rid of the zero bins (even though it's already folded -- this is just eliminating the extra zero bins), and calculate stats ############ 
df <- data.frame(population=character(),pi=numeric(),S=numeric(),Wattersons_theta=numeric(),stringsAsFactors = F)
for(pop in populations){
  file=paste(data.dir,pop,"_MAFpop0.obs",sep="")
  popdf <- data.frame(population=character(),pi=numeric(),S=numeric(),Wattersons_theta=numeric(),stringsAsFactors = F) # specific to model 
  # specific to populatio
  input <- read.table(file,header=T,stringsAsFactors = F,skip=1) # input in "R format" by my simulations
  # convert to R format:
 sfs <-  convertToRFormat(input)
  # get rid of monomorphic sites and fold SFS 
  # fold SFS -- note it's actually already folded, but there are zeroes in the unfolded bins
 # so we are 'folding' it to eliminate those zeroes, but it's not actually changing counts
 # and it's not double-folding *** be careful here -- if those zero bins aren't present you might accidentally double-fold! 
  sfsFolded <- foldSFS(sfs)
  nDip=max(sfsFolded$frequency)
  nHap=2*nDip # get unfolded sample size of total chromosomes from FOLDED SFS
  totalSites=sum(sfsFolded$count) # must include monomorphic in this total
  # need to exclude monomorphic!
  sfsFoldedExclMono <- sfsFolded[sfsFolded$frequency!=0,]
  pi=calculatePiFromSFS_empData(nHap,totalSites,sfsFoldedExclMono) # 
  # divide n by two to get diploid individuals (for input into function; gets converted to haploid within function) # no
  S=sum(sfsFoldedExclMono$count) #excludes monomorphic
  theta= wattersons_theta(nDip,S,totalSites) # note this uses Ndip in the function because the function itself converts nDip to nHap (this is a bit awkward but comes of merging functions across scripts)
  # add entry to data frame
  popdf[1,"population"] <- pop # need to include model
  popdf[1,"pi"] <- pi # and value of pi
  popdf[1,"S"] <- S
  popdf[1,"Wattersons_theta"] <- theta
  df <- rbind(df,popdf)
# combine all popdfs:
}

# write out the df in the approrpriate indir
write.table(df,paste(out.dir,"pi.S.Theta.CalculatedFromSFSes.allpops.txt",sep=""),quote = F,row.names=F)

