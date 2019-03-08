############ This calculates the SFS from fastsimcoal format
# need to figure out if input is folded or not
# look at other scripts to get file list by population 

data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/fastsimcoal2-plusMonomorphic/"
populations=c("AK","CA","AL","BER","MED","KUR")
########### THIS ONLY WORKS WITH FOLDED SFS (MAF) ###################
for(pop in populations){
  input <- read.table(paste(data.dir,pop,"_MAFpop0.obs",sep=""),skip=1,header=T,stringsAsFactors = F)# the SFS given by easy SFS (which is folded, but has zeroes in the latter half bins)
  sfs <- data.frame(t(input))
  colnames(sfs) <- c("count")
  # calculate pi: from sfs in fastsimcoal format (or dadi format?) or R format
  # input: frequency 
  # then make a function
  sfs$frequency <- as.numeric(lapply(strsplit(as.character(rownames(sfs)),"_"),"[",2))
  ###### SFS is to be in format count and frequency
  # get rid of monomorphic sites and make sure it's folded
  # skip 0 for pi calculation
  #n= # unfolded sample size, don't include zero bin
  # this has the unfolded bins as zeroes. but what if you don't have that?
  # need to standardize sfs format somehow. maybe go from dadi format since that's more usual? R format? folded/not folded? 
  # get n (the full unfolded sample size) with or without zero bins
  
  n=max(sfs$frequency) # get this from the d0_0 frequency; need a better way to get this
  # needs to be the full unfolded sample size. I want this to work even if the extra bins have been removed 
  # 
  totalSites=sum(sfs$count) # need to include monomorphic in the total sites for getting pi per site
  # set sum total as zero then go through sfs 
  sumTotal=0
  if(sum(sfs[sfs$frequency>n/2,]$count)==0){
    print("The SFS appears to be folded with zeroes in the final bins")
    for(i in seq(1,n/2)){
      # greek letter eta (pg 16 of Wakeley) is the folded count
      eta = sfs[sfs$frequency==i,]$count
      # from the equation 
      sumTotal=sumTotal+(i*(n-i)*eta)
    }
    # then multiply sumTotal by 1/(n choose 2)
    pi_overall = (1/(choose(n,2))) * sumTotal
    
    # then divide by total sites to get pi per site (note that monomorphic sites were also projected)
    pi_perSite = pi_overall/totalSites
  }
}


# write a calculation function that takes in n, and a vector of frequencies (excluding monomorphic)
# then do the processing upstream
# unfoldedHaploidSampleSize is n
calculatePiFromSFS <- function(unfoldedHaploidSampleSize,countsIncludingMonomorphic){
  n=as.numeric(unfoldedHaploidSampleSize)
  totalSites=sum(countsIncludingMonomorphic)
  sumTotal=0
  # need to skip over monomorphic so add 1 to each thing in the sequence
  for(i in seq(1+1,n/2+1)){
    # greek letter eta (pg 16 of Wakeley) is the folded count
    eta = countsIncludingMonomorphic[i]
    # from the equation 
    sumTotal=sumTotal+(i*(n-i)*eta)
  }
    # then multiply sumTotal by 1/(n choose 2)
    pi_overall = (1/(choose(n,2))) * sumTotal
    
    # then divide by total sites to get pi per site (note that monomorphic sites were also projected)
    pi_perSite = pi_overall/totalSites
    return(pi_perSite)
}

calculatePiFromSFS(20,sfs$count)
