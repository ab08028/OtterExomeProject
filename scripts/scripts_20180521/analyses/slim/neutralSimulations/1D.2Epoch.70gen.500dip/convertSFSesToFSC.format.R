######## convert SFS to FSC format ##########

require(ggplot2)
require(gridExtra)
####### sandbox trying to plot dadi results 
slimmodel="1D.2Epoch.70gen.500dip" # model you simulated under
slimdate=20190313 # date you ran slim
pop="generic" # what you named population in slim simulation
totalSitesSimulated=6000000 # 6Mb; so monomorphic will be 6000000 - total SFS sites 
# need to show a couple things
# 1) whether 1Epoch or 2Epoch fits better for the same replicate
# 2) if parameters converge for 2Epoch
# 3) what the range of those parameter estimates are 
dadimodel="1D.2Epoch" # model you inferred in dadi
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/neutralSimulations/",slimmodel,"/",slimdate,"/allSFSes/",sep = "")
data.dir

######## fold SFS function ######
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

############ generate all the folded fsc format sfses, with monomorphic sites ##########
### need to add in total sites simulated
# get R format SFSes (need to fold as well)
for(rep in seq(1,11)){
  results <- read.table(paste(data.dir,pop,".rep.",rep,".",slimmodel,".slim.output.unfolded.sfs.R.format.txt",sep=""),sep = "\t",header=T)
  # fold SFS:
  results = foldSFS(results)
  # select best LL
  # format for fsc
  # add to monomorphic bin
  totalSFSSites=sum(results$count)
  remainder=totalSitesSimulated - totalSFSSites 
  results[results$frequency==0,]$count = results[results$frequency==0,]$count + remainder # add remainder to monomorphic bin
  # add fsc ID
  results$fscID <- paste("d0_",as.character(results$frequency),sep="")
  if(sum(results$count) == totalSitesSimulated){
    outfile=paste(data.dir,pop,".simulation.rep.",rep,".",slimmodel,"_MAFpop0.obs",sep="")
    sink(outfile)
    cat("1 observations\n")
    cat(results$fscID,sep="\t")
    cat("\n")
    cat(results$count,sep="\t")
    cat("\n")
    sink()
  }
  else{
    print("Something is wrong with your counts")
    break
  }
}
