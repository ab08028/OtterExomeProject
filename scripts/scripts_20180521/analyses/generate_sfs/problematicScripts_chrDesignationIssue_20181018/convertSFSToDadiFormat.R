## Convert from Tanya SFS to dadi format (1D for now)

genotypeDate=20180806
sfs.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/neutralSFS/",sep="")
suffix="_all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled.filtered.sfs.20181009"
pop="CA"
# this is in Tanya's format:
# frequency	num_variants
# 1	489
# 2	231
# 3	248
# 4	278
# 5	240
sfs=read.table(paste(sfs.dir,pop,suffix,".out",sep=""),header=T)

# this assumes SFS is folded already, and will add in zeroes for the extra values that got folded 
Tanya2dadi_StartsFolded <- function(sfs,out.dir){
  todaysdate=format(Sys.Date(),format="%Y%m%d") # date you make plots
  dimSFS=(2*dim(sfs)[1])+1 # num haploids +1 - this is full unfolded dim.
  status="folded" # change if unfolded
  out.file=paste(out.dir,pop,suffix,".dadi.format.out",sep="")
  entries=sfs$num_variants 
  sink(out.file)
  # header:
  cat("# ",pop," folded SFS in dadi format. Generated: ",todaysdate,"\n",sep="")
  # give the dimensions, status (folde/unfolded), and pop name
  cat(dimSFS," ",status," \"",pop,"\"","\n",sep="")
  # have to add 0 for monomorphic (gets masked in dadi anyway)
  cat(0," ",sep="") # make it empty, with no \n so a newline doesn't start
  cat(entries,rep(0,dimSFS-dim(sfs)[1]-1))
  # put in mask (want first and last entry masked with "1", no others masked)
  cat("\n",1," ",sep="")
  cat(rep(0,dim(sfs)[1]),rep(1,dimSFS-dim(sfs)[1]-1),1,sep=" ")
  sink()
}

# do this for all populations *** question **** do we think Tanya SFS is affected by lack of chr designation too????