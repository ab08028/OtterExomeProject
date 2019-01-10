###################### Adding monomorphic sites to fastsimcoal SFSes ##########
# don't need to add to dadi because are masked anyway
# just need the total value of sites for dadi (L) that is sum of SFS + monomorphic. going to output that 
# make this work on Hoffman as part of easySFS script.
library("optparse")
option_list = list(
  make_option(c("-dataDir", "--dataDir"), type="character", default=NULL, 
              help="path to your data directory (projection output of easy SFS)", metavar="character"),
  make_option(c("-popFile", "--popFile"), type="character", default=NULL, 
              help="path to your data directory (projection output of easy SFS)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

data.dir=opt$dataDir
pop=opt$popFile
popFile=read.table(pop,header=F)

#genotypeDate=20181119
#projectionDate=20181212 # date projection was carried out
#data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/projection-",projectionDate,"/",sep="")
#popFile=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/easySFSPopMapFiles/samplesPop.Headers.forEasySFS.3.20181119.txt",header=F) # used in easySFS; update if used a different one. The order of populations in this script is the order that pops are assigned numbers in easy sfs (that's a bit hacky of easySFs). So get order of pops from this file for 0/1 (double check manually)

fsc.format.dir=paste(data.dir,"/fastsimcoal2/",sep="") # where easysfs output is
new.fsc.format.dir=paste(data.dir,"/fastsimcoal2-plusMonomorphic/",sep="") # where you'll put the new sfses that have monomorphic sites added in
dir.create(new.fsc.format.dir,showWarnings = F) 
dadi.format.dir = paste(data.dir,"/dadi/",sep="")
new.dadi.format.dir=paste(data.dir,"/dadi-plusMonomorphic/",sep="") # where you'll put the new sfses that have monomorphic sites added in
dir.create(new.dadi.format.dir,showWarnings = F) 
colnames(popFile) <- c("sample","population")
popOrder <- as.character(unique(popFile$population)) # this should be CA,AK,AL,COM,KUR  for my project

####################### read in file of monomorphic site counts that are monomorphic in all individuals and are called in at least [projection value] for each population ######
# note that sites that are variable in the whole sample, but monomorphic in one population, are already included in the zero bin of the sfs. # 

monomorph <- read.table(paste(data.dir,"/countsOfMonomorphicPassingProjectionThresholds.perPop.txt",sep=""),header = T)
monomorph2D <- read.table(paste(data.dir,"/countsOfMonomorphicPassingProjectionThresholds.perPair.txt",sep=""),header = T)


################################## fastsimcoal format ################################
#### go through SFSes and write them out with the monomorphic added in: #########
for(pop in popOrder) {
  input <- list.files(fsc.format.dir,pattern=pop,full.names = T)
  sfs <- read.table(input,skip = 1,header = T) # skip first line: "1 observation"
  monoCount <- monomorph[monomorph$population==pop,]$HomREFcount
  sfs$d0_0 <- sfs$d0_0+monoCount
  sink(paste(new.fsc.format.dir,pop,"_MAFpop0.obs",sep=""))
  cat("1 observation\n")
  write.table(sfs,quote=F,row.names = F) # writes new sfs to the table
  sink()
}

################################## dadi format ########################################

for(pop in popOrder) {
  # pattern is pop-##.sfs
  input <- list.files(dadi.format.dir,pattern=paste(pop,"-[0-9]+.sfs",sep=""),full.names = T)
  #print(input)
  header <- readLines(input,n=1) # get first line of file
  sfs <- read.table(input,skip = 1,header = F) # skip first line: "13 folded "KUR""
  monoCount <- monomorph[monomorph$population==pop,]$HomREFcount
  # two rows; first is counts, second is mask; don't change the mask
  sfs[1,"V1"] <- sfs[1,"V1"]+monoCount
  totalSites <- sum(sfs[1,])
  sink(paste(new.dadi.format.dir,pop,"-",as.character(length(sfs)),".plusMonomorphic.sfs",sep=""))
  cat(header,"\n") # put the header back in
  write.table(sfs,quote=F,row.names = F,col.names = F) # writes new sfs to the table
  sink()
  sink(paste(new.dadi.format.dir,pop,"-",as.character(length(sfs)),".totalSiteCount.L.withMonomorphic.txt",sep=""))
  cat("pop\ttotalSites\n")
  cat(pop, totalSites,"\n")
  sink()
}

####################################### 2D SFS ###################################
# fsc files are numbered by pop0_1 where pop numbers are from my pop order (check this carefully)
# order should be CA,AK,AL,COM,KUR  for my project (0,1,2,3,4)
#combos = combn(seq(0,length(popOrder)-1),2)
for(i in seq(0,length(popOrder)-1)) {
  for(j in seq(0,length(popOrder)-1)){
  # note counts are zero based in filename (eg CA is 0), but 1 based in R
  pop1=popOrder[i+1] # i and j are zero based so add 1 to get right index
  pop2=popOrder[j+1]
  input <- list.files(fsc.format.dir,pattern=paste("jointMAFpop",i,"_",j,sep=""),full.names = T)
  if(length(input)==1){
    sfs <- read.table(input,skip = 1,header = T) # skip first line: "1 observation"
    monoCount <- monomorph2D[monomorph2D$population1==pop1 & monomorph2D$population2==pop2,]$HomREFcountPassingBothProjThresholds
  # row and column names reflect pop numbers; first number in file name is colnames, second is rownames. so row is first and should be j, and col is second and should be i
    sfs[paste("d",j,"_0",sep=""),paste("d",i,"_0",sep = "")] <- sfs[paste("d",j,"_0",sep=""),paste("d",i,"_0",sep = "")]+monoCount
    sink(paste(new.fsc.format.dir,"/neutral.",pop1,".",pop2,"_jointMAFpop",i,"_",j,".obs",sep=""))
    cat("1 observation\n")
    write.table(sfs,quote=F,row.names = T) # writes new sfs to the table
    sink()
  }
  }
}
