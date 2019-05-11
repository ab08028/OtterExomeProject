###### need to fix naming with .bam in there 
###### where are the idx files? one level up -- make it all in per-individual
############# Get individual heterozyosity ###############
############# from the single individual sfs ##############
########### pi = sfs[2]/sum(sfs)
require(ggplot2)
angsdDate=20190508
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/SFS/",angsdDate,"/perIndividual/",sep="")
sfses <- list.files(data.dir,full.names = T)

allPis <- data.frame()
for(sfs in sfses){
  input <- read.table(sfs)
  piInfo= input[2]/sum(input)
  colnames(piInfo) <- "pi"
  # pull out info from the name
  piInfo$filename <- tail(unlist(strsplit(sfs,"/")),n=1)
  allPis = rbind(allPis,piInfo)
}


# want to label transversions
allPis$transv <- "Transitions+Transversions"
allPis[grep("TransversionsOnly",allPis$filename),]$transv <- "TransversionsOnly"
allPis$downsamp <- "Original"
allPis[grep("downsamp",allPis$filename),]$downsamp <- "Downsampled"
allPis$sampleID <- unlist(lapply(strsplit(allPis$filename,"\\."),"[",1))
allPis$reference <- NA
allPis[grep("Mustela",allPis$filename),]$reference <- "Ferret"
allPis[grep("sea_otter",allPis$filename),]$reference <- "Sea otter"

ggplot(allPis,aes(x=as.factor(sampleID),y=pi,fill=downsamp))+
  geom_bar(stat="identity",position="dodge",alpha=0.4)+
  coord_flip()+
  facet_wrap(reference~transv,scales="free")

########### SOMETHING HAS GONE QUITE WRONG ##############
#OKAY THESE ARE ALL 55. !!!!! boo.
# need to try script again
# ffff
# also want to plot the little SFSes
