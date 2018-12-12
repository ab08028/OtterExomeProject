########### Average FST between populations #############
genotypeDate="20181119"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/fst/",genotypeDate,"/sandbox",sep="")
## CA vs. BAJ ##
fileList = list.files(data.dir,pattern=".gz",full.names = T)
fileList

dat = read.table(fileList[1],header=T)
dim(dat) # 1587850 what are all these? are these sites iwth no-calls or something? figure that out
# get rid of NaNs :
dat <- na.omit(dat)
dim(dat)

avgFst <- mean(dat$WEIR_AND_COCKERHAM_FST)
avgFst

## CA vs. AK ##

dat2 = read.table(fileList[2],header = T)
dim(dat2) # 1587850 what are all these? are these sites iwth no-calls or something? figure that out
# get rid of NaNs :
dat2 <- na.omit(dat2)
dim(dat2)

avgFst2 <- mean(dat2$WEIR_AND_COCKERHAM_FST)
avgFst2