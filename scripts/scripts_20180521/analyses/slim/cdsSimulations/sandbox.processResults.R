#playing with plotting
require(ggplot2)
#require(ggforce)
#require(concaveman)
popsModelsRundates='COM/1D.3Epoch.1.5Mb.cds/20190404 AK/1D.2Epoch.1.5Mb.cds/20190404 AL/1D.2Epoch.1.5Mb.cds/20190404 CA/1D.2Epoch.1.5Mb.cds/20190404 KUR/1D.2Epoch.1.5Mb.cds/20190404 ' # maybe? -- this is kind of awkward, maybe have to deal with diff populations differently?
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/concattedSummaries/"
rep=1
state="Pre"
preContract <- read.table(paste(data.dir,"AK/1D.2Epoch.1.5Mb.cds/20190404/rep.",rep,".slim.output.",state,"Contraction.allConcatted.summary.txt.gz",sep=""),header=T,sep=",",stringsAsFactors = F)

state="Post"
postContract <- read.table(paste(data.dir,"AK/1D.2Epoch.1.5Mb.cds/20190404/rep.",rep,".slim.output.",state,"Contraction.allConcatted.summary.txt.gz",sep=""),header=T,sep=",",stringsAsFactors = F)


#postContract$popsizeDIP <- 30 # temp fix <-- temporary!!! fix this asap!


head(postContract)
postContract$label <- "post-contraction"
preContract$label <- "equilibrium"


# get frequency:
postContract$frequency <- (postContract$p1numhet+2*postContract$p1numhom) / (2*postContract$popsizeDIP)
preContract$frequency <- (preContract$p1numhet+2*preContract$p1numhom) / (2*preContract$popsizeDIP)

### combine:

all <- rbind(postContract,preContract)


ggplot(all,aes(x=s,fill=label))+
  geom_histogram(bins=3,position="dodge")+
  scale_y_log10()

ggplot(all,aes(x=p1numhet,fill=label))+
  geom_density()

all$desc <- c("long-term pop size of 4000","after contraction to 30 inds for 4 gen")[as.factor(all$label)]
### woah crazy discovery
# doing c("xx","yy")[as.factor(all$label)] will make a vector of the same length of all$label and have it change from xx to yy when the label changes . what happens if there's a mismatch in lens?
# this is crazy. how did I never know this? what's it actually doing?

ggplot(all,aes(x=s,y=frequency,fill=label,color=label))+
  geom_point(alpha=0.5)+
  theme_bw()+
  geom_mark_hull(aes(fill=label,label=label,description=desc))

all$h <- 0
# calculate genetic load
# sum of 2hspiqi * siqi
# qi is frequency
# should s be made positive?
postContract$h <- 0
postContract$Scomponent <- 2*postContract$h*-postContract$s*postContract$frequency*(1-postContract$frequency) + (-postContract$s * (postContract$frequency^2))
preContract$h <- 0
preContract$Scomponent <- 2*preContract$h*-preContract$s*preContract$frequency*(1-preContract$frequency) + -preContract$s * (preContract$frequency^2)

# I think s shouldn't be negative in these calculations perhaps?
sum(preContract$Scomponent)
sum(postContract$Scomponent)
Wbar_pre=exp(-sum(preContract$Scomponent))
Wbar_post=exp(-sum(postContract$Scomponent))
Wbar_pre
Wbar_post
Load_pre = 1 - Wbar_pre
Load_post = 1 - Wbar_post
Load_pre
Load_post
