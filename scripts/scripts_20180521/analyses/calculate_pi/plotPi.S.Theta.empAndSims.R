require(ggplot2)
require(reshape2)
######### Plot pi, theta and S for empirical and simulated data
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PI_THETA"
out.dir=paste(data.dir,"/plots/",sep="")
dir.create(out.dir,recursive = T)
# simulated models: 
models=c("1D.1Epoch/20190213/","1D.2Epoch/20190125/","1D.2Epoch.70gen.500dip/20190313/","1D.2Epoch.4gen/20190128/") # taking out ,"1D.2Epoch.30gen/20190227/" because it was not based on any inference and is really off -- dont use that model
genotypeDate="20181119" # for empirical
# input of all populations
emp <- read.table(paste(data.dir,"/empirical/",genotypeDate,"/pi.S.Theta.CalculatedFromSFSes.allpops.txt",sep=""),header=T)

# for each model read in simulated data:
simdata <- data.frame()
for(model in models){
  input <- read.table(paste(data.dir,"/simulated/",model,"pi.S.Theta.CalculatedFromSFSes.allreps.txt",sep=""),header=T)
  # combine each model together:
  simdata <- rbind(simdata,input)
}

# make colnames match:
simdata$population <- "generic"
simdata$label <- simdata$model
simdata$color <- "simulated"
emp$model <- "empirical"
emp$label <- emp$population
emp$color <- "empirical"
allData <- rbind(simdata,emp)
allData_melt <- melt(allData)
#################### Plot empirical and simulated data together ##########
require(ggplot2)
# exclude S, doesnt make much sense since number of called sites differ
p1 <- ggplot(allData_melt[allData_melt$variable!="S",],aes(x=label,y=value,color=color))+
  geom_boxplot()+
  facet_wrap(~variable)+
  theme_bw()+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  ggtitle("Note: 1D.2Epoch.4gen is based on AK parameters and has slightly higher Nanc than 1D.1Epoch\n1D.2Epoch model has T=10, seems too long ")
p1
ggsave(paste(out.dir,"/pi.theta.simulated.empirical.pdf",sep=""),p1,device="pdf",width=11,height=7)
