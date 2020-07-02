require(ggplot2)
require(reshape2)
######### Plot pi, theta and S for empirical and simulated data
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PI_THETA/"

out.dir=paste(data.dir,"/plots/",sep="")
dir.create(out.dir,recursive = T)
# simulated models: 
models=c("AK.1D.2Epoch.35Gen.250Inds/20200129/","CA.1D.2Epoch.25Gen.100Inds/20200129/")
genotypeDate="20181119" # for empirical
# input of all populations
emp <- read.table(paste(data.dir,"/empirical/",genotypeDate,"/pi.S.Theta.CalculatedFromSFSes.allpops.txt",sep=""),header=T)

# for each model read in simulated data:
simdata <- data.frame()
for(model in models){
  input <- read.table(paste(data.dir,"/simulated/",model,"/pi.S.Theta.CalculatedFromSFSes.allreps.PrePostContraction.txt",sep=""),header=T)
  # combine each model together:
  simdata <- rbind(simdata,input)
}

# make colnames match:
simdata$population <- unlist(lapply(strsplit(as.character(simdata$model),"\\."),"[",1))
simdata$label <- simdata$state
simdata$color <- "simulated"
emp$model <- "empirical"
emp$label <- "empirical"
emp$color <- "empirical"
emp$rep <- "empirical"
emp$state <- "postContraction"
# restric temp to just AK and CA:
emp <- emp[emp$population %in% c("AK","CA"),]
emp_melt <- melt(emp)

allData <- rbind(simdata,emp)
allData_melt <- melt(allData)
#################### Plot empirical and simulated data together ##########
require(ggplot2)
# exclude S, doesnt make much sense since number of called sites differ
# get an average across replicates 
simdata_avgAcrossReps <- simdata %>%
  group_by(model,state,population,label,color) %>%
  summarise(pi=mean(pi),Wattersons_theta=mean(Wattersons_theta),sdPi=sd(pi),sdTheta=sd(Wattersons_theta))

simdata_avgAcrossReps_melt <- melt(simdata_avgAcrossReps)
simdata_avgAcrossReps_melt$label <- factor(simdata_avgAcrossReps_melt$label,levels=c("preContraction","postContraction"))
p1 <- ggplot(simdata_avgAcrossReps_melt[simdata_avgAcrossReps_melt$variable %in% c("pi","Wattersons_theta"),],aes(x=label,y=value,fill=color))+
  geom_col()+
  geom_col(data=emp_melt[emp_melt$variable %in% c("pi","Wattersons_theta"),],aes(x=label,y=value))+
  facet_wrap(population~variable,scales="free")+
  theme_bw()
p1
ggsave(paste(out.dir,"NeutralPiFromSimulated.vsEmpirical.AK.CA.pdf",sep=""),p1,height=5,width=9)
