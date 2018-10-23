######### Dadi Results processing ###########
#populations=c("CA","AK","AL","COM","KUR") # I think make one script to explore each population separately
pop="CA"
# want to talk to Jim Estes about possible models.

# get all output:
genotypeDate=20180806
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/",genotypeDate,"/",sep="")
#models=c("1D.1Bottleneck","1D.2Bottlenecek","1D.2Epoch")
################ 2 Epoch ##################
model0="1D.2Epoch"
results0 <- read.table(paste(data.dir,"/",pop,"/",model0,"/dadi.inference.",model0,".all.output.concatted.txt",sep=""),sep="\t",header=T)
# calculate Nanc:
results0$Nanc <- results0$theta/(4*results0$mu*results0$L)
names(results0)
# scale results0 (eventually do in python)
results0$nu_scaled <- results0$nu * results0$Nanc
results0$T_scaled <- results0$T * 2 * results0$Nanc
################ 1 bottleneck ##################
model1="1D.1Bottleneck"
results1 <- read.table(paste(data.dir,"/",pop,"/",model1,"/dadi.inference.",model1,".all.output.concatted.txt",sep=""),sep="\t",header=T)
# calculate Nanc:
results1$Nanc <- results1$theta/(4*results1$mu*results1$L)

# scale results (eventually do in python)
results1$nuB_scaled <- results1$nuB * results1$Nanc
results1$nuF_scaled <- results1$nuF * results1$Nanc
results1$TB_scaled <- results1$TB * 2 * results1$Nanc
results1$TF_scaled <- results1$TF * 2 * results1$Nanc


################# 2 bottleneck ################

################ 1 bottleneck ##################
model2="1D.2Bottleneck"
results2 <- read.table(paste(data.dir,"/",pop,"/",model2,"/dadi.inference.",model2,".all.output.concatted.txt",sep=""),sep="\t",header=T)
# calculate Nanc:
results2$Nanc <- results2$theta/(4*results2$mu*results2$L)

# scale results (eventually do in python)
results2$nuB1_scaled <- results2$nuB1 * results2$Nanc
results2$nuF1_scaled <- results2$nuF1 * results2$Nanc
results2$TB1_scaled <- results2$TB1 * 2 * results2$Nanc
results2$TF1_scaled <- results2$TF1 * 2 * results2$Nanc
results2$nuB2_scaled <- results2$nuB2 * results2$Nanc
results2$nuF2_scaled <- results2$nuF2 * results2$Nanc
results2$TB2_scaled <- results2$TB2 * 2 * results2$Nanc
results2$TF2_scaled <- results2$TF2 * 2 * results2$Nanc

############# Convergence ###################
# plot iterations in order of inc. ll and show how parameters change
#subset results:
results0_sub <- results0[,c("runNumber","nu","T","LL")]  
results0_sub_melt <- melt(results0_sub,id.vars = c("runNumber","LL"))
ggplot(results0_sub_melt,aes(x=runNumber,y=value,fill=LL,color=LL))+
  geom_point()+
  facet_wrap(~variable,scales="free")

results1_sub <- results1[,c("runNumber","nuB","TB","LL")]  
results1_sub_melt <- melt(results1_sub,id.vars = c("runNumber","LL"))
ggplot(results1_sub_melt,aes(x=runNumber,y=value,fill=LL,color=LL))+
  geom_point()+
  facet_wrap(~variable,scales="free")

############# AIC COMPARISONS ##############
