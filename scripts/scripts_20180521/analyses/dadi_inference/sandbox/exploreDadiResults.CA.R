######### Dadi Results processing ###########
todaysdate=format(Sys.Date(),format="%Y%m%d") # date you make plots
pop="CA"
# want to talk to Jim Estes about possible models.
generationTime=6 # for now using 6 yr/ gen (Tinker says 6-7 is reasonable)
# get all output:
genotypeDate=20180806
dadiInferenceDate=
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/",genotypeDate,"/",sep="")
#models=c("1D.1Bottleneck","1D.2Bottlenecek","1D.2Epoch")
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/dadi_inference/",genotypeDate,"/",sep="")
dir.create(plot.dir,recursive = T)
#### Set up a little data.frame to keep track of best fit models ####
df = data.frame()
################ 2 Epoch ##################
model0="1D.2Epoch"
k0 = 2 # number of params
results0 <- read.table(paste(data.dir,"/",pop,"/",model0,"/dadi.inference.",model0,".all.output.concatted.txt",sep=""),sep="\t",header=T)
# calculate Nanc:
results0$Nanc <- results0$theta/(4*results0$mu*results0$L)
names(results0)
# scale results0 (eventually do in python)
results0$nu_scaled <- results0$nu * results0$Nanc
results0$T_scaled_gen <- results0$T * 2 * results0$Nanc
results0$T_scaled_yr <- results0$T * 2 * results0$Nanc * generationTime

results0$modelNumber <- 0 # internal to R
# get AIC score
# get best fit
model0_best = results0[results0$LL==-min(abs(results0$LL)),]
# get AIC
model0_best$AIC = 2*k0 - 2*model0_best$LL

# add model name and AIC to df (don't add params because will differ across models)
df = rbind(df,model0_best[,c("LL","AIC","modelFunction","rundate","runNumber","modelNumber")])
################ 1 bottleneck ##################
model1="1D.1Bottleneck"
k1 = 4
results1 <- read.table(paste(data.dir,"/",pop,"/",model1,"/dadi.inference.",model1,".all.output.concatted.txt",sep=""),sep="\t",header=T)
# calculate Nanc:
results1$Nanc <- results1$theta/(4*results1$mu*results1$L)

# scale results (eventually do in python)
results1$nuB_scaled <- results1$nuB * results1$Nanc
results1$nuF_scaled <- results1$nuF * results1$Nanc
results1$TB_scaled_gen <- results1$TB * 2 * results1$Nanc
results1$TF_scaled_gen <- results1$TF * 2 * results1$Nanc
results1$TB_scaled_yr <- results1$TB * 2 * results1$Nanc * generationTime
results1$TF_scaled_yr <- results1$TF * 2 * results1$Nanc * generationTime

results1$modelNumber <- 1 # internal to R

# get best fit
model1_best = results1[results1$LL==-min(abs(results1$LL)),]
# get AIC
model1_best$AIC = 2*k1 - 2*model1_best$LL
# add to tracking df:
df = rbind(df,model1_best[,c("LL","AIC","modelFunction","rundate","runNumber","modelNumber")])

################# 2 bottleneck ################
model2="1D.2Bottleneck"
k2 = 8
results2 <- read.table(paste(data.dir,"/",pop,"/",model2,"/dadi.inference.",model2,".all.output.concatted.txt",sep=""),sep="\t",header=T)
# calculate Nanc:
results2$Nanc <- results2$theta/(4*results2$mu*results2$L)

# scale results (eventually do in python)
results2$nuB1_scaled <- results2$nuB1 * results2$Nanc
results2$nuF1_scaled <- results2$nuF1 * results2$Nanc
results2$TB1_scaled_gen <- results2$TB1 * 2 * results2$Nanc
results2$TF1_scaled_gen <- results2$TF1 * 2 * results2$Nanc
results2$TB1_scaled_yr <- results2$TB1 * 2 * results2$Nanc * generationTime
results2$TF1_scaled_yr <- results2$TF1 * 2 * results2$Nanc * generationTime
results2$nuB2_scaled <- results2$nuB2 * results2$Nanc
results2$nuF2_scaled <- results2$nuF2 * results2$Nanc
results2$TB2_scaled_gen <- results2$TB2 * 2 * results2$Nanc
results2$TF2_scaled_gen <- results2$TF2 * 2 * results2$Nanc
results2$TB2_scaled_yr <- results2$TB2 * 2 * results2$Nanc * generationTime
results2$TF2_scaled_yr <- results2$TF2 * 2 * results2$Nanc * generationTime

results2$modelNumber <- 2 # internal to R

# get best fit
model2_best = results2[results2$LL==-min(abs(results2$LL)),]
# get AIC
model2_best$AIC = 2*k2 - 2*model2_best$LL

df = rbind(df,model2_best[,c("LL","AIC","modelFunction","rundate","runNumber","modelNumber")])

############## Pick best model #################


bestModel=df[df$AIC==min(df$AIC),]

# pull out those parameters: 
bestModelResults <- get(paste("results",bestModel$modelNumber,sep=""))
bestModelRunParams <- bestModelResults[bestModelResults$runNumber==bestModel$runNumber,]
bestModelName=bestModelRunParams$modelFunction

write.table(bestModelRunParams,paste(data.dir,"/",pop,"/",pop,".bestModelRunParams.AIC.",todaysdate,".txt",sep=""),quote=F,row.names=F)

# locate output files that correspond and copy them into a dir.:
# dadi.inference.1D.1Bottleneck.runNum.33.20181019.output
############# Plot Convergence ###################
# plot iterations in order of inc. ll and show how parameters change
#subset results:

results1_sub <- results1[,c("runNumber","Nanc","nuB_scaled","nuF_scaled","TB_scaled_gen","TF_scaled_gen","LL")]  
results1_sub_melt <- melt(results1_sub,id.vars = c("runNumber","LL"))
ggplot(results1_sub_melt,aes(x=runNumber,y=value,fill=LL,color=LL))+
  geom_point()+
  facet_wrap(~variable,scales="free")

convergePlot <- ggplot(bestModelResults,aes(x=reorder(runNumber,LL),y=LL))+
  geom_point()+
  ggtitle(paste(pop," ",bestModelName," LL Convergence",sep=""))
convergePlot
ggsave(paste(plot.dir,"/",pop,".",bestModelName,".LLConvergence.",todaysdate,".pdf",sep=""),convergePlot,height=5,width=7)