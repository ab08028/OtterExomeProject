########## Want to plot dadi and fsc2 results run on neutral data ###
######### dadi grid search ########
# I think I may skip grid search even though I ran it. Yeah I don't really need it.

################### Plot parameters ############
parametersInferredFromSims <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/SummaryOfDadiFSC_RunOnSimData/Summarize.Dadi.FSC.RunOnNeutralSimulations.ForR.2Epoch.txt",header=T)

simulatedTrueParams_melt=melt(data.frame(Nanc_FromTheta_scaled_dip=3500,nu_scaled_dip=200,T_scaled_gen=35))
simulatedTrueParams_melt$ParameterLabel <- NA
simulatedTrueParams_melt[simulatedTrueParams_melt$variable=="Nanc_FromTheta_scaled_dip",]$ParameterLabel <- "Nanc (diploids)"
simulatedTrueParams_melt[simulatedTrueParams_melt$variable=="nu_scaled_dip",]$ParameterLabel <- "Ncur (diploids)"
simulatedTrueParams_melt[simulatedTrueParams_melt$variable=="T_scaled_gen",]$ParameterLabel <- "T (gen.)"

require(reshape2)
require(ggplot2)
parametersInferredFromSims_melt <- melt(parametersInferredFromSims,id.vars = c("Model","Replicate","Label","Method"))
# make nice names for parameters:
parametersInferredFromSims_melt$ParameterLabel <- NA
parametersInferredFromSims_melt[parametersInferredFromSims_melt$variable=="Nanc_FromTheta_scaled_dip",]$ParameterLabel <- "Nanc (diploids)"
parametersInferredFromSims_melt[parametersInferredFromSims_melt$variable=="nu_scaled_dip",]$ParameterLabel <- "Ncur (diploids)"
parametersInferredFromSims_melt[parametersInferredFromSims_melt$variable=="T_scaled_gen",]$ParameterLabel <- "T (gen.)"

parametersInferredFromSims_melt$DataLabel2 <- NA
parametersInferredFromSims_melt[parametersInferredFromSims_melt$Label=="Empirical_CA",]$DataLabel2 <- "MLE based on empirical CA data"
parametersInferredFromSims_melt[parametersInferredFromSims_melt$Label=="Simulated",]$DataLabel2 <- "MLE based on simulated data"

plot1 <- ggplot(parametersInferredFromSims_melt,aes(x=Method,y=value,color=DataLabel2,shape=DataLabel2))+
  geom_point(size=4)+
  geom_hline(data=simulatedTrueParams_melt,aes(yintercept =value),linetype=2,color="purple")+
  facet_wrap(~ParameterLabel,scales="free_y")+
  theme_bw()+
  theme(legend.title = element_blank())+
  xlab("")+
  scale_shape_manual(values=c(8,1))+
  ylab("Parameter Estimate")+
  ggtitle("Comparing dadi and fastsimcoal inference\non simulated and empirical data\n(10 replicates)\nSimulated values marked by purple lines")

plot1
ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/SummaryOfDadiFSC_RunOnSimData/2Epoch.Dadi.Fsc.NeutralSimulations.CompareToEmpirical.pdf",plot1,device="pdf",height=5,width=9)

######################## LRT between simulated 1 and 2 epoch models #############
########### do likelihood ratio test #########
LLs <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/SummaryOfDadiFSC_RunOnSimData/Summarize.Dadi.FSC.RunOnNeutralSimulations.ForR.LLs.1Epoch.2Epoch.txt",header=T)
######## Now do Lhood ratio test for more complex to more simple models ##########
LRTResults <- data.frame()
############### 1Epoch vs 2 Epoch ############
model0="1D.1Epoch"
k0=1 # parameters for model (Nanc)
model1="1D.2Epoch"
k1=3 # parameters of 2 epoch model (Nanc,nu,T)

LLs$df <- k1-k0
LLs$LRT <- -2*(LLs$LL_1Epoch- LLs$LL_2Epoch)
LLs$pval <- pchisq(LLs$LRT,LLs$df,lower.tail = F)
alpha=0.05/20
LLs$alpha <- alpha
LLs$significantAtAlpha <- NA
LLs[LLs$pval<alpha,]$significantAtAlpha <- "Significant"
LLs[LLs$pval>=alpha,]$significantAtAlpha <- "Not Significant"

View(LLs)
write.table(LLs,"/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/SummaryOfDadiFSC_RunOnSimData/LLs.Dadi.Fastsimcoal.LikelihoodRatioTestComparison.txt",row.names=F,quote=F,sep="\t")


