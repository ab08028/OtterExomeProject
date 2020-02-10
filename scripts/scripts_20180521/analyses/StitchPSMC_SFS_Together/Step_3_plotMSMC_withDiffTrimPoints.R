require(ggplot2)
require(scales)
require(RColorBrewer)
require(scales)
elutCol="#377EB8"
pbraCol="#4DAF4A"
elut2Col="#C77CFF"
################### southern sea otter ############
################ read in scaled results #############

sso = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/MSMC_DemographicInference_WholeGenome/output_20180209_50_250DPFilter/elut.MSMC.results.allScalings.10_14_20_MyaDivergence.txt",header=T)


# northern sea otter (AK/AL)
nso=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/MSMC_50_250DPFilter_simsInSppDIrs/output_20190215/elut2.MSMC.results.allScalings.10_14_20_MyaDivergence.txt",header=T)
# trim: 
################## set parameters ############
mu=8.64e-09

################ PLOTS ###################

################# MAIN TEXT FIGURE ::: PLOT Inv Coal Rate and generations (reasonable mu), WITHTRIMMING ##########
# Trimming information: Chose index 33 for elut, 31 for elut2, 20 for pbra:
elutTimeIndices=c(39,33,27,22) # sso; 39 is  untrimmed - just taking off last inf one
elut2TimeIndices=c(39,31,27,20) # nso; 39 is  untrimmed -- just taking off last inf interval nso

# because you use left generatins and are setting that time index you remove everything before
# the time index -1 time index. (based on python script)
# so new Na is:
############# NOTE index -1: (spent a lot of time tirekicking; this is right)

############ NOTE A CHANGE IN DEPICTING CUTOFFS: need to draw a line at (index -1) because you are totally removing that index in python
for(elutTimeIndex in elutTimeIndices){
  
# get weighted average Ne (weighted by duration of each interval)
weightedNeAvg <- sum(sso[sso$time_index<elutTimeIndex,]$Ne.reasonable*(sso[sso$time_index<elutTimeIndex,]$Right_generations.reasonable-sso[sso$time_index<elutTimeIndex,]$Left_generations.reasonable)/sum(((sso[sso$time_index<elutTimeIndex,]$Right_generations.reasonable-sso[sso$time_index<elutTimeIndex,]$Left_generations.reasonable))))
plot1 <- ggplot(sso,aes(x=Left_generations.reasonable,y=Ne.reasonable))+
  geom_step(stat="identity",size=1,color=elutCol)+
  theme_bw()+
  theme(legend.title = element_blank())+
  xlab("Generations Ago")+
  ylab("IICR (~Ne)")+
  scale_y_log10(labels=comma)+
  scale_x_log10(labels=comma)+
  theme(legend.position= c(0.6,0.85),legend.background = element_rect(fill="transparent"))+
  geom_vline(xintercept = sso[sso$time_index==elutTimeIndex-1,]$Left_generations.reasonable,linetype=2,color="grey20")+
  theme(legend.direction=("vertical"),legend.position=c(0.5,0.8),legend.background = element_rect(fill="transparent"),legend.text=element_text(size=14),legend.key.size = unit(1,"cm"))+
  annotate("text",label=paste("Nanc=",round(sso[sso$time_index==elutTimeIndex-1,]$Ne.reasonable),"\nOldest gen = ",round(sso[sso$time_index==elutTimeIndex-1,]$Left_generations.reasonable),"\nWeighted (by interval time) Ne =", round(weightedNeAvg),sep=""),x=5000,y=12000)
plot1
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/StitchPSMC_SFS_Together/CA/trim.",elutTimeIndex,".msmcPlot.pdf",sep=""),plot1,height=5,width=7)
}


################# northern sea otter ###############################
nso = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/MSMC_50_250DPFilter_simsInSppDIrs/output_20190215/elut2.MSMC.results.allScalings.10_14_20_MyaDivergence.txt",header=T)
for(elutTimeIndex in elut2TimeIndices){
  
  # get weighted average Ne (weighted by duration of each interval)
  weightedNeAvg <- sum(nso[nso$time_index<elutTimeIndex,]$Ne.reasonable*(nso[nso$time_index<elutTimeIndex,]$Right_generations.reasonable-nso[nso$time_index<elutTimeIndex,]$Left_generations.reasonable)/sum(((nso[nso$time_index<elutTimeIndex,]$Right_generations.reasonable-nso[nso$time_index<elutTimeIndex,]$Left_generations.reasonable))))
  plot2 <- ggplot(nso,aes(x=Left_generations.reasonable,y=Ne.reasonable))+
    geom_step(stat="identity",size=1,color=elut2Col)+
    theme_bw()+
    theme(legend.title = element_blank())+
    xlab("Generations Ago")+
    ylab("IICR (~Ne)")+
    scale_y_log10(labels=comma)+
    scale_x_log10(labels=comma)+
    theme(legend.position= c(0.6,0.85),legend.background = element_rect(fill="transparent"))+
    geom_vline(xintercept = nso[nso$time_index==elutTimeIndex-1,]$Left_generations.reasonable,linetype=2,color="grey20")+
    theme(legend.direction=("vertical"),legend.position=c(0.5,0.8),legend.background = element_rect(fill="transparent"),legend.text=element_text(size=14),legend.key.size = unit(1,"cm"))+
    annotate("text",label=paste("Time Index=",elutTimeIndex,"\nNanc=",round(nso[nso$time_index==elutTimeIndex-1,]$Ne.reasonable),"\nOldest gen = ",round(nso[nso$time_index==elutTimeIndex-1,]$Left_generations.reasonable),"\nWeighted (by interval time) Ne =",round(weightedNeAvg),sep=""),x=5000,y=12000)
  plot2
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/StitchPSMC_SFS_Together/North-AK-AL/trim.",elutTimeIndex,".msmcPlot.pdf",sep=""),plot2,height=5,width=7)
}
  
