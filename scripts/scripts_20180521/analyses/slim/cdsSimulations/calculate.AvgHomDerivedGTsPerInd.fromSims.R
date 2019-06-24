####### something is weird wtih 20K model right now -- showing more hom and het derived than the small pop -should it be more hets and fewer homs? ######
# stay calm
# NOTE: genericPop simulations with 20K didn't have 20 finish due to Hoffman maintenance (24 and 25 didn't run) so eventually run more of these if you need them for hte paper; but not doing for now because you may do something else
################### get avg num of homs and hets ############
require(ggplot2)
require(dplyr)
todaysdate=format(Sys.Date(),format="%Y%m%d")

data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/concattedSummaries/"
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/loadCalcs/"
dir.create(plot.dir)
#pops=c("AK","AL","genericPop.LongerContract")
#models=c("1D.2Epoch.1.5Mb.cds")
#simdates=c(20190424,20190607)
# skipping AL "AL/1D.2Epoch.1.5Mb.cds/20190424/" and CA etc -- add those in next 
#popModDates=c("AK/1D.2Epoch.1.5Mb.cds/20190424/","AK/1D.2Epoch.1.5Mb.cds.LongerContract/20190607","genericPop/1D.2Epoch.1.5Mb.cds.20KAncSize/20190611") # AK and AL have dadi parameters, genericPop 
popModDates=c("AK/1D.2Epoch.1.5Mb.cds.LongerContract/20190607","genericPop/1D.2Epoch.1.5Mb.cds.20KAncSize/20190611")
#has parameters based on AK MLE grid that is fur-trade relevant. ### need to come up with better classification system for this. 
#reps=c(seq(1,23))
reps=c(seq(1,25)) # some reps don't make it through Hoffman; so I have a file.exists() test in the loop to skip reps that didn't yield output
hset=c(0,0.5)
states=c("PreContraction","PostContraction")
allLoads=data.frame()
allInputs = data.frame()
for(popModDate in popModDates){
  #for(model in models){
  #for(simdate in simdates){
  #  for(pop in pops){
  for(rep in reps){
    for(state in states){
      for(h in hset){
        # check if rep exists (some have random hoffman failures)
        infile=paste(data.dir,popModDate,"/h_",h,"/replicate_",rep,".slim.output.",state,".allConcatted.summary.txt.gz",sep="")
        if(file.exists(infile)){
          input = read.table(infile,sep=",",header=T)
          # want to categorize by s
          input$popModDate <- popModDate
          input$sCat <- NA
          # JAR categories from her 2018 paper
          input[input$s >= -1 & input$s < -0.01,]$sCat <- "strongly deleterious"
          input[input$s >= -0.01 & input$s < -0.001,]$sCat <- "moderately deleterious"
          input[input$s >= -0.001 & input$s < 0,]$sCat <- "weakly deleterious"
          input[input$s==0,]$sCat <- "neutral"
          # calculate q (alt allele frequency) and p (ref allele frequency) per site
          # be careful about which you use in equation! s*q^2 means that q is frequency of ALT allele with associated "s". so p is freq of ref allele
          #input$qFreq <- (input$p1numhet + (2*input$p1numhom)) / (2*input$popsizeDIP)
          #input$pFreq <- 1 - input$qFreq
          # equation is, per site: 2hspq + sq^2 is the contribution to load; p = 1-qFreq
          # my "s" is negative, so I want to absolute value s --> |s|
          #input$loadComponent <- (2*h*abs(input$s)*input$qFreq*input$pFreq) + (abs(input$s)*((input$qFreq)^2))
          # total sites:
          
          #S = sum(input$loadComponent)
          #W = exp(-S) # mean fitness e^-S
          #L  = 1 - W # mutation load 
          #### add to dataframe: #####
          # pull out population etc from popModDate
          pop= unlist(lapply(strsplit(popModDate,"/"),"[",1))
          model= unlist(lapply(strsplit(popModDate,"/"),"[",2))
          date= unlist(lapply(strsplit(popModDate,"/"),"[",3))
          input$population <- pop
          input$state <- state
          input$h <- h
          #loadDF <- data.frame(population=pop)
          input$model <- model
          input$date <- date
          #loadDF$rep <- rep
          #loadDF$state <- state
          #loadDF$h <- h
          #loadDF$da
          #loadDF$S_allsites <- S
          #loadDF$W_meanFitness <- W
          #loadDF$L_mutationLoad <- L
          #### combine with other reps: #####
          #allLoads = rbind(allLoads,loadDF)
          allInputs=rbind(allInputs,input)
          
        }}}}}

# label h:
allInputs$hLabel <- ""
allInputs[allInputs$h==0,]$hLabel <- "h = 0 (rec.)"
allInputs[allInputs$h==0.5,]$hLabel <- "h = 0.5 (add.)"

###### tallying up the number of homozygous genotypes per category divided by the population size before and after contraction for each replicate
# 
avgHomPerIndPersCat <- allInputs %>%
  group_by(state,h,hLabel,population,sCat,replicate,model,popModDate) %>%
  tally(p1numhom/popsizeDIP)

# order factors:
avgHomPerIndPersCat$state <- factor(avgHomPerIndPersCat$state,levels=c("PreContraction","PostContraction"))
avgHomPerIndPersCat$sCat <- factor(avgHomPerIndPersCat$sCat,levels=c("neutral","weakly deleterious","moderately deleterious","strongly deleterious"))
modelsUsed = paste(unique(avgHomPerIndPersCat$popModDate),collapse="; ")
# aha, the "collapse" option in paste works for this to collapse list elements and sep by something:

p1 <- ggplot(avgHomPerIndPersCat,aes(y=n,x=state,fill=model,group=interaction(popModDate,state)))+
  geom_violin(position=position_dodge(0.5))+
  geom_point(position=position_dodge(0.5),size=0.3)+
  theme_bw()+
  facet_grid(~sCat~hLabel,scales="free")+
  ggtitle(paste("Average number of homozygous derived genotypes per individual, per category\nneutral [0]; weak (0 to -0.001]; moderate (-0.001 to -0.01]; strong (-0.01 to -1]\nSimulations: ",modelsUsed,sep=""))+
  ylab("Average Genotypes Per Individual") +
  xlab("")
p1
ggsave(paste(plot.dir,"avgNumHomDerivedPerInd.pdf",sep = ""),p1,width=7,height=7)

############## get avg het per ind ############
avgHetPerIndPersCat <- allInputs %>%
  group_by(state,h,hLabel,population,model,sCat,replicate,popModDate) %>%
  tally(p1numhet/popsizeDIP)

# order factors:
avgHetPerIndPersCat$state <- factor(avgHetPerIndPersCat$state,levels=c("PreContraction","PostContraction"))
avgHetPerIndPersCat$sCat <- factor(avgHetPerIndPersCat$sCat,levels=c("neutral","weakly deleterious","moderately deleterious","strongly deleterious"))
# get list of models for the title:
modelsUsed = paste(unique(avgHetPerIndPersCat$popModDate),collapse="; ")
# aha, the "collapse" option in paste works for this to collapse list elements and sep by something:

p2 <- ggplot(avgHetPerIndPersCat,aes(y=n,x=state,fill=model,group=interaction(popModDate,state)))+
  geom_violin(position=position_dodge(0.5))+
  geom_point(position=position_dodge(0.5),size=0.3)+
  theme_bw()+
  facet_grid(~sCat~hLabel,scales="free")+
  theme(legend.position = "NA")+
  ggtitle(paste("Average number of heterozygous genotypes per individual, per category\nneutral [0]; weak (0 to -0.001]; moderate (-0.001 to -0.01]; strong (-0.01 to -1]\nSimulations: ",modelsUsed,sep=""))+
  ylab("Average Genotypes Per Individual") +
  xlab("")
p2
ggsave(paste(plot.dir,"avgNumHetPerInd.pdf",sep = ""),p2,width=7,height=7)

##################### what else would be good? total number of alleles? #######
totalDerivedAlleles <- allInputs %>%
  group_by(state,h,hLabel,population,sCat,replicate,popModDate) %>%
  tally((2*p1numhom+p1numhet)/(2*popsizeDIP))

p3 <- ggplot(totalDerivedAlleles,aes(y=n,x=state,fill=state,group=interaction(popModDate,state)))+
  geom_violin(position=position_dodge(0.5))+
  geom_point(position=position_dodge(0.5),size=0.3)+
  theme_bw()+
  facet_grid(~sCat~hLabel,scales="free")+
  theme(legend.position = "NA")+
  ggtitle("Total number of derived ALLELES per individual")
p3
# so maybe what we saw in sea otter vs giant otter is due to longer-term trends?
# try with 20K vs 5K vs 250- or run for a really long time or something? did I even plot this right?