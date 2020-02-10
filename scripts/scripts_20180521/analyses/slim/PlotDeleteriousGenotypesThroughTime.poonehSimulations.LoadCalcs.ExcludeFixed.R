colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[1],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5])
textsize=14
require(ggplot2)
require(dplyr)
todaysdate=format(Sys.Date(),format="%Y%m%d")

data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/concattedSummaries/"
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/loadPlots/"
dir.create(plot.dir)
popModDates=c("AK/1D.3Epoch.LongerRecovery/20191202/","CA/1D.3Epoch.LongerRecovery/20191013/","AK/1D.5Epoch/20190802/")
reps=c(seq(1,25)) # some reps don't make it through Hoffman; so I have a file.exists() test in the loop to skip reps that didn't yield output
hset=c(0,0.5)
# note these simulations start at gen 0, when pop is at 4500 (taht's the only gen with that pop size)
# so you want to only find things that are fixed at gen 0.
#states=c("PreContraction","PostContraction")
allAvgdInputs=data.frame()
allLoads=data.frame() 
# get avg homozygous derived per individual : (get hets too?)
for(popModDate in popModDates){
  for(rep in reps){
    print(rep)
    for(h in hset){
      # check if rep exists (some have random hoffman failures)
      infile=paste(data.dir,popModDate,"/h_",h,"/replicate_",rep,".slim.output.allConcatted.summary.txt.gz",sep="")
      if(file.exists(infile)){
        input = read.table(infile,sep=",",header=T)
        # want to exclude sites that are at frequency 1 *prior* to the bottleneck
        pop= unlist(lapply(strsplit(popModDate,"/"),"[",1))
        # any sites that are at frequency 1 prior to the bottleneck (at gen 0) should be excluded from load calcs. If they are at frequency 1 only after the bottleneck they can stay to be part of load. This will be if geneartion <= bneckGen (model specific) and if the numhom==popsizeDIP. can't just go by mutid because those can be duplicated across chunks. Must go by chunk AND mutID within a replicate. Then need to remove them at all other time points.... 
        # give each mutation a unique ID that is their chunk (ie chromosome #) and mutid
        # you need this because mutIDs can be duplicated between chunks, but not within a chunk
        input$chunk.mutID <- paste(input$chunk,".",input$mutid,sep="")
        # THIS DOESNT WORK AS EXPECTED
        fixedToRemove <- input[(input$gen==0 & input$numhom==input$popsizeDIP),]$chunk.mutID # ~4000 sites per replicate. cool
        # should each be a unique value:
        length(unique(fixedToRemove))==length(fixedToRemove)
        #exclude those sites:
        inputWithFixedRemoved <- input[!(input$chunk.mutID %in% fixedToRemove),]
        # see how this changes:
        length(unique(input$chunk.mutID)) - length(unique(inputWithFixedRemoved$chunk.mutID))==length(fixedToRemove)# should equal total to remove TRUE
        # should be none left:
        #inputWithFixedRemoved[inputWithFixedRemoved$generation==0 & inputWithFixedRemoved$numhom==inputWithFixedRemoved$popsizeDIP,] # there should be no sites left.
        # want to get load:
        inputWithFixedRemoved$qFreq <- (inputWithFixedRemoved$numhet + (2*inputWithFixedRemoved$numhom)) / (2*inputWithFixedRemoved$popsizeDIP)
        inputWithFixedRemoved$pFreq <- 1 - inputWithFixedRemoved$qFreq
        inputWithFixedRemoved$loadComponent <- (2*h*abs(inputWithFixedRemoved$s)*inputWithFixedRemoved$qFreq*inputWithFixedRemoved$pFreq) + (abs(inputWithFixedRemoved$s)*((inputWithFixedRemoved$qFreq)^2))
        # want to categorize by s
        inputWithFixedRemoved$popModDate <- popModDate
        inputWithFixedRemoved$sCat <- NA
        # JAR categories from her 2018 paper
        inputWithFixedRemoved[inputWithFixedRemoved$s >= -1 & inputWithFixedRemoved$s < -0.01,]$sCat <- "strongly deleterious"
        inputWithFixedRemoved[inputWithFixedRemoved$s >= -0.01 & inputWithFixedRemoved$s < -0.001,]$sCat <- "moderately deleterious"
        inputWithFixedRemoved[inputWithFixedRemoved$s >= -0.001 & inputWithFixedRemoved$s < 0,]$sCat <- "weakly deleterious"
        inputWithFixedRemoved[inputWithFixedRemoved$s==0,]$sCat <- "neutral"
        model= unlist(lapply(strsplit(popModDate,"/"),"[",2))
        date= unlist(lapply(strsplit(popModDate,"/"),"[",3))
        inputWithFixedRemoved$population <- pop
        #input$state <- state
        inputWithFixedRemoved$h <- h
        inputWithFixedRemoved$model <- model
        inputWithFixedRemoved$date <- date
        inputWithFixedRemoved$replicate <- rep
        # want to get total S per generation:
        # want to get totals and avgs across all chunks per generation
        avgHomPerIndPersCat <- inputWithFixedRemoved %>%
          group_by(generation,h,population,sCat,model,replicate,popModDate,popsizeDIP) %>%
          summarise(totalNumHom=sum(numhom),totalHet=sum(numhet)) %>%
          mutate(avgHomPerInd=totalNumHom/popsizeDIP) %>%
          mutate(avgHetPerInd=totalHet/popsizeDIP) %>%
          mutate(avgDerivedAllelesPerInd=((2*totalNumHom)+totalHet)/(2*popsizeDIP))
        # don't group by sCat for load calcs:
        LoadPerGeneration <- inputWithFixedRemoved %>%
          group_by(generation,h,population,model,replicate,popModDate,popsizeDIP) %>%
          summarise(totalS=sum(loadComponent))
        LoadPerGeneration$W <- exp(-LoadPerGeneration$totalS)
        LoadPerGeneration$L  = 1 - LoadPerGeneration$W # mutation load 
        
        allAvgdInputs=rbind(allAvgdInputs,data.frame(avgHomPerIndPersCat))
        allLoads = rbind(allLoads,data.frame(LoadPerGeneration))
      }}}}

# label h:
allAvgdInputs$hLabel <- ""
allAvgdInputs[allAvgdInputs$h==0,]$hLabel <- "h = 0 (rec.)"
allAvgdInputs[allAvgdInputs$h==0.5,]$hLabel <- "h = 0.5 (add.)"
allLoads$hLabel <- ""
allLoads[allLoads$h==0,]$hLabel <- "h = 0 (rec.)"
allLoads[allLoads$h==0.5,]$hLabel <- "h = 0.5 (add.)"
## want to write this out as a table so don't have to do it multiple times ###
write.table(allAvgdInputs,paste(data.dir,"AvgHomozygousDerivedGTs.PerInd.ThroughTime.AllReps.CA.AK.RemovedPreFurTradeFixedVar.txt",sep=""),row.names = F,col.names = T,quote=F,sep="\t")
write.table(allLoads,paste(data.dir,"LoadPerGeneration.RemovedPreFurTradeFixedVars.ThroughTime.AllReps.CA.AK.txt",sep=""),row.names = F,col.names = T,quote=F,sep="\t")

### need to get average across number of reps:
avgHomPerIndPersCat_AvgAcrossReps <- allAvgdInputs %>%
  group_by(generation,h,hLabel,population,sCat,model,popModDate,popsizeDIP) %>%
  summarise(avgHomPerInd_overAllReps=mean(avgHomPerInd))
# order factors:
#avgHomPerIndPersCat$state <- factor(avgHomPerIndPersCat$state,levels=c("PreContraction","PostContraction"))
allAvgdInputs$sCat <- factor(allAvgdInputs$sCat,levels=c("neutral","weakly deleterious","moderately deleterious","strongly deleterious"))
avgHomPerIndPersCat_AvgAcrossReps$sCat <- factor(avgHomPerIndPersCat_AvgAcrossReps$sCat,levels=c("neutral","weakly deleterious","moderately deleterious","strongly deleterious"))
modelsUsed = paste(unique(avgHomPerIndPersCat$popModDate),collapse="; ")
# aha, the "collapse" option in paste works for this to collapse list elements and sep by something:
pops=c("AK","CA")
########## do recessive and additive separately #######
for(pop in pops){
  for(h in c(0,0.5)){
    p1 <- ggplot(allAvgdInputs[allAvgdInputs$h==h & allAvgdInputs$population==pop,],aes(y=avgHomPerInd,x=generation))+
      geom_line(size=0.3,alpha=0.5,aes(group=replicate,color=population))+
      geom_line(data=avgHomPerIndPersCat_AvgAcrossReps[avgHomPerIndPersCat_AvgAcrossReps$h==h & avgHomPerIndPersCat_AvgAcrossReps$population==pop,],aes(y=avgHomPerInd_overAllReps,x=generation,color=population),size=1.5)+
      theme_bw()+
      facet_grid(~sCat~population,scales="free")+
      ggtitle(paste("Average number of homozygous derived genotypes per individual, per category\nneutral [0]; weak (0 to -0.001]; moderate (-0.001 to -0.01]; strong (-0.01 to -1]\nh=",h,sep=""))+
      ylab("Average Genotypes Per Individual") +
      xlab("")+
      scale_color_manual(values=unlist(colors))+
      theme(legend.position = 'none')
    p1
    ggsave(paste(plot.dir,"avgNumHomDerivedPerInd.ThroughTime.h.",h,".",pop,".ONLY.pdf",sep = ""),p1,width=7,height=7)
    
  }}
############## plot avg het per ind ############

avgHetPerIndPersCat_AvgAcrossReps <- allAvgdInputs %>%
  group_by(generation,h,hLabel,population,sCat,model,popModDate,popsizeDIP) %>%
  summarise(avgHetPerInd_overAllReps=mean(avgHetPerInd))
# order factors:
avgHetPerIndPersCat_AvgAcrossReps$sCat <- factor(avgHetPerIndPersCat_AvgAcrossReps$sCat,levels=c("neutral","weakly deleterious","moderately deleterious","strongly deleterious"))


pops=c("AK","CA")
for(pop in pops){
  for(h in c(0,0.5)){
    p2 <- ggplot(allAvgdInputs[allAvgdInputs$h==h & allAvgdInputs$population==pop,],aes(y=avgHetPerInd,x=generation))+
      geom_line(size=0.3,alpha=0.5,aes(group=replicate,color=population))+
      geom_line(data=avgHetPerIndPersCat_AvgAcrossReps[avgHetPerIndPersCat_AvgAcrossReps$h==h & avgHetPerIndPersCat_AvgAcrossReps$population==pop,],aes(y=avgHetPerInd_overAllReps,x=generation,color=population),size=1.5)+
      theme_bw()+
      facet_grid(~sCat~population,scales="free")+
      ggtitle(paste("Average number of heterozygous genotypes per individual, per category\nneutral [0]; weak (0 to -0.001]; moderate (-0.001 to -0.01]; strong (-0.01 to -1]\nh=",h,sep=""))+
      ylab("Average Genotypes Per Individual") +
      xlab("")+
      scale_color_manual(values=unlist(colors))+
      theme(legend.position = 'none')
    p2
    ggsave(paste(plot.dir,"avgNumHetDerivedPerInd.ThroughTime.h.",h,".",pop,".ONLY.pdf",sep = ""),p2,width=7,height=7)
    
  }}
##################### what else would be good? total number of alleles? #######
totalDerivedAlleles_AcrossReps <- allAvgdInputs %>%
  group_by(generation,h,hLabel,population,sCat,popModDate) %>%
  summarise(avgDerivedAcrossReps=mean(avgDerivedAllelesPerInd))

totalDerivedAlleles_AcrossReps$sCat <- factor(totalDerivedAlleles_AcrossReps$sCat,levels=c("neutral","weakly deleterious","moderately deleterious","strongly deleterious"))
for(pop in pops){
  for(h in c(0,0.5)){
    p3 <- ggplot(allAvgdInputs[allAvgdInputs$h==h & allAvgdInputs$population==pop,],aes(y=avgDerivedAllelesPerInd,x=generation,color=population))+
      geom_line(data=totalDerivedAlleles_AcrossReps[totalDerivedAlleles_AcrossReps$h==h & totalDerivedAlleles_AcrossReps$population==pop,],aes(y=avgDerivedAcrossReps,x=generation,color=population),size=1.5)+
      geom_line(size=0.2,alpha=0.5,aes(group=replicate,color=population))+
      theme_bw()+
      facet_grid(~sCat~population,scales="free")+
      xlab("")+
      scale_color_manual(values=unlist(colors))+
      theme(legend.position = 'none')+
      ggtitle(paste("Average number of derived alleles per individual, per category\nneutral [0]; weak (0 to -0.001]; moderate (-0.001 to -0.01]; strong (-0.01 to -1]\nh=",h,sep=""))+
      ylab("Average Derived Alleles Per Individual") 
    p3
    ggsave(paste(plot.dir,"avgDerivedAllelesPerInd.ThroughTime.h.",h,".",pop,".ONLY.pdf",sep = ""),p3,width=7,height=7)
  }}

