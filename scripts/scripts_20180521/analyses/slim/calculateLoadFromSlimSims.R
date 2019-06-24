require(ggplot2)
todaysdate=format(Sys.Date(),format="%Y%m%d")

data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/concattedSummaries/"
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/loadCalcs/"
dir.create(plot.dir)
#pops=c("AK","AL","genericPop.LongerContract")
#models=c("1D.2Epoch.1.5Mb.cds")
#simdates=c(20190424,20190607)
# skipping AL "AL/1D.2Epoch.1.5Mb.cds/20190424/" and CA etc -- add those in next 
popModDates=c("AK/1D.2Epoch.1.5Mb.cds/20190424/","AK/1D.2Epoch.1.5Mb.cds.LongerContract/20190607","genericPop/1D.2Epoch.1.5Mb.cds.20KAncSize/20190611") # AK and AL have dadi parameters, genericPop has parameters based on AK MLE grid that is fur-trade relevant. ### need to come up with better classification system for this. 
#reps=c(seq(1,23))
reps=c(seq(1,25)) # some reps don't make it through Hoffman; so I have a file.exists() test in the loop to skip reps that didn't yield output
hset=c(0,0.5)
states=c("PreContraction","PostContraction")
allLoads=data.frame()
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
            # calculate q (alt allele frequency) and p (ref allele frequency) per site
            # be careful about which you use in equation! s*q^2 means that q is frequency of ALT allele with associated "s". so p is freq of ref allele
            input$qFreq <- (input$p1numhet + (2*input$p1numhom)) / (2*input$popsizeDIP)
            input$pFreq <- 1 - input$qFreq
            # equation is, per site: 2hspq + sq^2 is the contribution to load; p = 1-qFreq
            # my "s" is negative, so I want to absolute value s --> |s|
            input$loadComponent <- (2*h*abs(input$s)*input$qFreq*input$pFreq) + (abs(input$s)*((input$qFreq)^2))
            # total sites:
            
            S = sum(input$loadComponent)
            W = exp(-S) # mean fitness e^-S
            L  = 1 - W # mutation load 
            #### add to dataframe: #####
            # pull out population etc from popModDate
            pop= unlist(lapply(strsplit(popModDate,"/"),"[",1))
            model= unlist(lapply(strsplit(popModDate,"/"),"[",2))
            date= unlist(lapply(strsplit(popModDate,"/"),"[",3))
            loadDF <- data.frame(population=pop)
            loadDF$model <- model
            loadDF$date <- date
            loadDF$rep <- rep
            loadDF$state <- state
            loadDF$h <- h
            loadDF$da
            loadDF$S_allsites <- S
            loadDF$W_meanFitness <- W
            loadDF$L_mutationLoad <- L
            #### combine with other reps: #####
            allLoads = rbind(allLoads,loadDF)
          }}}}}

# change order of factors:
allLoads$state <- factor(allLoads$state,levels=c("PreContraction","PostContraction"))
# label H:
allLoads$hLabel <- paste("h = ",allLoads$h)
allLoads$model <- factor(allLoads$model,levels=c("1D.2Epoch.1.5Mb.cds","1D.2Epoch.1.5Mb.cds.LongerContract", "1D.2Epoch.1.5Mb.cds.20KAncSize"))
p1 <- 
  ggplot(allLoads,aes(x=state,y=L_mutationLoad,fill=state))+
  geom_violin(position=position_dodge(.5))+
  geom_point(position=position_dodge(.5),size = 1,alpha=0.5)+
  theme_bw()+
  facet_grid(hLabel~interaction(population,model))+
  ylab("Genetic Load")+
  xlab("") +
  theme(legend.position = "none")
p1
ggsave(paste(plot.dir,"Load.perPop.PrePostContract.",todaysdate,".pdf",sep=""),p1,height=6,width=8)


##################### make a fitch figure of load : #############
forTalk <- allLoads[allLoads$population=="AK" & allLoads$model=="1D.2Epoch.1.5Mb.cds.LongerContract" & allLoads$date=="20190607",] # just want the 250 inds with 35 gens AK model
# use alaska color
colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[7],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors

head(forTalk)
forTalk$hLabel2 <- "Recessive"
forTalk[forTalk$h == "0.5",]$hLabel2 <- "Additive"
head(forTalk)

pTalk1 <- ggplot(forTalk, aes(x=state,y=L_mutationLoad,fill=state))+
  geom_violin(position=position_dodge(.5))+
  geom_point(position=position_dodge(.5),size = 1,alpha=0.5)+
  theme_bw()+
  facet_grid(~hLabel2)+
  ylab("Genetic Load")+
  xlab("") +
  theme(legend.position = "none",text=element_text(size=14),axis.text=element_text(size=14),strip.text = element_text(size=14))+
  ggtitle("Figure for talk: Simulated AK model with 250 inds for 35 gens")+
  scale_fill_manual(values=c("dodgerblue",colors$AK))
pTalk1
ggsave(paste(plot.dir,"AK.modelLongerContract.FigureForTALKS.Load.PrePostContract.",todaysdate,".pdf",sep=""),pTalk1,height=6,width=8)

# separate additive and recessive:
####### recessive:
pTalk2a <- ggplot(forTalk[forTalk$hLabel2=="Recessive",], aes(x=state,y=L_mutationLoad,fill=state))+
  geom_violin(position=position_dodge(.5))+
  geom_point(position=position_dodge(.5),size = 1,alpha=0.5)+
  theme_bw()+
  facet_grid(~hLabel2)+
  ylab("Genetic Load")+
  xlab("") +
  theme(legend.position = "none",text=element_text(size=14),axis.text=element_text(size=14),strip.text = element_text(size=14))+
  ggtitle("Figure for talk: Simulated AK model with 250 inds for 35 gens")+
  scale_fill_manual(values=c("dodgerblue",colors$AK))+
  scale_y_continuous(limits=c(0.19,0.34))
pTalk2a
ggsave(paste(plot.dir,"AK.modelLongerContract.FigureForTALKS.Load.PrePostContract.RecessiveOnly.",todaysdate,".pdf",sep=""),pTalk2a,height=6,width=4)

######## additive:
pTalk2b <- ggplot(forTalk[forTalk$hLabel2=="Additive",], aes(x=state,y=L_mutationLoad,fill=state))+
  geom_violin(position=position_dodge(.5))+
  geom_point(position=position_dodge(.5),size = 1,alpha=0.5)+
  theme_bw()+
  facet_grid(~hLabel2)+
  ylab("Genetic Load")+
  xlab("") +
  theme(legend.position = "none",text=element_text(size=14),axis.text=element_text(size=14),strip.text = element_text(size=14))+
  ggtitle("Figure for talk: Simulated AK model with 250 inds for 35 gens")+
  scale_fill_manual(values=c("dodgerblue",colors$AK))+
  scale_y_continuous(limits=c(0.19,0.34))
pTalk2b
ggsave(paste(plot.dir,"AK.modelLongerContract.FigureForTALKS.Load.PrePostContract.AdditiveOnly.",todaysdate,".pdf",sep=""),pTalk2b,height=6,width=4)

########## also make box plots because people understand them better ########
####### recessive:
pTalk3a <- ggplot(forTalk[forTalk$hLabel2=="Recessive",], aes(x=state,y=L_mutationLoad,fill=state))+
  geom_boxplot(position=position_dodge(.5))+
  #geom_point(position=position_dodge(.5),size = 1,alpha=0.5)+
  theme_bw()+
  facet_grid(~hLabel2)+
  ylab("Genetic Load")+
  xlab("") +
  theme(legend.position = "none",text=element_text(size=14),axis.text=element_text(size=14),strip.text = element_text(size=14))+
  ggtitle("Figure for talk: Simulated AK model with 250 inds for 35 gens")+
  scale_fill_manual(values=c("dodgerblue",colors$AK))
pTalk3a
ggsave(paste(plot.dir,"AK.modelLongerContract.FigureForTALKS.Load.PrePostContract.RecessiveOnly.BOXPLOT.",todaysdate,".pdf",sep=""),pTalk3a,height=6,width=4)

######## additive:
pTalk3b <- ggplot(forTalk[forTalk$hLabel2=="Additive",], aes(x=state,y=L_mutationLoad,fill=state))+
  geom_boxplot(position=position_dodge(.5))+
  #geom_point(position=position_dodge(.5),size = 1,alpha=0.5)+
  theme_bw()+
  facet_grid(~hLabel2)+
  ylab("Genetic Load")+
  xlab("") +
  theme(legend.position = "none",text=element_text(size=14),axis.text=element_text(size=14),strip.text = element_text(size=14))+
  ggtitle("Figure for talk: Simulated AK model with 250 inds for 35 gens")+
  scale_fill_manual(values=c("dodgerblue",colors$AK))
pTalk3b
ggsave(paste(plot.dir,"AK.modelLongerContract.FigureForTALKS.Load.PrePostContract.AdditiveOnly.BOXPLOT.",todaysdate,".pdf",sep=""),pTalk3b,height=6,width=4)
