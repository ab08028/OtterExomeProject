############## Plot zooROH stuff without running it again ###################
require(RZooRoH)
require(R.utils)
require(ggplot2)
require(scales)
require(RColorBrewer)

colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[1],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5])
wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/zooROH/"
### converted from vcf > ped > oxford Gen , which is the "GP" format in ZooROH
pops=c("CA","AK","AL","COM","KUR") # alreadydid AK 
modelName="mix10R"
### load in each population's results and assign to "pop_results":
for(pop in pops){
  data.dir=paste(wd,pop,"/",modelName,"/",sep="")
  load(paste(data.dir,pop,".ZooROH.Results.RData",sep="")) # load in the R data ; will be called "results" for now which isn't amazing.
  assign(paste(pop,"_results",sep=""), results)
}
popResultsList=list(CA=CA_results,AK=AK_results,AL=AL_results,COM=COM_results,KUR=KUR_results)
####### plots for all pops together ########
pdf(paste(wd,"allPopsTogether.zooROH.outputPlots.",modelName,".pdf",sep=""),height = 8,width=11)

zooplot_prophbd(popResultsList,style='barplot')

zooplot_prophbd(popResultsList,style='lines', cumulative = TRUE)

zooplot_individuals(popResultsList, cumulative = TRUE)

zooplot_partitioning(popResultsList, plotids = FALSE,ylim=c(0,0.5), nonhbd = FALSE)

dev.off()
###### Plot length dist per class:

AK_hbd <- AK_results@hbdseg
rates <- AK_results@krates[1,]
p1 <- ggplot(AK_hbd,aes(x=length/1e06,fill=as.factor(HBDclass),color=as.factor(HBDclass)))+
  geom_density(alpha=0.4)+
  theme_bw()+
  xlab("Length (Mb)")+
  ggtitle("Length distribution per class")+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,10,15,20))
ggsave(paste(wd,"AK.LengthDistPerClass.pdf",sep=""),height=4,width=7)
################# Plot actual segments for a scaffold or two ##########
zooplot_hbdseg(popResultsList, chr=1,coord = c(0,50e6))

############## Plot totals? #############
#AK_results@hbdseg
# want per individual and per class to be separated 
########## classify by size:
allPopLengthSummaries <- data.frame()
for(pop in pops){
  data.dir=paste(wd,pop,"/",modelName,"/",sep="")
  load(paste(data.dir,pop,".ZooROH.Results.RData",sep="")) # load in the R data ; will be called "results" for now which isn't amazing.
  results@hbdseg$category <- NA
  results@hbdseg[results@hbdseg$length < 5e6,]$category <- "< 5Mb"
  results@hbdseg[results@hbdseg$length >= 5e6 & results@hbdseg$length < 10e6,]$category <- "5-10Mb"
  results@hbdseg[results@hbdseg$length >= 10e6,]$category <- ">=10Mb"
  # another cateogry that just splits >< 5Mb:
  results@hbdseg$category2 <- NA
  perIndTotalBPInROH <- results@hbdseg %>%
    group_by(id,category) %>%
    summarise(sumLen=sum(length))
  perIndTotalBPInROH$pop <- pop
  allPopLengthSummaries <- rbind(allPopLengthSummaries,data.frame(perIndTotalBPInROH))

}

# order factors:
allPopLengthSummaries$category <- factor(allPopLengthSummaries$category,levels=c("< 5Mb","5-10Mb",">=10Mb"))
# ggplot(allPopLengthSummaries,aes(group=interaction(as.factor(id),pop),fill=pop,x=category,y=sumLen/1e06))+
#   geom_bar(stat="identity",position="dodge")+
#   ggtitle("MB of genome contained in runs of different size categories")+
#   ylab("MB")+
#   theme_bw()

lengthPlot <- ggplot(allPopLengthSummaries,aes(x=interaction(as.factor(id),pop),fill=category,y=sumLen/1e06))+
  geom_bar(stat="identity",position="stack")+
  ggtitle("MB of genome contained in runs of different size categories")+
  ylab("MB")+
  xlab("individuals")+
  theme_bw()+
  facet_wrap(~pop,nrow=1,scales="free_x")+
  theme(axis.text.x=element_blank())
lengthPlot
ggsave(paste(wd,"ROH.SeparatedByLengthCategories.pdf",sep=""),lengthPlot,height=6,width=9)

### get avg per pop: 
meanLengthsPerCategory <- allPopLengthSummaries %>%
  group_by(pop,category) %>%
  summarise(meanLen=mean(sumLen))



meanPlot <- ggplot(meanLengthsPerCategory,aes(x=pop,fill=category,y=meanLen/1e06))+
  geom_bar(stat="identity",position="stack")+
  ggtitle("mean MB of genome contained in runs of different size categories")+
  ylab("mean MB")+
  xlab("populations")+
  theme_bw()
meanPlot
ggsave(paste(wd,"meanMBROH.PerPop.SeparatedByLengthCategories.pdf",sep=""),meanPlot,height=6,width=9)

meanLengthsPerCategory$pop   <- factor(meanLengthsPerCategory$pop, levels=c("CA","AK","AL","COM","KUR"))
meanPlot2 <- ggplot(meanLengthsPerCategory,aes(fill=pop,x=category,y=meanLen/1e06))+
  geom_bar(stat="identity",position="dodge")+
  ggtitle("mean MB of genome contained in runs of different size categories")+
  ylab("mean MB")+
  xlab("ROH")+
  theme_bw()+
  scale_fill_manual(values=unlist(colors))
meanPlot2
ggsave(paste(wd,"meanMBROH.PerPop.SeparatedByCategories.2.pdf",sep=""),meanPlot2,height=6,width=9)

######### just do all ROH > 5Mb: ########
allPopLengthSummaries2 <- data.frame()
for(pop in pops){
  data.dir=paste(wd,pop,"/",modelName,"/",sep="")
  load(paste(data.dir,pop,".ZooROH.Results.RData",sep="")) # load in the R data ; will be called "results" for now which isn't amazing.
  results@hbdseg$category2 <- NA
  results@hbdseg[results@hbdseg$length < 5e6,]$category2 <- "< 5Mb"
  results@hbdseg[results@hbdseg$length >= 5e6,]$category2 <- ">=5Mb"

  perIndTotalBPInROH <- results@hbdseg %>%
    group_by(id,category2) %>%
    summarise(sumLen=sum(length))
  perIndTotalBPInROH$pop <- pop
  allPopLengthSummaries2 <- rbind(allPopLengthSummaries2,data.frame(perIndTotalBPInROH))
  
}

# order factors:
allPopLengthSummaries2$category2 <- factor(allPopLengthSummaries2$category2,levels=c("< 5Mb",">=5Mb"))


### get avg per pop: 
meanLengthsPerCategory2 <- allPopLengthSummaries2 %>%
  group_by(pop,category2) %>%
  summarise(meanLen=mean(sumLen))


meanPlot3 <- ggplot(meanLengthsPerCategory2,aes(fill=pop,x=category2,y=meanLen/1e06))+
  geom_bar(stat="identity",position="dodge")+
  ggtitle("mean MB of genome contained in runs of different size categories")+
  ylab("mean MB")+
  xlab("ROH")+
  theme_bw()+
  scale_fill_manual(values=unlist(colors))
meanPlot
ggsave(paste(wd,"meanMBROH.PerPop.SeparatedByCategories.2.pdf",sep=""),meanPlot2,height=6,width=9)