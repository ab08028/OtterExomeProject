######## Plink ROH's -- dubious if this will work on capture data. 
# how do I normalize for number of callable sites (do I?)
require(dplyr)
require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=8,name = "Dark2")
# make Baja brown
colors=list(CA=colorPal[1],BAJ=colorPal[7],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors
pops=c("CA","AK","AL","COM","KUR") # excluding BAJA because too low sample size
allPopsROHs <- data.frame()
wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/ROH_plink/"
for(pop in pops){
  data.dir=paste(wd,pop,"/",sep="")
  input <- read.table(paste(data.dir,"plink.ROH.",pop,".hom",sep=""),header=T)
  # want to focus just on large ROHs per individual
  input$pop <- pop
  allPopsROHs <- rbind(allPopsROHs,input)
}
p1 <- ggplot(allPopsROHs, aes(x=IID,y=KB,color=pop))+
  geom_point()+
  geom_hline(yintercept = 5000,color="blue")+
  geom_hline(yintercept = 10000,color="red")+
  theme_bw()+
  theme(axis.text.x=element_blank())
p1  
## do MB of genome contained in long ROH? or full on Froh? 
LongROHCutoffKB=5000 # in kb, so is 10Mb
allPopsROHs_long <- allPopsROHs[allPopsROHs$KB >=LongROHCutoffKB,]
# summarize KB contained in LongROH per ind:
allPopsROHS_long_summary <- allPopsROHs_long %>%
  group_by(IID) %>%
  mutate(totalMBInLongROH=sum(KB)/1000) # divide by 1000 to convert KB to MB
# e
p2 <- ggplot(allPopsROHS_long_summary,aes(x=pop,y=totalMBInLongROH,fill=pop))+
  geom_boxplot()+
  geom_point()+
  theme_bw()+
  ylab(paste("Total MB contained per ind.\nin ROHs longer than ",LongROHCutoffKB/1000," MB",sep=""))+
  xlab("")+
  scale_fill_manual(values=c(colors['AK'],colors['AL'],colors['CA'],colors['COM'],colors['KUR']))+
  theme(legend.position = "none")+
  ggtitle("MB contained in long ROH")

p2
# divide by genome size or not? 
ggsave(paste(wd,"MBContainedInLongROH.perpop.pdf",sep=""),p2,height=5,width=7)
####### try plotting by category of ROH? ###########
####### does Froh make any sense with capture data? ########
