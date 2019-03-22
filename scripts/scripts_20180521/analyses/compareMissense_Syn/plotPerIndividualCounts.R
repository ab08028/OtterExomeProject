# process each table
genotypeDate="20181119"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/compareMissense_Syn/",genotypeDate,"/countsOfGenotypesPerIndividual/",sep="")

categories=c("missense","syn")
df=data.frame()
# for(category in categories){
#   input = read.table(paste(data.dir,category,".countsPerIndividual.countOfHomAltRefHet.txt",sep=""),header=T)
#   input$category <- category
#   input$totalDerivedAlleles <- (input$HomAltCount*2) + input$HetCount # should I have excluded sites that are 1/1 for everybody?
#   # add to overall dataframe
#   df=rbind(input,df)
# }
missense=read.table(paste(data.dir,"missense.countsPerIndividual.countOfHomAltRefHet.txt",sep=""),header=T)
missense$category <- "missense"
missense$totalDerivedAlleles <- (missense$HomAltCount*2) + missense$HetCount

syn=read.table(paste(data.dir,"syn.countsPerIndividual.countOfHomAltRefHet.txt",sep=""),header=T)
syn$category <- "synonymous"
syn$totalDerivedAlleles <-  (syn$HomAltCount*2) + syn$HetCount


########## merge syn and missense ###########
both <- merge(syn,missense,by="individual",suffixes = c(".synonymous",".missense"))
both$NS_S <- both$totalDerivedAlleles.missense/both$totalDerivedAlleles.synonymous
# do I need to do some sort of called site correction??????? (or would it cancel out???)
# deal with neutral separately
# get populations
both$population <- unlist(lapply(strsplit(as.character(both$individual),"_"),"[",3))

##################### neutral ##################################
neutral = read.table(paste(data.dir,"neutral.countsPerIndividual.countOfHomAltRefHet.txt",sep=""),header=T)
neutral$category <- "neutral"
neutral$totalDerivedAlleles <-  (neutral$HomAltCount*2) + neutral$HetCount
neutral$totalCalledSites <- neutral$HomAltCount + neutral$HomRefCount + neutral$HetCount # MAKE SURE THE COUNTS INCLUDE *all* homRef sites including those that are totally monomorphic, otherwise the denominator is too small
neutral$heterozygosity <- neutral$HetCount/neutral$totalCalledSites

#### merge with both
all <- merge(both,neutral,by="individual",suffixes=c("",".neutral"))
### ############# try plotting ###########
require(ggplot2)
p1 <- ggplot(all,aes(x=heterozygosity,y=NS_S,color=population))+
  geom_point()+
  theme_bw()+
  ggtitle("Nonsynonymous to synonymous ratio vs. neutral heterozygosity per individual")+
  ylab("Missense / Synonymous Derived Alleles Ratio")+
  xlab("Neutral heterozygosity")
p1
ggsave(paste(data.dir,"NonSynToSynRatio.vs.NeutralHet.PerInd.pdf",sep=""),p1,device="pdf",width=7,height=5)

# based on clare's paper, trying to get at NS/S *heterozygosity*
# not entirely sure what to do for this (talk to Kirk)
# but going to try NS hets / S hets, with the idea that correcting for total called cds sites (denominator) would cancel out. However more proper denominator would be the total called sites of each type -- should I classify these? not sure...
#
all$NSHets_SHets <- all$HetCount.missense/all$HetCount.synonymous
p2 <- ggplot(all,aes(x=heterozygosity,y=NSHets_SHets,color=population))+
  geom_point()+
  theme_bw()+
  ggtitle("Nonsynonymous to synonymous ratio vs. neutral heterozygosity per individual")+
  ylab("Missense Hets / Synonymous Hets Ratio")+
  xlab("Neutral heterozygosity")
p2
ggsave(paste(data.dir,"NonSynToSynHetRatio.vs.NeutralHet.PerInd.pdf",sep=""),p2,device="pdf",width=7,height=5)

############# So what do I want to simulate?