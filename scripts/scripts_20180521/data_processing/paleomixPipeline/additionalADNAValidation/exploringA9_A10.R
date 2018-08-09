require(ggplot2)
# statsitics for Beth
# get this table from the file reads.settings in paleomix; pull out the distribution (can automate)
statsdate=20180730
plotoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/mappingStats/aDNA/",statsdate,"/",sep="")
todaysdate=format(Sys.Date(),format="%Y%m%d")

###### this isn't very well automated
A9len <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/plmx_mappingStats/A9_A10_addtlStats/A9_LengthDist_all.txt",header=T)
A10len <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/plmx_mappingStats/A9_A10_addtlStats/A10_LengthDist_all.txt",header=T)
# this is lenght dist of all reads
p1 <- ggplot(A9len,aes(x=Length,y=All-Discarded))+
  geom_bar(stat="identity")+
  scale_x_continuous(breaks=c(seq(0,300,by=25)))+
  ggtitle("A9 Read Length Distribution (All Reads)")+
  ylab("Occurences")+
  theme_bw()
p1
ggsave(paste(plotoutdir,"/A9_ReadLengthDist.all.pdf",todaysdate,".pdf",sep=""),p1,device="pdf",width = 8,height=5)

p2 <- ggplot(A10len,aes(x=Length,y=All-Discarded))+
  geom_bar(stat="identity")+
  scale_x_continuous(breaks=c(seq(0,300,by=25)))+
  ggtitle("A10 Read Length Distribution (All Reads)")+
  ylab("Occurences")+
  theme_bw()
p2
ggsave(paste(plotoutdir,"/A10_ReadLengthDist.all.pdf",todaysdate,".pdf",sep=""),p2,device="pdf",width = 8,height=5)

####################### mapped to sea otter ################
A9lenMapToElut <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/mapDamage/mapDamage_mappedToElut/A9_Elut_CA_AN_388_SN1_1a/lgdistribution.txt",header=T)
A10lenMapToElut <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/mapDamage/mapDamage_mappedToElut/A10_Elut_CA_NIC_1_SN1_1a//lgdistribution.txt",header=T)
p3 <- ggplot(A9lenMapToElut,aes(x=Length,y=Occurences))+
  geom_bar(stat="identity",fill="forestgreen")+
  scale_fill_manual(values="lightblue")+
  scale_x_continuous(breaks=c(seq(0,300,by=50)))+
  ggtitle("A9 Read Length Distribution (Mapped to Sea Otter)")+
  theme_bw()+
  ylab("Occurences")
p3
ggsave(paste(plotoutdir,"/A9_ReadLengthDist.mappedToSeaOtter.pdf",todaysdate,".pdf",sep=""),p3,device="pdf",width = 8,height=5)

p4 <- ggplot(A10lenMapToElut,aes(x=Length,y=Occurences))+
  geom_bar(stat="identity",fill="forestgreen")+
  scale_fill_manual(values="lightblue")+
  scale_x_continuous(breaks=c(seq(0,300,by=50)))+
  ggtitle("A10 Read Length Distribution (Mapped to Sea Otter)")+
  theme_bw()+
  ylab("Occurences")
p4
ggsave(paste(plotoutdir,"/A10_ReadLengthDist.mappedToSeaOtter.pdf",todaysdate,".pdf",sep=""),p4,device="pdf",width = 8,height=5)
