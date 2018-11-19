######### Explore mapping:
require(ggplot2)
data <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/plmx_mappingStats/20180730/plmx.mappingStats.manuallyCombined.keptMarked.20181115.txt",header=T)
dataSubset <- data[data$statistic %in% c("seq_retained_reads","hits_unique","hits_clonality"),]
ggplot(dataSubset,aes(x=sample,y=value,fill=reasonExcluded))+
  geom_bar(stat="identity")+
  facet_wrap(~statistic,scales="free")+
  coord_flip()

mean(data[data$statistic=="seq_retained_reads",]$value) # 26507598 avg reads per sample 
