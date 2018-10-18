require(ggplot2)
# these date are from the 20180608 gt calls that includes capture 2 (HISeq4000) and the medgenome samples (4 v low coverage already removed prior to gt calling)
gtDate=20180806
data <- read.table(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/missingnessPerInd/",gtDate,"/noCall_per_Ind_all_5_passingFilters_raw_variants.txt",sep=""),header=T)
plotoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/MissingnessPerInd/",gtDate,"/",sep="")
tableoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/missingnessPerInd/",gtDate,sep="")
dir.create(outdir)
head(data)

mn = mean(data$NoCallCount) # 51438600
std = sd(data$NoCallCount) # 11336774
data$status <- "NA"
data[data$NoCallCount <= mn,]$status <- "less than mean"
data[data$NoCallCount > mn,]$status <- "greater than mean"
data[data$NoCallCount >= mn+std,]$status <- "greater than mean + sd"
#data[data$NoCallCount > mn+2*std,]$color <- "greater than mean + 2sd"

# initial plot of all samples
p1 <-ggplot(data,aes(x=Sample,y=NoCallCount,fill=status))+
  geom_bar(stat="identity")+
  coord_flip()+
  geom_hline(yintercept = mn)+
  geom_hline(yintercept = (mn+std),linetype="dashed")
p1
ggsave(paste(plotoutdir,"/MissingnessPerInd.pdf",sep=""),device="pdf",height=20,width=10)

# exclude those with no-call GTs greater than 1std dev above the mean
exclude <- data[data$NoCallCount > mn + std,]
exclude
write.table(exclude,paste(tableoutdir,"/exclude.Missingness.GT.MeanPlus1SD.txt",sep=""),row.names = F,quote=F,sep="\t")
