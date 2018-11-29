# this is the result of 
# zcat all_1_TrimAlt_raw_variants.vcf.gz | grep -v "#" | awk '{print $6}' > all_1_TrimAlt_raw_variants.DP.dist.txt #### THIS ISN'T DP! THIS IS QUALITY! OOPS.
# going to plot DP and then delete the file because it is huge
########################## Plot DP #########################
genotypeDate="20180806"
scaffold="GL896899.1"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/DP_QD_DIST/",genotypeDate,"/",sep="")
DPFile=paste(data.dir,scaffold,".all_1_TrimAlt_raw_variants.DP.dist.txt.gz",sep="")
QDFile=paste(data.dir,scaffold,".all_1_TrimAlt_raw_variants.QD.dist.txt.gz",sep="")
QDFile_postFilter=paste(data.dir,scaffold,".snp_2_Filter_TrimAlt_raw_variants.QD.dist.txt.gz",sep="")
DPdata <- read.table(DPFile,stringsAsFactors = F)
colnames(DPdata) <- "DP"
QDdata1 <- read.table(QDFile,stringsAsFactors = F)
colnames(QDdata1) <- "QD"
QDdata2 <- read.table(QDFile_postFilter,stringsAsFactors = F)
colnames(QDdata2) <- "QD"
### too big!! need to subsample before hand # randomly subsample:
## what to do with . <- replace with 0
require(ggplot2)
length(data[data$DP>500 & data$DP <10000,])
######################## DP Plot #############################

p0 <- ggplot(DPdata, aes(x=as.numeric(DP)))+
  geom_histogram(stat="bin",bins=100)+
  geom_vline(xintercept = 500,color="red")+
  scale_x_continuous(breaks=c(seq(0,21000,1000)))+
  theme_bw()+
  xlab("overall DP per site")
p0
ggsave(paste(data.dir,"DP.dist.",scaffold,".TrimAlt.pdf",sep=""),p0,device="pdf",width=7,height=5)
######################## QD Plot #############################
p1 <- ggplot(QDdata1, aes(x=as.numeric(QD),color="preDPFiltering"))+
  geom_density()+
  geom_density(data=QDdata2,aes(x=as.numeric(QD),color="postFiltering"))+
  geom_vline(xintercept = 2,color="red")+
  geom_vline(xintercept = 10,color="dodgerblue")+
  theme_bw()+
  xlab("Quality by Depth")+
  ggtitle("QD Distribution, pre and post minDP500 filter")
p1
ggsave(paste(data.dir,"QD.dist.",scaffold,".TrimAltAndPostDP500.pdf",sep=""),p1,device="pdf",width=7,height=5)

######################## QD Plot #############################

