########### get scaff sizes and exclude <5Mb ###########
mustelaChrSizes <- (read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/SlidingWindowheterozygosity/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.224.ChrLengths.txt",sep=",",stringsAsFactors = F,strip.white = T))
mustelaChrSizes <- data.frame(t(mustelaChrSizes))
colnames(mustelaChrSizes) <- "record"
mustelaChrSizes$scaff <- unlist(lapply(strsplit(as.character(mustelaChrSizes$record),":"),"[",1))
mustelaChrSizes$size <- unlist(lapply(strsplit(as.character(mustelaChrSizes$record),":"),"[",2))
##### kirk's solution: do for one population, no fixed sites, and with projection's appropriate missingness threshold (done in vcftools)
# did it for alaska:
############ AK ##########
AK <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/basicSNPStats/DistOfNeutralSNPs/AK.Mac1.maxMiss4.recode.vcf.gz")
head(AK)
colnames(AK) <- c("Scaffold","Position")
p1 <- ggplot(AK,aes(x=Position/1e6,y=Scaffold,color=Scaffold))+
  geom_point(size=.5,alpha=0.5)+
  theme_bw()+
  theme(axis.text.y=element_blank(),legend.position = 'none')+
  ylab("Scaffold")+
  xlab("Pos (Mb)")+
  ggtitle("AK")+
  scale_color_manual(values=rep(c("tomato","dodgerblue"),length(unique(AK$Scaffold))+2/2))
p1
ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/basicSNPStats/DistOfNeutralSNPs/AK.NeutralSNPDist.pdf",p1,device="pdf",height=7,width=5)
#########
require(GenomicRanges)
grangeX <- makeGRangesFromDataFrame(AK,seqnames.field="Scaffold",start.field = "Position",end.field = "Position",ignore.strand = T)
# granges object
mustelaChrSizes$start <- as.character(1)
grangMustela <- makeGRangesFromDataFrame(mustelaChrSizes,seqnames.field="scaff",start.field="start",end.field="size")
tiles <- tile(grangMustela,width=1000) # tile the genome into 1000 bp intervals
length(unlist(tiles)) # need to unlist it or it remains grouped by scaffold!!!
########### REALLY IMPORTANT NOTE HERE: must "unlist" tiles!!! otherwise you just get counts per scaffold!!!!! ########## https://support.bioconductor.org/p/100222/
##### pay attention! make sure counts make sense!
overlaps <- data.frame(unlist(countOverlaps(query=unlist(tiles),subject=grangeX)))
colnames(overlaps) <- "snpsPer1KbBin"
overlapsNoZero <- overlaps[overlaps$snpsPer1KbBin>0,] # condition on at least one snp being present
overlapsNoZero <- data.frame(overlapsNoZero)
colnames(overlapsNoZero) <- "snpsPer1KbBin"

length(overlapsNoZero$snpsPer1KbBin) # 2139 1kb bins that contain >= 1 snp
avg = mean(overlapsNoZero$snpsPer1KbBin) # avg 1.2 snps/1 kb window (conditioning on there being a snp in a window)
avg # avg 1.2 snps per window
sd = sd(overlapsNoZero$snpsPer1KbBin) 
avg+1.96*sd
length(overlapsNoZero[overlapsNoZero$snpsPer1KbBin>3,])
#ggplot(overlapsNoZero,aes(x=snpsPer1KbBin))+
#  geom_density()

###### dist between hits? 
# distToNearest <- distanceToNearest(grangeX)
# length(distToNearest)
# distToNearest_df <- data.frame(distToNearest)
# mean(distToNearest_df$distance)
# ggplot(distToNearest_df,aes(x=distance))+
#   geom_density()

### Select the windows that have an overlap # this subsets tiles by the overlap hits
WindowsWithHits <- subsetByOverlaps(unlist(tiles),grangeX)
# then get distance between those windows
WindowsWithHitsDistanceBetween <- distanceToNearest(WindowsWithHits)
### this is good. There are some windows that are right next to each other (window = 0)
# how many?
WindowsWithHitsDistanceBetween_df <- data.frame(WindowsWithHitsDistanceBetween)
dim(WindowsWithHitsDistanceBetween_df[WindowsWithHitsDistanceBetween_df$distance == 0,])
# 329 windows. 
dim(WindowsWithHitsDistanceBetween_df[WindowsWithHitsDistanceBetween_df$distance == 0,])
median(WindowsWithHitsDistanceBetween_df$distance) # MEDIAN distance between windows 210,498 # means half are above that and half are below
