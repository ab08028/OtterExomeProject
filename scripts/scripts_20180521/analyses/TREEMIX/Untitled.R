require(GenomicRanges)
require(bootstrap)
binsize=100000 # start with 100kb
###### Idea
# get scaff lengths associated with cds coordinates
# get their full lengths and tile them with gRanges
############### ######################
cds=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/bootStrapRegions/test.out.bed",header=T,sep="\t",strip.white = T, comment.char = "") # converts # to X. in table
############# get cds superfile and make it a Granges object #######
mustelaChrSizes <- (read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/SlidingWindowheterozygosity/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.224.ChrLengths.txt",sep=",",stringsAsFactors = F,strip.white = T))
mustelaChrSizes <- data.frame(t(mustelaChrSizes))
colnames(mustelaChrSizes) <- "record"
mustelaChrSizes$scaff <- lapply(strsplit(as.character(mustelaChrSizes$record),":"),"[",1)
mustelaChrSizes$size <- lapply(strsplit(as.character(mustelaChrSizes$record),":"),"[",2)
mustelaChrSizes_AB <-  mustelaChrSizes[mustelaChrSizes$scaff %in% seqnames(abbaBaba_GRanges),] # 