########### Here I just want to plot the SNPs ###########
# installation; 
### note: these files contain monomorphic 1/1 sites, but those are automaticall
# excluded in the HMM method of ZooRoh. but you should remove them in other methods (Plink etc.)
# note that hbd refers to "homozygous by descent"
#install.packages("RZooRoH")
# install.packages("R.utils") # for reading gzip
require(RZooRoH)
require(R.utils)
require(ggplot2)
require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[1],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5])
modelName="mix10R"
wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/zooROH/"
### converted from vcf > ped > oxford Gen , which is the "GP" format in ZooROH
pops=c("CA","AK","AL","COM","KUR") # alreadydid AK 

########### get scaff sizes and exclude <2Mb ###########
mustelaChrSizes <- (read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/SlidingWindowheterozygosity/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.224.ChrLengths.txt",sep=",",stringsAsFactors = F,strip.white = T))
mustelaChrSizes <- data.frame(t(mustelaChrSizes))
colnames(mustelaChrSizes) <- "record"
mustelaChrSizes$scaff <- unlist(lapply(strsplit(as.character(mustelaChrSizes$record),":"),"[",1))
mustelaChrSizes$size <- unlist(lapply(strsplit(as.character(mustelaChrSizes$record),":"),"[",2))
mustelaChrSizes_gt5Mb <- mustelaChrSizes[as.numeric(mustelaChrSizes$size) > 5e6,]
################################### plots #################
allResults <- data.frame()
for(pop in pops){
  #pop="AK"
  print(paste("Starting pop",pop))
  data.dir=paste(wd,pop,"/",sep="")
  out.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PlottingSNPsPerScaffold/",pop,"/",sep="")
  dir.create(out.dir,showWarnings = F)
  ##################### read in the population specific data and sample IDs #############
  input <- read.table(paste(data.dir,"/",pop,".Oxford.maf.001.gen.gz",sep=""))
  colnames(input) <- c("Chr","ID","Pos")
  # this has no mono sites

  # this will include monomorphic -- want to remove somehow. 
  input_gt5Mb <- input[input$Chr %in% mustelaChrSizes_gt5Mb$scaff,] # restrict to chr htat are >2Mb
  input_gt5Mb$pop <- pop
  p1 <- ggplot(input_gt5Mb,aes(x=Pos/1e6,y=Chr,color=Chr))+
    geom_point(size=0.01)+
    theme_bw()+
    ggtitle(paste(pop,": ",dim(input)[1], " SNPs (not LD pruned; only scaffs >5Mb shown)",sep=""))+
    theme(axis.text.y=element_blank(),legend.position = 'none')+
    ylab("Scaffold")+
    xlab("Pos (Mb)")+
    scale_color_manual(values=rep(c("tomato","dodgerblue"),length(unique(input_gt5Mb$Chr))+2/2))
  p1
  ggsave(paste(out.dir,pop,".snpsPerScaffold.NoMonomorphic.notLDPruned.pdf",sep=""),p1,height = 5,width=7)
  
  # add together
  allResults <- rbind(allResults,input_gt5Mb[,c("Chr","Pos","pop")])
}

############# Plot just first scaffold (biggest) for comparison to ROHs #################
# order factors
allResults$pop <- factor(allResults$pop, levels=c("CA","AK","AL","COM","KUR"))
p2 <- ggplot(allResults[allResults$Chr=="GL896898.1",],aes(x=Pos/1e06,y=pop,color=pop))+
  geom_point(size=0.1)+
  xlab("Position (Mb) on scaffold GL896898.1")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_manual(values=unlist(colors))
p2
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PlottingSNPsPerScaffold/Scaffold1.EachPopulation.Snps.NoMono.pdf",sep=""),height=4,width=8)
