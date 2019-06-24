######################### PCA based on pseudohaploids ###################################
require(gdsfmt)
require(SNPRelate)
require(ggplot2)
date="20190612"
refs=c("mfur","elut")
#mafCutoff=0.3 # should be 1/nInds 
#missingCutoff=0
#conditions=c("highcov","lowcov")
wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/aDNA-ModernComparison/angsd-pseudoHaps/"
plot.dir=""

################################## LOW COVERAGE #################################
###### biallelic transversions only; filter on maf and missingness #########
cond="lowcov"
ref="mfur"
#for(ref in refs){
  ########### get sample info from low coverage bamlist #################
  bamList <- read.table(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/angsd.bamList.LowCoverageOnly.mappedto",ref,"fullpaths.txt",sep=""),stringsAsFactors = F)
  # this pulls bam info from bamList:
  # THIS IS A USEFUL BIT OF CODE THAT PULLS SAMPLE ID, REFERENCE AND WHETHER WAS DOWNSAMPLED
  # from list of bam names that was used in angsd
  # the order of these individuals will match the rows of the PCs 
  colnames(bamList) <- "bamPath"
  bamList$bamName <- lapply(strsplit(bamList$bamPath,"/"),tail,n=1) # tail pulls out last entry. useful because the different paths have different numbers of dirs so it's not always the same index # but it is the last entry always
  bamList$sampleID <- lapply(strsplit(as.character(bamList$bamName),"\\."),"[",1)
  bamList$plinkIndID <- paste("ind",seq(0,nrow(bamList)-1),sep="") # must be in same order as input and 0based
  bamList$category <- "modern"
  bamList[grep("downsamp",bamList$bamName),]$category <- "modern-downsampled"
  bamList[grep("^A",bamList$bamName),]$category <- "ancient"
  bamList$reference <- ref
  bamList$pop <- lapply(strsplit(as.character(bamList$sampleID),"_"),"[",3)
  # list of samples:
  sample.id = as.character(bamList$sampleID)
  pop1_code = as.character(bamList$pop)
  category_code = as.character(bamList$category)
  ########################## read in gds file ###########################
  data.dir=paste(wd,date,"-",cond,"-pseudoHaps/gdsFormat/",sep="")
  inGDS <- paste(data.dir,"angsdOut.mappedTo",ref,".BiallelicTransvOnly.noRefInfo.gds",sep="")
  genofile <- snpgdsOpen(inGDS) 
  # do LD pruning (r2 = 0.2) to get your snp set e.g. 7,322 markers are selected in total.
  snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F)
  snpset.id <- unlist(snpset)
  ############## Regular PCA ########################
  mafCutoff=0 # should be 1/nInds 
  missingCutoff=0.8
  pca <- snpgdsPCA(genofile,snp.id=snpset.id,autosome.only = F,maf = mafCutoff , missing.rate=missingCutoff)
  # update sampleids: 
  pca$sample.id <- unlist(bamList$sampleID)
  # get pc percents:
  pc.percent <- pca$varprop*100
  pc = head(round(pc.percent, 2))
  pc
  # make a dataframe
  tab1a <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], category = factor(category_code)[match(pca$sample.id, sample.id)], EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], EV3=pca$eigenvect[,3],stringsAsFactors = FALSE)
  # plot:
  p1a <- ggplot(tab1a,aes(x=EV1,y=EV2,shape=category,color=pop1))+
    geom_point(size=3)+
    theme_bw()+
    ylab(paste("PC2 (", pc[2],"%)")) +
    xlab(paste("PC1 (", pc[1],"%)"))+
    ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned Biallelic Transversions\nfrom random read sampling (angsd doHaploCall)\nMissingness Cutoff: ",missingCutoff,"; MAF Cutoff: ",mafCutoff,sep=""))+
    theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14)) +
    scale_shape_manual(values=c(1,16)) 
  p1a
  # pc 2 vs 3:
    p1b <- ggplot(tab1a,aes(x=EV2,y=EV3,shape=category,color=pop1))+
    geom_point(size=3)+
    theme_bw()+
    ylab(paste("PC3 (", pc[3],"%)")) +
    xlab(paste("PC2 (", pc[2],"%)"))+
    ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned Biallelic Transversions\nfrom random read sampling (angsd doHaploCall)\nMissingness Cutoff: ",missingCutoff,"; MAF Cutoff: ",mafCutoff,sep=""))+
    theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14)) +
    scale_shape_manual(values=c(1,16)) 
  p1b
  # pc 1 vs 3:
  p1c <- ggplot(tab1a,aes(x=EV1,y=EV3,shape=category,color=pop1))+
    geom_point(size=3)+
    theme_bw()+
    ylab(paste("PC3 (", pc[3],"%)")) +
    xlab(paste("PC1 (", pc[1],"%)"))+
    ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned Biallelic Transversions\nfrom random read sampling (angsd doHaploCall)\nMissingness Cutoff: ",missingCutoff,"; MAF Cutoff: ",mafCutoff,sep=""))+
    theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14)) +
    scale_shape_manual(values=c(1,16)) 
  p1c
  
  
  ############## PCA Projection ########################
  # first do pca only based on modern individuals both with their real names and their Ind0 IDs from the gds file:
  modernSamples = bamList[bamList$category!="ancient",c("sampleID","plinkIndID")]
  ancientSamples = bamList[bamList$category=="ancient",c("sampleID","plinkIndID")]
  # do pca based on modern samples only:
  MODERNONLY_pca <- snpgdsPCA(genofile,snp.id=snpset.id,autosome.only = F,maf = mafCutoff , missing.rate=missingCutoff,sample.id = modernSamples$plinkIndID)
  # get SNP loadings:
  modernSNPLoadings <- snpgdsPCASNPLoading(MODERNONLY_pca, genofile)
  # project ancient onto modern based on those snp loadings
  ancientProjectPCA <- snpgdsPCASampLoading(modernSNPLoadings, genofile, sample.id=ancientSamples$plinkIndID)
  ancientProjectPCA$sample.id <- unlist(ancientSamples$sampleID)
  MODERNONLY_pca$sample.id <- unlist(modernSamples$sampleID)
  
  #### put together ########
  modTab <- data.frame(sample.id = MODERNONLY_pca$sample.id, pop1 = factor(pop1_code)[match(MODERNONLY_pca$sample.id, sample.id)], category = factor(category_code)[match(MODERNONLY_pca$sample.id, sample.id)], EV1 = MODERNONLY_pca$eigenvect[,1], EV2 = MODERNONLY_pca$eigenvect[,2], stringsAsFactors = FALSE)
  ancProjTab <- data.frame(sample.id = ancientProjectPCA$sample.id, pop1 = factor(pop1_code)[match(ancientProjectPCA$sample.id, sample.id)], category = factor(category_code)[match(ancientProjectPCA$sample.id, sample.id)], EV1 = ancientProjectPCA$eigenvect[,1], EV2 = ancientProjectPCA$eigenvect[,2], stringsAsFactors = FALSE)
  # put together:
  modPlusAncProcTab <- rbind(modTab,ancProjTab)
  # plot:
  p2 <- ggplot(modPlusAncProcTab,aes(x=EV1,y=EV2,shape=category,color=pop1))+
    geom_point(size=3)+
    theme_bw()+
    ylab(paste("PC2 (", pc[2],"%)")) +
    xlab(paste("PC1 (", pc[1],"%)"))+
    ggtitle(paste("PROJECTED PCA based on ",as.character(length(MODERNONLY_pca$snp.id))," MODERN LD Pruned Biallelic Transversions\nfrom random read sampling (angsd doHaploCall)\nMissingness Cutoff: ",missingCutoff,"; MAF Cutoff: ",mafCutoff,sep=""))+
    theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14)) +
    scale_shape_manual(values=c(1,16)) 
  p2
  ggsave()
  #}
