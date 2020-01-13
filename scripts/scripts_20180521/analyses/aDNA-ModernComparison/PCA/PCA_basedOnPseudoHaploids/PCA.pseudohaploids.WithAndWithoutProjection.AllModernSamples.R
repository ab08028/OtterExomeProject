######################### PCA based on pseudohaploids ###################################
require(gdsfmt)
require(SNPRelate)
require(ggplot2)
date="20190612"
#refs=c("mfur","elut")
#mafCutoff=0.2 # should be 1/nInds  # to use the SNPs with ">= maf" only;
#missingCutoff=0.8 # to use the SNPs with "<= missing.rate" only
#conditions=c("highcov","lowcov")
wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/aDNA-ModernComparison/angsd-pseudoHaps/"
plot.dir=""

################################## HIGH COVERAGE #################################
###### biallelic transversions only; filter on maf and missingness #########
cond="highcov"
ref="mfur"
#for(ref in refs){
  ########### get sample info from low coverage bamlist #################
  bamList <- read.table(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/angsd.bamList.mappedto",ref,"fullpaths.HighCovPlusADNAOnly.txt",sep=""),stringsAsFactors = F)
  # this pulls bam info from bamList:
  # THIS IS A USEFUL BIT OF CODE THAT PULLS SAMPLE ID, REFERENCE AND WHETHER WAS DOWNSAMPLED
  # from list of bam names that was used in angsd
  # the order of these individuals will match the rows of the PCs 
  colnames(bamList) <- "bamPath"
  bamList$bamName <- lapply(strsplit(bamList$bamPath,"/"),tail,n=1) # tail pulls out last entry. useful because the different paths have different numbers of dirs so it's not always the same index # but it is the last entry always
  bamList$sampleID <- lapply(strsplit(as.character(bamList$bamName),"\\."),"[",1)
  bamList$plinkIndID <- paste("ind",seq(0,nrow(bamList)-1),sep="") # must be in same order as input and 0based
  bamList$category <- "modern"
  bamList[grep("^A",bamList$bamName),]$category <- "ancient"
  bamList$reference <- ref
  bamList$pop <- lapply(strsplit(as.character(bamList$sampleID),"_"),"[",3)
  # list of samples:
  sample.id = as.character(bamList$sampleID)
  pop1_code = as.character(bamList$pop)
  category_code = as.character(bamList$category)
  
   ########################## read in gds file ###########################
  data.dir=paste(wd,date,"-",cond,"-pseudoHaps/gdsFormat/",sep="")
  inGDSAncient <- paste(data.dir,"angsdOut.mappedTo",ref,".BiallelicTransvOnly.noRefInfo.gds",sep="")
  inGDSModern <- "/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/snps_gds/20181119/snp7/snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.gds" ###### 20191204 experiment: all modern data mapped to mfur
  #genofileModern <- snpgdsOpen(inGDSModern) 
  # make a new gds file that I can alter for Modern:
  genofileModern <- snpgdsOpen(inGDSModern,readonly = F) ## this allows you to modify CAUTION CAUTION (that's why we have a backup)
  #genofileModernNew <- createfn.gds(inGDSModern)
  #### replace the original ID with CHROM_POS (only do once!) 
  # ONLY ONCE ---> add.gdsn(genofileModern$root, "snp.id", paste(read.gdsn(index.gdsn(genofileModern, "snp.chromosome")),"_",read.gdsn(index.gdsn(genofileModern, "snp.position")),sep=""),replace = T) <---- ONLY ONCE #########
  # this is done now.
  ### want to restrict both to those snps ####
  ### adding a second snp ID that is based on location: ##### # only once
  # make modern geno file have snp ids consistent with ancient:
  snpsModern <- data.frame(snpID=read.gdsn(index.gdsn(genofileModern, "snp.id")),snpLocation=paste(read.gdsn(index.gdsn(genofileModern, "snp.chromosome")),"_",read.gdsn(index.gdsn(genofileModern, "snp.position")),sep=""))
  ######### ancient angsd genofile #######
  genofileAncient <- snpgdsOpen(inGDSAncient)

  snpsAncient <- data.frame(snpID=read.gdsn(index.gdsn(genofileAncient, "snp.id")),snpLocation=paste(read.gdsn(index.gdsn(genofileAncient, "snp.chromosome")),"_",read.gdsn(index.gdsn(genofileAncient, "snp.position")),sep=""))
  
  snpsIntersect <- merge(snpsAncient, snpsModern,by="snpLocation",all.x = F,all.y=F,suffixes=c(".ancient",".modern")) # 9673 cool. 
  
  #######
  # LD prune modern data set restricted to the intersected snps:
  snpSetModernLD <- snpgdsLDpruning(genofileModern, snp.id = snpsIntersect$snpLocation,ld.threshold=0.2,autosome.only = F) ## these are modern snps that intersect with the ancient set of snps and then you LD prune them 
  # do LD pruning (r2 = 0.2) to get your snp set e.g. 7,322 markers are selected in total.
  #snpsetAncient <- snpgdsLDpruning(genofileAncient, ld.threshold=0.2,autosome.only = F) # only use snps that are in ancient set
  snpset.id <- unlist(snpSetModernLD) # this is the set of snps to use in PCA: they exist in both ancient and modern and have been LD pruned!
  # 2,557 markers are selected in total.

  #mafCutoff=0.12
  #missingCutoff=0.2
  mafCutoffs=c(0,0.05,0.12,0.2)
  missingCutoffs=c(0,0.2,0.5,0.8)
  for(mafCutoff in mafCutoffs){
    for(missingCutoff in missingCutoffs){
  ############## PCA Projection ########################
  # first do pca only based on modern individuals both with their real names and their Ind0 IDs from the gds file:
  #modernSamples = bamList[bamList$category!="ancient",c("sampleID","plinkIndID")]
  #ancientSamples = bamList[bamList$category=="ancient",c("sampleID","plinkIndID")]
  # do pca based on modern samples only:
  #MODERNONLY_pca <- snpgdsPCA(genofile,snp.id=snpset.id,autosome.only = F,maf = mafCutoff , missing.rate=missingCutoff,sample.id = modernSamples$plinkIndID)
  MODERNONLY_pca <- snpgdsPCA(genofileModern,snp.id=snpset.id,autosome.only = F,maf = mafCutoff , missing.rate=missingCutoff) # all modern samples with PCA just based on the markers taht intersect with ancient dataset
  #MODERNONLY_pca$sample.id <- unlist(modernSamples$sampleID)
  # get SNP loadings:
  modernSNPLoadings <- snpgdsPCASNPLoading(MODERNONLY_pca, genofileModern)
  # project ancient onto modern based on those snp loadings
  #ancientProjectPCA <- snpgdsPCASampLoading(modernSNPLoadings, genofileAncient, sample.id=ancientSamples$plinkIndID)
  ancientProjectPCA <- snpgdsPCASampLoading(modernSNPLoadings, genofileAncient)
  ancientProjectPCA$sample.id <- unlist(bamList$sampleID) # for some reason the last eigenvector is Nas, not sure why
  #modTab <- data.frame(sample.id = MODERNONLY_pca$sample.id, pop1 = factor(pop1_code)[match(MODERNONLY_pca$sample.id, sample.id)], category = factor(category_code)[match(MODERNONLY_pca$sample.id, sample.id)], EV1 = MODERNONLY_pca$eigenvect[,1], EV2 = MODERNONLY_pca$eigenvect[,2],EV3 = MODERNONLY_pca$eigenvect[,3],EV4 = MODERNONLY_pca$eigenvect[,4], stringsAsFactors = FALSE)
  MODERNONLY_pca$pop1 <- unlist(lapply(strsplit(as.character(MODERNONLY_pca$sample.id),"_"),"[",3))
  
  modTab <- data.frame(sample.id = MODERNONLY_pca$sample.id, pop1=MODERNONLY_pca$pop1,category="modern", EV1 = MODERNONLY_pca$eigenvect[,1], EV2 = MODERNONLY_pca$eigenvect[,2],EV3 = MODERNONLY_pca$eigenvect[,3],EV4 = MODERNONLY_pca$eigenvect[,4], stringsAsFactors = FALSE)
  ancProjTab <- data.frame(sample.id = ancientProjectPCA$sample.id, pop1=ancProjTab$pop1,category="angsd", EV1 = ancientProjectPCA$eigenvect[,1], EV2 = ancientProjectPCA$eigenvect[,2], EV3 = ancientProjectPCA$eigenvect[,3], EV4 = ancientProjectPCA$eigenvect[,4],stringsAsFactors = FALSE)
  ancProjTab$pop1 <- unlist(lapply(strsplit(as.character(ancProjTab$sample.id),"_"),"[",3))
  
  # put together:
  modPlusAncProcTab <- rbind(modTab,ancProjTab)
  # plot:
  p2 <- ggplot(modPlusAncProcTab,aes(x=EV1,y=EV2,color=pop1,shape=category))+
    geom_point(size=3)+
    theme_bw()+
    #ylab(paste("PC2 (", pc[2],"%)")) +
    #xlab(paste("PC1 (", pc[1],"%)"))+
    ggtitle(paste("PROJECTED PCA based on ",as.character(length(MODERNONLY_pca$snp.id))," MODERN LD Pruned Biallelic Transversions\nfrom random read sampling (angsd doHaploCall)\nMissingness Cutoff: ",missingCutoff,"; MAF Cutoff: ",mafCutoff,sep=""))+
    theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14)) +
    scale_shape_manual(values=c(1,16)) +
    geom_text(aes(label=sample.id),stat="identity",size=3,position="jitter")
    
  p2
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/PCA/experiment.Projection/sandbox.ProjectOntoAllModernSamples.Pseudohaploids/sandbox.ProjectedSFS.LooksWeird.maf.",mafCutoff,".missing.",missingCutoff,".pdf",sep=""),height=5,width=8)
  
}}
  snpgdsClose(genofileAncient)
  snpgdsClose(genofileModern)
  
  


