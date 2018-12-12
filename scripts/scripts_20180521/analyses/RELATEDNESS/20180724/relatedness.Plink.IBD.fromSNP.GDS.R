#load R packages
library(gdsfmt)
library(SNPRelate)
# guide: https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#principal-component-analysis-pca
# For relatedness analysis, identity-by-descent (IBD) estimation in SNPRelate can be done by either the method of moments (MoM) (Purcell et al., 2007) or maximum likelihood estimation (MLE) (Milligan, 2003; Choi et al., 2009). For both of these methods it is preffered to use a LD pruned SNP set.

calldate=20181119 # date gt's were called in format YYYYMMDD
todaysdate=format(Sys.Date(),format="%Y%m%d")
# file locations:

indir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/snps_gds/",calldate,"/",sep="") # this is where your snp vcf file is and where you will save your gds file (downloaded from hoffman)
plotoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/Relatedness/",calldate,sep="")
fileoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/Relatedness/",calldate,sep="")
dir.create(plotoutdir,recursive = T)
dir.create(fileoutdir,recursive = T)

#open the gds file
genofile <- snpgdsOpen(paste(indir,"/snp_5_passingAllFilters_postMerge_raw_variants.gds",sep=""))

# LD snp pruning: (1min); r2 threshold : 0.2; recommended by SNPRelate tutorial
# https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#ld-based-snp-pruning
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F)
head(snpset)
# Get all selected snp id
snpset.id <- unlist(snpset)
head(snpset.id)

#population information
popmap = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/allSamples.Populations.txt",header=T)
head(popmap)
sample.id = as.character(popmap$Sample)
pop1_code = as.character(popmap$PrimaryPop)
pop2_code = as.character(popmap$SecondaryPop)

# separate by population:
AK.id <- popmap[popmap$PrimaryPop == "Alaska",]$Sample
CA.id <- popmap[popmap$PrimaryPop == "California",]$Sample
KUR.id <- popmap[popmap$PrimaryPop == "Kuril",]$Sample
COM.id <- popmap[popmap$PrimaryPop == "Commander",]$Sample
AL.id <- popmap[popmap$PrimaryPop == "Aleutian",]$Sample

# Estimate IBD coefficients (Estimating IBD Using PLINK method of moments (MoM))
ibd <- snpgdsIBDMoM(genofile, sample.id=AK.id, snp.id=snpset.id,
                    maf=0.05, missing.rate=0.2, num.thread=2,autosome.only = F)

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)
# k0 is the prb of sharing 0 ibd; k1 is prb of sharing 1 ibd (100% ibd?) (expect diagonal? why?)
p1 <- plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
           xlab="k0", ylab="k1", main="YRI samples (MLE)")
lines(c(0,1), c(1,0), col="red", lty=2)
p1
# plot kinship pairs 
p2 <- ggplot(ibd.coeff,aes(x=ID1,y=ID2,fill=kinship))+
  geom_tile()

p2
### how to interpret this?