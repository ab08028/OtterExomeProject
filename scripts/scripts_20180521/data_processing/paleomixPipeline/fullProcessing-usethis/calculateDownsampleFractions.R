require(dplyr)
numReps = 1 # how many times you want to resample the bam file (with replacement)
### Want to make a script that calculates the downsampling fraction for my modern samples
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/plmx_mappingStats/aDNA.modern.Comparison.SummaryStats-9samples-usethis/"
script.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/paleomixPipeline/fullProcessing-usethis/"
todaysdate=format(Sys.Date(),"%Y%m%d")
moderns=c("140_Elut_CA_403","141_Elut_CA_419","116_Elut_CA_307","126_Elut_AK_AF3394","55_Elut_AK_AF3736","129_Elut_AK_AL4660")
ancients=c("A30_Elut_CA_SM_35_SN1_CAP","A29_Elut_CA_SM_30_SN2_CAP","A13_Elut_CA_AN_388_SN1_2CAP_screen")
modern_ancient_pairs = 
elutref="sea_otter_23May2016_bS9RH.deduped.99"
mfurref="Mustela_putorius_furo.MusPutFur1.0.dna.toplevel"
# stats from plmx
stats <- read.table(paste(data.dir,"temp.incomplete.plmx.aDNA.summary.stats.fullProcessing.ancient.modern.comparison.20190422.txt",sep=""),header=T)
head(stats)

# so we want to summarize a few things
# we want to group by sample and reference and summarize the mean of hits_unique among the ancient samples; and also match 1:1
# should I resample a lot? that's an interesting idea. Resample 100 times and get a range? It's an option.

unique_hits <- stats[stats$statistic=="hits_unique",]

# label by ancient / modern
unique_hits$label <- NA
unique_hits[unique_hits$sample %in% ancients,]$label <- "ancient"
unique_hits[unique_hits$sample %in% moderns,]$label <- "modern"


# Want mean of Ancient reads separated by reference
unique_hits <- unique_hits %>%
  group_by(reference) %>% # group by reference genome
  mutate(meanAncient = mean(value[label=="ancient"])) %>% # calculates mean of ancient samples split by reference and adds it as a column (matches on reference)
  mutate(downSampleFrac_meanAncient = meanAncient/value )  # divide the ancient reads by the number of reads in each sample to get the downsample frac. (Maybe will downsample ancient samples to this number as well?) or maybe not. This will be the -s for the modern samples 

modern_unique_hits <- unique_hits[unique_hits$label=="modern",]
ancient_unique_hits <- unique_hits[unique_hits$label=="ancient",]

# Want to write a script:
# adding rep to fraction because rep will be the seed, frac will be the frac. the integer part could be anything (this is a funny samtools oddity)
# so 3.045 would get a fraction of 0.045 and 3 is the seed 

########### sample multiple times if desired (each replicate is a different seed)
# start your outfile
sink(paste(script.dir,"downsampleModernSamples.",todaysdate,".sh",sep=""))
cat("# Downsample modern bam files to match mean unique reads mapped to sea otter and ferret genomes in the best 3 ancient samples\n")
cat("# Mean reads mapped to sea otter: ",unique(unique_hits[unique_hits$reference==elutref,]$meanAncient))
cat("\n")
cat("# Mean reads mapped to ferret: ",unique(unique_hits[unique_hits$reference==mfurref,]$meanAncient))
cat("\n")
cat("wd=/u/flashscratch/a/ab08028/captures/paleomix/fullProcessing/\n")
cat("downsampledir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/downsampledBams/downsample_allEven/\n\n") 
# go through replicates if you want to resample
for(rep in seq(1,numReps)){
 # making them all even 
  cat("echo \"downsample to equal mean ancient\"\n\n" )
  cat(paste("samtools view -s ",rep+round(modern_unique_hits$downSampleFrac_meanAncient,4)," -b $wd/",modern_unique_hits$sample,"/",modern_unique_hits$sample,".",modern_unique_hits$reference,".bam > $downsampledir/",modern_unique_hits$sample,".",modern_unique_hits$reference,".downsamp.",rep+round(modern_unique_hits$downSampleFrac_meanAncient,4),".bam\n",sep=""))
  cat("\n")
  # and then count downsampled reads:
  # make a header
  cat("echo \"count downsampled reads\"\n\n" )
  cat("echo downsampledReadCounts > $downsampledir/downsampledReadCounts.txt\n")
  cat(paste("samtools flagstat ","$downsampledir/",modern_unique_hits$sample,".",modern_unique_hits$reference,".downsamp.",rep+round(modern_unique_hits$downSampleFrac_meanAncient,4),".bam | head -n1 | awk '{cat $1}' >> $downsampledir/downsampledReadCounts.txt\n",sep=""))
}
sink()

############################### Downsample to 1:1 match the ancient samples ######################
# need to make pairs where a CA and AK modern sample is matched to a specific ancient sample
head(unique_hits)

# create sets:
set1=c("55_Elut_AK_AF3736","116_Elut_CA_307") # "A13_Elut_CA_AN_388_SN1_2CAP_screen"
set2=c("126_Elut_AK_AF3394","141_Elut_CA_419") # "A29_Elut_CA_SM_30_SN2_CAP"
set3=c("129_Elut_AK_AL4660","140_Elut_CA_403") # "A30_Elut_CA_SM_35_SN1_CAP"

modern_unique_hits$downSampleFrac_matchSpecificAncient <- 1
modern_unique_hits$ancientToMatch <- "A13_Elut_CA_AN_388_SN1_2CAP_screen"

#### 
