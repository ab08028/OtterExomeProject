require(dplyr)
numReps = 1 # how many times you want to resample the bam file (with replacement)
### Want to make a script that calculates the downsampling fraction for my modern samples
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/plmx_mappingStats/aDNA.modern.Comparison.SummaryStats-9samples-usethis/"
script.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/paleomixPipeline/fullProcessing-usethis/"
todaysdate=format(Sys.Date(),"%Y%m%d")
moderns=c("140_Elut_CA_403","141_Elut_CA_419","116_Elut_CA_307","126_Elut_AK_AF3394","55_Elut_AK_AF3736","129_Elut_AK_AL4660")
ancients=c("A30_Elut_CA_SM_35_SN1_CAP","A29_Elut_CA_SM_30_SN2_CAP","A13_Elut_CA_AN_388_SN1_2CAP_screen")
# setting pairs manually: a CA and an AK for each ancient sample 
CA_ancientPairs = list(c("140_Elut_CA_403","A30_Elut_CA_SM_35_SN1_CAP"),c("141_Elut_CA_419","A29_Elut_CA_SM_30_SN2_CAP"),c("116_Elut_CA_307","A13_Elut_CA_AN_388_SN1_2CAP_screen")) 

AK_ancientPairs= list(c("126_Elut_AK_AF3394","A30_Elut_CA_SM_35_SN1_CAP"),c("55_Elut_AK_AF3736","A29_Elut_CA_SM_30_SN2_CAP"),c("129_Elut_AK_AL4660","A13_Elut_CA_AN_388_SN1_2CAP_screen")) 
all_modern_ancientPairs = c(CA_ancientPairs,AK_ancientPairs)

elutref="sea_otter_23May2016_bS9RH.deduped.99"
mfurref="Mustela_putorius_furo.MusPutFur1.0.dna.toplevel"
# stats from plmx
stats <- read.table(paste(data.dir,"plmx.aDNA.summary.stats.fullProcessing.ancient.modern.comparison.20190501.txt",sep=""),header=T)
head(stats)

# so we want to summarize a few things
# we want to group by sample and reference and summarize the mean of hits_unique among the ancient samples; and also match 1:1
# should I resample a lot? that's an interesting idea. Resample 100 times and get a range? It's an option.

unique_hits <- stats[stats$statistic=="hits_unique",]

# label by ancient / modern
unique_hits$label <- NA
unique_hits[unique_hits$sample %in% ancients,]$label <- "ancient"
unique_hits[unique_hits$sample %in% moderns,]$label <- "modern"

# want to plot them:
p0 <- ggplot(unique_hits,aes(x=sample,y=value,fill=reference))+
  geom_bar(stat="identity",position="dodge")+
  coord_flip()+
  theme_bw()+
  ggtitle("Unique reads mapped to reference")
p0
ggsave(paste(data.dir,"/uniqueHitsPerSample.fullProcessingModernAncient.pdf",sep=""),device="pdf",height=7,width=14)
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
############## skipping downsamplingt to the mean value, in favor of doing matched pairs with the each ancient sample #################### 
########### sample multiple times if desired (each replicate is a different seed) ###
# start your outfile
# sink(paste(script.dir,"downsampleModernSamples.ToMeanAncientReads.",todaysdate,".sh",sep=""))
# cat("# Downsample modern bam files to match mean unique reads mapped to sea otter and ferret genomes in the best 3 ancient samples\n")
# cat("# Mean reads mapped to sea otter: ",unique(unique_hits[unique_hits$reference==elutref,]$meanAncient))
# cat("\n")
# cat("# Mean reads mapped to ferret: ",unique(unique_hits[unique_hits$reference==mfurref,]$meanAncient))
# cat("\n")
# cat("wd=/u/flashscratch/a/ab08028/captures/paleomix/fullProcessing/\n")
# cat("downsampledir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/downsampledBams/downsample_allEven/\n\n") 
# # go through replicates if you want to resample
# for(rep in seq(1,numReps)){
#  # making them all even 
#   cat("echo \"downsample to equal mean ancient\"\n\n" )
#   cat(paste("samtools view -s ",rep+round(modern_unique_hits$downSampleFrac_meanAncient,4)," -b $wd/",modern_unique_hits$sample,"/",modern_unique_hits$sample,".",modern_unique_hits$reference,".bam > $downsampledir/",modern_unique_hits$sample,".",modern_unique_hits$reference,".downsamp.",rep+round(modern_unique_hits$downSampleFrac_meanAncient,4),".bam\n",sep=""))
#   cat("\n")
#   # and then count downsampled reads:
#   # make a header
#   cat("echo \"count downsampled reads\"\n\n" )
#   cat("echo downsampledReadCounts > $downsampledir/downsampledReadCounts.txt\n")
#   cat(paste("samtools flagstat ","$downsampledir/",modern_unique_hits$sample,".",modern_unique_hits$reference,".downsamp.",rep+round(modern_unique_hits$downSampleFrac_meanAncient,4),".bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt\n",sep=""))
# }
# sink()

############################### Make a script to downsample to 1:1 match the ancient samples ######################
# need to make pairs where a CA and AK modern sample is matched to a specific ancient sample
# now want to match number of mapped unique reads between pairs
sink(paste(script.dir,"downsampleModernSamples.MatchPairs.",todaysdate,".sh",sep=""))
cat("source /u/local/Modules/default/init/modules.sh\n")
cat("module load samtools\n")
cat("# Downsample modern bam files to match the number unique reads mapped to sea otter and ferret genomes in the best 3 ancient samples (in pairs -- each ancient sample matches 1 CA and 1 AK sample)\n")
cat("# note: the -s value has the replicate number as the integer which is the seed and the decimal part is the fraction. So 1.006 is replicate 1, fraction 0.006")
cat("\n")
cat("wd=/u/flashscratch/a/ab08028/captures/paleomix/fullProcessing/\n")
cat("downsampledir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/downsampledBams/downsample_Pairs/\n\n") 
# go through replicates if you want to resample
cat("echo downsampledReadCounts > $downsampledir/downsampledReadCounts.txt\n")
for(rep in seq(1,numReps)){
  for(ref in unlist(unique(unique_hits$reference))){
    for(pair in all_modern_ancientPairs){
      # first of pair: first entry is modern, second is ancient:
      modernID=unlist(pair)[1]
      ancID=unlist(pair)[2]
      modernCount=unique_hits[unique_hits$sample==modernID & unique_hits$reference==ref,]$value
      ancientCount=unique_hits[unique_hits$sample==ancID & unique_hits$reference==ref,]$value
      # check that mod > anc
      modernCount>=ancientCount
      # get ratio to get downsample rate
      downSampleFrac = ancientCount/modernCount
      cat("# downsample ",modernID," to equal ancient sample ",ancID,"\n" ,sep="")
      cat("# ",modernID," (modern) starting reads: ",modernCount,"\n",sep="")
      cat("# ",ancID," (ancient) starting reads: ",ancientCount,"\"\n\n",sep="")
      cat(paste("samtools view -s ",rep+round(downSampleFrac,4)," -b $wd/",modernID,"/",modernID,".",ref,".bam > $downsampledir/",modernID,".",ref,".downsamp.rep.",rep,".bam\n",sep=""))
      cat("# Count the resulting reads to make sure it downsampled properly\n")
      cat("echo \"",modernID,"\" >> $downsampledir/downsampledReadCounts.txt\n",sep="")
      cat("samtools flagstat ","$downsampledir/",modernID,".",ref,".downsamp.rep.",rep,".bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt\n",sep="")
      cat("\n\n")
    }
  }
}

sink()
