require(ggplot2)
require(scales)
#### Process pre-seq results

data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/preseq/"
headers=c("A29_Elut_CA_SM_30_SN2_CAP","A30_Elut_CA_SM_35_SN1_CAP","A13_Elut_CA_AN_388_SN1_2CAP_screen")
refs=c("sea_otter_23May2016_bS9RH.deduped.99","Mustela_putorius_furo.MusPutFur1.0.dna.toplevel")
sources=c("DupHist","Bam") # what it's calculated from: bam or dup histogram 
# want to apply this across a table of headers etc
processPreseq <- function(data.dir,header,REF,preseqMethod,source){
  input = read.table(paste(data.dir,header,".",REF,".preseq.",preseqMethod,".from",source,".txt",sep=""),header=T)
  input$header <- header
  input$REF <- REF
  input$preseqMethod <- preseqMethod
  input$source <- source
  return(input)
}
########## TEMPORARY PARAMETERS --- eventually make it a distribution ########
### parameters for extrapoloation: EVENTUALLY GET THESE MY CRIPT THAT PARSES SUMMARY STATS!! ####
desiredReadsPerSample=300000000 # how many sequencing reads I want to sequence
targetSize=57549858
FracBasesOnOrNearTarget=0.94
######### MAKE DATAFRAME OF ALL POSSIBLE COMBINATIONS OF HEADERS/SOURCES/ETC:
# dataframe of all combinations of references, headers, sources; exclude preseq method (do separately); mapply will work through all combos of these parameters to read in all files:

#df <- data.frame(header=headers,source=c("DupHist","Bam"),REF=refs)
# expand grid to get all combinations of parameters:
df_exp <- expand.grid(header=headers,REF=refs,source=sources) #### NOTE! this contains NO DATA -- it just shows the different analysis combos
head(df_exp)
# have to do c_curve and lc_extrap separately

c_curve_allResults <- do.call(rbind,(mapply(processPreseq,data.dir=data.dir,preseqMethod="c_curve",header=df_exp$header,source=df_exp$source,REF=df_exp$REF,SIMPLIFY=FALSE)))
c_curve_allResults$ID <- unlist(lapply(strsplit(as.character(c_curve_allResults$header),"_"),"[",1))

# get rid of row names
rownames(c_curve_allResults)<- NULL
head(c_curve_allResults)
lc_extrap_allResults <- do.call(rbind,(mapply(processPreseq,data.dir=data.dir,preseqMethod="lc_extrap",header=df_exp$header,source=df_exp$source,REF=df_exp$REF,SIMPLIFY=FALSE)))
# get rid of row names
rownames(lc_extrap_allResults)<- NULL
lc_extrap_allResults$ID <- unlist(lapply(strsplit(as.character(lc_extrap_allResults$header),"_"),"[",1))
# plot c-curves 
### WOAH DISCOVERY: facet_grid is way cool! can do 2 variables as row columns
c_curve_plot <- ggplot(c_curve_allResults,aes(x=total_reads,y=distinct_reads,color=source))+
  geom_line()+
  facet_grid(ID~REF,scales="free")+
  #ggtitle(paste("sample: ",header,"\nMapped to: ",REF,"\n",preseqMethod," from ",source,sep=""))+
  theme_bw()+
  scale_y_continuous(labels=comma)+
  scale_x_continuous(labels=comma)+
  geom_abline(slope=1,intercept = 0,linetype="dashed",alpha=0.5)+
  ggtitle("preseq c_curves\ncalculated from Bam or duplicate histogram\nReads mapped to different ref. genomes")
c_curve_plot
ggsave(paste(data.dir,"c_curve_plot.pdf"),c_curve_plot,height=5,width=7)

lc_extrap_plot <- ggplot(lc_extrap_allResults,aes(x=TOTAL_READS,y=EXPECTED_DISTINCT,color=source))+
  geom_line()+
  facet_grid(ID~REF,scales="free")+
  #ggtitle(paste("sample: ",header,"\nMapped to: ",REF,"\n",preseqMethod," from ",source,sep=""))+
  theme_bw()+
  scale_x_continuous(limits=c(0,1e7))+
  scale_y_continuous(limits=c(0,1e7))+
  geom_abline(slope=1,intercept = 0,linetype="dashed",alpha=0.5)+
  ggtitle("preseq lc_curves (extrapolation)\ncalculated from Bam or duplicate histogram\nReads mapped to different ref. genomes")

lc_extrap_plot
ggsave(paste(data.dir,"lc_extrap_plot.zoomedIn.pdf"),lc_extrap_plot,height=5,width=7)

# scale up TOTAL READS to account for the fact that only 3% of reads map to sea otter
# so you have to scale up by 1/0.03 to get the actual total reads
# and add line where 1 lane of Hiseq is ~300M reads 
# scale sequenced reads by actual frac of unique endogenous:
lc_extrap_allResults$FracEndogenousContent <- 0
lc_extrap_allResults[lc_extrap_allResults$ID=="A29",]$FracEndogenousContent <- 0.013 # gotten from paleomix summary stats 
lc_extrap_allResults[lc_extrap_allResults$ID=="A30",]$FracEndogenousContent <- 0.022 # gotten from paleomix summary stats 
lc_extrap_allResults[lc_extrap_allResults$ID=="A13",]$FracEndogenousContent <- 0.028 # gotten from paleomix summary stats 
# and specific info about the read length:
lc_extrap_allResults$AvgReadLen <- 0
lc_extrap_allResults[lc_extrap_allResults$ID=="A29",]$AvgReadLen <- 119 # gotten from paleomix summary stats 
lc_extrap_allResults[lc_extrap_allResults$ID=="A30",]$AvgReadLen <- 130 # gotten from paleomix summary stats 
lc_extrap_allResults[lc_extrap_allResults$ID=="A13",]$AvgReadLen <- 129 # gotten from paleomix summary stats 
# so for each one can get the probable perc Unique 
lc_extrap_plot2 <- ggplot(lc_extrap_allResults,aes(x=TOTAL_READS/FracEndogenousContent,y=EXPECTED_DISTINCT,color=source))+
  geom_line()+
  # add shadow of 95CI
  #geom_abline(slope=1,intercept = 0,linetype="dashed",alpha=0.5)+
  geom_ribbon(aes(ymin=LOWER_0.95CI,ymax=UPPER_0.95CI),alpha=0.2,linetype=0)+
  #ggtitle(paste("sample: ",header,"\nMapped to: ",REF,"\n",preseqMethod," from ",source,sep=""))+
  theme_bw()+
  #scale_y_log10()+
  #scale_x_log10()+
  scale_x_continuous(limits=c(0,5e8))+
  scale_y_continuous(limits=c(0,3e7))+
  #geom_abline(slope=1,intercept = 0,linetype="dashed",alpha=0.5)+
  #geom_vline(xintercept = 31000000)+
  geom_vline(xintercept = 300000000,linetype="dashed")+
  #geom_hline(yintercept = lc_extrap_allResults[lc_extrap_allResults$TOTAL_READS==1000000,]$EXPECTED_DISTINCT,aes(color=source,linetype=REF))+
  #geom_hline(yintercept = lc_extrap_allResults[lc_extrap_allResults$TOTAL_READS==9e6,]$EXPECTED_DISTINCT,aes(color=source,linetype=REF))+
  facet_grid(ID~REF,scales="free")+
  xlab(paste("Scaled up number of reads to sequence\ntotal reads predicted to map to otter scaled by unique endogenous fraction",sep=""))+
  ggtitle("preseq lc_curves (extrapolation)\ncalculated from Bam or duplicate histogram\nReads mapped to different ref. genomes\nWith Scaled up read count and dash for planned sequencing effort")


lc_extrap_plot2

ggsave(paste(data.dir,"lc_extrap_plot.zoomedInFurther.ScaledUpReadCount.pdf"),lc_extrap_plot2,height=5,width=7)

###################### make plots for medgenome: sea otter and dup histogram only ##########
### WOAH DISCOVERY: facet_grid is way cool! can do 2 variables as row columns
# sea otter ref only; source is DupHist only 
#### c curve ####
c_curve_plot_simple <- ggplot(c_curve_allResults[c_curve_allResults$REF=="sea_otter_23May2016_bS9RH.deduped.99" & c_curve_allResults$source=="DupHist",],aes(x=total_reads,y=distinct_reads))+
  geom_line(color="tomato3")+
  facet_grid(~ID)+
  #ggtitle(paste("sample: ",header,"\nMapped to: ",REF,"\n",preseqMethod," from ",source,sep=""))+
  theme_bw()+
  geom_abline(slope=1,intercept = 0,linetype="dashed",alpha=0.5)+
  ggtitle("preseq c_curves, sea otter ref\ncalculated from duplicate histogram")+
  xlab("Total reads mapped to sea otter")
c_curve_plot_simple
ggsave(paste(data.dir,"c_curve_plot.elutOnly.simplified.pdf"),c_curve_plot_simple,height=5,width=10)
#### extrapolation  ####
lc_extrap_plot_simple <- ggplot(lc_extrap_allResults[lc_extrap_allResults$REF=="sea_otter_23May2016_bS9RH.deduped.99" & lc_extrap_allResults$source=="DupHist",],aes(x=TOTAL_READS,y=EXPECTED_DISTINCT))+
  geom_line(color="darkorange")+
  facet_grid(~ID,scales="free")+
  #ggtitle(paste("sample: ",header,"\nMapped to: ",REF,"\n",preseqMethod," from ",source,sep=""))+
  theme_bw()+
  scale_x_continuous(limits=c(0,1e7))+
  scale_y_continuous(limits=c(0,1e7))+
  geom_abline(slope=1,intercept = 0,linetype="dashed",alpha=0.5)+
  ggtitle("preseq lc_curves (extrapolation)\ncalculated from duplicate histogram\nReads mapped to sea otter")+
  xlab("Total reads mapping to sea otter")

lc_extrap_plot_simple
ggsave(paste(data.dir,"lc_extrap_plot.zoomedIn.ElutOnly.pdf"),lc_extrap_plot_simple,height=5,width=10)

#### scale up number of reads by endogenous fraction: ####  
lc_extrap_plot2_simple <- ggplot(lc_extrap_allResults[lc_extrap_allResults$REF=="sea_otter_23May2016_bS9RH.deduped.99" & lc_extrap_allResults$source=="DupHist",],aes(x=TOTAL_READS/FracEndogenousContent,y=EXPECTED_DISTINCT))+
  geom_line(color="darkorange")+
  # add shadow of 95CI
  #geom_abline(slope=1,intercept = 0,linetype="dashed",alpha=0.5)+
  geom_ribbon(aes(ymin=LOWER_0.95CI,ymax=UPPER_0.95CI),alpha=0.2,linetype=0)+
  #ggtitle(paste("sample: ",header,"\nMapped to: ",REF,"\n",preseqMethod," from ",source,sep=""))+
  theme_bw()+
  #scale_y_log10()+
  #scale_x_log10()+
  scale_x_continuous(limits=c(0,5e8))+
  scale_y_continuous(limits=c(0,2e7))+
  #geom_abline(slope=1,intercept = 0,linetype="dashed",alpha=0.5)+
  #geom_vline(xintercept = 31000000)+
  geom_vline(xintercept = 450000000,linetype="dashed")+
  #geom_hline(yintercept = lc_extrap_allResults[lc_extrap_allResults$TOTAL_READS==1000000,]$EXPECTED_DISTINCT,aes(color=source,linetype=REF))+
  #geom_hline(yintercept = lc_extrap_allResults[lc_extrap_allResults$TOTAL_READS==9e6,]$EXPECTED_DISTINCT,aes(color=source,linetype=REF))+
  facet_grid(~ID,scales="free")+
  xlab(paste("Scaled up number of reads to sequence\n(total reads predicted to map to otter scaled by unique endogenous fraction)",sep=""))+
  ggtitle("preseq lc_curves (extrapolation)\ncalculated from duplicate histogram\nReads mapped to sea otter\nWith Scaled up read count and dash for planned sequencing effort")

lc_extrap_plot2_simple

ggsave(paste(data.dir,"lc_extrap_plot.zoomedInFurther.ScaledUpReadCount.MoreSequencing.ElutOnly.pdf"),lc_extrap_plot2_simple,height=5,width=10)

############ plot coverage ##############

# multiply the expected number of distinct reads by the *read len* on target frac and then divide by total bases
lc_extrap_plot4_simple <- ggplot(lc_extrap_allResults[lc_extrap_allResults$REF=="sea_otter_23May2016_bS9RH.deduped.99" & lc_extrap_allResults$source=="DupHist",],aes(x=TOTAL_READS/FracEndogenousContent,y=(EXPECTED_DISTINCT*AvgReadLen*FracBasesOnOrNearTarget)/targetSize))+
  geom_line(color="darkorange")+
  theme_bw()+
  scale_x_continuous(limits=c(0,6e8))+
  scale_y_continuous(limits=c(0,45))+
  geom_vline(xintercept = 450000000,linetype="dashed")+
  facet_grid(~ID,scales="free")+
  ylab("Approximate Coverage\n(expected reads * avg read length * near-target frac)/(target size)")+
  xlab("Total reads sequenced")

lc_extrap_plot4_simple

ggsave(paste(data.dir,"lc_extrap_plot.ExpCoverage.MoreSequencing.ElutOnly.pdf"),lc_extrap_plot4_simple,height=5,width=10)



###################### back of envelope calculations #################

# go from ~1M total reads mapping to otter which yields ~1M distinct reads (mean 2x coverage of targets)
unique(lc_extrap_allResults[lc_extrap_allResults$TOTAL_READS==1000000 & lc_extrap_allResults$source=="Bam" & lc_extrap_allResults$header=="A13_Elut_CA_AN_388_SN1_2CAP_screen" & lc_extrap_allResults$REF=="sea_otter_23May2016_bS9RH.deduped.99",]$EXPECTED_DISTINCT) # 963,334.5 distinct reads
# to 9M reads mapping to otter which yields 8M distinct reads (should hopefully be mean 18x coverage of targets)
unique(lc_extrap_allResults[lc_extrap_allResults$TOTAL_READS==9e6 & lc_extrap_allResults$source=="Bam" & lc_extrap_allResults$header=="A13_Elut_CA_AN_388_SN1_2CAP_screen" & lc_extrap_allResults$REF=="sea_otter_23May2016_bS9RH.deduped.99",]$EXPECTED_DISTINCT) # 7,266,957 distinct reads only (out of 300M)

# so with 7M reads, what sort of coverage would that give me? 60Mb of capture. 300M reads * 300bp per read pair (probably less than that honestly -- maybe set as 200?)
(unique(lc_extrap_allResults[lc_extrap_allResults$TOTAL_READS==9e6 & lc_extrap_allResults$source=="Bam" & lc_extrap_allResults$header=="A13_Elut_CA_AN_388_SN1_2CAP_screen" & lc_extrap_allResults$REF=="sea_otter_23May2016_bS9RH.deduped.99",]$EXPECTED_DISTINCT) * 200 ) / 60000000  # 20-30x coverage of targets . this would be great. let's do it.

# Kind of ridiculous to use a whole sequencing lane just to get an exome. but it's an ANCIENT exome! and at least won't oversaturate based on the curves.

# Plot expected coverage 
#desiredBP = probablePercUnique*avgBasesPerRead*desiredReadsPerSample # desired base pairs
#targetsize=57549858
#expCoverage=desiredBP/targetsize
# need to also scale by the lc_extrapolation factor for each individual 
# 

