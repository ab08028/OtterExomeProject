# practice assessing f-stats

# only care about ones that are *negative*

data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/snp7/f3f4Statistics/"
header="snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants"
marker="sepCA-BAJ.exclRelatives"
### standard error comes from jackknifing ### 
f3file <- paste(data.dir,header,".",marker,".f3.out.forR.txt",sep="")
f3 <- read.table(f3file,header=T)

f3$pop1 <- as.character(lapply(strsplit(as.character(f3$pops),";"),"[",1))
f3$sourcePops <- as.character(lapply(strsplit(as.character(f3$pops),";"),"[",2))
#f3$sourcePops <- as.factor(f3$sourcePops)

##### Find possibly significant results ########
# must be significantly NEGATIVE!
totalTests=dim(f3)[1]
alpha=0.05/totalTests # correct for total comparisons
f3_negative <- f3[f3$f3 < 0,]

# note if Z score is positive you would need to take 1-pnorm
# so we insist that it is negative using -abs()
f3_negative$pvalue <- pnorm(-abs(f3_negative$f3_zScore))
f3_negative$significance <- "NS"
f3_negative[f3_negative$pvalue <= alpha,]$significance <- "Significant"
f3_sig <- f3_negative[f3_negative$significance=="Significant",]
f3_sig <- as.data.frame(f3_sig)

require(ggplot2)
ggplot(f3,aes(x=as.factor(pop1),y=as.factor(sourcePops),fill=as.numeric(f3)))+
  geom_tile()+
  theme_bw()
### None of these make ANY sense. This is super weird. Step back and think about it all. Maybe population order got messed up or something. so weird.