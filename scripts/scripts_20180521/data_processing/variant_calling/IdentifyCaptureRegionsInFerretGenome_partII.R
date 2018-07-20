# these are blast results from this command:
#blastn -num_threads 10 -evalue 1e-10 -max_target_seqs 5 -query $wd/neutralPromoterCaptureRegions_noExon.seaOtterSequence.nonStranded.20180717.fasta -db $ferretDB -outfmt 6 > $wd/blastResults/blast.seaOtterNeutralPromoterCaptureRegions_vs_ferretGenome.top5.out
# they have been prefiltered to only exclude hits of min. 500bp in length (since each query is ~1kb)
blast <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/captureRegions_Ferret/blastResults/blast.seaOtterNeutralPromoterCaptureRegions_vs_ferretGenome.min500bp.out")
head(blast)
colnames(blast) <- c("query","scaffold","perc_ident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
View(blast)
dups <- blast[duplicated(blast$query),]
View(dups)
####### First choose the non duplicated ones
####### Then for the duplicated ones, only choose ones with high length
require(dplyr)
distinctBlast <- blast[!duplicated(blast$query) & !duplicated(blast$query,fromLast = T),]

dim(distinctBlast)
nonDistinctBlast <- blast[duplicated(blast$query,fromLast = T) | duplicated(blast$query,fromLast = F),]
dim(nonDistinctBlast)
nonDistinctBlast.900 <- nonDistinctBlast[nonDistinctBlast$length>900,]
dim(nonDistinctBlast.900) 
length(unique(nonDistinctBlast.900$query)) # 187/327 entries are present. Good enough? try to round up a few more? some are still duplicated though.

######## STUCK HERE #############
# look at tally for teach one
dupTally <- nonDistinctBlast %>%
  group_by(query) %>%
  tally()
# these add up to the total
# Now need to decide between the nonDistinct ones : 
### want to do something like: see if the two halves are near each other, or if one is much bigger than the other ** stuck here ** start here tomorrow! 