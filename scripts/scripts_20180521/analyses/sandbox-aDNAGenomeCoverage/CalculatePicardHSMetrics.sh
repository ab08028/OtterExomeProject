###### this should be part of paleomix pipeline

PICARD=/u/home/a/ab08028/klohmueldata/annabel_data/bin/Picard_2.8.1/picard.jar

# capture regions in sea otter genome:
# use targets and regions as the same (don't care about efficiency)
regions=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_CaptureScripts/captureRegions/byCategory_exon_neutral_promoter/dup99removed/finalFiles-0based-Correction_20170906_USEFROMNOWON/allCaptureRegions_IDCategoriesLengths_ForMyUse.DupScaff99Removed.0-based.20170906.PositionsOnlyForUseInBedtools.sorted.merged.noHeader.bed
targets=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_CaptureScripts/captureRegions/byCategory_exon_neutral_promoter/dup99removed/finalFiles-0based-Correction_20170906_USEFROMNOWON/allCaptureRegions_IDCategoriesLengths_ForMyUse.DupScaff99Removed.0-based.20170906.PositionsOnlyForUseInBedtools.sorted.merged.noHeader.bed
# convert
refDir=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/
refPrefix=sea_otter_23May2016_bS9RH.deduped.99
refDict=$refDir/${refPrefix}.dict 
ref=$refDir/${refPrefix}.fasta
java -jar $PICARD BedToIntervalList \
       I=$regions \
       O=${regions%.bed}.interval_list \
       SD=$refDict


header=A19_Elut_CA_SM_30_SN2_screen
plmxDir=/u/flashscratch/a/ab08028/captures/paleomix/testMapping/
plmxBam=$plmxDir/$header/${header}.${refPrefix}.bam
bedtools coverage -a $capRegions -b $plmxBam > $header.$ref.$bedtoolsGenomecov.bed
# bg: report in bedgraph format

# then want to intersect with capture regions? Or plot whole thing but mark capture regions in red
# to see if they stand out. Yes I like that plan. 
# I think Picard is better
java -jar $PICARD CollectHsMetrics \
      I=$plmxBam \
      O=${plmxBam%.bam}.hs_metrics.region2target.txt \
      R=$ref \
      BAIT_INTERVALS=${regions%.bed}.interval_list \
      TARGET_INTERVALS=${regions%.bed}.interval_list