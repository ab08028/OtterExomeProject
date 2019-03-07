#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=8G
#$ -N picardHSMetrics
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
# note you need java 1.8; if you used -V and have a different java installed you may have to unload it
###### this should be part of paleomix pipeline

#headers=/u/flashscratch/a/ab08028/captures/samples/ancientSamples.capture.screen.all.txt
headers=/u/flashscratch/a/ab08028/captures/samples/aDNA.Screens.2.txt # 20190219 new screens
PICARD=/u/home/a/ab08028/klohmueldata/annabel_data/bin/Picard_2.8.1/picard.jar

# capture regions in sea otter genome:
# use targets and regions as the same (don't care about efficiency)
regions=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_CaptureScripts/captureRegions/byCategory_exon_neutral_promoter/dup99removed/finalFiles-0based-Correction_20170906_USEFROMNOWON/allCaptureRegions_IDCategoriesLengths_ForMyUse.DupScaff99Removed.0-based.20170906.PositionsOnlyForUseInBedtools.sorted.merged.noHeader.bed
# convert
refDir=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/
refPrefix=sea_otter_23May2016_bS9RH.deduped.99
refDict=$refDir/${refPrefix}.dict 
ref=$refDir/${refPrefix}.fasta
# only do this once --> 
#java -jar $PICARD BedToIntervalList \
#       I=$regions \
#       O=${regions%.bed}.interval_list \
#       SD=$refDict
plmxDir=/u/flashscratch/a/ab08028/captures/paleomix/testMapping/
sumStatDir=/u/flashscratch/a/ab08028/captures/paleomix/summaryStats
todaysdate=`date +%Y%m%d`
# set up summary stats file; added "SAMPLE" column to the end (adding sample name in AWK), removed lib group field
# figured out that HS library size field is empty in picard output; so removed that from header and used tr -s to squeeze out the column in my awk command in the loop below
echo -e "BAIT_SET\tGENOME_SIZE\tBAIT_TERRITORY\tTARGET_TERRITORY\tBAIT_DESIGN_EFFICIENCY\tTOTAL_READS\tPF_READS\tPF_UNIQUE_READS\tPCT_PF_READS\tPCT_PF_UQ_READS\tPF_UQ_READS_ALIGNED\tPCT_PF_UQ_READS_ALIGNED\tPF_BASES_ALIGNED\tPF_UQ_BASES_ALIGNED\tON_BAIT_BASES\tNEAR_BAIT_BASES\tOFF_BAIT_BASES\tON_TARGET_BASES\tPCT_SELECTED_BASES\tPCT_OFF_BAIT\tON_BAIT_VS_SELECTED\tMEAN_BAIT_COVERAGE\tMEAN_TARGET_COVERAGE\tMEDIAN_TARGET_COVERAGE\tPCT_USABLE_BASES_ON_BAIT\tPCT_USABLE_BASES_ON_TARGET\tFOLD_ENRICHMENT\tZERO_CVG_TARGETS_PCT\tPCT_EXC_DUPE\tPCT_EXC_MAPQ\tPCT_EXC_BASEQ\tPCT_EXC_OVERLAP\tPCT_EXC_OFF_TARGET\tFOLD_80_BASE_PENALTY\tPCT_TARGET_BASES_1X\tPCT_TARGET_BASES_2X\tPCT_TARGET_BASES_10X\tPCT_TARGET_BASES_20X\tPCT_TARGET_BASES_30X\tPCT_TARGET_BASES_40X\tPCT_TARGET_BASES_50X\tPCT_TARGET_BASES_100X\tHS_PENALTY_10X\tHS_PENALTY_20X\tHS_PENALTY_30X\tHS_PENALTY_40X\tHS_PENALTY_50X\tHS_PENALTY_100X\tAT_DROPOUT\tGC_DROPOUT\tHET_SNP_SENSITIVITY\tHET_SNP_Q\tSAMPLE" > $sumStatDir/picard.HSMetrics.plmx.${todaysdate}.txt



cat $headers | while read header
do
echo $header
plmxBam=$plmxDir/$header/${header}.${refPrefix}.bam

java -jar $PICARD CollectHsMetrics \
      I=$plmxBam \
      O=${plmxBam%.bam}.hs_metrics.region2target.txt \
      R=$ref \
      BAIT_INTERVALS=${regions%.bed}.interval_list \
      TARGET_INTERVALS=${regions%.bed}.interval_list
done

# add to summary stats:
# use tr -s to squeeze out repeat tabs:
cat $headers | while read header
do
plmxBam=$plmxDir/$header/${header}.${refPrefix}.bam 
# using tr to squeeze out empty columns; fixed header above accordingly. this works now. Note that empty fields really messes stuff up.
grep -A1 "BAIT_SET" ${plmxBam%.bam}.hs_metrics.region2target.txt | grep -v "BAIT_SET" | awk -F '\t' -v header=$header '{OFS="\t";print $0, header}'  | tr -s '\t' '\t' >> $sumStatDir/picard.HSMetrics.plmx.${todaysdate}.txt
done

