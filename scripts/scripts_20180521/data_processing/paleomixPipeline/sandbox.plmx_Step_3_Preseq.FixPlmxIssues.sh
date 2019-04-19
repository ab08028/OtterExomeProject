### Want to look at library complexity of ANCIENT only 
## Issue: paleomix doesn't output the correct histogram if you have PCR duplicates set to "filter"
# doesn't correctly assess curve. Need to try with a different bam
# this runs fast; can be run in the shell with qrsh

module load samtools
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/paleomix/testMapping
#headers=$SCRATCH/captures/samples/ancientSamples.txt
headers=A13_Elut_CA_AN_388_SN1_2CAP_screen
REFPREFIX="Mustela_putorius_furo.MusPutFur1.0.dna.toplevel sea_otter_23May2016_bS9RH.deduped.99"

preseq=/u/home/a/ab08028/klohmueldata/annabel_data/bin/preseq_v2.0/preseq

mkdir -p $SCRATCH/captures/preseq
outdir=$SCRATCH/captures/preseq


### not sure why yet but preseq not working with dup histograms from paleomix
# figured out why: if < 1M reads c_curve will not output anything; have to change the step size to 10,000 or 100,000 reads

cat $headers | while read header
do 
for REF in $REFPREFIX
do
# get the actual curve of the library from subsampling:
$preseq c_curve -s 100000 -H ${wd}/${header}/${header}.${REF}.duphist/*txt -o $outdir/${header}.${REF}.preseq.c_curve.txt
# get extrapolated curve if you sequenced the library more
$preseq lc_extrap -H ${wd}/${header}/${header}.${REF}.duphist/*txt -o $outdir/${header}.${REF}.preseq.lc_extrap.fromDupHist.txt
# get extrapolated curve from bam file before rmdup:
## $preseq lc_extrap -B $wd/$header/$header/$REF/$header/${header}_1a/Lane_1/collapsed.bam -o $outdir/${header}.${REF}.preseq.lc_extrap.fromCollapsedBam.txt
# use main bam file. doesn't matter that pcr dups are removed, because this is not about pcr dups it's about
# oversequencing reads 
$preseq lc_extrap -B $wd/$header/${header}.${REF}.bam -o $outdir/$header.$REF.lc_extrap.fromBam.txt
done
done
# then am going to make plots in R (plot # of unique reads vs. # of reads)
# may need to adjust step size 
###########
# sandbox #
############

wd=/u/flashscratch/a/ab08028/captures/paleomix/testMapping/original.filter.A30_Elut_CA_SM_35_SN1_CAP/A30_Elut_CA_SM_35_SN1_CAP/sea_otter_23May2016_bS9RH.deduped.99/A30_Elut_CA_SM_35_SN1_CAP/A30_Elut_CA_SM_35_SN1_CAP_1a/Lane_1

# testing it with non-mark dup bam (have dups been marked???)
$preseq c_curve -s 100000 -B collapsed.bam -o testPreseq.1.collapsed.c_curve.txt

$preseq lc_extrap -B collapsed.bam -o testPreseq.1.collapsed.lc_extrap.txt

# need to intersect it with bed file
regions=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_CaptureScripts/captureRegions/byCategory_exon_neutral_promoter/dup99removed/finalFiles-0based-Correction_20170906_USEFROMNOWON/allCaptureRegions_IDCategoriesLengths_ForMyUse.DupScaff99Removed.0-based.20170906.PositionsOnlyForUseInBedtools.sorted.merged.noHeader.bed
bedtools intersect -a $wd/collapsed.bam -b $regions > test.collapsed.regionsOnly.bam

# now run preseq on this bam:
$preseq c_curve -s 100000 -B test.collapsed.regionsOnly.bam -o testPreseq.1.collapsed.c_curve.regionsOnly.txt

$preseq lc_extrap -B test.collapsed.regionsOnly.bam -o testPreseq.1.collapsed.lc_extrap.regionsOnly.txt
