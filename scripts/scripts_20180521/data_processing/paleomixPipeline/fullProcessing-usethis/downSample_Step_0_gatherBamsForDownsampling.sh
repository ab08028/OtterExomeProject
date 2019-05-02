###### Gather the relevant bam files #######
#### Make sure to use REALIGNED bams ####
SCRATCH=/u/flashscratch/a/ab08028/
captures=$SCRATCH/captures/
wd=$captures/aDNA-ModernComparison
plmxDir=$captures/paleomix/fullProcessing
headers=/u/flashscratch/a/ab08028/captures/samples/ancient.modern.comparison.txt



refs="sea_otter_23May2016_bS9RH.deduped.99 Mustela_putorius_furo.MusPutFur1.0.dna.toplevel"


cat $headers | while read header
do
for ref in $refs
do
####### very important to use REALIGNED bam file!!!! otherwise it hasn't gone through indel realignment #####
cp $plmxDir/$header/$header.$ref.realigned.ba* $wd/bams
# using *ba* to get the bam indices
done
done

# want to copy these into the correct dir