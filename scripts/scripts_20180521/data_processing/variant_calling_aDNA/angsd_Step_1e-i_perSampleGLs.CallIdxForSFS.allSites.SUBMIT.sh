gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA/
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison

script=angsd_Step_1e-i_perSampleGLs.CallIdxForSFS.allSites.sh
# step 1e-i submitter
elutBamList=$scriptDir/bamLists/angsd.bamList.mappedtoElutfullpaths.txt # list of bam files mapped to sea otter, including downsampled AND non-downsampled
mfurBamList=$scriptDir/bamLists/angsd.bamList.mappedtoMfurfullpaths.txt  # list of bam files mapped to ferret, including downsampled AND non-downsampled

# references:
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

# script usage:
# qsub script [ path to bam List for single individual ] [sampleID] [ref prefix (Mfur or Elut)] [ path to reference genome ]
# need to make a file for each bam
refPrefix=Mfur
cat $mfurBamList | while read bam
do
# get sample name from bam using basename:
filename=`basename $bam`
#sampleID=${filename%.Mustela*}
label=${filename%.newSampName.*} # gets info from filename 
### need to get downsampled label
echo $bam > $scriptDir/bamLists/$sampleID.$refPrefix.bamList.txt
qsub $scriptDir/${script} $scriptDir/bamLists/$sampleID.$refPrefix.bamList.txt $label $refPrefix $mfurRef
done

# elut-mapped:
refPrefix=Elut
cat $elutBamList | while read bam
do
filename=`basename $bam`
sampleID=${filename%.sea_otter*}
label=${filename%.newSampName.*} # gets info from filename 
echo $bam > $scriptDir/bamLists/$sampleID.$refPrefix.bamList.txt
qsub $scriptDir/${script} $scriptDir/bamLists/$sampleID.$refPrefix.bamList.txt $label $refPrefix $elutRef
done
