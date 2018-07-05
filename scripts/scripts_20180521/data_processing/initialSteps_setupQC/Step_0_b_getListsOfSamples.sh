#### Make list of files to process:

######## this script assumes fastq name is [SAMPLE ID]_SXX_R[12]_001.fastq.gz 
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
fastqs=$wd/fastqs
outdir=$wd/samples
mkdir -p $outdir

########## THIS DOESN"T DEAL WITH SN properly
# pull out sample IDs:

# ancient:
ls $fastqs | grep -E ^A[0-9]+.*gz | sed -e 's/_S[0-9]*_R.*_.*fastq.gz//g' | sort | uniq > $outdir/ancientSamples.txt

# modern:
ls $fastqs | grep -E ^[0-9]+_Elut_.*gz | sed -e 's/_S.*_R.*_.*fastq.gz//g' | sort | uniq > $outdir/modernSamples.txt

# dog
ls $fastqs | grep Cfam | sed -e 's/_S.*_R.*_.*fastq.gz//g' | sort | uniq > $outdir/dogSamples.txt

# all sea otter, plus blank, but no dog:
ls $fastqs | grep fastq.gz | grep -v Cfam | sed -e 's/_S.*_R.*_.*fastq.gz//g' | sort | uniq > $outdir/allElutSamples.txt

# then when submitting jobs you can do
# cat ancientSamples.txt | while read sample
# do
# ls $fastqs/$sample*
# done
# or whatever else you want to do 