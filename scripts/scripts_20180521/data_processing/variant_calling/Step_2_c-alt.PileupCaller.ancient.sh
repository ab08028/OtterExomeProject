######### samtools mpileup

source /u/local/Modules/default/init/modules.sh
module load samtools

REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
REFPREFIX=Mustela_putorius_furo.MusPutFur1.0.dna.toplevel
# header:
header=$1

# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
paleomixOutput=$wd/paleomix/$header
outdir=$wd/pileupCaller
mkdir -p $outdir

# these are covered intervals (result of previous step) -- min coverage of 1 read.
intervals=$wd/coveredIntervals
intervalFile=$intervals/${header}.coveredIntervals.list


samtools mpileup -R -B -q30 -Q30 -l $intervalFile \
	-f $REFERENCE \
    Sample1.bam Sample2.bam Sample3.bam > pileup.txt
    
# -R:
# -B: disables auto quality recalibration; creates ref bias in low coverage data
# -q30/Q30: mappping and base quality min. 30