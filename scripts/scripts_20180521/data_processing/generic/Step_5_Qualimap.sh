#! /bin/bash
#$ -cwd
#$ -l h_rt=20:00:00,h_data=2G,highp
#$ -pe shared 15
#$ -N qualimap
#$ -m abe

# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
fastqs=$wd/fastqs
bams=$wd/bams
outdir=$wd/qualimap
mkdir -p $outdir
header=$1

qmap=/u/home/a/ab08028/klohmueldata/annabel_data/bin/qualimap_v2.2.1/qualimap
source /u/local/Modules/default/init/modules.sh
module load R
module load java

$qmap bamqc --java-mem-size=24G -nt 15 -bam $bams/$header.bam -outdir $outdir -outfile ${header}.qualimap.out.pdf
# take out -c 

