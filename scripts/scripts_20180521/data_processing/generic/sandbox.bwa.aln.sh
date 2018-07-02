## sandbox of bwa aln stuff
source /u/local/Modules/default/init/modules.sh
module load samtools
module load bwa

numThreads=10
REFERENCE=$1
fileR1=$2
fileR2=$3
header=
# make sai files first: 
bwa aln -t $numThreads -n 0.01 -o 2 -l 100000 $REFERENCE $fileR1 > ${header}.R1.aln.sai
bwa aln -t $numThreads -n 0.01 -o 2 -l 100000 $REFERENCE $fileR2 > ${header}.R2.aln.sai
bwa sampe $REFERENCE ${header}.R1.aln.sai ${header}.R2.aln.sai $fileR1 $fileR2 | samtools view -h -bS > ${header}_aligned.bam
samtools flagstat ${header}_aligned.bam > ${header}_aligned.flagstat
              
# do I also need to do something to single-end reads that don't have a pair (kicked out of adapter removal?)
# what                   
## this lets you figure out %% mapping (endogenous)
# Then markduplicates step (will be same as before just change what bam file goes in?)
#  	-o Maximum number of gap opens [1] 
# 	-n The fraction of missing alignments given 2% uniform base error rate if FLOAT. The maximum edit distance is automatically chosen for different read lengths. [0.04] 
# 	-l Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. For long reads, this option is typically ranged from 25 to 35 for Ô-k 2Õ. [inf]
# make l very big to disable seed