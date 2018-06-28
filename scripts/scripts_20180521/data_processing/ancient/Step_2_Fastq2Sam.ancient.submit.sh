######### This script will submit a series of jobs that convert fastq to sam and adds readgroup info

# location of github:  which may be on remote server or local laptop
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/
QSUB=/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub
scriptname=Step_2_FastqToSam.sh
seqCenter=Medgenome
platform=illumina
errorLocation=/u/flashscratch/a/ab08028/captures/reports
user=ab08028 # where emails are sent
# all of these are Lib1 for now (could change later if I need to but shouldn't have to 
# because not trying to combine with my pilot data)
# starting prefix number
START=1
# ending prefix number 
END=2
# ancient DNA is treated separately because start with A1... # 
for (( c=$START; c<=$END; c++ ))
do
fileR1=`ls A${c}_Elut*R1*fastq.gz` # the R1 fastq file
fileR2=`ls A${c}_Elut*R2*fastq.gz` # the R2 fastq file
header=${fileR1%_S*_R*} # this is the header sample name
# example of fastq header: @A00354:22:HFCWVDMXX:1:1101:5358:1031 2:N:0:NCACAACA+NCACAACA
# want the HFCWVDMXX identifier (differs across fastqs; but is the same between R1 and R2)
UNIT=`zcat $fileR1 | head -n1 | awk -F ":" '{print $3}'` # pull out platform flowcell unit info
$QSUB -e $errorLocation -o $errorLocation -M $user -N fq2sam${c} $scriptDir/generic/$scriptname $fileR1 $fileR2 $header_FastqToSam.bam ${header}_1a ${header} Lib1 $UNIT $platform $seqCenter
# clear variables just in case:
fileR1=""
fileR2=""
header=""
UNIT=""
sleep 30m # sleep between job submissions
done
