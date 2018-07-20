module load samtools
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/paleomix/ancient_FirstRun_noIndelR_20180712
headers=$SCRATCH/captures/samples/ancientSamples.txt
REFPREFIX=Mustela_putorius_furo.MusPutFur1.0.dna.toplevel

mkdir -p $SCRATCH/captures/flagstat
> $SCRATCH/captures/flagstat/percentMapped.all.txt

cat $headers | while read header
do 
samtools flagstat ${wd}/${header}/${header}.${REFPREFIX}.bam > $SCRATCH/captures/flagstat/${header}.${REFPREFIX}.flagstat
# get mapped stats
stat=`grep "mapped (" $SCRATCH/captures/flagstat/${header}.${REFPREFIX}.flagstat`
echo $header $stat >> $SCRATCH/captures/flagstat/percentMapped.all.txt

done


