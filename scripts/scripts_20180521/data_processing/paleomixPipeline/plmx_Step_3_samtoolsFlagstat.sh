module load samtools
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/paleomix/ancient_FirstRun_noIndelR_20180712
makefileDir=$scriptDir/makefiles/ancientMakefiles
headers=$SCRATCH/captures/samples/ancientSamples.txt
REFPREFIX=Mustela_putorius_furo.MusPutFur1.0.dna.toplevel

cat $headers | while read header
do 
samtools flagstat ${wd}/${header}/${header}.${REFPREFIX}.bam > ${wd}/${header}/${header}.${REFPREFIX}.flagstat
done
