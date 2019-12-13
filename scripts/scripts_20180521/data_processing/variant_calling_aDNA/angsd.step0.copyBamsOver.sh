###### copy bams: formally: 
wd=$SCRATCH/captures
headers=$wd/samples/modern.SamplesToCompareToADNA.COM.AL.KUR.txt # adding 3 COM, KUR, AL samples
indir=$wd/paleomix/fullProcessing/
# 
# if you want to copy all the samples that you're going to use:
cat $headers | while read header
do
copy bams and bai files
cp $indir/$header/*realigned* /u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/bams
done

# just the new samples on 20191212:
# newSampleHeaders="100_Elut_KUR_24 101_Elut_KUR_3 102_Elut_KUR_4 111_Elut_AL_AT_GE91133 121_Elut_AL_AT_GE91135 124_Elut_AL_AT_GE91143 136_Elut_BER_46 137_Elut_BER_88 50_Elut_BER_100"	
# for header in $newSampleHeaders
# do
# echo $header
# ls -alh $indir/$header/*realigned*
# cp -i $indir/$header/*realigned* /u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/bams
# done

# just waiting on 50 and 137 then can get started with angsd
