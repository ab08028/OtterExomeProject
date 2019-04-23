######## Samtools downsampling #######

# Each modern sample is going to get a specific fraction to downsample by
# a) to match the mean(?) number of uniquely mapping reads to mfur or elut across the aDNA samples
# or
# b) to match a specific ancient sample (chosen at random of the three)

####### have to do for each sample separately (has its own fraction -- though maybe I could pull them ) ##########
wd=/u/flashscratch/a/ab08028/captures/paleomix/fullProcessing
downsampledir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/downsampledBams/sandbox

header=55_Elut_AK_AF3736
ref=sea_otter_23May2016_bS9RH.deduped.99
# example
input=$header.${ref}.bam 
frac=.058 # make sure you don't lead with a zero; this differs between refs too 
# count original reads:
samtools flagstat $wd/$header/$input | head -n1 | awk '{print $1}' > $downsampledir/$header.original.ReadCount.txt
for rep in {1..1}
do
# note: using rep as the seed 
# # kind of odd behavior -- when you input -s, the integer part is the seed (could be any integer)
# # and the decimal part is the fraction you want (what a bizarre approach)
# # so sampling -s 3.076 would set seed as 3 (so you can replicate the sampling)
# # and then samples 0.076 of the reads. That is super weird

# using rep as the seed so then you can replicate it 
echo "rep: $rep ; seed: $seed ; frac = $frac ; -s : $seed.$frac \n where integer part is the seed and after decimal is desired frac" > $downsampledir/${header}.downsample.parameters.txt
samtools view -s $rep$frac -b $wd/$header/$input > $downsampledir/${header}.downsample.$rep.0$frac.bam
# count number of before and after reads and write out: 
samtools flagstat $downsampledir/${header}.downsample.$rep.0$frac.bam | head -n1 | awk '{print $1}' > $downsampledir/$header.$rep.0$frac.ReadCount.txt
done
