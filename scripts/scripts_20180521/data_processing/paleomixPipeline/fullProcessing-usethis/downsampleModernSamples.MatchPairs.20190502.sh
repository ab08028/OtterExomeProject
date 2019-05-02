#! /bin/bash
#$ -l h_rt=5:00:00,h_data=6G
#$ -m abe
#$ -M ab08028
#$ -N downsamplesource /u/local/Modules/default/init/modules.sh
module load samtools
# Downsample modern bam files to match the number unique reads mapped to sea otter and ferret genomes in the best 3 ancient samples (in pairs -- each ancient sample matches 1 CA and 1 AK sample)
# note: the -s value has the replicate number as the integer which is the seed and the decimal part is the fraction. So 1.006 is replicate 1, fraction 0.006
wd=/u/flashscratch/a/ab08028/captures/paleomix/fullProcessing/
downsampledir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/downsampledBams/downsample_Pairs/

echo downsampledReadCounts > $downsampledir/downsampledReadCounts.txt
# downsample 140_Elut_CA_403 to equal ancient sample A30_Elut_CA_SM_35_SN1_CAP
# 140_Elut_CA_403 (modern) starting reads: 18089541
# A30_Elut_CA_SM_35_SN1_CAP (ancient) starting reads: 1195428"

samtools view -s 1.0661 -b $wd/140_Elut_CA_403/140_Elut_CA_403.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.bam > $downsampledir/140_Elut_CA_403.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam
# Count the resulting reads to make sure it downsampled properly
echo "140_Elut_CA_403" >> $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/140_Elut_CA_403.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt
# Rename sample

samtools view -H $downsampledir/140_Elut_CA_403.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam  | sed "s/SM:[^	]*/SM:140_Elut_CA_403_downsamp/g" | samtools reheader - $downsampledir/140_Elut_CA_403.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam > $downsampledir/140_Elut_CA_403.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.newSampName.bam


# downsample 141_Elut_CA_419 to equal ancient sample A29_Elut_CA_SM_30_SN2_CAP
# 141_Elut_CA_419 (modern) starting reads: 12818602
# A29_Elut_CA_SM_30_SN2_CAP (ancient) starting reads: 898036"

samtools view -s 1.0701 -b $wd/141_Elut_CA_419/141_Elut_CA_419.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.bam > $downsampledir/141_Elut_CA_419.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam
# Count the resulting reads to make sure it downsampled properly
echo "141_Elut_CA_419" >> $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/141_Elut_CA_419.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt
# Rename sample

samtools view -H $downsampledir/141_Elut_CA_419.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam  | sed "s/SM:[^	]*/SM:141_Elut_CA_419_downsamp/g" | samtools reheader - $downsampledir/141_Elut_CA_419.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam > $downsampledir/141_Elut_CA_419.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.newSampName.bam


# downsample 116_Elut_CA_307 to equal ancient sample A13_Elut_CA_AN_388_SN1_2CAP_screen
# 116_Elut_CA_307 (modern) starting reads: 17463663
# A13_Elut_CA_AN_388_SN1_2CAP_screen (ancient) starting reads: 505603"

samtools view -s 1.029 -b $wd/116_Elut_CA_307/116_Elut_CA_307.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.bam > $downsampledir/116_Elut_CA_307.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam
# Count the resulting reads to make sure it downsampled properly
echo "116_Elut_CA_307" >> $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/116_Elut_CA_307.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt
# Rename sample

samtools view -H $downsampledir/116_Elut_CA_307.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam  | sed "s/SM:[^	]*/SM:116_Elut_CA_307_downsamp/g" | samtools reheader - $downsampledir/116_Elut_CA_307.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam > $downsampledir/116_Elut_CA_307.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.newSampName.bam


# downsample 126_Elut_AK_AF3394 to equal ancient sample A30_Elut_CA_SM_35_SN1_CAP
# 126_Elut_AK_AF3394 (modern) starting reads: 11547167
# A30_Elut_CA_SM_35_SN1_CAP (ancient) starting reads: 1195428"

samtools view -s 1.1035 -b $wd/126_Elut_AK_AF3394/126_Elut_AK_AF3394.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.bam > $downsampledir/126_Elut_AK_AF3394.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam
# Count the resulting reads to make sure it downsampled properly
echo "126_Elut_AK_AF3394" >> $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/126_Elut_AK_AF3394.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt
# Rename sample

samtools view -H $downsampledir/126_Elut_AK_AF3394.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam  | sed "s/SM:[^	]*/SM:126_Elut_AK_AF3394_downsamp/g" | samtools reheader - $downsampledir/126_Elut_AK_AF3394.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam > $downsampledir/126_Elut_AK_AF3394.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.newSampName.bam


# downsample 55_Elut_AK_AF3736 to equal ancient sample A29_Elut_CA_SM_30_SN2_CAP
# 55_Elut_AK_AF3736 (modern) starting reads: 16637734
# A29_Elut_CA_SM_30_SN2_CAP (ancient) starting reads: 898036"

samtools view -s 1.054 -b $wd/55_Elut_AK_AF3736/55_Elut_AK_AF3736.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.bam > $downsampledir/55_Elut_AK_AF3736.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam
# Count the resulting reads to make sure it downsampled properly
echo "55_Elut_AK_AF3736" >> $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/55_Elut_AK_AF3736.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt
# Rename sample

samtools view -H $downsampledir/55_Elut_AK_AF3736.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam  | sed "s/SM:[^	]*/SM:55_Elut_AK_AF3736_downsamp/g" | samtools reheader - $downsampledir/55_Elut_AK_AF3736.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam > $downsampledir/55_Elut_AK_AF3736.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.newSampName.bam


# downsample 129_Elut_AK_AL4660 to equal ancient sample A13_Elut_CA_AN_388_SN1_2CAP_screen
# 129_Elut_AK_AL4660 (modern) starting reads: 18345387
# A13_Elut_CA_AN_388_SN1_2CAP_screen (ancient) starting reads: 505603"

samtools view -s 1.0276 -b $wd/129_Elut_AK_AL4660/129_Elut_AK_AL4660.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.bam > $downsampledir/129_Elut_AK_AL4660.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam
# Count the resulting reads to make sure it downsampled properly
echo "129_Elut_AK_AL4660" >> $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/129_Elut_AK_AL4660.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt
# Rename sample

samtools view -H $downsampledir/129_Elut_AK_AL4660.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam  | sed "s/SM:[^	]*/SM:129_Elut_AK_AL4660_downsamp/g" | samtools reheader - $downsampledir/129_Elut_AK_AL4660.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.bam > $downsampledir/129_Elut_AK_AL4660.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.rep.1.newSampName.bam


# downsample 140_Elut_CA_403 to equal ancient sample A30_Elut_CA_SM_35_SN1_CAP
# 140_Elut_CA_403 (modern) starting reads: 18221144
# A30_Elut_CA_SM_35_SN1_CAP (ancient) starting reads: 2141919"

samtools view -s 1.1176 -b $wd/140_Elut_CA_403/140_Elut_CA_403.sea_otter_23May2016_bS9RH.deduped.99.bam > $downsampledir/140_Elut_CA_403.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam
# Count the resulting reads to make sure it downsampled properly
echo "140_Elut_CA_403" >> $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/140_Elut_CA_403.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt
# Rename sample

samtools view -H $downsampledir/140_Elut_CA_403.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam  | sed "s/SM:[^	]*/SM:140_Elut_CA_403_downsamp/g" | samtools reheader - $downsampledir/140_Elut_CA_403.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam > $downsampledir/140_Elut_CA_403.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.newSampName.bam


# downsample 141_Elut_CA_419 to equal ancient sample A29_Elut_CA_SM_30_SN2_CAP
# 141_Elut_CA_419 (modern) starting reads: 13106833
# A29_Elut_CA_SM_30_SN2_CAP (ancient) starting reads: 1506414"

samtools view -s 1.1149 -b $wd/141_Elut_CA_419/141_Elut_CA_419.sea_otter_23May2016_bS9RH.deduped.99.bam > $downsampledir/141_Elut_CA_419.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam
# Count the resulting reads to make sure it downsampled properly
echo "141_Elut_CA_419" >> $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/141_Elut_CA_419.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt
# Rename sample

samtools view -H $downsampledir/141_Elut_CA_419.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam  | sed "s/SM:[^	]*/SM:141_Elut_CA_419_downsamp/g" | samtools reheader - $downsampledir/141_Elut_CA_419.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam > $downsampledir/141_Elut_CA_419.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.newSampName.bam


# downsample 116_Elut_CA_307 to equal ancient sample A13_Elut_CA_AN_388_SN1_2CAP_screen
# 116_Elut_CA_307 (modern) starting reads: 17696897
# A13_Elut_CA_AN_388_SN1_2CAP_screen (ancient) starting reads: 967603"

samtools view -s 1.0547 -b $wd/116_Elut_CA_307/116_Elut_CA_307.sea_otter_23May2016_bS9RH.deduped.99.bam > $downsampledir/116_Elut_CA_307.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam
# Count the resulting reads to make sure it downsampled properly
echo "116_Elut_CA_307" >> $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/116_Elut_CA_307.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt
# Rename sample

samtools view -H $downsampledir/116_Elut_CA_307.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam  | sed "s/SM:[^	]*/SM:116_Elut_CA_307_downsamp/g" | samtools reheader - $downsampledir/116_Elut_CA_307.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam > $downsampledir/116_Elut_CA_307.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.newSampName.bam


# downsample 126_Elut_AK_AF3394 to equal ancient sample A30_Elut_CA_SM_35_SN1_CAP
# 126_Elut_AK_AF3394 (modern) starting reads: 11454425
# A30_Elut_CA_SM_35_SN1_CAP (ancient) starting reads: 2141919"

samtools view -s 1.187 -b $wd/126_Elut_AK_AF3394/126_Elut_AK_AF3394.sea_otter_23May2016_bS9RH.deduped.99.bam > $downsampledir/126_Elut_AK_AF3394.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam
# Count the resulting reads to make sure it downsampled properly
echo "126_Elut_AK_AF3394" >> $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/126_Elut_AK_AF3394.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt
# Rename sample

samtools view -H $downsampledir/126_Elut_AK_AF3394.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam  | sed "s/SM:[^	]*/SM:126_Elut_AK_AF3394_downsamp/g" | samtools reheader - $downsampledir/126_Elut_AK_AF3394.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam > $downsampledir/126_Elut_AK_AF3394.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.newSampName.bam


# downsample 55_Elut_AK_AF3736 to equal ancient sample A29_Elut_CA_SM_30_SN2_CAP
# 55_Elut_AK_AF3736 (modern) starting reads: 16561142
# A29_Elut_CA_SM_30_SN2_CAP (ancient) starting reads: 1506414"

samtools view -s 1.091 -b $wd/55_Elut_AK_AF3736/55_Elut_AK_AF3736.sea_otter_23May2016_bS9RH.deduped.99.bam > $downsampledir/55_Elut_AK_AF3736.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam
# Count the resulting reads to make sure it downsampled properly
echo "55_Elut_AK_AF3736" >> $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/55_Elut_AK_AF3736.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt
# Rename sample

samtools view -H $downsampledir/55_Elut_AK_AF3736.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam  | sed "s/SM:[^	]*/SM:55_Elut_AK_AF3736_downsamp/g" | samtools reheader - $downsampledir/55_Elut_AK_AF3736.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam > $downsampledir/55_Elut_AK_AF3736.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.newSampName.bam


# downsample 129_Elut_AK_AL4660 to equal ancient sample A13_Elut_CA_AN_388_SN1_2CAP_screen
# 129_Elut_AK_AL4660 (modern) starting reads: 18052751
# A13_Elut_CA_AN_388_SN1_2CAP_screen (ancient) starting reads: 967603"

samtools view -s 1.0536 -b $wd/129_Elut_AK_AL4660/129_Elut_AK_AL4660.sea_otter_23May2016_bS9RH.deduped.99.bam > $downsampledir/129_Elut_AK_AL4660.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam
# Count the resulting reads to make sure it downsampled properly
echo "129_Elut_AK_AL4660" >> $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/129_Elut_AK_AL4660.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam | head -n1 | awk '{print $1}' >> $downsampledir/downsampledReadCounts.txt
# Rename sample

samtools view -H $downsampledir/129_Elut_AK_AL4660.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam  | sed "s/SM:[^	]*/SM:129_Elut_AK_AL4660_downsamp/g" | samtools reheader - $downsampledir/129_Elut_AK_AL4660.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.bam > $downsampledir/129_Elut_AK_AL4660.sea_otter_23May2016_bS9RH.deduped.99.downsamp.rep.1.newSampName.bam


