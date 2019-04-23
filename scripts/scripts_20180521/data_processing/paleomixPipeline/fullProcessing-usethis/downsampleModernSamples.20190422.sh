# Downsample modern bam files to match mean unique reads mapped to sea otter and ferret genomes in the best 3 ancient samples
# Mean reads mapped to sea otter:  967603
# Mean reads mapped to ferret:  505603
wd=/u/flashscratch/a/ab08028/captures/paleomix/fullProcessing/
downsampledir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/downsampledBams/downsample_allEven/

echo "downsample to equal mean ancient"

samtools view -s 1.0394 -b $wd/141_Elut_CA_419/141_Elut_CA_419.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.bam > $downsampledir/141_Elut_CA_419.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.1.0394.bam
 samtools view -s 1.0738 -b $wd/141_Elut_CA_419/141_Elut_CA_419.sea_otter_23May2016_bS9RH.deduped.99.bam > $downsampledir/141_Elut_CA_419.sea_otter_23May2016_bS9RH.deduped.99.downsamp.1.0738.bam
 samtools view -s 1.029 -b $wd/116_Elut_CA_307/116_Elut_CA_307.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.bam > $downsampledir/116_Elut_CA_307.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.1.029.bam
 samtools view -s 1.0547 -b $wd/116_Elut_CA_307/116_Elut_CA_307.sea_otter_23May2016_bS9RH.deduped.99.bam > $downsampledir/116_Elut_CA_307.sea_otter_23May2016_bS9RH.deduped.99.downsamp.1.0547.bam
 samtools view -s 1.0438 -b $wd/126_Elut_AK_AF3394/126_Elut_AK_AF3394.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.bam > $downsampledir/126_Elut_AK_AF3394.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.1.0438.bam
 samtools view -s 1.0845 -b $wd/126_Elut_AK_AF3394/126_Elut_AK_AF3394.sea_otter_23May2016_bS9RH.deduped.99.bam > $downsampledir/126_Elut_AK_AF3394.sea_otter_23May2016_bS9RH.deduped.99.downsamp.1.0845.bam
 samtools view -s 1.0304 -b $wd/55_Elut_AK_AF3736/55_Elut_AK_AF3736.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.bam > $downsampledir/55_Elut_AK_AF3736.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.1.0304.bam
 samtools view -s 1.0584 -b $wd/55_Elut_AK_AF3736/55_Elut_AK_AF3736.sea_otter_23May2016_bS9RH.deduped.99.bam > $downsampledir/55_Elut_AK_AF3736.sea_otter_23May2016_bS9RH.deduped.99.downsamp.1.0584.bam
 samtools view -s 1.0276 -b $wd/129_Elut_AK_AL4660/129_Elut_AK_AL4660.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.bam > $downsampledir/129_Elut_AK_AL4660.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.1.0276.bam
 samtools view -s 1.0536 -b $wd/129_Elut_AK_AL4660/129_Elut_AK_AL4660.sea_otter_23May2016_bS9RH.deduped.99.bam > $downsampledir/129_Elut_AK_AL4660.sea_otter_23May2016_bS9RH.deduped.99.downsamp.1.0536.bam

echo "count downsampled reads"

echo downsampledReadCounts > $downsampledir/downsampledReadCounts.txt
samtools flagstat $downsampledir/141_Elut_CA_419.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.1.0394.bam | head -n1 | awk '{cat $1}' >> $downsampledir/downsampledReadCounts.txt
 samtools flagstat $downsampledir/141_Elut_CA_419.sea_otter_23May2016_bS9RH.deduped.99.downsamp.1.0738.bam | head -n1 | awk '{cat $1}' >> $downsampledir/downsampledReadCounts.txt
 samtools flagstat $downsampledir/116_Elut_CA_307.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.1.029.bam | head -n1 | awk '{cat $1}' >> $downsampledir/downsampledReadCounts.txt
 samtools flagstat $downsampledir/116_Elut_CA_307.sea_otter_23May2016_bS9RH.deduped.99.downsamp.1.0547.bam | head -n1 | awk '{cat $1}' >> $downsampledir/downsampledReadCounts.txt
 samtools flagstat $downsampledir/126_Elut_AK_AF3394.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.1.0438.bam | head -n1 | awk '{cat $1}' >> $downsampledir/downsampledReadCounts.txt
 samtools flagstat $downsampledir/126_Elut_AK_AF3394.sea_otter_23May2016_bS9RH.deduped.99.downsamp.1.0845.bam | head -n1 | awk '{cat $1}' >> $downsampledir/downsampledReadCounts.txt
 samtools flagstat $downsampledir/55_Elut_AK_AF3736.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.1.0304.bam | head -n1 | awk '{cat $1}' >> $downsampledir/downsampledReadCounts.txt
 samtools flagstat $downsampledir/55_Elut_AK_AF3736.sea_otter_23May2016_bS9RH.deduped.99.downsamp.1.0584.bam | head -n1 | awk '{cat $1}' >> $downsampledir/downsampledReadCounts.txt
 samtools flagstat $downsampledir/129_Elut_AK_AL4660.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.downsamp.1.0276.bam | head -n1 | awk '{cat $1}' >> $downsampledir/downsampledReadCounts.txt
 samtools flagstat $downsampledir/129_Elut_AK_AL4660.sea_otter_23May2016_bS9RH.deduped.99.downsamp.1.0536.bam | head -n1 | awk '{cat $1}' >> $downsampledir/downsampledReadCounts.txt
