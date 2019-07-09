#! /bin/bash
#$ -cwd
#$ -l h_rt=40:00:00,h_data=2G,highp
#$ -m abe
#$ -pe shared 16
#$ -M ab08028
#$ -N neutralCDSBedtools
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd


source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env
module load vcftools
module load bedtools


########## script to get continuity input format ###########
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/

wd=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/

mkdir -p $wd

## modern GATK genotype file:
vcfDate=20181119
vcfDir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${vcfDate}_filtered
vcf=snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz
outdir=$wd/continuity
mac=1
########### California: get frequency information ###################
pop=CA
keep=/u/flashscratch/a/ab08028/captures/samples/keep.${pop}.${vcfDate}.txt # list of 7 individuals from /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/easySFSPopMapFiles/samplesPop.Headers.forEasySFS.3.20181119.txt
output=${pop}.freqs.fromModernData.${vcfDate}.mac.${mac}
# sets that there has to be at least one copy of a minor allele (excludes freq 0 and 1)
vcftools --gzvcf $vcfDir/$vcf --keep $keep --freq2 --out $outdir/modernDataGATK_Freqs/$output --mac ${mac}  # 
bedhead="#chrom\tstart0based\tend\tmarkerID\tempty5\tempty6\tempty7\tempty8\tempty9\tempty10\tempty11\tempty12"
frqhead="CHROM\tPOS\tN_ALLELES\tN_CHR\tREF_FREQ\tALT_FREQ\n"
comboheader1=`echo -e "$bedhead\t$frqhead"`
echo -e $comboheader1 > $outdir/${output}.0based.bed
grep -v "CHROM" $outdir/modernDataGATK_Freqs/${output}.frq | awk '{OFS="\t";print $1,$2-1,$2,$1"_"$2,".",".",".",".",".",".",".",".",$0}' >> $outdir/modernDataGATK_Freqs/${output}.0based.bed

# this gives output as :
#CHROM   POS     N_ALLELES       N_CHR   {ALLELE:FREQ}
#										# this is ref	# this is alt # I checked it against vcf file
# GL896898.1      23947   2       14      A:0.857143      G:0.142857
# GL896898.1      38867   2       14      T:0.857143      C:0.142857
# GL896898.1      7745485 2       14      T:0.0714286     C:0.928571 # here C is greater than T
# so that makes things easier! 

# you can tell it's ref/alt not major/minor
# because sometimes the second column is greater than the first 
# okay so that works
# then convert this to bed format (0 based remember) 
# and merge with dumpCounts where it has A T C G (see angsd doc)

########### Alaska: get frequency information ###################
pop=AK
keep=/u/flashscratch/a/ab08028/captures/samples/keep.${pop}.${vcfDate}.txt # list of 7 individuals from /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/easySFSPopMapFiles/samplesPop.Headers.forEasySFS.3.20181119.txt
output=${pop}.freqs.fromModernData.${vcfDate}.mac.${mac}
vcftools --gzvcf $vcfDir/$vcf --keep $keep --freq --out $outdir/modernDataGATK_Freqs/$output --mac ${mac} # 

# convert to bed format
# header:
bedhead="#chrom\tstart0based\tend\tmarkerID\tempty5\tempty6\tempty7\tempty8\tempty9\tempty10\tempty11\tempty12"
frqhead="CHROM\tPOS\tN_ALLELES\tN_CHR\tREF_FREQ\tALT_FREQ\n"
comboheader1=`echo -e "$bedhead\t$frqhead"`
echo -e $comboheader1 > $outdir/${output}.0based.bed
grep -v "CHROM" $outdir/modernDataGATK_Freqs/${output}.frq | awk '{OFS="\t";print $1,$2-1,$2,$1"_"$2,".",".",".",".",".",".",".",".",$0}' >> $outdir/modernDataGATK_Freqs/${output}.0based.bed

##################################################
############## angsd dump counts 4 ###############
##################################################

# want to ONLY use aDNA samples (don't care about high/low coverage)
# want to do elut and mfur though
# want to not do GLs or anything, just do counts
# so need bamList of just ancient 

### list of bam files to include: aDNA **only** ### 
elutBamList=$scriptDir/data_processing/variant_calling_aDNA/bamLists/angsd.ancient.bamList.mappedtoElutfullpaths.txt # ancient only 
mfurBamList=$scriptDir/data_processing/variant_calling_aDNA/bamLists/angsd.ancient.bamList.mappedtoMfurfullpaths.txt # ancient only 

# reference genomes:
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

snpCutoff=1E-06
trimValue=7
countsdir='ancientOnly-counts'
mkdir -p $outdir/$countsdir
# trying output in beagle format  doGlf 2
####### Mfur mapped bams ############
spp="mfur"
ref=$mfurRef
bamList=$mfurBamList
basename=angsdOut.mappedTo${spp}
# can't use rmTrans or snp value cutoff when just dumping counts. it will just dump counts for all sites
# 
angsd -nThreads 16 \
-ref $ref \
-bam $bamList \
-remove_bads 1 -uniqueOnly 1 \
-C 50 -baq 1 -trim $trimValue -minQ 20 -minMapQ 25 \
-out $outdir/$countsdir/$basename \
-doCounts 1 -dumpCounts 4

# need to combine .pos and .count 

angsdheaders=`paste <(zcat $outdir/$countsdir/${basename}.pos.gz | head -n1) <(zcat $outdir/$countsdir/${basename}.counts.gz | head -n1)`
comboheader2=`echo -e "$bedhead\t$angsdheaders"`
echo -e $comboheader2 > $outdir/$countsdir/${basename}.counts.0based.bed
paste <(zcat $outdir/$countsdir/${basename}.pos.gz) <(zcat $outdir/$countsdir/${basename}.counts.gz) | grep -v "totDepth" | awk '{OFS="\t";print $1,$2-1,$2,$1"_"$2,".",".",".",".",".",".",".",".",$0}' | sed 's/\t\t/\t/g' | sed 's/\t$//g' >> $outdir/$countsdir/${basename}.counts.0based.bed # go into awk and rearrange to make it bed format with extra columns 
################################## intersect #####################################################
############# then using bedtools intersect, want to combine the counts and frequencies ##########
# get a header:
comboheader3=`echo -e "$comboheader1\t$comboheader2"`
combodir=$wd/continuity/combinedCounts-Freqs
mkdir -p $combodir
echo -e $comboheader3 > $combodir/${basename}.ancient.counts.freqsFromModernGATK.superfile.0based.bed
bedtools intersect -a $outdir/modernDataGATK_Freqs/${output}.0based.bed -b $outdir/$countsdir/${basename}.counts.0based.bed -wa -wb >> $combodir/${basename}.ancient.counts.freqsFromModernGATK.superfile.0based.bed

####### Elut mapped bams ############
spp="elut"
ref=$elutRef
bamList=$elutBamList
basename=angsdOut.mappedTo${spp}
# can't use rmTrans or snp value cutoff when just dumping counts. it will just dump counts for all sites
# 
angsd -nThreads 16 \
-ref $ref \
-bam $bamList \
-remove_bads 1 -uniqueOnly 1 \
-C 50 -baq 1 -trim $trimValue -minQ 20 -minMapQ 25 \
-out $outdir/$countsdir/$basename \
-doCounts 1 -dumpCounts 4

# need to combine .pos and .count 

angsdheaders=`paste <(zcat $outdir/$countsdir/${basename}.pos.gz | head -n1) <(zcat $outdir/$countsdir/${basename}.counts.gz | head -n1)`
comboheader2=`echo -e "$bedhead\t$angsdheaders"`
echo -e $comboheader2 > $outdir/$countsdir/${basename}.counts.0based.bed
paste <(zcat $outdir/$countsdir/${basename}.pos.gz) <(zcat $outdir/$countsdir/${basename}.counts.gz) | grep -v "totDepth" | awk '{OFS="\t";print $1,$2-1,$2,$1"_"$2,".",".",".",".",".",".",".",".",$0}' | sed 's/\t\t/\t/g' | sed 's/\t$//g' >> $outdir/$countsdir/${basename}.counts.0based.bed # go into awk and rearrange to make it bed format with extra columns 
################################## intersect #####################################################
############# then using bedtools intersect, want to combine the counts and frequencies ##########
# get a header:
comboheader3=`echo -e "$comboheader1\t$comboheader2"`
combodir=$wd/continuity/combinedCounts-Freqs
mkdir -p $combodir
echo -e $comboheader3 > $combodir/${basename}.ancient.counts.freqsFromModernGATK.superfile.0based.bed
bedtools intersect -a $outdir/modernDataGATK_Freqs/${output}.0based.bed -b $outdir/$countsdir/${basename}.counts.0based.bed -wa -wb >> $combodir/${basename}.ancient.counts.freqsFromModernGATK.superfile.0based.bed

############### convert to bed format 
source deactivate

