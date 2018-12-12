################ TREEMIX ############# 
# this should be run *after* you've converted your vcf to plink format for faststructure
# maybe run with all_8? for now running with all-7...
module load treemix
module load plink

# converter script:
# wget https://bitbucket.org/nygcresearch/treemix/downloads/plink2treemix.py
# from treemix 3/12/12:
# Added a small script to convert stratified allele frequencies output from plink into TreeMix format. 
# This will be incorporated into the next release, but for the moment must be downloaded 
# separately. To run this, let's say you have data in plink format (e.g., data.bed, data.bim, 
# data.fam) and a plink cluster file matching each individual to a population (data.clust).
# data was formatted for PLINK to run faststructure (see those scripts)
## change third column of .fam to be the pop identifier, as save as population.clusters
# only want first three columns
#awk '{OFS="\t"; print $1,$2,$3}' population.clusters.temp > population.clusters

gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=$gitdir/analyses/TREEMIX/
genotypeDate=20181119
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/
plinkFileDir=$vcfdir/plinkFormat/
treeFileDir=$vcfdir/treemixFormat/
mkdir -p $treeFileDir
clusters=$plinkFileDir/population.clusters # 3 columns: 0 sampleID popID
header=snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants

# must add a snp name to my .bim file:
# only once!!! mv snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.bim original.bim
awk '{$2 = $1":"$4; print}' $plinkFileDir/original.bim >  $plinkFileDir/snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.bim


# get minor allele frequencies stratified by pop from plink
# https://www.cog-genomics.org/plink/1.9/input#within
plink --bfile $plinkFileDir/$header \
--freq \
--missing \
--within $clusters \
--allow-extra-chr \
--out $plinkFileDir/$header \
--nonfounders \
--keep-allele-order 

gzip -f $plinkFileDir/${header}.frq.strat


# output
# CHR  (scaffold) SNP (.)    CLST  (population) A1  (allele 1) A2 (allele 2)     MAF (minor alelle freq)   MAC (minor allele count)  NCHROBS (non missing allele count)
# need allow extra chr to have non standard chr names
# need nonfounders to not just calculate allelle freqs based on 'founders' (there are no founders in my dataset)
# need keep allele order so it doesn't switch A1 and A2 ref/alt alleles

# bfile: plink.bed + plink.bim + plink.fam  references
# --freq writes a MAF report
# --missing makes missing data report 
# modify the third column of the .fam file :
# https://www.cog-genomics.org/plink/1.9/input#within
# note:
# If you have binary plink files then open your FAM file on some editor (Excel?). Keep 1st and 2nd column and add 3rd column as your clusters: A, B, C etc. and save it as mycluster.dat. Use this file with --within option.

######## NOTE: the snps in your .bim file MUST BE NAMED. otherwise treemix thinks there's only one snp and stops.
# You can fix this by https://www.biostars.org/p/114616/
# adding a name to each snp in your bim file that is scaff:locus 
# put in code to do this above 
python $scriptdir/plink2treemix.py $plinkFileDir/${header}.frq.strat.gz $treeFileDir/${header}.frq.strat.treemixFormat.gz
