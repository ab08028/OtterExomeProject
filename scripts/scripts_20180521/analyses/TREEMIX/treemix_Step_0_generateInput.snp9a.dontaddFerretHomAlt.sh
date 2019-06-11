################ TREEMIX ############# 
# this should be run *after* you've converted your vcf to plink format
###### Involves some manual steps!!!!!!! ############ 

## It will make 3 versions of treemix input: 
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
############## get plink format files from fastStructure_Step_1_vcf2plinkBed.20181119.sh in the FASTSTRUCTURE analysis section ##########

gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=$gitdir/analyses/TREEMIX/
genotypeDate=20181119
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/
plinkFileDir=$vcfdir/plinkFormat/ 
treeFileDir=$vcfdir/treemixFormat/
mkdir -p $treeFileDir

#### DO FOR BOTH HEADERS (separate scripts) ##########

#header=snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants
header=snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants


# must add a snp name to my .bim file:
# make a backup
/bin/cp $plinkFileDir/$header.bim $plinkFileDir/${header}.original.bim
awk '{$2 = $1":"$4; print}' $plinkFileDir/${header}.original.bim >  $plinkFileDir/$header.bim
awk '{print $1,$2}' $plinkFileDir/$header.fam > $plinkFileDir/$header.samples
/bin/cp $plinkFileDir/$header.samples $plinkFileDir/$header.population.clusters 
clusters=$plinkFileDir/$header.population.clusters # 3 columns: 0 sampleID popID
# ########## MANUALLY ADD POP ASSIGNMENTS HERE ################
/bin/cp $clusters $clusters.backup


##################################### With separate CA/BAJA  ##############################
############# don't need to remove relatives, they are already gone from snp9a file. ######################

marker="sepCA-BAJ"
plink --bfile $plinkFileDir/$header \
--freq \
--missing \
--within $clusters \
--allow-extra-chr \
--out $plinkFileDir/$header.${marker} \
--nonfounders \
--keep-allele-order

gzip -f $plinkFileDir/${header}.${marker}.frq.strat

python $scriptdir/plink2treemix.py $plinkFileDir/${header}.${marker}.frq.strat.gz $treeFileDir/${header}.${marker}.frq.strat.treemixFormat.gz

######################## EXCLUDE BAJA ########################################
############# can use plink --filter ######################

# make a filter file:
# plink --file data --filter myfile.raw 1 --freq

#implies a file myfile.raw exists which has a similar format to phenotype and cluster files: that is, the first two columns are family and individual IDs; the third column is expected to be a numeric value (although the file can have more than 3 columns), and only individuals who have a value of 1 for this would be included in any subsequent analysis or file generation procedure. e.g. if myfile.raw were 
# copy cols 1, 2 and then a "1"

awk '{print $1,$2,1}' $plinkFileDir/$header.fam > $plinkFileDir/$header.exclList.allIndsIncluded
# 168_Elut_BAJ_TS2
# 169_Elut_BAJ_R1
echo " YOU MUST MANUALLY EXCLUDE BAJA  HERE" #### 

###### MANUALLY edit and rename as : $plinkFileDir/$header.exclList.rmBAJA
# then manually edit the file to set individuals you want to exclude to "0"
marker=noBAJA
plink --bfile $plinkFileDir/$header \
--freq \
--missing \
--within $clusters \
--allow-extra-chr \
--out $plinkFileDir/$header.${marker} \
--nonfounders \
--keep-allele-order \
--filter $plinkFileDir/$header.exclList.rmBAJA 1

gzip -f $plinkFileDir/${header}.${marker}.frq.strat

python $scriptdir/plink2treemix.py $plinkFileDir/${header}.${marker}.frq.strat.gz $treeFileDir/${header}.${marker}.frq.strat.treemixFormat.gz

###### cool this works.






########################## old portions of script, including adding ferret (dont use) ###############
########################### dont want: with CA and BAJA combined  (see below for them separated )####################################
### assign Alaska and AK-ADMIX to different pops, assigned CA-BAJ to same pop
# need to assign populations MANUALLY and need to be formatted 0 ID POP: 
#awk '{print $1,$2}' $plinkFileDir/$header.fam > $plinkFileDir/$header.samples
#/bin/cp $plinkFileDir/$header.samples $plinkFileDir/$header.population.clusters 
#echo " AT THIS STAGE YOU MUST ASSIGN POPULATIONS MANUALLY TO POPULATION.CLUSTERS"
#clusters=$plinkFileDir/$header.population.clusters # 3 columns: 0 sampleID popID
# make a back up once you've added populations
#/bin/cp $clusters $clusters.backup
# get minor allele frequencies stratified by pop from plink
# https://www.cog-genomics.org/plink/1.9/input#within
# also need to 
# I believe A2 is the reference allele based on Plink manual https://www.cog-genomics.org/plink/1.9/data#ax_allele
# but it's confusing. the distinction between minor and non-ref is not well documented
# plink --bfile $plinkFileDir/$header \
# --freq \
# --missing \
# --within $clusters \
# --allow-extra-chr \
# --out $plinkFileDir/$header \
# --nonfounders \
# --keep-allele-order 
# 
# gzip -f $plinkFileDir/${header}.frq.strat


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
#python $scriptdir/plink2treemix.py $plinkFileDir/${header}.frq.strat.gz $treeFileDir/${header}.frq.strat.treemixFormat.gz
# this generates a file called: snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.frq.strat.treemixFormat.gz


############################ DON'T DO THIS --- RESULTS IN STRANGE OUTPUT XXXadd in FERRET OUTGROUPXXX ####################
# can then add in the ferret as an outgroup ; it's a bit hacky but seems to work:
# in the strat file: A1 is the alternate allele, and A2 is the reference allele
# this is based on --keep-allele-order 
# and the "MAC" field is not in fact the minor allele count in this case
# it is the A1 (alternate) count 
# so then when it's converted into treemix format
# the format is "MAC",Total-"MAC" where here MAC is the MAC field from the .strat. file, which is not actually MAC but is the A1 count.
# so the treemix format is alt,ref (instead of minor, major; that is why you can have entries like 40,0)
# this is really annoying inconsistent nomenclature
# But it means I can easily add in ferret. going to have ferret be fixed for the reference in all cases which would be 0,2 (alt,ref) for all sites
# in real life, the ferret would be heterozygous at some sites, but we are treating it all as ancestral-fixed.

# make a header
#infile=${header}.frq.strat.treemixFormat.gz
#outfile=${header}.frq.strat.InclFerret.NewRefFmt.treemixFormat

# add ferret to the header as the last group  ; use ofs="" so two spaces don't get inserted between line ending and ferret
#zcat $treeFileDir/$infile | head -n1 | awk '{OFS="";print $0,"FERRET"}' > $treeFileDir/$outfile

# then add the 0,2 count to the end of each line
# going to exclude the header line using grep -v BAJ (doesn't have anything to do with Baja specifically, that's just the start of the header line)
# tried using 0,1 but got weird ferret migration
# instead doing 0,2 (as though ferret is hom ref for all sites)
#zcat $treeFileDir/$infile | grep -v BAJ | awk '{OFS="";print $0,"0,2"}' >> $treeFileDir/$outfile


# then gzip this file to be ready for Treemix
#gzip -f $treeFileDir/$outfile 