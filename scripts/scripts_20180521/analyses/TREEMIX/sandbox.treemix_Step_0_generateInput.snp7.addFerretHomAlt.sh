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

header=snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants
#header=snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants


# must add a snp name to my .bim file:
# make a backup
/bin/cp $plinkFileDir/$header.bim $plinkFileDir/${header}.original.bim
awk '{$2 = $1":"$4; print}' $plinkFileDir/${header}.original.bim >  $plinkFileDir/$header.bim
awk '{print $1,$2}' $plinkFileDir/$header.fam > $plinkFileDir/$header.samples
/bin/cp $plinkFileDir/$header.samples $plinkFileDir/$header.population.clusters 
clusters=$plinkFileDir/$header.population.clusters # 3 columns: 0 sampleID popID
# make a back up once you've added populations
/bin/cp $clusters $clusters.backup


##################################### With separate CA/BAJA and removing 3 relatives   ##############################
############# Excluding 3 Relatives: can use plink --filter ######################

# make a filter file:
# plink --file data --filter myfile.raw 1 --freq

#implies a file myfile.raw exists which has a similar format to phenotype and cluster files: that is, the first two columns are family and individual IDs; the third column is expected to be a numeric value (although the file can have more than 3 columns), and only individuals who have a value of 1 for this would be included in any subsequent analysis or file generation procedure. e.g. if myfile.raw were 
# copy cols 1, 2 and then a "1"

awk '{print $1,$2,1}' $plinkFileDir/$header.fam > $plinkFileDir/$header.exclList.allIndsIncluded
echo " YOU MUST MANUALLY EXCLUDE THE RELATIVES HERE" #### 
# excluding relatives from my excel sheet:
# 104_Elut_KUR_7 (keep 79_Elut_KUR_17)
# 77_Elut_KUR_14 (keep 81_Elut_KUR_2)
# 106_Elut_AL_AD_GE91109 (keep 118_Elut_AL_AD_GE91101)

###### MANUALLY edit and rename as : $plinkFileDir/$header.exclList.rmRelatives
# manually edit the file to set individuals you want to exclude to "0"
# tried: a2 vs a1 (doesn't make a diff); adding 4 inds (doesn't make diff); try no missing data
clusters=$plinkFileDir/$header.population.clusters.sepCA-BAJ # 3 columns: 0 sampleID popID
marker="sepCA-BAJ.exclRelatives-testADDINGMFUR-refA2"
plink --bfile $plinkFileDir/$header \
--freq \
--missing \
--within $clusters \
--allow-extra-chr \
--out $plinkFileDir/$header.${marker} \
--nonfounders \
--keep-allele-order \
--filter $plinkFileDir/${header}.exclList.rmRelatives 1 \
--a2-allele $plinkFileDir/${header}.ReferenceAlleles.txt

gzip -f $plinkFileDir/${header}.${marker}.frq.strat

python $scriptdir/plink2treemix.py $plinkFileDir/${header}.${marker}.frq.strat.gz $treeFileDir/${header}.${marker}.frq.strat.treemixFormat.gz

######## ADD MFUR ###########
# make a header
infile=${header}.${marker}.frq.strat.treemixFormat.gz
outfile=${header}.${marker}.frq.strat.treemixFormat.InclFERRET

# add ferret to the header as the last group  ; use ofs="" so two spaces don't get inserted between line ending and ferret
zcat $treeFileDir/$infile | head -n1 | awk '{OFS="";print $0,"FERRET"}' > $treeFileDir/$outfile

# I forced A1 to be ref allele using  # --reference-allele $plinkFileDir/${header}.ReferenceAllelesToSetAsA1.txt
# pattern in treemix for all is A1,A2 so want 2,0 for ferret: 
# so ferret should be all reference 2,0
zcat $treeFileDir/$infile | grep -v BAJ | awk '{OFS="";print $0,"0,8"}' >> $treeFileDir/$outfile

gzip $treeFileDir/$outfile
######## Try running tree mix on it : #######
######### CA Only (no Baja, no relatives) #############
k=500
wd=/u/flashscratch/a/ab08028/captures/analyses/TREEMIX/$genotypeDate/snp7 # update if using snp7 8 etc 
infile=$outfile.gz # this has CA-BAJ separated 
root='FERRET'
for m in {0..5}
do
outdir="root.${root}.mig.${m}.k.${k}.global.${marker}.treemix"
mkdir -p $wd/$outdir
treemix -i $treeFileDir/$infile -m ${m} -root ${root} -k ${k} -global -o $wd/$outdir/$outdir
done

######### experiment: try removing sites with missing data:
zcat $treeFileDir/$outfile | grep -v -w "0,0" > $treeFileDir/${outfile%.gz}.NOSITESMISSINGFROMWHOLEPOP
gzip $treeFileDir/${outfile%.gz}.NOSITESMISSINGFROMWHOLEPOP # removes about 400 sites

k=500
wd=/u/flashscratch/a/ab08028/captures/analyses/TREEMIX/$genotypeDate/snp7 # update if using snp7 8 etc 
infile=${outfile%.gz}.NOSITESMISSINGFROMWHOLEPOP.gz # this has CA-BAJ separated 
root='FERRET'
for m in {0..5}
do
outdir="root.${root}.mig.${m}.k.${k}.global.${marker}.MissingSitesRemoved.treemix"
mkdir -p $wd/$outdir
treemix -i $treeFileDir/$infile -m ${m} -root ${root} -k ${k} -global -o $wd/$outdir/$outdir
done

######## another experiment: less ld pruning #############

k=50
wd=/u/flashscratch/a/ab08028/captures/analyses/TREEMIX/$genotypeDate/snp7 # update if using snp7 8 etc 
infile=${outfile%.gz}.NOSITESMISSINGFROMWHOLEPOP.gz # this has CA-BAJ separated 
root='FERRET'
for m in {0..5}
do
outdir="root.${root}.mig.${m}.k.${k}.global.${marker}.MissingSitesRemoved.treemix"
mkdir -p $wd/$outdir
treemix -i $treeFileDir/$infile -m ${m} -root ${root} -k ${k} -global -o $wd/$outdir/$outdir
done