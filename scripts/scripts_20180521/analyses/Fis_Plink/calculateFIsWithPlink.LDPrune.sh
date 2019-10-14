# Experimenting with vcftools stats

# want to separate things by population

#module load vcftools 
module load plink
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/20181119_filtered
header=snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants
vcf=${header}.vcf.gz
# I think this is a good vcf. It doesn't have admixed or relatives present,
# has a 75% excess het filter already applied
# and has max no-call frac of 20% (snp9b has 100% max no call frac aka no missingness filter)
## need clusters (same as generating treemix input)


wd=/u/flashscratch/a/ab08028/captures/analyses/Fis_plink
cd $wd

# convert to plink format bed file : 
plink --vcf $vcfdir/$vcf \
--allow-extra-chr \
--const-fid \
--out $wd/$header \
--keep-allele-order \
--make-bed \
--set-missing-var-ids @-#

# or try with recoding it to ped/map files:
#plink --vcf $vcfdir/$vcf \
#--allow-extra-chr \
#--const-fid \
#--out $wd/$header \
#--keep-allele-order \
#--recode \
#--set-missing-var-ids @-#
# need --allow-extra-chr for non-standard chr names
# need --const-fid because there are underscores in my sample names (no family names)
# need keep-allele-order so that original vcf encoding (ref/alt) is kept
# need --set-missing-var-ids to make variant names instead of "." for LD pruning -- use @ to specify chr and # to specify bp position

#################### make clusters file (once) ############
# make clusters file (ONCE)
# format of cluster file is FID\tIndId\tClusterID
# but fam id is 0:
# awk '{OFS="\t";print 0, $2}' $wd/$header.fam > $wd/$header.clusters
# THEN MANUALLY ADD POPS! ### <---- 

clusters=$wd/$header.clusters

########### loop through clusters and separate infiles #############
for cluster in CA AL AK BAJ COM KUR
do
mkdir -p $cluster
plink --bfile $wd/$header \
--allow-extra-chr \
--const-fid \
--keep-allele-order \
--within $clusters \
--keep-cluster-names $cluster \
--out $wd/$cluster/plinkFormat.$cluster \
--make-bed \
--geno 0 ### trying this!!! removes 200K snps
done

# tried: filtering for missing data; didn't help. but when is it filtering? tried before and after; doesn't change
# what else could it be??? hwe? 
########### LD PRUNING ###############
# --indep-pairwise <window size>['kb'] <step size (variant ct)> <r^2 threshold>
# want like 50K variants left based on what I've read.
windowSize=100kb
variantStepSize=30 # not sure what this should be? starting with 30 because such low diversity? that is how much window shifts 
r2=0.2 # what I used when making the PCA
for cluster in CA AL AK BAJ COM KUR
do
# LD prune:
plink --bfile $wd/$cluster/plinkFormat.$cluster --within $clusters --allow-extra-chr --const-fid --indep-pairwise $windowSize $variantStepSize $r2 --out $wd/$cluster/plinkFormat.$cluster
# calculate F stat:
plink --bfile $wd/$cluster/plinkFormat.$cluster --het --within $clusters --allow-extra-chr --const-fid --extract $wd/$cluster/plinkFormat.$cluster.prune.in --out $wd/$cluster/plinkFormat.$cluster.ldpruned
### wtff now it's negative. Wahlund effect before with --within not actually calculating allele freqs per population??
done

### NOTE FOR LD PRUNING:
# MAKE SURE IDs of pruned  are real and not "." : 
# should look like:
# GL896898.1-23352
# GL896898.1-23947
# GL896898.1-38867
# GL896898.1-98062
# GL896898.1-99284
# Pruning complete.  1442723 of 1571275 variants removed. Leaves like 1.28K . See how this changes things.
# so the output is in: snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.LDPRUNED.prune.in
# then have the ".in" file be used to --extract to keep them in the analysis and exclude the others
# THIS DOESN"T WORK -- need variant id 



######## 

# try combining CA and BAJ to see if they still stand out 
# plink --bfile $header --het --within snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.BAJClassAsCA.clusters --allow-extra-chr --const-fid --out BAJ-CACombined

######## just pull out CA:
wd=
cluster=CA
mkdir $cluster
vcftools --gzvcf $vcfdir/$vcf --keep $cluster.txt --out $cluster/$cluster.only --recode

gzip $cluster/$cluster.only.recode.vcf


# [ DO REMOVE 75% HET SITES STEP]

### STILL NEED TO EVENTUALLY LD PRUNE, but for now just focus on basics since we saw LD pruning still yields neg Fis.
### download this vcf and play with it 
# Experimenting with vcftools:
#### just het filter (0.75)
vcftools --gzvcf testOut.maxHetFilterPerPop.0.75.vcf.gz --het
# doing it with all pops as one led to good spread of Fis values ; but that can just be Wahlund effect
### het filter + monomorphics filtered
vcftools --gzvcf testOut.maxHetFilterPerPop.0.75.vcf.gz --het --out hetFilterNoMonomorph

#vcftools --gzvcf $vcfdir/$vcf --het --keep CA.txt 
# hmmmm OOOH is it that monomorphic sites are appearing or something? could that be it? exclude those perhaps??? 

### Aha maf filter might be something important!! 
# maybe the monomorphs are being odd? try maf 0.05.
#vcftools --gzvcf $vcfdir/$vcf --maf 0.05 --keep CA.txt --recode
# but is that maf just within CA?? or is it among them all. figure that out. 