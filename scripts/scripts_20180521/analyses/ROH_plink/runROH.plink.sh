### sandbox, roh in plink
module load plink
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/20181119_filtered
header=snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants
vcf=${header}.vcf.gz
# I think this is a good vcf. It doesn't have admixed or relatives present,
# has a 75% excess het filter already applied
# and has max no-call frac of 20% (snp9b has 100% max no call frac aka no missingness filter)
## need clusters (same as generating treemix input)


wd=/u/flashscratch/a/ab08028/captures/analyses/ROH_plink
mkdir -p $wd
cd $wd
# need to convert bed to ped
# got bed (plink not UCSC) files in steps leading into Fis (merge those scripts if this ends up being useful)
# recoding vcf to ped/map files: (don't want bed)
plink --vcf $vcfdir/$vcf \
--allow-extra-chr \
--const-fid \
--out $wd/$header \
--keep-allele-order \
--recode \
--set-missing-var-ids @-#

# then want to separate into pops:
#################### make clusters file (once) ############
# make clusters file (ONCE)
# format of cluster file is FID\tIndId\tClusterID
# but fam id is 0:
# awk '{OFS="\t";print 0, $2}' $wd/$header.fam > $wd/$header.clusters
# THEN MANUALLY ADD POPS! ### <---- 

clusters=$wd/$header.clusters 

########### separate into clusters ##############

#for cluster in CA AL AK BAJ COM KUR
for cluster in AL BAJ COM KUR
do
mkdir -p $cluster
plink --file $wd/$header \
--allow-extra-chr \
--const-fid \
--keep-allele-order \
--within $clusters \
--keep-cluster-names $cluster \
--out $wd/$cluster/plinkFormat.$cluster \
--recode
done

####### run ROH ###########

for cluster in CA AL AK BAJ COM KUR
do

plink \
--file $wd/$cluster/plinkFormat.$cluster \
--allow-extra-chr \
--homozyg \
--homozyg-window-het 3 \
--homozyg-window-missing 5 \
--homozyg-kb 500 \
--homozyg-snp 50 \
--out $wd/$cluster/plink.ROH.$cluster
done