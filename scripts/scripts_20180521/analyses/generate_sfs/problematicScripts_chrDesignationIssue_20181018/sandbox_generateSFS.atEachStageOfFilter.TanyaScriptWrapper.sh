#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=8G,highp
#$ -N generate_sfs_allStages_allPops
#$ -o /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -e /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -m abe
#$ -M ab08028

### This is for the neutral SFS, but can modify choice of bed file to make coding SFS
############## THIS SCRIPT has issues
# problem: 
# during each of these filtering levels, not all sites have been called for all individuals
# instead, should make a miniature vcf file with some %% of the neutral bed file (1/4?)
# per population with all sites called
# and make an sfs from that
# that should be relatively easy to script (maybe haha)
# but will take some time to figure out/ run
source /u/local/Modules/default/init/modules.sh
module load python/2.7

rundate=20180806
### wrapper for Tanya's script
SCRATCH=/u/flashscratch/a/ab08028
popHeaders=$SCRATCH/captures/samples/popHeaderFiles/popHeaders_rmBadInd/

# location of per-population input VCFs (with NO no-call genotypes)
vcfdir=$SCRATCH/captures/vcf_filtering/${rundate}_filtered
neutralBed=${vcfdir}/bedCoords/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.useThis.bed

# output SFS location
SFSdir=$SCRATCH/captures/analyses/SFS/${rundate}
mkdir -p $SFSdir/neutralSFS

# location of tanya's scripts
# this is latest 20180822 where it runs faster (but need to give enough memory)
# and doesn't require pre-filtered SFS 
tanyaDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/scripts_from_others/tanya_scripts/

vcfs="all_1_TrimAlt_raw_variants.vcf.gz all_5_passingFilters_raw_variants.vcf.gz all_6_rmBadIndividuals_passingFilters_raw_variants.vcf.gz all_7_passingBespoke_maxNoCallFrac_0.2_rmBadIndividuals_passingFilters_raw_variants.vcf.gz"

# California:
echo "California"

CA_ind=$popHeaders/california.rmBad.txt
for vcf in $vcfs
do
echo $vcf
# if vcf file is not prefiltered to just contain neutral (or other) regions

python $tanyaDir/popgen_tools/popgen_tools.py \
--vcf_file $vcfdir/$vcf \
--sfs_out $SFSdir/neutralSFS/CA_${vcf%.vcf.gz}.sfs.out \
--no_pi \
--target_bed $neutralBed \
--names_list $CA_ind

done


# Kurils:
echo "Kuril"
KUR_ind=$popHeaders/kuril.rmBad.txt
for vcf in $vcfs
do
echo $vcf
# if vcf file is not prefiltered to just contain neutral (or other) regions

python $tanyaDir/popgen_tools/popgen_tools.py \
--vcf_file $vcfdir/$vcf \
--sfs_out $SFSdir/neutralSFS/KUR_${vcf%.vcf.gz}.sfs.out \
--no_pi \
--target_bed $neutralBed \
--names_list $KUR_ind

done

# Commanders:
echo "Commanders"
COM_ind=$popHeaders/commanders.rmBad.txt
for vcf in $vcfs
do
echo $vcf
# if vcf file is not prefiltered to just contain neutral (or other) regions

python $tanyaDir/popgen_tools/popgen_tools.py \
--vcf_file $vcfdir/$vcf \
--sfs_out $SFSdir/neutralSFS/COM_${vcf%.vcf.gz}.sfs.out \
--no_pi \
--target_bed $neutralBed \
--names_list $COM_ind

done

# Alaska:
echo "Alaska"
AK_ind=$popHeaders/alaska.rmBad.txt
for vcf in $vcfs
do
echo $vcf
# if vcf file is not prefiltered to just contain neutral (or other) regions

python $tanyaDir/popgen_tools/popgen_tools.py \
--vcf_file $vcfdir/$vcf \
--sfs_out $SFSdir/neutralSFS/AK_${vcf%.vcf.gz}.sfs.out \
--no_pi \
--target_bed $neutralBed \
--names_list $AK_ind

done

# Aleutian:
AL_ind=$popHeaders/aleutian.rmBad.txt
for vcf in $vcfs
do
echo $vcf
# if vcf file is not prefiltered to just contain neutral (or other) regions

python $tanyaDir/popgen_tools/popgen_tools.py \
--vcf_file $vcfdir/$vcf \
--sfs_out $SFSdir/neutralSFS/AL_${vcf%.vcf.gz}.sfs.out \
--no_pi \
--target_bed $neutralBed \
--names_list $AL_ind

done