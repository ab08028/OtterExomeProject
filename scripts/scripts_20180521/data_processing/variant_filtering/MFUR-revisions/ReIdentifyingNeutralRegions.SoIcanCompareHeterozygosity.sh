################################# Set up ###########################

###### NOTE for the future: bedtools intersect -v will get rid of the entire entry in A that intersects by even 1bp with a region in B
# for my capture data, this was okay because my regions in A were small and I didn't want them to be adjacent to any repeat or CpG region
# but for whole genome data this could perhaps be too stringent and result in too much lost data 
# so many adjust in future studies. But it's okay for the sea otters.

# modules
source /u/local/Modules/default/init/modules.sh
module load bedtools
module load blast

rundate=20181119
maxHetFilter=0.75 # het filter used across all samples (per population het filter occurs during easy sfs)

captures=~/klohmueldata/annabel_data/captures/
filterDir=$captures/vcf_filtering/${rundate}_filtered/
bedDir=$captures/vcf_filtering/${rundate}_filtered/bedCoords
wd=$captures/vcf_filtering/${rundate}_filtered/checkingNeutralSites-Revisions
mkdir -p $wd
mfurDir=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome
drerDir=/u/home/a/ab08028/klohmueldata/annabel_data/zebra_fish_genome
REFERENCE=$mfurDir/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
gff=$mfurDir/Mustela_putorius_furo.MusPutFur1.0.91.gff3
cdsRegions=$mfurDir/MusPutFuro1.0.91.cdsOnly.0based.sorted.merged.bed

exonicRegions=$mfurDir/MusPutFur1.0.91.exonsOnly.0based.sorted.bed # make this below if haven't already

mfurCpG=$mfurDir/CpG_Islands/cpgIslandExtUnmasked.already0based.bed 
mfurRepeat=$mfurDir/repeatMasker_UCSC/repeatRegions.0based.bed

# script to get gc content
getGC=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/scripts_from_others/GCContent/get_gc_content.pl

######## Filtering steps:
# 1. >10kb from exons
# 2. not inside CpG Island
# 3. not inside repeat region
# 4. normal GC content
# 5. doesn't blast to zebra fish
####### get exonic regions (once) and make sure is sorted: ############### 
### (do once) grep exon $gff | awk '{OFS="\t";print $1,$4-1,$5,$9}' | sort -k1,1 -k2,2n > $exonicRegions
# should be ~200,000 lines

############ HQ site coords ####################
# results of filtering snps (all populations; all nv and snps)
noCallFrac=1.0 # no filter at all. no limits on how many individuals must be called. but there is min dp 500 at site level
prefix=all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_${noCallFrac}_rmBadIndividuals_passingFilters
# all8 is fine for this -- all 9 just has gone through max het filter so a few sites missing, but as long as you work from all9 to pull the sites, using the all8 coords is no problem

# prefix=all_9_maxHetFilter_${maxHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters

hqSites=$bedDir/${prefix}.sorted.merged.coords.bed # bed coords (sorted, merged) of sites from step 7 of filtering (comes from filtering_Step_2.sh)

############# Get distance of every set of sites from exonic regions in ferret genome ################
mkdir -p $wd/distanceFromExons # this dir will have info on distance of sites from exons
mkdir -p $wd/CpG_Islands
mkdir -p $wd/repeatRegions
mkdir -p $wd/get_fasta
mkdir -p $wd/GC_Content
mkdir -p $wd/zebra_fish
mkdir -p $wd/passing_sites
### for revisions only:
mkdir -p $wd/sitesThatDidntPass

 # this dir will have neutral regions going through 3 checks: CpG Islands, GC content, and blast to fish
bedtools closest -d -a ${hqSites} -b ${exonicRegions} > $wd/distanceFromExons/${prefix}.distanceFromExons.0based.txt



#### NOTE: the output of this will be
# [Hq site info] [closest exon info] [distance between]; so I want the HQ sites that are >10000bp away from the closest exon.
# don't want it the other way around (getting info on each exon). want info on hq sites.

# last column (8) is the distance; want it to be at least 10,000, and want to keep
# track of the distance. Collect all that are >10,000 away. 
# pick the ones with high distance (awk) (get site totals below)
awk -F'\t' '{OFS="\t";if($8>10000)print $1,$2,$3}' $wd/distanceFromExons/${prefix}.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/distanceFromExons/${prefix}.min10kb.fromExon.0based.sorted.merged.bed

awk -F'\t' '{OFS="\t";if($8>100000)print $1,$2,$3}' $wd/distanceFromExons/${prefix}.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/distanceFromExons/${prefix}.min100kb.fromExon.0based.sorted.merged.bed

awk -F'\t' '{OFS="\t";if($8>1000000)print $1,$2,$3}' $wd/distanceFromExons/${prefix}.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/distanceFromExons/${prefix}.min1Mb.fromExon.0based.sorted.merged.bed


### for revision only: collect the ones that are <=10000 for heterozygosity measurement:
awk -F'\t' '{OFS="\t";if($8<=10000)print $1,$2,$3}' $wd/distanceFromExons/${prefix}.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/sitesThatDidntPass/${prefix}.lte10kb.fromExon.0based.sorted.merged.FAIL.bed

### Note: 1,2,3 columns are the HQ SITES position, NOT the position of the exon. (If you mess up what is a and b in bedtools closest this would be messed up)
######### get total amounts of sequence in each file: ########

> $filterDir/filteringStats/totalSequenceByDistanceFromExons.txt
for i in `ls $wd/distanceFromExons/*bed`
do
echo $i >> $filterDir/filteringStats/totalSequenceByDistanceFromExons.txt
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $i >> $filterDir/filteringStats/totalSequenceByDistanceFromExons.txt
done


###### for now use this 10kb bed file for making the SFS
###################### Check neutral regions (10kb)###################

######### Check to see if any regions intersect with CpG Islands ##############
# got CpG islands from UCSC browser;
# added .1 to each scaffold name so it matches my reference (see script in mfurDir/CpG_Islands)
# trying to intersect 
# get regions that DO NOT intersect with CpG islands
# *** -v **** this will output regions in "A" that ***DO NOT*** intersect with "B" (CpG Islands)
bedtools intersect -v -a $wd/distanceFromExons/${prefix}.min10kb.fromExon.0based.sorted.merged.bed -b $mfurCpG > $wd/CpG_Islands/${prefix}.min10kb.fromExon.noCpGIsland.0based.sorted.merged.bed
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' 

###### !!! for revisions: want to collect the coordinates that don't fit in:
# take out -v and so collect that sites that DO intersect with CpG island:
### SITES THAT ARE TOO CLOSE to CPG for reference DO NOT USE THESE SITES::::::
# note intersect and intersect -v don't really act as inverses of each other -- intersect will portion up regions, while intersect -v will remove the whole region -- so wont add up.
#bedtools intersect -a $wd/distanceFromExons/${prefix}.min10kb.fromExon.0based.sorted.merged.bed -b $mfurCpG > $wd/sitesThatDidntPass/${prefix}.TooCloseToCpG.0based.sorted.merged.FAIL.bed

########## 20181119 value:
#9,212,233 sites >10kb from genes (no max call frac) w no max no call frac filter
#9,062,287 after CpG filter (no max call frac filter)

########## 20180806 value
# total sequence before CpG filter: 9,030,072 w 0.9 max no call frac (7,818,882 when using 0.2 max no call frac)
# total sequence after CpG filter: 8,885,106 w 0.9 max no call frac  (6,816,729 when using 0.2 max no call frac)
# so lost ~100kb of sequence
# get CpG content
# overly conservative to filter >10kb from CpG Islands.
#bedtools closest -d -a $wd/distanceFromExons/all_7_passingBespoke.min10kb.fromExon.0based.sorted.merged.bed -b ${mfurCpG} > $wd/CpG_Islands/all_7_passingBespoke.min10kb.fromExon.distanceFromCpGIslands.0based.txt
#awk -F'\t' '{OFS="\t";if($8>10000)print $1,$2,$3}' $wd/CpG_Islands/all_7_passingBespoke.min10kb.fromExon.distanceFromCpGIslands.0based.txt > $wd/CpG_Islands/all_7_passingBespoke.min10kb.fromExon.min10kb.fromCpG.0based.bed
 

######### Check to see if any regions intersect with RepeatMasker Regions ##############
# got repeat mask from UCSC browser
# added .1 to each scaffold name so it matches my reference (see script in mfurDir/repeatMasker_UCSC)
# trying to intersect 
# get regions that DO NOT intersect with repeat regions

bedtools intersect -v -a $wd/CpG_Islands/${prefix}.min10kb.fromExon.noCpGIsland.0based.sorted.merged.bed -b $mfurRepeat > $wd/repeatRegions/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.merged.bed
# check amount  of sequence lost:
# awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'  $wd/repeatRegions/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.merged.bed
# 2018119 values: (no max call frac)
# 6,853,475 (so lost 2,208,812 sites (same as before))
# 20180806 values:
# 6,701,482 remain with 0.9 max no call frac (5,953,300 with 0.2 max no call frac) (lost 2Mb). Lost more with 0.9 setting at this stage: maybe because things that map poorly are in repeat regions more
#### !!! for revisions only: get list of sites that DO intersect with repeat regions -- DO NOT USE THESE SITES! Just getting for my recods
# remove -v to get the intersect:
# note intersect and intersect -v don't really act as inverses of each other -- intersect will portion up regions, while intersect -v will remove the whole region -- so wont add up.
#bedtools intersect -a $wd/CpG_Islands/${prefix}.min10kb.fromExon.noCpGIsland.0based.sorted.merged.bed -b $mfurRepeat > $wd/sitesThatDidntPass/${prefix}.IntersectWithRepeat.based.sorted.merged.FAIL.bed

### skipping drer as it removes so few loci.
 
########### So what I want to do is the following :
# just do intersect -v for the regions that are lte 10kb from genes and the regions in the final bed file
# from before:
finalBed=/u/home/a/ab08028/klohmueldata/annabel_data/captures/vcf_filtering/20181119_filtered/checkingNeutralSites/passing_sites/all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.useThis.bed
bedtools intersect -v -a $wd/distanceFromExons/${prefix}.min10kb.fromExon.0based.sorted.merged.bed -b $finalBed > $wd/RegionsThatDidntEndUpInFinalPassingSites.gt10kbFromExons.bed
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $wd/RegionsThatDidntEndUpInFinalPassingSites.gt10kbFromExons.bed
# total sequence that didn't make it in should be ~3Mb:
# 2,372,696
# final passing: 6839537
# sites gt 10kb from exons: 
# 9212233
### YES IT MATCHES! woohoo. 6,839,537