#! /bin/bash
#$ -cwd
#$ -l h_rt=40:00:00,h_data=8G,highp
#$ -m abe
#$ -M ab08028
#$ -N neutralCDSBedtools
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

######## get cds and neutral regions from super bed file that has all angsd output results in its extra columns ###########
source /u/local/Modules/default/init/modules.sh
module load bedtools
######################### MFUR coordinates ####################################
ref=mfur
# these ferret neutral regions are from the full modern dataset; eventually could detect neutral regions in my full angsd dataset as well; for now just want to compare to the modern dataset
neutBed=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/regionCoordinates/fromModernFullDataSet/bedCoords/all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.useThis.bed
# these ferret cds regions are from the genome annotation (gff)
cdsBed=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/MusPutFuro1.0.91.cdsOnly.0based.sorted.merged.bed
wd=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison

#dates="20190524-highcov-AFprior 20190524-lowcov-AFprior 20190524-highcov-UNIFprior 20190524-lowcov-UNIFprior" # set of angsdDates you want to process 
dates="20190612-highcov-pseudoHaps 20190612-lowcov-pseudoHaps"

for angsdDate in $dates
do
echo $angsdDate
indir=$wd/angsd-pseudoHaps/$angsdDate
basename=angsdOut.mappedTo${ref}
hapSuperfile=$indir/${basename}.pseudoHaps.superfile.0based.bed.gz
# bedtools intersect: 
# -wa says to print every entry of "a" that intersects with "b"; -header says to print the header of "a" before the results
# note that bed headers must start with "#"
### GPs:
echo "starting intersection"
bedtools intersect -a $hapSuperfile -wa -header -b $neutBed | gzip > ${hapSuperfile%.0based.bed.gz}.neutralOnly.0based.bed.gz

bedtools intersect -a $hapSuperfile -wa -header -b $cdsBed | gzip > ${hapSuperfile%.0based.bed.gz}.cdsOnly.0based.bed.gz


done


######################### ELUT coordinates ####################################
# if you want target beds they are in: /Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/ExomeCapture_INFO/FinalBEDS_round5ExomeCapture_Nov17_copiedfromAnnotationFolder/finalFiles-0based-Correction_20170906_USEFROMNOWON
# but note they are exons not CDS (can get cds from .gff file); and aren't based on 99-dedup genome
# So really need to do a bit of work to get them in shape that isn't worth it right now 
# because don't really want / need them for right now. don't have VEP database, phyloP scores, etc.
# so just subsetting ferret files, at least for now. 
