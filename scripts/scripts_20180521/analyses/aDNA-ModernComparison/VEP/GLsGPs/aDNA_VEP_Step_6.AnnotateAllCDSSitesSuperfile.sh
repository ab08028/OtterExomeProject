#!/bin/bash
#$ -l h_rt=24:00:00,h_data=3G
#$ -N VEPannotateSuperfile
#$ -cwd
#$ -m bea
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -M ab08028

# want to annotate my cds superfile 
source /u/local/Modules/default/init/modules.sh
module load bedtools
# directories: 
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/

SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison

ref=mfur # only work with MFUR when going into VEP!
#dates="20190524-highcov-AFprior 20190524-lowcov-AFprior"
dates="20190701-highcov-AFprior-MajorMinor4 20190701-lowcov-AFprior-MajorMinor4"

#categories="synonymous missense stopgained"

#  20190524-highcov-UNIFprior 20190524-lowcov-UNIFprior" # set of angsdDates you want to process 
basename=angsdOut.mappedTo${ref}
type="GPs" # for now
for angsdDate in $dates
do
indir=$wd/VEP/$angsdDate 
GLdir=$wd/angsd-GLs/$angsdDate 
mkdir -p $GLdir/cdsPerCategoryFromVEP
### convert vep output back to bed coords
# note vep is 1based, bed is 0based
## so you want to take the cds superfile (that you got in step XX -- look up!)
# and you want to annotate it based on the VEP output (but first convert to a bed )
cdsSuperfile=$GLdir/${basename}.superfile.${type}.mafs.counts.cdsOnly.0based.bed.gz
vepWholeOutput=$indir/angsdOut.mappedTomfur.superfile.GPs.mafs.counts.cdsOnly.1based.VEPInput.VEP.output.pick.tbl.gz

# convert vep output into bed format:
zcat $vepWholeOutput | grep -v "#"  | awk '{OFS="\t";split($2,pos,":");print pos[1],pos[2]-1,pos[2],$1,".",".",".",".",".",".",".",".",$0}' | sed 's/\t$//g' > ${vepWholeOutput%.tbl.gz}.0based.Annotated.bed
gzip -f ${vepWholeOutput%.tbl.gz}.0based.Annotated.bed
# then want to intersect it with cdsSuperfile but *keep* all superfile sites -- if no annotation just put no annotation (let's see how to do this....) -- might make an ugly table
# put in whole line of a and of b
# wao:  	Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0. Restricted by -f and -r.
# I think this is what I want; a null feature of B if there isn't an overlap (those will be classed as "other" when I sum stuff up)
# need to get a header
header1=`zcat $cdsSuperfile | grep "#" -m1`
header2=`zcat $vepWholeOutput | grep -v "##" | grep "#" -m1`

comboheader=`echo -e "${header1}\t#chrom\tstart0based\tend\tmarkerID\tempty5\tempty6\tempty7\tempty8\tempty9\tempty10\tempty11\tempty12\t${header2}\tOverlap"` # need the "" and -e to get the tabs in
echo -e "$comboheader" | sed 's/\t\t/\t/g'>  $GLdir/cdsPerCategoryFromVEP/${basename}.superfile.${type}.mafs.counts.0based.allCDSSites.AnnotatedWithVEP.bed
bedtools intersect -a $cdsSuperfile -b ${vepWholeOutput%.tbl.gz}.0based.Annotated.bed.gz -wao >> $GLdir/cdsPerCategoryFromVEP/${basename}.superfile.${type}.mafs.counts.0based.allCDSSites.AnnotatedWithVEP.bed
gzip -f $GLdir/cdsPerCategoryFromVEP/${basename}.superfile.${type}.mafs.counts.0based.allCDSSites.AnnotatedWithVEP.bed

done
# use awk to pull out the stuff I want 
# want:
# coordinates (0 and 1 based)
# marker ID
# GPs per individual
# Annotation label -- were they "picked"? figure out how to deal with multi annotations in python
# and whether it was a null annotation
# put them all in python and  then if it's syn or mis you count it (or figure out how to do multis)
# and then all the rest put as 'other'
# will have to check on CANONICAL=YES and all the other filter_vep stuff if I do stuff in python 
