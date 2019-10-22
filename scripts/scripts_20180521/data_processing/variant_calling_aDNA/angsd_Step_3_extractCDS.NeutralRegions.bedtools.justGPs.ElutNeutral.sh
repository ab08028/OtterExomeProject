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

######################### ELUT coordinates ####################################
# if you want target beds they are in: /Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/ExomeCapture_INFO/FinalBEDS_round5ExomeCapture_Nov17_copiedfromAnnotationFolder/finalFiles-0based-Correction_20170906_USEFROMNOWON
# but note they are exons not CDS (can get cds from .gff file); and aren't based on 99-dedup genome
# So really need to do a bit of work to get them in shape that isn't worth it right now 
# because don't really want / need them for right now. don't have VEP database, phyloP scores, etc.
# so just subsetting ferret files, at least for now.

# 20191022: wanting to do elut neutral regions at least; wait on CDS for now. So have neutral bed from my capture design - do I want any overhang? not for now.
ref=elut
neutBed=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/neutralCoordsElut/finalNeutralSet.10K.1kb.regions.NoOverhang.NoFishMatch.NoNNNs.DupScaff99Removed.0-based.20170906.bed
wd=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison
dates="20190701-highcov-AFprior-MajorMinor4 20190701-lowcov-AFprior-MajorMinor4" # set of angsdDates you want to process 

for angsdDate in $dates
do
echo $angsdDate
indir=$wd/angsd-GLs/$angsdDate
basename=angsdOut.mappedTo${ref}
GPSuperfile=$indir/${basename}.superfile.GPs.mafs.counts.0based.bed.gz
#GLSuperfile=$indir/${basename}.superfile.GLs.mafs.counts.0based.bed.gz
# bedtools intersect: 
# -wa says to print every entry of "a" that intersects with "b"; -header says to print the header of "a" before the results
# note that bed headers must start with "#"
### GPs:
echo "starting GPs"
bedtools intersect -a $GPSuperfile -wa -header -b $neutBed | gzip > ${GPSuperfile%.0based.bed.gz}.neutralOnly.0based.bed.gz

# not doing cds for now
#bedtools intersect -a $GPSuperfile -wa -header -b $cdsBed | gzip > ${GPSuperfile%.0based.bed.gz}.cdsOnly.0based.bed.gz
# not doing GLs
done
