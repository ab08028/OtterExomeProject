#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G
#$ -m abe
#$ -M ab08028
#$ -N CDSsuperfile2VEP
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

## goal of this script is to take the mafs output of angsd and convert it into a tab-delim table that can be input into VEP
# giving the scaff, position, ref/alt alleles, strand and marker name for your sites and then they will be annotated by VEP 
# you will then be able to figure out how to merge this with your beagle files and sum up GLs or GPs across different GP categories
# a goal will be to get the VEP output in the exact same order as the input, but I'm not sure that will work. Can merge across marker IDs instead -- figure out how to do efficiently

# directories: 
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/

SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
GLdir=$wd/angsd-GLs

ref=mfur # only work with MFUR when going into VEP!
#dates="20190524-highcov-AFprior 20190524-lowcov-AFprior 20190524-highcov-UNIFprior 20190524-lowcov-UNIFprior" # set of angsdDates you want to process 
dates="20190701-highcov-AFprior-MajorMinor4 20190701-lowcov-AFprior-MajorMinor4"
basename=angsdOut.mappedTo${ref}
for angsdDate in $dates
do
indir=$GLdir/$angsdDate # location of your posterior probs
outdir=$wd/VEP/$angsdDate
mkdir -p $outdir
# can use GP or GL superfile; they are in the same order and have same sites, so it doesn't matter, since GLs and GPs aren't taken along for the ride in the VEP input format
# I am using GPs, but again, could be either:
cdsSuperfile=${basename}.superfile.GPs.mafs.counts.cdsOnly.0based.bed.gz
output=${cdsSuperfile%.0based.bed.gz}.1based.VEPInput.txt


########### convert to VEP input format ###############
# VEP input format:
# Chr Start Stop Allele (ref/alt) Strand Identifier(optional)
# I think I want strand to always be +, otherwise it does rev comp of provided alleles.... 
# I want the identifier to match the "marker" of the beagle file, which is Chr_Position

# Okay so you don't want to use "major/minor" you want to do "ref/alt" so need to check which allele != the ref and use that one
# that's a little trickier; because based on Orlando command of doMajorMinor 1 I am inferring the major allele from the data, so it's not necessarily the ref allele

# converting from the section of the superfile that contains the angsd "mafs" output -- easiest conversion
# could also just work directly from "mafs" output itself (if you didn't have to subset your data with bedtools); see commented code below
# awk structure: if the major allele  ($3 in maf file; $15 in superfile)== reference allele ($5 in maf file; $17 in superfile), then minor allele is alt allele so do ref/alt = $5/$4 (or $17/$16)
# but if the minor allele ($4 in maf file; $16 in superfile) == ref allele ($5 ; $17 in superfile) then major allele is alt allele! and so you must do $5/$3 (or $17/$15) 
# THIS WORKS VERY NICELY :)
### --> this commented out code works if you're going straight from a mafs file from angsd (maybe once generated with -rf to already be in the cds region)
#zcat $postDir/$mafs | grep -v "chromo" | awk '{OFS="\t"; if($3==$5) print $1,$2,$2,$5"/"$4,"+",$1"_"$2; else if($4==$5) print $1,$2,$2,$5"/"$3,"+",$1"_"$2 }' > $outdir/$output
# but this won't work on the subsetted cds superfile
# need a different structure, but same principle. Want to work on the MAF columns which come directly after the bed columns --> code below
# So you just need to add "12" to all the numbers:
zcat $indir/$cdsSuperfile | grep -v "#chrom" | awk '{OFS="\t"; if($15==$17) print $13,$14,$14,$17"/"$16,"+",$13"_"$14; else if($16==$17) print $13,$14,$14,$17"/"$15,"+",$13"_"$14}' > $outdir/$output
# NOTE: this ordering of columns won't change based on the number of individuals
# it will always be 12 bed columns, followed by the MAF file columns. The GL or GP and counts files (which DO vary based on # of inds) come afterward, so don't affect the count

# NOTE: the coordinates in the MAFS section of the superfile are 1-based
# If you were using the first couple columns of the superfile (which is the bed-file part of the file)
# Then you'd need to fuss around with 0 vs 1 based
# But since we are using the mafs section, which is 1based, and going into VEP format, which is 1-based, we are good to go! No need to add or subtract from position.


# gzip output: 
gzip -f $outdir/$output
done
