######### can be run in the shell #############
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

# input beagle GL or GP file 
# don't want the genotypes (will izip them later)
# just want default VEP format which is:

# Chr Start Stop Allele (ref/alt) Strand Identifier(optional)

# I think I want strand to always be +, otherwise it does rev comp of provided alleles.... 
# I want the identifier to match the "marker" of the beagle file, which is Chr_Position
# note that in mafs file the sites are in the exact same order as the beagle files (checked it using number of lines -- same line count in both, and in same order)
# but the identifier will help as an internal check (in python script, insist that marker == marker otherwise break)
# So the columns I want from my maf file is:

# chromo	position	major	minor	ref	knownEM	nInd
# will be $1, $2, $2, $5"/"$

# note $2, $2 because it's a single snp, 1-based inclusive, DO NOT have it be $2+1, that would be an insertion), 

# Okay so you don't want to use "major/minor" you want to do "ref/alt" so need to check which allele != the ref and use that one
# that's a little trickier; because based on Orlando command of doMajorMinor 1 I am inferring the major allele from the data, so it's not necessarily the ref allele
# 
# will be easiest conversion from the MAF file from ANGSD, rather than a beagle format (because don't have to convert 0123 to ATGC)
ref=mfur # only works with MFUR for VEP !

for angsdDate in XX YY # fill in dates
do
postDir=$GLdir/$angsdDate/posteriorProbabilities # location of your posterior probs
outdir=$wd/VEP/$angsdDate
# work on mafs file (result of angsd -doMajorMinor 1, so note that major allele != referenence allele, it's okay because I account for that below with awk
mafs=angsdOut.mappedTo${ref}.OrlandoSettings.mafs.gz
output=angsdOut.mappedTo${ref}.OrlandoSettings.VEP.inputFormat.txt
# need to ignore header "chromo"
# awk structure: if the major allele  ($3)== reference allele ($5), then minor allele is alt allele so do ref/alt = $5/$4
# but if the minor allele ($4) == ref allele ($5) then major allele is alt allele! and so you must do $5/$3: 
# THIS WORKS VERY NICELY :)
zcat $postDir/$mafs | grep -v "chromo" | awk '{OFS="\t"; if($3==$5) print $1,$2,$2,$5"/"$4,"+",$1"_"$2; else if($4==$5) print $1,$2,$2,$5"/"$3,"+",$1"_"$2 }' > $outdir/$output

# gzip output: 
gzip -f $outdir/$output
done