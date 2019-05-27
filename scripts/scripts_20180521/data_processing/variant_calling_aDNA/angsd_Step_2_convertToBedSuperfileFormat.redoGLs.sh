#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G
#$ -m abe
#$ -M ab08028
#$ -N angsdStep2SuperFile
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

####### make a super angsd output file ##########
# this script takes about 16 hours to run  -- could submit each date separately to speed it up #
# want to combine beagle gprobs or gls with maf and count information
#then 
# Convert --> bed format (with all this info as extra columns) ############
# can then intersect with neutral bed file, cds bed file, convert to VEP format etc

######### dirs and files ###########
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
dates="20190524-highcov-AFprior 20190524-lowcov-AFprior 20190524-highcov-UNIFprior 20190524-lowcov-UNIFprior" # set of angsdDates you want to process 
for ref in mfur elut
do
echo $ref
for angsdDate in $dates
do
echo $angsdDate

indir=$wd/angsd-GLs/$angsdDate

basename=angsdOut.mappedTo${ref}
mafs=$indir/${basename}.mafs.gz
counts=$indir/${basename}.counts.gz
GPs=$indir/${basename}.beagle.gprobs.gz
GLs=$indir/${basename}.beagle.gz
GPoutput=$indir/${basename}.superfile.GPs.mafs.counts.0based.bed
GLoutput=$indir/${basename}.superfile.GLs.mafs.counts.0based.bed

####################################################
# okay so should paste mafs first because they are easier to deal with -- note that $GPs will differ in numbering of columns depending on number of individuals so that's not good
# so then you use awk to convert the first couple columns into bed format (scaff, 0based start,nonIncEnd,MarkerName (Scaff_Position))
# and then put in 8 "." empty columns for the rest of bed format
# bed columns > 12 are ignored, so can store all the rest of the info there
# so in the output, columns 1-12 are bed format; then come mafs information, then GP or GL information (# of cols will vary with individuals), then counts info (# of cols with vary with individuals)
# have to use sed to get rid of a terminal tab at the end of the line that bed tools isn't happy about (good tip for future: cat -t will show invisibles)
# THIS WORKS BECAUSE ALL 3 files are in the EXACT SAME ORDER when output from angsd (rows correspond to same set of sites, columns to same individuals)


# order of columns in the superfile coming out of paste:
# marker	allele1	allele2	Ind0	Ind0	Ind0	Ind1	Ind1	Ind1	Ind2	Ind2	Ind2	Ind3	Ind3	Ind3	Ind4	Ind4	Ind4	Ind5	Ind5	Ind5	Ind6	Ind6	Ind6	Ind7	Ind7	Ind7	Ind8	Ind8	Ind8	chromo	position	major	minor	ref	knownEM	nInd	ind0TotDepth	ind1TotDepth	ind2TotDepth	ind3TotDepth	ind4TotDepth	ind5TotDepth	ind6TotDepth	ind7TotDepth	ind8TotDepth
# want to make a bed format, which means we need specific information in the first few columns:
# dummy bed line:
# GL896898.1	752354	752355	GL896898.1_19406	.	.	.	.	.	.	.	1	0	0.999984	0.000016	0.000000	0.999984	0.000016	0.000000	0.999984	0.000016	0.000000	0.999984	0.000016	0.000000	0.999984	0.000016	0.000000	0.999984	0.000016	0.000000	0.999992	0.000008	0.000000	0.999984	0.000016	0.000000	0.999984	0.000016	0.000000

# GPs:
# Want to have the header starting with a "#" in bed file 
# note don't have to do for GPs and GLs separately; have same header
angsdheaders=`paste <(zcat $mafs | head -n1) <(zcat $GPs | head -n1) <(zcat $counts | head -n1)` # this is the headers of all the files pasted together
# so want to add a bed header too, which will be 
bedhead="#chrom\tstart0based\tend\tmarkerID\tempty5\tempty6\tempty7\tempty8\tempty9\tempty10\tempty11\tempty12"
comboheader=`echo -e "$bedhead\t$angsdheaders"` # need the "" and -e to get the tabs in
#echo -e "$comboheader" >  ${GPoutput} 
#paste <(zcat $mafs) <(zcat $GPs) <(zcat $counts) | grep -v "chromo" | awk '{OFS="\t";print $1,$2-1,$2,$1"_"$2,".",".",".",".",".",".",".",".",$0}' | sed 's/\t$//g' >> ${GPoutput} # go into awk and rearrange to make it bed format with extra columns 
#gzip -f ${GPoutput}
# GLs: 
echo -e "$comboheader" >  ${GLoutput} 
paste <(zcat $mafs) <(zcat $GLs) <(zcat $counts) | grep -v "chromo" | awk '{OFS="\t";print $1,$2-1,$2,$1"_"$2,".",".",".",".",".",".",".",".",$0}' | sed 's/\t$//g' >> ${GLoutput} # go into awk and rearrange to make it bed format with extra columns 
gzip -f ${GLoutput} 
done
done

############################ SANDBOX #######################################


#sandbox=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/sandbox/superfileSandbox
#mafs=$sandbox/mafs.txt
#counts=$sandbox/counts.txt
#GPs=$sandbox/grobs.txt