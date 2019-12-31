#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=2G
#$ -m abe
#$ -M ab08028
#$ -pe shared 16
#$ -N abbababaSingleInd
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd
######## ANGSD ABBA-BABA multipop
trimValue=7 # set value you want to trim from either end of read (looking at mapdamage plots)
#posterior=1 # setting for angsd -doPost : 1 for using allele frequencies as prior, 2 for using a uniform prior 
snpCutoff=1e-06
minInd=3 # no missing data (?)
#todaysdate=`date +%Y%m%d`'-highcov-AFprior-MajorMinor4'
#todaysdate='20191205-highcov-AFprior-MajorMinor4'
todaysdate=`date +%Y%m%d`'-highcov-singleInd-minInd-ancCA-AK-KUR'
#### ANGSD v 0.923 ####
source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env

######### dirs and files ###########
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
bamdir=$wd/bams/
ABBADIR=$wd/ABBA_BABA
mkdir -p $ABBADIR
mkdir -p $ABBADIR/$todaysdate
outdir=$ABBADIR/$todaysdate


# abcDstat2.cpp:
# 	-doAbbababa2	                0	run the abbababa analysis
# 	-rmTrans		        0       remove transitions
# 	-blockSize		       5000000	size of each block in bases
# 	-anc			       (null)	fasta file with outgroup
# 	-sample			        0	sample a single base in each individual
# 	-maxDepth		        100	max depth of each site allowed
# 	-sizeFile		       (null)   file with sizes of the populations	
# 	-enhance			0	only analyze sites where outgroup H4 is non poly
# 	-Aanc			        0	set H4 outgroup allele as A in each site
#         -useLast                        0       1=use the last group of bam files as outgroup

# bam list needs to be in order
# with a size file 
bamListDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/forABBABABA/singleInd
mfurBamList=$bamListDir/angsd.bamList.mappedtoMfurfullpaths.HighCovPlusADNAOnly.InPopOrder.3inds.ancCA.AK.KUR.txt

mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
### HAS TO BE MFUR TO HAVE ANCESTRAL!!! #########
############### mfur high coverage ##############
angsd -nThreads 16 \
-anc $mfurRef \
-ref $mfurRef \
-bam $mfurBamList \
-doAbbababa 1 \
-rmTrans 1 \
-blockSize 5000000 \
-doCounts 1 -remove_bads 1 -uniqueOnly 1 \
-C 50 -baq 1 -trim 7 -minQ 20 -minMapQ 25 \
-out $outdir/angsdOut.TransvOnly \
-minInd $minInd

source deactivate

