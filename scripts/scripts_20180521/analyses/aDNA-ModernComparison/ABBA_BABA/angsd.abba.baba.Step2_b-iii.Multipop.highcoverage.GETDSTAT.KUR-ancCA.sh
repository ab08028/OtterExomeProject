#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=2G
#$ -m abe
#$ -M ab08028
#$ -pe shared 16
#$ -N abbababa
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd
######## Then run DSTAT 

# start with output of abbababa2: 
#The output file is bam.Angsd.abbbababa2 (used for the 4-population test) Each line represents a block with a chromsome name (Column 1) for one of the possible 30 trees (so each block is written on 30 lines), a start position (Column 2), an end postion (Column 3). Columns 4,5 and 6 are the numerator, denominator and number of sites analyzed in the block. The next 256 columns are the counted patterns of alleles in the tree, starting from X0000=AAAA,X0001=AAAC,....,X3333=TTTT, with the correspondence 0=A,1=C,2=G,3=T. This file is used as input for the R script estAvgError.R.

#We run the R script specifying the error files for the population with 3 individuals. This is done defining the error files in each populations inside a text file (including a line for the outgroup population). If a population has no error file, it is sufficient to write NA. Create a file called errorList.error with written 

### ANGSD v 0.923 ####
source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
module load R
source activate angsd-conda-env # activate conda env

######### dirs and files ###########
abbadate="20191230-highcov-singleInd-minInd-ancCA-AK-KUR" # with min ind filter
# http://www.popgen.dk/angsd/index.php/Abbababa
# use a different script than the multipop case:
JACKKNIFE=/u/home/a/ab08028/klohmueldata/annabel_data/bin/angsd/R/jackKnife.R ### this script can be used for error est (not doing here) but also for getting the D statistics
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
bamdir=$wd/bams/
ABBADIR=$wd/ABBA_BABA
#mkdir -p $ABBADIR
#mkdir -p $ABBADIR/$todaysdate

outdir=$ABBADIR/$abbadate


bamListDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/forABBABABA/singleInd
IDfile=angsd.bamList.mappedtoMfurfullpaths.HighCovPlusADNAOnly.InPopOrder.3inds.ancCA.AK.KUR.txt
#mfurBamList=$bamListDir/angsd.bamList.mappedtoMfurfullpaths.HighCovPlusADNAOnly.InPopOrder.txt

#errorList=$bamListDir/errorList.txt # this should be NA for all pops because not doing error correction (at least for now)
#mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
# note: sizeList.pops.txt must include a 1 at the end for the ancestral(maybe?)
### HAS TO BE MFUR TO HAVE ANCESTRAL!!! #########
### in R need to install.packages("pracma") in the anaconda session
############ run DSTAT : ############
Rscript $JACKKNIFE file=$outdir/angsdOut.TransvOnly.abbababa outfile=$outdir/DStats indNames=$bamListDir/$IDfile


# results 1)result.Observed.txt

#D-statistic calculated WITHOUT Error Correction and WITHOUT Ancient Transition removal

#2) result.ErrorCorr.txt

#D-statistic calculated WITH Error Correction and WITHOUT Ancient Transition removal

#3) result.ErrorCorr.TransRem.txt 
#D-statistic calculated WITH Error Correction and WITH Ancient Transition removal

#4) result.TransRem.txt <---- this is what you want

#D-statistic calculated WITHOUT Error Correction and WITH Ancient Transition removal 


# Specifically, the values contained in the four files are: 
# mean(D)=average D-stat, 
# JK-D=jackknife estimate of the D-stat, 
# V(JK-D)=variance of the D-stat, 
# Z=Z score, 
# pvalue=pvalue from the Z score, 
# nABBA=number of ABBA patterns observed, 
# nBABA=number of BABA patterns observed, 
# nBlocks=number of blocks with observed data, 
# H*=the names of the four populations for the specific tree. 
# Note that the number of patterns might not be integer because of how ANGSD treats multiple genomes per populations. 
source deactivate

