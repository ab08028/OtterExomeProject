#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=25G,arch=intel*
#$ -m bea


## username, error locations, etc. are set by submission scripts.
adapterRemoval=/u/home/a/ab08028/klohmueldata/annabel_data/bin/AdapterRemoval

$adapterRemoval --file1 $1 --file2 $2 --basename $3 \
--trimns \
--trimqualities \
--collapse \
--identify-adapters  \
--combined-output \
--gzip
# trims qualities below 2; collapses reads (good for aDNA); trims NNs at ends of reads; 
# tries to identify adapters from PE reads (in case Illumina TruSeq not appropriate)
# combines the output
# gzips the output

#test:
$adapterRemoval --file1 A1_Elut_CA_AN_396_SN1_S61_R1_001.fastq.gz --file2 A1_Elut_CA_AN_396_SN1_S61_R1_001.fastq.gz --basename A1_Elut_CA_AN_396_SN1 \
--trimns \
--trimqualities \
--collapse \
--identify-adapters  \
--combined-output \
--gzip