### Set up file structure

SCRATCH=/u/flashscratch/a/ab08028/
wd=$SCRATCH/captures
fastqs=$wd/fastqs
bams=$wd/bams
reports=$wd/reports

# make dirs
mkdir -p $wd
mkdir -p $fastqs
mkdir -p $bams
mkdir -p $reports


# make sections for different reports:
mkdir -p $reports/step_1_fastqc
mkdir -p $reports/step_2_fq2sam
mkdir -p $reports/submissions 
### .. more 
