###### Need to convert vcf file to plink bed format 
# isn't there some issue with chromosomes?
# need to make sure it's only SNPs not invariant sites
# so want to convert my final SNP file

source /u/local/Modules/default/init/modules.sh
module load plink

indir=
infile=
outdir=$SCRATCH/captures/analysis/FASTSTRUCTURE
# you need to use const-fid 0 otherwise it thinks that family name_sample name is structure of ID and tries to split it (and fails)
# allow extra chromosomes: to get it to get over the fact that chr names are non standard (make sure these wont get ignored?)
plink --vcf $indir/$infile --make-bed --keep-allele-order --const-fid 0 --allow-extra-chr --maf 0.05 -out $indir/XX.test
### note for faststructure to work you have to filter on maf 0.05
#testing: infile=snp_5_passingAllFilters_postMerge_elut.raw_variants.20170914.vcf.gz
# see if faststructure works