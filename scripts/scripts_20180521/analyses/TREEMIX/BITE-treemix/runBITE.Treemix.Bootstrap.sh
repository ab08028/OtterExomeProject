######## submit it:
source /u/local/Modules/default/init/modules.sh
module load treemix
indir=/u/flashscratch/a/ab08028/captures/vcf_filtering/20181119_filtered/treemixFormat
infile=snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.noBAJA.exclRelatives.frq.strat.treemixFormat.gz		## treemix input file
# snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.noBAJA.exclRelatives.frq.strat.treemixFormat.gz
# snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.noBAJA.frq.strat.treemixFormat.gz


ncore=$3 		## max number of cores to use

blockk=100 		## block size

outgroup=NoOutgroup 	## name of the selected outgroup population (if you want to do an unrooted ML tree put here 'NoOutgroup' (without quotes))

nboot=100 # start with 100		## number of bootstrap replicates of the TREE

pathP=/u/home/a/ab08028/klohmueldata/annabel_data/bin/phylip-3.697/exe/consense		## path to Phylip consense program. Example: /biosoftware/phylip/phylip-3.696/exe/consense

outdir=/u/flashscratch/a/ab08028/captures/analyses/TREEMIX/20181119/BITE
outname=mig.${numk}.k.${blockk}.out.${outgroup}.out		## name for output file

script=Treemix_bootstrap.sh

for mig in 5
do
# can loop over migrations
$script $indir/$infile $mig $ncore $blockk $outgroup $nboot $pathP $outdir/$outname
do