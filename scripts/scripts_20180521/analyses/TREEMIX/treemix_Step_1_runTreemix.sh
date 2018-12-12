module load treemix
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=$gitdir/analyses/TREEMIX/
genotypeDate=20181119
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/
treeFileDir=$vcfdir/treemixFormat/
header=snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants
wd=/u/flashscratch/a/ab08028/captures/analyses/TREEMIX/$genotypeDate
# decisions: ld pruning or not? all individuals or not? 
# start with most basic treemix (no migration)
treemix -i $treeFileDir/${header}.frq.strat.treemixFormat.gz -o $wd/vanilla.treemix # change out stem to something more intelligent

# need to make a little wrapper to plot tree

# start with most basic treemix (no migration) + ld pruning (ish)
treemix -i $treeFileDir/${header}.frq.strat.treemixFormat.gz -k1000 -o $wd/vanilla.treemix.k1000 # change out stem to something more intelligent

# start with most basic treemix (no migration) + ld pruning (ish) + no sample size correction 
treemix -i $treeFileDir/${header}.frq.strat.treemixFormat.gz -k1000 -noss -o vanilla.treemix.k1000.noss # change out stem to something more intelligent


# start with most basic treemix + migration
treemix -i $treeFileDir/${header}.frq.strat.treemixFormat.gz -m -root KUR -o migration.treemix # change out stem to something more intelligent


# need some sort of LD pruning (something more sophisticated?)
treemix -i $treeFileDir/${header}.frq.strat.treemixFormat.gz -m -root KUR -o migration.treemix # change out stem to something more intelligent



# subset populations:
# columns:
cols=`zcat $treeFileDir/${header}.frq.strat.treemixFormat.gz | head -n1`
cols
# want to isolate AK, CA and BAJ (west coast pops):
zcat $treeFileDir/${header}.frq.strat.treemixFormat.gz | awk '{print $1,$3,$5}' > $treeFileDir/AK.CA.BAJ.only.{header}.frq.strat.treemixFormat.gz


# module load R
# R
# source("/u/flashscratch/a/ab08028/captures/analyses/TREEMIX/20181119/sandbox/plotting_func.R")
# plot_tree(migration.treemix)




# descriptoin of output from manual:
# TreeMix will output a number of files. If you have used the -o flag to designate the output stem outstem, these will be:
# 1. outstem.cov.gz. The covariance matrix (Wö in Pickrell and Pritchard [2012]) between pop- ulations estimated from the data
# 2. outstem.covse.gz. The standard errors for each entry in the covariance matrix
# 3. outstem.modelcov.gz. The fitted covariance (W in Pickrell and Pritchard [2012]) according
# to the model
# 4. outstem.treeout.gz. The fitted tree model and migration events
# 5. outstem.vertices.gz. This and the following file (outstem.edges.gz) contain the internal structure of the inferred graph. Modifying these files will cause issues if you try to read the graph back in, so we recommend against this.
# 6. outstem.edges.gz.
# The tree inferred from the data is in outstem.treeout.gz. The first line of this file is the Newick format ML tree, and the remaining lines contain the migration edges. The first column for these lines is the weight on the edge, followed (optionally) by the jackknife estimate of the weight, the jackknife estimate of the standard error, and the p-values. Then come the subtree below the origin of the migration edge, and the subtree below the destination of the migration edge.