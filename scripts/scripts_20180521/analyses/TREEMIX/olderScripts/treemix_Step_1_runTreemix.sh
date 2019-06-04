####### NOTES: the weird RWAB samples seemed to cause false migration edges to ferret; make sure they are not in your dataset
# on 20190307 I remade my plink files without the RWAB files (need to redo faststructure at some point too)
# going to redo treemix and see if those weird mig edges are gone (I hope so!)

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
##### aha! you need to specify "m" -- the number of migration events! and can ratchet it up and down 
infile=${header}.frq.strat.InclFerret.treemixFormat.gz
########### model 1: no ld pruning, no migration edges, ferret as root 

outdir=ferretOutgroup.noMig.treemix 
mkdir $wd/$outdir
echo "vanilla model: no migration, no ld pruning" > $wd/$outdir/modelInfo.txt
echo "model command: treemix -i $treeFileDir/$infile -root FERRET -o $wd/$outdir/ferretOutgroup.noMig.treemix" >> $wd/$outdir/modelInfo.txt
treemix -i $treeFileDir/$infile -root FERRET -o $wd/$outdir/ferretOutgroup.noMig.treemix 

########## model 2: adding in ld pruning 
# start with most basic treemix (no migration) + ld pruning (ish)
# i don't know how much the snp grouping makes sense, because my snps are across many scaffolds and are low density? 
# see if it makes a difference
outdir=ferretOutgroup.ldPruning.1000snps.treemix
mkdir $wd/$outdir
echo "model: no migration, ld pruning groups of 1000 snps" > $wd/$outdir/modelInfo.txt
echo "model command: treemix -i $treeFileDir/$infile -root FERRET -k1000 -o $wd/$outdir/ferretOutgroup.ldPruning.1000snps.treemix" >> $wd/$outdir/modelInfo.txt
treemix -i $treeFileDir/$infile -root FERRET -k1000 -o $wd/$outdir/ferretOutgroup.ldPruning.1000snps.treemix  # change out stem to something more intelligent

# skip: start with most basic treemix (no migration) + ld pruning (ish) + no sample size correction 
#treemix -i $treeFileDir/${header}.frq.strat.treemixFormat.gz -k1000 -noss -o vanilla.treemix.k1000.noss # change out stem to something more intelligent


# start with most basic treemix + migration
outdir=ferretOutgroup.migration.treemix
mkdir $wd/$outdir
echo "model: migration, no ld pruning" > $wd/$outdir/modelInfo.txt
echo "model command: treemix -i $treeFileDir/$infile -m -root FERRET -o $wd/$outdir/ferretOutgroup.migration.treemix" >> $wd/$outdir/modelInfo.txt
treemix -i $treeFileDir/$infile -m 10 -root FERRET -o $wd/$outdir/ferretOutgroup.migration.treemix  # change out stem to something more intelligent



# start with most basic treemix + migration + k1000 ld pruning
outdir=ferretOutgroup.migration.ldPruning.k1000.treemix
mkdir $wd/$outdir
echo "model: migration, with k=1000 ld pruning" > $wd/$outdir/modelInfo.txt
echo "model command: treemix -i $treeFileDir/$infile -m -root FERRET -k1000 -o $wd/$outdir/ferretOutgroup.migration.treemix" >> $wd/$outdir/modelInfo.txt
treemix -i $treeFileDir/$infile -m 10 -root FERRET -k1000 -o $wd/$outdir/ferretOutgroup.migration.ldPruning.k1000.treemix  # change out stem to something more intelligent



# subset populations:
# columns:
#cols=`zcat $treeFileDir/${header}.frq.strat.treemixFormat.gz | head -n1`
#cols
# want to isolate AK, CA and BAJ (west coast pops):
#zcat $treeFileDir/${header}.frq.strat.treemixFormat.gz | awk '{print $1,$3,$5}' > $treeFileDir/AK.CA.BAJ.only.{header}.frq.strat.treemixFormat.gz


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