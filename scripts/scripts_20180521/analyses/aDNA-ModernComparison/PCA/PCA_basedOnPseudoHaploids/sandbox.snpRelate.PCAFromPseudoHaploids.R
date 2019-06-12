##### pseudocode ######

# PED2GDS : need ped, map and can output a gds file # do this as its own step on hoffman

# ld prune (see other pca scripts)

# PCA with all individuals (with strict maf and missingness filters) --> plot 

# PCA with projections
## first do PCA based on modern individuals only (snp PCA with filters)
## snpgdspcasnploading -- need pca from previous step, gds file -- gets loading of snps
## snpgdspcasamploading -- need output form previous step, gds file, list of ancient samples; this will project the ancient samples onto the PCs --> plot

### results: for high and low cov datasets (transversions only, Ti+TV) should have PCAs either based on all individuals and based on projections #####
