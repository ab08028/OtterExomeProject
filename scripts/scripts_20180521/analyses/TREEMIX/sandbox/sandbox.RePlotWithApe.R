### make a nicer tree in ape
require(ape)
# read in the tree -- well that was easyish
test <- read.tree("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/snp7/goodParams-USETHESE/root.CA.mig.4.k.500.global.noBAJA.exclRelatives.treemix/root.CA.mig.4.k.500.global.noBAJA.exclRelatives.treemix.treeout.gz")
plot.phylo(test,use.edge.length = T,no.margin = T,type = "cladogram")
axisPhylo(side=1)
