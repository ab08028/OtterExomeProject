genotypeDate=20181119
results.dir = paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/",genotypeDate,"/",sep="")

# want to sort by LL
# want to rescale results by theta 
# and write out as a table that can be easily examined in Excel
mu= # can test different mus 
test$Nanc <- test$theta / (4*test$mu*test$L)

# scale:
test$nuB_scale <- test$nuB * test$Nanc
test$nuF_scale <- test$nuF * test$Nanc
test$TB_scale <- test$TB * test$Nanc * 2
test$TF_scale <- test$TF * test$Nanc * 2


test[test$LL==max(test$LL),]


