test <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20180806/Bottleneck_1D/all.output.concatted.txt",header=T)

test$Nanc <- test$theta / (4*test$mu*test$L)

# scale:
test$nuB_scale <- test$nuB * test$Nanc
test$nuF_scale <- test$nuF * test$Nanc
test$TB_scale <- test$TB * test$Nanc * 2
test$TF_scale <- test$TF * test$Nanc * 2


test[test$LL==max(test$LL),]


