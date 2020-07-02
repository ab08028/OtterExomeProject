### note: note that mapped to mfur had minMaf 0.01 filter applied to get rid of fixed sites because of long divergence
indir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/basicSNPStats/aDNA_TVSNPStats/20191212-highcov-AFprior-MajorMinor4-plusCOM-KUR-AL-RedoneToReplaceDeletedFiles/"
indIDs=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/angsd.bamList.mappedtoMfurfullpaths.HighCovPlusADNA.PlusCOM.KUR.AL.IDsOnly.txt")
# ancient=c("ind6TotDepth","ind7TotDepth","ind8TotDepth")
ancient=c("A13_Elut_CA_AN_388_SN1_2CAP_screen","A29_Elut_CA_SM_30_SN2_CAP","A30_Elut_CA_SM_35_SN1_CAP")
infodf <- data.frame()
refs=c("elut","mfur")
for(ref in refs){
input <- read.table(paste(indir,"/angsdOut.mappedTo",ref,".1e-06.snpsOnly.TransvOnly.counts.gz",sep=""),header=T,stringsAsFactors = F)

# want to count up non zero entries in each column
totalNonZeroEntriesPerInd <- data.frame(colSums=(colSums(input != 0)))
# bring in individual IDs (they are in same order so that is why cbind works -- be careful if your ids aren't in same order)
totalNonZeroEntriesPerInd$IDs <- indIDs$V1
# need to get individual mappings 
totalNonZeroEntriesPerInd[totalNonZeroEntriesPerInd$IDs %in% ancient,]
avgAncTVSnps = sum(totalNonZeroEntriesPerInd[totalNonZeroEntriesPerInd$IDs %in% ancient,]$colSums)/length(ancient)
avgAncTVSnps

avgModTVSnps = sum(totalNonZeroEntriesPerInd[!(totalNonZeroEntriesPerInd$IDs %in% ancient),]$colSums)/(length(totalNonZeroEntriesPerInd$IDs) - length(ancient))
avgModTVSnps

df = data.frame(ref=ref,avgModTVSNPs=avgModTVSnps,avgAncTVSNPs=avgAncTVSnps)
infodf <- rbind(infodf,df)
write.table(totalNonZeroEntriesPerInd,paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/basicSNPStats/aDNA_TVSNPStats/20191212-highcov-AFprior-MajorMinor4-plusCOM-KUR-AL-RedoneToReplaceDeletedFiles/PerIndividual.mappedTo.",ref,".AvgAncModernTVSNPs.txt",sep=""),sep="\t",row.names=F,quote=F)
}
write.table(infodf,"/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/basicSNPStats/aDNA_TVSNPStats/20191212-highcov-AFprior-MajorMinor4-plusCOM-KUR-AL-RedoneToReplaceDeletedFiles/Summary.AvgAncModernTVSNPs.txt",sep="\t",row.names=F,quote=F)
