######## Get average called sites across all individuals: 
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/sumGPsGLsPerVEPCategory/"
dates=c("20190701-highcov-AFprior-MajorMinor4","20190701-lowcov-AFprior-MajorMinor4")
minGP=0.95
minDepths=c(1,2,4)
minInd=1
# get counts of all CDS sites that pass filters:
# start a file:
for(minDepth in minDepths){
sink(paste(data.dir,"AVERAGECALLEDSITES.allInds.HighCov.LowCov.minDepth.",minDepth,".minInd.",minInd,".minGP.",minGP,".txt",sep=""))
cat("date\tminDepth\tavgCalledCDSSites\n")
for(date in dates){
  input <- read.table(paste(data.dir,date,"/angsdOut.mappedTomfur.hetHomTotals.GPs.ProbCutoff.",minGP,".DepthCutoff.",minDepth,".minInd.",minInd,".",date,".CDS.txt",sep=""),header = T)
  avgCalledSites <- mean(input$callableSites)
  #cat("date: ",date," minDepth: ",minDepth,"avgCalledCDSSites: ",avgCalledSites)
  cat(date,"\t",minDepth,"\t",avgCalledSites,"\n")
  #cat("\n")
  
}

sink()
}
