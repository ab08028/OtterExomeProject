####### Plot treemix output ##########
# source the treemix plotting funcs

source("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/analyses/TREEMIX/plotting_funcs.R") 
genotypeDate="20181119"
#pops="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/poporder.2.higherOrder.txt" # for plotting residuals
#pops="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/poporder.5.sepCA-BAJ.noFerret.txt" # for plotting residuals
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/",genotypeDate,"/snp7/properFerretExpt/",sep="")
######################## First do models that include BAJA ######################
#models=c("ferretOutgroup.noMig.treemix") # can specify by hand or get from list files in the dir
markers=c("sepCA-BAJ.exclRelatives-testADDINGMFUR")
for(marker in markers){
  models=list.files(data.dir,pattern = marker)
  # pop file for models incld BAJA:
  pops=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/poporder.",marker,".txt",sep="") # for plotting residuals
  
  # loop through models and save plots
  for(modelName in models){
    # go to wd of that model
    setwd(paste(data.dir,modelName,sep=""))
    # tree plot:
    png(paste("plot.",modelName,".treePlot.png",sep=""),width=2400,height=2000,res=300)
    plot_tree(modelName)
    dev.off()
    # residuals: 
    png(paste("plot.",modelName,".residualsPlot.png",sep=""),width=2400,height=2000,res=300)
    plot_resid(modelName,pop_order = pops) #poporder is a text file with populations in the order you want them
    dev.off()
    
  }
}

