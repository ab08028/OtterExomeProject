####### Plot treemix output ##########
# source the treemix plotting funcs

source("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/analyses/TREEMIX/plotting_funcs.R") 
genotypeDate="20181119"
pops="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/poporder.txt" # for plotting residuals
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/",genotypeDate,"/snp9a/",sep="")
#models=c("ferretOutgroup.noMig.treemix") # can specify by hand or get from list files in the dir
models=list.files(data.dir,pattern = "treemix")
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

# residuals explanation : "Positive residuals indicate pairs of populations where the model underestimates the observed covariance, and thus populations where the fit might be improved by adding additional edges. Negative residuals indicate pairs of populations where the model overestimates the observed covariance; these are a necessary outcome of having positive residuals, but can also sometimes be interpreted as populations that are forced too close together due to unmodeled migration elsewhere in the graph. "

