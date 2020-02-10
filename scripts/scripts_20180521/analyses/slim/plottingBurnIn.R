########## Plotting burn in ################
require(ggplot2)
require(reshape2)
rep=1
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/burnin/h_0/"
input <-  read.table(paste(data.dir,"replicate_",rep,".slim.output.BurnIn.allConcatted.summary.txt.gz",sep=""),header=T,sep=",")
###### take sum of chunks #####

input_melt <- melt(input,id.vars = c("generation","chunk","replicate","popSize"))
input_melt_sum <- input_melt %>%
  group_by(generation,replicate,popSize,variable) %>%
  summarise(sumOfAvgs=sum(value))
# don't sum up het - that is per bp so should be averaged 

ggplot(input_melt_sum[input_melt_sum$variable!="meanHet",],aes(x=generation,y=sumOfAvgs))+
  facet_wrap(~variable,scales="free_y")+
  geom_line()+
  theme_bw()
