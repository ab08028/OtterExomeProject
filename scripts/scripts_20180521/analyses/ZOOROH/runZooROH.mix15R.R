# installation; 
### note: these files contain monomorphic 1/1 sites, but those are automaticall
# excluded in the HMM method of ZooRoh. but you should remove them in other methods (Plink etc.)
# note that hbd refers to "homozygous by descent"
#install.packages("RZooRoH")
# install.packages("R.utils") # for reading gzip
require(RZooRoH)
require(R.utils)
require(ggplot2)
require(ggplot)
require(RColorBrewer)
colors=c(brewer.pal(n=8,"Spectral"),brewer.pal(n=8,"Paired"))
modelName="mix15R"
wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/zooROH/"
### converted from vcf > ped > oxford Gen , which is the "GP" format in ZooROH
pops=c("CA","AK","AL","COM","KUR") # alreadydid AK 
for(pop in pops){
#pop="AK"
print(paste("Starting pop",pop))
data.dir=paste(wd,pop,"/",sep="")
out.dir=paste(wd,pop,"/",modelName,"/",sep="")
dir.create(out.dir,showWarnings = F)
##################### read in the population specific data and sample IDs #############
zooInput <- zoodata(paste(data.dir,"/",pop,".Oxford.gen.gz",sep=""),zformat="gp",samplefile = paste(data.dir,"/",pop,".Oxford.SampleIDs",sep=""),min_maf =0.05) # R.utils lets you read in a gzipped file. 
#[1] "Number of positions in original file ::" "1571275"                                
#[1] "Number of positions after MAF filtering  ::" "1571275" 
# apparently ZooROH is not sensitive to maf filtering, so you don't need to do it
# can access various parts of zooInput using @ symbol
# eg zooInput@sample_ids
# try with maf 0.05 just to make it smaller (monos will be ignored no matter what)
# 
##################### set up the model #############
# want to use the fixed Rx model (mix model0 so that each individual has same set of rates for comparison
# the default is good: 10 rates each a factor of 2 away (2,4,8,16,32, ) -- 9 for hbd, and one for non-hbd
model <- zoomodel(K=15,base=2) # mix15R
model

#tempModel <- zoomodel(K=5,base=2) # temporary small model
#tempModel
# default values:
# type of model ; mixkr
# error rates 0.001
# krates:
# 2   4   8  16  32  64 128 256 512 512 # note last one repeats and refers to non-hbd segments (I think; not that clear in manual)
# mix coeff :
#0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.91
# rate info from manual
# The rate of the class is approximately equal to the size of the inbreeding loop in generations (Druet and Gautier, 2017). So, the rate of the class is approximately equal to twice the number of generations to the common ancestor. For instance, a class with a rate Rk equal to 10 would correspond approximately to ancestors five generations ago whereas a class with a rate Rk equal to 100 would correspond approximately to ancestors fifty generations ago. These are off course no precise estimations of age of HBD segments but rather qualitative measures.
############################### run the model ##########################
results <- zoorun(model, zooInput) # don't need to specify ids because this has already been separated by individual. Am assigning to pop_results eg AL_results
# final  value 14502.636500 
# converged -- want to log this somewhere
## save results:
#Save and restore multiple R objects: 
save(results, file = paste(out.dir,pop,".ZooROH.Results.RData",sep=""))
# can then load with: load(“my_data.RData”)
############################## can start here if you don't want to reload everything; just comment out previous loop ############
#load(paste(data.dir,pop,".ZooROH.Results.RData",sep="")) ## this works!!! So don't have to run everything above this once you've done it once

################## examine the results ##############
# what i want to write out: Krates,ll, bic, ...
# rows are per-individual
# mixing coefficients: results@mixc
# krates: results@krates
# get LL or BIC:
#results@modlik
#results@modbic
# results@realized #realized autozygosity per HBD class, or the partitioning of the genome is different classes is one of the main outputs of the model.

####### functions to summarize some things #########
###### contributions of diff classes (realized autozygosity) ###### 
#realAutozygosity <- realized(results) # The function returns a data frame with one row per individual and one column per extracted classes.shows partitioning of genome into different classes (more sophisticated than just summing up different bins and dividing nby genome size apparently)

##### inbreeding coefficient: needs a threshold T after which you say 'not inbred' #####
#thresholds=c(10,20,50,100)
#inbreedingCoefficients <- data.frame(pop=pop)
#for(threshold in thresholds){
#  F <- data.frame(cumhbd(results,T=threshold))
#  colnames(F) <- paste("F",threshold,sep="")
#  inbreedingCoefficients <- cbind(inbreedingCoefficients,F)
#}

######## skipping for now; exracting ROH for a particular genomic region ####
# rohbd  pg 31 of manual
# probhbd to get individual prob of HBD at a particular region pg 32 of manual
##################### Plotting ##############################

########## zooplot_prophbd: prob of genome assoc with diff. classes POP-WIDE #######
pdf(paste(out.dir,pop,".zooROH.outputPlots.pdf",sep=""),height = 8,width=11)

zooplot_prophbd(list(pop = results), cols = 'tomato', style = 'boxplot')


# multiple populations:
#zooplot_prophbd(list(Soay=soay_mix10r,Wiltshire=wilt_mix10r, RasaAragonesa=rara_mix10r),style='barplot')

#### inbreeding coeffs: use cumulative dist from prophbd
#To plot the average inbreeding coefficients (cumulative values) for three populations with lines:
zooplot_prophbd(list(pop=results),style='lines', cumulative = TRUE)

#################### plot individual-wide inbreeding coeffs; can encompass multiple pops too. #########
zooplot_individuals(results,cumulative=TRUE)

################# partition individual genomes into different classes #######
#par(cols=colors)
#zooplot_partitioning(list(pop=results), ylim = c(0,0.5), nonhbd = FALSE,col = c("#D53E4F", "#F46D43", "#FDAE61" ,"#FEE08B", "#E6F598", "#ABDDA4" ,"#66C2A5" ,"#3288BD", "#A6CEE3", "#1F78B4", "#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F"))
### see manual here -- a lot of options. 

##### plot distribution of lengths
# ggplot(CA_hbd,aes(x=length/1e06,fill=as.factor(HBDclass),color=as.factor(HBDclass)))+
#   geom_density(alpha=0.4)+
#   theme_bw()+
#   scale_x_continuous(label=comma)+
#   xlab("Length (Mb)")+
#   ggtitle("Length distribution per class")
dev.off()
}


