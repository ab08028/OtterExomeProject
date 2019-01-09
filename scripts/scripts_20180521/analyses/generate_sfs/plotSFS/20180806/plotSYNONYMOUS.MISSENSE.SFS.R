require(ggplot2)
require(dplyr)
require(grid)
require(purrr) # for map, reduce
require(readr) # for read_csv
require(tidyr)# for unnest 

todaysdate=format(Sys.Date(),format="%Y%m%d")
genotypeDate=20180806
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/SFS/",genotypeDate,"/cdsSFS/",sep="")
dir.create(plot.dir,recursive = T)
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/cdsSFS/",sep="")
neut.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/neutralSFS/",sep="")
################ Plot intermediate filtered SFSes #################
####################### fold sfs function ###############
foldSFS <- function(sfs){
  foldedSFS <- data.frame()
  ss=length(sfs$frequency) - 1 # this is the ss in chromosomes
  foldedBin=ss/2  # half as many as ss ; do seq from 0 --> foldedBins to include monomorphics
  # note you have to start seq at 0 to include monomorphic bin 
  for(i in seq(0,foldedBin)){
    # if it's not the center point (which doesn't get added together)
    # see wakeley coalescent eq 1.2
    if(i==ss-i){
      foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count)/2)) # if it's add midpoint (foldedLen, including 0 monomorphic bin), add together and divide by two (equivalent of not adding together)
    }
    else{
      foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count))) # if not at mid point, just add together like normal
    }
  }
  return(foldedSFS)
}



########################### get all files ###############
# list of syn/mis files:
fileList=list.files(pattern=paste("R.format",sep=""),path = data.dir,full.names = F )
#"One limitation of the previous approach is that we don’t keep any auxilliary information #we may want to, such as the filenames of the files read. To keep the filename alongside #the data, we can read the data into a nested dataframe rather than a list, using the #mutate() function from dplyr. This gives us the following result:"
# https://serialmentor.com/blog/2016/6/13/reading-and-combining-many-tidy-data-files-in-R

# list of neut files:
neutFileList=list.files(pattern=paste("R.format",sep=""),path = neut.dir,full.names = F )

# combine file lists:

data1 <- data_frame(filename = fileList) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_delim(delim="\t",file.path(data.dir, .)))# a new data column 
         
  )  
# neutral data: 
data2 <- data_frame(filename = neutFileList) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_delim(delim="\t",file.path(neut.dir, .)))# a new data column 
         
  )
# this tell you in red that is "Parsed with column specification" - this isn't an error
# combine data:
data <- rbind(data1,data2)
data <- unnest(data) # to make it useful to use, will put the filename as a column
# fold each SFS:
foldedData <- data %>%
  group_by(filename) %>%
  do(foldSFS(.))
#### omg this actually worked! 

############# get populations ###########
foldedData$population <- sapply(strsplit(as.character(foldedData$filename),"\\."),'[',1)
foldedData$label <- sapply(strsplit(as.character(foldedData$filename),"\\."),'[',2)
# update labels:
foldedData$newLabel <- foldedData$label
foldedData[foldedData$label=="all_9",]$newLabel <- "all_9_neutral"
# exclude monomorphic before getting proportions:
foldedData_exclMono <- foldedData[foldedData$frequency!=0,]
########## get proportional SFS ########
foldedData_exclMono_prop <- foldedData_exclMono %>%
  group_by(filename) %>%
  mutate(proportion=count/sum(count))
############# plot ###########
##### Plot counts: ###########
plot <- ggplot(foldedData_exclMono,aes(x=frequency,y=count,fill=newLabel))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~population,scales="free")+
  theme_bw()
plot

ggsave(paste(plot.dir,"allPops.",todaysdate,"syn.mis.neut.count.SFS.pdf",sep=""),plot,device="pdf",width=7,height=5)

pops=c("CA","AK","AL","COM","KUR")
allPlots_list=list()
for(i in (1:length(pops))){
  pop=pops[i]
  plot <- ggplot(foldedData_exclMono[foldedData_exclMono$population==pop,],aes(x=frequency,y=count,fill=newLabel))+
    geom_bar(stat="identity",position="dodge")+
    theme_bw()+
    ggtitle(pop)+
    theme(legend.position = c(0.5,0.6))
  ggsave(paste(plot.dir,pop,".",todaysdate,"syn.mis.count.SFS.pdf",sep=""),plot,device="pdf",width=10,height=5)
  allPlots_list[[i]]=plot
}
allPlots_list

##### Plot proportions: ##########
plot <- ggplot(foldedData_exclMono_prop,aes(x=frequency,y=proportion,fill=newLabel))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~population,scales="free")+
  theme_bw()
plot

ggsave(paste(plot.dir,"allPops.",todaysdate,"syn.mis.proportion.SFS.pdf",sep=""),plot,device="pdf",width=7,height=5)

pops=c("CA","AK","AL","COM","KUR")
allPlots_list=list()
for(i in (1:length(pops))){
  pop=pops[i]
  plot <- ggplot(foldedData_exclMono_prop[foldedData_exclMono_prop$population==pop,],aes(x=frequency,y=proportion,fill=newLabel))+
    geom_bar(stat="identity",position="dodge")+
    theme_bw()+
    ggtitle(pop)+
    theme(legend.position = c(0.5,0.6))
  ggsave(paste(plot.dir,pop,".",todaysdate,"syn.mis.proportion.SFS.pdf",sep=""),plot,device="pdf",width=10,height=5)
  allPlots_list[[i]]=plot
}
allPlots_list
