## these files are the result of easySFS_2_parseProjection
genotypeDate="20180806"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/",sep="")
plot.dir=data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/SFS/",genotypeDate,"/easySFS_projection/",sep="")
#dir.create(plot.dir)
# get list of files:
fileList=list.files(pattern=paste("R.format",sep=""),path = data.dir,full.names = F )

data1 <- data_frame(filename = fileList) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_delim(delim=",",file.path(data.dir, .)))# a new data column 
         
  )  


data <- unnest(data1)
data$population <-  sapply(strsplit(as.character(data$filename),"\\."),'[',1)

# summarise maxima: 
data %>% 
  group_by(population) %>%
  summarise(max(snps))

p0 <- ggplot(data,aes(x=projection,y=snps))+
  geom_point()+
  facet_wrap(~population,scales="free_x")+
  theme_bw()+
  scale_x_continuous(breaks=seq(2,80,2))+
  geom_vline(xintercept = 20,color="red")
p0
# want to find the maximum # of snps that still has sample size >10
ggsave(paste(plot.dir,"/easySFS.projections.allPops.pdf",sep=""),p0,height=5,width=10)
  
##### Want to plot difference in snp number #######
before <- data %>% 
  group_by(population) %>%
  filter(projection == max(projection)) 
before$label <- "original"

after <- data %>% 
  group_by(population) %>%
  filter(projection == 20) 
after$label <- "projection"
both <- rbind(before,after)
p1 <- ggplot(both,aes(x=population,y=snps,fill=label))+
  geom_bar(stat="identity",position="dodge")+
  theme_bw()+
  ggtitle("Increase in SNPs (of all types) per population predicted after projection")
p1
ggsave(paste(plot.dir,"/easySFS.predictedGainInSNPs.postProj.allPops.pdf",sep=""),p1,height=5,width=10)
