##### Read angsd sfses in R #########

angsdDate=20190506
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/SFS/",angsdDate,sep="")
SFSes <- list.files(data.dir,pattern="SFS.txt",full.names = T)

##################### functions #######################
calculatePiFromSFS_empData <- function(unfoldedHaploidSampleSize,totalSimulatedSites,sfsFoldedExcludingMonomorphic){
  n=as.numeric(unfoldedHaploidSampleSize)
  #print(n)
  totalSites=totalSimulatedSites
  sumTotal=0
  # need to skip over monomorphic so add 1 to each thing in the sequence
  # make sure monomorphic is excluded: (should already be, but this makes sure)
  sfsExclMono = sfsFoldedExcludingMonomorphic[sfsFoldedExcludingMonomorphic$frequency!=0,]
  #print(sfsExclMono)
  for(i in seq(1,n/2)){
    # greek letter eta (pg 16 of Wakeley) is the folded count
    eta = sfsExclMono[i,2] # second column is count 
    #print(eta)
    # from the equation 
    sumTotal=sumTotal+(i*(n-i)*eta)
    #print(sumTotal)
  }
  # then multiply sumTotal by 1/(n choose 2)
  pi_overall = (1/(choose(n,2))) * sumTotal
  #print(pi_overall)
  # then divide by total sites to get pi per site (note that monomorphic sites were also projected)
  pi_perSite = pi_overall/totalSites
  return(pi_perSite)
}

####################### functions to calculate S and Wattersons Theta ############
# watterson's theta uses n-1 haploid individuals for harmonic number
# my harmonicNumber_2nMinus1 function takes in the number of *diploid* individuals
# and converts it to haploid 
harmonicNumber_2nMinus1 <- function(num_indv){
  n_dip=num_indv # the numver of diploid individuals
  n_hap=2*n_dip # number of haploids
  n_harm=n_hap - 1 # number of haps minus 1 for harmonic number
  harm_number=0
  for(i in 1:(n_hap-1)){
    harm_number <- harm_number + 1/i
  }
  return(harm_number)
}

# example for 26 individuals:
#harmonicNumber_2nMinus1(26) # is the harmonic # for 2*26-1 = 51 -->  4.518813 # checked it here for 51: https://www.math.utah.edu/~carlson/teaching/calculus/harmonic.html

#### WATTERSON'S THETA (uses harmonic number function)
# this function will calculate watterson's theta using the number of diploid inds, 
# number of SNPs, and the number of callable (HQ) sites used for the calculation of # of SNPs. 
wattersons_theta <- function(num_indv,numSNPs,callablesites){
  harm_num <- harmonicNumber_2nMinus1(num_indv)
  theta <- numSNPs / harm_num
  theta_divbyL <- theta/callablesites
  return(theta_divbyL)
}

#################### read in SFSes, get Pi #############

allSFSes <- data.frame()
allPiTheta <- data.frame()
######### THESE SFSes ARE FOLDED. DO NOT USE UNFOLDED SFSES IN THIS SCRIPT ########### 
for(SFS in SFSes){
  input = read.table(SFS,header=F)
  colnames(input) <- seq(0,length(input)-1) # label frequencies; first bin is 0
  input = data.frame(t(input)) # transpose
  colnames(input) <- "count"
  input$frequency <- row.names(input) # set frequencies as a column
  # get metadata from file name:
  info <- tail(unlist(strsplit(SFS,"/")),n=1) # pull out file name, not path
  input$filename <- info
  input$label <- unlist(strsplit(info,"\\."))[2]
  input$reference <- unlist(strsplit(info,"\\."))[3]
  input$proportion <- input$count / sum(input[input$frequency!=0,]$count)
  input[input$frequency==0,]$proportion <- NA # make proportion for 0 bin 
  # a short hand label:
  input$label2 <- unlist(strsplit(unlist(strsplit(info,"angsd."))[2],".saf"))[1]
  allSFSes = rbind(allSFSes,input)
  # input: unfolded diploid sample size, 
  unfoldedHaploidSampleSize=as.numeric(max(input$frequency))*2 # highest frequency * 2 is haploid sample size (number of unfolded bins)
  # when you do total sites here it should INCLUDE the monomorphic bin (input$count incl 0 bin)
  # then when you input SFS it should NOT include the 0 bin:  input[input$frequency!=0,c("frequency","count")]
  pi=calculatePiFromSFS_empData(unfoldedHaploidSampleSize,sum(input$count),input[input$frequency!=0,c("frequency","count")])
  piInfo=data.frame(pi=pi)
  piInfo$filename <- info
  piInfo$label <- unlist(strsplit(info,"\\."))[2]
  piInfo$reference <- unlist(strsplit(info,"\\."))[3]
  piInfo$label2 <- unlist(strsplit(unlist(strsplit(info,"angsd."))[2],".saf"))[1]
  # calculate watterson's theta:
  piInfo$S = sum(input[input$frequency!=0,]$count)
  # watterson's theta function takes  function(num_indv,numSNPs,callablesites)
  # this takes the DIPLOID sample size and converts to haploid!!!!!!! don't make a mistake here! 
  # so DON'T multiply the diploid size by two in your function input (different from pi above)
  piInfo$WattersonsTheta=wattersons_theta(as.numeric(max(input$frequency)),sum(input[input$frequency!=0,]$count),sum(input$count))
  allPiTheta=rbind(allPiTheta,piInfo)
}
  
# plot count SFSes
ggplot(allSFSes[allSFSes$frequency!=0,],aes(x=as.numeric(frequency),y=count,fill=reference))+
  geom_bar(stat="identity",position="dodge")+
  theme_bw()+
  facet_wrap(label~TransRemoved,scales="free")


# plot proportional sfses:
ggplot(allSFSes[allSFSes$frequency!=0,],aes(x=as.numeric(frequency),y=proportion,fill=reference))+
  geom_bar(stat="identity",position="dodge")+
  theme_bw()+
  facet_wrap(~label2)

##### just plot a few SFSes with minInd3 and mapped to Mfur ####
allSFSes_select <- allSFSes[allSFSes$reference=="mappedToMfur" & grepl("minInd",allSFSes$label2) & allSFSes$frequency!=0,] # no monomorphic, only mapped to ferret, minInd3
allSFSes_select$transv <- "Transitions+Transversions"
allSFSes_select[grep("TransversionsOnly",allSFSes_select$label2),]$transv <- "TransversionsOnly"
p2a <- ggplot(allSFSes_select,aes(x=frequency,y=count,fill=transv))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~label,scales="free")+
  theme_bw()+
  ggtitle("SFSes -- Counts \nminInd=3")+
  theme(legend.title = element_blank())
p2a
ggsave(paste(data.dir,"/SFSes.allGroups.TiTv.minInd3.Counts.pdf",sep=""),p2a,height=5,width=9)


p2b <- ggplot(allSFSes_select,aes(x=frequency,y=proportion,fill=transv))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~label,scales="free")+
  theme_bw()+
  ggtitle("SFSes -- Proportions\nminInd=3")+
  theme(legend.title = element_blank())
p2b
ggsave(paste(data.dir,"/SFSes.allGroups.TiTv.minInd3.Proportions.pdf",sep=""),p2b,height=5,width=9)

############## plot pi: #############
# just doing mapped to MFUR -- need to troubleshoot elut
p3a <- ggplot(allPiTheta[allPiTheta$reference=="mappedToMfur",],aes(y=pi,x=label2,fill=label))+
  geom_bar(stat="identity",position="dodge",alpha=0.6)+
  coord_flip()+
  theme_bw()+
  xlab("Filters and data info")
p3a
ggsave(paste(data.dir,"/allPis.mappedToMfur.pdf",sep=""),p3a,height=5,width=9)
# just doing mapped to Elut -- need to troubleshoot elut
p3b <- ggplot(allPiTheta[allPiTheta$reference=="mappedToElut",],aes(y=pi,x=label2,fill=label))+
  geom_bar(stat="identity",position="dodge",alpha=0.6)+
  coord_flip()+
  theme_bw()+
  xlab("Filters and data info")
p3b
ggsave(paste(data.dir,"/allPis.mappedToElut-WEIRD.pdf",sep=""),p3b,height=5,width=9)


####### Pull out cool ones for wayne lab meeting : minInd3 with and without transversions and ONLY mapped to Mfur: 
allPiTheta_MfurOnly <- allPiTheta[allPiTheta$reference=="mappedToMfur",]
allPiTheta_MfurOnly$transv <- "Transitions+Transversions"
allPiTheta_MfurOnly[grep("TransversionsOnly",allPiTheta_MfurOnly$filename),]$transv <- "TransversionsOnly"

p4a <- ggplot(allPiTheta_MfurOnly[grep("minInd",allPiTheta_MfurOnly$filename),],aes(y=pi,x=label2,fill=label))+
  geom_bar(stat="identity",position="dodge",alpha=0.6)+
  coord_flip()+
  theme_bw()+
  xlab("Filters and data info")+
  theme(legend.title = element_blank(),legend.position=c(0.5,0.5),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))
p4a
ggsave(paste(data.dir,"/minInd3.Pis.mappedToMfur.pdf",sep=""),p4a,height=8,width=14)

# make another plot that is just transversions:

p4b <- ggplot(allPiTheta_MfurOnly[allPiTheta_MfurOnly$transv=="TransversionsOnly" & grepl("minInd",allPiTheta_MfurOnly$filename),],aes(y=pi,x=label2,fill=label))+
  geom_bar(stat="identity",position="dodge",alpha=0.6)+
  coord_flip()+
  theme_bw()+
  xlab("Filters and data info")+
  theme(legend.title = element_blank(),legend.position=c(0.5,0.5),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  ggtitle("Transversions Only")
p4b
ggsave(paste(data.dir,"/minInd3.Pis.mappedToMfur.TransversionsOnly.pdf",sep=""),p4b,height=8,width=14)

p4c <- ggplot(allPiTheta_MfurOnly[allPiTheta_MfurOnly$transv=="Transitions+Transversions" & grepl("minInd",allPiTheta_MfurOnly$filename),],aes(y=pi,x=label2,fill=label))+
  geom_bar(stat="identity",position="dodge",alpha=0.6)+
  coord_flip()+
  theme_bw()+
  xlab("Filters and data info")+
  theme(legend.title = element_blank(),legend.position=c(0.5,0.5),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  ggtitle("Transitions+Transversions")
p4c
ggsave(paste(data.dir,"/minInd3.Pis.mappedToMfur.TransitionsAndTransversions.pdf",sep=""),p4c,height=8,width=14)

############# Plot watterson's theta with Pi ##############
require(reshape2)
allPiTheta_MfurOnly_melt <- melt(allPiTheta_MfurOnly)
# exclude S as a variable

p5a <- ggplot(allPiTheta_MfurOnly_melt[allPiTheta_MfurOnly_melt$variable!="S" & grepl("minInd",allPiTheta_MfurOnly_melt$filename),],aes(y=value,x=label2,fill=variable))+
  geom_bar(stat="identity",position="dodge",alpha=0.6)+
  coord_flip()+
  theme_bw()+
  xlab("Filters and data info")+
  theme(legend.title = element_blank(),legend.position=c(0.5,0.5),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))
p5a
ggsave(paste(data.dir,"/minInd3.PisPlusWattersonsTheta.mappedToMfur.pdf",sep=""),p5a,height=8,width=14)


p5b <- ggplot(allPiTheta_MfurOnly_melt[allPiTheta_MfurOnly_melt$variable!="S" & grepl("minInd",allPiTheta_MfurOnly_melt$filename) & allPiTheta_MfurOnly_melt$transv=="TransversionsOnly",],aes(y=value,x=label2,fill=variable))+
  geom_bar(stat="identity",position="dodge",alpha=0.6)+
  coord_flip()+
  theme_bw()+
  xlab("Filters and data info")+
  theme(legend.title = element_blank(),legend.position=c(0.5,0.5),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  ggtitle("TransversionsOnly")
p5b
ggsave(paste(data.dir,"/minInd3.PisPlusWattersonsTheta.mappedToMfur.TransversionsOnly.pdf",sep=""),p5b,height=8,width=14)


p5c <- ggplot(allPiTheta_MfurOnly_melt[allPiTheta_MfurOnly_melt$variable!="S" & grepl("minInd",allPiTheta_MfurOnly_melt$filename) & allPiTheta_MfurOnly_melt$transv=="Transitions+Transversions",],aes(y=value,x=label2,fill=variable))+
  geom_bar(stat="identity",position="dodge",alpha=0.6)+
  coord_flip()+
  theme_bw()+
  xlab("Filters and data info")+
  theme(legend.title = element_blank(),legend.position=c(0.5,0.5),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  ggtitle("Transitions+Transversions")
p5c
ggsave(paste(data.dir,"/minInd3.PisPlusWattersonsTheta.mappedToMfur.TransitionsAndTransversions.pdf",sep=""),p5c,height=8,width=14)
