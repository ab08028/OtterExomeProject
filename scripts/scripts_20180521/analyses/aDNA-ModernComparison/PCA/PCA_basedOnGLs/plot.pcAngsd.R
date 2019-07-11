library(RcppCNPy)
#comp=c(1,2) # set which two PCs you want to plot ; for now let's go with 1 and 2


###### Doing two sets: all individuals or just low coverage individuals ##########
HCSampleList=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt")
LCSampleList=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.LowCoverageOnly.BeCarefulOfOrder.txt")
################################# high coverage ###############################
angsddates=c("20190701-highcov-AFprior-MajorMinor4") # high cov only
refs=c("elut","mfur") # set what reference you mapped to 
states=c("1e-06.snpsOnly", "1e-06.snpsOnly.TransvOnly")
minMafs=c("0.12","0.2","0.05")
# see below for downsampled inds only 
# with and without transversions
for(angsddate in angsddates){
  data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/PCA/",angsddate,"/",sep="")
  
  for(ref in refs){
    # ooh individuals may be in different orders! 
    
    #bamList <- read.table(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/angsd.bamList.mappedto",ref,"fullpaths.txt",sep=""),stringsAsFactors = F)
    for(state in states){
      for(minMaf in minMafs){
      ####### Plot PCA from covariance matrix #########
      # following example in https://github.com/mfumagalli/ngsPopGen/blob/master/scripts/plotPCA.R
      ### test: read in binary covariance matrix:
      # read in covariance matrix
      covar <- npyLoad(paste(data.dir,"pcAngsd.",ref,".",state,".minMaf.",minMaf,".cov.npy",sep="")) # ready from binary ## NOTE; this test dataset is only based on a tiny part of a chromosome. Only 4 snps! I think this test was mapped to elut. let's see
      # read numpy binary arrays  in R 
      # http://dirk.eddelbuettel.com/code/rcpp.cnpy.html
      
      #dim(covar) # 15 x 15. so it's all 15 individuals with a covariance matrix between them 
      
      ####### assign information from list of bamfiles ########
      # THIS IS A USEFUL BIT OF CODE THAT PULLS SAMPLE ID, REFERENCE AND WHETHER WAS DOWNSAMPLED
      # from list of bam names that was used in angsd
      # the order of these individuals will match the rows of the PCs 
      bamList=HCSampleList # HIGH COVERAGE
      colnames(bamList) <- "sampleID"
      #bamList$bamName <- lapply(strsplit(bamList$bamPath,"/"),tail,n=1) # tail pulls out last entry. useful because the different paths have different numbers of dirs so it's not always the same index # but it is the last entry always
      #bamList$sampleID <- lapply(strsplit(as.character(bamList$bamName),"\\."),"[",1)
      bamList$category <- "modern"
      #bamList[grep("downsamp",bamList$bamName),]$category <- "modern-downsampled"
      bamList[grep("^A",bamList$sampleID),]$category <- "ancient"
      # fill in reference:
      bamList$reference <- ref
      bamList$pop <- lapply(strsplit(as.character(bamList$sampleID),"_"),"[",3)
      #######################################################################
      # get eigen values from covariance matrix
      # Eigenvalues
      eig <- eigen(covar, symm=TRUE);
      # contains eigenvectors and eigenvalues
      eig$val <- eig$val/sum(eig$val); # rescale to be proportional out of 1
      #cat(signif(eig$val, digits=3)*100,"\n"); # rounds eigenvalues and makes percentages
      # eigen$val is a set of numInd values that are how much variance is explained (?)
      # eigen$vectors is an numInd x numInd matrix that has all your PCs in it! Each PC is a column and
      # rows are individuals
      # so for each pair of PCs you can pull out the 
      # Plot
      PC <- as.data.frame(eig$vectors) # the PCs are the eigenvectors; the eigenvalues are how much each of those vectors explains of the variance, I think. So if you want PC3 and PC4 then it the x,y coordinates for each individual found in those two columns 
      colnames(PC) <- gsub("V", "PC", colnames(PC)) # rename columns to be PC1 instead of V1
      # add metadata. This assumes that the individuals in your bamlist used in angsd are in the same order as the rows of the PCs
      PC$sampleID <- factor(unlist(bamList$sampleID))
      PC$reference <- factor(unlist(bamList$reference))
      PC$category <- factor(unlist(bamList$category))
      PC$population <- factor(unlist(bamList$pop))
      
      comp=c(1,2)
      x_axis = paste("PC",comp[1],sep="")
      y_axis = paste("PC",comp[2],sep="")
      
      title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)\nMapped to ",ref,"\n",state,sep="",collapse="")
      
      p12 <- ggplot(PC, aes_string(x=x_axis, y=y_axis)) + 
        geom_point(aes(color=population,shape=category),size=8) +
        ggtitle(title)+
        theme_bw()+
        geom_text(aes(label=sampleID),stat="identity",size=3,position="jitter")+
        scale_shape_manual(values=c(8,16,1))
      p12
      ggsave(paste(data.dir,"PCA.mappedTo",ref,".",state,".minMaf",minMaf,".PC1.PC2.pdf",sep=""),p12,device="pdf",height=7, width=10)
      
      comp=c(1,3)
      x_axis = paste("PC",comp[1],sep="")
      y_axis = paste("PC",comp[2],sep="")
      
      title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)\nMapped to ",ref,"\n",state,sep="",collapse="")
      
      p13 <- ggplot(PC, aes_string(x=x_axis, y=y_axis)) + 
        geom_point(aes(color=population,shape=category),size=8) +
        ggtitle(title)+
        theme_bw()+
        geom_text(aes(label=sampleID),stat="identity",size=3)+
        scale_shape_manual(values=c(8,16,1))
      p13
      ggsave(paste(data.dir,"PCA.mappedTo",ref,".",state,".minMaf",minMaf,".PC1.PC3.pdf",sep=""),p13,device="pdf",height=7, width=10)
      
      
      
      comp=c(2,3)
      x_axis = paste("PC",comp[1],sep="")
      y_axis = paste("PC",comp[2],sep="")
      
      title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)\nMapped to ",ref,"\n",state,sep="",collapse="")
      
      p23 <- ggplot(PC, aes_string(x=x_axis, y=y_axis)) + 
        geom_point(aes(color=population,shape=category),size=8) +
        ggtitle(title)+
        theme_bw()+
        geom_text(aes(label=sampleID),stat="identity",size=3)+
        scale_shape_manual(values=c(8,16,1))
      p23
      ggsave(paste(data.dir,"PCA.mappedTo",ref,".",state,".minMaf",minMaf,".PC2.PC3.pdf",sep=""),p23,device="pdf",height=7, width=10)
    }
  }
}
}

################################# high coverage ###############################
angsddates=c("20190701-lowcov-AFprior-MajorMinor4") # high cov only
refs=c("elut","mfur") # set what reference you mapped to 
states=c("1e-06.snpsOnly", "1e-06.snpsOnly.TransvOnly")
minMafs=c("0.12","0.2","0.05")
# see below for downsampled inds only 
# with and without transversions
for(angsddate in angsddates){
  data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/PCA/",angsddate,"/",sep="")
  
  for(ref in refs){
    # ooh individuals may be in different orders! 
    
    #bamList <- read.table(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/angsd.bamList.mappedto",ref,"fullpaths.txt",sep=""),stringsAsFactors = F)
    for(state in states){
      for(minMaf in minMafs){
        ####### Plot PCA from covariance matrix #########
        # following example in https://github.com/mfumagalli/ngsPopGen/blob/master/scripts/plotPCA.R
        ### test: read in binary covariance matrix:
        # read in covariance matrix
        covar <- npyLoad(paste(data.dir,"pcAngsd.",ref,".",state,".minMaf.",minMaf,".cov.npy",sep="")) # ready from binary ## NOTE; this test dataset is only based on a tiny part of a chromosome. Only 4 snps! I think this test was mapped to elut. let's see
        # read numpy binary arrays  in R 
        # http://dirk.eddelbuettel.com/code/rcpp.cnpy.html
        
        #dim(covar) # 15 x 15. so it's all 15 individuals with a covariance matrix between them 
        
        ####### assign information from list of bamfiles ########
        # THIS IS A USEFUL BIT OF CODE THAT PULLS SAMPLE ID, REFERENCE AND WHETHER WAS DOWNSAMPLED
        # from list of bam names that was used in angsd
        # the order of these individuals will match the rows of the PCs 
        bamList=LCSampleList # LOW COVERAGE
        colnames(bamList) <- "sampleID"
        #bamList$bamName <- lapply(strsplit(bamList$bamPath,"/"),tail,n=1) # tail pulls out last entry. useful because the different paths have different numbers of dirs so it's not always the same index # but it is the last entry always
        #bamList$sampleID <- lapply(strsplit(as.character(bamList$bamName),"\\."),"[",1)
        bamList$category <- "modern"
        #bamList[grep("downsamp",bamList$bamName),]$category <- "modern-downsampled"
        bamList[grep("^A",bamList$sampleID),]$category <- "ancient"
        # fill in reference:
        bamList$reference <- ref
        bamList$pop <- lapply(strsplit(as.character(bamList$sampleID),"_"),"[",3)
        #######################################################################
        # get eigen values from covariance matrix
        # Eigenvalues
        eig <- eigen(covar, symm=TRUE);
        # contains eigenvectors and eigenvalues
        eig$val <- eig$val/sum(eig$val); # rescale to be proportional out of 1
        #cat(signif(eig$val, digits=3)*100,"\n"); # rounds eigenvalues and makes percentages
        # eigen$val is a set of numInd values that are how much variance is explained (?)
        # eigen$vectors is an numInd x numInd matrix that has all your PCs in it! Each PC is a column and
        # rows are individuals
        # so for each pair of PCs you can pull out the 
        # Plot
        PC <- as.data.frame(eig$vectors) # the PCs are the eigenvectors; the eigenvalues are how much each of those vectors explains of the variance, I think. So if you want PC3 and PC4 then it the x,y coordinates for each individual found in those two columns 
        colnames(PC) <- gsub("V", "PC", colnames(PC)) # rename columns to be PC1 instead of V1
        # add metadata. This assumes that the individuals in your bamlist used in angsd are in the same order as the rows of the PCs
        PC$sampleID <- factor(unlist(bamList$sampleID))
        PC$reference <- factor(unlist(bamList$reference))
        PC$category <- factor(unlist(bamList$category))
        PC$population <- factor(unlist(bamList$pop))
        
        comp=c(1,2)
        x_axis = paste("PC",comp[1],sep="")
        y_axis = paste("PC",comp[2],sep="")
        
        title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)\nMapped to ",ref,"\n",state,sep="",collapse="")
        
        p12 <- ggplot(PC, aes_string(x=x_axis, y=y_axis)) + 
          geom_point(aes(color=population,shape=category),size=8) +
          ggtitle(title)+
          theme_bw()+
          geom_text(aes(label=sampleID),stat="identity",size=3,position="jitter")+
          scale_shape_manual(values=c(8,16,1))
        p12
        ggsave(paste(data.dir,"PCA.mappedTo",ref,".",state,".minMaf",minMaf,".PC1.PC2.pdf",sep=""),p12,device="pdf",height=7, width=10)
        
        comp=c(1,3)
        x_axis = paste("PC",comp[1],sep="")
        y_axis = paste("PC",comp[2],sep="")
        
        title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)\nMapped to ",ref,"\n",state,sep="",collapse="")
        
        p13 <- ggplot(PC, aes_string(x=x_axis, y=y_axis)) + 
          geom_point(aes(color=population,shape=category),size=8) +
          ggtitle(title)+
          theme_bw()+
          geom_text(aes(label=sampleID),stat="identity",size=3)+
          scale_shape_manual(values=c(8,16,1))
        p13
        ggsave(paste(data.dir,"PCA.mappedTo",ref,".",state,".minMaf",minMaf,".PC1.PC3.pdf",sep=""),p13,device="pdf",height=7, width=10)
        
        
        
        comp=c(2,3)
        x_axis = paste("PC",comp[1],sep="")
        y_axis = paste("PC",comp[2],sep="")
        
        title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)\nMapped to ",ref,"\n",state,sep="",collapse="")
        
        p23 <- ggplot(PC, aes_string(x=x_axis, y=y_axis)) + 
          geom_point(aes(color=population,shape=category),size=8) +
          ggtitle(title)+
          theme_bw()+
          geom_text(aes(label=sampleID),stat="identity",size=3)+
          scale_shape_manual(values=c(8,16,1))
        p23
        ggsave(paste(data.dir,"PCA.mappedTo",ref,".",state,".minMaf",minMaf,".PC2.PC3.pdf",sep=""),p23,device="pdf",height=7, width=10)
      }
    }
  }
}

