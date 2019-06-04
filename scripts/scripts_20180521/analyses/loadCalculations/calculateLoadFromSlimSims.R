
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/concattedSummaries/"
pops=c("AK","AL")
models=c("1D.2Epoch.1.5Mb.cds")
simdates=c(20190424)
#reps=c(seq(1,23))
reps=c(seq(1,25)) # some reps don't make it through Hoffman; so I have a file.exists() test in the loop to skip reps that didn't yield output
hset=c(0,0.5)
states=c("PreContraction","PostContraction")
allLoads=data.frame()

for(date in simdates){
  for(pop in pops){
    for(rep in reps){
      for(state in states){
        for(h in hset){
          # check if rep exists (some have random hoffman failures)
          infile=paste(data.dir,pop,"/",model,"/",simdate,"/h_",h,"/replicate_",rep,".slim.output.",state,".allConcatted.summary.txt.gz",sep="")
          if(file.exists(infile)){
            input = read.table(paste(data.dir,pop,"/",model,"/",simdate,"/h_",h,"/replicate_",rep,".slim.output.",state,".allConcatted.summary.txt.gz",sep=""),sep=",",header=T)
            # calculate q (alt allele frequency) and p (ref allele frequency) per site
            # be careful about which you use in equation! s*q^2 means that q is frequency of ALT allele with associated "s". so p is freq of ref allele
            input$qFreq <- (input$p1numhet + (2*input$p1numhom)) / (2*input$popsizeDIP)
            input$pFreq <- 1 - input$qFreq
            # equation is, per site: 2hspq + sq^2 is the contribution to load; p = 1-qFreq
            # my "s" is negative, so I want to absolute value s --> |s|
            input$loadComponent <- (2*h*abs(input$s)*input$qFreq*input$pFreq) + (abs(input$s)*((input$qFreq)^2))
            # total sites:
            
            S = sum(input$loadComponent)
            W = exp(-S) # mean fitness e^-S
            L  = 1 - W # mutation load 
            #### add to dataframe: #####
            loadDF <- data.frame(population=pop)
            loadDF$rep <- rep
            loadDF$state <- state
            loadDF$h <- h
            loadDF$da
            loadDF$S_allsites <- S
            loadDF$W_meanFitness <- W
            loadDF$L_mutationLoad <- L
            #### combine with other reps: #####
            allLoads = rbind(allLoads,loadDF)
          }}}}}}

# change order of factors:
allLoads$state <- factor(allLoads$state,levels=c("PreContraction","PostContraction"))
# label H:
allLoads$hLabel <- paste("h = ",allLoads$h)
ggplot(allLoads,aes(x=state,y=L_mutationLoad,fill=state))+
  geom_boxplot(position=position_dodge(.5))+
  #geom_point(position=position_dodge(.5),size = 1,alpha=0.75)+
  theme_bw()+
  facet_wrap(population~hLabel)
