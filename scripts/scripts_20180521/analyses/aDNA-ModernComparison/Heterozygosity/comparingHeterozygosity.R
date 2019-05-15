### R plot pis from parsing Beagle:
# this is only based on 200K sites
require(ggplot2)
require(reshape2)
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/Heterozygosity"
input <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/parseBeagleFile/testout.txt",header=T,sep="\t")
input$label <- "modern"
input[grep("^A",input$sample),]$label <- "ancient"
input[grep("downsamp",input$sample),]$label <- "modern-downsampled"

input_melt <- melt(input,measure.vars = c("HetPerSite","HetPerSite_transversionsOnly"),id.vars = c("sample","label"))

# label transitions+transv vs transv only
input_melt$label2 <- "Transitions+Transversions"
input_melt[input_melt$variable=="HetPerSite_transversionsOnly",]$label2 <- "TransversionsOnly"

p0 <- ggplot(input_melt,aes(x=sample,y=value,fill=label2))+
  geom_bar(stat="identity",position="dodge")+
  coord_flip()+
  theme_bw()
p0
ggsave(paste(data.dir,""))
p1 <- ggplot(input_melt,aes(x=label,y=value,fill=label2))+
  geom_violin(position=position_dodge(.5))+
  geom_point(position=position_dodge(.5),size = 1)+
  theme_bw()+
  ggtitle("Comparing Heterozygosity")+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),legend.text = element_text(size=14),legend.position = c(0.6,0.8),legend.background = element_rect(fill = "transparent"))+
  xlab("")
p1

p2 <- ggplot(input_melt[input_melt$variable=="HetPerSite_transversionsOnly",],aes(x=label,y=value,fill=label))+
  geom_violin(position=position_dodge(.5))+
  geom_point(position=position_dodge(.5),size = 1)+
  theme_bw()+
  ggtitle("Comparing Heterozygosity -- Transversions Only")+
  theme(legend.title = element_blank(),axis.text = element_text(size=14))+
  xlab("")
p2
