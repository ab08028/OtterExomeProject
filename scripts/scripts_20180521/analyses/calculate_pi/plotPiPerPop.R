input=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PI_THETA/empirical/20181119/pi.S.Theta.CalculatedFromSFSes.allpops.TajimasD.txt",header=T)
input

require(ggplot2)
require(reshape2)
genotypeDate="20181119"
model="1D.2Epoch"
mu=8.64e-09
require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[1],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5])
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PI_THETA/empirical/20181119/"
# arrange levels:
input$population <- factor(input$population,levels=c("CA","AK","AL","COM","KUR"))

p1 <- ggplot(input,aes(x=population,y=pi_perSite,fill=population))+
  geom_col(size=5)+
  theme_bw()+
  xlab("")+
  ylab(expression("neutral "*pi*" per site"))+
  scale_color_manual(values=unlist(colors))+
  theme(legend.position="none",text=element_text(size=14))+
  scale_y_continuous(labels=function(x) format(x,scientific=T))

p1
ggsave(paste(plot.dir,"neutral.pi.per.site.per.pop.pdf",sep=""),height=3,width=4)
