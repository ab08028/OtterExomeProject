##### Faststructure outputs SVG plots 
# download them and copy to : /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/FASTSTRUCTURE
wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/FASTSTRUCTURE/"
calldate=20180806 # date genotypes were called
# install.packages("rsvg")
# install.packages("svglite")
require(rsvg)
#require(svglite)
# convert to PDF
for (i in seq(1,10)){
  print(i)
  svg1File <- paste(wd,calldate,"/plots/snp_5_passingAllFilters_postMerge_raw_variants.faststructure_plot.",i,".svg",sep="")
  rsvg_pdf(svg1File, paste(wd,calldate,"/plots/snp_5_passingAllFilters_postMerge_raw_variants.faststructure_plot",i,".pdf",sep=""))
}
