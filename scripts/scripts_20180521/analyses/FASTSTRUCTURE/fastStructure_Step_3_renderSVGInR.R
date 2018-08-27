##### Faststructure outputs SVG plots 
# download them and copy to : /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/FASTSTRUCTURE
<<<<<<< HEAD
calldate=20180806 # date genotypes were called
wd=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/FASTSTRUCTURE/",calldate,"/",sep="")
outdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/FASTSTRUCTURE/",calldate,"/",sep="")
dir.create(outdir)
=======
#wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/FASTSTRUCTURE/"
calldate=20180806 # date genotypes were called
wd="/u/flashscratch/a/ab08028/captures/analyses/FASTSTRUCTURE/${calldate}_filtered/"
>>>>>>> 71dab9a1455625836fb70d98d0b8e66695b36603
# install.packages("rsvg")
# install.packages("svglite")
require(rsvg)
#require(svglite)
# convert to PDF
for (i in seq(1,10)){
  print(i)
  svg1File <- paste(wd,"/svg/snp_7_downSampCOM_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.faststructure_plot.",i,".svg",sep="")
  rsvg_pdf(svg1File, paste(outdir,"/snp_7_downSampCOM_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.faststructure_plot.",i,".pdf",sep=""))
}
