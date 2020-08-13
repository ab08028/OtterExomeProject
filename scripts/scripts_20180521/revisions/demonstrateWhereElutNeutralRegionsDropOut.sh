############## # goal is to figure out why there is a 3Mb difference between
# elut target regions and mfur neutral regions
# MK thinks its because of non-mapping
# I think they map, but due to superior annotation of mfur they aren't passing neutral criteria based on that annotation

# a few ways to determine what's going on
# in the past, I blasted the elut regions against the mfur genome. I required a min overlap of 500bp
# and I got a bed file. 
# let's see what's going on with that it's on sirius:
# extract a single file:
tar -zxvf round5ExomeCaptureDesign_Oct2016.20200512.tar.gz round5ExomeCaptureDesign_Oct2016/round5_AED1.0_newnames/findFerretCaptureRegions_20180717/blastResults/blast.seaOtterNeutralCaptureRegions_vs_ferretGenome.min500bp.bed