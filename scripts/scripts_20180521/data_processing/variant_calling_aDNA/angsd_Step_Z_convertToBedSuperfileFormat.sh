####### make a super angsd output file ##########
# want to combine beagle gprobs or gls with maf and count information
#then 
# Convert --> bed format (with all this info as extra columns) ############
# can then intersect with neutral bed file, cds bed file, convert to VEP format etc



# want to do this as part of running angsd ; just automatically do a conversion and gzip it 
# do I want to convert mafs too? maybe. need to figure out how do this for counts and mafs as well
# maybe can do it all? 

# maybe make a mega file in python that izips maf and counts together?

ref=mfur
basename=angsdOut.mappedTo${ref}
mafs=${basename}.mafs.gz
counts=${basename}.counts.gz
gprobs=${basename}.beagle.gprobs.gz
GLs=

# probs first:
# I think I can do it without intermediate files, just pipe to awk and convert to bed format

paste $gprobs $mafs $counts | # go into awk and rearrange to make it bed format with extra columns 
 