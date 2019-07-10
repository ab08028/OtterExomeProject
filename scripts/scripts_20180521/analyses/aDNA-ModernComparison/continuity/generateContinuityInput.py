# -*- coding: utf-8 -*-
import gzip
import sys
"""
Created on Mon Jul  8 14:46:47 2019

@author: annabelbeichman

Want to take in a superfile that combines a bed formatted file of frequencies from modern data, and counts for each basepair from angsd for ancient data
and combined them

This assumes ***3*** ancient individuals -- you have to change the script if there are more/fewer
"""
filepath = sys.argv[1] #path to input file
#sampleIDFile=sys.argv[2] # path to file with list of names in SAME ORDER as bamList you used for angsd
outfilepath= sys.argv[2] # output file
outfilepathTV=sys.argv[3] # output file for transversions

####### for testing: #######
#filepath = "/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/generateContinuityInput/dummy.freqs.counts.gz"
#outfilepath= "/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/generateContinuityInput/test.out.txt"
#outfilepathTV= "/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/generateContinuityInput/test.out.TV.txt"
superfile = open(filepath,"r")
outfile = open(outfilepath,"w")
outfileTV = open(outfilepathTV,"w") # for transversions
# get header:
# for now spaces is okay -- eventually want to make it tabs
header=[]
for line in superfile:
    if "#" in line:
        header=line.strip().split('\t')

        break
### things we care about:
# REF_FREQ
        
markerIndex=header.index("markerID") # there are two of these, just case about the first one
# index of modern ref and alt freqs:
refFreqIndex=header.index("REF_FREQ")
altFreqIndex=header.index("ALT_FREQ")

scaffIndex=header.index("chr")
posIndex=header.index("pos")
countIndex=header.index("ind0_A") # where counts start
###### 3 ancient individuals --- change this if you have more/fewer #####
totInds=3
# eventually want to enter this better
# len(header[countIndex:])/4 # count up total counts columns and divide by 4 (A C T G bases) - that gives total inds. should be 3 ancient individuals

# order of bases :  'ind2_A', 'ind2_C', 'ind2_G', 'ind2_T']
bases=['A','C','G','T'] # this is the order of bases in dumpCounts 4
transversions=[('A','C'),('C','A'),('A','T'),('T','A'),('C','G'),('G','C'),('G','T'),('T','G')]

# need to get the A C T G stuff 
# set up header
outfile.write("Chrom\tPos\tAF\tind0_der\tind0_anc\tind0_other\tind1_der\tind1_anc\tind1_other\tind2_der\tind2_anc\tind2_other\n")
outfileTV.write("Chrom\tPos\tAF\tind0_der\tind0_anc\tind0_other\tind1_der\tind1_anc\tind1_other\tind2_der\tind2_anc\tind2_other\n")
# go through file
for line0 in superfile:
    # skip header 
    if line0.startswith("#"):
        continue
    # process line, split by tabs:
    line=line0.strip().split('\t')
    # First want to check if the line passes the nInd filter
    # But want that nInd filter to also account for the minDepth filter
    # so need a couple things:
    scaff=line[scaffIndex]
    pos=line[posIndex]
    ##### get ref and alt info (from modern data) ######
    refAllele=line[refFreqIndex].split(":")[0] # A T C or G
    refFreq=line[refFreqIndex].split(":")[1] # frequency (must be between 0 and 1)
    altAllele=line[altFreqIndex].split(":")[0] # A T C or G
    altFreq=line[altFreqIndex].split(":")[1] # frequency (must be between 0 and 1)
    #### get counts of ref and alt alleles and sum up the others #####
    # write these out:
    outfile.write("\t".join([scaff,pos,altFreq]))
    if tuple([refAllele,altAllele]) in transversions:
        outfileTV.write("\t".join([scaff,pos,altFreq]))
    #outfile.write("\t")
    counts=line[countIndex:]
    for ind in range(0,totInds):
        # get the 4 counts for a given individual from 0 --> 2 (0, 1, 2); 4 columsn per individual (so you do ind *4 = start, and ind*4  + 4 = end)
        indCounts = counts[ind*int(4):ind*int(4)+int(4)]
        indCountsDict = dict(zip(bases,indCounts))
        refCount=indCountsDict[refAllele]
        altCount=indCountsDict[altAllele]
        # sum up all counts
        totalCount=sum([int(x) for x in indCountsDict.values()])
        # to get other count, subtract refCount and altCount from totalCounts -- this is the count of non ref / non alt reads:
        otherCount=totalCount - int(refCount) - int(altCount)

        # make sure you FIRST write out *alt* (derived) count THEN *ref* (ancestral) allele (according to continuity guidelines)
        out=["\t",str(altCount),"\t",str(refCount),"\t",str(otherCount)]
        outfile.write("".join(out))
        if tuple([refAllele,altAllele]) in transversions:
            outfileTV.write("".join(out))
            
    outfile.write("\n")
    if tuple([refAllele,altAllele]) in transversions:
        outfileTV.write("\n")
        
    # okay so have the counts, but need to figure out the 

        
# reset file:
superfile.seek(0)
superfile.close()
outfile.close()
outfileTV.close()
# did some checks on a couple sites, seems to work. 