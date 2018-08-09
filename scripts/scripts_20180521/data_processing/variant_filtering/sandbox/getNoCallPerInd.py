# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:45:50 2018

@author: annabelbeichman

"""
import sys
import gzip

filepath = sys.argv[1] #input file
outname= sys.argv[2] # output file
#errorname= sys.argv[3] # error file
# these are filepaths for dummy testing:
#filepath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/dummyVCF.forSandbox.allSites_5_passingFilters.vcf.gz" # this was my dummy file that I used to test (had some artifically bad sites for testing)
outname="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/noCallPerInd.getFromFilteredVCF.txt"
# do i want that from filtered vcf? then I have to refilter a bunch.
# What about also getting DP dist at same time
# outfile
# sample meanDP totalNoCall
#errorname="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/sitesFailing_bespokeFilters.vcf"

# useful commands that are now part of function so dont' need to use:
# inVCF = gzip.open(filepath, 'r')

################ this function will count the number of no-call genotypes per individual for a VCF ##############
def getNoCallPerInd(inputvcfilename,outfilename):
    # open your files:
    inVCF = gzip.open(inputvcfilename, 'r')
    outList = open(outfilename, 'w')
    # get sample names
    samples=[]
    for line in inVCF:
	if line.startswith('##'):
		pass
	else:
		for i in line.split()[9:]: samples.append(i)
		break
    # reset vcf
    inVCF.seek(0)
    # set up an empty dictionary filled with zeroes to be a counter
    noCallDict=dict()
    for sample in samples:
        noCallDict[sample]=0
    # skip the header lines
    for line0 in inVCF:
        if line0.startswith('#'): 
            continue
        ### For all other non-header lines:
        line=line0.strip().split('\t') # this splits line by tabs
        #CHROM0	POS1	ID2	REF3	ALT4	QUAL5	FILTER	INFO	FORMAT	[indivudals]
        mygenoinfo=line[9:]
        # zip it together with samples
        # get genotype info
        ####### This is really useful: makes queriable dictionary 
        # of sample names and their genotypes
        allCalls=[i.split(":")[0] for i in mygenoinfo] # get genotype calls
        callDict = dict(zip(samples,allCalls))
        # now want to detect any no calls and add to a counter for each individual
        # if a sample has a no call genotype, want to add 1 to its counter
        for sample in samples:
            if callDict[sample]=="./.":
                noCallDict[sample]+=1
    inVCF.close()
    outList.write('Sample\tNoCallCount\n')
    for sample, nocall in noCallDict.items():
        outList.write('{}\t{}\n'.format(sample, nocall))





#### run the function ##########
getNoCallPerInd(filepath,outname)

# skipping a few things from clare's script
# not bothering to update PASS flag for no call genotypes
# not bothering to update alt allele for sites that become invariant to "." -- leaving as is
# but making sure the AC= . 
