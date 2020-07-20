# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:45:50 2018

@author: annabelbeichman

usage: python 2.7
python filtering_getNoCallPerInd.py [infile full path] [outfile full path]
"""
import sys
import gzip
### making this on 20200720 to get total het sites per individual
# want to output a table with sample id, total number of non missing sites, total hets, total monomorphic, total hom-alt and total hom-ref
# this should input a all-sites vcf file (snps and monomorphic)
filepath = sys.argv[1] #input file
outname= sys.argv[2] # output file
#errorname= sys.argv[3] # error file
# these are filepaths for dummy testing:
#filepath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/dummyVCF.forSandbox.allSites_5_passingFilters.vcf.gz" # this was my dummy file that I used to test (had some artifically bad sites for testing) <-- NOTE THIS HAS ARTIFICALLY BAD SITES FOR TESTING PURPOSES (!!!) and is NOT representative.
#outname="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/noCallPerInd.getFromFilteredVCF.txt"
#outfilename="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/hetPerIndsandbox.output.txt"
# do i want that from filtered vcf? then I have to refilter a bunch.
# What about also getting DP dist at same time
# outfile
# sample meanDP totalNoCall
#errorname="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/sitesFailing_bespokeFilters.vcf"

# useful commands that are now part of function so dont' need to use:
# inVCF = gzip.open(filepath, 'r')

################ this function will count the number of no-call genotypes per individual for a VCF ##############
def getHetPerInd(inputvcfilename,outfilename):
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
    calledSiteDict=dict()
    hetSiteDict=dict()
    monomorphicSiteDict=dict()
    homAltDict=dict()
    homRefDict=dict()
    # note homRef +homALt = monomorphic and monomorphic + het = all called sites
    for sample in samples:
        calledSiteDict[sample]=0
        hetSiteDict[sample]=0
        monomorphicSiteDict[sample]=0
        homAltDict[sample]=0
        homRefDict[sample]=0
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
            # count past missing data
            if callDict[sample]!="./.":
                calledSiteDict[sample]+=1
            # count hets
            if callDict[sample]=="0/1" or callDict[sample]=="1/0" or callDict[sample]=="0|1" or callDict[sample]=="1|0":
                hetSiteDict[sample]+=1
            # count hom ref
            if callDict[sample]=="0/0" or callDict[sample]=="0|0":
                # add to monomorphic and to hom ref categories
                monomorphicSiteDict[sample]+=1
                homRefDict[sample]+=1
            if callDict[sample]=="1/1" or callDict[sample]=="1|1":
                # add to monomorphic and to hom alt categories                
                monomorphicSiteDict[sample]+=1
                homAltDict[sample]+=1
   # (calledSiteDict)
   # print(hetSiteDict)
            
    inVCF.close()
    outList.write('Sample\tTotalCallCount\tHetSiteCount\tTotalMonomorphicCount\tHomRefCount\tHomAltCount\tHeterozygosityPerSite\n')
    #for sample, nocall in noCallDict.items():
    #    outList.write('{}\t{}\n'.format(sample, nocall))
    for sample in samples:
        outList.write('{}\t'.format(sample))
        outList.write('{}\t'.format(calledSiteDict[sample]))
        outList.write('{}\t'.format(hetSiteDict[sample]))
        outList.write('{}\t'.format(monomorphicSiteDict[sample]))
        outList.write('{}\t'.format(homRefDict[sample]))
        outList.write('{}\t'.format(homAltDict[sample]))
        # het is hets / total called sites
        heterozygosity=float(hetSiteDict[sample])/float(calledSiteDict[sample])
        outList.write('{}\n'.format(heterozygosity))
    
    outList.close()



#### run the function ##########
getHetPerInd(filepath,outfilename)


sys.exit()
