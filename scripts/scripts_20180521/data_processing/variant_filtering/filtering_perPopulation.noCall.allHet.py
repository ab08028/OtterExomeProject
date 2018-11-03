# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:45:50 2018

@author: annabelbeichman
# These filters should be carried out on the per-population VCF and will check for:

1. sites where there are any no-call genotypes (removed); this stringency can be adjusted
2. sites where all called genotypes are heterozygotes (removed) (this is no longer done across all populations at once)

usage: use python 2.7
python filtering_bespokeFiltersAndChecks.py [infile full path] [outfile full path] [error file full path]
"""
import sys
import gzip
#import re
#import time
#from collections import Counter
# sys argv 1 = filename
filepath = sys.argv[1] #input file
outname= sys.argv[2] # output file
errorname= sys.argv[3] # error file
maxNoCallFrac=sys.argv[4] # this is fraction of no-call genotypes you'll permit per site (anything greater than that will be excluded)
# example: 118 samples, 118*.2 = 23.6 so any site with > 24 no-call sites will be excluded. A site with 23 no-call will be kept.
#
#filepath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/dummyVCF.forSandbox.allSites_5_passingFilters.vcf.gz" # this was my dummy file that I used to test (had some artifically bad sites for testing)
#outname="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/6_allSites_bespokeFilters_passingFilters_80percCall_elut.raw_variants.20170914.vcf"
#errorname="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/sitesFailing_bespokeFilters.vcf"

# useful commands that are now part of function so dont' need to use:
#inVCF = gzip.open(filepath, 'r')
#outVCF = open(outname, 'w')
#errorVCF = open(errorname, 'w')


################ write header lines to new vcf file #############
# this then prints out the rest of the header lines
def main_vcf_check(inputvcfilename,outfilename,errorfilename):
    # open your files:
    inVCF = gzip.open(inputvcfilename, 'r')
    outVCF = open(outfilename, 'w')
    errorVCF = open(errorfilename, 'w')
    # set up your possible ref and alt bases
    #refbases=set(['A', 'C', 'G', 'T'])
    #altbases=set(['A', 'C', 'G', 'T', '.'])
    # set up your expected genotypes (change if phased!!!!)
    #expected_genos=set(['0/0', '0/1', '1/1', './.']) # this doesn't expect phase
    # set up a counter of no-pass (nop) or pass (p) sites
    counter_nop=0 # count no pass
    counter_p=0 # count pass
    # set up the header of the new vcf:
    for line0 in inVCF:
        if line0.startswith('#'):
            outVCF.write(line0)
            continue
        ### For all other non-header lines:
        line=line0.strip().split('\t') # this splits line by tabs
        #CHROM0	POS1	ID2	REF3	ALT4	QUAL5	FILTER	INFO	FORMAT	[indivudals]
        #myref=line[3]
        #myalt=line[4]
        #myqual=line[5]
        #myfilter=line[6]
        #myinfo=line[7]
        #infoFields=dict(s.split('=') for s in myinfo.split(";")) # this makes my info into a dictionary that you can query; e.g. infoFields['AN'] will give you the AN value
        #myformat=line[8]
        mygenoinfo=line[9:]
        # get genotype info
        allCalls=[i.split(":")[0] for i in mygenoinfo] # get genotype calls
        myHomRef=allCalls.count("0/0") + allCalls.count("0|0") # number of hom ref gts
        myHet=allCalls.count("0/1") + allCalls.count("1/0") + allCalls.count("0|1")+ allCalls.count("1|0") # number of het gts
        myHomAlt=allCalls.count("1/1") + allCalls.count("1|1") # num of hom alt gts
        myCalled=myHomRef+myHet+myHomAlt # num of called genotypes
        myNoCalled=allCalls.count("./.") # num of no called genotypes
        # get AC and AN values:
        #myAC=myHet+(2*myHomAlt) # alternate alleles
        #myAN=myCalled *2 # all called alleles
        # do a series of checks
        # check if there are more than X% no-called genotypes
        if float(myNoCalled) > round(float(len(allCalls))* float(maxNoCallFrac)):
            errorVCF.write('# More than ' + str(maxNoCallFrac) + ' fraction of genotypes are  no-call at this site\n')
            errorVCF.write(line0)
            counter_nop+=1
        # check if all calls are heterozygous at the population level
        elif myHet==myCalled:
            errorVCF.write('# All genotypes are heterozygous\n')
            errorVCF.write(line0)
            counter_nop+=1
        else:
            outVCF.write(line0)
            counter_p+=1
    # write out final information
    errorVCF.write('## This many lines failed a filter and werent printed' + '\t' + str(counter_nop) + '\n')
    errorVCF.write('## This many lines PASSED and were printed' + '\t' + str(counter_p) + '\n')
    errorVCF.close()
    outVCF.close()
    inVCF.close()
    
#### run the function ##########
main_vcf_check(filepath,outname,errorname)
# exit
sys.exit()
