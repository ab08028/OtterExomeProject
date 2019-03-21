# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 15:27:15 2018

@author: annabelbeichman
"""
import sys
import gzip
#import datetime
import argparse
# Goal of the script:
# To count up the number of heterozygous and homozygous sites per individual
# can then be run on the missense or synonymous vcfs 
# want to output:
# ind (pop) 11_count 01_count 00_count noCall_count
############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Total up the number of homozygous ref, homozygous alt, heterozygous and no-call genotypes per individual')
parser.add_argument("--vcf",required=True,help="path to vcf file")
# maybe don't need pop
#parser.add_argument("--pop",required=True,help="population identifier, e.g. 'CA'")
parser.add_argument("--outdir",required=True,help="path to output directory")
parser.add_argument("--outPREFIX",required=False,help="output file prefix (optional)",default="")

args = parser.parse_args()
vcfFile=args.vcf
#pop=str(args.pop)
outdir=str(args.outdir)
prefix=str(args.outPREFIX)

# dummy files for testing:
#vcfFile="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/dummyVCF.forSandbox.allSites_5_passingFilters.vcf.gz"
#outdir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/troubleshooting/gettingNS_SynCounts"
#prefix="test"
#### OPEN VCF TO READ ######### 
# first check if it's gzipped or not (from Tanya Phung)
def file_test(vcf_file):
	"""
    This function checks if the input VCF file is gzip or not.
    """
	if vcf_file.endswith('.gz'):
		return gzip.open, 'rb'
	else:
           return open, 'rb'
           
open_func, mode = file_test(vcf_file=vcfFile)
#inVCF = gzip.open(vcfFile, 'r')
# open the VCF using the appropriate function (open for not gzipped; gzip.open for gzipped)
inVCF=open_func(vcfFile,mode)
############# reset vcf to make sure no lines are missed #########
inVCF.seek(0)

########### GET SAMPLE NAMES #############
# get sample names
samples=[]
for line in inVCF:
	if line.startswith('##'):
		pass
	else:
		for i in line.split()[9:]: samples.append(i)
		break
######### set up empty dicts with samples as the keys and 0 as starting values ############
HomAltCount=dict.fromkeys(samples,0)
HomRefCount=dict.fromkeys(samples,0)
HetCount=dict.fromkeys(samples,0)
noCallCount=dict.fromkeys(samples,0)
linesProcessed=0
# so for each individual, am going to count the 1/1 0/0 and 0/1 values
###### READ THROUGH VCF AND EXTRACT INFO LINE BY LINE #######
# first read the header lines ("#") 
for line0 in inVCF:
    if line0.startswith('#'):
        continue

### For all other non-header lines, split along tabs to get each entry as a seprate entry in the list "line"
    line=line0.strip().split('\t') # this splits line by tabs
#CHROM0	POS1	ID2	REF3	ALT4	QUAL5	FILTER	INFO	FORMAT	[indivudals]
# get all the genotypes for all individuals, then split each individual's genotype info by ":" to get all the calls
    mygenoinfo=line[9:]
    allCalls=[i.split(":")[0] for i in mygenoinfo] # get genotype calls
# Get the counts of HomozygousREference, Heterozygous and Homozygous Alternate alleles (for now has all combos of genotypes; though if unphased most likely will only see 0/0 0/1 and 1/1)
    allCallsSamples=dict(zip(samples,allCalls))
    for sample in allCallsSamples.keys():
        if allCallsSamples[sample]=="0/0" or allCallsSamples[sample]=="0|0":
            HomRefCount[sample]+=1
            linesProcessed+=1
            print(linesProcessed)
        elif allCallsSamples[sample]=="1/1" or allCallsSamples[sample]=="1|1":
            HomAltCount[sample]+=1
            linesProcessed+=1
            print(linesProcessed)
        elif allCallsSamples[sample]=="0/1" or allCallsSamples[sample]=="1/0" or allCallsSamples[sample]=="1|0" or allCallsSamples[sample]=="0|1":
            HetCount[sample]+=1
            linesProcessed+=1
            print(linesProcessed)
        elif allCallsSamples[sample]=="./.":
            noCallCount[sample]+=1
            linesProcessed+=1
            print(linesProcessed)
        else:
            sys.exit("There are weird genotypes present!")

############################ Write out a table of counts ################

outputFile=open(str(outdir)+"/"+prefix+".countOfHomAltRefHet.txt","w")
outputFile.write("individual\tHomRefCount\tHomAltCount\tHetCount\tNoCallCount\n")
    
for ind in samples:
    dataLine=('\t'.join(str(x) for x in (ind,HomRefCount[ind],HomAltCount[ind],HetCount[ind],noCallCount[ind])))
    outputFile.write(dataLine + "\n")

sys.exit()
