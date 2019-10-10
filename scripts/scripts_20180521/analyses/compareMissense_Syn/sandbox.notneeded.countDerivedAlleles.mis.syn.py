# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 13:34:45 2019

@author: annabelbeichman
"""
import argparse
import gzip

############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Count number of derived alleles per individual from vcf file (missense or syn sites only in file)')
parser.add_argument("--vcf",required=True,help="path to vcf file containing only missense or syn. variants")
parser.add_argument("--outfile",required=True,help="path to output file")
################ testing params ############
# hoffman vcf location /u/flashscratch/a/ab08028/captures/vcf_filtering/20181119_filtered/neutral_and_cds_VCFs/cdsVCFs

vcf="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/troubleshooting/sandboxCountDerived/test.syn.vcf.gz"
outfilename="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/troubleshooting/sandboxCountDerived/syn.count.txt"
# want to count per individuals
#### OPEN A VCF TO READ ######### 
inVCF = gzip.open(vcf, 'r')
############# reset a vcf that you've read through #########
inVCF.seek(0)
#### OPEN A VCF TO WRITE OUT TO ######### 
outfile = open(outfilename, 'w')

########### GET SAMPLE NAMES #############
# get sample names
samples=[]
for line in inVCF:
	if line.startswith('##'):
		pass
	else:
		for i in line.split()[9:]: samples.append(i)
		break
####### set up a dictionary of samples #####
derivedAllelesPerSample = dict()
for sample in samples:
    derivedAllelesPerSample[sample]=0
###### READ THROUGH VCF AND EXTRACT INFO LINE BY LINE #######
# first read the header lines ("#") and write out as the new header
inVCF.seek(0)
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
    # convert to count of derived alleles -- 1/1 = 2; 0/1 = 1; 0/0 = 0 ; ./.=0 
    
    allCallsDict = dict(zip(samples,allCalls))
    for sample in samples:
        if allCallsDict[sample]=="0/0" or allCallsDict[sample]=="./.":
            continue
        elif allCallsDict[sample]=="0/1": # add one for hets
            derivedAllelesPerSample[sample] = derivedAllelesPerSample[sample] +1
        elif allCallsDict[sample]=="1/1": # add two for hom-alt
            derivedAllelesPerSample[sample] = derivedAllelesPerSample[sample] +2
## write out ###
outfile.write("id\tcount\n")
for id, value in derivedAllelesPerSample.iteritems():
    outfile.write(('{0}\t{1}\t').format(id,value))
    outfile.write("\n")
outfile.close()
inVCF.close()