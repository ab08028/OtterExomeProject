
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 15:27:15 2018

@author: annabelbeichman
"""
import sys
import gzip
import datetime
import argparse
todaysdate=datetime.datetime.today().strftime('%Y%m%d')

############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Generate a 1D SFS in dadi format and in easy input format for R based on a pre-filtered VCF that contains ONLY the regions you of interest')
parser.add_argument("--vcf",required=True,help="path to vcf file")
parser.add_argument("--pop",required=True,help="population identifier, e.g. 'CA'")
parser.add_argument("--outdir",required=True,help="path to output directory")
parser.add_argument("--outPREFIX",required=False,help="output file prefix (optional)",default="")

args = parser.parse_args()
vcfFile=args.vcf
pop=str(args.pop)
outdir=str(args.outdir)
prefix=str(args.outPREFIX)

#### OPEN VCF TO READ ######### 
inVCF = gzip.open(vcfFile, 'r')
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

##### get sample size and set up sfs ###########
totalChr=len(samples)*2
unfoldedSFSlen=totalChr+1 # adds in a monomorphic category because could show up on 0,1,2...,2n
sfs=dict()
for i in range(0,unfoldedSFSlen):
    sfs[i]=0


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
    myHomRef=allCalls.count("0/0") + allCalls.count("0|0") # number of hom ref gts
    myHet=allCalls.count("0/1") + allCalls.count("1/0") + allCalls.count("0|1")+ allCalls.count("1|0") # number of het gts
    myHomAlt=allCalls.count("1/1") + allCalls.count("1|1") # num of hom alt gts
# get counts of all called genotypes:    
    myCalled=myHomRef+myHet+myHomAlt # num of called genotypes
# get the counts of no-call genotypes:
    myNoCalled=allCalls.count("./.") # num of no called genotypes
# set AC and AN values: (note that GATK doesn't always reliably update AN and AC, so we do it ourselves. AC is the count of alternate alleles, AN is the total called alleles
    myAC=myHet+(2*myHomAlt) # alternate alleles
    myAN=myCalled *2 # all called alleles
    ############ Fill in sfs ################
    # make sure that there are no "./." genotypes:
    if myNoCalled==0:
        sfs[myAC]+=1 # add the snp to the appropriate bin of SFS; this works by accessing the key of the dict that corresponds to the allele frequency (AC), e.g. myAC of 30 is the bin in the unfolded
        # this version of the script allows there to be no-call sites in the sfs, but they won't go into the SFS (contrast to generateSFS.py which will error out if there are no-call because that means something has gone wrong with your filtering)
    else:
        continue
    # SFS to add +1 to.
    # print(sum(sfs.values()))
############################ write out SFS in dadi format ################
# write it out in a couple formats (dadi, fsc, R) #### dadi format : 
outputFile=open(str(outdir)+"/"+pop+"."+prefix+".unfolded.sfs.dadi.format."+todaysdate+".txt","w")
outputFile.write(str(unfoldedSFSlen)+" unfolded "+"\""+str(pop)+"\"\n")
#freqs='\t'.join(str(x) for x in sfs) # set up first row
values='\t'.join(str(x) for x in sfs.values()) # set up second row 
outputFile.write(values)
outputFile.close()
# skipping mask for now 

############################ eventual functionality: write out SFS in fastsimcoal format (or convert it in R) ################
#outputFile=open(str(outdir)+"/"+pop+".sfs_jointMAFpop1.obs","w")
# header=
#outputFile.close()
############################ write out SFS for easy plotting in R ################
outputFile=open(str(outdir)+"/"+pop+"."+prefix+".unfolded.sfs.R.format."+todaysdate+".txt","w")
outputFile.write("frequency\tcount\n")
[outputFile.write('{0}\t{1}\n'.format(key,value)) for key, value in sfs.items()]
outputFile.close()

sys.exit()