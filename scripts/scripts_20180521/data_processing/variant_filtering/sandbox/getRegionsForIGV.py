# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:45:50 2018

@author: annabelbeichman
# This script will pull out singletons +- 50bp from a population level filtered vcf file
and make a samtools script so that you can pull them out of the apporpriate bam file
for igv visualization

usage: use python 2.7
python getRegionsForIGV.py [infile full path] [outfile script path] [path to paleomix dir where bams are] [path to where you want new bams to go] 
"""
import sys
import gzip
#import re
#import time
#from collections import Counter
# sys argv 1 = filename
filepath = sys.argv[1] #input file
outname= sys.argv[2] # output file
bamdir= sys.argv[3] # where bam  are located (upper level)
igvdir = sys.argv[4] # where you want new bams to go
#filepath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/dummyVCF.forSandbox.allSites_5_passingFilters.vcf.gz" # this was my dummy file that I used to test (had some artifically bad sites for testing)
# useful commands that are now part of function so dont' need to use:
inVCF = gzip.open(filepath, 'r')
#putVCF = open(outname, 'w')
#errorVCF = open(errorname, 'w')

# get sample names
samples=[]
for line in inVCF:
	if line.startswith('##'):
		pass
	else:
		for i in line.split()[9:]: samples.append(i)
		break
################ write header lines to new vcf file #############
# this then prints out the rest of the header lines
  # you must give input vcf file, output script file name, bamdir and igvdir where you want new files to go
def selectForIGV(inputvcfilename,outfilename,bamdir,igvdir):
    # open your files:
    # these can be changed:
    #bamdir="/u/flashscratch/a/ab08028/captures/paleomix/"
    #outdir="/u/flashscratch/a/ab08028/captures/IGV/"
    inVCF = gzip.open(inputvcfilename, 'r')
    script = open(outfilename, 'w')
    # set up script :
    script.write("module load samtools\n")
    script.write("bamDir="+bamdir+"\n")
    script.write("igvdir="+igvdir+"\n")
    for line0 in inVCF:
        if line0.startswith('#'):
            continue
        ### For all other non-header lines:
        line=line0.strip().split('\t') # this splits line by tabs
        mygenoinfo=line[9:]
        # get genotype info
        allCalls=[i.split(":")[0] for i in mygenoinfo] # get genotype calls
        allAD=[i.split(":")[1] for i in mygenoinfo] # get genotype depths should be in format X,Y; not using this for now, but will in the future so leaving it in
        #myHomRef=allCalls.count("0/0") + allCalls.count("0|0") # number of hom ref gts
        myHet=allCalls.count("0/1") + allCalls.count("1/0") + allCalls.count("0|1")+ allCalls.count("1|0") # number of het gts
        myHomAlt=allCalls.count("1/1") + allCalls.count("1|1") # num of hom alt gts
       # get AC and AN values:
        myAC=myHet+(2*myHomAlt) # alternate alleles
        #myAN=myCalled *2 # all called alleles
        # if it's a singleton (assumes no no-calls)
        if myAC==1 or myAC==(len(allCalls)-1):
            print("found a singleton!"+"AC = " + myAC)
            samplesCalls=dict(zip(samples,allCalls))
            # gets you name of heterozygous individual(s)
            keys=[key for key,value in samplesCalls.items() if value=="0/1"]
            for key in keys:
                bampath=str(bamdir+key+"/*bam")
                samtoolsEntry=("samtools view -b " + bampath+" \""+str(line[0])+":"+str(int(line[1])-50)+"-"+str(int(line[1])+50)+"\" >  "+igvdir+key+"_"+str(line[0])+"-"+str(line[1])+".bam"+"\n")
                script.write(samtoolsEntry)
        else:
            continue
        ##### Set some conditions
    script.close()
    inVCF.close()


#### run the function ##########
selectForIGV(filepath,outname,bamdir,igvdir)

sys.exit()
