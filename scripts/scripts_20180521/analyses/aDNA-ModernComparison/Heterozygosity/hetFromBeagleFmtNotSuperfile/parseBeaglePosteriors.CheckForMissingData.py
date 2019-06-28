"""
Created on Fri May 10 12:14:23 2019

@author: annabelbeichman

In this script, we want to parse a beagle file and count the number of callable sites per individual that pass posterior prob threshold

The script will calculate heterozygosity for all sites and just for transversions 

It sums up heterozygosity posterior probability for each individual and divides by callable GTs (following methods of Fages et al. 2019 in Cell)

Explantion: the script sets up a dictionary for your list of individuals (order of list MUST match the input bam file list in ANGSD)

It then goes through a BEAGLE file site by site and for each individual it checks if the maximum posterior for the individual's three genotypes is >= some threshold (e.g. 0.5) .

If it is, it adds the heterozygosity posterior probabilty into that individual's dictionary entry (heterozygosity numerator), and tallies the GT as a callable Site for that individual (heterozygosity denominator)
# It will also calculate heterozygosity for transversions only
# and then reports the total sites passing the threshold per individual, the total het probabilities, and divides the two (hets/total sites) for all hets or just transversions


usage: python script.py inputFilepath countsFile sampleIDFile outputFile MaxProbCutoff PerIndividualDepthMinimum
"""
import gzip
import sys
from itertools import izip

filepath = sys.argv[1] #path to input file, beagle posterior GLs from angsd; should contain transitions and transversions (.gprobs.gz)
countsFile=sys.argv[2] ## path to counts file from angsd; sites should be in exact same order as beagle posterior grobs file (results of -doCount 1 -dumpCounts 2)
sampleIDFile=sys.argv[3] # path to file with list of names in SAME ORDER as bamList you used for angsd
outname= sys.argv[4] # output file
MaxProbCutoff=float(sys.argv[5]) # 20180626 setting to float -- it's also a float further in the script from before; this the cutoff for a genotype to be considered callable
PerIndividualDepthMinimum=float(sys.argv[6]) # per individual min depth

################# list of possible transversions ###############
# In beagle format, nucleotides are labeled as numbers
# # Beagle codes: the allele codes as 0=A, 1=C, 2=G, 3=T
#  so the 8 possible transversions are:
            # 0-1 : A-C  
            # 1-0 : C-A 
            # 0-3 : A-T 
            # 3-0: T-A 
            # 1-2 : C-G 
            # 2-1 : G-C 
            # 2-3 : G-T 
            # 3-2: T-G 
transversions=[('0','1'),('1','0'),('0','3'),('3','0'),('1','2'),('2','1'),('2','3'),('3','2')]

######### for testing only (set these as arguments eventually) #######################
#filepath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/parseBeagleFile/miniSample-testPosteriors.beagle.gprobs.gz"

# list of samples   # check order super carefully!!! must be in same order as input bam list for angsd!!!!!!!
#sampleIDFile="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.BeCarefulOfOrder.txt" 

#outname="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/parseBeagleFile/testout.txt"
#MissingnessMaxProbCutoff=0.5 # if the max of the 3 probs is below this, discard; keep if >= to the cutoff

# new files to test:
#filepath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/Heterozygosity/experiments/compareGL-GP-MAF/20190523-experiment-GLvsGP-highcov/posteriorProbabilities/angsdOut.mappedTomfur.OrlandoSettings.beagle.gprobs.gz"
#sampleIDFile="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt"
#outname="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/parseBeagleFile/testout.CheckForMissing.txt"
#MissingnessMaxProbCutoff=0.5
#CountsFile="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/Heterozygosity/experiments/compareGL-GP-MAF/20190523-experiment-GLvsGP-highcov/posteriorProbabilities/angsdOut.mappedTomfur.OrlandoSettings.counts.gz"
# from dumpCounts 2 (gives per individual counts per site)


######################################################

# read samples into a list (strip \n from end of each one)
sampList = [line.rstrip('\n') for line in open(sampleIDFile)]
numInd=len(sampList)
# these are now in the proper order ### MAKE SURE IT'S ORDER YOUR BAMLIST WAS IN FOR ANGSD!!!!! OTHERWISE INDS WILL BE ASSIGNED INCORRECTLY

###### make empty dictionaries : ##############
#### counts of missing GTs: 
missingDict=dict()
for sample in sampList:
    missingDict[sample]=0
#### counts of called GTs: 
calledDict=dict()
for sample in sampList:
    calledDict[sample]=0
#### sums of heterozygosous posterior Probs
hetProbSumDict=dict()
for sample in sampList:
    hetProbSumDict[sample]=0
#### sums of transversions heterozygous posterior Probs
#### sums of heterozygosous posterior Probs
TransvOnly_HetProbSumDict=dict()
for sample in sampList:
    TransvOnly_HetProbSumDict[sample]=0
########### Open beagle GL posteriors file #############

beagle = gzip.open(filepath,"rb")
counts = gzip.open(countsFile,"rb")

# get beagle header:
header=[]
for line in beagle:
    if "marker" in line:
        header=line.strip().split('\t')
        # check that length of header is numInd*3 + 3 (marker allele1 allele2):
        len(header)==(numInd*3) + 3 # should be TRUE
        break
# reset file:
beagle.seek(0)
counts.seek(0) # reset for good measure
# reading through the beagle and counts files at the same time; they are the same number of lines in the same order, which is why this works:
# using izip to zip together each line of each file into a pair 
# so that the two files can be read simultaneously
# https://stackoverflow.com/questions/8461154/use-izip-to-read-lines-from-two-files-simultaneously-in-python
for beagleline0,countline0 in izip(beagle,counts):
    # skip header and process things directly
    if beagleline0.startswith("marker"):
        continue
    if countline0.startswith("ind"):
        continue

    # process beagle line, split by tabs:
    beagleline=beagleline0.strip().split('\t')
    scaff=beagleline[0] # scaffold name
    allele1=beagleline[1] # allele 1 (reference)
    allele2=beagleline[2] # allele 2 alternate
    GT_Probs=beagleline[3:] # genotype posterior probabilities for all GTs/ inds
    # Split GTs into per-individual ; note that the second [1] of each set is the het GT
    GTs_perInd = [GT_Probs[i:i+3] for i in range(0,len(GT_Probs),3)] # groups together every set of three GTs (3 per individual) # be careful here; checked it carefully

    # make a dictionary:
    GTs_perInd_Dict = dict(zip(sampList,GTs_perInd)) # note that zip maintains the respective orders, but then dict orders alphabetically. this is okay as long as zip happens before dict 
    # iterate through the individuals
    # get max per set:
    # process count line, split by tabs:
    countline=countline0.strip().split('\t')
   # make a dictionary of the read counts per individual:
    counts_perInd_Dict = dict(zip(sampList,countline))
    for sample in sampList:
        # first check if there isn't missing data, otherwise skip it
        # careful, was treating counts as strings for some reason
        if float(counts_perInd_Dict[sample]) < float(PerIndividualDepthMinimum):
            missingDict[sample]+=1 # add to the missing count 
        # but if it's greater than 1 (non missing data),
            # first check to see if it passes missingness threshold,
            # and if it does, keep and add to het sum, if not, discard.
        elif float(counts_perInd_Dict[sample])>= float(PerIndividualDepthMinimum):
            GTs=GTs_perInd_Dict[sample]
        # first check if the maximum post prob in the GT is >= cutoff (e.g. 50%)
        # If not, then call genotype as missing
            if float(max(GTs)) >= float(MaxProbCutoff):
                hetProb=GTs[1] # the middle value is the heterozygosity posterior value
                hetProbSumDict[sample]+=float(hetProb) # add it to the total probability 
            # then add 1 to the called site dictionary for that individual:
                calledDict[sample]+=1 
            # check if it's a transversion, if yes, add to the transvHetDict, if not don't
                if (allele1,allele2) in transversions:
                    TransvOnly_HetProbSumDict[sample]+=float(hetProb)
            
            elif float(max(GTs)) < float(MaxProbCutoff):
            # if it doesn't pass the cutoff, skip it and add a 1 to the missing dict for bookkeeping purposes
                missingDict[sample]+=1

heterozygosityDict=dict()
for sample in sampList:
    heterozygosityDict[sample] = hetProbSumDict[sample]/calledDict[sample]

TransvOnly_heterozygosityDict=dict()
for sample in sampList:
    TransvOnly_heterozygosityDict[sample] = TransvOnly_HetProbSumDict[sample]/calledDict[sample]

# note that transv and regular het have same denominator (all callable sites)
## want to save this information

# missingDict, calledDict, hetProbDict, heterozygosityDict
outfile=open(outname, "w")
outheader="sample\tuncallableSites\tcallableSites\tsumHetProb\tsumHetProb_TransversionsOnly\tHetPerSite\tHetPerSite_transversionsOnly\tPerIndividualDepthMinimum\tProbThresholdForCallableSite\n"
outfile.write(outheader)
for sample in sampList:
    out=[sample,str(missingDict[sample]),str(calledDict[sample]),str(hetProbSumDict[sample]),str(TransvOnly_HetProbSumDict[sample]),str(heterozygosityDict[sample]),str(TransvOnly_heterozygosityDict[sample]),str(PerIndividualDepthMinimum),str(MaxProbCutoff)]
    outfile.write("\t".join(out))
    outfile.write("\n")
outfile.close()
beagle.close()
counts.close()
