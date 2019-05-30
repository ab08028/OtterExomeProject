"""
Created on Fri May 10 12:14:23 2019

@author: annabelbeichman

In this script, we want to parse a angsd "superfile" that I generate that has concatenated bed-format (0based) coordinates in the first 12 columns, then the angsd maf output, then genotype probabilities or genotype lhoods (can be either), and finally count data per individual/per site

We can parse that file and count the number of callable sites per individual that pass posterior prob threshold, depth filter, and minimum number of individuals with data at that site (affects prior)

The script will calculate heterozygosity for all sites and just for transversions 

It sums up heterozygosity posterior probability for each individual and divides by callable GTs (following methods of Fages et al. 2019 in Cell)

Explantion: the script sets up a dictionary for your list of individuals (order of list MUST match the input bam file list in ANGSD)

It then goes through the file site by site and for each individual it checks if the maximum posterior for the individual's three genotypes is >= some threshold (e.g. 0.5 or 0.95) .

If it is, it adds the heterozygosity posterior probabilty into that individual's dictionary entry (heterozygosity numerator), and tallies the GT as a callable Site for that individual (heterozygosity denominator)
# It will also calculate heterozygosity for transversions only
# and then reports the total sites passing the threshold per individual, the total het probabilities, and divides the two (hets/total sites) for all hets or just transversions


usage: python script.py inputFilepath sampleIDFile outputFile MaxProbCutoff PerIndividualDepthMinimum minIndsPerSite
"""
import gzip
import sys

filepath = sys.argv[1] #path to input superfile file, should contain transitions and transversions (a concatenation of angsd results in bed format, then mafs, then GPs or GLs, then counts -- generated from previous script)
print(filepath)
sampleIDFile=sys.argv[2] # path to file with list of names in SAME ORDER as bamList you used for angsd
outname= sys.argv[3] # output file
MaxProbCutoff=sys.argv[4] # # if the max of the 3 probs is below this, discard; keep if >= to the cutoff
PerIndividualDepthMinimum=sys.argv[5]
minIndsPerSite=sys.argv[6] # min number of individuals that have data at a site (note: won't be counted unless they also pass the PerIndividualDepthMinimum; so if you require 2 inds at 2 read depth (not recommended, just an example) then if you had an individual with 1 read, and an individual with 2 reads, it wouldn't pass.)

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
#filepath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/parseBeagleSuperfile/test.superfileForPythonTests.bed.gz"
# this was drawn from high coverage, AF prior 20190524
# list of samples   # check order super carefully!!! must be in same order as input bam list for angsd!!!!!!!
#sampleIDFile="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt"

#outname="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/parseBeagleSuperfile/testout.txt"
#MaxProbCutoff=0.95 # if the max of the 3 probs is below this, discard; keep if >= to the cutoff
#PerIndividualDepthMinimum=1
#minIndsPerSite=2


######################################################

# read samples into a list (strip \n from end of each one)
sampList = [line.rstrip('\n') for line in open(sampleIDFile)]
numInd=len(sampList)
# these are now in the proper order ### MAKE SURE IT'S ORDER YOUR BAMLIST WAS IN FOR ANGSD!!!!! OTHERWISE INDS WILL BE ASSIGNED INCORRECTLY

###### make empty dictionaries : ##############
#### counts of missing GTs: 
missingDict=dict()
#### counts of called GTs: 
calledDict=dict()
#### sums of het/homAlt/homRef GPs or GLs passing filters
hetProbSumDict=dict()
homAltProbSumDict=dict()
homRefProbSumDict=dict()
####  transversions:
TransvOnly_HetProbSumDict=dict()
TransvOnly_HomAltProbSumDict=dict()
TransvOnly_HomRefProbSumDict=dict()
# populate all the dicts with sample IDs and 0s:
for sample in sampList:
    calledDict[sample]=0
    missingDict[sample]=0
    hetProbSumDict[sample]=0
    homAltProbSumDict[sample]=0
    homRefProbSumDict[sample]=0
    TransvOnly_HetProbSumDict[sample]=0
    TransvOnly_HomAltProbSumDict[sample]=0
    #TransvOnly_HomRefProbSumDict[sample]=0 # don't need HomRef transversions -- they are HomRef so aren't elevated by transversions

########### Open beagle GL posteriors file #############

superfile = gzip.open(filepath,"rb")
#counts = gzip.open(countsFile,"rb") # don't need counts file because it's all in super file

# get beagle header:
header=[]
for line in superfile:
    if "#" in line:
        header=line.strip().split('\t')
        # check that length of header is numInd*3 + 3 (marker allele1 allele2):
        #len(header)==(numInd*3) + 3 # should be TRUE
        break
# need to figure out positions within the header of certain things
# know header 0-11 will be bed fmt stuff
# and 12-18 (?) will be mafs
# note that .index("foo") will only return teh FIRST entry and will throw and error if it isn't there. in this case that is okay, but just be aware.
# things I need for filtering: nInd; 
#perSitenIndIndex=header.index("nInd") # so not really using this information from MAFs; actually just getting it from the "counts" part of the file. (counting number of inds that pass depth filter so both filters are applied together to the whole line; then again perGT am applying the depth filter)
allele1Index=header.index("allele1")
allele2Index=header.index("allele2")
beagleIndex=header.index("Ind0") # where beagle GPs or GLs start
countIndex=header.index("ind0TotDepth") # where per-ind read counts start
markerIndex=header.index("marker")
# reset file:
superfile.seek(0)


for line0 in superfile:
    # skip header and process things directly
    if line0.startswith("#"):
        continue
    # process beagle line, split by tabs:
    line=line0.strip().split('\t')
    # First want to check if the line passes the nInd filter
    # But want that nInd filter to also account for the minDepth filter
    # so need a couple things:
    counts=line[countIndex:]
    # need to check if a minimum number of counts pass the count threshold
    indsPassingThreshold = sum(float(i) >= float(PerIndividualDepthMinimum) for i in counts) # want to check if this value is >= the min number of individiduals at the site
    if float(indsPassingThreshold) < float(minIndsPerSite):
        continue
    # dont need to add this line to failed GTs dict, because it fails for all individuals. just don't include it in counts at all.
    # if the number of individuals with depth above the depth cutoff passes the minInd threshold, then we are off to the races (but want to record nInd passing threshold)
    # but note that if an individual GT still doesn't pass depth cutoff, it won't be included in the calculation for that individual which is good
    elif float(indsPassingThreshold) >= float(minIndsPerSite):
        marker=line[markerIndex] # marker name
        allele1=line[allele1Index] # allele 1 (reference)
        allele2=line[allele2Index] # allele 2 alternate
        # GT_Probs (or GLs) are located from the first "Ind0" until that index + 2*Number of individuals because there are 3 GTs per individual; note python range is non-inclusive of last number so adding these together works as 22:49, where 48 is the last one you want. I checked that this works and it does.
        GT_Probs=line[beagleIndex:(beagleIndex+(numInd*3))] # genotype posterior probabilities or Lhoods for all GTs/ inds
# Split GTs into per-individual ; note that the second [1] of each set is the het GT
        GTs_perInd = [GT_Probs[i:i+3] for i in range(0,len(GT_Probs),3)] # groups together every set of three GTs (3 per individual) # be careful here; checked it carefully
# make a dictionary:
        GTs_perInd_Dict = dict(zip(sampList,GTs_perInd)) # note that zip maintains the respective orders, but then dict orders alphabetically. this is okay as long as zip happens before dict 
        # iterate through the individuals
        # get max per set:
        # want to get count information:
        # make a dictionary of the read counts per individual:
        counts_perInd_Dict = dict(zip(sampList,counts))
        for sample in sampList:
            # first check if there isn't missing data, otherwise skip it
            # careful, was treating counts as strings for some reason
            # if either of these things are true ()
            GTs=GTs_perInd_Dict[sample]
            if float(counts_perInd_Dict[sample]) < float(PerIndividualDepthMinimum) or float(max(GTs)) < float(MaxProbCutoff):
                missingDict[sample]+=1 # add to the missing count 
                # check if it passes both filters (read counts per individual and the maxProbCutoff); must pass both to keep add to het sum
            elif float(counts_perInd_Dict[sample])>=    float(PerIndividualDepthMinimum) and float(max(GTs)) >=    float(MaxProbCutoff):
                # count the site as callable:
                calledDict[sample]+=1 
                # get het/homAlt/homRef GPs or GLs:
                # don't need to define these as variables, just plunk them in (make sure indices are right); should be homRef 0, het 1, homAlt 2
                #homRefProb=GTs[0] # the first value is homRef
                #hetProb=GTs[1] # the middle value is the heterozygosity posterior value
                #homAltProb=GTs[2] # the third value is homALt
                homRefProbSumDict[sample]+=float(GTs[0]) # add it to the total probability 
                hetProbSumDict[sample]+=float(GTs[1]) # add it to the total probability 
                homAltProbSumDict[sample]+=float(GTs[2]) # add it to the total probability 
    # check if it's a transversion, if yes, add to the transvHetDict, if not don't
                if (allele1,allele2) in transversions:
                    # note transversions only for HomRef doesn't make sense, so not tracking.
                    #TransvOnly_HomRefProbSumDict[sample]+=float(GTs[0])
                    TransvOnly_HetProbSumDict[sample]+=float(GTs[1])
                    TransvOnly_HomAltProbSumDict[sample]+=float(GTs[2])


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
# also want to calculate total hom-alt GTs (GPs or GLs -- not totally sure. try it both ways.)
outheader="sample\tuncallableSites\tcallableSites\tsumHetGLsOrGPs\tsumHetGLsOrGPs_TransversionsOnly\tsumHomAltGLsOrGPs\tsumHomAltGLsOrGPs_TransversionsOnly\tsumHomRefGLsOrGPs\tHetPerSite\tHetPerSite_TransversionsOnly\tFilter_PerIndividualDepthMinimum\tFilter_minIndsPerSite\tFilter_ProbThresholdForCallableSite\n"
outfile.write(outheader)
for sample in sampList:
    out=[sample,str(missingDict[sample]),str(calledDict[sample]),str(hetProbSumDict[sample]),str(TransvOnly_HetProbSumDict[sample]),str(homAltProbSumDict[sample]),str(TransvOnly_HomAltProbSumDict[sample]),str(homRefProbSumDict[sample]),str(heterozygosityDict[sample]),str(TransvOnly_heterozygosityDict[sample]),str(PerIndividualDepthMinimum),str(minIndsPerSite),str(MaxProbCutoff)]
    outfile.write("\t".join(out))
    outfile.write("\n")
outfile.close()
superfile.close()
