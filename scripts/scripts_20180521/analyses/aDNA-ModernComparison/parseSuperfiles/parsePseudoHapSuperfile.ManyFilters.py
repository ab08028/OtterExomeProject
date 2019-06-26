"""
Created on Fri May 10 12:14:23 2019

@author: annabelbeichman

In this script, we want to parse a angsd pseudoHaploid "superfile" that I generate that has concatenated bed-format (0based) coordinates in the first 12 columns, then the angsd maf output, then pseudohaploid calls from randomly sampling reads, and finally count data per individual/per site

We can parse that file and count the number of callable sites per individual that pass a depth filter and minimum number of individuals with data at that site (optional)

With this script we count up the number reference and non-ref calls, as well as number of missing lines (overall) and missing GTs per ind

usage: python script.py inputFilepath sampleIDFile outputFile MaxProbCutoff PerIndividualDepthMinimum minIndsPerSite
"""
import gzip
import sys
filepath = sys.argv[1] #path to input superfile file, should contain transitions and transversions (a concatenation of angsd results in bed format, then mafs, then GPs or GLs, then counts -- generated from previous script)
sampleIDFile=sys.argv[2] # path to file with list of names in SAME ORDER as bamList you used for angsd
outname= sys.argv[3] # output file
# it's very important for these to be floats
#MaxProbCutoff=float(sys.argv[4]) # # if the max of the 3 probs is below this, discard; keep if >= to the cutoff
PerIndividualDepthMinimum=float(sys.argv[4])
minIndsPerSite=float(sys.argv[5]) # min number of individuals that have data at a site (note: won't be counted unless they also pass the PerIndividualDepthMinimum; so if you require 2 inds at 2 read depth (not recommended, just an example) then if you had an individual with 1 read, and an individual with 2 reads, it wouldn't pass.)
#sys.stdout.write("helloworld")
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
transversions=[('A','C'),('C','A'),('A','T'),('T','A'),('C','G'),('G','C'),('G','T'),('T','G')]

######### for testing only (set these as arguments eventually) #######################
#filepath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/parseHaploFile/dummy.pseudohapsuperfile.0based.bed.gz"

# this was drawn from high coverage, AF prior 20190524
# list of samples   # check order super carefully!!! must be in same order as input bam list for angsd!!!!!!!
#sampleIDFile="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt"

#outname="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/parseHaploFile/testout.superfile.txt"
#MaxProbCutoff=0.95 # if the max of the 3 probs is below this, discard; keep if >= to the cutoff
#PerIndividualDepthMinimum=1
#minIndsPerSite=2


######################################################

# read samples into a list (strip \n from end of each one)
sampList = [line.rstrip('\n') for line in open(sampleIDFile)]
numInd=len(sampList)
# these are now in the proper order ### MAKE SURE IT'S ORDER YOUR BAMLIST WAS IN FOR ANGSD!!!!! OTHERWISE INDS WILL BE ASSIGNED INCORRECTLY

###### make empty dictionaries : ##############
#### counts of missing sites: 
missingDict=dict()
#### counts of called sites: 
calledDict=dict()
#### count of non-ref calls
refCallCountDict=dict()
nonRefCallCountDict=dict()
#### count of transversions non-ref calls
nonRefCallCountDict_TvOnly=dict()

# populate all the dicts with sample IDs and 0s:
for sample in sampList:
    calledDict[sample]=0
    missingDict[sample]=0
    nonRefCallCountDict[sample]=0
    nonRefCallCountDict_TvOnly[sample]=0
    refCallCountDict[sample]=0
# want to count up number of TV non-ref calls per individual from pseudohaps
# Can use this script on the overall cds file to get number of callable sites per individual and then can use it on the syn, mis, sg category files to get those counts per individual.

########### Open beagle GL posteriors file #############

superfile = gzip.open(filepath,"r")

# get header:
header=[]
for line in superfile:
    if "#" in line:
        header=line.strip().split('\t')

        break

markerIndex=header.index("markerID")
scaffIndex=header.index("chr")
posIndex=header.index("pos")
refIndex=header.index("ref")
hapIndex=header.index("ind0") # where beagle pseudohaps start
countIndex=header.index("ind0TotDepth") # where per-ind read counts start
# reset file:
superfile.seek(0)

# count lines that you skip because they are triallelic:
skippedLines=0
triallelic=0
for line0 in superfile:
    # skip header 
    if line0.startswith("#"):
        continue
    # process line, split by tabs:
    line=line0.strip().split('\t')
    # First want to check if the line passes the nInd filter
    # But want that nInd filter to also account for the minDepth filter
    # so need a couple things:
    counts=line[countIndex:]
    # need to check if a minimum number of counts pass the count threshold
    indsPassingThreshold = sum(float(i) >= float(PerIndividualDepthMinimum) for i in counts) # want to check if this value is >= the min number of individiduals at the site
    if float(indsPassingThreshold) < float(minIndsPerSite):
        skippedLines+=1
        continue # skip lines that don't pass
      # if the number of individuals with depth above the depth cutoff passes the minInd threshold, then we are off to the races (but want to record nInd passing threshold)
    # but note that if an individual GT still doesn't pass depth cutoff, it won't be included in the calculation for that individual which is good
    elif float(indsPassingThreshold) >= float(minIndsPerSite):
        marker=line[markerIndex] # marker name
        scaff=line[scaffIndex]
        pos=line[posIndex]
    # make sure that marker and chrom/pos are the same
        if marker != scaff+"_"+pos:
            print("You have concatenated the wrong files together into a superfile!!")
            break
        ref=line[refIndex]
        HapsList=line[hapIndex:(hapIndex+(numInd))] # pseudoHaploid calls
        HapsSetList=list(set(HapsList))
        HapsSetListPlusRef=list(set(HapsList+list(ref))) # add in the reference alelle to the Haps list so it's part of the set, will help with identifying triallelic sites speedily

        # want to combine this with the reference allele to make a full set of possible haps
        # GET the alleles including ref allele: (and make sure it's not triallelic)
        HapsPlusRefNoNs = [x for x in HapsSetListPlusRef if x!="N" ]

        if len(HapsPlusRefNoNs)>2: # skip if triallelic when you include the ref allele
            triallelic+=1
            continue

        # if it's monomorphic for ref (len = 1) or biallelic (len =2, one of which is ref allele because you added it in)    
        elif len(HapsPlusRefNoNs) <=2:

            Haps_perInd_Dict = dict(zip(sampList,HapsList)) # note that zip maintains the respective orders, but then dict orders alphabetically. this is okay as long as zip happens before dict 
        # iterate through the individuals

        # want to get count information:
        # make a dictionary of the read counts per individual:
            counts_perInd_Dict = dict(zip(sampList,counts))
            for sample in sampList:
                Hap=Haps_perInd_Dict[sample]
                # skip if it's uncalled 
                if Hap=="N":
                    missingDict[sample]+=1
                    continue
            # check if too few reads per individual: 
                elif float(counts_perInd_Dict[sample]) < float(PerIndividualDepthMinimum):
                    missingDict[sample]+=1 # add to the missing count 
            # check if it passes read counts per individual
                elif float(counts_perInd_Dict[sample])>=    float(PerIndividualDepthMinimum):
                # count the site as callable:
                    calledDict[sample]+=1 
                    # and if it's a ref allele count it
                    if Hap == ref:
                        refCallCountDict[sample]+=1
                    # if the alle is not an 
                    elif Hap != ref:
                        nonRefCallCountDict[sample]+=1 # no matter what
                        # and if it's a transversion add it in

                        if tuple(HapsPlusRefNoNs) in transversions:
                            nonRefCallCountDict_TvOnly[sample]+=1
# how I checked this script (have to be very careful, while writing it I had some logic errors that are fixed now.)
# did checks in bash based on awk for ind0:
# count up sites where ind0 has a non-N GT and it is not equal to the ref allele
# zcat dummy.pseudohapsuperfile.0based.bed.gz | awk '{if($23!="N" && $17!=$23)print}' | wc -l --> this equaled 2672
# and then if you subtract the 3 triallelic sites you get 2669
# which is what this python script calculates. So it's correct now! hooray

# glad I checked, I wasn't adding the Ti+TV together but now I am




# missingDict, calledDict, hetProbDict, heterozygosityDict
outfile=open(outname, "w")
# also want to calculate total hom-alt GTs (GPs or GLs -- not totally sure. try it both ways.)
outheader="sample\tuncallableSites\tcallableSites\trefCalls\tnonRefCalls_all\tnonRefCalls_TransversionsOnly\tFilter_PerIndividualDepthMinimum\tFilter_minIndsPerSite\tSkippedLinesOverallNotPerInd\tTriallelicOverallNotPerInd\n"
outfile.write(outheader)
for sample in sampList:
    out=[sample,str(missingDict[sample]),str(calledDict[sample]),str(refCallCountDict[sample]),str(nonRefCallCountDict[sample]),str(nonRefCallCountDict_TvOnly[sample]),str(PerIndividualDepthMinimum),str(minIndsPerSite),str(skippedLines),str(triallelic)]
    outfile.write("\t".join(out))
    outfile.write("\n")
outfile.close()
superfile.close()
