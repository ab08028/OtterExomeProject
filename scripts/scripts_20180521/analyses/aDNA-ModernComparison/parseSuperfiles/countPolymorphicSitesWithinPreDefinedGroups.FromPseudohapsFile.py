# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 13:35:27 2019

@author: annabelbeichman
This script is not optimized for other projects -- is exploratory
Idea: to go through the pseudohaploid file
and for each set of three inds (CA, AK, ancient)
want to count the total number of sites where at least 2/3 are called
and then of those, how many are polymorphic
(see if ancient is more polymorphic)

# this script requires that you pre-set the group membership of the individuals, though the order of individuals in the input file comes from a different file -- this isn't ideal for using in other projects, but works for me for now # 
"""
import gzip
import sys
# usage: script [input] [output] [sampleIDFile] [maxMissingInds]
filepath = sys.argv[1]  #path to input haplofile
outname= sys.argv[2] # output file
sampleIDFile=sys.argv[3] # path to file with list of names in SAME ORDER as bamList you used for angsd
maxMissingInd = sys.argv[4] # max missing inds per population for a site to be considered 'called' for that population ; count of "Ns" within a popualtion must be <= this number. 
# I think order of individuals is: (TEMP)
#maxMissingInd= 1 
sampList = [line.rstrip('\n') for line in open(sampleIDFile)]
#indOrder=["116_Elut_CA_307","126_Elut_AK_AF3394","129_Elut_AK_AL4660","140_Elut_CA_403","141_Elut_CA_419","55_Elut_AK_AF3736","A13_Elut_CA_AN_388_SN1_2CAP_screen","A29_Elut_CA_SM_30_SN2_CAP","A30_Elut_CA_SM_35_SN1_CAP"]

CAInds=["116_Elut_CA_307","140_Elut_CA_403","141_Elut_CA_419"]
AKInds=["126_Elut_AK_AF3394","129_Elut_AK_AL4660","55_Elut_AK_AF3736"]
ancInds=["A13_Elut_CA_AN_388_SN1_2CAP_screen","A29_Elut_CA_SM_30_SN2_CAP","A30_Elut_CA_SM_35_SN1_CAP"]

# populate a dict with individuals per pop:

indDict=dict()
indDict["CA"] = CAInds
indDict["AK"] = AKInds
indDict["ancient"] = ancInds

# populate empty dicts to keep track of monomorphic and polymorphic sites
popList=["CA","AK","ancient"]
monoDictAll=dict() # for all sites
polyDictAll=dict() # for all sites
polyDictTv=dict() # for transversions only
sitesCalled=dict() # keep track of sites that are called in 2/3 inds
for pop in popList:
    monoDictAll[pop]=0 
    polyDictAll[pop]=0
    sitesCalled[pop]=0
    polyDictTv[pop]=0
# for testing:
#filepath='/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/parseHaploFile/dummy.haplo.gz'
#outname='/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/testout.polymorphic.txt'
#badOut='/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/badout.txt'
haploFile =gzip.open(filepath,"r")
outfile=open(outname,"w")
#badFile=open(badOut,"w")
############### TRANSVERSIONS #############  so the 8 possible transversions are:
            # 0-1 : A-C  
            # 1-0 : C-A 
            # 0-3 : A-T 
            # 3-0: T-A 
            # 1-2 : C-G 
            # 2-1 : G-C 
            # 2-3 : G-T 
            # 3-2: T-G 
transversions=[('A','C'),('C','A'),('A','T'),('T','A'),('C','G'),('G','C'),('G','T'),('T','G')]
haploFile.seek(0)
for line0 in haploFile:
    # write out header:
    if 'major' in line0:
        outfile.write(line0)
        continue
    line=line0.strip().split('\t')
    scaff=line[0]
    pos=line[1]
    major=line[2]
    indGTs=line[3:]
    indGTdict = dict(zip(sampList,indGTs)) # this only works because order is the same between sampList and indGTs; so you have to figure out how to make sure that's true! 
    for pop in popList:
        popGTs = [indGTdict[key] for key in indDict[pop]] # get GTs for the individuals that are in the pop
        # check if it's polymorphic within a population (?) or fixed in others?
        # must be at least 2 non-N entries
        #if popGTs.count("N") < (len(popGTs) -1): # if there are more Ns than 
        if popGTs.count("N") <= maxMissingInd:
            sitesCalled[pop]+=1
            # and figure out if it's polymorphic:
            GTSetList=list(set(popGTs))
            GTSetListNoNs = [x for x in GTSetList if x!="N" ]
            numAlleles=len(GTSetListNoNs) # if numALleles =1 then it's monomorphic for these 3 individuals
            # if numAlleles > 2 they are disregarded
            if numAlleles ==1:
                monoDictAll[pop]+=1
            elif numAlleles ==2:
                if tuple(GTSetListNoNs) in transversions:
                    # add to both dicts:
                    polyDictTv[pop]+=1
                    polyDictAll[pop]+=1
                else:
                    polyDictAll[pop]+=1 # just add to all sites dict
outfile=open(outname, "w")

outheader="pop\tsitesCalled\tPolymorphicSitesAll\tPolymorphicSitesTV\tMaxMissingInd\n"
outfile.write(outheader)
for pop in popList:
    out=[pop,str(sitesCalled[pop]),str(polyDictAll[pop]),str(polyDictTv[pop]),str(maxMissingInd)]
    outfile.write("\t".join(out))
    outfile.write("\n")
outfile.close()
haploFile.close()
