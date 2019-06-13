# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 11:14:57 2019

@author: annabelbeichman

In this script, we want to get rid of transitions (big deal) and multiallelic sites (not super important, but good to do) by filtering the .haplo. file from angsd (-doHaploCall 1 )

You should output a file that has has transitions and multiallelic sites removed

Want to filter:
1a) monomorphic sites (uninformative for PCA); these are kept in my GL files because I'm interested in fixation of derived alleles, but for PCA we remove monomorphic anyway and they just make the file big here. 
1b) multi allelic sites (note that from this kind of file you dont' know what the reference allele is. so it's possible that if CC and GG are observed that the reference is T. Not worrying about this small category of sites for now, because at least it's biallelic within the otter samples, but be aware -- Fages doesn't deal with these at all, so I'm not *too* worried about it)
2) remove transitions (note in the same way as #1, since you don't have the reference allele, you are basing transitions as what alleles are present in your sample, not relative to the reference. this is slightly different from what I do with GLs which is relative to reference. Should affect only a very small portion of triallelic sites (which many studies don't even account for, eg Fages))
Result: transversions only (no monomorphic, multiallelic)

## NOTE: I am *not* filtering on allele frequency here (doing that as part of snpRelate PCA)  ##


USAGE: python <script.py> [input haplo file ] [path to output file]
"""
import gzip
import sys

filepath = sys.argv[1] #path to input superfile file, should contain transitions and transversions (a concatenation of angsd results in bed format, then mafs, then GPs or GLs, then counts -- generated from previous script)
outname= sys.argv[2] # output file


# for testing:
#filepath='/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/parseHaploFile/dummy.haplo.gz'
#outname='/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/testout.txt'
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

for line0 in haploFile:
    line=line0.strip().split('\t')
    scaff=line[0]
    pos=line[1]
    major=line[2]
    indGTs=line[3:]
    # first want to exclude any multiallelic site (don't count Ns)
    # make GTs a set; a set will only list each unique entry once so it's a fast way to condense (like unique)
    GTSetList=list(set(indGTs)) # then can use this to test for things
    GTSetListNoNs = [x for x in GTSetList if x!="N" ]
    # a few checks
    # first if the length of the set after Ns removed is 1, then it's monomorphic and you can skip (or write it out if you uncomment the code) and move on
    # then check if it's > 2, in which case you want to skip it because it's triallelic+
    # finally if it == 2 whether it's in transversions or not
    numAlleles=len(GTSetListNoNs)
    ###### don't want to write out monomorphic #############
   # if numAlleles == 1:
    #    outfile.write(line0)
    if numAlleles != 2: # so if it's < 2 (monomorph) or > 2 (multiallelic)
        #badFile.write(line0)
        continue
    elif numAlleles == 2:
        # check if it's a transversion
        if tuple(GTSetListNoNs) in transversions:
            outfile.write(line0)
        else:
            #badFile.write(line0)
            continue

haploFile.close()
outfile.close()
#badFile.close()