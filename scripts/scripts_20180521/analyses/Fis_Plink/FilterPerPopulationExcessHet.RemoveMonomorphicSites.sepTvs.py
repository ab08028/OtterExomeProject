# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 09:43:11 2019

@author: annabelbeichman

Goal: to filter per-population vcf file to avoid sites with too many hets which may be biasing Fis calculations. 
"""
import gzip
######### want to filter excess het sites out #############
### if going to use, need to 
VCF="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/troubleshooting/troubleshootFis/CA.only.recode.vcf.gz"
inVCF = gzip.open(VCF, 'r')
inVCF.seek(0)
maxHetFilter=0.70
outfileTv=open("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/troubleshooting/troubleshootFis/testOut.maxHetFilterPerPop.0.70.FilterOutMonomorphic.TVOnly.vcf","w")
outfileTi=open("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/troubleshooting/troubleshootFis/testOut.maxHetFilterPerPop.0.70.FilterOutMonomorphic.TiOnly.vcf","w")
## want to only output lines with < 75% het
hetFailSiteCounter=0 
monoSitesCount=0
transversions=[('A','C'),('C','A'),('A','T'),('T','A'),('C','G'),('G','C'),('G','T'),('T','G')]

for line0 in inVCF:
    if line0.startswith('#'):
        outfileTi.write(line0)
        outfileTv.write(line0)
        continue
    else:
        # check if it's a Ti or Tv
        ### For all other non-header lines, set up a dicionary of genotype calls
        line=line0.strip().split('\t') # this splits line by tabs
        alleles=tuple([line[3],line[4]]) # get ref and alt alleles
        if alleles in transversions:
    
            #CHROM0    POS1    ID2    REF3    ALT4    QUAL5    FILTER    INFO    FORMAT    [indivudals]
            mygenoinfo=line[9:]
            allCalls=[i.split(":")[0] for i in mygenoinfo] # get genotype calls
            # want to know how many of those are not missing dataa;
            called_gts=sum([x!="./." for x in allCalls])
            het_count = sum([x == "1/0" or x == "0/1" or x == "1|0" or x == "0|1" for x in allCalls])
            # skip sites that have no derived alleles:
            alt_count = sum([x == "1/1" or x == "1|1" for x in allCalls]) * 2
            ref_count = sum([x == "0/0" or x == "0|0" for x in allCalls]) * 2
            #  need to add het count in!
            alt_count = alt_count + het_count
            ref_count = ref_count + het_count
            # if monomorphic 0/0 or 1/1 across all inds don't use it 
            # must be some alts and refs present (variable site)
            if alt_count==0 or ref_count==0:
                print("found monomorphic site (all 0/0 or all 1/1s) in population")
                monoSitesCount=monoSitesCount+1
            # if too many hets, don't write out the line
            elif het_count !=0 and het_count >= called_gts*float(maxHetFilter):
                    print("found a site with >="+str(float(maxHetFilter)*100)+"% of all calls hets. het count = "+str(het_count)+" genotypes: "+ str(allCalls))
                    hetFailSiteCounter = hetFailSiteCounter+1
            else:
                outfileTv.write(line0)
        elif alleles not in transversions:
                        #CHROM0    POS1    ID2    REF3    ALT4    QUAL5    FILTER    INFO    FORMAT    [indivudals]
            mygenoinfo=line[9:]
            allCalls=[i.split(":")[0] for i in mygenoinfo] # get genotype calls
            # want to know how many of those are not missing dataa;
            called_gts=sum([x!="./." for x in allCalls])
            het_count = sum([x == "1/0" or x == "0/1" or x == "1|0" or x == "0|1" for x in allCalls])
            # skip sites that have no derived alleles:
            alt_count = sum([x == "1/1" or x == "1|1" for x in allCalls]) * 2
            ref_count = sum([x == "0/0" or x == "0|0" for x in allCalls]) * 2
            #  need to add het count in!
            alt_count = alt_count + het_count
            ref_count = ref_count + het_count
            # if monomorphic 0/0 or 1/1 across all inds don't use it 
            # must be some alts and refs present (variable site)
            if alt_count==0 or ref_count==0:
                print("found monomorphic site (all 0/0 or all 1/1s) in population")
                monoSitesCount=monoSitesCount+1
            # if too many hets, don't write out the line
            elif het_count !=0 and het_count >= called_gts*float(maxHetFilter):
                    print("found a site with >="+str(float(maxHetFilter)*100)+"% of all calls hets. het count = "+str(het_count)+" genotypes: "+ str(allCalls))
                    hetFailSiteCounter = hetFailSiteCounter+1
            else:
                outfileTi.write(line0)
            
            
outfileTi.close()
outfileTv.close()
inVCF.close()