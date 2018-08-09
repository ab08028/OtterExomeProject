# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 10:49:04 2018

@author: annabelbeichman
"""
######## this is a sandbox to play with pyvcf ################

import sys
import gzip
import re
import time
import vcf # this is pyvcf

filepath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/dummyVCF.forSandbox.allSites_5_passingFilters.vcf.gz" # eventually want this to just be snps
outfname="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/6_allSites_bespokeFilters_passingFilters_80percCall_elut.raw_variants.20170914.vcf"
errorname="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/sitesFailing_bespokeFilters.vcf"
inVCF = gzip.open(filepath, 'r')

outVCF = open(outfname, 'w')
errorVCF = open(errorname, 'w')

# https://github.com/jamescasbon/PyVCF
# reinitialize each time
vcf_reader = vcf.Reader(open(filepath, 'r'))
vcf_writer = vcf.Writer(open(outfname, 'w'), vcf_reader)
#record=next(vcf_reader)
# this is very handy
# Practice tasks:
    # 1. eliminate lines with >1bp ref allele
    # 2. 
test=[]

def findAllHetLines(vcf_reader):
    test=[]
    for record in vcf_reader:
        if record.num_het==record.num_called:
            errorVCF.write("This line is all heterozygotes / no-calls)
            errorVCF.write(record)
        else:
            outVCF.write(record) 
            return record
test = findAllHetLines(vcf_reader)


# idea: get het count
# get all count count 
def countNoCallGTs(vcf_reader):
    samples=vcf_reader.samples
    for record in vcf_reader:
        for sample in samples:
            if record.genotype(sample)['GT']=="./."
            return record.genotype
            #noCallPerSample[samples[1]]
            
        return record.num_called
    
### useful items:
samples=vcf_reader.samples # get list of samples *** keep **** 
test=[]
# first test: get records where all are hets:
allHets=[]
while vcf_reader:
    record=next(vcf_reader)
    print record
    if record.num_het == record.num_called:
        allHets.append(record)
    else:
        continue
    
print record.num_called, record.call_rate, record.num_unknown
# 26 1.0 0 # this is useful; gives how many are called, and the call rate (100%), and # unknown
# can use this to get an/ac maybe?

# print number of hom ref, hets and num_hom_alt (so can do if num_hom alt is num_called?) -- FAIL het. 
#        Record.FORMAT
#        Record.samples
 #       Record.genotype
print record.num_hom_ref, record.num_het, record.num_hom_alt
# 26 0 0 
print record.get_hets()
print record.num_called, record.call_rate, record.num_unknown
print record.num_hom_ref, record.num_het, record.num_hom_alt
print record.nucl_diversity, record.aaf, record.heterozygosity
print record.get_hets()
#[Call(sample=NA00002, CallData(GT=1|0, GQ=48, DP=8, HQ=[51, 51]))]
print record.is_snp, record.is_indel, record.is_transition, record.is_deletion
print record.var_type, record.var_subtype
print record.is_monomorphic
print record.samples # nice way to explore
# then you can look at each entry
for sample in record.samples:
    print sample['RGQ']
    

