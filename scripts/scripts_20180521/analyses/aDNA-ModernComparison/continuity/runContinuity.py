# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 15:59:52 2019

@author: annabelbeichman
"""
import sys
import numpy as np
import scipy
# path to continuity:
#sys.path.append('/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/scripts_from_others/continuity/')
#ancient_genotypes = imp.load_source('ancient_genotypes', '/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/generateContinuityInput/continuity/ancient_genotypes.py')


pathToContinuity= sys.argv[1]
readsFilePath = sys.argv[2] #path to input file in continuity format
indsFilePath=sys.argv[3] # path to file specifying group membership of ancient samples
modernReference = sys.argv[4]# name of modern reference 
outfilepath= sys.argv[5] # output file
sys.path.append(pathToContinuity) # add continuity
from ancient_genotypes import *

#readsFilePath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/continuity/continuityInput/CA.angsdOut.mappedTomfur.ancient.counts.freqsFromModernGATK.superfile.1based.contInput.TVOnly.txt"
#indsFilePath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/generateContinuityInput/ancInds.txt"
print("******** reading in data **********")
unique_pops, inds, label, pops, freqs, read_lists = parse_reads_by_pop(readsFilePath,indsFilePath,cutoff=0)

####### could filter coverage here ##########
####### may want to do alpha and beta potentiall, but would need allele counts but can account for variable coverage #####

# Okay so I think this is sort of working ; the freqs are unique frequencies (representing multiple sites?), and then they are part of a dictionary with the corresponding ancient reads for each of those sites -- I think? confused here.

# from continuity docs: "For various reasons, such as cryptic structural variation, you may have some sites with unusually high/low coverage in your ancient samples. These sites can severely mislead analyses. A common approach is to remove sites that fall in the tails of the coverage distribution. You can do that once you've read in your data by using the function coverage_filter(). It will set any sites that violate a coverage filter in an individual to have 0 coverage. It takes 3 arguments:"

    #read_lists: The read_lists object that comes out of parse_reads_by_pop()
    #min_cutoff: sites with coverage less than this percentile will be removed (default 2.5)
    #max_cutoff: sites with coverage more than this percentile will be removed (default 97.5)
# So cut off the bottom and top 2.5% percentiles -- sure, I guess. Might cut off too many sites though. Try with and without
# This function operates on read_lists IN PLACE, meaning that read_lists will be modified. The function returns a single value, which is a list of the coverage cutoffs for each individual in each population.
print("******** filtering ancient coverage -- removing < 2.5% percentile, >97.5% percentile. The percentiles are: **********")
coverage_filter(read_lists,2.5,97.5)

######## optimize #######
print("********* beginning optimization ***********")
num_core=2
opts_cont_false=optimize_pop_params_error_parallel(freqs,read_lists,num_core,detail=True,continuity=False)
opts_cont_true=optimize_pop_params_error_parallel(freqs,read_lists,num_core,detail=True,continuity=True)
### a simple demographic model: from cont paper: "To calculate the sampling probability, we assume a simple demographic model, in which the ancient individual belongs to a population that split off from the modern population τ1 generations ago, and subsequently existed as an isolated population for τ2 generations."
# so opts_cont_false[0][0] will be for population 0 and give t1, t2 and error for ind1, ind2, ...

# from continuity docs:
likelihood_false = np.array([-x[1] for x in opts_cont_false]) #minus sign is because scipy.optimize minimizes the negative log likelihood
likelihood_true = np.array([-x[1] for x in opts_cont_true])
LRT = 2*(likelihood_false - likelihood_true)
log_p_vals = scipy.stats.chi2.logsf(LRT,1) #returns the LOG p-values

############ write out ################
# lhood false, true, log_p_vals, params
outfile=open(outfilepath,"w")
outfile.write("ModernRefPopulation\tancientPopulation\tlhood_false\tlhood_true\tLRT\tlog_p_vals\tparams_false\tparams_true\n")
# for each ancient population go through and write out a line:
for i in range(len(unique_pops)):
    pop=unique_pops[i]
    lhoodFalse_pop=likelihood_false[i]
    lhoodTrue_pop=likelihood_true[i]
    LRT_pop=LRT[i]
    pvalPop=log_p_vals[i]
    paramsFalse_pop=opts_cont_false[i][0]
    paramsTrue_pop=opts_cont_true[i][0]
    out="\t".join([modernReference,pop,str(lhoodFalse_pop),str(lhoodTrue_pop),str(LRT_pop),str(pvalPop),str(paramsFalse_pop),str(paramsTrue_pop)])
    outfile.write(out)
    outfile.write("\n")

outfile.close()