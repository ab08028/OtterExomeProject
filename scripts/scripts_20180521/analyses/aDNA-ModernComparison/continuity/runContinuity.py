# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 15:59:52 2019

@author: annabelbeichman
"""
import sys
import numpy as np
import scipy
# path to continuity:
sys.path.append('/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/generateContinuityInput/continuity/')
#ancient_genotypes = imp.load_source('ancient_genotypes', '/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/generateContinuityInput/continuity/ancient_genotypes.py')


from ancient_genotypes import *

readsFilePath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/generateContinuityInput/test.out.txt"
indsFilePath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/generateContinuityInput/ancInds.txt"
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
num_core=2
opts_cont_false=optimize_pop_params_error_parallel(freqs,read_lists,num_core,detail=True,continuity=False)
opts_cont_true=optimize_pop_params_error_parallel(freqs,read_lists,num_core,detail=True,continuity=True)

# from continuity docs:
likelihood_false = np.array([-x[1] for x in opts_cont_false]) #minus sign is because scipy.optimize minimizes the negative log likelihood
likelihood_true = np.array([-x[1] for x in opts_cont_true])
LRT = 2*(likelihood_false - likelihood_true)
log_p_vals = scipy.stats.chi2.logsf(LRT,1) #returns the LOG p-values