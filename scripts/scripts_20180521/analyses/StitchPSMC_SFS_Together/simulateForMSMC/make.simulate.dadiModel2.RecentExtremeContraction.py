# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 13:40:10 2017

@author: annabelbeichman
"""
# This is to set up a MACS script to simulate the full demographic history of the sea otter or giant otter
# want to make generic
# input
# spp
# mu
# Len

import pandas as pd

mu= float(8.64e-09)


r = float(1e-08)
Len = 30000000
num=60
blocksPerGroup=10
groups=num/blocksPerGroup

Na = float(3500) # set based on dadi model
nu = float(0.06) # set at 195/3500
T_macs = float(0.0025) # set at 35/(2*3500) /2 to get it into macs format (scaled by 4N not 2N)

ss=2 # two haplotypes (one genome)
theta = 4*Na*mu
rho = 4 * Na*r

# Note that theta will be the same across different values of Mu for the same
# models, because mutation rate changes scaling of MSMC and so alters Na, but then
# theta = 4Namu = 4 (1/(lamba*2mu))*mu = 2/lamba [is this right?]
############################## WRITE SCRIPT ################
print("#!/bin/bash")
print("#$ -cwd")
print("#$ -l h_rt=5:00:00,h_data=28G,highp")
print("#$ -N simDadiModel2")
print("#$ -m abe")
print("#$ -M ab08028")

print("module load python/3.7") # >3.2 to have argparse
print("#may need to pip install argparse for ms2multihetsep.py")
print("rundate=`date +%Y%m%d`")
print("replicate=$SGE_TASK_ID")

print("model=dadiModel2.RecentExtremeContraction")
print("wd=/u/flashscratch/a/ab08028/captures/analyses/simulateForMSMC")
print("cd $wd")
print("macsFile=/u/home/a/ab08028/klohmueldata/annabel_data/bin/macs")
print("msformatterFile=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msformatter")
print("ms2multiFile=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msmc-tools/ms2multihetsep.py")
#print("cp $macsFile $msformatter $ms2multiFile $wd")
print("mkdir -p ${model}")
print("for j in {1.."+str(groups)+"}")
# simulate slightly more than you need
print("do")
print("mkdir -p ${model}/rep_${replicate}/group_$j.${model}")

print("cd ${model}/rep_${replicate}/group_$j.${model}")
print("cp -n $macsFile ./")
print("cp -n $msformatterFile ./")
print("cp -n $ms2multiFile ./")
print("for i in {1.."+str(blocksPerGroup)+"}")
print("do")

print("# dadi model 2 for msmc")
print("mu="+str(mu))
print("r="+str(r))
print("Na="+str(Na))
print("rho=" +str(rho))
print("theta="+str(theta))
print("date=`date +%Y%m%d`")
print("SEED=$((date+$RANDOM+((j-1)*"+str(blocksPerGroup)+")+i))") #
print("# this is a new addition! need to have a different random seed for each simulation; if they start within a second of each other, they will have the same seed. not an issue for big simulations of 30Mb because those are slow, but 100kb can start within a second of each other!")
print("./macs " +str(ss) +" "+str(Len)+" -t "+str(theta)+" -r "+str(rho)+" -s $SEED"+" -eN 0.0 "+str(nu)+" -eN "+str(T_macs)+" 1"),
#for x, y in zip(times_gen_trimancient_4Na,diploids_trimancient_Na):
#    print("-eN " + str(x)+" "+str(y)),
print(" > group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt")

print("#convert to ms format")
print("./msformatter < group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt > group_${j}_block_${i}.${model}.msFormat.OutputFile.${rundate}.txt")
print("#convert to msmc input format")
print("python3 ms2multihetsep.py $i "+ str(Len) +" < group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt > group_${j}_block_${i}.${model}.MSMCFormat.OutputFile.${rundate}.txt")

###################################################
print("done")
print("cd $wd")
print("done")
