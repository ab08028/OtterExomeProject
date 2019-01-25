# Script to make SLIM job scripts and submission scripts for wolf simulations
# To test without running the simulations, comment out the qsub line

# Set g, number of genes
g=1000

# Set t, number of burn-in generations
t=450000

# Do each dominance coefficient
for h in 0.0 0.5
do

# Write submission scripts
cat > submit_slim_wolf_${g}genes_${t}burn_h${h}_101517.job.sh << EOM
#! /bin/bash
#$ -wd /u/scratch/b/bkim331/jacqueline/wolf
#$ -l h_rt=300:00:00,h_data=12G,highp
#$ -o /u/scratch/b/bkim331/jacqueline/wolf
#$ -e /u/scratch/b/bkim331/jacqueline/wolf
#$ -N wolfsim
#$ -t 1-250:1

source /u/local/Modules/default/init/modules.sh
module load gcc/6.3.0

SLIMDIR=/u/scratch2/b/bkim331/jacqueline/SLiM/bin

/u/scratch/b/bkim331/jacqueline/wolf/make_slim_wolf_101517.job.sh ${g} ${t} ${h} \$(printf %03d \${SGE_TASK_ID})

\${SLIMDIR}/slim -d "mu=1e-8" /u/scratch/b/bkim331/jacqueline/wolf/slim_wolf_${g}genes_${t}burn_h${h}_101517.job

EOM

# Submit
qsub submit_slim_wolf_${g}genes_${t}burn_h${h}_101517.job.sh

done
