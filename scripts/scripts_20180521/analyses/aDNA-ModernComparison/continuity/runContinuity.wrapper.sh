###### wrapper to run continuity: can be run in the shell on home computer ########

#source /u/local/Modules/default/init/modules.sh
#module load python

#gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
gitDir=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/
pathToContinuity=/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/scripts_from_others/continuity
script=$scriptDir/analyses/aDNA-ModernComparison/continuity/runContinuity.py
ref=mfur # only mfur for continuity
#indir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/continuity/combinedCounts-Freqs # on hoffman
indir=$gitDir/results/analysisResults/aDNA-ModernComparison/continuity/continuityInput
IDdir=$gitDir/information/samples/continuityIDFiles # files that give population of 'ancient' (or modern) angsd samples, labeled ind0 ind1 ind2
outdir=$gitDir/results/analysisResults/aDNA-ModernComparison/continuity/continuityOutput
########### modern reference population: california #######################
modernPop=CA

for category in ancient lowcovCA highcovCA
do
input=${modernPop}.angsdOut.mappedTo${ref}.${category}.counts.freqsFromModernGATK.superfile.1based.contInput.txt
# transversions only:
tvInput=${modernPop}.angsdOut.mappedTo${ref}.${category}.counts.freqsFromModernGATK.superfile.1based.contInput.TVOnly.txt
ancIDs=${category}.IDS.txt
outfile=${modernPop}.mappedTo${ref}.${category}.continuity.output.txt
tvOutfile=${modernPop}.mappedTo${ref}.${category}.continuity.output.TVOnly.txt

python $script $pathToContinuity $indir/$input $IDdir/${category}.IDs.txt $modernPop $outdir/$outfile
python $script $pathToContinuity $indir/$tvInput $IDdir/${category}.IDs.txt $modernPop $outdir/$tvOutfile

done

########### modern reference population: alaska #######################
modernPop=AK

for category in ancient
do
input=${modernPop}.angsdOut.mappedTo${ref}.${category}.counts.freqsFromModernGATK.superfile.1based.contInput.txt
# transversions only:
tvInput=${modernPop}.angsdOut.mappedTo${ref}.${category}.counts.freqsFromModernGATK.superfile.1based.contInput.TVOnly.txt
ancIDs=${category}.IDS.txt
outfile=${modernPop}.mappedTo${ref}.${category}.continuity.output.txt
tvOutfile=${modernPop}.mappedTo${ref}.${category}.continuity.output.TVOnly.txt

python $script $pathToContinuity $indir/$input $IDdir/${category}.IDs.txt $modernPop $outdir/$outfile
python $script $pathToContinuity $indir/$tvInput $IDdir/${category}.IDs.txt $modernPop $outdir/$tvOutfile

done

