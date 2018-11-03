#! /bin/bash
#$ -cwd
#$ -l h_rt=100:00:00,h_data=16G,highp,arch=intel*
#$ -N vcf1e_filtering
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028

# modules
source /u/local/Modules/default/init/modules.sh
module load java
module load python/2.7
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix
bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip

#### parameters:
rundate=20180806 # date genotypes were called (vcf_20180806 includes capture 02)
noCallFrac=0 # maximum fraction of genotypes that can be "no call" (./.) # not allowing any no-call 


#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
mkdir -p $wd
indir=$SCRATCH/captures/vcfs/vcf_${rundate}
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering
vcfdir=$wd/${rundate}_filtered # date you called genotypes

bespokeFilterScript=$scriptdir/filtering_perPopulation.noCall.allHet.py

for pop in CA KUR COM AK AL
do
echo $pop
# usage: python [script] [full path to invcf] [full path to out vcf] [full path to error file] [max no call fraction]
# 20181102 updated it to 
# note that my bespoke script does the nocallfrac so you don't have to do it in step 8! 
python $bespokeFilterScript ${vcfdir}/populationVCFs/${pop}_'all_8_rmRelativesAdmixed_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
${vcfdir}/populationVCFs/${pop}_'all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled.vcf' \
${vcfdir}/populationVCFs/${pop}_'fail_all_9_FAILING_perPopulationBespoke_Filters_'${infile%.vcf.gz}.txt \
$noCallFrac
# bgzip the result: (note: must use bgzip not gzip)
$bgzip ${vcfdir}/populationVCFs/${pop}_'all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled.vcf'
$tabix -p vcf ${vcfdir}/populationVCFs/${pop}_'all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz' # index the vcf

done


# also filter the admixed vcfs:
for pop in KUR AK
do
echo "admixed " $pop
# usage: python [script] [full path to invcf] [full path to out vcf] [full path to error file] [max no call fraction]

python $bespokeFilterScript ${vcfdir}/populationVCFs/admixedVCFs/admixIndOnly_${pop}_'all_8_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
${vcfdir}/populationVCFs/admixedVCFs/admixIndOnly_${pop}_'all_9_rmAllHet_passingAllFilters_allCalled.vcf' \
${vcfdir}/populationVCFs/admixedVCFs/admixIndOnly_${pop}_'fail_all_9_FAILING_perPopulationBespoke_Filters_'${infile%.vcf.gz}.txt \
$noCallFrac
# bgzip the result: (note: must use bgzip not gzip)
$bgzip ${vcfdir}/populationVCFs/admixedVCFs/admixIndOnly_${pop}_'all_9_rmAllHet_passingAllFilters_allCalled.vcf'
$tabix -p vcf ${vcfdir}/populationVCFs/admixedVCFs/admixIndOnly_${pop}_'all_9_rmAllHet_passingAllFilters_allCalled.vcf.gz' # index the vcf

done