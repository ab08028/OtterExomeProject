#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=20G,highp
#$ -N sift
#$ -o /u/flashscratch/a/ab08028/captures/reports/sift
#$ -e /u/flashscratch/a/ab08028/captures/reports/sift
#$ -m abe
#$ -M ab08028

source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111


####  Run sift
# sift jar:
sift=/u/home/a/ab08028/klohmueldata/annabel_data/bin/sift/SIFT4G_Annotator_v2.3.jar
# mustela db:
mustelaDB=/u/home/a/ab08028/klohmueldata/annabel_data/bin/sift-MustelaDatabase/MusPutFur1.0.83
# variants (both homozygous and heterozygous, relative to ferret)
bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip
# must use bgzip to zip vcfs

# all sites, each pop separately, filtered out all het. 
wd=/u/flashscratch/a/ab08028/captures/
genotypeDate=20180806
vcfdir=$wd/vcf_filtering/${genotypeDate}_filtered/populationVCFs
suffix=all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled # this vcf is per-population but has sites that are all-het removed, relatives removed, admixed removed, passing all filters
outdir=$wd/analyses/sift/${genotypeDate}
mkdir -p $outdir

pops="CA AK AL COM KUR"

for pop in $pops
do
vcf=${pop}_${suffix}.vcf
gunzip $vcf.gz
# http://sift.bii.a-star.edu.sg/sift4g/Commandline.html
java -jar -Xmx10G $sift -c -i $vcfdir/$vcf -d $mustelaDB -r $outdir -t
$bgzip $vcf
done