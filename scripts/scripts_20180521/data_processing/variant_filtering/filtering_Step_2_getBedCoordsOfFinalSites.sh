source /u/local/Modules/default/init/modules.sh
module load bedtools
#### parameters:
rundate=20180806 # date genotypes were called
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering/
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
outdir=$wd/${rundate}_filtered
noCallFrac=0.2 # this is the tolerated no call fraction
# if you want to use this with your new 90% filtering, just have to change this to 0.9? or use the final all sites file (all_9)?
sitesPassingAllFilters=${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile}
mkdir -p $outdir/bedCoords

# want to get bed coords of all passing sites:
zcat $sitesPassingAllFilters | grep -v "#" | awk '{OFS="\t";print $1,$2-1,$2}' | bedtools sort -i stdin | bedtools merge -i stdin > ${outdir}/bedCoords/all_7_passingBespoke.sorted.merged.coords.bed

# get total sum of covered sites:
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' ${outdir}/bedCoords/all_7_passingBespoke.sorted.merged.coords.bed > ${outdir}/filteringStats/all_7_passingBespoke.TOTALSITES.txt
