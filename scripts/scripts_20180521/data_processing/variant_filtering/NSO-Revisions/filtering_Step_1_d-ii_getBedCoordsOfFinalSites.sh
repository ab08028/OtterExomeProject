source /u/local/Modules/default/init/modules.sh
module load bedtools
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering/
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/northern_sea_otter_genome/northern_sea_otter_genome.fasta
REFSHORTCODE=NSO # new thing ; a short code for the reference (if no code, it's mfur; otherwise NSO or SSO)

#### parameters:
rundate=20200719_${REFSHORTCODE} # date genotypes were called and ref code 20200719_NSO 
outdir=$wd/${rundate}_filtered
noCallFrac=1.0 # no filter # this is the tolerated no call fraction 20181008: doing this with the liberal filter instead of the stringent 0.2 filter
#sitesPassingAllFilters=${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} # old filtering scheme 
prefix='all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters'
sitesPassingAllFilters=${outdir}/${prefix}'_raw_variants.vcf.gz' # new filteirng scheme Oct 2018
mkdir -p $outdir/bedCoords

# want to get bed coords of all passing sites:
zcat $sitesPassingAllFilters | grep -v "#" | awk '{OFS="\t";print $1,$2-1,$2}' | bedtools sort -i stdin | bedtools merge -i stdin > ${outdir}/bedCoords/${prefix}.sorted.merged.coords.bed

# get total sum of covered sites:
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' ${outdir}/bedCoords/${prefix}.sorted.merged.coords.bed > ${outdir}/filteringStats/${prefix}.TOTALSITES.txt
