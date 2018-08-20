# modules
source /u/local/Modules/default/init/modules.sh
module load bedtools
module load blast

rundate=20180806

SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering/${rundate}_filtered/
mfurDir=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome
drerDir=/u/home/a/ab08028/klohmueldata/annabel_data/zebra_fish_genome
REFERENCE=$mfurDir/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
gff=$mfurDir/Mustela_putorius_furo.MusPutFur1.0.91.gff3
cdsRegions=$mfurDir/MusPutFuro1.0.91.cdsOnly.0based.sorted.merged.bed

exonicRegions=$mfurDir/MusPutFur1.0.91.exonsOnly.0based.sorted.bed # make this below if haven't already

####### get exonic regions (once) and make sure is sorted:
### (do once) grep exon $gff | awk '{OFS="\t";print $1,$4-1,$5,$9}' | sort -k1,1 -k2,2n > $exonicRegions
# should be ~200,000 lines

# results of filtering snps (all populations; all nv and snps)
hqSites=$wd/bedCoords/all_7_passingBespoke.sorted.merged.coords.bed # bed coords (sorted, merged) of sites from step 7 of filtering (comes from filtering_Step_2.sh)

############# Get distance of every set of sites from exonic regions in ferret genome ################
mkdir -p $wd/distanceFromExons
bedtools closest -d -a ${hqSites} -b ${exonicRegions} > $wd/distanceFromExons/all_7_passingBespoke.distanceFromExons.0based.txt
#### NOTE: the output of this will be
# [Hq site info] [closest exon info] [distance between]; so I want the HQ sites that are >10000bp away from the closest exon.
# don't want it the other way around (getting info on each exon). want info on hq sites.

# last column (8) is the distance; want it to be at least 10,000, and want to keep
# track of the distance. Collect all that are >10,000 away. 
# pick the ones with high distance (awk)
awk '{OFS="\t";if($8>10000)print $1,$2,$3}' $wd/distanceFromExons/all_7_passingBespoke.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/distanceFromExons/all_7_passingBespoke.min10kb.fromExon.0based.sorted.merged.bed

awk '{OFS="\t";if($8>100000)print $1,$2,$3}' $wd/distanceFromExons/all_7_passingBespoke.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/distanceFromExons/all_7_passingBespoke.min100kb.fromExon.0based.sorted.merged.bed

awk '{OFS="\t";if($8>1000000)print $1,$2,$3}' $wd/distanceFromExons/all_7_passingBespoke.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/distanceFromExons/all_7_passingBespoke.min1Mb.fromExon.0based.sorted.merged.bed

### Note: 1,2,3 columns are the HQ SITES position, NOT the position of the exon. (If you mess up what is a and b in bedtools closest this would be messed up)
# get amounts in each file:
> $wd/filteringStats/totalSequenceByDistanceFromExons.txt
for i in `ls $wd/distanceFromExons/*bed`
do
echo $i >> $wd/filteringStats/totalSequenceByDistanceFromExons.txt
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $i >> $wd/filteringStats/totalSequenceByDistanceFromExons.txt
done


###### for now use this 10kb bed file for making the SFS

##### re-merge the bed file with -d 10 setting (if things are 10bp apart you can still merge them)
# use this to get fasta
# doing this because if there is a single isolated base,
# it becomes its own entry in the fasta, which is going to make BLASTING a pain
# for now looking at overall region (allowing gaps of up to 10bp), not just called sites.
bedtools merge -d 10 $wd/distanceFromExons/all_7_passingBespoke.min10kb.fromExon.0based.sorted.merged.bed > $wd/distanceFromExons/all_7_passingBespoke.min10kb.fromExon.0based.sorted.mergedMaxDistance10.forFasta.notForSFS.bed

bedtools getfasta -fi $REFERENCE -bed $wd/distanceFromExons/all_7_passingBespoke.min10kb.fromExon.0based.sorted.merged.bed -fo $wd/distanceFromExons/putative.neutral.seqs.min10kb.fromExon.fasta
#  ********* YOU ARE HERE ************ 
# to do:  get GC content 
# get CpG content
# blast to fish
### you are here. it's an awkward fasta because some sites are only a single bp . need to merge them further?
# Could blast against sea otter and make sure don't land in cds region?
# Could blast against fish again  ???

############# Blast against zebra fish genome to look for conservation #############
# retrieved 20180620
# do this once:
# cd $drerDir
# wget ftp://ftp.ensembl.org/pub/release-86/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz
# gunzip Danio_rerio.GRCz10.dna.toplevel.fa.gz
# makeblastdb -in Danio_rerio.GRCz10.dna.toplevel.fa -out Drer_blastdb -dbtype nucl
blastn -query putative.neutral.seqs.fasta -db ../Drer_blastdb -outfmt 7 > neutralBlast_ZebraFish_blastn.out
# maybe also blast against sea otter cds regions?

# then want to choose the regions that are far (>10kb) from exonic regions; should hopefully be ~10,000
# May also want to exclude CpG islands, etc. from them (ask Tanya)
# then you can make a vcf based on this bed file using bedtools intersect
sleep 10m