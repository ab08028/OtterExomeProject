###### Gather up summary stats
# I want the fraction of aligned (hits) reads out of total retained reads
wd=/u/flashscratch/a/ab08028/captures/paleomix/testMapping
#headers=$SCRATCH/captures/samples/ancientSamples.txt
# 20181005: headers for library screens
headers=$SCRATCH/captures/samples/ancientLibs.Screen.txt
outdir=$SCRATCH/captures/paleomix/summaryStats
mkdir $outdir
outfile=$outdir/plmx.aDNA.summary.stats.screening.txt
REF=sea_otter_23May2016_bS9RH.deduped.99 # note that you'll collect ALL ref genomes for most stats; you're
#just setting the REF here for the MD coverage stats, which you're only collecting for the sea otter
# This script will still collect all other stats for all other ref genomes beause you aren't specifying 
# when you grep; so it'll grep all instances of seq_retained_reads for every genome from the summary file.
echo "sample statistic reference value" > $outfile
cat $headers | while read header
do
grep -w seq_retained_reads $wd/$header/*summary | awk '{print($1,$4,"noREF",$5)}' | sort | uniq >> $outfile
grep -w hits_raw $wd/$header/*summary | awk '{print($1,$4,$5)}' | sed 's/(/ /g' | sed 's/)//g' | sort | uniq >> $outfile
grep -w hits_raw_frac $wd/$header/*summary | awk '{print($1,$4,$5)}' | sed 's/(/ /g' | sed 's/)//g' | sort | uniq >> $outfile
grep -w hits_clonality $wd/$header/*summary | awk '{print($1,$4,$5)}' | sed 's/(/ /g' | sed 's/)//g' | sort | uniq >> $outfile
grep -w hits_unique $wd/$header/*summary | awk '{print($1,$4,$5)}' | sed 's/(/ /g' | sed 's/)//g' | sort | uniq >> $outfile
grep -w hits_unique_frac $wd/$header/*summary | awk '{print($1,$4,$5)}' | sed 's/(/ /g' | sed 's/)//g' | sort | uniq >> $outfile
grep -w hits_length $wd/$header/*summary | awk '{print($1,$4,$5)}' | sed 's/(/ /g' | sed 's/)//g' | sort | uniq >> $outfile
grep MD_001 -m1 -A1 $wd/$header/$header.$REF.depths | tail -n1 | awk -v REF=$REF '{print $1,"MD_001",REF,$7}' >> $outfile
grep MD_002 -m1 -A1 $wd/$header/$header.$REF.depths | tail -n1 | awk -v REF=$REF '{print $1,"MD_002",REF,$8}' >> $outfile
grep MD_003 -m1 -A1 $wd/$header/$header.$REF.depths | tail -n1 | awk -v REF=$REF '{print $1,"MD_003",REF,$9}' >> $outfile

done

# also want to gather MD_100; this is fraction of genome (not capture regions) with >=1x coverage; so multiply by genome length
# to get >=1x sites.
# note that if you don't use -w then hits_raw and hits_raw_frac will both be extracted (same with unique); this can mess up your
# bar graphs in R due to duplicate entries (smart to remove dups in R just in case anyway)
