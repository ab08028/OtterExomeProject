###### Gather up summary stats
# I want the fraction of aligned (hits) reads out of total retained reads
wd=/u/flashscratch/a/ab08028/captures/paleomix/
headers=$SCRATCH/captures/samples/modernSamples.txt

outdir=$SCRATCH/captures/paleomix/summaryStats
mkdir $outdir
outfile=$outdir/plmx.modernDNA.summary.stats.txt
echo "sample statistic reference value" > $outfile

cat $headers | while read header
do
# need to make spaces explicit 
grep -w seq_retained_reads $wd/$header/*summary | awk '{print($1,$4,"noREF",$5)}' | sort | uniq >> $outfile
grep -w hits_raw $wd/$header/*summary | awk '{print($1,$4,$5)}' | sed 's/(/ /g' | sed 's/)//g' | sort | uniq >> $outfile
grep -w hits_raw_frac $wd/$header/*summary | awk '{print($1,$4,$5)}' | sed 's/(/ /g' | sed 's/)//g' | sort | uniq >> $outfile
grep -w hits_clonality $wd/$header/*summary | awk '{print($1,$4,$5)}' | sed 's/(/ /g' | sed 's/)//g' | sort | uniq >> $outfile
grep -w hits_unique $wd/$header/*summary | awk '{print($1,$4,$5)}' | sed 's/(/ /g' | sed 's/)//g' | sort | uniq >> $outfile
grep -w hits_unique_frac $wd/$header/*summary | awk '{print($1,$4,$5)}' | sed 's/(/ /g' | sed 's/)//g' | sort | uniq >> $outfile
grep -w hits_length $wd/$header/*summary | awk '{print($1,$4,$5)}' | sed 's/(/ /g' | sed 's/)//g' | sort | uniq >> $outfile
# don't need these stats for modern
#grep MD_001 -m1 -A1 $wd/$header/$header.$REF.depths | tail -n1 | awk -v REF=$REF '{print $1,"MD_001",REF,$7}' >> $outfile
#grep MD_002 -m1 -A1 $wd/$header/$header.$REF.depths | tail -n1 | awk -v REF=$REF '{print $1,"MD_002",REF,$8}' >> $outfile
#grep MD_003 -m1 -A1 $wd/$header/$header.$REF.depths | tail -n1 | awk -v REF=$REF '{print $1,"MD_003",REF,$9}' >> $outfile
done

# note that if you don't use -w then hits_raw and hits_raw_frac will both be extracted (same with unique); this can mess up your
# bar graphs in R due to duplicate entries (smart to remove dups in R just in case anyway)