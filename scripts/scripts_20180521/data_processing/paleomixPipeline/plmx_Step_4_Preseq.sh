### Want to look at library complexity of ANCIENT only 

# this runs fast; can be run in the shell with qrsh

module load samtools
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/paleomix/
headers=$SCRATCH/captures/samples/ancientSamples.txt
REFPREFIX=Mustela_putorius_furo.MusPutFur1.0.dna.toplevel

preseq=/u/home/a/ab08028/klohmueldata/annabel_data/bin/preseq_v2.0/preseq

mkdir -p $SCRATCH/captures/preseq
outdir=$SCRATCH/captures/preseq

cat $headers | while read header
do 
# get the actual curve of the library from subsampling:
preseq c_curve -H ${wd}/${header}/${header}.${REFPREFIX}.duphist -o $outdir/${header}.preseq.cCurve.list
# get extrapolated curve if you sequenced the library more
preseq lc_extrap -H ${wd}/${header}/${header}.${REFPREFIX}.duphist -o $outdir/${header}.preseq.lcExtrap.list
done
# then am going to make plots in R (plot # of unique reads vs. # of reads)
# may need to adjust step size 