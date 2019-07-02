#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G
#$ -m abe
#$ -M ab08028
#$ -N HapFile2VEP
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

## goal of this script is to convert the pseudohaploids file into a bed file ####
# directories: 
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/

SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
hapDir=$wd/angsd-pseudoHaps

ref=mfur # only work with MFUR when going into VEP!
#dates="20190524-highcov-AFprior 20190524-lowcov-AFprior 20190524-highcov-UNIFprior 20190524-lowcov-UNIFprior" # set of angsdDates you want to process 
#dates="20190612-highcov-pseudoHaps 20190612-lowcov-pseudoHaps"
# have corresponding maf and count files that have exact same filters (eventually want to make a script that makes these at the same time)
# these have EXACT same filters and same # of sites -- eventually want to combine the forming of these files so they can't be offset
############################################### high coverage ################################################
###### these were hand checked that they have the same sites contained; if diff filters had been used these wouldn't be safe to combine
# but it's okay because I check # of sites in this script, and I check that the marker IDs are the same using my python parsing script
# so if you accidentally concat things that shouldn't, it should be caught (but don't ;) ) 
#hapDate="20190612-highcov-pseudoHaps"
#mafCountDate="20190524-highcov-AFprior"
#basename=angsdOut.mappedTo${ref}

#indir=$hapDir/$hapDate # location of your posterior probs
#GLDir=$wd/angsd-GLs/$mafCountDate
#outdir=$wd/VEP/pseduoHaps/$hapDate
#mkdir -p $outdir
# can use GP or GL superfile; they are in the same order and have same sites, so it doesn't matter, since GLs and GPs aren't taken along for the ride in the VEP input format
# I am using GPs, but again, could be either:
#cdsSuperfile=${basename}.superfile.GPs.mafs.counts.cdsOnly.0based.bed.gz
#hapFile=${basename}.haplo.gz
#hapoutput=${hapFile%.haplo.gz}.pseudoHaps.superfile.0based.bed
#hapVEPOutput=${hapFile%.haplo.gz}.1based.VEP.input.txt
#mafs=${basename}.mafs.gz
#counts=${basename}.counts.gz
#hapSites=`zcat $indir/$hapFile | wc -l`
#mafSites=`zcat $GLDir/$mafs | wc -l`
#if [ $hapSites -eq $mafSites ]
#then
#echo "they have the same number of sites"
#else 
#echo "watch out -- differing number of sites"
#fi

# w
# so want to add a bed header too, which will be 
#angsdheaders=`paste <(zcat $GLDir/$mafs | head -n1) <(zcat $indir/$hapFile | head -n1) <(zcat $GLDir/$counts | head -n1)` # this is the headers of all the files pasted together
# so want to add a bed header too, which will be 
#bedhead="#chrom\tstart0based\tend\tmarkerID\tempty5\tempty6\tempty7\tempty8\tempty9\tempty10\tempty11\tempty12"

#comboheader=`echo -e "$bedhead\t$angsdheaders"` # need the "" and -e to get the tabs in
#echo -e "$comboheader" >  $indir/${hapoutput} 


#paste <(zcat $GLDir/$mafs) <(zcat $indir/$hapFile) <(zcat $GLDir/$counts) | grep -v "chromo" | awk '{OFS="\t";print $1,$2-1,$2,$1"_"$2,".",".",".",".",".",".",".",".",$0}' | sed 's/\t\t/\t/g' | sed 's/\t$//g' >> $indir/${hapoutput} # go into awk and rearrange to make it bed format with extra columns 
#gzip -f $indir/${hapoutput}

############################################### low coverage ################################################
###### these were hand checked that they have the same sites contained; if diff filters had been used these wouldn't be safe to combine
# but it's okay because I check # of sites in this script, and I check that the marker IDs are the same using my python parsing script
# so if you accidentally concat things that shouldn't, it should be caught (but don't ;) ) 
hapDate="20190612-lowcov-pseudoHaps"
mafCountDate="20190524-lowcov-AFprior"
basename=angsdOut.mappedTo${ref}

indir=$hapDir/$hapDate # location of your posterior probs
GLDir=$wd/angsd-GLs/$mafCountDate
outdir=$wd/VEP/pseduoHaps/$hapDate
mkdir -p $outdir
# can use GP or GL superfile; they are in the same order and have same sites, so it doesn't matter, since GLs and GPs aren't taken along for the ride in the VEP input format
# I am using GPs, but again, could be either:
#cdsSuperfile=${basename}.superfile.GPs.mafs.counts.cdsOnly.0based.bed.gz
hapFile=${basename}.haplo.gz
hapoutput=${hapFile%.haplo.gz}.pseudoHaps.superfile.0based.bed
hapVEPOutput=${hapFile%.haplo.gz}.1based.VEP.input.txt
mafs=${basename}.mafs.gz
counts=${basename}.counts.gz
hapSites=`zcat $indir/$hapFile | wc -l`
mafSites=`zcat $GLDir/$mafs | wc -l`
if [ $hapSites -eq $mafSites ]
then
echo "they have the same number of sites"
else 
echo "watch out -- differing number of sites"
fi

# w
# so want to add a bed header too, which will be 
angsdheaders=`paste <(zcat $GLDir/$mafs | head -n1) <(zcat $indir/$hapFile | head -n1) <(zcat $GLDir/$counts | head -n1)` # this is the headers of all the files pasted together
# so want to add a bed header too, which will be 
bedhead="#chrom\tstart0based\tend\tmarkerID\tempty5\tempty6\tempty7\tempty8\tempty9\tempty10\tempty11\tempty12"

comboheader=`echo -e "$bedhead\t$angsdheaders"` # need the "" and -e to get the tabs in
echo -e "$comboheader" >  $indir/${hapoutput} 


paste <(zcat $GLDir/$mafs) <(zcat $indir/$hapFile) <(zcat $GLDir/$counts) | grep -v "chromo" | awk '{OFS="\t";print $1,$2-1,$2,$1"_"$2,".",".",".",".",".",".",".",".",$0}' | sed 's/\t\t/\t/g' | sed 's/\t$//g' >> $indir/${hapoutput} # go into awk and rearrange to make it bed format with extra columns 
gzip -f $indir/${hapoutput}


