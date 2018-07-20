############## Find my sea otter capture targets in the ferret genome ##########

# I have the bed tools file with labels:
## note this is the bed file of the regions *I want* -- the actual deliverable from Nimblegen has buffers, combined regions, etc.
## ** on sirius ** ##
# temporarily move regions, and ferret genome to work dir (remove them after I'm done)
wd=/work2/abeichman/findFerretCaptureRegions_20180717
cp -r /data3/abeichman/ferret_genome_mfur/ $wd
cp /data3/abeichman/round5ExomeCaptureDesign_Oct2016/round5_AED1.0_newnames/FinalBEDS_round5ExomeCapture_Nov17_copiedfromAnnotationFolder/allCaptureRegions_IDCategoriesLengths_ForMyUse.Nov17.bed $wd
cp  /data3/abeichman/gidgetDeNovoGenome/DovetailAssembly/sea_otter_23May2016_bS9RH.fasta $wd # use original (not dedupped) genome
cp  /data3/abeichman/gidgetDeNovoGenome/DovetailAssembly/sea_otter_23May2016_bS9RH.fasta.fai $wd # use original (not dedupped) genome

# okay now I have everything
# So I want to pull these regions out of sea otter genome into a fasta file using bedtools
# then I want to blast them against ferret genome
# and be correctly stranded (luckily blast looks at both strands)
cd $wd
elutRef=$wd/sea_otter_23May2016_bS9RH.fasta # use full genome because that is what was used to make probe design; but be careful (but Nimblegen checked for stuff mapping multiple places)
ferretDB=$wd/ferret_genome_mfur/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa
regions=$wd/allCaptureRegions_IDCategoriesLengths_ForMyUse.Nov17.bed # not in perfect bed format -- it is 1-based and has extra columns.


# Need to make a better bed. Make it zero based and with better name column.
# this makes it zero based ($2-1, don't change $3 bc is noninclusive)
# and creates a name column that is 1-based (don't subtract from $2), that gives the category.gene.Scaffold:start-stop
# This will become the name in the fasta
# skip the header line: 
# need to not have 0-starts become negative
grep -v "scaffold" $regions | awk '{OFS="\t"; print($1,$2-1,$3,$5"."$4"."$1":"$2"-"$3)}' > $wd/allCaptureRegions_0Based_Names_UseThis.ProperFormat.20180717.bed
# but also just want the neutral/promoter regions (exclude exon regions)
grep -v "exon" allCaptureRegions_0Based_Names_UseThis.ProperFormat.20180717.bed > neutralPromoterCaptureRegions_noExon_0Based_Names_UseThis.ProperFormat.20180717.bed
# want to check that nothing got made negative 1 (if it started at 0)
allbed=$wd/allCaptureRegions_0Based_Names_UseThis.ProperFormat.20180717.bed
neutralbed=$wd/neutralPromoterCaptureRegions_noExon_0Based_Names_UseThis.ProperFormat.20180717.bed

sed -i'' 's/\t-1\t/\t0\t/g' $allbed
sed -i'' 's/\t-1\t/\t0\t/g' $neutralbed

bedtools getfasta -name -fi $elutRef -bed $allbed -fo $wd/allCaptureRegions.seaOtterSequence.nonStranded.20180717.fasta
bedtools getfasta -name -fi $elutRef -bed $neutralbed -fo $wd/neutralPromoterCaptureRegions_noExon.seaOtterSequence.nonStranded.20180717.fasta

mkdir -p $wd/blastResults
# then blast only neutral promoter sequences against ferret
# this is too slow: blastn -num_threads 30 -evalue 1e-10 -query $wd/neutralPromoterCaptureRegions_noExon.seaOtterSequence.nonStranded.20180717.fasta -db $ferretDB -outfmt 6 > $wd/blastResults/blast.seaOtterNeutralPromoterCaptureRegions_vs_ferretGenome.out

# note that max target seqs is part of the search algorithm, it doesn't filter at the end 
blastn -num_threads 10 -evalue 1e-10 -max_target_seqs 5 -query $wd/neutralPromoterCaptureRegions_noExon.seaOtterSequence.nonStranded.20180717.fasta -db $ferretDB -outfmt 6 > $wd/blastResults/blast.seaOtterNeutralPromoterCaptureRegions_vs_ferretGenome.top5.out

# do some prefiltering: select region lengths that are 500 bp (below that, seems like probably not mapping well)
awk '{if($4>=500) print}' $wd/blastResults/blast.seaOtterNeutralPromoterCaptureRegions_vs_ferretGenome.top5.out > $wd/blastResults/blast.seaOtterNeutralPromoterCaptureRegions_vs_ferretGenome.min500bp.out
# this retains 20,473 queries out of 21,383. Good enough. ~9K neutral regions
# and makes the file much more manageable.
# still need to get rid of duplicates. Look at in R
# turn neutral min 500 bp regions into bed file 
# blast output doesn't always orient start - stop such taht stop > start due to strandedness
# so if you want it to always have start < stop you need to use the below awk command 
grep neutral $wd/blastResults/blast.seaOtterNeutralPromoterCaptureRegions_vs_ferretGenome.min500bp.out | awk '{OFS="\t"; if($10>$9) print $2,$9-1,$10,$1; else print $2,$10-1,$9,$1}' > $wd/blastResults/blast.seaOtterNeutralCaptureRegions_vs_ferretGenome.min500bp.bed

# seeing how far neutral regions are from exons: must be at least 10kb
#bedtools closest -d -a ferret.Exon.Coordinates.0based.bed -b blastResults/blast.seaOtterNeutralCaptureRegions_vs_ferretGenome.min500bp.bed > neutralRegions.DistanceFromExons.ferret.bed
#sort -k9n neutralRegions.DistanceFromExons.ferret.bed | awk '{if($7!=-1)print $0;}' | head 

bedtools closest -d -b ferret.Exon.Coordinates.0based.bed -a blastResults/blast.seaOtterNeutralCaptureRegions_vs_ferretGenome.min500bp.bed > closestExonToNeutral.ferret.bed
sort -k9n closestExonToNeutral.ferret.bed | awk '{if($7!=-1)print $0;}' | less -S # check output
awk '{OFS="\t";if($9>=10000)print $1,$2,$3,$4}' closestExonToNeutral.ferret.bed | sort | uniq > ferret.NeutralRegions.gt10kbfromExon.minLen500bp.bed

### make some nice final table for these 
########## Ferret Exon Intervals: ###########
# get ferret gff file:
wget ftp://ftp.ensembl.org/pub/release-93/gff3/mustela_putorius_furo/Mustela_putorius_furo.MusPutFur1.0.93.gff3.gz
gunzip ftp://ftp.ensembl.org/pub/release-93/gff3/mustela_putorius_furo/Mustela_putorius_furo.MusPutFur1.0.93.gff3.gz
gff=$wd/Mustela_putorius_furo.MusPutFur1.0.93.gff3
# get exon coordinates to add to bed file (doing exon, rather than CDS, because I want UTRs)
# using Info column (9) as bed name 
grep exon $gff | awk '{OFS="\t";print $1,$4-1,$5,$9}' > ferretCoords.Exons.0based.bed
# make sure nothing got made -1
sed -i'' 's/\t-1\t/\t0\t/g' ferretCoords.Exons.0based.bed

# manually organized stuff and deleted ref files from work dir.
########## Transfer to Hoffman ###########
# neutral bed:
rsync ferret.NeutralRegions.gt10kbfromExon.minLen500bp.bed ab08028@hoffman2.idre.ucla.edu:/u/flashscratch/a/ab08028/captures/captureRegions
# exon bed:
rsync ferretCoords/ferret.Exon.Coordinates.0based.bed ab08028@hoffman2.idre.ucla.edu:/u/flashscratch/a/ab08028/captures/captureRegions
# distances bed for neutral regions (so I can categorize)
rsync bedtoolsClosest/closestExonToNeutral.ferret.bed ab08028@hoffman2.idre.ucla.edu:/u/flashscratch/a/ab08028/captures/captureRegions

# 7339 regions. Wish there were more...
# want closestExonToNeutral.ferret.bed this file and this file: ferret.NeutralRegions.gt10kbfromExon.minLen500bp.bed
# to work with
# talk to tanya about sorting neutral regions by distnace
# put onto Hoffman

#### Want to get coverage for: (note that I genotype everything with cov = 5 or greater)
# 1. ferretCoords.Exons.0based.bed # exons
# 2. ferret.NeutralRegions.gt10kbfromExon.minLen500bp.bed
# (promoters can wait)
