######## get cds and neutral regions from super bed file ###########

module load bedtools
# so far only works for mfur; need to get coordinates for elut (later)
ref=mfur
neutBed=
cdsBed=
GPSuperfile
GLSuperfile
basename=
# -wa says to print every entry of "a" that intersects with "b"; -header says to print the header of "a" before the results
# note that bed headers must start with "#"
bedtools intersect -a $GPSuperfile -wa -header -b $neutBed | gzip > ${basename}.??.neutral.bed.gz

bedtools intersect -a $GPSuperfile -wa -header -b $cdsBed | gzip > ${basename}.??.cds.bed.gz

