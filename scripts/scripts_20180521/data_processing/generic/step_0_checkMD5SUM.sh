#### Little script to check md5sum for all

DATE=`date +%Y%m%d`
> md5sum.${DATE}.txt
for i in `ls *fastq.gz`
do
echo $i
md5sum $i >> md5sum.${DATE}.txt
done

# then compare md5sums with R
# this checks if the lines are the same and prints the ones that aren't

grep -vFxf md5sum.${DATE}.txt md5sum.txt > non-matching.md5sum.${DATE}.txt

# should be empty!

