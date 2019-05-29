### convert vep output back to bed coords
# note vep is 1based, bed is 0based


input=#vep file filtered for NS, S or LOF

grep -v "#" $input | awk '{OFS="\t";split($2,pos,":");print pos[1],pos[2]-1,pos[2],$1}' > ${input%.tbl}.0based.bed


# or possibly can just izip it with my superfile and work with it that way
