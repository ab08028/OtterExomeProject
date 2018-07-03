### do something like what I did for the exome capture;


# closest exonic region to neutral region is 10001 bp ## need the -d option!! ## 
#bedtools closest -a NeutralRegions_tab.bed -b ../Post_doc_ChevironLab/exome_array/Peromyscus_array_design/files_from_Gideon/pman_excap_sorted_for_order_tab.bed -d > closest.bed
exonicRegions=[make a bed file of ferret exonic regions ; maybe +- 1000bp?]
coveredRegions=[some sort of consensus set of regions where variants were called for most/all individuals]

bedtools closest -d -a ${exonicRegions} -b ${coveredRegions} > closest.bed

sort -k8n closest.bed | awk '{if($7!=-1)print $0;}' | head 
# last column (8) is the distance; want it to be at least 10,000, and want to keep
# track of the distance. Collect all that are >10,000 away. 
#NW_006501281.1	4653851	4654851	NW_006501281.1	4643138	4643851	10001
#ScbS9RH_95815   316271  316294  ELUT_00004602   ScbS9RH_95815   326201  327201  9908
# pick the ones with high distance (awk)
awk '{if($8>=10000)print $0}' # ... something like this to select 'neutral' regions

# But maybe have to do some other checks? Because they could have been annotated in sea otter?
# Could blast against sea otter and make sure don't land in cds region?
# Could blast against fish again 

bedtools getfasta -fi $REFERENCE -bed $neutralBed -fo putative.neutral.seqs.fasta

# Could blast against sea otter and make sure don't land in cds region?
# Could blast against fish again 
blastn -query putative.neutral.seqs.Oct03.noSnapMask.WITH2ZEBRAFISHTESTSEQS.fasta -db ../Drer_blastdb -outfmt 7 > neutralBlast_ZebraFish_blastn.out
# maybe also blast against sea otter cds regions?

# then want to choose the regions that are far (>10kb) from exonic regions; should hopefully be ~10,000
# May also want to exclude CpG islands, etc. from them (ask Tanya)
# then you can make a vcf based on this bed file using bedtools intersect
sleep 10m