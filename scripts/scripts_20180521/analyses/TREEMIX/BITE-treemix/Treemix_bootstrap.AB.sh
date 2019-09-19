###########################################
##### Treemix pipeline with bootstrap #####
#License: GPLv3 or later
#Modification date: 2017-05-10
#Written by: Elia Vajana, Marco Milanesi
#Contact: Elia Vajana <vajana.elia@gmail.com>, Marco Milanesi <marco.milanesi.mm@gmail.com>
#Description: Run Treemix software for a specific migration number previously chosen applying a bootstrap procedure 
###########################################
infile=$1 		## treemix input file
numk=$2 		## the migration number you choose and you want test
ncore=$3 		## max number of cores to use
blockk=$4 		## block size
outgroup=$5 	## name of the selected outgroup population (if you want to do an unrooted ML tree put here 'NoOutgroup' (without quotes))
nboot=$6		## number of bootstrap replicates of the TREE
pathP=$7		## path to Phylip consense program. Example: /biosoftware/phylip/phylip-3.696/exe/consense
outname=$8		## name for output file

##########################################################################################
#########################
##### Settings file #####
#########################
echo "Input file name = $1" > $outname"_Settings.txt"
echo "Output file name = $8" >> $outname"_Settings.txt"
echo "Number of migration modeled = $2" >> $outname"_Settings.txt"
echo "Number of SNPS per block = $4" >> $outname"_Settings.txt"
if [ $outgroup = "NoOutgroup" ]; then
	echo "Unrooted ML tree" >> $outname"_Settings.txt"
else
	echo "Outgroup = $5" >> $outname"_Settings.txt"
fi
echo "Number of replications = $6" >> $outname"_Settings.txt"

mkdir bootstrap



#############################
#### BOOTSTRAP PROCEDURE ####
#############################

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### N.B. If you need to modify the treemix parameters please do it here ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

### Seed number settings ###
echo "**** RUNNING TREEMIX ****"
if [ $outgroup = "NoOutgroup" ]; then
	dowork() { 
		a=$RANDOM
		b=1`date +%N`
		let "c=$a+$b"
		treemix -i $2 -k $3 -se -m $4 -seed $c -o "bootstrap/"$5"_treemix_bootrep_"$1
	}
	export -f dowork
	seq 1 $nboot | dowork {} $infile $blockk $numk $outname > $outname"_logfile_treemix_bootrep.log"
else
	dowork() { 
		a=$RANDOM
		b=1`date +%N`
		let "c=$a+$b"
		treemix -i $2 -k $3 -se -m $4 -root $6 -seed $c -o "bootstrap/"$5"_treemix_bootrep_"$1
	}
	export -f dowork
	seq 1 $nboot | dowork {} $infile $blockk $numk $outname $outgroup > $outname"_logfile_treemix_bootrep.log"
fi
echo "**** RUNNING TREEMIX: DONE ****"


### Create a file with all the bootstrapped trees
for a in `seq 1 $nboot`
do
	bootf="bootstrap/"$outname"_treemix_bootrep_"$a".treeout.gz"
	gunzip -c $bootf | head -1 >> $outname"_boottree.tre"
done

echo "   ***** Boostrap procedure: DONE *****"



#########################################################################
#### Run Phylip on the bootstrapped trees to obtain a consensus tree ####
#########################################################################
echo "   ***** Phylip - consensus tree construction: START *****"
### Clean the environment
rm -rf outfile outtree screanout

# Create parameters file
if [ $outgroup = "NoOutgroup" ]; then
	echo $outname"_boottree.tre" > $outname"_PhylipInputFile"
	echo "Y" >> $outname"_PhylipInputFile"
else
	# Find the position of Outgroup population
	posOutgroup=`head -1 $outname"_boottree.tre" | tr "," "\n" | grep $outgroup -n | cut -d":" -f1`
	# echo $posOutgroup
	echo $outname"_boottree.tre" > $outname"_PhylipInputFile"
	echo "O" >> $outname"_PhylipInputFile"
	echo $posOutgroup >> $outname"_PhylipInputFile"
	echo "Y" >> $outname"_PhylipInputFile"
fi

# Run Phylip
$pathP < $outname"_PhylipInputFile" > screanout

### The output from Phylip will be modified because:
### 1) Treemix accept only one line tree 
### 2) Treemix accept newick format file 
#sed ':a;N;$!ba;s/\n//g' outtree > $outname"_outtree.newick"
cat outtree | tr -d "\n" > $outname"_outtree.newick"
echo >> $outname"_outtree.newick"
echo "   ***** Phylip - consensus tree construction: DONE *****"


######################################################################################
### Run TreeMix with the chosen number of migrations by loading the consensus tree ###
######################################################################################

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### N.B. If you need to modify the treemix parameters please do it here ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

echo "**** RUNNING TREEMIX with Consensus tree ****"
if [ $outgroup = "NoOutgroup" ]; then
	treemix -i $infile -m $numk -k $blockk -se -tf $outname"_outtree.newick" -o $outname > $outname"_logfile_treemix_boot.log"
else
	treemix -i $infile -m $numk -k $blockk -root $outgroup -se -tf $outname"_outtree.newick" -o $outname > $outname"_logfile_treemix_boot.log"
fi
echo "**** RUNNING TREEMIX: DONE ****"

echo "TreeMix - Bootstrap Analysis --> DONE"

