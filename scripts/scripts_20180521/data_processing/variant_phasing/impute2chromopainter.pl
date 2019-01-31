# this is from https://people.maths.bris.ac.uk/~madjl/finestructure-old/chromopainter_info.html
#!/usr/bin/perl
## CONVERTS PHASED IMPUTE2 OUTPUT TO CHROMOPAINTER-STYLE INPUT FILES

sub help {
print("CONVERTS PHASED SHAPEIT/IMPUTE2 OUTPUT TO CHROMOPAINTER-STYLE INPUT FILES\n");

print("usage:   perl impute2chromopainter.pl <options> <-r recommap_filein> impute_output_file.haps output_filename_prefix\n");

print("where:\n");
print("        (i) impute_output_file.haps = filename of IMPUTE2 output file with suffix \".haps\" that contains phased haplotypes\n");
print("        (ii) output_filename_prefix = filename prefix for chromopainter input file(s). The suffix \".phase\" is added\n\n");
print("The output, by default, is in CHROMOPAINTER v2 input format.\n");

print("<options>:\n");
print("-J:                 Jitter (add 1) snp locations if snps are not strictly ascending. Otherwise an error is produced.\n");
print("-r recommap_filein: Filename of genetic map used to phase haplotypes in IMPUTE2/shapeit. If provided, a \".recomrates\" file is also produced.\n");

print("<further options>   NOTE: YOU ONLY NEED THESE OPTIONS FOR BACKWARDS COMPATABILITY!\n");
print("-v1:                Produce output compatible with CHROMOPAINTER v1, i.e. include the line of \"S\" for each SNP. \n");
print("-f:                 By default, this script produces PHASE-style output, which differs from \n");
print("			   ChromoPainter input which requires an additional first line.  This option creates the correct\n");
print("		           first line for standard fineSTRUCTURE usage (i.e. the first line is \"0\", all other lines are appended)\n\n");

print("NOTE: TO USE IN CHROMOPAINTER: You also need a recombination map. Create this with the \"-r\" option, or use the \"convertrecfile.pl\" or \"makeuniformrecfile.pl\" scripts provided.\n\n");
print(" !!! WARNING:  THIS PROGRAM DOES NOT SUFFICIENTLY CHECK FOR MISSPECIFIED FILES. WE ARE NOT ACCOUNTABLE FOR THIS RUNNING INCORRECTLY !!!\n");
die "\n";
}

use Switch;


###############################
## INPUT:

$numindsmaxsize=500;     ## ONLY READ IN THIS NUMBER OF IMPUTE2 INDS AT A TIME (TO MINIMIZE RAM REQUIREMENTS -- MAKES PROGRAM RUN A FACTOR OF (N/$numindsmaxsize) SLOWER THAN IF YOU SET $numindsmaxsize=N, where N is the total number of individuals in your ".haps" IMPUTE2 output file, but also uses a factor of (N/$numindsmaxsize) less RAM

###############################
## ARGUMENT PROCESSING

$v1=0; ## version 1 mode
$fsmode=0; ## finestructure mode (i.e. start with an additional line containing 0)
$jitter=0; ## whether we jitter snp locations

$Mb = 1000000.0;
$IMPUTEinfile="";
$recommapinfile="";
$outfilePRE="";

$argon=0;
for (my $i = 0; $i < scalar(@ARGV); ++$i){
	if(@ARGV[$i] eq "-f"){
	    $fsmode=1;
	}elsif(@ARGV[$i] eq "-r"){
	    $recommapinfile=$ARGV[$i++];
	}elsif(@ARGV[$i] eq "-v1"){
	    $v1=1;
	}elsif(@ARGV[$i] eq "-J"){
	    $jitter=1;
	}else{
	    switch($argon){
		case 0 {$IMPUTEinfile="$ARGV[$i]";}
		case 1 {$outfilePRE="$ARGV[$i]";}
		else {
		    help();
		}
	    }
	    $argon++;
	}
}

if($outfilePRE eq "" || $argon != 2) {help();}
$outfilePRE =~ s/.phase$//;

##############################
## PROGRAM:

if("$recommapinfile" ne ""){
    ## (I) GET RECOM-RATE INFO:
    open(IN,"$recommapinfile");
    $line=<IN>;        ## header
    @posvecFULL=();
    @recomrateFULL=();
    $numfoundind=0;
    $prevpos=0;
    while(<IN>)
    {
	$line=$_;
	@linearray=split(/\s+/,$line);
	push(@posvecFULL,$linearray[0]);
	push(@recomrateFULL,$linearray[1]);
	if ($linearray[0] < $prevpos)
	{
	    print "genetic map error -- positions out of order ($linearray[0] < ${prevpos})\n";
	    die;
	}
	$prevpos = $linearray[0];
    }
    $nsitesFULL=@posvecFULL;
}

## (II) GET NUMBER OF SITES AND INDS: 
open(IN,"$IMPUTEinfile");
$line=<IN>;
@linearray=split(/\s+/,$line);
$totalINDS=(@linearray-5)/2;
$totalhaps=@linearray-5;
$numsplits=int($totalINDS/$numindsmaxsize);
if (($numsplits*$numindsmaxsize)<$totalINDS)
{
    $numsplits=$numsplits+1;
}
$nsites=1;
while(<IN>)
{
    $line=$_;
    $nsites=$nsites+1;
}
#print "$numsplits $nsitesFULL $totalINDS\n";

              ## (III) READ IN IMPUTE2 HAPLOTYPES AND MAKE CHROMOPAINTER HAPLOTYPE INPUT FILE:
open(OUT,">${outfilePRE}.phase");
if($fsmode==1) {
	print OUT "0\n";
}
for ($a=0; $a < $numsplits; $a+=1)
{
    $startIND=$a*$numindsmaxsize;
    $endIND=($a+1)*$numindsmaxsize;
    if ($endIND > $totalINDS)
    {
	$endIND=$totalINDS;
    }

                              ## read in:
    open(IN,"$IMPUTEinfile");
    @rsvec=();
    @posvec=();
    @genomat=();
    $snpcount=0;
    while(<IN>)
    {
	$line=$_;
	@linearray=split(/\s+/,$line);
	push(@rsvec,$linearray[1]);
#	print("BEFORE: SNP location $linearray[2] lup $lastuniquesnp\n");
	if(scalar(@posvec)>0 &&($linearray[2] <= $posvec[-1]) && $linearray[2]>=0){
	    if(!$jitter){
		die("ERROR: SNPs are not strictly ascending, exiting. Rerun with -J to jitter the SNP locations.\n");
	    }
#	    print("Duplication found: setting $linearray[2] to $posvec[-1]+1\n");
	    $linearray[2]=$posvec[-1]+1;
	}
#	print("AFTER:  SNP location $linearray[2] lup $lastuniquesnp\n");
	if( $linearray[2] == $posvec[-1] ){
	    die("Strange error due to jittering?\n");
	}
	push(@posvec,$linearray[2]);
	shift(@linearray);
	shift(@linearray);
	shift(@linearray);
	shift(@linearray);
	shift(@linearray);

	for ($i=$startIND; $i < $endIND; $i+=1)
	{
	    $genomat[(($i-$startIND)*2)][$snpcount]=$linearray[(2*$i)];
	    $genomat[(($i-$startIND)*2+1)][$snpcount]=$linearray[(2*$i+1)];
	}

	$snpcount=$snpcount+1;
    }

                                ## print out:	
    if ($a==0)
    {
	if($v1){
	    print OUT "$totalINDS\n";
	}else {
	    print OUT "$totalhaps\n";
	}
	print OUT "$nsites\n";
	print OUT "P @posvec\n";
	if($v1){
	    for ($j=0; $j < $nsites; $j+=1)
	    {
		print OUT "S";
	    }
	    print OUT "\n";
	}
    }
    for ($i=0; $i < (2*($endIND-$startIND)); $i+=1)
    {
	for ($j=0; $j < $nsites; $j+=1)
	{
	    print OUT "$genomat[$i][$j]";
	}
	print OUT "\n";
    }
}

if($recommapinfile ne ""){
    ## (IV) MAKE CHROMOPAINTER RECOM-MAP INPUT FILE:
    $start = 0;
    @recomrate=();
    if ($posvec[0] < $posvecFULL[0])
    {
	#print "ERROR - first basepair of file is less than first basepair of genetic map!\n";
    }
    for ($i=0; $i < ($nsites-1); $i+=1)
    {
	push(@recomrate,0);
	if ($posvec[$i] >= $posvecFULL[($nsitesFULL-1)])
	{
	    $recomrate[$i] = $recomrateFULL[($nsitesFULL-2)]/$Mb;
	}
	if ($posvec[$i] < $posvecFULL[($nsitesFULL-1)])
	{
	    for ($j=$start; $j < ($nsitesFULL-1); $j+=1)
	    {
		if (($posvec[$i] >= $posvecFULL[$j]) && ($posvec[$i] < $posvecFULL[($j+1)]) && ($posvec[($i+1)] <= $posvecFULL[($j+1)]))
		{
		    $recomrate[$i] = $recomrateFULL[$j]/$Mb;  
		    $start = $j;
		    last;
		}
		if (($posvec[$i] >= $posvecFULL[$j]) && ($posvec[$i] < $posvecFULL[($j+1)]) && ($posvec[($i+1)] > $posvecFULL[($j+1)]))
		{
		    $recomcurrent = $recomrateFULL[$j]*($posvecFULL[($j+1)]-$posvec[$i]);
		    $endspot = $j+1;
		    if ($endspot == ($nsitesFULL-1))
		    {
			$recomcurrent = $recomcurrent + $recomrateFULL[$j]*($posvec[($i+1)]-$posvecFULL[($j+1)]);
			$recomrate[$i] = ($recomcurrent/($posvec[($i+1)]-$posvec[$i]))/$Mb;  
			last;
		    }
		    while($posvec[($i+1)] > $posvecFULL[($endspot+1)])
		    {
			$recomcurrent = $recomcurrent + $recomrateFULL[$endspot]*($posvecFULL[($endspot+1)]-$posvecFULL[$endspot]);
			$endspot = $endspot+1;
			if ($endspot == ($nsitesFULL-1))
			{
			    last;
			}
		    }
		    $recomcurrent = $recomcurrent + $recomrateFULL[$endspot]*($posvec[($i+1)]-$posvecFULL[$endspot]);
		    $recomrate[$i] = ($recomcurrent/($posvec[($i+1)]-$posvec[$i]))/$Mb;  
		    $start = $j;
		    last;
		}
	    }
	}
	if ($recomrate[$i] == 0)
	{
	    #print "error - $i\n";
	    if(($posvec[$i] < $posvecFULL[0]) && ($posvec[($i+1)] <= $posvecFULL[1]))
	    {
		$recomrate[$i] = $recomrateFULL[0]/$Mb;
	    }
	    if(($posvec[$i] < $posvecFULL[0]) && ($posvec[($i+1)] > $posvecFULL[1]))
	    {
		$recomcurrent = $recomrateFULL[0]*($posvecFULL[1]-$posvec[$i]);
		$endspot = 1;
		while($posvec[($i+1)] > $posvecFULL[($endspot+1)])
		{
		    $recomcurrent = $recomcurrent + $recomrateFULL[$endspot]*($posvecFULL[($endspot+1)]-$posvecFULL[$endspot]);
		    $endspot = $endspot+1;
		    if ($endspot == ($nsitesFULL-1))
		    {
			last;
		    }
		}
		$recomcurrent = $recomcurrent + $recomrateFULL[$endspot]*($posvec[($i+1)]-$posvecFULL[$endspot]);
		$recomrate[$i] = ($recomcurrent/($posvec[($i+1)]-$posvec[$i]))/$Mb;  
	    }
	}
    }
    ## print:
    open(OUT2,">${outfilePRE}.recomrates");
    print OUT2 "start.pos recom.rate.perbp\n";
    for ($i=0; $i < ($nsites-1); $i+=1)
    {
	$recomI = sprintf("%.15f",$recomrate[$i]/100);      ## as should be a prob
	print OUT2 "$posvec[$i] $recomI\n";
    }
    print OUT2 "$posvec[($nsites-1)] 0\n";
}
