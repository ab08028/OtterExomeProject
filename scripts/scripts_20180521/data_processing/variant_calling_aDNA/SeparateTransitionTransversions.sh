# Remove transitions from beagle format

# Transitions
# A <-> G
# C <-> T

# Ref  = A , alt = G
# Ref = G, alt = A
# Ref = C, alt = T
# Ref = T, alt = C
# Beagle codes: the allele codes as 0=A, 1=C, 2=G, 3=T

# and these appear in columns 2 and 3 of the beagle file
# So if there's 0 <-> 2 and 1 <-> 3

angsdDate=20190503
wd=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/angsd-GLs/$angsdDate

for ref in Elut Mfur
do
for state in 1e-6.snpsOnly 1e-6.snpsOnly.downSampOnly.minInd.5 # skipping allsites 
do
input=$wd/angsdOut.mappedTo${ref}.${state}.beagle.gz
 
# Want to exclude  0 <-> 2 and 1 <-> 3
# So want to *KEEP* 
# 0 <-> 1
# 0 <-> 3
# 1 <-> 2
# 2 <-> 3
 
# 8 possible transversion genotypes
# 0-1 : A-C ; 
# 1-0 : C-A ;
# 0-3 : A-T ;
# 3-0: T-A ;
# 1-2 : C-G ;
# 2-1 : G-C ;
# 2-3 : G-T ;
# 3-2: T-G

#### pull transversions (leave transitions behind)
zcat $input | awk '{

if(($2==0 && $3==1) || ($2==1 && $3==0) || ($2==0 && $3==3) || ($2==3 && $3==0) || ($2==1 && $3==2) || ($2==2 && $3==1) || ($2==2 && $3==3) || ($2==3 && $3==2)) print}' | gzip -f > ${input%.beagle.gz}.transversionsOnly.beagle.gz
# 128662 transversions 
# checking:
#zcat $input | grep -c -v marker 
# 379972 sites in the 1e-6.snpsOnly

# get transitions

zcat $input | awk '{

if(($2==0 && $3==2) || ($2==2 && $3==0) || ($2==1 && $3==3) || ($2==3 && $3==1)) print}' | gzip -f > ${input%.beagle.gz}.transitionsOnly.beagle.gz

# number of transitions: 251310

# they add up: 128662 (Transversions)+ 251310 (transitions) = 379972 total 
done
done


