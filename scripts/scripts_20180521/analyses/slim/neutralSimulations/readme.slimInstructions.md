#readme for slim --> dadi

## order of events

### 1 slim replicate:

An array job (slim.1D.2Epoch.array.sh, with modifiable parameters for the simulation) is submitted using slim.1D.2Epoch.submitter.sh. 

Each task in the array simulates a chunk that is made up of 100x1kb independent blocks 

(1e-08 recomb rate within block; 0.5 between blocks)

so there is a total of 6Mb (6000x1kb blocks) simulated


##### current slim model: based on the results of dadi, the population has a long-term Ne of 4000, followed by a rapid crash to 30 which lasts for 10 generations, and then 7 individuals are sampled for the VCF file
##### only neutral mutations are simulated. mu is from otter genome paper (8.64e-9)

### 100 slim replicates:
The above process is repeated so that there are 100 replicates of the 6Mb, yielding 100 different SFSes (one for each replicate), each based on 6Mb of neutral sequence


### Concatenation 
After all chunks of all replicates finish, you run slim.1D.2Epoch.concatVCF.generateSFS.sh to concatenate results into one VCF (with artificial chromosome numbers to signify each chunk) and generate an SFS based on the concatenated VCF 


### Dadi inference
Then dadi inference is carried out on each of those 100 SFSes, with 50 replicates per SFS. The results are found in the replicate's directory in the slim/ dir. 