## Data Processing: Fastqs --> Genotypes
### These scripts will take you through initial QC, the paleomix mapping pipeline, variant calling and filtering

#### 0. initialSteps_setupQC: Initial directory setup, fastq download and QC

	a. Step_0_a_checkMD5SUM.sh: check MD5Sum on Fastq files to make sure they downloaded correctly (note: make sure to use "binary" setting in FTP)

	b. Step_0_b_getListsOfSamples.sh: get lists of sample prefixes (referred to as headers) -- these will be what all other scripts iterate through

	c. Step_0_c_FastQC.Hoffman.sh: carry out fastqc (submit with Step_0_c_FastQC.Hoffman.submit.sh)

	d. (not a script) use multiqc to gather fastqc reports

#### 1. paleomixPipeline: maps reads (ancient and modern)

	plmx.0.a plmx_Step_0_a_makeDefaultMakefile.sh: Generate a default paleomix pipeline; modify by hand in text editor to have the settings you want (different for ancient and modern)

	plmx.0.b plmx_Step_0_b_MakeAllMakefiles.sh: Using the templates you made in (a) with you unique settings, make 1 makefile for every sample (diff. anc/mod)

	plmx.1 plmx_Step_1_runPlmx.sh: run paleomix! Submit ancient and modern jobs separately with plmx_Step_1_runPlmx.submit.[ancient/modern].sh (different memory requirements, etc.)

	plmx.2 plmx_Step_2_DownloadMapDamageReports.sh: gather mapDamage plots to transfer to home computer to go into SI materials

##### **Paleomix Notes:** 
	Paleomix Documentation: https://paleomix.readthedocs.io/en/latest/
	Paleomix carries out AdapterRemoval, read mapping, mapping DNA damage and correcting it, and validation.
	It can be used with multiple reference genomes and multiple libraries per sample.
	First, AdapterRemoval trims adapters and low quality bases at the ends of reads. It also collapses overlapping reads.
	For ancient samples, I only use these collapsed reads downstream. Modern samples retain all reads.
	For ancient samples, bwa backtrack with the seed disabled and settings from Kircher's protocol is used to map (-n = 0.01; -o=2; min qual 30)
	For modern samples, bwa mem is used to map (min qual 30)
	All samples have mapDamage plots produced, but the mapDamage is only rescaled for ancient samples.
	Indel realignment is not carried out because GATK Haplotype Caller does it internally

#### 2. variant_calling: find covered variants and call genotypes (GATK)

	a. Step_2_a_FindCoveredIntervals.sh (submit with wrapper Step_2_a_FindCoveredIntervals.[modern/ancient].submit.sh): Detect covered intervals using GATK's FindCoveredIntervals to use downstream (min cov. = 1 read; MAPQ min. 30; min base qual. 20)
	
	b. Step_2_b_QualimapIntervals.sh (submit with wrapper Step_2_b_QualimapIntervals.submit.sh): Run Qualimap on covered regions (use multiqc to gather reports)  

	c. [modern samples only] Step_2_c_HaplotypeCaller.sh (submit with Step_2_c_HaplotypeCaller.submit.sh): call individual variants using GATK's HaplotypeCaller to generate one g.vcf file per sample
	c-alt. [aDNA samples only] -- still determining best way to call genotypes; will likely use pileupcaller --
		
	d. [modern samples only] Step_3_genotypeGVCFS.modern.sh: call genotypes jointly for all populations together (exclude some extremely low coverage individuals)


#### 3.  variant_filtering: filter variants on depth and quality (GATK)

	a. filtering_Step_1_a_filterVariantAndInvariantSites.sh: trim alternate alleles, select biallelic SNPs, filter using GATK HF (QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0), GQ < 20, individual genotype DP > 1000, individual genotype DP < 12, clustered SNPs (3/10)
		
