# Script to make SLIM job script
# USAGE: ./make_slim_bott_4319.sh [Na] [Nb] [nF] [h]


# Set Na, the ancestral population size
Na=${1}

# Set Nb, the bottleneck population size
Nb=${2}

# Set nF, the number of founders
nF=${3}

# Set h, dominance coefficient
h=${4}

# Make script
cat > slim_otter_bottleneck_ageStruc_${Na}Na_${Nb}Nb_${nF}nF_h${h}_82319.slim << EOM

initialize() {
	
	initializeSLiMModelType("nonWF");
	defineConstant("K1", ${Na});
	defineConstant("K3", ${Nb});
	defineConstant("num_founders", ${nF});
	defineConstant("sampleSize", 30);
	defineConstant("g",20000); //number of genes
	defineConstant("ROHcutoff", 1000000);
	defineConstant("geneLength", 1500);
	defineConstant("seqLength", g*geneLength);
	//cat("Genome length:"+seqLength+"\n");
	initializeSex("A");
	
	
	initializeMutationRate(8.64e-9);
	initializeMutationType("m1", ${h}, "g",-0.01314833, 0.186);
	initializeMutationType("m2", 0.5, "f", 0.0);
	defineConstant("L", c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0));
	
	initializeGenomicElementType("g1", c(m1,m2), c(2.31,1.0));
	
	
	// approach for setting up genes on different chromosomes adopted from Jacqueline's wolf scripts 
	
	// vector of # genes on 38 different dog chromosomes (scaled by Jacqueline to 1000 total genes) - need to rescale according to number of desired genes
	gene_nums=c(56,39,42,40,40,35,37,34,28,31,34,33,29,28,29,27,29,25,24,26,23,28,24,21,23,18,21,19,19,18,18,18,14,19,12,14,14,11);
	gene_nums = gene_nums*g/1000; //need to scale to number of desired genes since above array was originally set up for 1000 genes
	
	
	for (i in 1:g){
		initializeGenomicElement(g1, ((i-1)*geneLength)+(i-1), (i*geneLength)+(i-2) );
	}
		
	
	rates=NULL;
	
	// Multiple chromosomes:
	for (i in 1:(size(gene_nums)-1)){
		rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_nums[i-1]-1)), 0.5);
	}
	rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_nums[size(gene_nums)-1]-1)));
	
	ends=NULL;
	for (i in 1:g){
		ends=c(ends, (i*geneLength)+(i-2), (i*geneLength)+(i-1));
	}
	ends=ends[0:(size(ends)-2)];
	
	initializeRecombinationRate(rates, ends);

}



// define function getStats that randomly samples a subpopulation for sampSize # of inds and outputs a string of: 
// pop size, mean fitness, heterozygosity, mean Froh, and avg num of variants of different classes per individual (very str del, str del, mod del, wk del)

function (s) getStats(o pop, i sampSize)
{
	i = sample(pop.individuals, sampSize, F);
	
	m = sortBy(i.genomes.mutations, "position"); //get all mutations in sample
	m_uniq = unique(m); // get rid of redundant muts
	DAF = sapply(m_uniq, "sum(m == applyValue);"); // count number of each mut in pop
	m_uniq_polym = m_uniq[DAF != i.genomes.size()]; //remove fixed mutations??
	
	//initialize vectors
	ROH_length_sumPerInd = c();
	Num_VstrDel_muts = c();
	Num_strDel_muts = c();
	Num_modDel_muts = c();
	Num_wkDel_muts = c();
	ind_het = c();
	fitness_population = c();
	
	for (individual in i) {
		
		indm = sortBy(individual.genomes.mutations, "position");
		indm = indm[match(indm, m_uniq_polym) >= 0];   // Check that individual mutations are not fixed 
		indm_uniq = unique(indm);
		
		genotype = sapply(indm_uniq, "sum(indm == applyValue);");
		
		// tally number of mutations for different classes of selection coefficient per individual
		s = individual.genomes.mutations.selectionCoeff;
		
		Num_VstrDel_muts = c(Num_VstrDel_muts, sum(s<=-0.05));
		Num_strDel_muts = c(Num_strDel_muts, sum(s<=-0.01));
		Num_modDel_muts = c(Num_modDel_muts, sum(s<=-0.001 & s > -0.01));
		Num_wkDel_muts = c(Num_wkDel_muts, sum(s<=-0.00001 & s > -0.001));
		
		if (isNULL(genotype)) {
			ind_het = c(ind_het, 0); //putting this here to avoid error when trying to sum null vector
			next;
		}
		
		ind_het = c(ind_het, sum(genotype==1)/(seqLength));
		
		//code for getting ROHs
		
		ID_het = (genotype == 1); //outputs T/F for genotypes if they are het or homDer
		ID_homDer = (genotype == 2);
		pos_het = indm_uniq.position[ID_het]; //outputs positions of heterozgoys genotypes
			
		startpos = c(0, pos_het); //adds 0 to beggining of vector of hets
		endpos = c(pos_het, sim.chromosome.lastPosition); //adds last position in genome to vector of hets
		pos_het_diff = endpos - startpos;
		ROH_startpos = startpos[pos_het_diff > ROHcutoff]; //filter out startpos that dont correspond to ROH > 1Mb
		ROH_endpos = endpos[pos_het_diff > ROHcutoff];
		ROH_length = pos_het_diff[pos_het_diff > ROHcutoff]; //vector of ROHs for each individual	
		ROH_length_sum = sum(ROH_length);
		ROH_length_sumPerInd = c(ROH_length_sumPerInd, ROH_length_sum); // add sum of ROHs for each individual to vector of ROHs for all individuals
		
		// calculate individual fitness - code from Bernard	
		allmuts = c(individual.genomes[0].mutationsOfType(m1), individual.genomes[1].mutationsOfType(m1));
		uniquemuts = individual.uniqueMutationsOfType(m1);
		
		fitness_individual = c();
		
		if (size(uniquemuts) > 0){
			for (u in uniquemuts){
				places = (allmuts.id == u.id);
				uu = allmuts[places];
				if ((m1.dominanceCoeff == 0.0) & (size(uu) == 2)) {
					fitness = 1 + sum(uu.selectionCoeff)/2;
				} else if ((m1.dominanceCoeff == 0.0) & (size(uu) == 1)) {
					fitness = 1;
				}
				fitness_individual = c(fitness_individual, fitness);
			}
			fitness_individual = product(fitness_individual);
			fitness_population = c(fitness_population, fitness_individual);
		} else {
			fitness_population = c(fitness_population, 1);
		}
	}
	
	return(pop.individuals.size() + "," + mean(fitness_population) + "," + mean(ind_het) + "," + mean(ROH_length_sumPerInd)/seqLength + "," + mean(Num_VstrDel_muts) + "," + mean(Num_strDel_muts)+ "," + mean(Num_modDel_muts) + "," + mean(Num_wkDel_muts));
}

1:$((${Na}*30+1)) reproduction() {
	males = p1.individuals[p1.individuals.sex=='M'];
	females = p1.individuals[p1.individuals.sex=='F'];
	
	old_females = females[females.age > 5];
	old_males = males[males.age > 10];
	
	for(mom in old_females){
		dad = sample(old_males,1);
		p1.addCrossed(mom, dad);
	}
	
	self.active = 0;
	//subpop.addCrossed(individual, subpop.sampleIndividuals(1));

}



1 early() {
	cat("gen,popSize,meanFitness,meanHet,FROH,avgVstrDel,avgStrDel,avgModDel,avgWkDel" + "\n");
	sim.addSubpop("p1", 100);
	p1.individuals.age = rdunif(100, min=0, max=20);

}

1:$((${Na}*30)) early() {
	// life table based individual mortality
	inds = p1.individuals;
	ages = inds.age;
	mortality = L[ages];
	//survival = 1;
	survival = 1 - mortality;
	inds.fitnessScaling = survival;

	p1.fitnessScaling = K1 /( p1.individualCount * mean(survival));
}

//track statistics pre-bottleneck every 1000 generations
1:$((${Na}*30)) late() {
		if (sim.generation % 1000 == 0) {
			stats = getStats(p1, sampleSize);
			cat(sim.generation + "," + stats + "\n");
		}
}



$((${Na}*30+2)):$((${Na}*30+5000)) reproduction() {
        males = p3.individuals[p3.individuals.sex=='M'];
        females = p3.individuals[p3.individuals.sex=='F'];
        
        old_females = females[females.age > 5];
        old_males = males[males.age > 10];
        
        for(mom in old_females){
                dad = sample(old_males,1);
                p3.addCrossed(mom, dad);
        }
        
        self.active = 0;
        //subpop.addCrossed(individual, subpop.sampleIndividuals(1));

}

// bottleneck to p3
$((${Na}*30+1)) early(){
	sim.addSubpop("p3",0);
	migrants = sample(p1.individuals, num_founders);
	p3.takeMigrants(migrants);
	cat("gen,K3,p_death,popSizeP3,meanFitness,meanHet,FROH,avgVStrDel,avgStrDel,avgModDel,avgWkDel" + "\n");
	
	sim.tag = K3; // use sim.tag to keep track of K3 from one generation to the next
}



// fitness scaling for p3

$((${Na}*30+1)):$((${Na}*30+5000)) early() {
	p1.fitnessScaling = 0; // kill off p1


	// life table based individual mortality
	inds = p3.individuals;
	ages = inds.age;
	mortality = L[ages];
	survival = 1 - mortality;
	inds.fitnessScaling = survival;
	
	// kill off individuals at random - not sure if I should then adjust the individualCount
	inds = p3.individuals;
	
	//simulate beta distribution
	alpha = 0.5;
	beta = 8;
	x1 = rgamma(1, mean = alpha, shape=alpha);
	x2 = rgamma(1, mean = beta, shape=beta);
	beta = x1/(x1+x2); //probability of stochastic mortality this generation
	
	
	//set probability of death for each generation equal to outcome of beta 	
	for(i in inds){
		kill = rbinom(1,1,beta);
		if(kill==1){
			i.fitnessScaling = 0.0;
		}
	}
	
	//OU model with phi=0.9
	sim.tag = asInteger(exp((1-0.9)*log(K3)+0.9*log(sim.tag)+rnorm(n = 1, mean = 0, sd = log10(1.3))));
	p3.fitnessScaling = sim.tag / (p3.individualCount * mean(survival));
	
	cat(sim.generation + "," + sim.tag + "," + beta + ",");

}



// track statistics for P3 every generation and terminate when the population goes to 1 individual or after 5000 generations
$((${Na}*30+1)):$((${Na}*30+5000)) late() {
	if(p3.individuals.size() < 2){
		stats_P3 = c("NA,NA,NA,NA,NA,NA,NA,NA"); //cant get stats from just one individual
	}
	if(p3.individuals.size() < sampleSize & p3.individuals.size() > 1){	// case when p3 size is less than sample size but greater than 1
		stats_P3 = getStats(p3, p3.individuals.size());
	}
	if(p3.individuals.size() >= sampleSize){ //case when p3 size is greater than or equal to sample size
		stats_P3 = getStats(p3, sampleSize);
	}
		
	cat(stats_P3 + "\n");	

	if(p3.individuals.size() < 2){
			sim.simulationFinished();
			cat("The population has gone extinct");
	}
}
EOM
