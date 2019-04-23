# Script to make SLIM job script
# USAGE: ./make_slim_wolf_101517.job.sh [g] [t] [h] [j]

# Set g, number of genes
g=${1}

# Set t, number of burn-in generations
t=${2}

# Set h, dominance coefficient
h=${3}

# Set j, the job number
j=${4}

# Make script
cat > slim_wolf_${g}genes_${t}burn_h${h}_101517.job << EOM

// Wolf Sim
// slim_wolf_${g}genes_${t}burn_h${h}_101517.job

initialize() {
initializeMutationRate(mu); 
initializeMutationType("m1", ${h}, "g", -0.01314833, 0.186); 
initializeMutationType("m2", 0.5, "f", 0.0);

initializeGenomicElementType("g1", c(m1,m2), c(2.31,1.0)); 

for (i in 1:${g}){
	initializeGenomicElement(g1, ((i-1)*1000)+(i-1), (i*1000)+(i-2) );
}
// AB: don't want to use this --> because sea otters have 19 chr pairs
// what to do instead? how bad is it to make it uniform? can I split them up if I'm not worried about extinction modeling?
// see if my existing simulations work or make any sense maybe? 
// # of genes per chromosome for 1000 genes:
gene_nums=c(56,39,42,40,40,35,37,34,28,31,34,33,29,28,29,27,29,25,24,26,23,28,24,21,23,18,21,19,19,18,18,18,14,19,12,14,14,11);

rates=NULL;

// Multiple chromosomes:
for (i in 1:(size(gene_nums)-1)){
	rates=c(rates, 0, rep(c(1e-3, 0), gene_nums[i-1]-1), 0.5);
}
rates=c(rates, 0, rep(c(1e-3, 0), gene_nums[size(gene_nums)-1]-1));

ends=NULL;
for (i in 1:${g}){
	ends=c(ends, (i*1000)+(i-2), (i*1000)+(i-1));
}
ends=ends[0:(size(ends)-2)];

initializeRecombinationRate(rates, ends);

}

1 { 
	sim.addSubpop("p1", 45000); 
} 

${t} {
	sim.addSubpopSplit("p2", 2500, p1);
	p1.setSubpopulationSize(17350);
}

1:$((${t} + 4167)) late() {
	if (sim.generation % 1000 == 0){
		writeFile("/u/scratch/b/bkim331/jacqueline/wolf/slim_wolf_${g}genes_${t}burn_h${h}_sim${j}.gen", paste(sim.generation));
	}
}

$((${t} + 10000)) late() {
  //output population sizes for easy reference
  cat("#p1: " + size(p1.individuals) + " individuals; p2: " + size(p2.individuals) + " individuals.\n");
  
  //file header
  //mutation id
  //mutation type
  //selection coefficient
  //age of mutation in generations
  //subpopulation it arose in
  //number of heterozygote derived in p1
  //number of homozygote derived in p1
  //number of heterozygote derived in p2
  //number of homozygote derived in p2
  //these are genotype counts not allele counts
  cat("mutid,type,s,age,subpop,p1numhet,p1numhom,p2numhet,p2numhom\n");
  
  //for every mutation in the simulation
  for (mut in sim.mutations){
    id = mut.id;
    s = mut.selectionCoeff;
    subpop = mut.subpopID;
    age = sim.generation - mut.originGeneration;
    type = mut.mutationType;
    
    //initialize genotype counts
    p1numhet = 0;
    p1numhom = 0;
    p2numhet = 0;
    p2numhom = 0;
    
    //count hom and het derived in p1
    for (p1i in p1.individuals){
      gt = sum(c(p1i.genomes[1].containsMutations(mut), p1i.genomes[0].containsMutations(mut)));
      if (gt == 1){
        p1numhet = p1numhet + 1;
      } else if (gt == 2){
        p1numhom = p1numhom + 1;
      }
    }
    
    //count hom and het derived in p2
    for (p2i in p2.individuals){
      gt = sum(c(p2i.genomes[1].containsMutations(mut), p2i.genomes[0].containsMutations(mut)));
      if (gt == 1){
        p2numhet = p2numhet + 1;
      } else if (gt == 2){
        p2numhom = p2numhom + 1;
      }
    }
    
    // string for mutation type. add m3, m4, etc. if you have multiple types
    if (type == m1){
      type = "m1";
    } else if (type == m2){
      type = "m2";
    }
    
    //print results
    cat(id + ",");
    cat(type);
    cat("," + s + "," + age + "," + subpop + "," + p1numhet + "," + p1numhom + "," + p2numhet + "," + p2numhom + "\n");
  }
}

EOM
