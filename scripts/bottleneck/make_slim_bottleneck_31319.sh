# Script to make SLIM job script
# USAGE: ./make_slim_bott_31319.sh [Na] [Nb] [nF] [h]


# Set Na, the ancestral population size
Na=${1}

# Set Nb, the bottleneck population size
Nb=${2}

# Set nF, the bottleneck population size
nF=${3}

# Set h, dominance coefficient
h=${4}

# Make script
cat > slim_bottleneck_${Na}Na_${Nb}Nb_${nF}nF_h${h}_31319.slim << EOM

initialize() {
	
	initializeSLiMModelType("nonWF");
	defineConstant("K1", ${Na});
	defineConstant("K3", ${Nb});
	defineConstant("num_founders", ${nF});
	defineConstant("phi", 0.9); //parameter for OU model used for stochastic demography
	defineConstant("sampleSize", 60);
	defineConstant("g",20000); //number of genes
	defineConstant("ROHcutoff", 1000000);
	defineConstant("geneLength", 1500);
	defineConstant("seqLength", g*geneLength);
	//cat("Genome length:"+seqLength+"\n");	
	
	initializeMutationRate(1e-8);
	initializeMutationType("m1", 0.0, "g",-0.01314833, 0.186);
	initializeMutationType("m2", 0.5, "f", 0.0);
	
	initializeGenomicElementType("g1", c(m1,m2), c(2.31,1.0));
	
	
	// approach for setting up genes on different chromosomes adopted from Jacqueline's wolf scripts 
	
	gene_nums=c(56,39,42,40,40,35,37,34,28,31,34,33,29,28,29,27,29,25,24,26,23,28,24,21,23,18,21,19,19,18,18,18,14,19,12,14,14,11);
	gene_nums = gene_nums*g/1000; //need to scale to number of desired genes since above array was originally set up for 1000 genes
	
	
	for (i in 1:g){
		initializeGenomicElement(g1, ((i-1)*geneLength)+(i-1), (i*geneLength)+(i-2) );
	}
	
	// # of genes per chromosome for 1000 genes:
	
	
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


reproduction() {
	subpop.addCrossed(individual, subpop.sampleIndividuals(1));
}


1 early() {
	cat("gen,K,popSize,meanFitness,avgVstrDel,avgStrDel,avgModDel,avgWkDel" + "\n");
	sim.addSubpop("p1", 10);

}
1:$((${Na}*10)) early() {
	p1.fitnessScaling = K1 / p1.individualCount;
}

//track statistics pre-bottleneck every 500 generations
1:$((${Na}*10)) late() {
	if (sim.generation % 500 == 0) {
		
		// Sample scale statistics
		i = sample(p1.individuals, sampleSize, F);
		
		Num_VstrDel_muts = c();
		Num_strDel_muts = c();
		Num_modDel_muts = c();
		Num_wkDel_muts = c();
		fitness_population = c();
		
		for (individual in i) {
		
			// get selection coefficients for all muts and sum avg num per individual for different types (weak, moderate, strong, very strong)
			s = individual.genomes.mutations.selectionCoeff;
			
			Num_VstrDel_muts = c(Num_VstrDel_muts, sum(s<=-0.05));
			Num_strDel_muts = c(Num_strDel_muts, sum(s<=-0.01));
			Num_modDel_muts = c(Num_modDel_muts, sum(s<=-0.001 & s > -0.01));
			Num_wkDel_muts = c(Num_wkDel_muts, sum(s<=-0.00001 & s > -0.001));
			
			
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
		
		cat(sim.generation + "," + K1 + "," + p1.individuals.size() + "," + mean(fitness_population) + "," + sum(Num_VstrDel_muts)/i.size() + "," + sum(Num_strDel_muts)/i.size() + "," + sum(Num_modDel_muts)/i.size() + "," + sum(Num_wkDel_muts)/i.size() + ","+ "\n");
	
	}
}


// bottleneck to p3
$((${Na}*10+1)) early(){
	sim.addSubpop("p3",0);
	migrants = sample(p1.individuals, num_founders);
	p3.takeMigrants(migrants);
	cat("gen,K3,p_death,popSizeP3,meanFitness,meanHet,FROH,meanDelMut,avgVStrDel,avgStrDel,avgModDel,avgWkDel," + "\n");
	sim.tag = K3;
}




// fitness scaling for p3

$((${Na}*10+1)):$((${Na}*10+5000)) early() {
	p1.fitnessScaling = 0;
	
	
	
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
	
	
	sim.tag = asInteger(exp((1-phi)*log(K3)+phi*log(sim.tag)+rnorm(n = 1, mean = 0, sd = log10(1.3))));
	p3.fitnessScaling = sim.tag / p3.individualCount;
	
	cat(sim.generation + "," + sim.tag + "," + beta + ",");

}



// track statistics for P3 every generation
$((${Na}*10+1)):$((${Na}*10+5000)) late() {
	if (p3.individuals.size() > 1){
		if(p3.individuals.size() > 60){
			i = sample(p3.individuals, sampleSize, F);
		}
		else{
			i=p3.individuals;
		}
		
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
		
		ind_heteroz = c();
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
				ind_heteroz = c(ind_heteroz, 0); //putting this here to avoid error when trying to sum null vector
				next;
			}
			
			ind_heteroz = c(ind_heteroz, sum(genotype==1)/(seqLength));
			
			
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
		
		cat(p3.individuals.size() + "," + mean(fitness_population) + "," + mean(ind_heteroz) + "," + mean(ROH_length_sumPerInd)/seqLength + ","  +  mean(p3.individuals.sumOfMutationsOfType(m1)) + "," + sum(Num_VstrDel_muts)/i.size() + "," + sum(Num_strDel_muts)/i.size() + "," + sum(Num_modDel_muts)/i.size() + "," + sum(Num_wkDel_muts)/i.size() + ","+ "\n");
	
	}
	
	if(p3.individuals.size() < 2){
		sim.simulationFinished();
		cat("The population has gone extinct" + "\n");
	}
}
EOM