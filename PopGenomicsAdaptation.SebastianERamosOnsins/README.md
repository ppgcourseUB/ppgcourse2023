# Population Genomics to Adaptation
## The Detection of the Proportion of Beneficial Substitutions

Instructor: **Sebastian E. Ramos-Onsins**

In this practical session we will focus on the effect of positive and negative selection on the polymorphism and the divergence at functional and neutral positions across the genome. We will simulate a number of different evolutionary scenarios and analyze the different observed patterns.

The objectives of this practical session are:

* To familiarize with the forward-simulator *Slim*. Understanding how it works and be able to construct your own code.
* Run different scenarios containing positive and/or negative selection, plus some demographic effects that can affect the detection of beneficial mutations.
* Estimate the proportion of beneficial substitutions using different methods, and contrast with the real value per scenario.

## Forward Simulator: *Slim*
1.	We will simulate a number of different scenarios using [*Slim*](https://messerlab.org/slim/) (Messer, Genetics 2013, Haller and Messer, MBE 2017, [Haller and Messer, MBE 2019](https://academic.oup.com/mbe/article/36/3/632/5229931?login=true)). Slim is a forward simulator that allows to simulate many selective positions at the same time in complex demographic patterns. Slim has a **graphical interface** (we will see an example on the practical class) but to speed up the simulation analysis we will use the **command line** program in the practical session. 
2. The simulator contains an extended manual and many "recipes" or examples for many different user, including complex metapopulation models in spatial context, non-Wright-Fischer models and also allows the use of phenotypic traits and QTLs in relation to the genotypes. 
3. The simulator is designed in a very versatile way, which allows to include new functions, debugging code and controling variable using *Eidos*. *Slim* also allows to output trees, genomes, mutations and substitutions at any step in the simulation.

If your computer allows it, download [*Slim*](https://messerlab.org/slim/). Start the graphical interface application and open the initial recipe 4.1 - A basic neutral simulation. Understand the basic commands included. Use the Simulation panel to "Dump Population State" and see the output of all genomes, mutations and substitutions at the current population.

```
// set up a simple neutral simulation
initialize()
{
	// set the overall mutation rate
	initializeMutationRate(1e-7);
	
	// m1 mutation type: neutral
	initializeMutationType("m1", 0.5, "f", 0.0);
	
	// g1 genomic element type: uses m1 for all mutations
	initializeGenomicElementType("g1", m1, 1.0);
	
	// uniform chromosome of length 100 kb
	initializeGenomicElement(g1, 0, 99999);
	
	// uniform recombination along the chromosome
	initializeRecombinationRate(1e-8);
}
```
The simulation needs to initialize some parameters. Here it is defined the mutation and the recombination rate per position, the type of mutation used (named "m1", with dominance 0.5 and a fixed selective value of s=0 (neutral)). The mutation "m1" is included in a genomic element named "g1", which contains 1e5 positions.

```
// create a population of 500 individuals
1
{
	sim.addSubpop("p1", 500);
}

// run to generation 10000
10000
{
	sim.simulationFinished();
}
```
The simulation starts at generation 1, creating a population named "p1" with 500 diploid individuals. Finally, the simulation finish at generation 10000, following a number of steps (which can also be modified) that includes migration (if defined), generation of offspring by choosing parents based on fitness, mutation, recombination and recalculation of fitness for the new individual.

### slim_template.slim ###
Open the file located in this Github folder named "slim_template.slim", using the graphical interface or a Text Editor to visualize. Explanation of the commands included in the file.

```
// command line including necessary parameters:
// slim -t -m -d "Ne=500" -d "L=30000" -d "Neb=500" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=2" -d "rate_ben=0.1,1" -d "rate_del=8" -d "s_backg_ben=0.005" -d "s_backg_del=-1" -d "nsweeps=1" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.1" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='~/slim.test_file.txt'" ./slim_template.slim


```
slim can run in command line. Here is one example where a number of parameters are externally defined by the user at the command.

```
// set up a simple neutral simulation
initialize() {

	if (exists("slimgui")) {
		defineConstant("Ne",500); //number of diploid infdividuals
		defineConstant("L",500000); //length of the genome		
		defineConstant("Neb",500); //number of diploid infdividuals given sudden change pop size
		
		defineConstant("mut_rate",1e-6);	
		defineConstant("rec_rate",1e-4);
		defineConstant("ngenes",100); //number of independent genes

		defineConstant("rate_ben",0.5); //rate beneficial mutations versus neutral (1) and deleterious
		defineConstant("rate_del",8); //rate deleterious mutations versus neutral (1) annd benficial
		defineConstant("s_backg_ben",+0.005);	//s of beneficial background 
		defineConstant("s_backg_del",-0.05);	//s of deleteriious background 

		defineConstant("nsweeps",0); //if >0, define number of sweeps frequencies to start and end plus strength 
		defineConstant("freq_sel_init",0.05);
		defineConstant("freq_sel_end",0.5);
		defineConstant("s_beneficial",0.1);
		
		defineConstant("ind_sample_size",25); //number of samples to use
		defineConstant("out_sample_size",1); //number of outgroup samples to use
		defineConstant("file_output1","~/slim.test_file.txt"); //name of output file
	}

	//define a fixed demographic process plus a possible selective sweep before sampling
	defineConstant("tsplit",5*Ne); //time split outgroup
	defineConstant("tne",14*Ne); //time change Ne
	defineConstant("tsweep",asInteger(14*Ne + Ne/2)); //time sweep
	defineConstant("tend",15*Ne); //time end

	initializeSLiMOptions(nucleotideBased=T);
	initializeAncestralNucleotides(randomNucleotides(L));
	
	//separate each gene using recombination 0.5
	rrates = NULL;
	ends = NULL;
	len = asInteger(L/ngenes-1);
	for(i in 1:ngenes) {
		rrates = c(rrates,rec_rate,0.5);
		ends = c(ends,len,len+1);
		len = len + asInteger(L/ngenes); 
	}	
	rrates = rrates[0:(2*ngenes-2)];
	ends = ends[0:(2*ngenes-2)];
	initializeRecombinationRate(rrates,ends); 

	// m1 mutation type: (neutral)
	initializeMutationTypeNuc("m1", 0.5, "f", 0.0);
	// m2 mutation type: (deleterious)
	if(s_backg_del==0) {initializeMutationTypeNuc("m2", 1.0, "f", s_backg_del);}
	else { initializeMutationTypeNuc("m2", 1.0, "g", s_backg_del, 0.2);}
	// m3 mutation type: (beneficial)
	initializeMutationTypeNuc("m3", 0.5, "f", s_backg_ben);

	// g1 genomic element type: (synonymous) uses m1 for all mutations
	initializeGenomicElementType("g1", m1, 1.0,mmJukesCantor(mut_rate/3));
	// g2 genomic element type: (nonsynonymous) uses all mutations
	initializeGenomicElementType("g2", c(m1,m2,m3), c(1,rate_del,rate_ben),mmJukesCantor(mut_rate/3));

	//the chromosome is made with L/3 codons having two nonsyn + 1 syn position
	types = rep(c(g2,g2,g1), asInteger(L/3));
	if(L%3==1) {types=c(types,g2);}
	if(L%3==2) {types=c(types,g2,g2);}
	position = seq(0,L-1);
	initializeGenomicElement(types, position, position); //each codon contains g2,g2,g1
}
```
The first step (initialize) is here more complex. Include the definition of parameters (*e.g.,* Ne, ngenes, selective effect ad proportion of beneficial and deleterious mutations, or the sample size collected at the end) in case the code is used with the graphical interface. Then, a number of parameters related to the demographic scenario are also included. In this example, we use the option to simulate nucleotide sequences (A,C,G,T) instead of only mutations. The recombination rate is defined within and between genes (being free recombination between them). 

It is very important to understand the definition of two different genomic elements: g1 (synonymous) and g2 (nonsynonymous). g1 contains only "m1" (neutral) mutations, while g2 contains "m1" (neutral), "m2" (deleterious) and "m3" (beneficial) in defined proportions for each. g1 and g2 are defined as (g2, g2, g1) per codon and consecutively in the whole genome.

```
//START SIMULATION
1{
	sim.addSubpop("p1", Ne); //create initial population and run up to equilibrium for 10Ne gen.
	sim.rescheduleScriptBlock(s1, start=tsplit, end=tsplit); //split OUTG
	sim.rescheduleScriptBlock(s2, start=tne, end=tne); //Sudden change Ne
	sim.rescheduleScriptBlock(s3, start=tsweep, end=tsweep); //possible sweeep
	sim.rescheduleScriptBlock(s4, start=tend, end=tend); //end simulation. Print data
}

// Split de p1 in the generation 5000
s1 50000 { sim.addSubpopSplit("p2", Ne, p1); }

// Sudden Change in population size
s2 97500 { p2.setSubpopulationSize(Neb); }

//if required, force a strong selective sweep in n different positions. Can be successful or not
s3 98000 {	
	if(nsweeps) {
		muts = sim.mutations;
		muts = muts[muts.mutationType==m3]; //added: look only at positive mutations
		mutsp2b = sim.mutationFrequencies(p2, muts);
		muts = muts[mutsp2b >= freq_sel_init & mutsp2b <= freq_sel_init + (freq_sel_end - freq_sel_init)/5];	
		if (size(muts)>=nsweeps)
		{
			mut = sample(muts, nsweeps);
			mut.setSelectionCoeff(s_beneficial);
			print("Sweep at positions: ");
			print(mut.position);
		}
		else {
			cat("No contender of sufficient frequency found.\n");
		}
	}
}

```
Here, the simulation is separated in blocks from the generation 1 (s1, s2,s3,s4). These blocks can be defined at the desired time by user. The blocks indicate a split of ancestral population in two populations, outgroup and target (s1), the change in Ne at the target population (s2) and the possibility to have strong selective sweeps (s3).

```

// Run up 10000 generations
s4 100000 late() {
	//OUTPUT
	// Select no samples from the outgroup and ni samples of the target population and output
	// obtain random samples of genomes from the three subpopulations
	g_1 = sample(p2.genomes,2*ind_sample_size,F);
	g_2 = sample(p1.genomes,2*out_sample_size,F);
	//Concatenate the 2 samples
	g_12=c(g_1,g_2);
	//Get the unique mutations in the sample, sorted by position
	m = sortBy(unique(g_12.mutations),"position");
	
	//separate mutations of element g2 or/and in codon position 1 and 2 (nonsyn) 
	nonsyn_m = m[(m.position+1) % 3 != 0];
	//separate mutations of element g1 or/and in codon position 3 (syn)
	syn_m = m[(m.position+1) % 3 == 0];
	//look for polymorphisms and substitutions in target population that are fixed in outgroup
	sfs = array(rep(0,2*(2*ind_sample_size+3)),dim=c(2,2*ind_sample_size+3));
	//SFS for nonsyn
	for(position in nonsyn_m.position) {
		fr = 0; //frequency of each variant
		for(genome in g_1) {
			fr = fr + sum(match(position,genome.mutations.position) >= 0);
		}
		//print("\nposition="+position);
		//print("fr="+fr);
		fro = 0;
		for(genome in g_2) {
			fro = fro + sum(match(position,genome.mutations.position) >= 0);
		}
		//print("fro="+fro);
		if(fr>0 & fro==0) {
			sfs[0,fr-1] = sfs[0,fr-1] + 1;
		}
		if(fr==0 & fro==2*out_sample_size) {
			sfs[0,2*ind_sample_size-1] = sfs[0,2*ind_sample_size-1] + 1;
			//number of beneficial substitutions (m3)
			if(nonsyn_m.mutationType[which(nonsyn_m.position==position)[0]]==m3) { 
				sfs[0,2*ind_sample_size+2]	= sfs[0,2*ind_sample_size+2] + 1; 
			}
		}
	}
	//add the mutations from added selective sweeps (assumed fixed)
	//sfs[0,2*ind_sample_size+2]	= sfs[0,2*ind_sample_size+2] + nsweeps;
	
	//SFS for nsyn
	for(position in syn_m.position) {
		fr = 0;
		for(genome in g_1) {
			fr = fr + sum(match(position,genome.mutations.position) >= 0);
		}
		fro = 0;
		for(genome in g_2) {
			fro = fro + sum(match(position,genome.mutations.position) >= 0);
		}
		if(fr>0 & fro==0) {
			sfs[1,fr-1] = sfs[1,fr-1] + 1;
		}
		if(fr==0 & fro==2*out_sample_size) {
			sfs[1,2*ind_sample_size-1] = sfs[1,2*ind_sample_size-1] + 1;
		}
	}
	//include length of sequences and move the divegence like in polyDFE output
	sfs[0,2*ind_sample_size+0] = sfs[0,2*ind_sample_size-1];
	sfs[0,2*ind_sample_size-1] = asInteger(L*2/3);
	sfs[0,2*ind_sample_size+1] = asInteger(L*2/3);
	sfs[1,2*ind_sample_size+0] = sfs[1,2*ind_sample_size-1];
	sfs[1,2*ind_sample_size-1] = asInteger(L*1/3);
	sfs[1,2*ind_sample_size+1] = asInteger(L*1/3);
	
	//manipulate for easy print of sfs matrix
	print_header = "SFS";
	print_header = print_header + paste("	fr"+c(1:(2*ind_sample_size-1)));
	print_header = print_header + "	" + "PosP" + "	" + "Fixed" + "	" + "PosF" + "	" + "FixBen";
	print_sfs_nsyn = "nsyn" + "	";
	print_sfs_syn  = "syn" + "	";
	for(i in c(0:(2*ind_sample_size+2))) {
		print_sfs_nsyn = print_sfs_nsyn + sfs[0,i] + "	"; 
		print_sfs_syn  = print_sfs_syn  + sfs[1,i] + "	"; 
	}
	
	//OUTPUT:
	//for syn and for nonsyn:
	// print the sfs of target, the substitutions vs outgroup and the total sites
	print("Saving results to file " + file_output1);
	writeFile(filePath=file_output1,contents=(print_header),append=F);
	writeFile(filePath=file_output1,contents=(print_sfs_nsyn),append=T);
	writeFile(filePath=file_output1,contents=(print_sfs_syn),append=T);
				
	print("Simulation finished");
}

```
The step s4 finish the simulation and collects the Site Frequency Spectrum (SFS) for a sample in synonymous and nonsynonymous positions, plus the number of fixations and the number of true beneficial fixations.

##Run Simulations under Different Selective Scenarios

Here, we will simulate different selective scenarios in order to evaluate the ability to detect the proportion of beneficial substitutions under the defined conditions. Nine different scenarios are defined in the script "**run\_construct\_slim\_conditions.sh**", although the user can modify the conditions if desired (but be careful to not include unrealistic or never-finish conditions!). 

Please modify the name of the job to identify yours.

```
#header for run in slurm
echo #!/bin/bash > ./run_slim_conditions.sh
echo #SBATCH -p normal > ./run_slim_conditions.sh
echo #SBATCH --job-name-USER > ./run_slim_conditions.sh 
echo #SBATCH -o %j.out > ./run_slim_conditions.sh
echo #SBATCH -e %j.err > ./run_slim_conditions.sh
echo module load slim > ./run_slim_conditions.sh

#fixed paraneters
Ne=500; L=500000; ngenes=100;
mut_rate=1e-6;
ind_sample_size=25; out_sample_size=1;

# CONDITION 0:
#Neutral. No change Ne.
FILEOUT="'./00_slim_SFS_SNM.txt'"
Neb=500; nsweeps=0;
rec_rate=1e-4;
rate_ben=0; s_backg_ben=0;
rate_del=0; s_backg_del=0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim > ./run_slim_conditions.sh

# CONDITION 1:
#Strong BACKGROUND SELECTION. No beneficial selection. No change Ne. No sweep
FILEOUT="'./01_slim_SFS_BGS.txt'"
Neb=500; nsweeps=0;
rec_rate=1e-4
rate_ben=0; s_backg_ben=0;
rate_del=8; s_backg_del=-0.1;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 2:
#No background selection. BENEFICIAL SELECTION. No change Ne. No sweep
FILEOUT="'./02_slim_SFS_WBGS.txt'"
Neb=500; nsweeps=0;
rec_rate=1e-4
rate_ben=0.005; s_backg_ben=0.005;
rate_del=0; s_backg_del=0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 3:
#Strong BACKGROUND SELECTION. BENEFICIAL SELECTION. POPULATION REDUCTION. No sweep.
FILEOUT="'./03_slim_SFS_BGS_RED.txt'"
Neb=200; nsweeps=0;
rec_rate=1e-4
rate_ben=0.05; s_backg_ben=0;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 4:
#Strong BACKGROUND SELECTION. BENEFICIAL SELECTION. POPULATION EXPANSION. No sweep.
FILEOUT="'./04_slim_SFS_BGS_EXP.txt'"
Neb=2500; nsweeps=0;
rec_rate=1e-4
rate_ben=0.05; s_backg_ben=0;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 5:
#Strong BACKGROUND SELECTION. SMALL PROPORTION BENEFICIAL SELECTION. No change Ne. No sweep.
FILEOUT="'./05_slim_SFS_BGS_PSL.txt'"
Neb=500; nsweeps=0;
rec_rate=1e-4
rate_ben=0.05; s_backg_ben=0.005;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 6:
#Strong BACKGROUND SELECTION. MIDDLE PROPORTION BENEFICIAL SELECTION. No change Ne. No sweep.
FILEOUT="'./06_slim_SFS_BGS_PSM.txt'"
Neb=500; nsweeps=0;
rec_rate=1e-4
rate_ben=0.5; s_backg_ben=0.005;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 7:
#Strong BACKGROUND SELECTION. HIGH PROPORTION BENEFICIAL SELECTION. No change Ne. No sweep.
FILEOUT="'./07_slim_SFS_BGS_PSH.txt'"
Neb=500; nsweeps=0;
rec_rate=1e-4
rate_ben=2; s_backg_ben=0.005;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 8:
#Strong BACKGROUND SELECTION. BENEFICIAL SELECTION. No change Ne. SWEEPS.
FILEOUT="'./08_slim_SFS_BGS_PSM_SW.txt'"
Neb=500; nsweeps=10;
rec_rate=1e-4
rate_ben=0.5; s_backg_ben=0.005;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

```
The script run in shell directly:

```
sh  ./run\_construct\_slim\_conditions.sh
```
This command creates a new script named "*run_slim_conditions.sh*".This file will be run in the cluster using slurm.
 
```
sbatch ./run\_construct\_slim\_conditions.sh
```
and check the progress with:

```
squeue
```
The results will be separated in nine different files, one for each simulation condition. The files will start with the name "[number]\_slim\_SFS\_\*.txt". You can see using the command *more*, *less*, *cat* or using a text editor such as *nano*.

The results will be the SFS of nonsynonymous and synonymous sites. The output will be similar to this:

```
SFS	fr1 	fr2 	fr3 	fr4 	fr5 	fr6 	fr7 	fr8 	fr9 	fr10 	fr11 	fr12 	fr13 	fr14 	fr15 	fr16 	fr17 	fr18 	fr19 	fr20 	fr21 	fr22 	fr23 	fr24 	fr25 	fr26 	fr27 	fr28 	fr29 	fr30 	fr31 	fr32 	fr33 	fr34 	fr35 	fr36 	fr37 	fr38 	fr39 	fr40 	fr41 	fr42 	fr43 	fr44 	fr45 	fr46 	fr47 	fr48 	fr49	PosP	Fixed	PosF	FixBen
nsyn	321	152	88	72	40	40	28	25	26	32	17	22	23	21	19	12	13	12	16	7	10	19	10	14	9	13	12	14	11	15	6	6	12	6	16	11	10	8	10	12	8	13	11	7	11	12	9	8	8	333333	2743	333333	1064
syn	263	116	79	78	71	35	32	25	42	23	21	20	12	19	15	18	16	11	16	20	12	15	19	7	9	17	6	11	10	7	9	8	10	10	8	2	7	12	12	6	9	4	8	7	6	6	5	5	6	166666	1552	166666	0	
```

##Estimating the proportion of Beneficial Substitutions (alpha) from Simulation Data




```
module load r-mass r-proto
```
