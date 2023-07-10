# Bayesian inference of species divergences

## 1. Setting the file structure to run `MCMCtree`

Given that we have already generated the calibrated trees, we have the alignment files, and we have just generated the `in.BV` files with the gradient and the Hessian for each alignment... We have everything we need to run `MCMCtree`!

However, we first need to create the corresponding file structure, which will follow the same structure as the directories where the Hessian and the gradient were calculated. You will need to run the following code snippet:

```sh
# Run from `main` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
# Note that, with large datasets, I
# recommend you run at least
# 16 chains and perhaps a maximum of 32 
# chains based on personal experience.
# Nevertheless, things may change 
# depending on your data! 
# This bear dataset is quite small,
# so we will just run 5 chains.
num_chains=5
for j in `seq 1 $num_chains`
do
mkdir -p MCMCtree/$j/{part,conc}/{fosscal,seccal}/{GBM,ILN}
done
mkdir -p pipelines_MCMCtree/{part,conc}/{fosscal,seccal}/{GBM,ILN}
```

----

>**IMPORTANT NOTE**: When sampling from the posterior, the likelihood is being calculated or approximated, depending on the `userdata` option you set in the control file to run `MCMCtree`. In other words, the larger the dataset, the more time it will take for `MCMCtree` to finish. For this bear dataset, five chains are more than enough (we could have even had enough samples with 2 chains, which is the minimum number of chains you must run to check for chain convergence). Nevertheless, sometimes you might find that five chains may not be enough once you run the MCMC diagnostics (we will see that later). If that's the case, you may want to increase the number of chains you are to run. However, that's a decision that you need to take based on the following:
>
> * **Types of computational resources you can use**: can you use as many nodes and partitions as you like (always within reason, please be mindful of your carbon footprint and the amount of CO2 levels you will release due to your analyses!) in your HPCs or are you sharing and need to be cautious of not sending too many jobs so that other people can also run their analyses?
> * **Type of analysis you are running**: are you just exploring your data and want a rough result or have you already done that and want results to add in your manuscripts?
> * **Type of data**: which kind of dataset are you working with? Are there just 3 or 5 genes for 10 taxa or are you working with a huge phylogenomic alignment with 6K taxa and billions of base pairs?
> You might need to find a "sweet spot" for the number of chains you are to run when sampling from the posterior in case you have computational limitations given your resources. You can always start by running fewer chains and, if you need more samples because most of your chains have not passed filters or seem to not have converged or the ESS is too low, you can always send more jobs to your HPCs to run more chains.

----

Now, we will use the same approach we used to generate the pipelines to run `BASEML`. In this case, however, I have modified the bash scripts accordingly. You will find them in the [`scripts` directory](scripts), from where you will copy them onto the `main/scripts` directory that you will have already generated.

```sh
# Run from the `01_PAML/01_MCMCtree/scripts`
# dir on your local PC. Please change
# directories until you are there. Then run
# the following commands.
cp generate_job_MCMCtree.sh ../../../main/scripts/
cp pipeline_MCMCtree_template.sh ../../../main/scripts/
```

The `main` directory will now have these two extra directories with the corresponding subdirectories:

```text
main
  |- MCMCtree
  |    |- [1-5]
  |         |- [GBM|ILN]
  |              |- [0-9] # As many as available datasets
  |                
  |- pipelines_MCMCtree
       |- [GBM|ILN]/
```

To prepare the script to run `MCMCtree`, you can run the bash scripts mentioned above using the commands below:

```sh
# Run from the `main` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
home_dir=$( pwd )
cd scripts
chmod 775 *sh
num_aln=2
num_chains=5
for i in `seq 1 $num_aln`
do 
printf "Generating job array for dir "$i" for both clocks, both sets of calibrations, and both partitioning schemes ... ...\n\n"
# The arguments are the following
# arg1   Alignment #1, #2, #3... The alignment being parsed at the moment.
# arg2   Clock model: e.g., "GBM" or "ILN".
# arg3   Number of partitions: 1, 2, 3... As many blocks as partitions in the alignment.
# arg4   Path to MCMCtree pipeline dir.
# arg5   Command to execute MCMCtree. E.g. "mcmctree", "mcmctree_4.10.7", etc.
# arg6   Number of MCMCs to run.
# arg7   Type of calibration: "seccal" (secondary calibrations) or "fosscal" (fossil calibrations)
./generate_job_MCMCtree.sh $i GBM 1 $home_dir/pipelines_MCMCtree mcmctree $num_chains seccal
./generate_job_MCMCtree.sh $i GBM 1 $home_dir/pipelines_MCMCtree mcmctree $num_chains fosscal
./generate_job_MCMCtree.sh $i GBM 1 $home_dir/pipelines_MCMCtree mcmctree $num_chains seccal
./generate_job_MCMCtree.sh $i ILN 1 $home_dir/pipelines_MCMCtree mcmctree $num_chains fosscal
./generate_job_MCMCtree.sh $i ILN 1 $home_dir/pipelines_MCMCtree mcmctree $num_chains seccal
done
```

Now, before running `MCMCtree` to sample from the posterior, we will run `MCMCtree` but sampling from the prior. In this way, we can verify whether the user-specified prior (the probability distributions we have used to calibrate the nodes in our phylogeny, and thus the ones we want the software to use) and the effective prior (the probability distribution the software will actually use when mathematical constraints are applied). Oftentimes, to deal with truncation issues, the effective prior might differ from the user-specified prior and not all possible ages framed by the distribution specified by the user will be included (see an extensive analysis about this effect in [dos Reis et al. 2015](https://pubmed.ncbi.nlm.nih.gov/26603774/)). To assess whether this occurs in our analysis, we will run first `MCMCtree` without the data (i.e., `MCMCtree` will be sampling from the prior distribution and will only use the information provided by the calibrations included in the user tree, not the sequence alignment).

First, we will generate a directory where `MCMCtree` will run and where all the results when sampling from the prior will be saved. Given that we are not using the data and only relying on additional information such as the calibrations, it does not matter whether the data are partitioned or not. In that way, we will keep only one directory that tests both fossil calibrations:

```sh
# Run from `main` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
num_chains=5
for i in `seq 1 $num_chains`
do
mkdir -p MCMCtree_prior/$i/{fosscal,seccal}
done
```

>**IMPORTANT NOTE**: When sampling from the prior, the likelihood is not being calculated or estimated. In other words, the most time-consuming part of the MCMC does not take place. In that way, you should be able to gather enough samples with fewer runs than those needed when sampling from the posterior. E.g., at least 2 runs and, perhaps, no more than 5 chains if you keep using the settings we have specified in the template control files in this tutorial.

Then, we will copy the directory `pipelines_MCMCtree` and will generate a copy called `pipelines_MCMCtree_prior`, which file structure we will modify:

```sh
# Run from `main` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
cp -R pipelines_MCMCtree pipelines_MCMCtree_prior
cd pipelines_MCMCtree_prior
```

We will modify the bash script that will be submitted as a job array so the `userdata` option in the control file is equal to `0`, which enables `MCMCtree` to sample from the prior instead of sampling from the posterior (i.e., the data alignment is ignored). We will also set the nucleotide model to `0` (i.e., JC69) given that no data are being used and hence this option is not enabled when sampling from the prior. We will specify the strict clock (i.e., `clock = 1`, the default option in the template control file) given that we are not using the data to estimate evolutionary rates -- in fact, that option is not even enable, so the option chosen would not really matter! Last, we will change the path to where the results will be stored:

```sh
# Run from `pipelines_MCMCtree_prior` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
# 
# Prepare directories to sample from the prior,
# only one needed as nucleotide subsitution models 
# are not used.
rm -r conc
mv part/fosscal .
mv part/seccal .
rm -r part
rm -r */ILN
mv fosscal/GBM/*sh fosscal/
mv seccal/GBM/*sh seccal
for i in */*sh
do
name=$( echo $i | sed 's/GBM\_//' | sed 's/\_part//' )
mv $i $name
done
rm -r */GBM

# Modify path to save results 
sed -i 's/MCMCtree/MCMCtree\_prior/' */*sh
sed -i 's/\/GBM//' */*sh
# Comment soft link
sed -i 's/ln \-s/\#ln \-s/' */*sh

# Modify bash script: options `usedata` and `model`
sed -i 's/part\/fosscal/fosscal/g' fosscal/*sh
sed -i 's/part\/seccal/seccal/g' seccal/*sh
sed -i 's/dir\$dir\"\_r\"/fosscal\_r/g' fosscal/*sh
sed -i 's/dir\$dir\"\_r\"/seccal\_r/g' seccal/*sh
sed -i 's/\_part\_\$dir\"\_r\"/\_r/' */*sh

# Modify bash script: options `usedata` and `model`
sed -i "s/..*usedata..*/sed \-i \'s\/usedata\.\.\*\/usedata \= 0\/\' \$home\_dir\/mcmctree\_fosscal_r\$TASK\_ID\"\.ctl\"\nsed \-i \'s\/model\.\.\*\/model \= 0\/\' \$home\_dir\/mcmctree\_fosscal_r\$TASK\_ID\"\.ctl\"/" fosscal/*sh
sed -i "s/..*usedata..*/sed \-i \'s\/usedata\.\.\*\/usedata \= 0\/\' \$home\_dir\/mcmctree\_seccal_r\$TASK\_ID\"\.ctl\"\nsed \-i \'s\/model\.\.\*\/model \= 0\/\' \$home\_dir\/mcmctree\_seccal_r\$TASK\_ID\"\.ctl\"/" seccal/*sh
# Modify clock model 
sed -i "s/^fi/fi\nsed \-i \'s\/clock\.\.\*\/clock \= 1\/\' \$home\_dir\/mcmctree\_fosscal_r\$TASK\_ID\"\.ctl\"/" fosscal/*sh
sed -i "s/^fi/fi\nsed \-i \'s\/clock\.\.\*\/clock \= 1\/\' \$home\_dir\/mcmctree\_seccal_r\$TASK\_ID\"\.ctl\"/" seccal/*sh
```

Now, you can check that the lines have been correctly modified:

```sh
# Run from `pipelines_MCMCtree_prior` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
grep 'usedata' */*sh
grep 'model'  */*sh
grep 'MCMCtree_prior' */*sh
```

## 2. Analyses with `MCMCtree` when sampling from the prior

Now, we will be able to run `MCMCtree` first when sampling from the prior (i.e., no data used!) using the code snippet below:

```sh
# The command below will start a process in which `MCMCtree`
# will run when sampling from the prior under the strict clock. 
# E.g., you will see a line in the format of `[1] 902` will
# be printed after the command `./pipeline_CLK_fosscal.sh &` or
## `./pipeline_CLK_seccal.sh`. The `&`  
# allows you to keep using the same terminal to run more commands.
# You can use the same terminal now to run the next commands 
# or open a new terminal to continue.

# Run from `pipelines_MCMCtree_prior/` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
chmod 775 */*sh
cd seccal
./pipeline_seccal.sh &
cd ../fosscal
./pipeline_fosscal.sh &
```

### 2.1. Setting the file structure

We will now create a `sum_analyses` directory to analyse the `MCMCtree` output.
First, we need to come back to the `01_analyses/01_MCMCtree` directory and run the following code snippet:

```sh
# Run from `01_PAML/01_MCMCtree` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands. If you are running this code with your
# own analyses, make sure that you have correctly
# defined `num_aln` and `num_chains` variables with
# the correct values!
# Note that we will generate some directories for
# when the analyses when sampling from the posterior
# are ready!
mkdir sum_analyses
cd sum_analyses
num_chains=5
for j in `seq 1 $num_chains`
do
mkdir -p posterior/{conc,part}/{fosscal,seccal}/{GBM,ILN}/$j
mkdir -p prior/{fosscal,seccal}/$j
done

# Now, you can start copying the
# data we only need!
# Run the commands below as they are written if you
# are still inside `01_PAML/01_MCMCtree/sum_analyses`.
# If not, run them once you are in this directory.
for i in `seq 1 $num_chains` 
do
printf "\n[[ Copying run "$i" into the corresponding directories ]]\n\n"
cp ../../../main/MCMCtree_prior/$i/fosscal/mcmc.txt prior/fosscal/$i
cp ../../../main/MCMCtree_prior/$i/fosscal/*ctl prior/fosscal/$i
cp ../../../main/MCMCtree_prior/$i/seccal/mcmc.txt prior/seccal/$i
cp ../../../main/MCMCtree_prior/$i/seccal/*ctl prior/seccal/$i
done
```

### 2.2. MCMC diagnostics

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to run the MCMC diagnostics!

We are going to run an R script I wrote, [`MCMC_diagnostics_prior.R`](scripts/MCMC_diagnostics_prior.R), and follow the detailed step-by-step instructions detailed in the script. In a nutshell, the protocol we will follow is the following:

1. Load the `mcmc.txt` files generated after each run.
2. Generate a convergence plot with the unfiltered chains.
3. Find whether there are major differences between the time estimates sampled across the chains for the same nodes in the 97.5% and the 2.5% quantiles. If so, flag and delete said chains.
4. If some chains have not passed the filters mentioned above, create an object with the chains that have passed the filter.
5. Generate a new convergence plot with those chains that passed filters.
6. Calculate the effective sample size for all model parameters and check whether chain convergence has been reached with the chains that have passed filters.

When you finish to run the R script, answer the following questions:

----

#### **EXERCISE 1**

* Write down the tail-ESS, the bulk-ESS, and the Rhat. What do these stats suggest with regards to chain convergence?
* What other tests could you use to assess chain autocorrelation and chain efficiency?

----

If the MCMC diagnostics did not find any of the chains problematic, then you can run the R script [`MCMC_diagnostics_prior.R`](scripts/MCMC_diagnostics_prior.R). Before that, however, we need to use an in-house bash I wrote, [`Combine_MCMC.sh`](scripts/Combine_MCMC.sh), to concatenate in a unique file all the `mcmc.txt` files that correspond to the analysis where chains passed the filters:

```sh
# Run from `01_MCMCtree/scripts`
cp Combine_MCMC.sh ../sum_analyses/prior
# One argument taken: number of chains
cd ../sum_analyses/prior
## Variables needed
## arg1   Mame of directory where analyses have taken place (e.g., `fosscal`, `seccal`)
## arg2   Name of the output directory: mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3   Filtered chains. E.g., `seq 1 36`, "1 2 5", etc. 
## arg4   Clock model: ILN, GBM, CLK
./Combine_MCMC.sh fosscal mcmc_files_fosscal "`seq 1 5`" CLK
./Combine_MCMC.sh seccal mcmc_files_seccal "`seq 1 5`" CLK
```

The script above will generate two directories called `mcmc_files_fosscal` and `mcmc_files_seccal` inside the `prior` directory, where the `mcmc.txt` with the concatenated samples will be saved. A template script to generate the `FigTree.tre` file with this `mcmc.txt` has been saved inside the [`dummy_ctl_files`](dummy_ctl_files) directory. We now will create a dummy alignment with only 2 nucleotides to generate the `FigTree` files using the concatenated `mcmc.txt` files. In order to do that, we can run the [`Generate_dummy_aln.R`](scripts/Generate_dummy_aln.R). Once you run it, a new directory called `dummy_aln` will be created, which will contain the dummy alignment. We have also generated a dummy control file with option `print = -1`, which will not run an MCMC but, instead, will use the input files (file with the dummy alignment, calibrated tree file, and concatenated `mcmc.txt` file) to generate a `FigTree.tre` file with the mean estimated divergence times and the corresponding mean CIs using all the samples collected during all the MCMCs.

```sh
# Run from `sum_analyses/prior`
base_dir=$( pwd )
cd ../../dummy_ctl_files
ctl_dir=$( pwd )
cd ../../../00_data_formatting/01_inp_data
tt_dir=$( pwd )
name_tt=`ls *_fosscal.tree`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_fosscal
printf "[[ Generating tree file for concatenated \"mcmc.txt\"  ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
printf "\n"
mv FigTree.tre FigTree_fosscal.tree
cd $tt_dir
name_tt=`ls *_seccal.tree`
cd $base_dir/mcmc_files_seccal
printf "[[ Generating tree file for concatenated \"mcmc.txt\"  ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
printf "\n"
mv FigTree.tre FigTree_seccal.tree
cd $base_dir
```

The next step is to plot the user-specified prior VS the effective prior. I have written the R script [`Check_priors_effVSuser.R`](scripts/Check_priors_effVSuser.R) to generate these plots. Once this script has finished, you will see that a new directory `plots/effVSuser` will have been created. Inside this directory, you will find one directory for each individual dataset with individual plots for each node. In addition, all these plots have been merged into a unique document as well (note: some plots may be too small to see for each node, hence why you will also generate individual plots).

Now, once the MCMC diagnostics have finished, you can extract the final results when sampling from the posterior so that you can go through them in a much easier manner:

```sh
# Run from `01_MCMCtree`
mkdir sum_files_prior
cp -R sum_analyses/prior/mcmc_files_*/*tree sum_files_prior/
cp -R sum_analyses/prior/*cal/*cal/*all_mean*tsv sum_files_prior/
cp -R plots/ESS_and_chains_convergence/*prior*pdf sum_files_prior/
cp -R plots/effVSuser sum_files_prior/
```

## 3. Analyses with `MCMCtree` when sampling from the posterior

Now that we have verified that there are no issues between the user-specified prior and the effective prior, we can run `MCMCtree` when sampling from the posterior. We will do these analyses under the GBM and ILN relaxed-clock models using the code snippet below. Please run only those commands that correspond to the datasets your team has been assigned to run:

```sh
# Go to directory `main/pipelines_MCMCtree/` dir on your local
# PC and run the following commands. Please change directories until
# you are there.
chmod 775 */*/*/*sh
cd conc/fosscal/GBM/
./pipeline_GBM_fosscal_conc.sh &
cd ../ILN
./pipeline_ILN_fosscal_conc.sh &
cd ../../seccal/GBM
./pipeline_GBM_seccal_conc.sh &
cd ../ILN
./pipeline_ILN_seccal_conc.sh &
cd ../../../part/fosscal/GBM/
./pipeline_GBM_fosscal_part.sh &
cd ../ILN
./pipeline_ILN_fosscal_part.sh &
cd ../../seccal/GBM
./pipeline_GBM_seccal_part.sh &
cd ../ILN
./pipeline_ILN_seccal_part.sh &
```

### 3.1. Setting the file structure

We will now go to the previously created `sum_analyses` directory to analyse the `MCMCtree` output.
First, we need to come back to the `01_analyses/01_MCMCtree` directory and run the following code snippet:

```sh
# It is assumed that you have already
# generated the `posterior` directory 
# inside the `01_MCMCtree/sum_analyses`
# directory following the commands
# detailed in `section 2.1`.
# Run the commands below 
# inside `01_PAML/01_MCMCtree/sum_analyses`.
## NOTE: If you have not run all the analyses with your group,
## please modify the content of the `for` loop so only
## those directories in `main/MCMCtree` that have `MCMCtree`
## output files are parsed!
num_chains=5
for j in `seq 1 $num_chains` 
do
printf "\n[[ Copying run "$j" into the corresponding directories ]]\n\n"
cp ../../../main/MCMCtree/$j/conc/fosscal/GBM/mcmc.txt posterior/conc/fosscal/GBM/$j
cp ../../../main/MCMCtree/$j/conc/fosscal/GBM/*ctl posterior/conc/fosscal/GBM/$j
cp ../../../main/MCMCtree/$j/conc/seccal/GBM/mcmc.txt posterior/conc/seccal/GBM/$j
cp ../../../main/MCMCtree/$j/conc/seccal/GBM/*ctl posterior/conc/seccal/GBM/$j

cp ../../../main/MCMCtree/$j/conc/fosscal/ILN/mcmc.txt posterior/conc/fosscal/ILN/$j
cp ../../../main/MCMCtree/$j/conc/fosscal/ILN/*ctl posterior/conc/fosscal/ILN/$j
cp ../../../main/MCMCtree/$j/conc/seccal/ILN/mcmc.txt posterior/conc/seccal/ILN/$j
cp ../../../main/MCMCtree/$j/conc/seccal/ILN/*ctl posterior/conc/seccal/ILN/$j

cp ../../../main/MCMCtree/$j/part/fosscal/GBM/mcmc.txt posterior/part/fosscal/GBM/$j
cp ../../../main/MCMCtree/$j/part/fosscal/GBM/*ctl posterior/part/fosscal/GBM/$j
cp ../../../main/MCMCtree/$j/part/seccal/GBM/mcmc.txt posterior/part/seccal/GBM/$j
cp ../../../main/MCMCtree/$j/part/seccal/GBM/*ctl posterior/part/seccal/GBM/$j

cp ../../../main/MCMCtree/$j/part/fosscal/ILN/mcmc.txt posterior/part/fosscal/ILN/$j
cp ../../../main/MCMCtree/$j/part/fosscal/ILN/*ctl posterior/part/fosscal/ILN/$j
cp ../../../main/MCMCtree/$j/part/seccal/ILN/mcmc.txt posterior/part/seccal/ILN/$j
cp ../../../main/MCMCtree/$j/part/seccal/ILN/*ctl posterior/part/seccal/ILN/$j
done
```

### 3.2. MCMC diagnostics

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to check for chain convergence!

Depending on the group you are part of, please run the corresponding `MCMC_diagnostic_posterior*.R` file inside the [`scripts` directory](scripts) and follow the detailed step-by-step instructions written in the script -- they are essentially the same ones used when analysing the chains when sampling from the prior. If all the checks above have been passed, then you are ready to generate the summarised `mcmc.txt` file with the samples collected by the chains that passed the filters:

```sh
# Run from `01_MCMCtree/scripts`
cp Combine_MCMC.sh ../sum_analyses/posterior
# One argument taken: number of chains
cd ../sum_analyses/posterior
## Variables needed
## arg1   Name of directory where analyses have taken place (e.g., conc/fosscal/GBM).
## arg2   Name of the output directory: mcmc_files_conc_fosscal_GBM, mcmc_files_conc_fosscal_ILN, etc.
## arg3   Depends on whether some chains were filtered out or not. E.g., `seq 1 36`, "1 2 5", etc.
## arg4   Clock model: ILN, GBM, CLK.
##
## NOTE: Please only run those commands that correspond to the analyses
## you have run!
./Combine_MCMC.sh conc/fosscal/GBM mcmc_files_conc_fosscal_GBM "`seq 1 5`" GBM
./Combine_MCMC.sh conc/fosscal/ILN mcmc_files_conc_fosscal_ILN "`seq 1 5`" ILN
./Combine_MCMC.sh part/fosscal/GBM mcmc_files_part_fosscal_GBM "`seq 1 5`" GBM
./Combine_MCMC.sh part/fosscal/ILN mcmc_files_part_fosscal_ILN "`seq 1 5`" ILN
./Combine_MCMC.sh conc/seccal/GBM mcmc_files_conc_seccal_GBM "`seq 1 5`" GBM
./Combine_MCMC.sh conc/seccal/ILN mcmc_files_conc_seccal_ILN "`seq 1 5`" ILN
./Combine_MCMC.sh part/seccal/GBM mcmc_files_part_seccal_GBM "`seq 1 5`" GBM
./Combine_MCMC.sh part/seccal/ILN mcmc_files_part_seccal_ILN "`seq 1 5`" ILN
```

Once the scripts above have finished, eight new directories called `mcmc_files_*` will be created inside the `posterior` directory. To infer the final timetrees with the mean time estimates using the samples collected by the filtered chains, we need a control file, the calibrated Newick tree, and the dummy alignment we previously generated in section 2 inside this directory:

```sh
# Run from `sum_analyses/posterior` directory.
# Please change directories until
# you are there. 
# Depending on the dataset you are to analyse,
# please change the following variables called 
# `dat_set` and `calib` accoridngly:
dat_set=$( echo conc_fosscal_GBM ) # conc_fosscal_GBM | conc_fosscal_ILN | conc_seccal_GBM | conc_seccal_ILN 
                                   # part_fosscal_GBM | part_fosscal_ILN | part_seccal_GBM | part_seccal_ILN 
calib=$( echo fosscal ) # fosscal | seccal
# Now, run the following commands
# to summarise the samples collected
base_dir=$( pwd )
cd ../../dummy_ctl_files
ctl_dir=$( pwd )
cd ../../../00_data_formatting/01_inp_data
tt_dir=$( pwd )
name_tt=`ls *_$calib".tree"`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_$dat_set
printf "[[ Generating tree file for concatenated \"mcmc.txt\" ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_$dat_set.tree
printf "\n"
cd $base_dir
```

Now, once the MCMC diagnostics have finished, you can extract the results generated when sampling from the posterior so that you can go through them in a much easier way:

```sh
# Run from `01_MCMCtree`
mkdir sum_files_post
cp -R sum_analyses/posterior/mcmc_files_*/*tree sum_files_post/
cp -R sum_analyses/posterior/*/*/*/*/*all_mean*tsv sum_files_post/
cp -R plots/ESS_and_chains_convergence/*post*pdf sum_files_post/
```

----

#### **EXERCISE 2**

* Compare the time estimates you have obtained under the different clock models, partitioning schemes, and sets of calibrations. Do you observe any clear pattern?
* Plot the divergence times estimated when sampling from the prior vs those from the posterior (e.g., write your own R code, use [`Tracer`](https://github.com/beast-dev/tracer/releases/tag/v1.7.2), etc.). What can you tell?

> TIP: You may one to read [Nascimento et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5624502/) to help you with such discussions!

----

If you have any more questions, let's just discuss them at the end of the session :)
