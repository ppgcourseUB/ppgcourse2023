# Hessian and gradient calculation

Before running `MCMCtree`, we need to calculate the gradient and the Hessian so we can use the approximate likelihood during timetree inference to save computational time ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)).

## 1. Pick rate prior

We will use a vague gamma distribution for the dataset considering the tree height (molecular distance in substitutions per site) and the divergence time at the root of the phylogeny (in time unit). As the [tree file we inferred with BAli-Phy](../../00_data_formatting/00_raw_data/cytb_rooted_bl.tree) has information about the branch lengths, we can load this file in `R` to estimate the tree height. We also have the root calibration in the calibrated tree file that we just generated, which is a rough idea of the age of the root of the phylogeny based on the fossil record. Following the first set of calibrations, the mean root age is 0.519x100 Ma (i.e., average between the minimum, 37.71 Ma, and the maximum, 66.09 Ma, soft bounds defined in the uniform distribution), while the one using the second set of calibrations is 0.2338x100 Ma (i.e., minimum age 20.16 Ma and maximum age is 26.6 Ma). We will estimate the mean evolutionary rate under both scenarios.

By setting a vague shape ($\alpha=2$) for the gamma distribution, we can account for the uncertainty on the mean rate. If we had more knowledge on the mean rate, however, we should use a narrower prior with a larger $\alpha$ that better represents our prior information.

Now, we have all the information we need to calculate the $\beta$ parameter for the Gamma distribution that will be the prior on the estimated mean evolutionary rate. I have written the [R script `calculate_rateprior.R`](scripts/calculate_rateprior.R) to carry out all the tasks mentioned above. You can open this file in RStudio to find out the value of $\beta$ for each analysis depending on the root calibration being used, and plot the final prior on the rates. A summary of what you will find in the script is described below:

```text
First, we know that the molecular distance (tree height, distance from the root to present time) is equal to the mean evolutionary rate (in substitutions per site per year) times the age of the divergence time at the root (in time unit, which we can define later). If we have estimated our phylogeny, and therefore have estimated the branch lengths, we will be able to estimate the tree height. The units of the tree height will be the following:

tree_height = rate * root_age --> units_tree_height = subst/site/y * y = subst/site

One way of estimating the tree height is by using the R function `phytools::nodeHeights`. The maximum height calculated by this function corresponds to the length from the root to the heighest tip. 

After estimating the tree height of our phylogeny (in subst/site) and considering the age of the root based on fossils (time unit = 100 Ma), we can get a rough estimate of the mean rate depending. We will calculate the mean rate using two different time units:

Time unit = 100 Ma (mean root age in Ma) --> mean_rate = tree_height / root_age = (subst/site) / (Ma) = subst/site per time unit (time unit = 100 Ma = 10^8 years)  --> (subst/site)/10^8 years

We also know that the mean of the gamma distribution is our parameter of interest: the mean evolutionary rate. Therefore:

mean_G = mean_rate = alpha / beta 
Time unit = 100 Ma: mean_rate = alpha / beta --> beta = alpha / mean_rate = 2 / mean_rate

In that way, the calibrated tree should be in this time unit too (i.e., do not forget to scale the calibrations accordingly if needed!). 
```

If you run the [R script `calculate_rateprior.R`](scripts/calculate_rateprior.R), you will see how all the steps described above take place and a new PDF file with the prior distribution to be used in each analysis will be generated in a new directory called `out_RData`.

Two [template control files](control_files) with the $\alpha$ and $\beta$ parameters (as defined using the R script above) for the gamma distribution as a prior on the rates have also been generated. Note that several options will be subsequently modified to fit the analysis with this dataset (i.e., you will see some options that have flags in capital letters, which will be replaced with the correct value for said option). Given the tree is not too shallow, the clock may not be seriously violated, and thus we have fixed a mean for the `sigma2` parameter (i.e., variation in the clock) as 0.05 using a gamma prior with $\alpha=2$ and $\beta=40$: `sigma2_gamma 2 40​`.

## 2. Set up the file structure

Before running `MCMCtree` using the approximate likelihood calculation to speed up timetree inference, we first need to calculate the gradient and the Hessian. We will use `BASEML` for this purpose!

The file structure we will use is the following:

```text
main/
  |
  |- alignments/
  |    |- X/ # Directory for alignment X -- have as many directories as alignments
  |       
  |- control_files/ # Pre-defined control files with flags to be later replaced with specific settings
  |
  |- Hessian/
  |    |- X # Directory for alignment X -- have as many directories as alignments
  |          
  |- pipelines_Hessian # Directory where the pipeline to run `BASEML` will be executed
  |
  |- scripts # Scripts used to prepare control files to run `BASEML
  |
  |- trees
      |- calibrated   # Directory with the calibrated tree for `MCMCtree`
      |- uncalibrated # Directory with the uncalibrated tree for `BASEML`
```

To create the `main` file structure in your PC, we run the following commands:

```sh
# Run the following commands from the
# main working directory (i.e., where 
# directories `00_data_formatting`,
# `01_PAML`, etc. can be found)
mkdir main 
cd main
# We have two alignments, and hence we
# type 2
num_aln=2
for i in `seq 1 $num_aln`
do
mkdir -p alignments/$i
mkdir -p Hessian/$i/prepare_baseml
done
mkdir pipelines_Hessian
mkdir control_files
mkdir -p trees/{uncalibrated,calibrated}
mkdir scripts
```

Once the file structure is created, we can now populate it with the input files we have generated some minutes ago: alignment files, tree files, and control files. E.g.:

```sh
# Run from `00_data_formatting/01_inp_data`
# Please change directories until 
# you are there. Then, run the 
# following commands.
#
# Copy input alignment and tree files 
cp aln.phy ../../main/alignments/1/
cp aln_12CP.phy ../../main/alignments/2/
cp *_uncalib.tree ../../main/trees/uncalibrated/
cp *_fosscal.tree ../../main/trees/calibrated/
cp *_seccal.tree ../../main/trees/calibrated/
# Change to where the control file is and transfer the 
# template control files
cd ../../01_PAML/00_Hessian/control_files/
cp *ctl ../../../main/control_files/
```

While the alignment files, the tree files, and the control files (with the already correct prior on the rates, MCMC settings, and flags to later replace with the corresponding parameter values) have already been generated and only need to be copied onto their corresponding directories as shown above, we need to generate other input files to estimate the Hessian and the gradient: the input control files for `BASEML`.

To do this in a reproducible manner, you can use the [script `generate_prepbaseml.sh`](scripts/generate_prepbaseml.sh), which you can find in the [`01_PAML/00_Hessian/scripts` directory](scripts). You should copy this bash script in the `main/scripts` directory previously created:

```sh
# Run from the `01_PAML/00_Hessian/scripts`
# dir on your local PC. Please change
# directories until you are there. Then,
# run the following commands.
cp generate_prepbaseml.sh ../../../main/scripts
```

The [`generate_prepbaseml.sh` script](scripts/generate_prepbaseml.sh) needs two arguments: the directory name (i.e., `1`, `2`, etc.) and the name of the template control file. As we are using two alignment files, we will use `2` as the argument. Also, we are not going to be using any of the other settings specified in the control file now (e.g., we are not going to be using any of the settings that enable priors or MCMC settings, they are just pre-defined for later when using `MCMCtree` and save us some time). In that way, we can just randomly pick one of the template control files now so it gets basic information to run `BASEML` (i.e., input file names, nucleotide substitution model chosen):

```sh
# Run from `main/scripts` on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
chmod 775 *sh
./generate_prepbaseml.sh 1 tmp_fosscal.ctl
./generate_prepbaseml.sh 2 tmp_fosscal.ctl
```

To make sure that all the paths have been properly extracted, you can run the following code snippet:

```sh
# Run from `main/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */prepare_baseml/*ctl
grep 'treefile' */prepare_baseml/*ctl
```

## 2. `BASEML`

### Preparing input files

Now that we have the input files (alignment and tree files) and the instructions to run `BASEML` (control file), we will be manually running `MCMCtree` inside each `prepare_baseml` directory (see file structure above) in a special mode that launches `BASEML` for the purpose we want: calculating the gradient and the Hessian.

----

Before you proceed with this analysis, it is important to note a few things that you need to do while running `MCMCtree` in this type of analysis when using large phylogenomic data. First, you will see that `MCMCtree` starts parsing the first locus. Then, you will see something like the following printed on your screen (some sections may change depending on the PAML version you have installed on your PC!):

```text
*** Locus 1 ***
running baseml tmp0001.ctl
BASEML in paml version 4.10.6, November 2022
ns = 4          ls = 48
Reading sequences, sequential format..
Reading seq # 4: sp4
Sequences read..

48 site patterns read, 609 sites
Counting frequencies..
```

As soon as you see the last line, you will see that the `tmp000X*` files will have been created, and hence you can stop this run by typing `ctrl+C` on the terminal you used to run such command. In that way, the `for` loop will continue and proceed with the next alignment, and you should do exactly the same. Given that this dataset is small and will finish quick, you will see that you may not even have time to stop the run as it will end before you even have the chance of doing that! Nevertheless, you need to do this with larger datasets. Therefore, please do not stop the run until the `tmp000X*` files are created in case you are using large datasets.

**NOTE**: Remember that we do not want to run `MCMCtree` now, we just want to obtain the `tmp000X*` files to run `BASEML` to estimate the gradient and the Hessian.

----

Now, let's run the following code snippet, even though you will not even have time to kill the jobs as aforementioned!

```sh
# Run the next commands from `main/Hessian`
home_dir=$( pwd )
for i in `seq 1 2`
do
cd $home_dir/$i/prepare_baseml
mcmctree *ctl
done
cd $home_dir
```

You can check that all the needed files have been created by running the following command:

```sh
# Run from the `main/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */*/tmp0001.ctl | wc -l # You should get as many datasets as you have, in this case 2
```

In addition, you need to make sure that option `method = 1` is enabled, which will speed up the computation of the Hessian and the gradient:

```sh
# Run from the `main/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
sed -i 's/method\ \=\ 0/method\ \=\ 1/' */*/tmp0001.ctl
grep 'method = 1' */*/tmp0001.ctl | wc -l # You should get as many as datasets you have, in this case 2
```

### Running `BASEML`

Now that we have the control file ready to run `BASEML` as well as the required input files, we can run `BASEML`!

I have created a template bash script with flags
(i.e., see script `pipeline_Hessian_template.sh` in the [`scripts` directory](scripts)) that will be replaced with the appropriate values by another bash script (`generate_job_BASEML.sh`, also saved in the [`scripts` directory](scripts)). Please note that the second bash script will edit the template bash script according to the data alignment/s that will be analysed. Now, we just need to copy them onto the `main/scripts` directory:

```sh
# Run from `main` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
home_dir=$( pwd )
cp ../01_PAML/00_Hessian/scripts/generate_job_BASEML.sh $home_dir/scripts
cp ../01_PAML/00_Hessian/scripts/pipeline_Hessian_template.sh $home_dir/scripts
cd $home_dir/scripts
chmod 775 *sh
num_aln=2
for i in `seq 1 $num_aln`
do
# The first argument corresponds to the number of 
# gene alignments and the second to the path to 
# where the pipeline will be executed: `pipelines_Hessian`
./generate_job_BASEML.sh $num_aln $home_dir/pipelines_Hessian
done
```

Now, we just need to go to the `pipelines_Hessian` directory and run the script that you will now see there by using the commands below:

```sh
# Run from `main/pipelines_Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
#
# If you list the content of this directory,
# you will see the pipeline you will need 
# to execute in a bash script called
# `pipelines_BASEML.sh`
ls *
# Now, execute this bash script
chmod 775 *sh
./pipeline_BASEML.sh &
```

Once `BASEML` finishes, we are ready to generate the `in.BV` file that we will later use when running `MCMCtree` to approximate the likelihood calculation:

```sh
# Run from dir `main/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
num_aln=2
for i in `seq 1 $num_aln`
do 
printf "\nGenerating in.BV files for dir "$i"  ... ...\n\n"
mv $i/rst2 $i/in.BV
done
```

You can now [start the next tutorial to run `MCMCtree`](../01_MCMCtree/README.md)!
