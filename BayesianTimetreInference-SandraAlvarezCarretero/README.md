# Bayesian timetree inference

**DISCLAIMER**: This tutorial is based on a phylogenetics tool that I am working on at the moment, which I am still developing and is yet to be published. While some of the scripts/tools that you will find here have been validated and used in published research ([Álvarez-Carretero et al., 2022](https://doi.org/10.1038/s41586-021-04341-1)), I am actively implementing new features as part of the current workflow of this pipeline as well as developing new scripts/tools. In other words, the code is not stable and I am still validating the new features. If you want to use the tools that I have developed as part of this tutorial/pipeline, please first contact me at <a href="mailto:sandra.ac93@gmail.com"><b>sandra.ac93@gmail.com</b></a>. Thank you :)

**NOTE**: In case you stop having access to this repository, and thus to this tutorial, please note that I have created a repository in my personal GitHub account to host these resources: [**MolClockDat-tutorial**](https://github.com/sabifo4/MolClockDat-tutorial).

## Overview

During this practical session, we will use various in-house scripts and tools together with `BASEML` and `MCMCtree`, two programs that are part of the `PAML` software ([Yang 2007](https://pubmed.ncbi.nlm.nih.gov/17483113/)), to run a **Bayesian analysis for timetree inference** using a **node-dating approach** and an **approximation to the likelihood calculation** implemented in `MCMCtree` ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) to speed up such analysis.

Before you get started, let's look at the data that you have managed to collect so far:

* [**Molecular alignment**](raw_data/aln_seq/baliphy/P1-max.fasta): alignments inferred with `BAli-Phy` ([Redelings, 2021](https://academic.oup.com/bioinformatics/article/37/18/3032/6156619)) using the CYTB sequences of 8 mammals downloaded from the NCBI web server. For more information on how this molecular alignment has been generated, you can always [read the step-by-step tutorial in the `raw_data` directory](raw_data/README.md). To make sure you have time to go through this tutorial, however, please do not go through this dataset assembly pipeline until the end of this session.
* [**Set of calibrations #1**](00_data_formatting/calibs/Calib_converter_1.txt): calibrations based on the fossil record. More information about how they have been established later in the tutorial.
* [**Set of calibrations #2**](00_data_formatting/calibs/Calib_converter_2.txt): secondary calibrations for the same nodes in set #1 based on the divergence times inferred by [Álvarez-Carretero et al., 2022](https://doi.org/10.1038/s41586-021-04341-1) using the Bayesian sequential-subtree (BSS) approach. More information about how they have been established later in the tutorial.
* [**Phylogeny**](raw_data/aln_seq/baliphy/cytb_rooted_bl.tree): tree topology with branch lengths inferred with `BAli-Phy` ([Redelings, 2021](https://academic.oup.com/bioinformatics/article/37/18/3032/6156619)). For more information about how this phylogeny has been inferred, you can always [read the step-by-step tutorial in the `raw_data` directory](raw_data/README.md). To make sure you have time to go through this tutorial, however, please do not go through this dataset assembly pipeline until the end of this session.

Remember that, before proceeding with timetree inference, you need to get familiar with your dataset! For instance, you may ask yourselves questions such as "how were the data collected?", "how were the alignments generated?", or "how are the files going to be organised?". In this practical session, we are not going to address such questions as we will only have time to focus on the subsequent steps. Nevertheless, at the end of the session, you can always go through the step-by-step data assembly workflow available in the [`raw_data` directory](raw_data) to understand how the data were collected and processed. In addition, to learn more about the different steps that the phylogenomic workflow consists of, you may also want to read [Álvarez-Carretero & dos Reis, 2022](https://link.springer.com/chapter/10.1007/978-3-030-60181-2_13).

### Goals

At the end of this practical session, you should...

* ... be mindful about how important it is to be familiar with your dataset before proceeding with timetree inference.
* ... understand how to parse and format input data.
* ... understand how to run `PAML` software for timetree inference analysis. E.g., specifying substitution models, selecting the most adequate priors according to your dataset, specifying MCMC settings, etc.
* ... understand how to run MCMC diagnostics to confidently filter the chains that were run under each type of analysis and assess chain convergence.
* ... be able to critically discuss the results you have obtained according to your prior hypotheses and the settings under which phylogenetic software have been executed.

### Workflow

The summary of the workflow that you will follow during this practical session is the following:

* Parsing and formatting the input data required to run `PAML` software: files with the sequence alignment and a fixed tree topology.
* Inferring the mean evolutionary rate to specify a sensible rate prior.
* Running `PAML` software for timetree inference:
  * Using various in-house pipelines to set up the working environment, the file structure, and the control files required to run `PAML` software - feel free to go back to the lecture content if you have questions about the options that are to be specified in `PAML` control files (e.g., nucleotide substitution model, tree prior, rate prior, MCMC settings, etc.).
  * Running `BASEML` to calculate the gradient and the Hessian, which will be required by `MCMCtree` to enable the approximate likelihood calculation.
  * Running `MCMCtree` with the approximate likelihood calculation enabled for timetree inference. We will run the following analyses to assess the impact that (i) different sets of calibrations, (ii) relaxed-clock models (i.e., autocorrelated-rates model or geometric Brownian diffusion model [[Thorne et al. 1998](http://www.ncbi.nlm.nih.gov/pubmed/9866200), [Yang and Rannala 2006](http://www.ncbi.nlm.nih.gov/pubmed/16177230)] and independent log-normal rate model [[Rannala and Yang 2007](http://www.ncbi.nlm.nih.gov/pubmed/17558967), [Lemey et al. 2010](http://www.ncbi.nlm.nih.gov/pubmed/20203288)]), and (iii) partitioning schemes (i.e., all codon positions or only the first and the second codon positions) can have on species divergence times estimation:
    * **GBM_CAL1_ALLCP**: autocorrelated-rates model + set of calibrations #1 + 1st/2nd/3rd codon positions.
    * **GBM_CAL1_12CP**: autocorrelated-rates model + set of calibrations #1 + 1st/2nd codon positions.
    * **GBM_CAL2_ALLCP**: autocorrelated-rates model + set of calibrations #2 + concatenated codon positions.
    * **GBM_CAL2_12CP**: autocorrelated-rates model + set of calibrations #2 + 1st/2nd codon positions.
    * **ILN_CAL1_ALLCP**: independent-rates model + set of calibrations #2 + 1st/2nd/3rd codon positions.
    * **ILN_CAL1_12CP**: independent-rates model + set of calibrations #2 + 1st/2nd codon positions.
    * **ILN_CAL2_ALLCP**: independent-rates model + set of calibrations #2 + 1st/2nd/3rd codon positions.
    * **ILN_CAL2_12CP**: independent-rates model + set of calibrations #2 + 1st/2nd/3rd codon positions.
* Running MCMC diagnostics for all the chains under each analysis.
* General discussion.

> **NOTE 1**: we will go through the first two steps together so that everyone starts the `PAML` analyses once all the details with regards to data formatting and the rate prior have been understood.
>
> **NOTE 2**: to speed things up, not everyone will individually run the 8 `MCMCtree` analyses on their own when sampling from the posterior. Depending on how many people are attending this session, we will work in either 8 or 4 groups. If the former, each group of participants will run one analysis when sampling from the posterior while, if the latter, each group will run two.
>
> **NOTE 3**: we will stop 15 minutes before the end of the practical session to discuss the results that you have obtained by that time. If you have not finished the analyses, please do not worry! I will share with you the results under the 8 different scenarios, and we can all have a general discussion about how different calibrations, partitioning schemes, and relaxed-clock models have affected the estimated divergence times in our bear phylogeny.

## Software

Before you start this practical session, please make sure you have the following software installed on your PCs:

* **`PAML`**: you will be using the latest `PAML` release, v4.10.7, available from the [`PAML` GitHub repository](https://github.com/abacus-gene/paml). If you do not want to install the software from the source code, then follow (A). If you want to install `PAML` from the source code, then follow (B):

  * Installation (A): if you have problems installing `PAML` from the source code or you do not have the tools required to compile the source code, then you can [download the pre-compiled binaries available from the latest release by following this link](https://github.com/abacus-gene/paml/releases/tag/4.10.7). Please choose the pre-compiled binaries you need according to your OS (i.e., Windows, Mac OSX, Linux), download the corresponding compressed file and save it in your preferred directory. Then, after decompressing the file, please give executable permissions, export the path to this binary file so you can execute it from a terminal, and you should be ready to go!
  * Installation (B): to install `PAML v4.10.7` from the source code, please follow the instructions given in the code snippet below:

    ```sh
    # Clone to the `PAML` GitHub repository to get the latest `PAML` version
    # You can go to "https://github.com/abacus-gene/paml" and manually clone
    # the repository or continue below from the command line
    git clone https://github.com/abacus-gene/paml
    # Change name of cloned directory to keep track of version
    mv paml paml4.10.7
    # Move to `src` directory and compile programs
    cd paml4.10.7/src
    make -f Makefile
    rm *o
    # Move the new executable files to the `bin` directory and give executable
    # permissions
    mkdir ../bin
    mv baseml basemlg chi2 codeml evolver infinitesites mcmctree pamp yn00 ../bin
    chmod 775 ../bin/*
    ```
  
    Now, you just need to export the path to the `bin` directory where you have saved the executable file. If you want to automatically export this path to your `./bashrc` or your `~/.bash_profile`, you can run the following commands **AFTER ADAPTING** the absolute paths written in the code snippet below to those in your filesystem:

    ```sh
    # Run from any location. Change `~/.bashrc` if you are 
    # using another file
    printf "\n# Export path to PAML\n" >> ~/.bashrc
    # Replace "/c/Users/Bioinfor_tools/" with the path
    # that leads to the location where you have saved the
    # `paml4.10.7` directory. Modify any other part of the
    # absolute path if you have made other changes to the 
    # name of the directory where you have downloaded `PAML`
    printf "export PATH=/c/usr/bioinfo_tools/paml4.10.7/bin:\$PATH\n" >>  ~/.bashrc
    # Now, source the `~/.bashrc` file (or the file you are 
    # using) to update the changes
    source ~/.bashrc
    ```

    Alternatively, you can edit this file using your preferred text editor (e.g., `vim`, `nano`, etc.).

* **`R`** and **`RStudio`**: please download [R](https://cran.r-project.org/) and [RStudio](https://posit.co/download/rstudio-desktop/) as we will be using them throughout the practical. The packages we will be using should work with R versions that are either newer than or equal to v4.1.2. If you are a Windows user, please make sure that you have the correct version of `RTools` installed, which will allow you to install packages from the source code if required. For instance, if you have R v4.1.2, then installing `RTools4.0` shall be fine. If you have another R version installed on your PC, please check whether you need to install `RTools 4.2` or `RTools 4.3`. For more information on which version you should download, [please go to the CRAN website by following this link and download the version you need](https://cran.r-project.org/bin/windows/Rtools/).

    Before you proceed, however, please make sure that you install the following packages too:

    ```R
    # Run from the R console in RStudio
    # Check that you have at least R v4.1.2
    version$version.string
    # Now, install the packages we will be using
    # Note that it may take a while if you have not 
    # installed all these software before
    install.packages( c('rstudioapi', 'ape', 'phytools', 'sn', 'stringr', 'rstan', 'colorBlindness'), dep = TRUE )
    ## NOTE: If you are a Windows user and see the message "Do you want to install from sources the 
    ## packages which need compilarion?", please make sure that you have installed the `RTools`
    ## aforementioned.
    ```

* **`FigTree`**: you can use this graphical interface to display tree topologies with/without branch lengths and with/without additional labels. You can then decide what you want to be displayed by selecting the buttons and options that you require for that to happen. You can [download the latest pre-compiled binaries, `FigTree v1.4.4`, from the `FigTree` GitHub repository](https://github.com/rambaut/figtree/releases).

* **`Tracer`**: you can use this graphical interface to visually assess the MCMCs you have run during your analyses (e.g., chain efficiency, chain convergence, autocorrelation, etc.). You can [download the latest pre-compiled binaries, `Tracer v1.7.2`, from the `Tracer` GitHub repository](https://github.com/beast-dev/tracer/releases/tag/v1.7.2).

## Data analysis

If you have gone through the previous sections and have a clear understanding of the dataset you will be using, the workflow of the analyses you will be running, and have installed the required software to do so... Then you are ready to go!

You can start formatting the input data needed to infer the divergence times of our bear phylogeny [by following this link](00_data_formatting/README.md).

Happy timetree inference! :)

----

ⓒ Dr Sandra Álvarez-Carretero | [`@sabifo4`](https://github.com/sabifo4/)
