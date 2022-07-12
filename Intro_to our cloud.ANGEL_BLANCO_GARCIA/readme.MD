# Intro to our cloud

Instructor: Àngel Blanco García


In this session we will get you started on how to use the cloud computing infrastructure available for this course.

## Computer clusters

In a very simplistic way, a cluster of computers is a group of **connected** computers execute tasks in a integrated way, as a single general system. Beyond this generic definition, there are different approaches to define the structure and behavior of a cluster, depending on its purpose. 

One of the most common purposes for a cluster is to enable a number of users to execute tasks that would be too resource intensive or time consuming to run on their personal computers. Users' tasks are distributed over the computers composing the cluster, based on priorities, resource availability, etc. Clusters of this kind are often called **computing clusters**. This way, many user tasks can be managed automatically by the system to decide when and how to be run, and can use combined resources from several computers of the cluster, which are often called computing nodes. The number of computing nodes could be from 1 up to thousands, and the cluster can increase its computing power just by incorporating new computing nodes. Most of the [World's top ranking HPC](https://www.top500.org/lists/top500/list/2022/06/) (High Performance Computing) systems follow this basic structure.


## Our cluster

> **The infrastructure** you will be using runs on a public cloud provider, and has been set up as a small computing cluster as described above. We will use it in our course mainly **to exemplify a computing cluster environment and allow the participants to familiarize with it**, given that many of the computer programs that will be used during the course are often run on computing clusters in real-case scenarios. Because of this, please keep in mind that ours is just a '***lab***' computing cluster. This means that **the computing resources are limited and the environment is simplified in some aspects compared to a typical real production-level computing cluster**. 

    
## Accessing the cluster

You can access the cluster by connecting to the login node.You can only access the cluster using SSH and other related secure protocols. You will receive by email the information about the name/address of this login node and your login credentials (username and password).

* **Command line access**: Needed for launching jobs on the cluster. Use any relatively recent ssh client:
    - **Linux**: open any terminal program and type (replace '*user*' and '*loginnodeaddress*' accordingly):
 
    `> ssh user@loginnodeaddress`
    
    - **Windows**:
        - Using **[WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux)** (Windows Subsystem for Linux): open the command line application an type the previous ssh command)
        - Using alternative free programs:
            - **[PuTTy](https://www.putty.org/)**
            - **[MobaXTerm](https://mobaxterm.mobatek.net/download-home-edition.html)**
    - **Mac**: open the terminal app and type the same ssh command
            
* **File transfer**:
    - **Graphical applications**:
        - **[FileZilla client](https://filezilla-project.org/download.php?type=client)** is one of the easiest to install and run, and works on Linux, Mac and Windows (it's free too).
        - **[MobaXTerm](https://mobaxterm.mobatek.net/download-home-edition.html)**
    - **Command line**:
        - Once you have a working ssh connection established (see above), there are several commands you can use for file transfer with similar syntax to ssh (scp for instance). Your instructor will give you specific indications if needed. 
        

## Using the cluster

After accessing the system, we will be submitting our computing tasks to a specific software that manages the cluster, schedules our tasks to be run on specific nodes according to several criteria (like the resources available on the computing nodes, permissions, priorities, etc), supervises its execution on the allocated nodes once there are enough available resources, and manages its termination (either succesful or not). In our case we will be using **[SLURM](https://slurm.schedmd.com/slurm.html)**, and the computing tasks that we submit to it are called jobs. 

Before submitting our jobs to SLURM, there may be previous steps for setting up the data we will use and the environment of the software we need to run in the jobs. Some of them are quite common, like loading software modules if our cluster supports it (although they are not mandatory and will depend on any specific case).

SLURM is a powerful software that gives users and administrators a lot capabilities and options, so there are many ways of setting up and controlling our jobs in SLURM and several of them can get rather complex. However, since in our particular case the goal is to exemplify the process in a relatively straightforward fashion, we will focus on a simple process, that we describe next:

### Getting the data for your job

Depending on each case, the data needed for your jobs may need to be downloaded from another internet site (like our GitHub repository), or uploaded from your personal computer, or it may already be present in a specific location of the cluster. Your instructor will specify how to get or reach the data you will need for each of your jobs. 

### Loading environment modules

**[Environment Modules](https://en.wikipedia.org/wiki/Environment_Modules_(software))** are tools to automate the configuration of your environment for running a program or programs. Keep in mind that on clusters and other shared systems, the amount of software installed can increase quickly, and each software may need specific environment variables to run correctly. The modules tools make the retrieval and loadind/unloading of these settings simpler for the users.

Each specific software that you run in a SLURM job may need one or mode modules to be loaded. Use the '*module load*' command to **load a module**:

```bash
> module load fastqc
```

You can **search for modules** with the '*module spider*' command (you can search for partial names):

```bash
> module spider trimmo

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  trimmomatic: trimmomatic/0.39-zeq7myp
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    This module can be loaded directly: module load trimmomatic/0.39-zeq7myp

    Help:
      A flexible read trimming tool for Illumina NGS data.

```


You can use '*module list*' to *get a list of your currently loaded modules*:

```bash
> module list

Currently Loaded Modules:
  1) autotools                   5) libiconv/1.16-de7moqn    9) ncurses/6.2-jbju4f4  13) openmpi4/4.1.1               17) bzip2/1.0.8-lswvzm3   21) fastqc/0.11.9-o7cpx62
  2) prun/2.2                    6) xz/5.2.5-3jrgbz6        10) hwloc/2.7.1-ie472nf  14) ohpc                         18) readline/8.1-if5qfy5
  3) gnu9/9.4.0                  7) zlib/1.2.12-ljlym6n     11) ucx/1.11.2           15) openjdk/11.0.15_10-lctc4i2   19) gdbm/1.19-3tckisz
  4) libpciaccess/0.16-kufodxm   8) libxml2/2.9.13-sgjxnxo  12) libfabric/1.13.0     16) berkeley-db/18.1.40-sxhhmz4  20) perl/5.34.1-dpunoxt
```

To **unload a module** use '*module unload*':

```bash
> module unload fastqc
```

Sometimes you would need to unload all your loaded modules to ensure there are no module incompatibilities before loading a specific module. You can **unload all your loaded modules** with '*module purge*'.

Please keep in mind that loading modules is a common step before submitting your jobs to SLURM, but NOT mandatory. Some software available on the cluster may not need this previous configuration. Your instructors will specify what to do for any specific case.

### Submitting jobs to SLURM

You can submit jobs to SLURM directly with a command, or use the batch mode which needs a batch file (a shell script actually) where you specify all the SLURM options, variables, and the specific command (or commands) for the job you are submitting. We will use this second method, as it's easier for specifying all the needed components for a job, and thus less prone to errors.

#### **Batch (run) files**

Next you can see an example of a batch file for submitting a job that submits *Orthofinder* to SLURM:

```bash
!/bin/bash                                                                                                             

##This is a script to run orthofinder                                                        

#SBATCH -p normal                                                                                      
#SBATCH -c 8                                                                                                     
#SBATCH --mem=6GB                                                                                       
#SBATCH --job-name orthofinder-job01                                                                            
#SBATCH -o %j.out                                                          
#SBATCH -e %j.err                                                      

#module loadding. Check available modules with `module avail` 
module load orthofinder

#running orthofinder
orthofinder -f proteomes -og -t 8

```

Some relevant bits in this script:

* '**-p**' : SLURM option to set the partition (roughly equivalent to queue with other similar software)
* '**-c**' :  SLURM option to set the CPUs (in this case threads, simultaneous tasks inside the same program process) that SLURM will allocate. It must be bigger or equal to the particular value set for the specific command line of the program you are going to run in this job.
* '**--job-name**' : SLURM option to set a name for the job, so you can identify it more easily.
* '**--mem**' : SLURM option to allocate an amout of RAM for the job.
* '**-o**' : File for the the job's output messages. Note that '***%***' translates to the number identifier of the job, which is an integer that will be showed to you when you submit the job
* '**-e**' : File for the job's error messages (if any). Note that '%' is used again.

Next, you have other commands that may be needed to run correctly the program specified below. In this case we load the needed modules directly from this script. For instance. We could also set other needed environment variables, create directories for the results (**mkdir**), create files with some content for input or configuration of the program, remove temporary files that could have been left from other runs...

And finally, we have the actual main command we want to run, in this case:

```bash
orthofinder -f proteomes -og -t 8
```

Note that the specific parameter of this program to specify the number of execution threads to run (in this case **-t**), matches the value of the SLURM option '**-c**'.

This is just a simple example, but there are many other SLURM options that can be set to modify job submissions, and several of them depend on the kind of the program to be run (what kind of parallelization it supports for example). Again, your instructors will explain other relevant aspects about SLURM job files, on a case by case basis.


#### **Submitting your job**:

To **submit a job** based on a batch file, the '***sbatch***' command is used, followed by the name (including path if needed) of the run script:

```bash
> sbatch orthofinder.run
```

After hitting enter SLURM will show you the job ID.

#### **Checking your submitted jobs**

Once the job is submitted, you **can check its state in several ways**. Just 2 examples:

* **default simple *squeue* query**:

```bash
> squeue

JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
   92    normal orthofin     alex  R       0:14      1 ip-192-168-1-101
```
If there are many jobs on this list, you can show only your user's jobs adding the '***-u***' parameter followed by your user name.

* **complex *squeue* query** (we can specify a list of individual information fields and format to be displayed), including the **CPUs allocated** for the job:

```bash
>squeue -o"%.7i %.9P %.8j %.8u %.2t %.10M %.6D %C"

JOBID PARTITION     NAME     USER ST       TIME  NODES CPUS
    94    normal orthofin     alex  R       0:02      1 8
```

Note that the state of the jobs in these cases is 'R', which means running. If the sate is 'PD', it would mean that there are no resources available on any node to run the job. SLURM will keep it waiting until there are resources available, and then run it.

You can **check the state of the computing nodes** with a query like this for example:

```bash
> sinfo -lNe

NODELIST          NODES PARTITION       STATE CPUS    S:C:T MEMORY TMP_DISK WEIGHT AVAIL_FE REASON              
ip-192-168-1-101      1   normal*       mixed 96     96:1:1 185900        0      1   (null) none                
ip-192-168-1-102      1   normal*        idle 96     96:1:1 185900        0      1   (null) none
```


## General Usage Remarks:

* Please **don't run computing tasks directly on the login node**. Once you log in to this node, use it only for ordinary commands (*ls*, *cd*, *cat*, *cp*, ...), and for other tasks that your instructor may indicate.
* Computing clusters are not storage-dedicated systems:
    - Use only the **storage space amount you need for your currently running and recently finished computing tasks**. Please don't leave data that you don't need for current computations for a long time in the cluster storage (transfer whatever data you need to other storage outside the cluster and then delete the data in the cluster).
    - Since this is a '*lab*' cluster, **there is no backup system of any kind**. You are responsible for keeping your data safe.
* The cluster **can only run command line tasks** (no graphical interface programs)
* Please don't use the cluster for running tasks that can be run easily and/or conveniently on your regular computer. A cluster is not a substitute for your computer, it's a computing dedicated system. 

