# Bayesian Inference in RevBayes
# Lab 1: Introduction to RevBayes


Instructor: **Isabel Sanmartin**

In this session, I will introduce you to [RevBayes](https://revbayes.github.io), an interactive environment for all types of Bayesian Inference analyses.

You can download RevBayes for your own use from [Downloads](https://revbayes.github.io/download). This includes compiled executables for Windows and MacOSX, and a singularity image for Linux. Users can also compile from the source code, which might be a good idea since this is the most up to date.
The RevBayes website also posts incredibly detailed [Tutorials](https://revbayes.github.io/tutorials/), including the theory behind the methods. READ THEM!

## Software description and hints

* RevBayes is a "command line" program, which means that it does not have a graphic user interface, and we must interact with it through the Terminal/Console.
* Mathematical and operational functions in RevBayes are written in the object-oriented programming language C++. This makes RevBayes computationally efficient for even the most complex analyses.
* Additionally, RevBayes uses a scripting language to call functions and provide arguments: the `Rev` language. It is similar to `R` but more heavily scripted.
* RevBayes is a multi-purpose software, meaning that, in addition to phylogenetics, it can be used for molecular dating, trait evolution, diversification analysis, biogeography, epidemiology, etc.


## Data description and access from the cluster

We will be working in the cluster *interactively*

First, access the data files and scripts for the practice. Access the lab cluster via the Terminal in MacOSX and Linux, or the Console in Windows machine, with the `ssh`, as you learnt at the *Intro to our cloud* practice.

```
ssh user@ec2-34-242-61-70.eu-west-1.compute.amazonaws.com
```
where the *user* are your credentials. It will ask for your password. Once inside, clone the entire folder *Bayesian_Inference_ISABEL_SANMARTIN* with the `svn` command. This will copy all files and folders in your home directory within the cluster.

```
svn export https://github.com/ppgcourseUB/ppgcourse2023/trunk//Bayesian_Inference_ISABEL_SANMARTIN
```

Check the contents of the folder. There are two subfolders named `Lab`. We will start with `Lab_1`, which introduces the RevBayes software and the `Rev` language. To see the contents of the folder, you need to move with the `cd` command

```
cd Lab_1
ls
```
You would see that there are several files, including scripts and data we need for the practice: `Intro.Rev`, `myScript.sh`, and `myScript-cluster.sh`

## Data description and access from your computer

Alternatively, if you want to run RevBayes from your own computer, first, you need to download the software in the version that is appropriate for your operating system (see *Downloads* above).
Next, access the pggcourse lab materials in Github [ppgcourse2022] (https://github.com/ppgcourseUB/ppgcourse2023), and download the scripts and files for the practice.
OBS!! The scripts and  files should be downloaded into the folder containing the *RevBayes* executable (the `rb` binary in the Mac). Inside this folder, create a new folder `lab_1` and copy the `Intro.Rev` and `myScript.sh` files into this folder.
 
```
mkdir lab_1
cp -p Intro.Rev lab_1/
cp -p myScript.sh lab1/
```

## Launching RevBayes from the cluster
To launch RevBayes using the version already installed in the cluster, first load the module, using the commands you learnt in the *Intro* class:

```
module load revbayes/1.1.1-zjzfb6s
```
Then, launch RevBayes, type:

```
rb-mpi
```

(notice that the version installed in the cluster is the MPI version (message passage interface), which allows running the software using multiple threads or processors)

This should launch RevBayes and give you a command prompt (the `>` character); this means RevBayes is waiting for input. 


## Launching RevBayes from your computer
Navigate to your `lab_1` directory. In Windows, simply double click on the executable. In MacOSX, launch RevBayes by typing `rb` into the command line (OBS: Because `rb` is a program (and unless you have set up the path), you need to launch it using the command below, which will ask for your administrator password:
```
./rb
```
This should launch RevBayes and give you a command prompt (the `>` character); this means RevBayes is waiting for input. 
The _working directory_ is the directory that RevBayes is currently working in. When you tell RevBayes to look up a file in a particular path, the path you provide is interpreted _relative_ to the working directory. You can print the current working directory using the `getwd()` command.


## Learning the basic structure of the `Rev` language.

## Mathematical operations
`Rev` is an interpreted language for statistical computing and phylogenetic analysis. Therefore, the basics are simple mathematical operations. Entering each of the following lines into the RevBayes prompt will execute these operations. You can also execute multiple operations in one line if you separate them with a semicolon.

```
# Simple mathematical operators:
1 + 1             # Addition
5 * 5             # Multiplication
2^3               # Exponentiation

1 + 3; 4 / 2 # Multiple statements in one line
```

## Functions

Functions are commands that perform more complex procedures than the above operations. Notice that RevBayes is case-sensitive, so `exp(1)` will work but `Exp(1)` will give you an error.
```
# Math functions:
exp(1)            # exponential function
ln(1)             # natural logarithmic function
sqrt(16)          # square root function
power(2,2)        # power function: power(a,b) = a^b
```
## Variables
One of the most important features of RevBayes is the ability to declare and assign variables. There are three types of variables, called _constant_, _deterministic_, and _stochastic_ variables. Variables are also the "nodes" in the directed acyclical graphs that are used to create the RevBayes models.

_Constant variables_ contain values that adopt fixed values. The left arrow (`<-`) creates a constant variable and automatically assign the following value to it. Here, we create the constant node `a` and assign a value of `1`. We can print the value by typing `a` and pressing enter.

```
# variable assignment constant
a <- 1 # assignment of constant variable 'a'
a      # printing the value of 'a'
```

_Deterministic variables_  are variables whose value depends on another random variable: the value changes when the variable they depend on change via a function or transformation. They are created with the colon-equal assignment (`:=`). Here, we create a deterministic variable, `b`, whose values are dependent on the value assigned to `a` via the exponential function `exp(a)`.

```
# Variable assignment: deterministic
b := exp(a)
b # printing the value of 'b' ("2.718282")
# Assigning another value to 'a' changes the value of 'b'
a <- 2
b # printing the new value of 'b' ("7.389056")
```

_Stochastic variables_ are random variables whose values are drawn from a statistical distribution with its own parameters. Because they are random, values will change during the analysis. They are created with the tilde assignment (`~`). Here, we create a stochastic variable `x` that is drawn from an exponential distribution with rate parameter (`lambda`).

```
# Variable assignment: stochastic
# First, we assign a value to the lambda parameter governing the exponential distribution
lambda <- 1.0
# Next, we create the stochastic variable with values drawn from this exponential distribution
x ~ dnExponential(lambda)
x # print value of stochastic node 'x' ("1.256852") 
```

## Vectors

Vectors are containers that contain multiple variables of the same type. To create a vector with values

```
Z <- v(1.0,2.0,3.0) # create a vector "implicitly"
Z # print the vector ("[ 1.000, 2.000, 3.000 ]")

# or alternatively fill an empty vector one by one "explicitly
Z[1] <- 1.0 # make the first element
Z[2] <- 2.0 # make the second element
Z[3] <- 3.0 # make the third element
Z # print the vector ("[ 1.000, 2.000, 3.000 ]")
```

## `for` Loops
`for` loops are important programming structures that allow you to repeat the same statement a number of times on different variables. This simple `for` loop creates the variable `i`, and for each value of `i` from 1 to 100, prints the value of `i` to the screen ("1, 2, 3, 4... 100").

```
for (i in 1:100) {
  i
}
```

## Quitting RevBayes
When we're done with RevBayes, or want to relaunch the program, we can quit using the `q()` command:

```
q() # quitting RevBayes

```

## Scripts

So far we've been using RevBayes *interactively*: by typing commands in line-by-line. Most often, however, we use *scripts*: a text file that contains a sequence of commands for the program to execute.

You can *source* the contents of a script from RevBayes using the `source("name of file")` command (the quotation marks are critical!). Source the `Intro.Rev` script that contains all commands above (OBS: First, you need to launch RevBayes again with the command `rb-mpi` (cluster) or `rb` (your Terminal).

```
source("Intro.Rev")
```

Alternatively, you can run the script file from the Terminal directly (outside RevBayes) using the command `rb`. OBS: In this case, we don't need the quotation marks!)

```
rb-mpi Intro.Rev
```

Or you can run it from a bash script from your Terminal. You can find an example of such a script on myScript.sh

```
#!/bin/bash
rb-mpi Intro.Rev
exit
```
And you can run it, using
```
bash myScript.sh
```

## Launching RevBayes with a bash file

Finally, we can use a bash file if we are running the script within a cluster. Below, is an example of a bash file to run the `Intro.Rev` script in RevBayes. Because we are running in a cluster, we need to include a command for the output to be saved. You can find this script `myScript-cluster.sh` in the `lab_1` folder.

```
#!/bin/bash                                                                                                             

#SBATCH -p normal                                                                                                       
#SBATCH -n 8                                                                                                            
#SBATCH -c 1                                                                                                            
#SBATCH --mem=6GB                                                                                                       
#SBATCH --job-name orthofinder-job01                                                                           \        
#SBATCH -o %j.out                                                                                                       
#SBATCH -e %j.err                                                                                                       

module load revbayes

mpirun -np 8 rb-mpi Intro.Rev
```
And to submit the job to SLURM

```
sbatch myScript-cluster.sh
```
And to check the progress

```
squeue
```
The results are inside the *[filename].out* file. You can read it with the command `less`, or with a text processor program such as VIM or AWK

## Exercises

Exit and restart RevBayes. Create a fresh, blank script in a text editor. You can use Vim or AWK, as you learnt the first day. Alternatively, you can type directly at the RevBayes prompt.

Create a variable called `z` with the value `10`. What kind of variable is this?

Create a second variable `y` which is `y := ln(z)`. What kind of variable is `y`? What is its value?

Change the value of `z` to `100`. Before printing `y`, can you guess if the value will be lower or higher?

Write a `for` loop that creates  a variable `i`, and for each value of `i` from 1 to 100, creates a second variable `z`, which is a deterministic function of `f(i) = 2 * i`, and then prints all the values of `z` to the screen. Then, using the `mean()` function, calculate the mean of those numbers.
