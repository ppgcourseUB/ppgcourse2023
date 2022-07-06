# Divergence Time Estimation based on Molecular clocks
# Tutorial in a nutshell


Instructor: **Jesus Lozano-Fernandez**

We are going to implement a **NODE CALIBRATED** divergence time estimation and a **TOTAL EVIDENCE ANALYSIS** in the context of bear evolution; and to compare both outputs to infer the age of the MRCA of Ursinae, and assess the impact of these different methods, models and priors. For coherence with the Bayesian class, we are going to use also [**RevBayes**](https://revbayes.github.io), which is a flexible software with his own programming language similar to R, named Rev. Is out of the scope of this short tutorial to be able to fully understand the language, the models or the parameters. At the end of the class, you are expected to have a grasp of the main difference between methods, the utility of fossil and stratigraphic calibrations and the relevance of the prior assumptions. This tutorial is based on material available for learning divergence time estimatin using RevBayes, available at [https://revbayes.github.io/tutorials/dating/](https://revbayes.github.io/tutorials/dating/). All the scripting to run the analyss and the datasets has been already built and are in place, so you won't need to modify them, just running the analyses. Morevoer, it has been uploaded also the output of the analyses which were previously ran, on case you face some difficulties. We will, though, examine the inpt files and run the analyes, and learn about the interpretation of results (using the softwares [**Tracer**](http://tree.bio.ed.ac.uk/software/tracer/), and [**FigTree**](http://tree.bio.ed.ac.uk/software/figtree/)).

 ![question](img/MolecularDating.jpg)

***

## NODE DATING

* *Based on the 3rd exercise ‘Estimating speciation times using node dating’*
https://revbayes.github.io/tutorials/dating/nodedate
* Additionally, RevBayes uses a scripting language to call functions and provide arguments: the `Rev` language. It is similar to `R` but more heavily scripted.


## Background

First, create a folder on your home directory and name it `lab_1`
Copy the `Intro.Rev` file from the pggcourse lab materials folder into the `lab_1` folder. 
Copy the `myScript.sh` file from the pggcourse lab materials folder into the `lab_1` folder.
 
```
mkdir lab_1
cp -p Intro.Rev lab_1/
cp -p myScript.sh lab1/
```
