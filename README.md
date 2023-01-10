# HIVtree
HIVtree is Bayesian phylogenetic inference program which infers latent integration times of HIV sequences and internal node ages, conditional on a tree topology.
It was originally modified from PAML version 4.9 by Ziheng Yang.

## Download and Compile the Program
To clone the program from GitHub 

```
git clone https://github.com/nage0178/HIVtree.git
```

To compile the program, navigate to the directory and use the make command.

```
cd HIVtree/
make
```
You must have a C compiler, such as gcc or clang, installed to compile the program. 
To use the make command, you must have the make system installed.

## Running HIVtree
HIVtree is run in the same way as mcmctree. 
For a full explanation of options and the control file, please see the manual. 

```
./HIVtree control.ctl
```
Examples are provided in the example directory which includes one example data set with a single gene and another with multiple genes illustrating how to combine inferences across multiple genes.

## R packages
A recent version of R must be installed with the following packages which are necessary to run the combined analysis. 

```
kdensity
GoFKernel
getopt
```
