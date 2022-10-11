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
cd ~/HIVtree/
make
```
You must have a C compiler, such as gcc or clang, installed to compile the program. 
To use the make command, you must have the make system.

## Running the HIVtree
HIVtree is run the same way mcmctree. 
For a full explanation of options and the control file, please see the manual. 

```
./HIVtree control.ctl
```
Examples are provided in the example directory.
There is one example data set with a single gene and an example of how to combine inferences across multiple genes.

## R packages
The following packages are necessary to run the combined analysis. 

```
kdensity
GoFKernel
```
