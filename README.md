# HIVtree
HIVtree is Bayesian phylogenetic inference program which infers latent integration times of HIV sequences and internal node ages, conditional on a tree topology.
It was originally modified from PAML version 4.9 by Ziheng Yang.

## Download and Compile the Program
To clone the program from GitHub 

```
git clone https://github.com/nage0178/HIVtree.git
```

To compile the program, 

```
cd ~/HIVtree/
make
```

## Running the HIVtree
HIVtree is run the same way mcmctree. 
For a full explanation of options and the control file, please see the manual. 

```
./mcmctree control.ctl
```

