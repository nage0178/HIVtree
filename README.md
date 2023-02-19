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
## LICENSE
    Copyright (C) 2022 Anna Nagel
    Originally modified from PAML version 4.9 by Ziheng Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
