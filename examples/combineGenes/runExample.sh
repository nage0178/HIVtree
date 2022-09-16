#!/bin/bash
# Data is originally from Abrahams et al. 2019

mkdir 257_ENV2_d10_0.75Prior_1 257_ENV2_d10_0.75_1 257_ENV3_d10_0.75Prior_1 257_ENV3_d10_0.75_1

cd 257_ENV2_d10_0.75Prior_1
../.././HIVtree 257_ENV2_d10_0.75Prior_1.ctl &> output & 


cd ../257_ENV3_d10_0.75Prior_1
../.././HIVtree 257_ENV3_d10_0.75Prior_1.ctl &> output & 

cd ../257_ENV2_d10_0.75_1
../.././HIVtree 257_ENV2_d10_0.75_1.ctl &> output & 

cd ../257_ENV2_d10_0.75_1
../.././HIVtree 257_ENV2_d10_0.75_1.ctl &> output & 

cd ../
../.././parseMCMC.sh sequences257.txt

# Example of how to run one of the output files
Rscript ../../combineEstimates.R CAP257_ENV_C1C2_W9_QVOA_1_3921.txt 0 3.921 1000 3921 2

wait;
