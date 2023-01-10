#!/bin/bash
# Data is originally from Abrahams et al. 2019

mkdir ENV2_Prior ENV2 ENV3_Prior ENV3

echo Running HIVtree. The MCMCs may take a few minutes to run.
cd ENV2_Prior
../../../HIVtree ../ENV2_Prior.ctl &> output & 

cd ../ENV3_Prior
../../../HIVtree ../ENV3_Prior.ctl &> output & 

cd ../ENV2
../../../HIVtree ../ENV2.ctl &> output & 

cd ../ENV3
../../../HIVtree ../ENV3.ctl &> output & 

cd ../

wait;

echo Finished running MCMCs. 

echo Parsing MCMCs and preparing files to combine estimates.
../.././parseMCMC.sh sequences.csv

echo 
echo Combining estimates across regions.
# Example of how to run the output files

echo Running combineEstimates.R with C1C2_W19_QVOA_3921.txt
Rscript ../../combineEstimates.R -m C1C2_W19_QVOA_3921.txt -s 0 -b 3.921 -t 1000 -l 3921 -g 2

echo
echo Running combineEstimates.R with C1C2_W14_QVOA_3921.txt
Rscript ../../combineEstimates.R -m C1C2_W14_QVOA_3921.txt -s 0 -b 3.921 -t 1000 -l 3921 -g 2
wait;
