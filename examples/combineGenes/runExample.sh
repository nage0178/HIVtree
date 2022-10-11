#!/bin/bash
# Data is originally from Abrahams et al. 2019

echo Parsing MCMCs and preparing files to combine estimates.
../.././parseMCMC.sh sequences.csv

echo 
echo Combining estimates across regions
# Example of how to run one of the output files
Rscript ../../combineEstimates.R C1C2_W19_QVOA_3921.txt 0 3.921 1000 3921 2

wait;
