#!/bin/bash
# Data is originally from Abrahams et al. 2019

../.././parseMCMC.sh sequences257Short.txt

# Example of how to run one of the output files
Rscript ../../combineEstimates.R CAP257_ENV_C1C2_W9_QVOA_1_3921.txt 0 3.921 1000 3921 2

wait;
