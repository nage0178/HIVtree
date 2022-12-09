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

wait;
