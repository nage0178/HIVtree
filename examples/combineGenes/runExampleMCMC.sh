#!/bin/bash
# Data is originally from Abrahams et al. 2019

mkdir 257_ENV2_d10_0.75Prior_1 257_ENV2_d10_0.75_1 257_ENV3_d10_0.75Prior_1 257_ENV3_d10_0.75_1

cd 257_ENV2_d10_0.75Prior_1
../../.././HIVLateTree ../257_ENV2_d10_0.75Prior_1.ctl &> output & 


cd ../257_ENV3_d10_0.75Prior_1
../../.././HIVLateTree ../257_ENV3_d10_0.75Prior_1.ctl &> output & 

cd ../257_ENV2_d10_0.75_1
../../.././HIVLateTree ../257_ENV2_d10_0.75_1.ctl &> output & 

cd ../257_ENV3_d10_0.75_1
../../.././HIVLateTree ../257_ENV3_d10_0.75_1.ctl &> output & 

wait;
