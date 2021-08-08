#!/usr/bin/env bash

gfortran -o xpara parametersGenerator.F90
./xpara
rm xpara

################# model ##############
cd ../model/
gfortran -o xvel velocity3d_true.F90
./xvel
rm xvel
cd ../commandCenter/

################## data ##############
cd ../inversion/
gfortran -o xsyn syntheticTraveltimes.f90
./xsyn
rm xsyn velmodel.mod

cd ../data
cp traveltimeReceiverGathers_synthetic traveltimeReceiverGathers
cd ../commandCenter/
