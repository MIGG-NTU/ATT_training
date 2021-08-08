#!/usr/bin/env bash

gfortran -o xpara parametersGenerator.F90
./xpara
rm xpara

##################### model ##################
cd ../model/
gfortran -o xvel velocity3d.F90
./xvel
rm xvel
cd ../commandCenter/

###################### mesh ##################
cd ../mesh/
gfortran -o xmesh meshgenerator.F90
./xmesh
rm xmesh regmesh.mod
cd ../commandCenter/

###################### data ###################
cd ../inversion/
gfortran -o xinv inversion.F90
./xinv
rm xinv parameter_inversion mesh.mod velmodel.mod

cd ../commandCenter/
rm gridnumber
