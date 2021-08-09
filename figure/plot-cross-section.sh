#!/usr/bin/env bash

gfortran -o zCutVelocity zCutVelocity.f90
./zCutVelocity
rm zCutVelocity

gfortran -o zCutVelocity_true zCutVelocity_true.f90
./zCutVelocity_true
rm zCutVelocity_true

./vcut.gmt
