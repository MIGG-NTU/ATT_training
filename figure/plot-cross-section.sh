#!/usr/bin/env bash

gfortran -o zCutVelocity zCutVelocity.f90
./zCutVelocity
rm zCutVelocity

./vcut.gmt
