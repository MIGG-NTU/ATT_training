#!/usr/bin/env bash

rm -vf ./region

cd ../data/
rm -vf evtstaminmax nevtstaray receivers sources traveltimeReceiverGathers traveltimeReceiverGathers_synthetic

cd ./dataSelection/1_picks/
rm -vf largeEarthquakes  obstime-rec-gather  obstime-src-gather  phaseDataInResearchRegion
cd ../../

cd ../figure/
rm -vf dvelCut* dveltrueCut001 mesh.mod tomo.grd vz1d vcut-example.pdf

cd ../inversion/
rm -vf adjtfield* kernel residual0* velocity3d0*
rm -vf objectiveFunction

cd ../mesh/
rm -vf mesh_cube0* multiple-grid

cd ../model/
rm -vf velocity3d velocity3d_true velocity2d_true
