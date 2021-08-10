#!/usr/bin/env bash

# 
# Prepare cross-section data
#
gfortran -o zCutVelocity zCutVelocity.f90
./zCutVelocity
rm zCutVelocity

gfortran -o zCutVelocity_true zCutVelocity_true.f90
./zCutVelocity_true
rm zCutVelocity_true

#
# Plot velocity model and kernel
#

gmt begin vcut-example pdf
gmt set FONT_LABEL 10p
gmt subplot begin 3x1 -Fs12c/2.0c -A -SCb -Bxa20.0f10.0+l"Distance (km)" -Bya10.0f5.0+l"Depth (km)" -BWSen

# Plot velocity model (a)
gmt subplot set 0
gmt basemap -R0.0/240.0/0.0/36.0 -JX?/-? 
awk ' { print $1,$2,$5 }' dvelCut001 | gmt surface -R0.0/240.0/0.0/36.0 -I0.2 -Gtomo.grd
T=$(gmt grdinfo -T+s tomo.grd)
gmt makecpt -Cpsrgb.cpt $T -I
gmt grdimage tomo.grd
gmt colorbar -DJRM+w1.8c/0.2c+e -Ba+u"%"


# Plot kernel (b)
gmt subplot set 1
gmt basemap -R0.0/240.0/0.0/36.0 -JX?/-?
awk ' { print $1,-1.0*$3,$4 }' ../inversion/kernel | gmt surface -R0.0/240.0/0.0/36.0 -I0.2 -Gtomo.grd
T=$(gmt grdinfo -T+a5+s tomo.grd)
gmt makecpt -Cpsrgb.cpt $T -I
gmt grdimage tomo.grd
gmt colorbar -DJRM+w1.8c/0.2c+e -Ba

# Plot true velocity model (a)
gmt subplot set 2
gmt basemap -R0.0/240.0/0.0/36.0 -JX?/-? 
awk ' { print $1,$2,$5 }' dveltrueCut001 | gmt surface -R0.0/240.0/0.0/36.0 -I0.2 -Gtomo.grd
T=$(gmt grdinfo -T+s tomo.grd)
gmt makecpt -Cpsrgb.cpt $T -I
gmt grdimage tomo.grd
gmt colorbar -DJRM+w1.8c/0.2c+e -Ba+u"%"
awk '{print $3, $2}' ../data/receivers | gmt plot -Si0.3c -Gblue -N
awk '{print $3, -$2}' ../data/sources | gmt plot -Sa0.2c -Gred -W0.5p,white

gmt subplot end
gmt end show
